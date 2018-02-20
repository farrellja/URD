#' Build Tree
#' 
#' @importFrom stats weighted.mean
#' 
#' @param object An URD object
#' @param pseudotime (Character) Pseudotime to use for building tree
#' @param tips.use (Character vector) Name of tips to use in the final tree (Default is NULL, which uses all tips.)
#' @param divergence.method (Character: "ks" or "preference") Test to use to determine whether visitation has diverged for each pseudotime window.
#' @param cells.per.pseudotime.bin (Numeric) Approximate number of cells to assign to each pseudotime bin for branchpoint finding.
#' @param bins.per.pseudotime.window (Numeric) Width of moving window in pseudotime used for branchpoint finding, in terms of bins.
#' @param minimum.visits (Numeric) Minimum number of random walk visits to a cell to retain it in the tree
#' @param visit.threshold (Numeric) Cells are considered potential members for segments/tips from which random walks visited them 
#' at least this fraction of their maximum visitation from a single tip
#' @param save.breakpoint.plots (Path) Path to save plots summarizing (default is NULL, which does not save plots as they are somewhat slow)
#' @param p.thresh (Numeric) p-value threshold to use in determining whether visitation is significantly different from pairs of tips
#' @param min.cells.per.segment (Numeric) Segments with fewer assigned cells weill be collapsed during tree construction
#' @param min.pseudotime.per.segment (Numeric) Segments shorter than this in pseudotime will be collapsed during tree construction
#' @param dendro.node.size (Numeric) Number of cells to assign per node (used for averaging expression for tree dendrogram branch coloring)
#' @param dendro.cell.jitter (Numeric) For the dendrogram tree layout, how much jitter to apply to cells. This can be revised after building the tree by re-generating the cell layout using \code{\link{treeLayoutCells}} with a different parameter.
#' @param dendro.cell.dist.to.tree (Numeric) For the dendrogram tree layout, how far to push cells away from the dendrogram branches. This can be revised after building the tree by re-generating the cell layout using \code{\link{treeLayoutCells}} with a different parameter.
#' @param verbose (Logical) Report on progress?
#' @return An URD object with an URD-recovered tree structure stored in \code{@@tree}
#' @export
buildTree <- function(object, pseudotime, tips.use=NULL, divergence.method=c("ks", "preference"), weighted.fusion=T, use.only.original.tips=T, cells.per.pseudotime.bin=80, bins.per.pseudotime.window=5, minimum.visits=10, visit.threshold=0.7, save.breakpoint.plots=NULL, p.thresh=.01, min.cells.per.segment=1, min.pseudotime.per.segment=.01, dendro.node.size=100, dendro.cell.jitter=0.15, dendro.cell.dist.to.tree=0.05, verbose=T) {
  # Check divergence.method parameter
  if (length(divergence.method) > 1) divergence.method <- divergence.method[1]
  if (!(divergence.method %in% c("ks", "preference"))) stop("Divergence method must be 'ks' or 'preference'.")
  # Make sure tips is a character vector
  tips <- as.character(tips.use)
  # Record tips for posterity
  object@tree$tips <- tips
  # If you're going to do balanced fusion, calculate number of cells per tip.
  if (weighted.fusion) {
    tip.size <- unlist(lapply(object@tree$cells.in.tip, length))
    if (length(tip.size)==0) stop("Either weighted.fusion must be false, or loadTipCells must be run prior.")
  }
  # Load the pseudotime used into the tree object for use later with plotting
  object@tree$pseudotime <- object@pseudotime[,pseudotime]
  names(object@tree$pseudotime) <- rownames(object@pseudotime)
  # Initialize segment to add
  seg.add <- as.character(max(as.numeric(tips)) + 1)
  # Grab cells in segments
  cells.in.segments <- putativeCellsInSegment(object, tips, minimum.visits=minimum.visits, visit.threshold=visit.threshold) 
  object@tree$cells.in.segments <- cells.in.segments
  # Initialize pseudotime start & end data.frame
  object@tree$segment.pseudotime.limits <- data.frame(
    start=sapply(tips, function(segment) min(object@pseudotime[rownames(cells.in.segments)[which(cells.in.segments[,segment])],pseudotime])),
    end=sapply(tips, function(segment) max(object@pseudotime[rownames(cells.in.segments)[which(cells.in.segments[,segment])],pseudotime])),
    stringsAsFactors = F, row.names = tips
  )
  # Calculate segment pseudotime divergence, WITHOUT cache in case this function is being re-run
  object <- allSegmentDivergenceByPseudotime(object, pseudotime=pseudotime, divergence.method=divergence.method, segments=tips, pseudotime.cuts=cells.per.pseudotime.bin, window.size=bins.per.pseudotime.window, minimum.visits=minimum.visits, visit.threshold=visit.threshold, p.thresh=p.thresh, breakpoint.decision.plots=save.breakpoint.plots, cache=F, verbose=verbose)
  # Loop through until everything has joined into a single root.
  while(length(tips) >= 2) {
    
    # Find the oldest pseudotime breakpoint, and thereby which branches to join at which pseudotime
    fuse.id <- which.max(object@tree$segment.divergence$pseudotime.breakpoint)
    seg.1 <- object@tree$segment.divergence[fuse.id, "seg.1"]
    seg.2 <- object@tree$segment.divergence[fuse.id, "seg.2"]
    pt.break <- object@tree$segment.divergence[fuse.id, "pseudotime.breakpoint"]
    
    # Record the fusion
    if (verbose) print(paste("Joining segments", seg.1, "and", seg.2, "at pseudotime", round(pt.break, digits=3), "to create segment", seg.add))
    if (is.null(object@tree$segment.joins)) {
      # Make a tree structure data frame if it doesn't exist
      object@tree$segment.joins <- data.frame(child.1=seg.1, child.2=seg.2, parent=as.character(seg.add), pseudotime=pt.break, stringsAsFactors = F)
    } else {
      # Or add the row to a previously existing data frame
      object@tree$segment.joins <- rbind(object@tree$segment.joins, c(seg.1, seg.2, as.character(seg.add), pt.break))
    }
    
    # Alter the pseudotime start/end for these segments
    object@tree$segment.pseudotime.limits[seg.add,] <- c(min(object@tree$segment.pseudotime.limits[c(seg.1, seg.2),"start"]), pt.break)
    object@tree$segment.pseudotime.limits[c(seg.1,seg.2),"start"] <- pt.break
    
    # Create the new visitation data as the mean of all of the children of the two segments being joined.
    children.of.branchpoint <- unique(c(segChildrenAll(object, seg.1, include.self=T, original.joins = T, point.at.segment.joins.initial = F), segChildrenAll(object, seg.2, include.self=T, original.joins = T, point.at.segment.joins.initial = F)))
    
    # Limit children of branchpoints to original tips only.
    if (use.only.original.tips) children.of.branchpoint <- intersect(children.of.branchpoint, object@tree$tips)
    
    # Do fusion to generate new raw visitation data
    if (weighted.fusion) {
      weights <- tip.size[children.of.branchpoint]
      weights <- weights/sum(weights)
      object@diff.data[,paste0("visitfreq.raw.", seg.add)] <- apply(object@diff.data[,paste0("visitfreq.raw.", children.of.branchpoint)], 1, weighted.mean, w=weights)
    } else {
      object@diff.data[,paste0("visitfreq.raw.", seg.add)] <- apply(object@diff.data[,paste0("visitfreq.raw.", children.of.branchpoint)], 1, mean)
    }
    # Create log-transformed version also.
    object@diff.data[,paste0("visitfreq.log.", seg.add)] <- log10(object@diff.data[,paste0("visitfreq.raw.", seg.add)] + 1)
    
    # Create cells in segment membership for the new segment
    object@tree$cells.in.segments[,seg.add] <- apply(object@tree$cells.in.segments[,c(seg.1, seg.2)], 1, any)
    
    # Update tips under consideration
    tips <- setdiff(c(tips, as.character(seg.add)), c(seg.1, seg.2))
    
    # Update the divergence data
    if (length(tips) >= 2) {
      object <- allSegmentDivergenceByPseudotime(object, pseudotime=pseudotime, segments=tips, pseudotime.cuts=cells.per.pseudotime.bin, window.size=bins.per.pseudotime.window, minimum.visits=minimum.visits, visit.threshold=visit.threshold, p.thresh=p.thresh, breakpoint.decision.plots=save.breakpoint.plots, cache=T, verbose=verbose)
    }
    
    # Prepare for next join
    seg.add <- as.character(as.numeric(seg.add) + 1)
  }
  
  # Determine complete list of all segments
  object@tree$segments <- as.character(sort(as.numeric(unique(unlist(object@tree$segment.joins[,c("child.1", "child.2", "parent")])))))
  
  # Make a backup copy of segment.joins in case segments are later deleted and you need this history
  object@tree$segment.joins.initial <- object@tree$segment.joins
  
  # Assign cells to segments (for tree structure)
  if (verbose) print("Assigning cells to segments.")
  object <- assignCellsToSegments(object, pseudotime, verbose)
  
  # Delete overly short segments from the tree, and reassign segments
  object <- reformatSegmentJoins(object)
  if (verbose) print ("Collapsing short segments.")
  object <- collapseShortSegments(object, min.cells.per.segment = min.cells.per.segment, min.pseudotime.per.segment = min.pseudotime.per.segment)
  if (verbose) print ("Removing singleton segments.")
  object <- removeUnitarySegments(object)
  if (verbose) print ("Reassigning cells to segments.")
  object <- assignCellsToSegments(object, pseudotime, verbose)
  
  # Throw out cells that weren't assigned.
  object@diff.data <- object@diff.data[unlist(object@tree$cells.in.segment),]
  
  # Assign cells to nodes for plotting along trajectories
  if (verbose) print ("Assigning cells to nodes.")
  object <- assignCellsToNodes(object, node.size=dendro.node.size, pseudotime=pseudotime)
  
  # Generate tree layout
  if (verbose) print ("Laying out tree.")
  object <- treeLayoutDendrogram(object)
  object <- treeLayoutElaborate(object)
  
  # Generate cell layout
  if (verbose) print ("Adding cells to tree.")
  object <- treeLayoutCells(object=object, pseudotime=pseudotime, jitter=dendro.cell.jitter, jitter.push=dendro.cell.dist.to.tree)
  
  # Return the finished object
  return(object)
}