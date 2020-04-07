#' Build Tree
#' 
#' This function recovers the branching dendrogram structure from the developmental trajectories
#' found by the random walks. IMPORTANT: CURRENTLY THIS FUNCTION IS DESTRUCTIVE. Output the
#' results to a new URD object (i.e. \code{object2 <- buildTree(object1, ...)} or save your object
#' (using \code{\link{saveRDS}}) before running.
#' 
#' This function starts from each tip and agglomeratively joins trajectories that
#' visited the same cells (which indicates a common progenitor type). It works by
#' first comparing all tips in a pair-wise fashion. For each pair, the cells
#' visited by either tip are divided up by a moving window through pseudotime. Then,
#' a test is used to determine whether the cells in each window were visited significantly
#' differently by walks from the two tips. A putative branchpoint is chosen when
#' the test starts being significant. After comparing all tips, the latest branchpoint
#' is chosen, and the two segments are combined upstream of the branchpoint into a new
#' segment. This process is repeated iteratively until only one trajectory remains.
#' Then, the branching tree is refined to remove overly short segments, thereby
#' allowing for multi-furcating branchpoints. Finally, the dendrogram layouts are
#' generated.
#' 
#' There are several parameters that can be modified that affect the resultant tree
#' structure slightly. (1) The test used for determining whether visitation of cells
#' is different from two tips (\code{divergence.method}), the p-value threshold used (\code{p.thresh}),
#' and the number of cells in each pseudotime window (\code{cells.per.pseudotime.bin}
#' and \code{bins.per.pseudotime.window}). In general, adjusting the p-value threshold
#' will make all branchpoints slightly earlier or later. Adjusting the number of cells
#' in each window may be important to make sure that the procedure is buffered
#' from noise (too small windows can cause bizarre fusion results), but if it is
#' too large, cell populations that split near the end of your timecourse may immediately fuse.
#' 
#' A number of plots can help shed light on the tree-building process. If \code{save.breakpoint.plots}
#' points to a path, then a plot is saved every time two branches fuse showing the
#' pseudotime of the putative branchpoint between every pair of segments. The plots are
#' named by the number of remaining segments (so the largest number is the first
#' fusion event). Within a plot, each bar represents the comparison between a pair
#' of segments with pseudotime along the x-axis, red represents pseudotime windows
#' that are not significantly different and blue represents pseudotime windows that
#' are significantly different, and the orange dotted line represents the putative
#' pseudotime of the breakpoint. Additionally, branchpoint detail plots can be very
#' useful if \code{save.all.breakpoint.info=T} -- see \code{\link{branchpointDetailsVisitTsne}},
#' \code{\link{branchpointDetailsVisitDist}}, and \code{\link{branchpointDetailsPreferenceDist}}.
#' 
#' @importFrom stats weighted.mean
#' @importFrom gridExtra grid.arrange
#' 
#' @param object An URD object
#' @param pseudotime (Character) Pseudotime to use for building tree
#' @param tips.use (Character vector) Name of tips to use in the final tree (Default is NULL, which uses all tips.) Currently, these should all be numeric -- the tree building function will fail otherwise.
#' @param divergence.method (Character: "ks" or "preference") Test to use to determine whether visitation has diverged for each pseudotime window.
#' @param cells.per.pseudotime.bin (Numeric) Approximate number of cells to assign to each pseudotime bin for branchpoint finding.
#' @param bins.per.pseudotime.window (Numeric) Width of moving window in pseudotime used for branchpoint finding, in terms of bins.
#' @param minimum.visits (Numeric) Minimum number of random walk visits to a cell to consider it for the tree
#' @param visit.threshold (Numeric) Cells are considered potential members for segments/tips from which random walks visited them 
#' at least this fraction of their maximum visitation from a single tip (i.e. if \code{visit.treshold=0.7} and a cell was visited
#' most heavily from tip X with 10,000 visits, then all tips that visited that cell at least 7,000 times will include it when
#' determining branchpoints.)
#' @param save.breakpoint.plots (Path) Path to save plots summarizing (default is \code{NULL}, which does not save plots as they are somewhat slow)
#' @param save.all.breakpoint.info (Logical) Should all information about breakpoints be stored in the object for use in diagnostic plots? (Can add several hundred MB to object size, but enables branchpointDetails plots.)
#' @param p.thresh (Numeric) p-value threshold to use in determining whether visitation is significantly different from pairs of tips
#' @param min.cells.per.segment (Numeric) Segments with fewer assigned cells will be collapsed during tree construction
#' @param min.pseudotime.per.segment (Numeric) Segments shorter than this in pseudotime will be collapsed during tree construction
#' @param dendro.node.size (Numeric) Number of cells to assign per node (used for averaging expression for tree dendrogram branch coloring)
#' @param dendro.cell.jitter (Numeric) For the dendrogram tree layout, how much jitter to apply to cells. This can be revised after building the tree by re-generating the cell layout using \code{\link{treeLayoutCells}} with a different parameter.
#' @param dendro.cell.dist.to.tree (Numeric) For the dendrogram tree layout, how far to push cells away from the dendrogram branches. This can be revised after building the tree by re-generating the cell layout using \code{\link{treeLayoutCells}} with a different parameter.
#' @param verbose (Logical) Report on progress?
#' 
#' @return An URD object with an URD-recovered tree structure stored in \code{@@tree}.
#' 
#' @examples 
#' # Load the cells used for each tip into the URD object
#' axial.tree <- loadTipCells(axial, "tip.clusters")
#' 
#' # Build the tree
#' axial.tree <- buildTree(axial.tree, pseudotime = "pseudotime", tips.use=1:2, divergence.method = "preference", cells.per.pseudotime.bin = 25, bins.per.pseudotime.window = 8, save.all.breakpoint.info = T, p.thresh=0.001)
#' 
#' object.built <- buildTree(object = object, pseudotime="pseudotime", divergence.method = "ks", tips.use=tips.to.use, weighted.fusion = T, use.only.original.tips = T, cells.per.pseudotime.bin=80, bins.per.pseudotime.window = 5, minimum.visits = 1, visit.threshold = 0.7, p.thresh = 0.025, save.breakpoint.plots = NULL, dendro.node.size = 100, min.cells.per.segment = 10, min.pseudotime.per.segment = .01, verbose = F)
#' 
#' @export
buildTree <- function(object, pseudotime, tips.use=NULL, divergence.method=c("ks", "preference"), weighted.fusion=T, use.only.original.tips=T, cells.per.pseudotime.bin=80, bins.per.pseudotime.window=5, minimum.visits=10, visit.threshold=0.7, save.breakpoint.plots=NULL, save.all.breakpoint.info=F, p.thresh=.01, min.cells.per.segment=1, min.pseudotime.per.segment=.01, dendro.node.size=100, dendro.cell.jitter=0.15, dendro.cell.dist.to.tree=0.05, verbose=T) {
  # Check divergence.method parameter
  if (length(divergence.method) > 1) divergence.method <- divergence.method[1]
  if (!(divergence.method %in% c("ks", "preference"))) stop("Divergence method must be 'ks' or 'preference'.")
  # If tips were not specified, define them.
  if (is.null(tips.use)) {
    tips.use <- suppressWarnings(as.character(setdiff(as.numeric(gsub("visitfreq.raw.", "", grep("visitfreq.raw.", colnames(object@diff.data), value=T))), NA)))
    if (verbose) message(paste0("Tips not provided, so using: ", paste0(tips.use, collapse=", ")))
  }
  # Make sure tips is a character vector
  tips <- as.character(tips.use)
  # Record tips for posterity
  object@tree$tips <- tips
  # If you're going to do balanced fusion, calculate number of cells per tip.
  if (weighted.fusion) {
    tip.size <- unlist(lapply(object@tree$cells.in.tip, length))
    if (length(tip.size)==0) stop("Either weighted.fusion must be false, or loadTipCells must be run prior.")
    tip.size <- tip.size[tip.size > 0]
    if !all(tips.use %in% names(tip.size)) {
      stop(paste0("tips.use included tips (", paste0(setdiff(tips.use, names(tip.size)), collapse=", "), ") for which no cells were defined by loadTipCells. Does that tip exist in your tip clustering?"))
    }
  }
  # Load the pseudotime used into the tree object for use later with plotting
  object@tree$pseudotime <- object@pseudotime[,pseudotime]
  names(object@tree$pseudotime) <- rownames(object@pseudotime)
  # Initialize segment to add
  seg.add <- as.character(max(suppressWarnings(as.numeric(tips))) + 1)
  if (is.na(seg.add)) seg.add <- "1" # If none of the tip names were numeric.
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
  # If storing all pseudotime breakpoint info, keep track of it.
  if (save.all.breakpoint.info) {
    all.pseudotime.breakpoint.details <- object@tree$pseudotime.breakpoint.details
    ptbreak.stored <- names(all.pseudotime.breakpoint.details)
  }
  # Loop through until everything has joined into a single root.
  while(length(tips) >= 2) {
    
    # Find the oldest pseudotime breakpoint, and thereby which branches to join at which pseudotime
    # If every breakpoint has received NA, set them to 0 to close up the top of the tree.
    if (all(is.na(object@tree$segment.divergence$pseudotime.breakpoint))) {
      object@tree$segment.divergence$pseudotime.breakpoint <- 0
    }
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
    children.of.branchpoint <- unique(c(segChildrenAll(object, seg.1, include.self=T, original.joins = F, format="binary"), segChildrenAll(object, seg.2, include.self=T, original.joins = F, format="binary")))
    
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
      object <- allSegmentDivergenceByPseudotime(object, pseudotime=pseudotime, divergence.method=divergence.method, segments=tips, pseudotime.cuts=cells.per.pseudotime.bin, window.size=bins.per.pseudotime.window, minimum.visits=minimum.visits, visit.threshold=visit.threshold, p.thresh=p.thresh, breakpoint.decision.plots=save.breakpoint.plots, cache=T, verbose=verbose)
      # If storing all pseudotime breakpoint info, add the new breakpoint info.
      if (save.all.breakpoint.info) {
        new.breakpoints <- setdiff(names(object@tree$pseudotime.breakpoint.details), ptbreak.stored)
        all.pseudotime.breakpoint.details <- unlist(list(all.pseudotime.breakpoint.details, object@tree$pseudotime.breakpoint.details[new.breakpoints]), recursive = F)
        ptbreak.stored <- c(ptbreak.stored, new.breakpoints)
      }
    }
    
    # Prepare for next join
    seg.add <- as.character(as.numeric(seg.add) + 1)
  }
  
  # If storing all breakpoints place them in the tree.
  if (save.all.breakpoint.info) object@tree$pseudotime.breakpoint.details <- all.pseudotime.breakpoint.details else object@tree$pseudotime.breakpoint.details <- NULL
  
  # Determine complete list of all segments
  object@tree$segments <- as.character(sort(as.numeric(unique(unlist(object@tree$segment.joins[,c("child.1", "child.2", "parent")])))))
  
  # Make a backup copy of segment.joins in case segments are later deleted and you need this history
  object@tree$segment.joins.initial <- object@tree$segment.joins
  
  # Assign cells to segments (for tree structure)
  if (verbose) print("Assigning cells to segments.")
  object <- assignCellsToSegments(object, pseudotime, verbose)
  
  # Delete overly short segments from the tree, and reassign segments
  object <- reformatSegmentJoins(object)
  object <- reformatSegmentJoins(object, segment.joins.initial = T)
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
