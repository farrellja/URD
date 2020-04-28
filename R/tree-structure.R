#' Putative Cells in Segment
#' 
#' @param object An URD object
#' @param segments (Character vector) Segment names to potentially assign cells to 
#' @param minimum.visits (Numeric) Minimum number of times a cell was visited 
#' @param visit.threshold (Numeric) Proportion of maximum visit frequency to use as cut-off for assignment to a segment
#' 
#' @return data.frame Rows are cells, columns are segments, values are logical and encode whether a cell might be part of a given segment
#' 
#' @export
#' @keywords internal
putativeCellsInSegment <- function(object, segments, minimum.visits, visit.threshold) {
  # Get visit data
  visit.data <- object@diff.data[, paste0("visitfreq.raw.", segments)]
  # Figure out the maximum visit frequency for each cell
  max.visit <- apply(visit.data, 1, max)
  # Figure out the threshold for each cell's visitation based on its maximum visit frequency
  enough.visits <- visit.threshold * max.visit
  # Generate matrix of whether each visitation frequency passes that cell's threshold.
  cell.in.segment <- sweep(visit.data, 1, enough.visits, ">=")
  # Which cells pass the minimum.visits threshold?
  pass.min <- max.visit >= minimum.visits
  # For cells that didn't pass, set them as not belonging to any segments.
  cell.in.segment <- sweep(cell.in.segment, 1, pass.min, "&")
  # Rename columns and return
  colnames(cell.in.segment) <- segments
  return(as.data.frame(cell.in.segment))
}

#' All Segment Divergence By Pseudotime
#' 
#' @importFrom utils combn
#' 
#' @param object An URD object
#' @param pseudotime (Character) Name of pseudotime to use for determining pseudotime of breakpoints
#' @param segments (Character vector) Segments to consider fot determining visit divergence
#' @param divergence.method (Character: "ks" or "preference") Test to use to determine whether visitation has diverged for each pseudotime window.
#' @param pseudotime.cuts (Numeric) Approximate number of cells to assign to each pseudotime bin for branchpoint finding.
#' @param window.size (Numeric) Width of moving window in pseudotime used for branchpoint finding, in terms of bins.
#' @param minimum.visits (Numeric) Minimum number of random walk visits to a cell to retain it in the tree
#' @param visit.threshold (Numeric) Cells are considered potential members for segments/tips from which random walks visited them 
#' at least this fraction of their maximum visitation from a single tip
#' @param p.thresh (Numeric) p-value threshold to use in determining whether visitation is significantly different from pairs of tips
#' @param pref.thresh (Numeric) Maximum mean preference within a window to be considered 'different' if \code{divergence.method=="preference"}
#' @param breakpoint.decision.plots (Path) Path to save plots summarizing (default is NULL, which does not save plots as they are somewhat slow)
#' @param cache (Logical) Used cached values? This will check \code{object@@tree$segment.divergence} and only calculate values for new segments.
#' @param verbose (Logical) Report on progress?
#' @keywords internal
allSegmentDivergenceByPseudotime <- function(object, pseudotime, segments, divergence.method=c("ks","preference"), pseudotime.cuts=80, window.size=5, minimum.visits=10, visit.threshold=0.7, p.thresh=.01, pref.thresh=0.5, breakpoint.decision.plots=NULL, cache=T, verbose=F) {
  if (length(divergence.method) > 1) divergence.method=divergence.method[1]
  if (!(divergence.method %in% c("ks", "preference"))) stop("Divergence method must be 'ks' or 'preference'.")
  # Make sure this is properly used as names, not indices
  segments <- as.character(segments)
  # Figure out which cells are in each segment
  cells.in.segments <- putativeCellsInSegment(object, segments, minimum.visits=minimum.visits, visit.threshold=visit.threshold)
  segment.overlaps <- as.data.frame(t(combn(segments, 2, simplify=T)), stringsAsFactors=F)
  names(segment.overlaps) <- c("seg.1","seg.2")
  # If you've already calculated many of these, don't recalculate them. Just update existing DF.
  if (!is.null(object@tree$segment.divergence) & cache) {
    # Figure out which combos need to be deleted
    old.segments <- unique(unlist(object@tree$segment.divergence[,c("seg.1", "seg.2")]))
    segs.to.delete <- setdiff(old.segments, segments)
    # Trim the existing data frame to remove anything you no longer want.
    divergences.to.keep <- which(!(object@tree$segment.divergence$seg.1 %in% segs.to.delete | object@tree$segment.divergence$seg.2 %in% segs.to.delete))
    object@tree$segment.divergence <- object@tree$segment.divergence[divergences.to.keep,]
    # Trim the combinations to generate to only those that are not in the existing data.frame
    divergence.cached <- rbind(segment.overlaps, object@tree$segment.divergence[,c("seg.1","seg.2")])
    segment.overlaps.not.cached <- !duplicated(divergence.cached, fromLast=T)[1:dim(segment.overlaps)[1]]
    segment.overlaps <- segment.overlaps[segment.overlaps.not.cached,]
    # Trim pseudotime.breakpoint.details also
    object@tree$pseudotime.breakpoint.details <- object@tree$pseudotime.breakpoint.details[divergences.to.keep]
  }
  # Calculate the pseudotime divergence between every pair of tips
  pseudotime.divergences <- lapply(1:dim(segment.overlaps)[1], function(overlap) {
    trim.before <- max(object@tree$segment.pseudotime.limits[unlist(segment.overlaps[overlap,c("seg.1","seg.2")]), "start"])
    trim.after <- min(object@tree$segment.pseudotime.limits[unlist(segment.overlaps[overlap,c("seg.1","seg.2")]), "end"])
    if (verbose) print(paste0("Calculating divergence between ", segment.overlaps[overlap, "seg.1"], " and ", segment.overlaps[overlap, "seg.2"], " (Pseudotime ", round(trim.before, digits=3), " to ", round(trim.after, digits=3), ")"))
    return(visitDivergenceByPseudotime(object, pseudotime, segment.1 = segment.overlaps[overlap, "seg.1"], segment.2=segment.overlaps[overlap, "seg.2"], cells.in.segments = cells.in.segments, pseudotime.cuts = pseudotime.cuts, pseudotime.min = trim.before, pseudotime.max = trim.after, window.size = window.size, p.thresh = p.thresh, divergence.method = divergence.method, verbose=verbose))
  })
  names(pseudotime.divergences) <- apply(segment.overlaps, 1, function(x) paste0(x[1], "-", x[2]))
  # Extract pseudotime breakpoints
  segment.overlaps$pseudotime.breakpoint <- unlist(lapply(pseudotime.divergences, function(x) x$breakpoint))
  # If this was a cache job, then add in the cached data.
  if (!is.null(object@tree$segment.divergence) & cache) {
    segment.overlaps <- rbind(object@tree$segment.divergence, segment.overlaps)
    pseudotime.divergences <- unlist(list(object@tree$pseudotime.breakpoint.details, pseudotime.divergences), recursive=F)
  }
  # If desired, save a plot to document the breakpoint decisions
  if (!is.null(breakpoint.decision.plots)) {
    if (verbose) print("Saving breakpoint decision plot.")
    grobs <- lapply(1:dim(segment.overlaps)[1], function(x) {
      this.decision <- pseudotime.divergences[[x]]
      this.decision$details$y <- 0
      this.decision$details$y1 <- 1
      gg <- ggplot(data=this.decision$details, aes(xmin=min.pseudotime, xmax=max.pseudotime, ymin=y, ymax=y1, fill=different)) + geom_rect(alpha=0.25) + guides(fill=F) + labs(y=paste0(segment.overlaps[x,"seg.1"], " vs. ", segment.overlaps[x,"seg.2"])) + theme_bw() + theme(axis.line=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank()) + scale_x_continuous(expand=c(0,0), limits=c(0,1)) + scale_y_continuous(expand=c(0,0), limits=c(0,1)) + scale_fill_manual(values=c("TRUE"="cyan3", "FALSE"="red"))
      if (!is.na(this.decision$breakpoint)) gg <- gg + geom_vline(aes(xintercept=max(this.decision$breakpoint, 0.0025)), color='orange', lty=2, lwd=0.8)
      return(gg)
    })
    pdf(file=paste0(breakpoint.decision.plots, "-", dim(segment.overlaps)[1], ".pdf"), width=3, height=(0.75*dim(segment.overlaps)[1]))
    grid.arrange(grobs=grobs, ncol=1)
    dev.off()
  }
  # Record the results in the object and return it
  object@tree$segment.divergence <- segment.overlaps
  object@tree$pseudotime.breakpoint.details <- pseudotime.divergences
  return(object)
}

#' Visitation Divergence for Pseudotime Windows
#' 
#' @importFrom stats embed p.adjust
#' 
#' @param object An URD object
#' @param pseudotime (Character) Name of pseudotime to use for determining pseudotime of breakpoints
#' @param segment.1 (Character) First segment to compare (Not necessary if \code{cells.in.segment.1} is set)
#' @param segment.2 (Character) Second segment to compare (Not necessary if \code{cells.in.segment.2} is set)
#' @param cells.in.segments (data.frame) Output from \code{\link{putativeCellsInSegments}} that describes which cells can belong to which segments. (Not necessary if \code{cells.in.segment.1} and \code{cells.in.segment.2} are both set)
#' @param cells.segment.1 (Character Vector) List of cells that belong to segment 1 for comparison (not necessary if \code{segment.1} and \code{cells.in.segments} are provided, in which case set to NULL)
#' @param cells.segment.2 (Character Vector) List of cells that belong to segment 2 for comparison (not necessary if \code{segment.2} and \code{cells.in.segments} are provided, in which case set to NULL)
#' @param divergence.method (Character: "ks" or "preference") Test to use to determine whether visitation has diverged for each pseudotime window.
#' @param pseudotime.cuts (Numeric) Approximate number of cells to assign to each pseudotime bin for branchpoint finding.
#' @param window.size (Numeric) Width of moving window in pseudotime used for branchpoint finding, in terms of bins.
#' @param pseudotime.min (Numeric) Minimum pseudotime of cells to compare
#' @param pseudotime.max (Numeric) Maximum pseudotime of cells to compare
#' @param p.thresh (Numeric) p-value threshold to use in determining whether visitation is significantly different from pairs of tips
#' @param pref.thresh (Numeric) Maximum mean preference within a window to be considered 'different' if \code{divergence.method=="preference"}
#' 
#' @return (List) Number of cells considered ("cells.considered"), the determined pseudotime breakpoint ("breakpoint"), and the
#' details of the calculation for each window ("details"), which is a data.frame: Rows are pseudotime windows, columns are KS-test p-value adjusted for multiple hypotheses ("p"), pseudotime of cells in the window ("mean.pseudotime", "min.pseudotime", "max.pseudotime"), number of cells considered from each segment ("cells.visited.seg1", "cells.visited.seg2"), and whether the window passed the p-value threshold for significance ("different")
#' 
#' @keywords internal
visitDivergenceByPseudotime <- function(object, pseudotime, segment.1, segment.2, cells.in.segments=NULL, cells.segment.1=NULL, cells.segment.2=NULL, divergence.method=c("ks","preference"), pseudotime.cuts=80, window.size=5, pseudotime.min=NULL, pseudotime.max=NULL, p.thresh=.01, pref.thresh=0.5, verbose=T) {
  # If method left as default, cull down.
  if (length(divergence.method) > 1) divergence.method <- divergence.method[1]
  # Make sure that cells.to.operate on are provided.
  if (is.null(cells.in.segments) && (is.null(cells.segment.1) || is.null(cells.segment.2))) {
    stop("Please provide either the cells.in.segments data.frame or cells.segment.1 and cells.segment.2 as vectors of cell names.")
  }
  # Extract cell lists from data.frame if necessary
  if (is.null(cells.segment.1)) cells.segment.1 <- rownames(cells.in.segments)[which(cells.in.segments[,segment.1])]
  if (is.null(cells.segment.2)) cells.segment.2 <- rownames(cells.in.segments)[which(cells.in.segments[,segment.2])]
  cells.to.use <- unique(c(cells.segment.1, cells.segment.2))
  # Get visit data and pseudotime for those two segments
  visit.data <- object@diff.data[cells.to.use, paste0("visitfreq.raw.", c(segment.1, segment.2))]
  names(visit.data) <- c("segment.1", "segment.2")
  visit.data$pseudotime <- object@pseudotime[rownames(visit.data), pseudotime]
  # Eliminate cells by pseudotime if trim parameters are provided
  if (!is.null(pseudotime.min)) visit.data <- visit.data[which(visit.data$pseudotime >= pseudotime.min),]
  if (!is.null(pseudotime.max)) visit.data <- visit.data[which(visit.data$pseudotime <= pseudotime.max),]
  # Cut pseudotime into groups
  visit.data <- visit.data[order(visit.data$pseudotime),]
  n.pseudotime.groups <- max(round(dim(visit.data)[1] / pseudotime.cuts),1)
  visit.data$pseudotime.group <- rep(1:n.pseudotime.groups, diff(round((dim(visit.data)[1]/n.pseudotime.groups)*0:n.pseudotime.groups)))
  # Turn groups into windows, based on window.size
  pseudotime.windows <- embed(1:n.pseudotime.groups, dimension=min(window.size,length(unique(visit.data$pseudotime.group))))
  # Convert pseudotime windows to cell lists
  cells.in.windows <- lapply(1:nrow(pseudotime.windows), function(window) {
    groups <- as.vector(as.numeric(pseudotime.windows[window,]))
    cells <- rownames(visit.data)[which(visit.data$pseudotime.group %in% groups)]
    return(cells)
  })
  
  # Calculate divergence
  if (divergence.method=="ks") {
    # Calculate divergence by KS test in each window
    div.pseudotime <- divergenceKSVisitation(visit.data=visit.data, pseudotime.windows=pseudotime.windows, cells.segment.1=cells.segment.1, cells.segment.2=cells.segment.2)
    # Multiple hypothesis correction because ran several tests
    div.pseudotime$p <- p.adjust(div.pseudotime$p, method="holm")
    # Determine whether the visit distributions are 'different' in each pseudotime window by p-value
    div.pseudotime$different <- div.pseudotime$p <= p.thresh
  } else if (divergence.method=="preference") {
    # Calculate divergence by dip test on preference in each window
    div.pseudotime <- divergencePreferenceDip(visit.data=visit.data, cells.in.windows=cells.in.windows, cells.segment.1=cells.segment.1, cells.segment.2=cells.segment.2)
    # Multiple hypothesis correction because ran several tests
    div.pseudotime$p <- p.adjust(div.pseudotime$p, method="holm")
    # Determine whether the the visit distributions are 'different': must have significant p-value for bimodality OR be unimodal, but with mean preference far from 0. (Sometimes can get a unimodal distribution, but it's not near 0, so should not fuse there)
    div.pseudotime$different <- div.pseudotime$p <= p.thresh | div.pseudotime$mean.preference > pref.thresh
  } else {
    stop("Preference must be either 'ks' or 'preference'.")
  }
  
  # Determine pseudotime breakpoint
  pt.break <- pseudotimeBreakpointByStretch(div.pseudotime, segment.1, segment.2, visit.data, pseudotime.windows, verbose=verbose)
  
  # Return it
  to.return <- list(cells.considered=length(cells.to.use), breakpoint=pt.break, cells.in.windows=cells.in.windows, details=div.pseudotime)
  return(to.return)
}

#' Calculate visitation divergence using KS-test
#' 
#' @importFrom stats ks.test
#' 
#' @param visit.data (data.frame) Rows are cells, columns are: "segment.1" and "segment.2" (visitation frequency from random walks
#' initiated in segment 1 and 2), "pseudotime" (the pseudotime of each cell), and "pseudotime.group" (the number of the pseudotime
#' bin each cell belongs to)
#' @param pseudotime.windows (matrix) Rows are pseudotime windows, and the columns in each row represent the pseudotime bins that belong to that window.
#' @param cells.segment.1 (Character vector) Cells in segment 1
#' @param cells.segment.2 (Character vector) Cells in segment 2
#'
#' @return (data.frame) Rows are pseudotime windows, columns are KS-test p-value ("p"), pseudotime of cells in the window ("mean.pseudotime", "min.pseudotime", "max.pseudotime"), and number of cells considered from each segment ("cells.visited.seg1", "cells.visited.seg2")
#' @keywords internal
divergenceKSVisitation <- function(visit.data, pseudotime.windows, cells.segment.1, cells.segment.2) {
  # Kolmogorov-Smirnov Test to determine likelihood that visitation distribution is different in each pseudotime group
  div.pseudotime <- as.data.frame(t(as.data.frame(lapply(1:dim(pseudotime.windows)[1], function(pseudotime.window) {
    # Get cells for this pseudotime moving window
    cells.in.pt.group <- rownames(visit.data)[which(visit.data$pseudotime.group %in% pseudotime.windows[pseudotime.window,])]
    # Do K-S test
    ks.p <- suppressWarnings(ks.test(x=visit.data[cells.in.pt.group, "segment.1"], y=visit.data[cells.in.pt.group, "segment.2"], exact=F)$p.value)
    # Calculate mean pseudotime of this window
    mean.pt <- mean(visit.data[cells.in.pt.group, "pseudotime"])
    min.pt <- min(visit.data[cells.in.pt.group, "pseudotime"])
    max.pt <- max(visit.data[cells.in.pt.group, "pseudotime"])
    # Cells in each segment
    cells.seg1.pt.group <- length(intersect(cells.segment.1, cells.in.pt.group))
    cells.seg2.pt.group <- length(intersect(cells.segment.2, cells.in.pt.group))
    return(c(ks.p, mean.pt, min.pt, max.pt, cells.seg1.pt.group, cells.seg2.pt.group))
  }))))
  colnames(div.pseudotime) <- c("p", "mean.pseudotime", "min.pseudotime", "max.pseudotime", "cells.visited.seg1", "cells.visited.seg2")
  rownames(div.pseudotime) <- NULL
  return(div.pseudotime)
}

#' Calculate visitation divergence based on bimodality of visitation preference
#' 
#' Calculates the preference of each cell for the two tips, tests for unimodality using Hartigan's dip test, and
#' if unimodal, tests the mean to make sure that it's 
#' 
#' @importFrom diptest dip.test
#' 
#' @param visit.data (data.frame) Rows are cells, columns are: "segment.1" and "segment.2" (visitation frequency from random walks
#' initiated in segment 1 and 2), "pseudotime" (the pseudotime of each cell), and "pseudotime.group" (the number of the pseudotime
#' bin each cell belongs to)
#' @param cells.in.windows (List) List of character vectors, giving names of cells in each pseudotime window.
#' @param cells.segment.1 (Character vector) Cells in segment 1
#' @param cells.segment.2 (Character vector) Cells in segment 2
#' 
#' @return (data.frame) Rows are pseudotime windows, columns are a numeric representation of whether (0 or 1) ("p"), pseudotime of cells in the window ("mean.pseudotime", "min.pseudotime", "max.pseudotime"), and number of cells considered from each segment ("cells.visited.seg1", "cells.visited.seg2")
#' @keywords internal
divergencePreferenceDip <- function(visit.data, cells.in.windows, cells.segment.1, cells.segment.2) {
  
  # Calculate visitation preference for each cell in visit.data
  visit.data$preference <- apply(visit.data, 1, function(x) preference(x[1], x[2], signed=T))
  
  # Cycle through pseudotime moving windows to build div.pseudotime data frame
  div.pseudotime <- as.data.frame(t(as.data.frame(lapply(cells.in.windows, function(cells.in.pt.group) {
    # Calculate mean pseudotime of this window
    mean.pt <- mean(visit.data[cells.in.pt.group, "pseudotime"])
    min.pt <- min(visit.data[cells.in.pt.group, "pseudotime"])
    max.pt <- max(visit.data[cells.in.pt.group, "pseudotime"])
    # Cells in each segment
    cells.seg1.pt.group <- length(intersect(cells.segment.1, cells.in.pt.group))
    cells.seg2.pt.group <- length(intersect(cells.segment.2, cells.in.pt.group))
    # Test for unimodality with Hartigan's diptest
    dip <- diptest::dip.test(visit.data[cells.in.pt.group,"preference"])
    p <- dip$p.value
    # Determine mean of preference
    mean.pref <- mean(abs(visit.data[cells.in.pt.group,"preference"]), na.rm=T)
    return(c(mean.pt, min.pt, max.pt, cells.seg1.pt.group, cells.seg2.pt.group, p, mean.pref))
  }))))
  colnames(div.pseudotime) <- c("mean.pseudotime", "min.pseudotime", "max.pseudotime", "cells.visited.seg1", "cells.visited.seg2", "p.value", "mean.preference")
  rownames(div.pseudotime) <- NULL
  
  return(div.pseudotime)
}

#' Find Pseudotime Breakpoint
#' 
#' 
pseudotimeBreakpointByStretchV1 <- function(div.pseudotime, segment.1, segment.2, visit.data, pseudotime.windows, verbose=T) {
  # Find the longest stretches of all-different and all-not-significantly-different pseudotime windows, and find the area in between
  pt.rle <- rle(div.pseudotime$different)
  pt.rle <- data.frame(lengths=pt.rle$lengths,
                       values=pt.rle$values)
  pt.rle$end <- cumsum(pt.rle$lengths)
  pt.rle$start <- head(c(0,pt.rle$end) + 1, -1)
  # If the comparison is not-significantly different the entire way across, then it means that you should set the breakpoint of these to
  # the very end of the pseudotime range they occupy
  if (dim(pt.rle)[1] == 1 && pt.rle[1,"values"] == FALSE) {
    pt.break <- tail(div.pseudotime, 1)$max.pseudotime
  }
  # These are constantly different, which suggests that there may be a disconnection in the graph.
  # Set their breakpoint to be at the very beginning. (You'll eventually need this to hook in the
  # PGCs.)
  else if (dim(pt.rle[1]) == 1 && pt.rle[1,"values"] == TRUE) {
    if (verbose) warning(paste("No obvious breakpoint between", segment.1, "and", segment.2, ". Difference constantly", pt.rle[1, "values"], "\n"))
    pt.break <- min(div.pseudotime$min.pseudotime)
  }
  # If you have perfect stretches where it switches, just take the time in between
  else if (dim(pt.rle)[1] == 2) {
    windows.to.use.for.pt.break <- table(unlist(pseudotime.windows[c(pt.rle[1,"end"], pt.rle[2,"start"]),]))
    groups.to.use.for.pt.break <- names(windows.to.use.for.pt.break)[windows.to.use.for.pt.break == 2]
    pt.break <- mean(visit.data[visit.data$pseudotime.group %in% groups.to.use.for.pt.break, "pseudotime"])
  }
  # Divergence fluctuates back and forth between being significant, so find longest stretch
  #   of significantly different & not-significantly-different
  else {
    pt.rle <- pt.rle[order(pt.rle$lengths, decreasing = T),]
    true.starts <- pt.rle[min(which(pt.rle$values)), "start"]
    false.ends <- pt.rle[min(which(!pt.rle$values)), "end"]
    # If true and false appear in wrong order, then just warn that there's no sensical branchpoint and don't calculate one.
    if (false.ends > true.starts) {
      if (verbose) warning(paste("No obvious breakpoint between", segment.1, "and", segment.2, ". Longest stretch of difference is upstream of longest stretch of non-different.\n"))
      groups.to.use.for.pt.break <- NULL
      pt.break <- NA
    } else {
      # Use only those pseudotime breaks that appear in all windows in the uncertain region
      uncertain.breaks <- table(unlist(pseudotime.windows[(false.ends+1):(true.starts-1),]))
      groups.to.use.for.pt.break <- names(uncertain.breaks)[uncertain.breaks == max(uncertain.breaks)]
      pt.break <- mean(visit.data[visit.data$pseudotime.group %in% groups.to.use.for.pt.break, "pseudotime"])
    }
  }
  return(pt.break)
}

#' Pseudotime Breakpoint By Stretch
#' 
#' @keywords internal
pseudotimeBreakpointByStretch <- function(div.pseudotime, segment.1, segment.2, visit.data, pseudotime.windows, verbose=T) {
  # Find the longest stretches of all-different and all-not-significantly-different pseudotime windows, and find the area in between
  pt.rle <- rle(div.pseudotime$different)
  pt.rle <- data.frame(lengths=pt.rle$lengths,
                       values=pt.rle$values)
  pt.rle$end <- cumsum(pt.rle$lengths)
  pt.rle$start <- head(c(0,pt.rle$end) + 1, -1)
  # If the two populations are never significantly different, set the breakpoint to the end -- these populations will fuse immediately.
  if (nrow(pt.rle) == 1 && pt.rle[1,"values"] == FALSE) {
    if (verbose) message(paste("Difference between", segment.1, "and", segment.2, "always FALSE -- setting breakpoint to end."))
    pt.break <- tail(div.pseudotime, 1)$max.pseudotime
  }
  # If the two populations are always significantly different, set the breakpoint to the very beginning -- these populations were different at the root, so they should "connect" at the very root.
  else if (nrow(pt.rle) == 1 && pt.rle[1,"values"] == TRUE) {
    if (verbose) message(paste("Difference between", segment.1, "and", segment.2, "always TRUE -- setting breakpoint to beginning."))
    pt.break <- min(div.pseudotime$min.pseudotime)
  }
  # If you have two stretches: not-different (FALSE) then different (TRUE), assign pseudotime as the boundary between the two windows.
  else if (nrow(pt.rle) == 2 && pt.rle[1,"values"] == FALSE) {
    windows.to.use.for.pt.break <- table(unlist(pseudotime.windows[c(pt.rle[1,"end"], pt.rle[2,"start"]),]))
    groups.to.use.for.pt.break <- names(windows.to.use.for.pt.break)[windows.to.use.for.pt.break == 2]
    pt.break <- mean(visit.data[visit.data$pseudotime.group %in% groups.to.use.for.pt.break, "pseudotime"])
  } 
  # If you have two stretches: different (TRUE) then not-different (FALSE), which doesn't make any biological sense, assign as NA.
  else if (nrow(pt.rle) == 2 && pt.rle[1,"values"] == TRUE) {
    if (verbose) message(paste("No obvious breakpoint between", segment.1, "and", segment.2, ". Longest stretch of difference is upstream of longest stretch of non-different."))
    groups.to.use.for.pt.break <- NULL
    pt.break <- NA
  }
    
  # Divergence fluctuates back and forth between being significant, so find longest stretch
  #   of significantly different & not-significantly-different
  else {
    pt.rle <- pt.rle[order(pt.rle$lengths, decreasing = T),]
    true.starts <- pt.rle[min(which(pt.rle$values)), "start"]
    false.ends <- pt.rle[min(which(!pt.rle$values)), "end"]
    # If true and false appear in wrong order, then just warn that there's no sensical branchpoint and don't calculate one.
    if (false.ends > true.starts) {
      if (verbose) warning(paste("No obvious breakpoint between", segment.1, "and", segment.2, ". Longest stretch of difference is upstream of longest stretch of non-different.\n"))
      groups.to.use.for.pt.break <- NULL
      pt.break <- NA
    } else {
      # Use only those pseudotime breaks that appear in all windows in the uncertain region
      uncertain.breaks <- table(unlist(pseudotime.windows[(false.ends+1):(true.starts-1),]))
      groups.to.use.for.pt.break <- names(uncertain.breaks)[uncertain.breaks == max(uncertain.breaks)]
      pt.break <- mean(visit.data[visit.data$pseudotime.group %in% groups.to.use.for.pt.break, "pseudotime"])
    }
  }
  return(pt.break)
}

#' Assign Cells to Segments
#' 
#' Assigns cells get assigned to segments based on whichever segment visited them
#' the most often, out of those segments that exist at that cell's pseudotime. This
#' function is called automatically by \code{\link{buildTree}}, but can be re-run
#' if necessary.
#' 
#' @param object An URD object
#' @param pseudotime (Character) Pseudotime to use (i.e. a column name of \code{@@pseudotime})
#' 
#' @return An URD object with segment assignments in \code{object@tree$cells.in.segment}, as a
#' list by segments, and in \code{object@diff.data$segment}, which allows them \code{"segment"}
#' to be used as a plotting label.
#' 
#' @keywords internal
assignCellsToSegments <- function(object, pseudotime, verbose=T) {
  # Get rid of the old 'cells.in.segments' data.frame now that you're done
  if (!is.null(object@tree$cells.in.segments)) object@tree$cells.in.segments <- NULL
  # Get visit data
  visit.data <- object@diff.data[,paste0("visitfreq.raw.", object@tree$segments)]
  segments <- object@tree$segments
  # For cells in tip.clusters, make sure they have high visitation from their own tip
  # to compensate for the fact that tip cells have unusual visitation parameters,
  # because they are randomly chosen for starting walks with a uniform distribution.
  # Set to max visitation from any segment + 1, so that cells in tips will be assigned
  # to the segment which they define regardless of their visitation.
  for (tip in object@tree$tips) {
    tip.to.assign <- tip
    # Get tip cells
    cells.in.tip <- object@tree$cells.in.tip[[tip]]
    # Check whether segment exists in the tree
    while(!(tip.to.assign %in% rownames(object@tree$segment.pseudotime.limits))) {
      tip.to.assign <- segParent(object, tip.to.assign, original.joins=T)
    }
    # Set visitation frequency to max + 1
    visit.data[cells.in.tip, paste0("visitfreq.raw.", tip.to.assign)] <- apply(visit.data[cells.in.tip,], 1, max) + 1
  }
  # Zero out visitation outside of the segment limits for each segment
  for (segment in segments) {
    visit.data[which(object@pseudotime[,pseudotime] > object@tree$segment.pseudotime.limits[segment, "end"] | object@pseudotime[,pseudotime] < object@tree$segment.pseudotime.limits[segment, "start"]),paste0("visitfreq.raw.", segment)] <- 0
  }
  # Figure out which cells are now not visited enough
  cells.still.visited <- rownames(visit.data)[which(apply(visit.data, 1, max) > 0)]
  cells.removed <- rownames(visit.data)[which(apply(visit.data, 1, max) == 0)]
  if ((length(cells.removed) > 0) & verbose) warning(paste(length(cells.removed), "cells were not visited by a branch that exists at their pseudotime and were not assigned."))
  # Cells get assigned to the branch for which they have the maximum visit.frequency
  vd.segs <- gsub("visitfreq.raw.", "", names(visit.data))
  cell.assignments <- apply(visit.data, 1, function(x) vd.segs[which.max(x)])
  cell.assignments[cells.removed] <- NA
  # Record the information
  object@diff.data$segment <- cell.assignments
  object@group.ids[rownames(object@diff.data),"segment"] <- object@diff.data$segment
  object@tree$cells.in.segment <- lapply(segments, function(segment) names(which(cell.assignments == segment)))
  names(object@tree$cells.in.segment) <- segments
  
  return(object)
}

# Restructure segment.joins data.frame so that it could potentially deal with something non-binary before collapsing short segments
#' Reformat Segment Joins
#' 
#' @importFrom gdata interleave
#' 
#' @keywords internal
reformatSegmentJoins <- function(object, segment.joins.initial=F) {
  # Get proper segment joins 
  if (segment.joins.initial) sj <- object@tree$segment.joins.initial else sj <- object@tree$segment.joins
  # Reformat it
  sj.c1 <- sj[,c("parent","child.1","pseudotime")]
  sj.c2 <- sj[,c("parent","child.2","pseudotime")]
  names(sj.c1) <- c("parent", "child", "pseudotime"); names(sj.c2) <- c("parent", "child", "pseudotime")
  sj <- gdata::interleave(sj.c1, sj.c2)
  rownames(sj) <- NULL
  # Place it back in the right place in the object
  if (segment.joins.initial) object@tree$segment.joins.initial <- sj else object@tree$segment.joins <- sj
  # Return the object
  return(object)
}

#' Collapse Short Segments
#' 
#' This refines the dendrogram structure found by URD by removing overly short segments.
#' It is automatically called by \code{\link{buildTree}}, but can be re-run if needed.
#' 
#' @param object An URD object
#' @param min.cells.per.segment (Numeric) Remove segments with fewer cells assigned
#' @param min.pseudotime.per.segment (Numeric) Remove segments shorter than this in pseudotime
#' @param collapse.root (Logical) Is the very root most segment OK to delete? (Sometimes it would have no segments, but be important not to collapse)
#' 
#' @keywords internal
collapseShortSegments <- function(object, min.cells.per.segment=1, min.pseudotime.per.segment=0.01, collapse.root=F) {
  # Figure out how many cells and how long (in pseudotime) each segment is
  n.cells.segment <- unlist(lapply(object@tree$cells.in.segment, length))
  pt.length.segment <- abs(apply(object@tree$segment.pseudotime.limits, 1, function(x)diff(as.numeric(x))))
  # Which segments fail to meet the conditions
  segments.gotta.go <- unique(c(names(which(n.cells.segment < min.cells.per.segment)), names(which(pt.length.segment < min.pseudotime.per.segment))))
  # Collapse root = F: don't collapse the very top segment no matter what (prevents PGC-somatic breaking)
  root.segment <- tail(object@tree$segments, 1)
  if (!collapse.root) segments.gotta.go <- setdiff(segments.gotta.go, root.segment)
  # Loop through those segments
  for (segment in segments.gotta.go) {
    ## object@tree$segment.joins
    # Which segment joins did the segment you're deleting belong to?
    seg.is.child.id <- which(object@tree$segment.joins$child == segment)
    seg.is.parent.id <- which(object@tree$segment.joins$parent == segment)
    # Replace any connections where this segment was the parent
    new.parent <- object@tree$segment.joins[seg.is.child.id, "parent"]
    object@tree$segment.joins[seg.is.parent.id, "parent"] <- new.parent
    # Sync the pseudotime breakpoint between these so all segments with the new parent have the same breakpoint
    adjusted.pt.breakpoint <- max(object@tree$segment.joins[c(seg.is.parent.id,seg.is.child.id), "pseudotime"])
    object@tree$segment.joins[which(object@tree$segment.joins$parent == new.parent),"pseudotime"] <- adjusted.pt.breakpoint
    # Delete the rows where the deleted segment was a child
    object@tree$segment.joins <- object@tree$segment.joins[-seg.is.child.id,]
    ## object@tree$segments
    # Delete from the list of segments
    object@tree$segments <- setdiff(object@tree$segments, segment)
    ## object@tree$pseudotime.limits
    # Delete the row for the segment you're deleting
    object@tree$segment.pseudotime.limits <- object@tree$segment.pseudotime.limits[object@tree$segments,]
    # Adjust end pseudotime for the new parent
    object@tree$segment.pseudotime.limits[new.parent,"end"] <- adjusted.pt.breakpoint
    # Adjust start pseudotime for all children of the new parent
    new.children <- object@tree$segment.joins[which(object@tree$segment.joins$parent == new.parent),"child"]
    object@tree$segment.pseudotime.limits[new.children,"start"] <- adjusted.pt.breakpoint
    
  }
  # Return the cleaned up object
  return(object)
}

# Function to combine segments when a parent has only one child
# (i.e. "middleman" segments)
#' Remove Unitary Segments
#' @keywords internal
removeUnitarySegments <- function(object) {
  # Which segments have only one child
  unitary.segments <- names(which(table(object@tree$segment.joins$parent) == 1))
  # What are those children?
  unitary.children <- object@tree$segment.joins[which(object@tree$segment.joins$parent %in% unitary.segments), "child"]
  # Remove those singleton children & make them part of their parent
  for (segment in unitary.children) {
    ## object@tree$segment.joins
    # Which segment joins did the segment you're deleting belong to?
    seg.is.child.id <- which(object@tree$segment.joins$child == segment)
    seg.is.parent.id <- which(object@tree$segment.joins$parent == segment)
    # Replace any connections where this segment was the parent
    new.parent <- object@tree$segment.joins[seg.is.child.id, "parent"]
    object@tree$segment.joins[seg.is.parent.id, "parent"] <- new.parent
    # Delete the rows where the deleted segment was a child
    object@tree$segment.joins <- object@tree$segment.joins[-seg.is.child.id,]
    ## object@tree$segments
    # Delete from the list of segments
    object@tree$segments <- setdiff(object@tree$segments, segment)
    ## object@tree$pseudotime.limits
    # Adjust end pseudotime for the new parent
    object@tree$segment.pseudotime.limits[new.parent,"end"] <- object@tree$segment.pseudotime.limits[segment,"end"]
    # Delete the row for the segment you're deleting
    object@tree$segment.pseudotime.limits <- object@tree$segment.pseudotime.limits[object@tree$segments,]
    
  }
  # Return the cleaned up object
  return(object)
}



# Segments get divided up into nodes of approximately node.size
# Node membership is stored in object@tree$cells.in.node
# Node membership is also stored in the "node" column object@diff.data
# Mean pseudotime of each node is calculated and stored in object@tree$node.pseudotime
# Also builds an edge list for building the igraph plots, which is stored in object@tree$edge.list
#' Assign Cells to Nodes
#' 
#' @param object An URD object
#' @param node.size (Numeric) Number of cells per node
#' @param pseudotime (Character) Pseudotime to use (i.e. a column name of \code{@@pseudotime})
#' 
#' @keywords internal
assignCellsToNodes <- function(object, node.size=100, pseudotime) {
  # Initialize edgelist and node.pseudotime
  edgelist <- NULL
  # Root segment
  root.segment <- tail(object@tree$segments, 1)
  # Process all segments
  cells.in.nodes <- unlist(lapply(object@tree$segments, function(segment) {
    # Get cells from this segment and their pseudotime, then order by pseudotime
    seg.cells <- data.frame(cell=object@tree$cells.in.segment[[segment]], stringsAsFactors=F)
    seg.cells$pseudotime <- object@pseudotime[seg.cells$cell, pseudotime]
    seg.cells <- seg.cells[order(seg.cells$pseudotime),]
    # Adjust number and size of nodes to be best fit node.size and number of cells in segment
    n.nodes <- max(round(dim(seg.cells)[1] / node.size), 1)
    s.nodes <- dim(seg.cells)[1] / n.nodes
    # Add connections between all nodes in this segment to edgelist (if there are multiple nodes in this segment)
    if (n.nodes > 1) {
      edgelist <<- rbind(edgelist, t(sapply(1:(n.nodes-1), function(node) paste(segment, c(node, node+1), sep="-"))))
    }
    # Add connections between the final node of this segment and first of children segments (if there are any)
    segment.joins.thisparent <- unique(unlist(object@tree$segment.joins[which(object@tree$segment.joins$parent==segment),"child"]))
    if (length(segment.joins.thisparent) > 0) {
      edgelist <<- rbind(edgelist, t(sapply(segment.joins.thisparent, function(child) c(paste(segment, n.nodes, sep="-"), paste0(child, "-1")))))
    }
    # Assign cells to a node, and generate a list of cells in each node, return it
    seg.cells$node.membership <- rep(1:n.nodes, times=diff(round(0:n.nodes*s.nodes)))
    object@diff.data[seg.cells$cell, "node"] <<- paste(segment, seg.cells$node.membership, sep="-")
    cells.in.node.thisseg <- lapply(1:n.nodes, function(node) seg.cells[which(seg.cells$node.membership==node), "cell"])
    names(cells.in.node.thisseg) <- paste(segment, 1:n.nodes, sep="-")
    return(cells.in.node.thisseg)
  }), recursive = F)
  # Calculate node mean and max pseudotime
  node.mean.pseudotime <- unlist(lapply(cells.in.nodes, function(node.cells) {
    if (length(node.cells) > 0) mean(object@pseudotime[node.cells, pseudotime]) else NA
  }))
  node.max.pseudotime <- unlist(lapply(cells.in.nodes, function(node.cells) {
    if (length(node.cells) > 0) max(object@pseudotime[node.cells, pseudotime]) else NA
  }))
  # Modify and return object
  object@tree$cells.in.nodes <- cells.in.nodes
  object@tree$node.mean.pseudotime <- node.mean.pseudotime
  object@tree$node.max.pseudotime <- node.max.pseudotime
  object@tree$edge.list <- edgelist
  object@group.ids[rownames(object@diff.data),"node"] <- object@diff.data$node
  return(object)
}
