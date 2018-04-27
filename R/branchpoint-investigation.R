#' Investigate Branchpoint - tSNE Visitation Groups
#' 
#' For a given pair of segments, this shows the branchpoint fusion details in the form of a series of tSNE plots.
#' Each moving window through pseudotime that was considered is rendered as a plot, in order of pseudotime. (So,
#' the first plots are near the root.) Cells within the window are
#' colored according to their visitation frequency from the two segments in a dual color plot. Red is visitation
#' from one segment, green is visitation from the other segment. Cells that are visited equally will be yellow.
#' 
#' @importFrom gridExtra grid.arrange
#' 
#' @param object An URD object
#' @param seg.1 (Character) Segment 1 of fusion
#' @param seg.2 (Character) Segment 2 of fusion
#' @param file (Character) Path to save plot to. (If \code{NULL}, will display live)
#' @param file.width (Numeric) Width per plot (in pixels) if saving file.
#' @param file.height (Numeric) Height per plot (in pixels) if saving file.
#' 
#' @return Nothing. If \code{file==NULL}, will display a plot, otherwise will produce a file.
#' 
#' @export
branchpointDetailsVisitTsne <- function(object, seg.1, seg.2, file=NULL, file.width=750, file.height=600) {
  # Check whether seg.1 and seg.2 need to be switched
  if (paste0(seg.1, "-", seg.2) %in% names(object@tree$pseudotime.breakpoint.details)) {
    pbd.name <- paste0(seg.1, "-", seg.2)
  } else if (paste0(seg.2, "-", seg.2) %in% names(object@tree$pseudotime.breakpoint.details)) {
    pbd.name <- paste0(seg.2, "-", seg.1)
    seg.switch <- seg.1; seg.1 <- seg.2; seg.2 <- seg.switch
  } else {
    stop("A comparison between segments ", seg.1, " and ", seg.2, " was not found.")
  }
  pbd <- object@tree$pseudotime.breakpoint.details[[pbd.name]]
  plots <- lapply(1:length(pbd$cells.in.windows), function(pt.window) {
    cells <- pbd$cells.in.windows[[pt.window]]
    object@diff.data$visitfreq.log.seg1 <- NA
    object@diff.data$visitfreq.log.seg2 <- NA
    object@diff.data[cells, "visitfreq.log.seg1"] <- object@diff.data[cells, paste0("visitfreq.log.", seg.1)]
    object@diff.data[cells, "visitfreq.log.seg2"] <- object@diff.data[cells, paste0("visitfreq.log.", seg.2)]
    plotDimDual(object, "visitfreq.log.seg1", "visitfreq.log.seg2", na.rm=F, na.alpha=0.1, plot.title=paste(pt.window, pbd$details[pt.window,"different"]))
  })
  # Either return the plot or save it directly to a PNG
  if (is.null(file)) {
    return(gridExtra::grid.arrange(grobs=plots, ncol=ceiling(sqrt(length(pbd$cells.in.windows))), top=paste0(seg.1, " vs ", seg.2)))
  } else {
    n.cols=ceiling(sqrt(length(pbd$cells.in.windows)))
    n.rows=ceiling(length(pbd$cells.in.windows) / n.cols)
    png(file=file, width=file.width*n.cols, height=file.height*n.rows)
    gridExtra::grid.arrange(grobs=plots, ncol=ceiling(sqrt(length(pbd$cells.in.windows))), top=paste0(seg.1, " vs ", seg.2))
    dev.off()
  }
}

#' Investigate Branchpoint - Visit Frequency Distributions
#' 
#' For a given pair of segments, this shows the branchpoint fusion details in the form of a series of distributionplots.
#' Each moving window through pseudotime that was considered is rendered as a plot, in order of pseudotime. (So,
#' the first plots are near the root.) The distribution of visitation frequency from the two segments is plotted in
#' red (seg.1) and green (seg.2). Plots are titled with the window number, p-value (drawn from the test used during
#' construction of the tree), and the determination of whether cell visitation within the window is significantly
#' different (TRUE) or not (FALSE).
#' 
#' @importFrom gridExtra grid.arrange
#' 
#' @param object An URD object
#' @param seg.1 (Character) Segment 1 of fusion
#' @param seg.2 (Character) Segment 2 of fusion
#' @param file (Character) Path to save plot to. (If \code{NULL}, will display live)
#' @param file.width (Numeric) Width per plot (in pixels) if saving file.
#' @param file.height (Numeric) Height per plot (in pixels) if saving file.
#' 
#' @return Nothing. If \code{file==NULL}, will display a plot, otherwise will produce a file.
#' 
#' @export
branchpointDetailsVisitDist <- function(object, seg.1, seg.2, file=NULL, file.width=375, file.height=300) {
  # Check whether seg.1 and seg.2 need to be switched
  if (paste0(seg.1, "-", seg.2) %in% names(object@tree$pseudotime.breakpoint.details)) {
    pbd.name <- paste0(seg.1, "-", seg.2)
  } else if (paste0(seg.2, "-", seg.2) %in% names(object@tree$pseudotime.breakpoint.details)) {
    pbd.name <- paste0(seg.2, "-", seg.1)
    seg.switch <- seg.1; seg.1 <- seg.2; seg.2 <- seg.switch
  } else {
    stop("A comparison between segments ", seg.1, " and ", seg.2, " was not found.")
  }
  # Get breakpoint details out of object
  pbd <- object@tree$pseudotime.breakpoint.details[[pbd.name]]
  # Make a plot for each window
  plots <- lapply(1:length(pbd$cells.in.windows), function(pt.window) {
    cells <- pbd$cells.in.windows[[pt.window]]
    gg.data <- data.frame(
      cell=c(cells,cells),
      seg=rep(c(seg.1, seg.2), each=length(cells)),
      visit=c(object@diff.data[cells,paste0("visitfreq.log.", seg.1)],object@diff.data[cells,paste0("visitfreq.log.", seg.2)]),
      stringsAsFactors=F
    )
    color.vec <- c("red", "green"); names(color.vec) <- c(seg.1, seg.2)
    return(ggplot(gg.data, aes(x=visit, fill=seg)) + geom_density(alpha=0.5, color=NA) + theme_bw() + scale_fill_manual(values=color.vec) + guides(fill=F) + labs(x="Visit Frequency (log10)", y="", title=paste0(pt.window, ": p=", round(pbd$details[pt.window,"p"],3), " (", pbd$details[pt.window,"different"], ")")))
  })
  # Either return the plot or save it directly to a PNG
  if (is.null(file)) {
    return(gridExtra::grid.arrange(grobs=plots, ncol=ceiling(sqrt(length(pbd$cells.in.windows))), top=paste0(seg.1, "(Red) vs ", seg.2, "(Green)")))
  } else {
    n.cols=ceiling(sqrt(length(pbd$cells.in.windows)))
    n.rows=ceiling(length(pbd$cells.in.windows) / n.cols)
    png(file=file, width=file.width*n.cols, height=file.height*n.rows)
    gridExtra::grid.arrange(grobs=plots, ncol=ceiling(sqrt(length(pbd$cells.in.windows))), top=paste0(seg.1, "(Red) vs ", seg.2, "(Green)"))
    dev.off()
  }
}


#' Investigate Branchpoint - Visit Preference Distributions
#' 
#' For a given pair of segments, this shows the branchpoint fusion details in the form of a series of distributionplots.
#' Each moving window through pseudotime that was considered is rendered as a plot, in order of pseudotime. (So,
#' the first plots are near the root.) The distribution of visitation preference for cells within the windows is
#' plotted (-1 is cells visited exclusively by seg.1, 0 is cells visited equally, and 1 is cells visited exclusively
#' by seg.2). Plots are titled with the window number, p-value (drawn from the test used during
#' construction of the tree), and the determination of whether cell visitation within the window is significantly
#' different (TRUE) or not (FALSE).
#' 
#' @importFrom gridExtra grid.arrange
#' 
#' @param object An URD object
#' @param seg.1 (Character) Segment 1 of fusion
#' @param seg.2 (Character) Segment 2 of fusion
#' @param file (Character) Path to save plot to. (If \code{NULL}, will display live)
#' @param file.width (Numeric) Width per plot (in pixels) if saving file.
#' @param file.height (Numeric) Height per plot (in pixels) if saving file.
#' 
#' @return Nothing. If \code{file==NULL}, will display a plot, otherwise will produce a file.
#' 
#' @export
branchpointDetailsPreferenceDist <- function(object, seg.1, seg.2, file=NULL, file.width=375, file.height=300) {
  # Check whether seg.1 and seg.2 need to be switched
  if (paste0(seg.1, "-", seg.2) %in% names(object@tree$pseudotime.breakpoint.details)) {
    pbd.name <- paste0(seg.1, "-", seg.2)
  } else if (paste0(seg.2, "-", seg.2) %in% names(object@tree$pseudotime.breakpoint.details)) {
    pbd.name <- paste0(seg.2, "-", seg.1)
  } else {
    stop("A comparison between segments ", seg.1, " and ", seg.2, " was not found.")
  }
  # Get breakpoint details out of object
  pbd <- object@tree$pseudotime.breakpoint.details[[pbd.name]]
  # Make a plot for each window
  plots <- lapply(1:length(pbd$cells.in.windows), function(pt.window) {
    cells <- pbd$cells.in.windows[[pt.window]]
    visit.data <- object@diff.data[cells, paste0("visitfreq.raw.", c(seg.1, seg.2))]
    visit.data$preference <- apply(visit.data, 1, function(x) preference(x[1], x[2], signed=T))
    title <- paste0(pt.window, ": p=", round(pbd$details[pt.window,"p"],3), " (", pbd$details[pt.window,"different"], ")")
    return(ggplot(visit.data, aes(x=preference)) + geom_density(alpha=0.5, fill='black') + theme_bw() + labs(x="Preference", y="", title=title))
  })
  # Either return the plot or save it directly to a PNG
  if (is.null(file)) {
    return(gridExtra::grid.arrange(grobs=plots, ncol=ceiling(sqrt(length(pbd$cells.in.windows))), top=paste0(seg.1, " vs ", seg.2)))
  } else {
    n.cols=ceiling(sqrt(length(pbd$cells.in.windows)))
    n.rows=ceiling(length(pbd$cells.in.windows) / n.cols)
    png(file=file, width=file.width*n.cols, height=file.height*n.rows)
    gridExtra::grid.arrange(grobs=plots, ncol=ceiling(sqrt(length(pbd$cells.in.windows))), top=pbd.name)
    dev.off()
  }
}