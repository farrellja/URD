#' Plot Pseudotime Stability (Individual Cells)
#' 
#' @param object An URD object
#' @param cells.plot (Character vector or Numeric) Either names of cells to plot, or a number of cells to plot (chosen randomly).
#' @param log.visits (Logical) Whether to log10 transform visitation frequency
#' @return A ggplot2 object
#' @export
pseudotimePlotStabilityCells <- function(object, cells.plot=9, log.visits=T) {
  # Get cells out
  if (is.numeric(cells.plot)) cells.plot <- sample(rownames(object@pseudotime.stability$pseudotime), cells.plot)
  # Get pseudotime values for ggplot
  pseudotime.stability.melt <- melt(as.matrix(object@pseudotime.stability$pseudotime[cells.plot,]),stringsAsFactors=F)
  pseudotime.stability.melt <- pseudotime.stability.melt[complete.cases(pseudotime.stability.melt),]
  names(pseudotime.stability.melt) <- c("Cell", "Walks", "Pseudotime")
  # Get visit values for ggplot
  pseudotime.stability.melt$Visits <- object@pseudotime.stability$walks.per.cell[cbind(as.character(pseudotime.stability.melt$Cell), as.character(pseudotime.stability.melt$Walks))]
  # Do the plot
  gg <- ggplot(pseudotime.stability.melt, aes(x=Walks, y=Pseudotime)) + facet_wrap(~Cell) + ylim(c(0,1)) + geom_line() + geom_point(aes(color=Visits)) + labs(title="Pseudotime Stability Per Cell", x="Simulations", y="Pseudotime")
  if (log.visits) {
    gg <- gg + scale_color_gradientn(trans='log2', colors=defaultURDContinuousColors())  
  } else {
    gg <- gg + scale_color_gradientn(colors=defaultURDContinuousColors())
  }
  return(gg)
}

#' Plot Pseudotime Stability
#' 
#' This plots the overall change in pseudotimes (across all cells), as the number
#' of simulations is increased. This plot should become asymptotic if enough
#' simulations have been run. The number of data points shown is determined by
#' the \code{stability.div} parameter when processing floods or walks.
#' 
#' This relies on the data stored in \code{@@pseudotime.stability}, which is
#' written by both \code{\link{floodPseudotimeProcess}} and
#' \code{\link{processRandomWalksFromTips}}. Each call to either function overwrites 
#' this slot, so it will display plots from the most recent time either function
#' has been called.
#' 
#' @param object An URD object
#' 
#' @return A ggplot2 object
#' 
#' @examples
#' pseudotimePlotStabilityOverall(object)
#' 
#' @export
pseudotimePlotStabilityOverall <- function(object) {
  # Plot overall change in cell pseudotimes as number of walks increases
  pseudotime.stability.rate <- t(diff(t(object@pseudotime.stability$pseudotime)))
  pseudotime.overall.change <- apply(pseudotime.stability.rate, 2, function(pseudotime.diffs) sum(abs(pseudotime.diffs), na.rm=T))
  gg.data <- data.frame(
    Walks=as.numeric(names(pseudotime.overall.change)),
    Change=pseudotime.overall.change,
    stringsAsFactors=F
  )
  return(ggplot(gg.data, aes(x=Walks, y=Change)) + geom_line() + labs(x="Simulations", y="Total change in pseudotime", title="Overall Pseudotime Stability"))
}

# Plot histograms of visit frequences as number of walks increases
# xlim: (vector) x-axis limits (e.g. c(0,30)) to apply.
#' Plot Cell Visitation By Random Walks or Flood Pseudotime
#' 
#' @param object An URD object
#' @param xlim (Numeric vector, length 2) Limits of histogram x-axis. (NULL allows ggplot to determine automatically.)
#' @return A ggplot2 object
#' @export
pseudotimePlotVisits <- function(object, xlim=NULL) {
  # Melt data in preparation
  walks.per.cell.melt <- melt(object@pseudotime.stability$walks.per.cell)
  names(walks.per.cell.melt) <- c("Cell", "Walks", "Visits")
  walks.per.cell.melt$Walks <- as.numeric(walks.per.cell.melt$Walks)
  gg <- ggplot(walks.per.cell.melt, aes(x=Visits, group=Walks, color=Walks, fill=Walks)) + geom_histogram(binwidth = 1) + scale_color_gradientn(colors=defaultURDContinuousColors()) + scale_fill_gradientn(colors=defaultURDContinuousColors()) + labs(y="Cells", title="Cell Visit Frequencies As Random Walks Increase")
  if (!is.null(xlim)) gg <- gg + xlim(xlim)
  return(gg)
}