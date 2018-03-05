#' Dot Plot
#' 
#' importFrom reshape2 melt
#' 
#' @param object An URD object
#' @param genes (Character Vector) Genes to plot
#' @param clustering (Character) Name of clustering to use (i.e. a column name of \code{@@group.ids})
#' @param clusters.use (Character vector) Names of specific clusters to plot (default: all clusters)
#' @param min.exp (Numeric) Minimum proportion of expressing cells (0-1) to be shown on the plot
#' @param mean.expressing.only (Logical) Should mean expression value exclude cells with no expression
#' @param title (Character) How should the plot be titled? (Default: no title)
#' 
#' @return A ggplot2 object
#' 
#' @export
plotDot <- function(object, genes, clustering, clusters.use=NULL, min.exp=.05, mean.expressing.only=F, title="") {
  # Get data out that you want to plot
  cluster.ids <- object@group.ids[colnames(object@logupx.data), clustering]
  if (is.null(clusters.use)) clusters.use <- sort(unique(cluster.ids))
  data.use <- as.data.frame(as.matrix(t(object@logupx.data[genes, which(cluster.ids %in% clusters.use)])))
  cluster.ids <- cluster.ids[cluster.ids %in% clusters.use]
  
  # Calculate Stats
  if (mean.expressing.only) mean.fun <- "mean.of.logs.pos" else mean.fun <- "mean.of.logs"
  mean.expression <- aggregate(data.use, by=list(cluster.ids), FUN=mean.fun)
  percent.expressed <- aggregate(data.use, by=list(cluster.ids), FUN="prop.exp")
  
  # Make ggplot ready
  mean.melt <- reshape2::melt(data=mean.expression, id.vars=c("Group.1"))
  names(mean.melt) <- c("Cluster", "Gene", "Mean")
  expr.melt <- reshape2::melt(data=percent.expressed, id.vars=c("Group.1"))
  names(expr.melt) <- c("Cluster", "Gene", "Prop.Exp")
  gg.data <- merge(mean.melt, expr.melt)
  
  # Remove values where proportion expressing is too small and should not be plotted.
  gg.data <- gg.data[which(gg.data$Prop.Exp >= min.exp),]
  
  # Go ahead and ggplot
  the.plot <- ggplot(data=gg.data, aes(x=Cluster, y=Gene, color=Mean, size=Prop.Exp)) + geom_point() + scale_color_gradientn(colours = defaultURDContinuousColors()) + ggtitle(title)
  if (!is.null(clusters.use)) the.plot <- the.plot + scale_x_discrete(limits=clusters.use)
  return(the.plot)
}

#' Violin plot of gene expression
#' 
#' @param object An URD object
#' @param labels.plot (Character vector)
#' @param clustering (Character) Name of clustering to use (i.e. a column name of \code{@@group.ids})
#' @param clusters.use (Character vector) Names of clusters to display (default NULL uses all clusters)
#' @param legend (Logical) Display a legend?
#' @param free.axes (Logical) Should the scale of each label be individually determined?
#' 
#' @return A ggplot2 object
#' 
#' @export
plotViolin <- function(object, labels.plot, clustering, clusters.use=NULL, legend=T, free.axes=F) {
  # Get data.to.plot
  data.plot <- as.data.frame(lapply(labels.plot, function (label) {
    data.for.plot(object = object, label=label)
  }))
  names(data.plot) <- labels.plot
  data.plot$Cluster <- object@group.ids[rownames(data.plot), clustering]
  if (!any(is.na(as.numeric(data.plot$Cluster)))) {
    cluster.names <- sort(as.numeric(unique(data.plot$Cluster)))
  } else {
    cluster.names <- sort(unique(data.plot$Cluster))
  }
  if (!is.null(clusters.use)) data.plot <- data.plot[data.plot$Cluster %in% clusters.use,]
  data.plot$Cell <- rownames(data.plot)
  data.plot.melt <- melt(data.plot, id.vars = c("Cell", "Cluster"))
  data.plot.melt$Cluster <- factor(data.plot.melt$Cluster, levels=cluster.names, ordered=T)
  if (free.axes) {free="free"} else {free="fixed"}
  the.plot <- ggplot(data.plot.melt, aes(x=Cluster, y=value, fill=Cluster)) + geom_violin() + facet_wrap(~variable, scales=free) + geom_jitter(size=0.5) + ylab("Expression (log2)")
  if (!legend) the.plot <- the.plot + guides(fill=F)
  return(the.plot)
}

#' Gene Expression Scatterplot
#' 
#' @importFrom reshape2 melt
#' @importFrom KernSmooth bkde2D
#' 
#' @param object An URD object
#' @param label.x (Character) Value to plot on the x-axis
#' @param label.y (Character) Value to plot on the y-axis
#' @param label.color (Character) Value to plot as point color
#' @param label.x.type (Character) Where to look for the label.x; Default is "search" which looks in order: "meta", "group", "sig", "gene", "pseudotime", "pc", "diff.data"
#' @param label.y.type (Character) Where to look for the label.y; Default is "search" which looks in order: "meta", "group", "sig", "gene", "pseudotime", "pc", "diff.data"
#' @param label.color.type (Character) Where to look for the label.color; Default is "search" which looks in order: "meta", "group", "sig", "gene", "pseudotime", "pc", "diff.data"
#' @param cells (Character vector) Names of cells to include in the plot (Default, \code{NULL} is all cells.)
#' @param point.size (Numeric) Adjust size of points on the plot
#' @param point.alpha (Numeric) Adjust transparency of points on the plot
#' @param title (Character) Title of plot
#' @param xlim (Numeric vector) X-axis limits (Default \code{NULL} lets ggplot2 automatically decide.)
#' @param ylim (Numeric vector) Y-axis limits (Default \code{NULL} lets ggplot2 automatically decide.)
#' @param density.background (Logical) Should 
#' @param density.color (Character) Darkest color of the density scale
#' @param xlab (Character) X-axis label
#' @param ylab (Character) Y-axis label
#' @param colorlab (Character) Color guide label
#' @param legend (Logical) Display a legend
#' 
#' @export

plotScatter <- function(object, label.x, label.y, label.color=NULL, label.x.type="search", label.y.type="search", label.color.type="search", cells=NULL, point.size=1, point.alpha=1, xlim=NULL, ylim=NULL, density.background=T, density.color="#888888", xlab=label.x, ylab=label.y, colorlab=label.color, title="", legend=T) {
  # Get data for plot
  x <- data.for.plot(object, label=label.x, label.type=label.x.type, as.discrete.list=T, cells.use=cells)
  y <- data.for.plot(object, label=label.y, label.type=label.y.type, as.discrete.list=T, cells.use=cells)
  if (!is.null(label.color)) col <- data.for.plot(object, label=label.color, label.type=label.color.type, as.discrete.list=T, cells.use=cells)
  
  # Check for improper data choices
  if (x$discrete | y$discrete) stop ("label.x and label.y must point to continuous data types.")
  
  # Determine x and y limits if not provided
  if (is.null(xlim)) xlim <- range(x$data)
  if (is.null(ylim)) ylim <- range(y$data)
  
  # Build data frame for ggplot
  gg.data <- data.frame(
    x=x$data,
    y=y$data,
    stringsAsFactors=F
  )
  if (!is.null(label.color)) gg.data$color <- col$data
  
  # Build basic ggplot
  the.plot <- ggplot(data=gg.data, aes(x=x, y=y)) + theme_bw() + labs(x=xlab, y=ylab, title=title, color=colorlab) + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  
  # Add density
  if (density.background) {
    # Auto-estimate bandwidth to use, lifted from graphics::smoothScater
    bandwidth <- diff(apply(gg.data[,1:2], 2, stats::quantile, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE))/25
    bandwidth[bandwidth == 0] <- 1
    # Do 2D kernel-density estimate
    kde <- KernSmooth::bkde2D(gg.data[,1:2], bandwidth=bandwidth, gridsize=c(128,128))
    # Transform into gg-ready data.frame
    rownames(kde$fhat) <- kde$x1
    colnames(kde$fhat) <- kde$x2
    kde.gg <- reshape2::melt(kde$fhat)
    names(kde.gg) <- c("x","y","density")
    kde.gg$density <- kde.gg$density ^ 0.25
    # Add it to ggplot
    the.plot <- the.plot + geom_raster(data=kde.gg, aes(x=x,y=y,fill=density)) + scale_fill_gradient(low="#FFFFFF", high=density.color)
  }
  
  # Add points
  if (!is.null(label.color)) {
    the.plot <- the.plot + geom_point(aes(color=color), alpha=point.alpha, size=point.size)
  } else {
    the.plot <- the.plot + geom_point(alpha=point.alpha, size=point.size)
  }
  
  # Legend
  if (legend) {
    the.plot <- the.plot + guides(fill=FALSE)
  } else {
    the.plot <- the.plot + guides(fill=FALSE, color=FALSE)
  }
  
  return(the.plot)

}

