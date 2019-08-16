#' Dot Plot
#' 
#' importFrom reshape2 melt
#' 
#' @param object An URD object
#' @param genes (Character Vector) Genes to plot
#' @param clustering (Character) Name of clustering to use (i.e. a column name of \code{@@group.ids})
#' @param clusters.use (Character vector) Names of specific clusters to plot (default: all clusters)
#' @param min.exp (Numeric) Minimum proportion of expressing cells (0-1) to be shown on the plot
#' @param size.min (Numeric) Point size for clusters with no cells expressing
#' @param size.max (Numeric) Point size for clusters with all cells expressing
#' @param scale.by (Character) Should point size scale by \code{radius} or \code{area}? See \code{\link[ggplot2]{scale_radius}} and \code{\link[ggplot2]{scale_size}}.
#' @param colors (Character vector) Vector of colors to use for a palette.
#' @param mean.expressing.only (Logical) Should mean expression value exclude cells with no expression
#' @param title (Character) How should the plot be titled? (Default: no title)
#' 
#' @return A ggplot2 object
#' 
#' @export
plotDot <- function(object, genes, clustering, clusters.use=NULL, min.exp=.05, size.min=0, size.max=5, scale.by=c("radius", "area"), colors=NULL, mean.expressing.only=F, title="") {
  # Decide on scale function
  if (length(scale.by) > 1) scale.by <- scale.by[1]
  if (tolower(scale.by) == "radius") {
    scale.by <- ggplot2::scale_radius
  } else if (tolower(scale.by) == "area") {
    scale.by <- ggplot2::scale_size
  } else {
    stop("scale.by must be 'radius' or 'area'.")
  }
  
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
  
  # Remove values where proportion expressing is too small and should not be plotted
  gg.data <- gg.data[which(gg.data$Prop.Exp >= min.exp),]
  
  # Turn proportion into %
  gg.data$Prop.Exp <- 100 * gg.data$Prop.Exp

  if (is.null(colors)) colors <- defaultURDContinuousColors(evenly.spaced=T)
  
  # Go ahead and ggplot
    the.plot <- ggplot(data=gg.data, aes(x=Cluster, y=Gene, color=Mean, size=Prop.Exp)) + geom_point() + scale_color_gradientn(colours = colors) + ggtitle(title) + scale_y_discrete(limits=rev(genes)) + scale.by(range=c(size.min,size.max), name="% Exp", breaks=c(20,40,60,80,100), limits=c(0,100))
  if (!is.null(clusters.use)) the.plot <- the.plot + scale_x_discrete(limits=clusters.use)
  return(the.plot)
}

#' Violin plot of gene expression
#' 
#' @importFrom reshape2 melt
#' 
#' @param object An URD object
#' @param labels.plot (Character vector) Data to use for coloring points (e.g. a metadata name, group ID from clustering, or a gene name)
#' @param clustering (Character) Name of clustering to use (i.e. a column name of \code{@@group.ids})
#' @param clusters.use (Character vector) Names of clusters to display (default NULL uses all clusters)
#' @param legend (Logical) Display a legend?
#' @param free.axes (Logical) Should the scale of each label be individually determined?
#' 
#' @return A ggplot2 object
#' 
#' @export
plotViolin <- function(object, labels.plot, clustering, clusters.use=NULL, legend=T, free.axes=F, point.size=0.2, point.color='black', point.alpha=0.5) {
  # Get data.to.plot
  data.plot <- as.data.frame(lapply(labels.plot, function (label) {
    data.for.plot(object = object, label=label)
  }))
  names(data.plot) <- labels.plot
  data.plot$Cluster <- object@group.ids[rownames(data.plot), clustering]
  if (suppressWarnings(!any(is.na(as.numeric(as.character(data.plot$Cluster)))))) {
    cluster.names <- sort(as.numeric(unique(data.plot$Cluster)))
  } else {
    cluster.names <- sort(unique(as.character(data.plot$Cluster)))
  }
  if (!is.null(clusters.use)) data.plot <- data.plot[data.plot$Cluster %in% clusters.use,]
  data.plot$Cell <- rownames(data.plot)
  data.plot.melt <- reshape2::melt(data.plot, id.vars = c("Cell", "Cluster"))
  data.plot.melt$Cluster <- factor(data.plot.melt$Cluster, levels=cluster.names, ordered=T)
  if (free.axes) {free="free"} else {free="fixed"}
  the.plot <- ggplot(data.plot.melt, aes(x=Cluster, y=value, fill=Cluster)) + geom_violin() + facet_wrap(~variable, scales=free) + geom_jitter(size=point.size, color=point.color, alpha=point.alpha) + ylab("Expression (log2)") + theme_bw()
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
#' @param density.bandwidth (Numeric) The kernel bandwidth used in calculating the density background is multiplied by this factor.
#' @param density.bins (Numeric) Number of bins to use for 2D kernel estimation.
#' @param density.exclude.origin (Logical) Whether values with expression (0,0) should be excluded from the density calculation. (The number of non-expressing cells can often overwhelm the density calculation and make it useless if the genes plotted are specific to a particular lineage and \code{cells} is not used to restrict the cells on the plot.)
#' @param xlab (Character) X-axis label
#' @param ylab (Character) Y-axis label
#' @param colorlab (Character) Color guide label
#' @param legend (Logical) Display a legend
#' 
#' @export

plotScatter <- function(object, label.x, label.y, label.color=NULL, label.x.type="search", label.y.type="search", label.color.type="search", cells=NULL, point.size=1, point.alpha=1, xlim=NULL, ylim=NULL, density.background=T, density.color="#888888", density.bandwidth=1, density.bins=128, density.exclude.origin=T, xlab=label.x, ylab=label.y, colorlab=label.color, title="", legend=T) {
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
    if (density.exclude.origin) {
      gg.density <- gg.data[gg.data$x > 0 | gg.data$y > 0,]
    } else {
      gg.density <- gg.data
    }
    # Auto-estimate bandwidth to use, lifted from graphics::smoothScater
    bandwidth <- diff(apply(gg.density[,1:2], 2, stats::quantile, probs = c(0.05, 0.95), na.rm = TRUE, names = FALSE))/30*density.bandwidth
    bandwidth[bandwidth == 0] <- 1
    # Do 2D kernel-density estimate
    kde <- KernSmooth::bkde2D(gg.density[,1:2], bandwidth=bandwidth, gridsize=c(density.bins,density.bins))
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

#' Plot distributions
#' 
#' This will plot the distribution of a continuous label across discrete classes.
#' For instance, it can be used to plot the distribution of pseudotime in each
#' developmental stage.
#' 
#' @param object An URD object
#' @param label (Character) Data to use for distributions (e.g. "pseudotime", a metadata name, or a gene name). Must be continuous.
#' @param category.label (Character) Data to use for dividing into separate curves (e.g. developmental stage or a clustering). Must be discrete.
#' @param label.type (Character) Type of data to search for the label. Default is "search" which checks several data types in order. For more information: \code{\link{data.for.plot}}
#' @param category.label.type (Character) Type of data to search for category.label. Default is "search" which checks several data types in order. For more information: \code{\link{data.for.plot}}
#' @param legend (Logical) Show a legend?
#' @param plot.title (Character) Title of the plot
#' 
#' @examples 
#' 
#' @return A ggplot2 object. If used in a loop, plot(plotDists(...)) must be used.
#' 
#' @export
plotDists <- function(object, label, category.label, label.type="search", category.label.type="search", legend=T, plot.title="") {
  # Get the data
  label.data <- data.for.plot(object, label = label, label.type = label.type, as.discrete.list = T)
  category.data <- data.for.plot(object, label = category.label, label.type=category.label.type, as.discrete.list=T)

  # Make sure label is continuous and category is discrete
  if (label.data$discrete | !category.data$discrete) stop("label must be a continuous variable and category.label must be a discrete variable.")
  
  # Make a data frame for ggplot.
  gg.data <- data.frame(
    label=label.data$data,
    category=category.data$data,
    stringsAsFactors=F,
    row.names = names(label.data$data)
  )
  
  # Make the plot
  the.plot <- ggplot(data=gg.data, aes(x=label, color=category, fill=category)) + geom_density(alpha=0.4) + theme_bw() + labs(x=label, fill=category.label, color=category.label, title=plot.title)
  
  # Get rid of legends if desired
  if (!legend) the.plot <- the.plot + guides(color=F, fill=F)
  
  # Return the plot
  return(the.plot)
}


