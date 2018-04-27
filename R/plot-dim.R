#' Dimensionality Reduction Plot
#' 
#' Plots cells according to their coordinates in a dimensionality reduction (tSNE by default,
#' but also PCA or diffusion map). Cells are colored according to a user set \code{label} that
#' can range from gene expression, to metadata values, or cluster identity. See 
#' \code{\link{data.for.plot}} for more information about labels that can be chosen. Additionally, 
#' see \code{\link{plotDimDual}} to plot two
#' continuous variables simultaneously, \code{\link{plotDimHighlight}} to highlight one group
#' from a discrete label, \code{\link{plotDimArray}} to repeat the same plot across several
#' sets of dimensions, and \code{\link{plotDim3D}} to plot in three dimensions.. Additionally, 
#' \code{transitions.plot} can plot the connections from the diffusion map onto the plot.
#' 
#' @param object An URD object
#' @param label (Character) Data to use for coloring points (e.g. a metadata name, group ID from clustering, or a gene name)
#' @param label.type (Character) Type of data to search for the label. Default is "search" which checks several data types in order. For more information: \code{\link{data.for.plot}}
#' @param reduction.use (Character) Dimensionality reduction to use (tSNE, PCA, or Diffusion Map)
#' @param dim.x (Numeric) Component to use on x-axis
#' @param dim.y (Numeric) Component to use on y-axis
#' @param colors (Character vector) Colors to use to generate a gradient scale to color continuous data
#' @param discrete.colors (Character vector) Colors to use to color 
#' @param point.size (Numeric) Size of points on plot
#' @param alpha (Numeric) Transparency of points on plot: 0 (Transparent) - 1 (Opaque)
#' @param point.shapes (Logical) Should point shapes vary? This is useful in plots of discrete data types with many categories to help differentiate between similar colors.
#' @param plot.title (Character) Title of the plot
#' @param legend (Logical) Show a legend?
#' @param legend.title (Character) Should the legend be titled?
#' @param legend.point.size (Numeric) How big should points be in the legend?
#' @param label.clusters (Logical) Label centroids of a discrete label?
#' @param cells (Character vector) Cells to show on the plot (Default \code{NULL} is all cells.)
#' @param x.lim (Numeric) Limits of x-axis (NULL autodetects)
#' @param y.lim (Numeric) Limits of y-axis (NULL autodetects)
#' @param color.lim (Numeric) Limits of the point color scale (NULL autodetects)
#' @param na.rm (Logical) Should points with value NA for the desired data be removed from the plot?
#' @param transitions.plot (Numeric or NULL) Number of transition matrix connections to add to the plot. \code{NULL} will plot all connections. (WARNING: Too many connections will produce an unreadable plot that takes a long time to plot. Start with 10,000.)
#' @param transitions.alpha (Numeric) Maximum transparency of line segments representing transitions. (They are scaled based on their transition probability).
#' @param transitions.df (data.frame) Output from \link{edgesFromDM} (potentially further curated) to display on the plot. If provided, \code{transitions.plot} is ignored and all transitions in the provided data.frame are plotted.
#' 
#' @return A ggplot2 object
#' 
#' @export
plotDim <- function(object, label, label.type="search", reduction.use=c("tsne", "pca", "dm"), dim.x=1, dim.y=2, colors=NULL, discrete.colors=NULL, point.size=1, alpha=1, point.shapes=F, plot.title=label, legend=T, legend.title="", legend.point.size=3*point.size, label.clusters=F, cells=NULL, x.lim=NULL, y.lim=NULL, color.lim=NULL, na.rm=F, transitions.plot=0, transitions.alpha=0.5, transitions.df=NULL) {
  
  # Get the data to plot
  if (length(reduction.use) > 1) reduction.use <- reduction.use[1]
  if (tolower(reduction.use)=="tsne") {
    data.plot <- object@tsne.y
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2]) stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("tSNE", dim.x)
    dim.y <- paste0("tSNE", dim.y)
  } else if (tolower(reduction.use)=="pca") {
    data.plot <- object@pca.scores
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2]) stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("PC", dim.x)
    dim.y <- paste0("PC", dim.y)
    data.plot <- data.plot[,c(dim.x, dim.y)]
  } else if (tolower(reduction.use)=="dm") {
    data.plot <- object@dm@eigenvectors
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2]) stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("DC", dim.x)
    dim.y <- paste0("DC", dim.y)
    data.plot <- as.data.frame(data.plot[,c(dim.x, dim.y)])
  } else {
    stop("The reduction provided is invalid.")
  }
  
  # Get the info to color by
  sig.score <- data.for.plot(object, label = label, label.type = label.type, as.color = F, as.discrete.list = T)
  data.plot$SIG <- sig.score[[2]][rownames(data.plot)]
  
  # Remove NAs if desired
  if (na.rm) {
    data.plot <- data.plot[complete.cases(data.plot),]
  }
  
  # Limit cells if desired
  if (!is.null(cells)) {
    cells <- intersect(cells, rownames(data.plot))
    data.plot <- data.plot[cells,]
  }
  
  # Get transitions if desired
  if (is.null(transitions.plot) || transitions.plot > 0 || !is.null(transitions.df)) {
    # If transitions aren't provided, get edge list
    if (is.null(transitions.df)) transitions.df <- edgesFromDM(object, cells=rownames(data.plot), edges.return=transitions.plot)
    # Add coordinates
    transitions.df$x1 <- data.plot[transitions.df$from, dim.x]
    transitions.df$x2 <- data.plot[transitions.df$to, dim.x]
    transitions.df$y1 <- data.plot[transitions.df$from, dim.y]
    transitions.df$y2 <- data.plot[transitions.df$to, dim.y]
    # Normalize alpha
    transitions.df$alpha <- transitions.df$weight / max(transitions.df$weight) * transitions.alpha
  }
  
  # Start the plot
  this.plot <- ggplot(data=data.plot, aes_string(x=dim.x, y=dim.y))
  
  # Add the transitions if desired
  if (!is.null(transitions.df)) this.plot <- this.plot + geom_segment(inherit.aes=F, data=transitions.df, aes(x=x1, y=y1, xend=x2, yend=y2, alpha=alpha))
  
  # Add the points (color based on whether or not label is discrete)
  if (sig.score[[1]]) {
    # Discrete
    if (point.shapes) {
      shape.rep <- ceiling(length(unique(data.plot$SIG)) / 4) + 1
      #this.plot <- this.plot + geom_point(aes(color=SIG, shape=SIG), size=point.size, alpha=alpha) + scale_shape_manual(values=rep(c(15, 17, 18, 19), shape.rep))
      this.plot <- this.plot + geom_point(aes(color=SIG, shape=SIG), size=point.size, alpha=alpha) + scale_shape_manual(values=rep(c(0, 2, 8, 9), shape.rep))
    } else {
      this.plot <- this.plot + geom_point(aes(color=SIG), size=point.size, alpha=alpha, stroke=0)
    }
    if (!is.null(discrete.colors)) {
      this.plot <- this.plot + scale_color_manual(values=discrete.colors)
    }
  } else {
    # Continuous
    if (is.null(colors)) colors <- defaultURDContinuousColors()
    this.plot <- this.plot + geom_point(aes(color=SIG), size=point.size) + 
      scale_color_gradientn(colors=colors, limits=color.lim)
  }
  
  # Label/title things appropriately
  this.plot <- this.plot + labs(title=plot.title, color=legend.title, shape=legend.title)
  
  # Format it to your liking.
  this.plot <- this.plot + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title=element_text(face="bold"))
  
  # Label clusters/groups if desired
  if (label.clusters && sig.score[[1]]) {
    # Get info about clusters
    data.plot$CLUSTER <- data.plot$SIG
    # Calculate center of each cluster
    k.centers <- aggregate(data.plot[,1:2], by=list(data.plot$CLUSTER), FUN="mean")
    # Add labels
    this.plot <- this.plot + geom_label(data=k.centers, aes_string(x=dim.x, y=dim.y, label="Group.1"), color="black", alpha=0.6, show.legend = F)
  }
  
  # Remove legend if desired
  if (!legend) {
    this.plot <- this.plot + guides(color=FALSE, shape=FALSE)
  } else if (sig.score[[1]]) {
    # Otherwise, make the legend points bigger if coloring by a discrete value
    this.plot <- this.plot + guides(color=guide_legend(override.aes = list(size=legend.point.size)))
  }
  # Add limits if desired
  this.plot <- this.plot + guides(alpha=F)
  if (!is.null(x.lim)) this.plot <- this.plot + xlim(x.lim[1],x.lim[2])
  if (!is.null(y.lim)) this.plot <- this.plot + ylim(y.lim[1],y.lim[2])
  return(this.plot)
}

#' Dimensionality Reduction Plot (Dual Color)
#' 
#' Plots cells according to their coordinates in a dimensionality reduction (tSNE by default).
#' Cells are colored according to two user set labels (\code{label.red} and \code{label.green})
#' in a mode that simulates dual-color microscopy - cells that express neither label are black,
#' cells that express one label are red or green, and cells that express both labels are yellow.
#' Both labels must be continuous variables (i.e. not cluster identities).
#' 
#' @importFrom scales squish rescale
#' @importFrom stats quantile
#' 
#' @param object An URD object
#' @param label.red (Character) Data to use for coloring points for the red channel
#' @param label.green (Character) Data to use for coloring points for the green channel
#' @param label.red.type (Character) Type of data to search for the label for the red channel. Default is "search" which checks several data types in order. For more information: \code{\link{data.for.plot}}
#' @param label.green.type (Character) Type of data to search for the label for the green channel. Default is "search" which checks several data types in order. For more information: \code{\link{data.for.plot}}
#' @param reduction.use (Character) Dimensionality reduction to use (tSNE, PCA, or Diffusion Map)
#' @param dim.x (Numeric) Component to use on x-axis
#' @param dim.y (Numeric) Component to use on y-axis
#' @param colors (Character vector) Colors to use to generate a gradient scale to color continuous data
#' @param discrete.colors (Character vector) Colors to use to color 
#' @param point.size (Numeric) Size of points on plot
#' @param alpha (Numeric) Transparency of points on plot: 0 (Transparent) - 1 (Opaque)
#' @param plot.title (Character) Title of the plot
#' @param legend (Logical) Show a legend?
#' @param legend.size (Numeric) Adjusts the size of the legend.
#' @param legend.offset.x (Numeric) Adjust the legend position (in terms of dimensionality reduction coordinates)
#' @param legend.offset.y (Numeric) Adjust the legend position (in terms of dimensionality reduction coordinates)
#' @param label.clusters (Logical) Label centroids of a discrete label?
#' @param x.lim (Numeric) Limits of x-axis (NULL autodetects)
#' @param y.lim (Numeric) Limits of y-axis (NULL autodetects)
#' @param na.rm (Logical) If \code{TRUE}, points with an NA value for either label are displayed as transparent grey. If \code{FALSE}, they are removed from the plot.
#' @param na.alpha (Numeric) If \code{na.rm=FALSE}, thae alpha value that should be used for NA points
#' 
#' @return A ggplot2 object
#' 
#' @export
plotDimDual <- function(object, label.red, label.green, label.red.type="search", label.green.type="search", reduction.use=c("tsne", "pca", "dm"), dim.x=1, dim.y=2, point.size=1, alpha=1, plot.title="", legend=T, legend.size=1/5.5, legend.offset.x=0, legend.offset.y=0, label.clusters=F, x.lim=NULL, y.lim=NULL, na.rm=F, na.alpha=0.4 * alpha) {
  
  # Get the data to plot
  if (length(reduction.use) > 1) reduction.use <- reduction.use[1]
  if (tolower(reduction.use)=="tsne") {
    data.plot <- object@tsne.y
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2]) stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("tSNE", dim.x)
    dim.y <- paste0("tSNE", dim.y)
  } else if (tolower(reduction.use)=="pca") {
    data.plot <- object@pca.scores
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2]) stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("PC", dim.x)
    dim.y <- paste0("PC", dim.y)
    data.plot <- data.plot[,c(dim.x, dim.y)]
  } else if (tolower(reduction.use)=="dm") {
    data.plot <- object@dm@eigenvectors
    if (dim.x > dim(data.plot)[2] | dim.y > dim(data.plot)[2]) stop("Dimensions requested were not previously calculated.")
    dim.x <- paste0("DC", dim.x)
    dim.y <- paste0("DC", dim.y)
    data.plot <- as.data.frame(data.plot[,c(dim.x, dim.y)])
  } else {
    stop("The reduction provided is invalid.")
  }
  
  # Get label expression
  plot.red <- data.for.plot(object, label=label.red, label.type = label.red.type, as.color = F, as.discrete.list = T)
  plot.green <- data.for.plot(object, label=label.green, label.type = label.green.type, as.color = F, as.discrete.list = T)
  if (plot.red[[1]] | plot.green[[1]]) stop("Cannot use discrete labels in dual-color plots.")
  data.plot$gene.red <- plot.red[[2]][rownames(data.plot)]
  data.plot$gene.green <- plot.green[[2]][rownames(data.plot)]
  
  # Scale gene expression and generate colors
  gene.red.max <- quantile(data.plot$gene.red[data.plot$gene.red > 0], prob=0.975, na.rm=T)
  gene.green.max <- quantile(data.plot$gene.green[data.plot$gene.green > 0], prob=0.975, na.rm=T)
  data.plot$gene.red.scaled <- squish(rescale(data.plot$gene.red, from=c(0,gene.red.max)), c(0,1))
  data.plot$gene.green.scaled <- squish(rescale(data.plot$gene.green, from=c(0,gene.green.max)), c(0,1))
  cc <- which(complete.cases(data.plot))
  data.plot[cc,"color.plot"] <- rgb(data.plot[cc,"gene.red.scaled"], data.plot[cc,"gene.green.scaled"], 0)
  
  # Add alpha
  data.plot$alpha <- alpha
  
  # Deal with NAs.
  if (!na.rm) {
    noncc <- setdiff(1:nrow(data.plot), cc)
    data.plot[noncc, "color.plot"] <- "#CECECE"
    data.plot[noncc, "alpha"] <- na.alpha
  }
  
  # Do the plot
  this.plot <- ggplot()
  
  # Add gene signature points
  this.plot <- this.plot + geom_point(data=data.plot, aes_string(x=dim.x, y=dim.y), color=data.plot$color.plot, size=point.size, alpha=data.plot$alpha) + guides(color=FALSE)
  
  # Label/title things appropriately
  this.plot <- this.plot + labs(title=plot.title)
  
  # Format it to your liking.
  this.plot <- this.plot + theme_bw() + theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), plot.title=element_text(face="bold"))
  
  # Label clusters if desired
  if (label.clusters) {
    # Get info about clusters
    data.plot$CLUSTER <- object@group[rownames(data.plot)]
    # Calculate center of each cluster
    k.centers <- aggregate(data.plot[,1:2], by=list(data.plot$CLUSTER), FUN="mean")
    # Add labels
    this.plot <- this.plot + geom_label(data=k.centers, aes_string(x=dim.x, y=dim.y, label="Group.1"), color="black", alpha=0.6, show.legend = F)
  }
  
  # Add color legend
  if (legend) {
    leg.x.max <- max(data.plot[,1]) + legend.offset.x
    leg.y.max <- max(data.plot[,2]) + legend.offset.y
    leg.x.min <- leg.x.max - (diff(range(data.plot[,1])) * legend.size) + legend.offset.x
    leg.y.min <- leg.y.max - (diff(range(data.plot[,2])) * legend.size) + legend.offset.y
    greens <- round(seq(0,gene.green.max, length.out = 6), digits=2)
    reds <- round(seq(0,gene.red.max, length.out = 6), digits=2)
    leg.x.breaks <- seq(leg.x.min, leg.x.max, length.out=10)
    leg.y.breaks <- seq(leg.y.min, leg.y.max, length.out=10)
    legend.squares <- data.frame(stringsAsFactors=F,
                                 r=rep(c(0,0.2,0.4,0.6,0.8,1.0), each=6),
                                 g=rep(c(0,0.2,0.4,0.6,0.8,1.0), 6),
                                 x.1=rep(leg.x.breaks[2:7], each=6),
                                 x.2=rep(leg.x.breaks[3:8], each=6),
                                 y.1=rep(leg.y.breaks[2:7], 6),
                                 y.2=rep(leg.y.breaks[3:8], 6)
    )
    legend.squares$color.plot <- rgb(legend.squares$r, legend.squares$g, 0)
    this.plot <- this.plot + geom_rect(data=legend.squares, aes(xmin=legend.squares$x.1, xmax=legend.squares$x.2, ymin=legend.squares$y.1, max=legend.squares$y.2), fill=legend.squares$color.plot)
    side.labels <- data.frame(stringsAsFactors=F,
                              label.x=c(leg.x.breaks[2:7],rep(leg.x.breaks[8], 6)),
                              label.y=c(rep(leg.y.breaks[8], 6), leg.y.breaks[2:7]),
                              label.text=c(reds,greens)
    )
    side.labels$label.x <- side.labels$label.x + (0.5*(leg.x.breaks[2]-leg.x.breaks[1]))
    side.labels$label.y <- side.labels$label.y + (0.5*(leg.y.breaks[2]-leg.y.breaks[1]))
    side.labels$angle <- c(rep(90,6), rep(0,6)); side.labels$hjust <- 0; side.labels$vjust <- 0.5; side.labels$size <- 2
    side.labels <- rbind(side.labels, c(leg.x.breaks[5], leg.y.breaks[1], 0, 0, 0.5, 0.5, 3))
    side.labels <- rbind(side.labels, c(leg.x.breaks[1], leg.y.breaks[5], 0, 90, 0.5, 0.5, 3))
    side.labels[13:14, "label.text"] <- c(label.red, label.green)
    this.plot <- this.plot + geom_text(data=side.labels, aes(label=label.text, x=label.x, y=label.y), angle=side.labels$angle, hjust=side.labels$hjust, vjust=side.labels$vjust, size=side.labels$size)
    
    
  }
  # Add limits if desired
  if (!is.null(x.lim)) this.plot <- this.plot + xlim(x.lim[1],x.lim[2])
  if (!is.null(y.lim)) this.plot <- this.plot + ylim(y.lim[1],y.lim[2])
  return(this.plot)
}

#' Dimensionality Reduction Plot With Highlighted Clusters
#' 
#' Produces a plot with \code{\link{plotDim}} with cluster colors dimmed, and
#' a single cluster highlighted in a bright color (by default, red).
#' 
#' @param object An URD object
#' @param clustering (Character) Name of column in \code{@@group.ids} that identifies the clustering to pull from
#' @param cluster (Character vector) A cluster name that you want to highlight.
#' @param highlight.color (Character) A color name to use for the highlighted cluster.
#' @param ... Additional parameters to pass to \code{\link{plotDim}}
#' 
#' @return A ggplot2 object
#' 
#' @export
plotDimHighlight <- function(object, clustering, cluster, highlight.color='red', ...) {
  cells <- cellsInCluster(object, clustering = clustering, cluster=cluster)
  base.plot <- plotDim(object, label=clustering, ...) + scale_color_discrete(c=40)
  base.plot <- base.plot + geom_point(data=base.plot$data[which(base.plot$data$SIG %in% cluster),], aes_string(x=names(base.plot$data)[1], y=names(base.plot$data)[2]), color=highlight.color)
  base.plot$labels$title <- paste0(base.plot$labels$title, " (Highlight ", paste0(cluster, collapse=", "), ")")
  return(base.plot)
}

identify.cluster <- function(object, clustering, n=20, point.size=0.5, color.palette=rainbow) {
  color.factor <- as.factor(object@group.ids[rownames(object@tsne.y),clustering])
  color.values <- color.palette(length(levels(color.factor)))
  names(color.values) <- levels(color.factor)
  color.plot <- color.values[color.factor]
  par(mar=c(0,0,0,0))
  plot(x=object@tsne.y$tSNE1, y=object@tsne.y$tSNE2, pch=16, cex=point.size, col=color.plot, xaxt='n', yaxt='n', ann=F)
  print("Click on the plot to identify cluster IDs of nearby cells! Esc to abort.")
  click <- locator(n = 1)
  par(mar=c(5.1,4.1,4.1,2.1))
  tsne.dist <- object@tsne.y
  names(tsne.dist) <- c("x","y")
  tsne.dist$x <- abs(tsne.dist$x - click$x) 
  tsne.dist$y <- abs(tsne.dist$y - click$y)
  tsne.dist$d <- sqrt(tsne.dist$x^2 + tsne.dist$y^2)
  tsne.dist$c <- object@group.ids[rownames(tsne.dist),clustering]
  sort(table(tsne.dist[order(tsne.dist$d)[1:n], "c"]), decreasing=T)
}

#' Dimensionality Reduction Plot Array
#' 
#' For surveying several dimensions of a dimensionality reduction (as in
#' \code{\link{plotDim}}). This will produce the same plot, but varying
#' the pair of dimensions shown. It optionally saves directly to a file,
#' since this can produce large plots.
#' 
#' @param object An URD object
#' @param reduction.use (Character) Dimensionality reduction to use
#' @param dims.to.plot (Numeric vector) Dimensions to plot. Will be plotted as pairs, with odd indices on the x-axes and even indices on the y-axes.
#' @param outer.title (Charater) A title to place outside the entire array of plots.
#' @param file (Character) Path to plot to save (if NULL, plot is returned)
#' @return If \code{file=NULL}, returns a grid.array, otherwise returns nothing.
#' @export
plotDimArray <- function(object, reduction.use=c("dm", "pca"), dims.to.plot, outer.title=NULL, file=NULL, file.width=750, file.height=600, ...) {
  # Need an even number of dimensions
  if (length(dims.to.plot) %% 2 != 0) dims.to.plot <- head(dims.to.plot, -1)
  # Check that all of the relevant dimensions were calculated.
  if (tolower(reduction.use)=="pca") {
    if(!all(dims.to.plot < ncol(object@pca.scores))) stop("dims.to.plot referenced PCs that were not calculated.")
  } else if (tolower(reduction.use)=="dm") {
    if(!all(dims.to.plot < ncol(object@dm@eigenvectors))) stop("dims.to.plot referenced DCs that were not calculated.")
  } else {
    stop("reduction.use must be either 'pca' or 'dm'.")
  }
  # Calculate grid layout
  n.cols <- ceiling(sqrt(length(dims.to.plot)/2))
  n.rows <- ceiling(length(dims.to.plot)/(2*n.cols))
  # Do pair plots
  the.plots <- lapply(seq(1, length(dims.to.plot), 2), function(dim.n) {
    plotDim(object, reduction.use = reduction.use, dim.x = dims.to.plot[dim.n], dim.y=dims.to.plot[dim.n+1], ...)
  })
  # Either return the plot or save it directly to a PNG
  if (is.null(file)) {
    return(grid.arrange(grobs=the.plots, n.cols=n.cols, top=outer.title))
  } else {
    png(file=file, width=file.width*n.cols, height=file.height*n.rows)
    grid.arrange(grobs=the.plots, n.cols=n.cols, top=outer.title)
    dev.off()
  }
}