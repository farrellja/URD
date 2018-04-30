#' Generate force-directed layout using tip-walk data.
#' 
#' This function generates a k-nearest neighbor network for cells, based on their visitation
#' frequency by the biased random walks from different tips (and optionally pseudotime), and
#' uses it as input into a force-directed layout (powered by igraph). The force-directed
#' layout is generated in 2 dimensions, and pseudotime is used as a third dimension.
#' 
#' Several settings adjust the k-nearest neighbor network to assist in the layout. Outlier
#' cells (that are only closely connected to a single other cells) or outlier edges (with
#' unusually long distances) can be eliminated, though these parameters are disabled by
#' default. However, the dendrogram structure recovered by URD is used by default to
#' refine the k-nearest neighbor network, breaking links between cells that are distant
#' in the dendrogram. This emphasizes producing a parseable layout at the expense of
#' representing rare transitions in the data.
#' 
#' @importFrom RANN nn2
#' @importFrom stats quantile
#' @importFrom reshape2 melt
#' @importFrom igraph graph_from_data_frame layout_with_fr layout_with_drl layout_with_kk vcount V distances
#' 
#' @param object An URD object
#' @param num.nn (Numeric) Number of nearest neighbors to use. (\code{NULL} will use the square root of the number of cells as default.)
#' @param pseudotime (Character) Pseudotime to use in nearest-neighbor calculation (i.e. a column name of \code{@@pseudotime}). (\code{NULL} omits pseudotime from the distance calculation of the nearest neighbor network.)
#' @param method (Character: "fr", "drl", or "kk") Which force-directed layout algorithm to use (\link[igraph:layout_with_fr]{Fruchterman-Reingold}, \link[igraph:layout_with_drl]{DrL}, or \link[igraph:layout_with_kk]{Kamada-Kawai})
#' @param cells.to.do (Character vector) Cells to use in the layout (default \code{NULL} is all cells in the tree.)
#' @param cut.outlier.cells (Numeric) If desired, omit cells with unusual second nearest neighbor distances (i.e. those that are likely outliers only well connected to one other cell). Parameter is given as a factor of the interquartile range calculated across all cells' distance to their second nearest neighbor. (Default is not to omit any cells.)
#' @param cut.outlier.edges (Numeric) If desired, cut edges in the nearest neighbor graph with unusually long distances. Parameters is given as a factor of the interquartile range calculated across all edges in the graph. (Default is not to cut any edges based on their length.)
#' @param max.pseudotime.diff (Numeric) If desired, cut edges in the nearest neighbor graph between cells with longer difference in pseudotime. (Default is not to cut any edges based on their pseudotime.)
#' @param cut.unconnected.segments (Numeric) Cut connections in the nearest-neighbor graph to cells that are 
#' more segments away in the dendrogram structure. For instance, the most aggressive setting (1) will only 
#' maintain that are at most 1 segment away (so it will maintain connections within a segment and to that 
#' segment's parent or children. Higher values permit more distant connections. The default value is 2, 
#' which permits connections up to two segments away (i.e. connections within a segment, to that segment's 
#' parent, grandparent, children, grandchildren, and siblings.) \code{NA} or \code{NULL} disable this setting
#' and permit all connections.
#' @param min.final.neighbors (Numeric) After trimming outlier and unconnected connections in the nearest
#' neighbor graph, remove any cells that remain connected to fewer than this many other cells.
#' @param tips (Character vector) Tips for which walk visitation data should be used in the construction of the nearest neighbor graph. (Default is all tips )
#' @param coords (Matrix: Cells as rows, 2 columns) Starting coordinates for the force directed layout. Default (\code{"auto"}) takes them from the cell layout of the dendrogram.
#' @param start.temp (Numeric) Starting temperature for the force-directed layout (if \code{method="fr"}), which controls how much cells can move in the initial iterations of the algorithm. Default (\code{NULL}) is the square root of the number of cells.
#' @param density.neighbors (Numeric) Distance to this nearest neighbor (default is 10th nearest neighbor) is used as a proxy for local density in the force-directed layout. This is used by \code{\link{plotTreeForce}} if \code{density.alpha=T} for increasing transparency in more high density regions of the layout. 
#' @param plot.outlier.cuts (Logical) If \code{cut.outlier.cells=T} or \code{cut.outlier.edges=T}, this displays the 
#' @param verbose (Logical) Print progress and time stamps?
#' 
#' @export

treeForceDirectedLayout <- function(object, num.nn=NULL, pseudotime, method=c("fr", "drl", "kk"), cells.to.do=NULL, cut.outlier.cells=NULL, cut.outlier.edges=NULL, max.pseudotime.diff=NULL, cut.unconnected.segments=2, min.final.neighbors=2, tips=object@tree$tips, coords="auto", start.temp=NULL, density.neighbors=10, plot.outlier.cuts=F, verbose=F) {
  # Params
  if (length(method) > 1) method <- method[1]
  if (is.null(cells.to.do)) cells.to.do <- rownames(object@diff.data)
  starting.cells <- length(cells.to.do)
  if (is.null(num.nn)) num.nn <- ceiling(sqrt(length(cells.to.do)))
  if (is.null(tips)) tips <- object@tree$tips
  if (coords=="auto") {
    coords <- as.matrix(object@tree$cell.layout[cells.to.do,c("x","y")])
  } else if (!is.null(coords) && class(coords) != "matrix") {
    stop("coords must either be 'auto', NULL, or a matrix.")
  }
  dim=2
  
  if (verbose) print(paste(Sys.time(), ": Starting with parameters", method, num.nn, "NN", dim, "D", length(cells.to.do), "cells"))
  
  # Get and normalize walk data, then add pseudotime (unless NULL in which case don't use)
  if (verbose) print(paste0(Sys.time(), ": Preparing walk data."))
  walk.data <- object@diff.data[cells.to.do,paste0("visitfreq.raw.", object@tree$tips)]
  walk.total <- apply(walk.data, 1, sum)
  walk.data <- sweep(walk.data, 1, walk.total, "/")
  if (!is.null(pseudotime)) walk.data$pseudotime <- object@pseudotime[cells.to.do, pseudotime]
  walk.data <- as.matrix(walk.data)
  
  # Calculate nearest neighbor graph
  if (verbose) print(paste0(Sys.time(), ": Calculating nearest neighbor graph."))
  walk.nn <- RANN::nn2(data=walk.data,query=walk.data,k=max(num.nn)+1, treetype = "kd", searchtype="priority")

  rownames(walk.nn$nn.idx) <- rownames(walk.data)
  walk.nn$nn.idx <- walk.nn$nn.idx[,2:(num.nn+1)]
  rownames(walk.nn$nn.dists) <- rownames(walk.data)
  walk.nn$nn.dists <- walk.nn$nn.dists[,2:(num.nn+1)]
  walk.nn$nn.label <- apply(walk.nn$nn.idx, 2, function(y) rownames(walk.data)[y])
  rownames(walk.nn$nn.label) <- rownames(walk.data)
  
  # Cull outlier cells with unusual second nearest neighbor distances
  if (!is.null(cut.outlier.cells)) {
    q <- quantile(walk.nn$nn.dists[,2])
    outer.fence <- q[4] + cut.outlier.cells*(q[4]-q[1])
    if (plot.outlier.cuts) {
      hist(walk.nn$nn.dists[,2], breaks=100, main="Outlier Removal: Second Nearest Neighbor", xlab="Distance to 2nd Nearest Neighbor")
      abline(v=outer.fence, col='red')
    }
    if (verbose) print(paste0(Sys.time(), ": Removing ", round(length(which(walk.nn$nn.dists[,2] > outer.fence)) / starting.cells, digits=2), "% of cells as outliers."))
    cells.keep <- which(walk.nn$nn.dists[,2] <= outer.fence)
  } else {
    cells.keep <- 1:dim(walk.data)[1]
  }
  
  # Turn nearest neighbors into edge list
  if (verbose) print(paste0(Sys.time(), ": Preparing edge list."))
  edges <- melt(walk.nn$nn.label[cells.keep,], stringsAsFactors=F)
  dists <- melt(walk.nn$nn.dists[cells.keep,])
  edges <- edges[,c(1,3)]
  names(edges) <- c("V1", "V2")
  edges$dists <- dists$value
  edges$V1 <- as.character(edges$V1)
  edges$V2 <- as.character(edges$V2)
  
  starting.edges <- dim(edges)[1]
  
  # Cut connections with unusually long distances
  if (!is.null(cut.outlier.edges)) {
    q <- quantile(edges$dists)
    outer.fence <- q[4] + cut.outlier.edges*(q[4]-q[1])
    if (verbose) print(paste0(Sys.time(), ": Removing ", round(length(which(edges$dists > outer.fence)) / starting.edges * 100, digits=3), "% of edges as outliers."))
    if (plot.outlier.cuts) {
      hist(edges$dists, breaks=100, main="Outlier Removal: Edge distances", xlab="Edge distance")
      abline(v=outer.fence, col='red')
    }
    edges <- edges[edges$dists <= outer.fence,]
  }
  
  # Cut connections with overly large pseudotime differences
  if (!is.null(max.pseudotime.diff)) {
    edges$pt1 <- object@pseudotime[edges$V1,pseudotime]
    edges$pt2 <- object@pseudotime[edges$V2,pseudotime]
    edges$dpt <- abs(edges$pt1 - edges$pt2)
    if (verbose) print(paste0(Sys.time(), ": Removing ", round(length(which(edges$dpt > max.pseudotime.diff)) / starting.edges * 100, digits=3), "% of edges with too large pseudotime difference."))
    edges <- edges[which(edges$dpt <= max.pseudotime.diff),]
  }
  
  # Cut edges between overly distant (e.g. non-connected) segments
  if (!is.na(cut.unconnected.segments) & !is.null(cut.unconnected.segments) & cut.unconnected.segments > 0) {
    # Determine which segments each connected cell belongs to
    edges$seg1 <- object@diff.data[edges$V1, "segment"]
    edges$seg2 <- object@diff.data[edges$V2, "segment"]
    # Determine distance in the tree between each segment
    seg.dist <- igraph::distances(object@tree$tree.igraph, mode="all")
    edges$seg.dist <- apply(edges[,c("seg1","seg2")], 1, function(segs) seg.dist[segs[1],segs[2]])
    
    if (verbose) print(paste0(Sys.time(), ": Removing ", round(length(which(edges$seg.dist > cut.unconnected.segments)) / starting.edges * 100, digits=2), "% of edges that are between segments with distance > ", cut.unconnected.segments))
    edges <- edges[edges$seg.dist <= cut.unconnected.segments,]
  }
  
  # Trim cells that are no longer well connected to the graph
  if (verbose) print(paste0(Sys.time(), ": Trimming cells that are no longer well connected."))
  connections.remaining <- table(edges$V1)
  cells.without.enough.connections <- sum(connections.remaining < min.final.neighbors)
  cells.with.enough.connections <- names(which(connections.remaining >= min.final.neighbors))
  while(cells.without.enough.connections > 0) {
    cells.with.enough.connections <- names(which(connections.remaining >= min.final.neighbors))
    edges <- edges[which(edges$V1 %in% cells.with.enough.connections & edges$V2 %in% cells.with.enough.connections),]
    connections.remaining <- table(edges$V1)
    cells.without.enough.connections <- sum(connections.remaining < min.final.neighbors)
  }
  if (verbose) print(paste0(Sys.time(), ": ", round(length(cells.with.enough.connections)/starting.cells*100, digits=2), "% of starting cells preserved."))
  
  # Create igraph representation
  if (verbose) print(paste0(Sys.time(), ": Preparing igraph object."))
  if (method=="kk") {
    # kk requires that weights be higher for points that should be far apart -- use distance directly
    edges$weight <- edges$dists
  } else {
    # everything else places points with high weight nearby -- use inverse of distance as graph weight
    edges$weight <- 1/edges$dists
  }
  object@tree$walks.force.edges <- edges
  igraph.walk.weights <- graph_from_data_frame(edges, directed=F)
  
  if (!is.null(coords)) {
    cells.remain <- unique(c(edges$V1, edges$V2))
    coords.remain <- rownames(coords) %in% cells.remain
    coords <- coords[coords.remain,]
  }
  
  if (is.null(start.temp)) start.temp <- sqrt(vcount(igraph.walk.weights))
  
  # Do force-directed layout
  if (verbose) print(paste0(Sys.time(), ": Doing force-directed layout."))
  if (method=="fr") igraph.walk.layout <- layout_with_fr(igraph.walk.weights, dim=dim, coords=coords, start.temp=start.temp)
  if (method=="drl") igraph.walk.layout <- layout_with_drl(igraph.walk.weights, dim=dim, options=list(edge.cut=0))
  if (method=="kk") igraph.walk.layout <- layout_with_kk(igraph.walk.weights, dim=dim, coords=coords)
  
  # Store the layout
  if (dim==2) {
    # Calculate weighted pseudotime of cells' neighbors for tree layout.
    if (verbose) print(paste0(Sys.time(), ": Calculating Z."))
    neighbor.pt <- unlist(lapply(1:dim(walk.nn$nn.label)[1], function(i) {
      weighted.mean(x=object@pseudotime[walk.nn$nn.label[i,],pseudotime], w=1/walk.nn$nn.dists[i,])
    }))
    names(neighbor.pt) <- rownames(walk.nn$nn.label)
    # Store the layout
    object@tree$walks.force.layout <- as.data.frame(igraph.walk.layout, stringsAsFactors=F)
    names(object@tree$walks.force.layout) <- c("x","y")
    rownames(object@tree$walks.force.layout) <- V(igraph.walk.weights)$name
    # Normalize neighbor.pt range
    neighbor.pt.factor <- mean(apply(object@tree$walks.force.layout[,c("x","y")], 2, max)) / max(neighbor.pt)
    object@tree$walks.force.layout$telescope.pt <- neighbor.pt[rownames(object@tree$walks.force.layout)] * neighbor.pt.factor
    
    # Calculate local density for adjusting alpha
    if (verbose) print(paste0(Sys.time(), ": Calculating neighbor distance."))
    object <- fdl.density(object, neighbor=density.neighbors)
    
  } else if (dim==3) {
    object@tree$walks.force.layout.3d <- as.data.frame(igraph.walk.layout, stringsAsFactors=F)
    names(object@tree$walks.force.layout.3d) <- c("x","y","z")
    rownames(object@tree$walks.force.layout.3d) <- V(igraph.walk.weights)$name
  }
  
  object@tree$walks.force.labels <- treeForcePositionLabels(object)
  
  # Return the object
  if (verbose) print(paste0(Sys.time(), ": Finished."))
  return(object)
}

#' Plot tree force-directed layout in 2D
#' 
#' @import ggplot2
#' 
#' @param object An URD object
#' 
#' @export
plotTreeForce2D <- function(object, label=NULL, label.type="search", title=label, show.points=T, point.alpha=1, point.size=1, show.neighbors=F, neighbors.max=10000, colors=NULL) {
  # Load colors if not specified
  if (is.null(colors)) colors <- defaultURDContinuousColors()
  
  gg.data <- object@tree$walks.force.layout
  if (!is.null(label)) color.data <- data.for.plot(object = object, label = label, label.type = label.type, as.color = F, as.discrete.list=T, cells.use = rownames(gg.data))
  gg.data$expression <- color.data$data
  gg <- ggplot() + theme_bw() + labs(x="", y="", color="", title=title)
  if (show.neighbors) {
    edge.data <- object@tree$walks.force.edges
    if (dim(edge.data)[1] > neighbors.max) edge.data <- edge.data[sample(1:dim(edge.data)[1], neighbors.max, replace=F),]
    edge.data[,c("x1","y1")] <- gg.data[edge.data$V1,c("x","y")]
    edge.data[,c("x2","y2")] <- gg.data[edge.data$V2,c("x","y")]
    if (show.points) gg <- gg + geom_segment(data=edge.data, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=0.02)
    if (!show.points) gg <- gg + geom_segment(data=edge.data, aes(x=x1, y=y1, xend=x2, yend=y2, color=weight), alpha=0.04)
  }
  if (show.points) gg <- gg + geom_point(data=gg.data, aes(x=x, y=y, color=expression), size=point.size, alpha=point.alpha) 
  if (!color.data$discrete) gg <- gg + scale_color_gradientn(colors = colors)
  return(gg)
}

#' Plot force-directed layout of tree
#' 
#' This plots data on the force-directed layout of an URD tree, gerenated by \code{\link{treeForceDirectedLayout}}.
#' 
#' @param object An URD object
#' @param label (Character) Data to use for coloring tree (see \link{data.for.plot})
#' @param label.type (Character) Where to find data for label (default \code{"search"} auto-detects, see \link{data.for.plot})
#' @param view (Character) Name of view to use (i.e. the name of an entry in \code{@@tree$force.view.list}) created by \link{plotTreeForceStore3DView}. Default is last view created.
#' @param alpha (Numeric) Maximum transparency of points
#' @param alpha.fade (Numeric) Minimum transparency of points
#' @param size (Numeric) Size of points in the plot
#' @param title (Character) Title to add to the plot. (This is sensitive to resizing the window after plotting, but if a view is stored, the window will be resized before the title is added, so it will be acceptable resolution for figures.)
#' @param title.cex (Numeric) Adjust the title font size
#' @param title.line (Numeric) Adjust the position of the title. Positive numbers move the title upward.
#' @param label.tips (Logical) Should text identifying the tips of the tree be placed in the force directed layout? Defaults to TRUE if \code{@@tree$segment.names} has been set.
#' @param use.short.names (Logical) Should short names from \code{@@tree$segment.names.short} be used? Defaults to TRUE if those values have been set.
#' @param seg.show (Character vector) Segments of the tree to put on the plot (Default \code{NULL} is all segments)
#' @param cells.show (Character vector) Cells of the tree to show (Default \code{NULL} is all cells)
#' @param fade.below (Numeric) If desired, transparency can be lowered for points with low expression values (thereby highlighting the portions of the tree with expression). This value is the portion of the range of expression values to fade. (The default, 2/9, fades the bottom two-ninths of the expression values.) A fully 'faded' point will have alpha of \code{alpha.fade}.
#' @param density.alpha (Logical) Should points with short distance to nearest neighbors (i.e. in denser regions of the plot) be more transparent?
#' @param label.spacing (Numeric) How far should labels be from the detected tips
#' @param text.cex (Numeric) Size of the label text
#' @param colors (Character vector) Vector of colors to use if plotting continuous data
#' @param discrete.colors (Character vector) Vector of colors to use if plotting discrete data
#' 
#' @return Nothing. Produces a plot using the \code{rgl} package, displayed in an X11 window.
#' 
#' @export
plotTreeForce <- function(object, label, label.type="search", view="default", alpha=0.8, alpha.fade=0.1, size=5, title=NULL, title.cex=3, title.line=0, label.tips=(!is.null(object@tree$segment.names) | !is.null(object@tree$segment.names.short)), use.short.names=!is.null(object@tree$segment.names.short), seg.show=NULL, cells.show=NULL, fade.below=(2/9), density.alpha=T, label.spacing=5, text.cex=0.8, colors=NULL, discrete.colors=NULL) {
  
  if (requireNamespace("rgl", quietly = TRUE)) {
    
    # Figure out which view to use
    if (!is.null(view)) {
      if (view=="default") {
        if ("force.view.default" %in% names(object@tree)) view <- object@tree$force.view.default
        else view <- NULL
      }
      view <- object@tree$force.view.list[[view]]
    }
    
    # Get default colors
    if (is.null(colors)) colors <- defaultURDContinuousColors()
    
    # Get layout data and expression data
    gg.data <- object@tree$walks.force.layout
    gg.data$segment <- object@diff.data[rownames(gg.data),"segment"]
    if (!is.null(label)) {
      color.data <- data.for.plot(object = object, label = label, label.type = label.type, as.color = T, as.discrete.list=T, cells.use = rownames(gg.data), continuous.colors=colors, colors.use = discrete.colors)
      gg.data$expression <- color.data$data
    }
    gg.data$alpha <- alpha
    gg.data$size <- size
    
    # Focus on expression
    # For this, grey out any thing below focus, set alpha low for those, and increase progressively until 2*focus.
    if (fade.below > 0 & !color.data$discrete) {
      # Get the actual expression data
      expression.data <- data.for.plot(object=object, label=label, label.type=label.type, as.color=F, cells.use=rownames(gg.data))
      # Figure out ranges
      expression.range <- range(expression.data)
      fade.range <- diff(expression.range) * fade.below + expression.range[1]
      # Adjust gg.data
      fade <- (fade.range-expression.data+expression.range[1])/fade.range
      fade[fade < 0] <- 0
      gg.data$alpha <- alpha - ((alpha-alpha.fade) * fade)
    }
    
    # Open 3D view
    if (!is.null(view)) {
      rgl::open3d(
        zoom=view$rgl.setting$zoom, 
        scale=view$rgl.setting$scale, 
        userMatrix=view$rgl.setting$userMatrix, 
        windowRect=view$rgl.setting$windowRect
      )
    } else {
      rgl::open3d()
    }
    
    # Adjust transparency for local density
    if (density.alpha) {
      dr <- range(object@tree$walks.force.layout$n.dist)
      n.dist <- object@tree$walks.force.layout$n.dist / dr[2]
      density.adj <- sqrt(n.dist)
      density.adj <- density.adj / median(density.adj)
      gg.data$alpha <- gg.data$alpha * density.adj
    }
    
    # Convert alpha to 0 on cells not to show (so the plot size doesn't get changed)
    if (!is.null(seg.show) & is.null(cells.show)) {
      segs.show <- segChildrenAll(object, seg.show, include.self=T)
      cells.show <- rownames(gg.data)[which(gg.data$segment %in% segs.show)]
    }
    if (!is.null(cells.show)) {
      cells.hide <- setdiff(rownames(gg.data), cells.show)
      gg.data[cells.hide,"alpha"] <- 0
    }
    
    # Add the points to the plot
    rgl::points3d(x=gg.data$x, y=gg.data$y, z=gg.data$telescope.pt, col=gg.data$expression, alpha=gg.data$alpha, size=size)
    
    # Add title to plot
    if (!is.null(title)) {
      Sys.sleep(0.1)
      rgl::bgplot3d({plot.new(); title(main=title, line=title.line, cex.main=title.cex)})
      Sys.sleep(0.1)
    }
    
    # Add labels to tips
    if (label.tips) {
      
      # Check for existing tree labels
      if (!is.null(object@tree$walks.force.labels)) {
        tip.labels <- object@tree$walks.force.labels
      } else {
        # If there aren't any, generate some.
        tip.labels <- treeForcePositionLabels(object)
      }
      
      # Choose short or long labels
      if (use.short.names) tip.labels$label <- tip.labels$name.short else tip.labels$label <- tip.labels$name
      
      # Limit to labels of cells that appear in the plot
      segs.to.label <- intersect(unique(gg.data[which(gg.data$alpha > 0),"segment"]), tip.labels$seg)
      tip.labels <- tip.labels[which(tip.labels$seg %in% segs.to.label),]
      
      # Add the text
      rgl::text3d(x=tip.labels$x, y=tip.labels$y, z=tip.labels$z, text=tip.labels$label, adj=0.5, cex=text.cex)
    }
    
    # Add a brief pause to ensure that rendering completes before moving on to next function.
    Sys.sleep(0.1)
  } else {
    stop("Package rgl is required for this function. To install: install.packages('rgl')\n")
  }
}

#' Plot force-directed layout of tree
#' 
#' This plots data on the force-directed layout of an URD tree, gerenated by \code{\link{treeForceDirectedLayout}}.
#' 
#' @param object An URD object
#' @param label.red (Character) Data to use for coloring tree (red channel) (see \link{data.for.plot})
#' @param label.green (Character) Data to use for coloring tree (green channel) (see \link{data.for.plot})
#' @param label.type.red (Character) Where to find data for red channel label (default \code{"search"} auto-detects, see \link{data.for.plot})
#' @param label.type.green (Character) Where to find data for green channel label (default \code{"search"} auto-detects, see \link{data.for.plot})
#' @param view (Character) Name of view to use (i.e. the name of an entry in \code{@@tree$force.view.list}) created by \link{plotTreeForceStore3DView}. Default is last view created.
#' @param alpha (Numeric) Maximum transparency of points
#' @param alpha.fade (Numeric) Minimum transparency of points
#' @param size (Numeric) Size of points in the plot
#' @param title (Character) Title to add to the plot. (This is sensitive to resizing the window after plotting, but if a view is stored, the window will be resized before the title is added, so it will be acceptable resolution for figures.)
#' @param title.cex (Numeric) Adjust the title font size
#' @param title.line (Numeric) Adjust the position of the title. Positive numbers move the title upward.
#' @param label.tips (Logical) Should text identifying the tips of the tree be placed in the force directed layout? Defaults to TRUE if \code{@@tree$segment.names} has been set.
#' @param use.short.names (Logical) Should short names from \code{@@tree$segment.names.short} be used? Defaults to TRUE if those values have been set.
#' @param seg.show (Character vector) Segments of the tree to put on the plot (Default \code{NULL} is all segments)
#' @param cells.show (Character vector) Cells of the tree to show (Default \code{NULL} is all cells)
#' @param fade.below (Numeric) If desired, transparency can be lowered for points with low expression values (thereby highlighting the portions of the tree with expression). This value is the portion of the range of expression values to fade. (The default, 1/3, fades the bottom third of the expression values.) A fully 'faded' point will have alpha of \code{alpha.fade}.
#' @param fade.method (Character: "max", "min", or "mean") Should fade be calculated on max expression of the two labels (emphasize any cells expressing either label), min expression of the two labels (emphasize only cells co-expressing both labels), or mean expression
#' @param density.alpha (Logical) Should points with short distance to nearest neighbors (i.e. in denser regions of the plot) be more transparent?
#' @param label.spacing (Numeric) How far should labels be from the detected tips
#' @param text.cex (Numeric) Size of the label text
#' 
#' @return Nothing. Produces a plot using the \code{rgl} package, displayed in an X11 window.
#' 
#' @export
plotTreeForceDual <- function(object, label.red, label.green, label.red.type="search", label.green.type="search", view="default", alpha=0.8, alpha.fade=0.025, size=5, title=NULL, title.cex=3, title.line=0, label.tips=(!is.null(object@tree$segment.names) | !is.null(object@tree$segment.names.short)), use.short.names=!is.null(object@tree$segment.names.short), seg.show=NULL, cells.show=NULL, fade.below=(1/3), fade.method=c("max", "min", "mean"), density.alpha=T, label.spacing=5, text.cex=0.8) {
  
  if (requireNamespace("rgl", quietly = TRUE)) {
    
    # Figure out which view to use
    if (!is.null(view)) {
      if (view=="default") {
        if ("force.view.default" %in% names(object@tree)) view <- object@tree$force.view.default
        else view <- NULL
      }
      view <- object@tree$force.view.list[[view]]
    }
    
    # Get layout data
    gg.data <- object@tree$walks.force.layout
    gg.data$segment <- object@diff.data[rownames(gg.data),"segment"]
    
    # Get expression data
    red.data <- data.for.plot(object=object, label=label.red, label.type=label.red.type, as.single.color = T, as.discrete.list=T, cells.use=rownames(gg.data))
    green.data <- data.for.plot(object=object, label=label.green, label.type=label.green.type, as.single.color = T, as.discrete.list=T, cells.use=rownames(gg.data))
    if (red.data$discrete | green.data$discrete) stop("For dual-color plots, both labels must be continuous values.")
    gg.data$expression <- rgb(red.data$data, green.data$data, 0)
    gg.data$alpha <- alpha
    gg.data$size <- size
    
    # Focus on expression
    # For this, grey out any thing below focus, set alpha low for those, and increase progressively until 2*focus.
    if (fade.below > 0) {
      # Figure out value for fading
      if (length(fade.method) > 1) fade.method <- fade.method[1]
      if (tolower(fade.method) == "max") {
        # Use max of red or green for fading
        expression.data <- pmax(red.data$data, green.data$data)
      } else if (tolower(fade.method) == "min") {
        # Use min of red or green for fading
        expression.data <- pmin(red.data$data, green.data$data)
      } else if (tolower(fade.method) == "mean") {
        # Use min of red or green for fading
        expression.data <- apply(data.frame(red.data$data, green.data$data), 1, mean)
      } else {
        stop("fade.method must be max, min, or mean.")
      }
      # Adjust gg.data
      fade <- (fade.below - expression.data)/fade.below
      fade[fade < 0] <- 0
      gg.data$alpha <- alpha - ((alpha-alpha.fade) * fade)
    }
    
    # Open 3D view
    if (!is.null(view)) {
      rgl::open3d(
        zoom=view$rgl.setting$zoom, 
        scale=view$rgl.setting$scale, 
        userMatrix=view$rgl.setting$userMatrix, 
        windowRect=view$rgl.setting$windowRect
      )
    } else {
      rgl::open3d()
    }
    
    # Adjust transparency for local density
    if (density.alpha) {
      dr <- range(object@tree$walks.force.layout$n.dist)
      n.dist <- object@tree$walks.force.layout$n.dist / dr[2]
      density.adj <- sqrt(n.dist)
      density.adj <- density.adj / median(density.adj)
      gg.data$alpha <- gg.data$alpha * density.adj
    }
    
    # Convert alpha to 0 on cells not to show (so the plot size doesn't get changed)
    if (!is.null(seg.show) & is.null(cells.show)) {
      segs.show <- segChildrenAll(object, seg.show, include.self=T)
      cells.show <- rownames(gg.data)[which(gg.data$segment %in% segs.show)]
    }
    if (!is.null(cells.show)) {
      cells.hide <- setdiff(rownames(gg.data), cells.show)
      gg.data[cells.hide,"alpha"] <- 0
    }
    
    # Add the points to the plot
    rgl::points3d(x=gg.data$x, y=gg.data$y, z=gg.data$telescope.pt, col=gg.data$expression, alpha=gg.data$alpha, size=size)
    
    # Add title to plot
    if (!is.null(title)) {
      Sys.sleep(0.1)
      rgl::bgplot3d({plot.new(); title(main=title, line=title.line, cex.main=title.cex)})
      Sys.sleep(0.1)
    }
    
    # Add labels to tips
    if (label.tips) {
      
      # Check for existing tree labels
      if (!is.null(object@tree$walks.force.labels)) {
        tip.labels <- object@tree$walks.force.labels
      } else {
        # If there aren't any, generate some.
        tip.labels <- treeForcePositionLabels(object)
      }
      
      # Choose short or long labels
      if (use.short.names) tip.labels$label <- tip.labels$name.short else tip.labels$label <- tip.labels$name
      
      # Limit to labels of cells that appear in the plot
      segs.to.label <- intersect(unique(gg.data[which(gg.data$alpha > 0),"segment"]), tip.labels$seg)
      tip.labels <- tip.labels[which(tip.labels$seg %in% segs.to.label),]
      
      # Add the text
      rgl::text3d(x=tip.labels$x, y=tip.labels$y, z=tip.labels$z, text=tip.labels$label, adj=0.5, cex=text.cex)
    }
    
    # Add a brief pause to ensure that rendering completes before moving on to next function.
    Sys.sleep(0.1)
  } else {
    stop("Package rgl is required for this function. To install: install.packages('rgl')\n")
  }
}

#' Determine position for segment labels in force-directed layout
#' 
#' Finds cells near the end of each tip, and attempts to calculate a vector from
#' them to determine a good location for a label near each tip. Called automatically
#' by \code{\link{treeForceDirectedLayout}}, but can be re-run if needed.
#' 
#' @param object An URD object
#' @param label.spacing (Numeric) Length of vector from final cell in each tip to label
#' 
#' @return A data.frame with rows as tips, and columns containing labels and coordinates.
#' This is normally stored in \code{@@tree$walks.force.labels}.
#' 
#' @examples 
#' object@tree$walks.force.labels <- treeForcePositionLabels(object, label.spacing=2) # Move labels closer
#' 
#' @export
treeForcePositionLabels <- function(object, label.spacing=5) {
  
  # Grab the tips and their names
  tip.labels <- data.frame(
    seg=segTerminal(object),
    name=object@tree$segment.names[segTerminal(object)],
    name.short=object@tree$segment.names.short[segTerminal(object)],
    stringsAsFactors=F, row.names=segTerminal(object)
  )
  
  # Anything that wasn't assigned a name previously: NA -> blank
  tip.labels[is.na(tip.labels)] <- ""
  
  # Get force-directed layout and add segment and pseudotime information
  fdl.data <- object@tree$walks.force.layout
  fdl.data$segment <- object@diff.data[rownames(fdl.data),"segment"]
  fdl.data$pseudotime <- object@tree$pseudotime[rownames(fdl.data)]
  
  # Order fdl data by pseudotime
  fdl.data <- fdl.data[order(fdl.data$pseudotime, decreasing=T),]
  
  # Figure out the location for each label by projecting out the vector at the end of each tip
  tip.labels[,c("x","y","z")] <- t(as.data.frame(lapply(tip.labels$seg, function(tip) {
    # Get cells from this tip
    this.fdl <- fdl.data[which(fdl.data$segment==tip),]
    # Find the vector describing direction that the cells at the tip are moving
    tip.vec <- as.numeric(apply(this.fdl[c(1,5),c("x","y","telescope.pt")], 2, diff))
    # Normalize the length of vector
    s <- sqrt(label.spacing^2/(tip.vec[1]^2+tip.vec[2]^2+tip.vec[3]^2))
    tip.vec <- s*tip.vec
    # Place the label
    label.coords <- apply(rbind(this.fdl[1,c("x","y","telescope.pt")], tip.vec),2,sum)
    return(label.coords)
  })))
  
  return(tip.labels)
}

#' Calculate force-directed layout local density
#' 
#' Add distance to nth nearest neighbor to force directed layout so that you can later
#' calculate local density and adjust alpha values.
#' 
#' @importFrom RANN nn2
#' 
#' @param object An URD object
#' @param neighbor (Numeric) Distance to \code{neighbor}th nearest neighbor is used for
#' density adjustment of point transparency on force-directed layout plots.
fdl.density <- function(object, neighbor=10) {
  nn <- RANN::nn2(object@tree$walks.force.layout, k = neighbor, treetype = 'bd')
  object@tree$walks.force.layout$n.dist <- nn$nn.dists[,neighbor]
  return(object)
}

#' Store a 3D view for force-directed layouts
#' 
#' The last stored view will be the "default."
#' 
#' @param object An URD object
#' @param view.name (Character) Name to use for this view
#' 
#' @return An URD object with the parameters for that view stored in \code{@@tree$force.view.list}.
#' 
#' @export
plotTreeForceStore3DView <- function(object, view.name) {
  if (requireNamespace("rgl", quietly = TRUE)) {
    
    # If not appending to a views list, make one
    if (is.null(object@tree$force.view.list)) object@tree$force.view.list <- list()
    # Grab rgl plot settings, if desired
    rgl.settings <- rgl::par3d(no.readonly=F)
    # Make a list of everything
    this.view <- list(rgl.settings=rgl.settings)
    # Store if in the views list
    object@tree$force.view.list[[view.name]] <- this.view
    # Store 'default' view
    object@tree$force.view.default <- view.name
    # Return the final views list
    return(object)
  } else {
    stop("Package rgl is required for this function. To install: install.packages('rgl')\n")
  }
}

#' Rotate coordinates in 3D
#' 
#' @param m (Matrix) n cells x 3 matrix of coordinates
#' @param a (Numeric) Angle to rotate
#' @param axis (Character) Which axis to rotate around
#' 
#' @return n x 3 matrix of rotated coordinates
rotateCoords3d <- function(m, a, axis=c("x","y","z")) {
  if (length(axis) > 1) axis <- axis[1]
  if (dim(m)[2] != 3) stop("Supply an n x 3 matrix.")
  if (axis == "x") r.mat <- matrix(c(1, 0, 0, 0, cos(a), -1*sin(a), 0, sin(a), cos(a)), ncol=3, byrow=T)
  if (axis == "y") r.mat <- matrix(c(cos(a), 0, sin(a), 0, 1, 0, -1*sin(a), 0, cos(a)), ncol=3, byrow=T)
  if (axis == "z") r.mat <- matrix(c(cos(a), -1*sin(a), 0, sin(a), cos(a), 0, 0, 0, 1), ncol=3, byrow=T)
  m %*% r.mat
}

#' Translate coordinates in 3D
#' 
#' @param m (Matrix) n cells x 3 coordinates matrix
#' @param x (Numeric) Distance to move in x
#' @param y (Numeric) Distance to move in y
#' @param z (Numeric) Distance to move in z
#' 
#' @return n x 3 matrix of translated coordinates
translateCoords3d <- function(m, x, y, z) {
  if (dim(m)[2] != 3) stop("Supply an n x 3 matrix.")
  sweep(m, 2, c(x, y, z), "+")
}

#' Rotate points in force-directed layout
#' 
#' This function can be used in addition to \code{\link{treeForceTranslateCoords}}
#' in order to fine-tune the presentation of a force-directed layout. For instance,
#' this can be used to achieve an improved 2D visualization by modifying overlapping
#' portions of the layout.
#' 
#' @param object An URD object
#' @param cells (Character vector) Cells to modify (Default \code{NULL} is all cells in the force-directed layout)
#' @param seg (Character) Instead of specifying cells, just grab all cells from this segment and downstream. Ignored if \code{cells} is specified.
#' @param angle (Numeric) Angle to rotate
#' @param axis (Character: "x", "y", "z") Axis around which to rotate
#' @param around.cell (Character or Numeric) If a character, then the name of a cell to rotate around. If numeric, then will choose the Nth cell by pseudotime. (This can be helpful for rotating a segment around its base to spread segments of the tree apart more. If \code{NULL}, then rotates around (0,0,0).
#' @param throw.out.cells (Numeric) Throw out the youngest N cells in pseudotime
#' @param pseudotime (Character) Pseudotime (i.e. column name of \code{@@pseudotime}) to use for determining cells to throw out.
#' 
#' @return An URD object with the coordinates of some cells in \code{@@tree$walks.force.layout} modified.
#' 
#' @export
treeForceRotateCoords <- function(object, cells=NULL, seg=NULL, angle, axis=c("x", "y", "z"), around.cell=NULL, throw.out.cells=0, pseudotime=NULL) {
  # If segment specified, grab the cells from that segment
  if (!is.null(seg) & is.null(cells)) cells <- cellsInCluster(object, "segment", segChildrenAll(object, seg, include.self=T))
  # If cells not specified, do all cells in the layout.
  if (is.null(cells)) cells <- rownames(object@tree$walks.force.layout)
  # Make sure all cells specified are actually in the force-directed layout
  cells <- intersect(cells, rownames(object@tree$walks.force.layout))
  # If throw.out.cells is on, by removing n cells with the youngest pseudotime.
  if ((throw.out.cells > 0)) {
    if (is.null(pseudotime)) stop ("To use throw.out.cells, must specify pseudotime.")
    cells <- setdiff(cells, cells[order(object@pseudotime[cells,pseudotime])[1:throw.out.cells]])
  }
  # Grab the appropriate chunk of the layout.
  m <- as.matrix(object@tree$walks.force.layout[cells, c("x","y","telescope.pt")])
  # If a cell is specified to rotate around, translate coords so that cell is at 0, 0, 0
  if (!is.null(around.cell)) {
    if (is.numeric(around.cell)) {
      #around.cell <- rownames(m)[order(m[,3])[around.cell]]
      around.cell <- cells[order(object@pseudotime[cells,pseudotime])[around.cell]]
    }
    around.coords <- m[around.cell,]
    m <- translateCoords3d(m=m, x=-1*around.coords[1], y=-1*around.coords[2], z=-1*around.coords[3])
  }
  # Do the rotation
  mprime <- rotateCoords3d(m=m, a=angle, axis=axis)
  # If you were rotating around an arbitrary cell, put them back
  if (!is.null(around.cell)) {
    mprime <- translateCoords3d(m=mprime, x=around.coords[1], y=around.coords[2], z=around.coords[3])
  }
  object@tree$walks.force.layout[cells, c("x", "y", "telescope.pt")] <- mprime
  return(object)
}

#' Translate points in force-directed layout
#' 
#' This function can be used in addition to \code{\link{treeForceRotateCoords}}
#' in order to fine-tune the presentation of a force-directed layout. For instance,
#' this can be used to achieve an improved 2D visualization by modifying overlapping
#' portions of the layout.
#' 
#' @param object An URD object
#' @param cells (Character vector) Cells to modify (Default \code{NULL} is all cells in the force-directed layout)
#' @param seg (Character) Instead of specifying cells, just grab all cells from this segment and downstream. Ignored if \code{cells} is specified.
#' @param x (Numeric) Distance to move cells along x-axis
#' @param y (Numeric) Distance to move cells along y-axis
#' @param z (Numeric) Distance to move cells along z-axis
#' 
#' @return An URD object with the coordinates of some cells in \code{@@tree$walks.force.layout} modified.
#' 
#' @export
treeForceTranslateCoords <- function(object, cells=NULL, seg=NULL, x, y, z) {
  if (!is.null(seg) & is.null(cells)) cells <- cellsInCluster(object, "segment", segChildrenAll(object, seg, include.self=T))
  if (is.null(cells)) cells <- rownames(object@tree$walks.force.layout)
  cells <- intersect(cells, rownames(object@tree$walks.force.layout))
  m <- as.matrix(object@tree$walks.force.layout[cells, c("x","y","telescope.pt")])
  mprime <- translateCoords3d(m=m, x=x, y=y, z=z)
  object@tree$walks.force.layout[cells, c("x", "y", "telescope.pt")] <- mprime
  return(object)
}