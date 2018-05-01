#' Dimensionality reduction plot in 3D.
#' 
#' Plots cells according to their coordinates in a dimensionality reduction (tSNE by default,
#' but also PCA or diffusion map). Cells are colored according to a user set \code{label} that
#' can range from gene expression, to metadata values, or cluster identity. See 
#' \code{\link{data.for.plot}} for more information about labels that can be chosen. Views can
#' be stored using \code{\link{plotDim3DStoreView}} that can reproduce the same orientation
#' plot repeatedly.
#' 
#' @param object An URD object
#' @param label (Character) Data to use for coloring points (e.g. a metadata name, group ID from clustering, or a gene name)
#' @param label.type (Character) Type of data to search for the label. Default is "search" which checks several data types in order. For more information: \code{\link{data.for.plot}}
#' @param reduction.use (Character) Dimensionality reduction to use (Diffusion Map ("dm") or PCA ("pca"))
#' @param view (Character) Stored view to use. These are created with \code{\link{plotDim3DStoreView}} and contain the dimensions and weights to use, as well as the stored orientation of the plot. Should be the name of an entry in the list \code{@@plot.3d}.
#' @param dim.1 (Numeric vector) Component(s) to use on x-axis. If > 1, they will be interpolated according to weights in \code{w.1}.
#' @param dim.2 (Numeric vector) Component(s) to use on y-axis. If > 1, they will be interpolated according to weights in \code{w.2}.
#' @param dim.3 (Numeric vector) Component(s) to use on z-axis. If > 1, they will be interpolated according to weights in \code{w.3}.
#' @param w.1 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param w.2 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param w.3 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param cells (Character vector) Cells to plot. (Default \code{NULL} plots all cells.)
#' @param alpha (Numeric) Transparency of points on plot: 0 (Transparent) - 1 (Opaque)
#' @param size (Numeric) Size of points on plot
#' @param title (Character) Title for plot. (This is sensitive to resizing the window after plotting, but if a view is stored, the window will be resized before the title is added, so it will be acceptable resolution for figures.)
#' @param title.cex (Numeric) Adjust the title font size
#' @param title.line (Numeric) Adjust the position of the title. Positive numbers move the title upward.
#' @param bounding.box (Logical) Should a bounding box with axis labels be displayed on the plot?
#' @param xlab (Character) Label for x-axis (if \code{bounding.box=T})
#' @param ylab (Character) Label for y-axis (if \code{bounding.box=T})
#' @param zlab (Character) Label for z-axis (if \code{bounding.box=T})
#' @param continuous.colors (Character vector) Vector of colors to use if data is continuous
#' @param continuous.color.limits (Numeric vector, length 2) Limits of data for color scale. Data outside these limits will be squished. \code{NULL} uses the original range of the data.
#' @param discrete.colors (Character vector) Vector of colors to use if data is discrete
#'
#' @return Returns nothing, but produces an output through rgl.
#' 
#' @export
plotDim3D <- function(object, label, label.type="search", reduction.use=c("dm", "pca"), view=NULL, dim.1=NULL, dim.2=NULL, dim.3=NULL, w.1=NULL, w.2=NULL, w.3=NULL, cells=NULL, alpha=0.2, size=4, title=NULL, title.cex=3, title.line=0, bounding.box=T, xlab=NULL, ylab=NULL, zlab=NULL, continuous.colors=NULL, continuous.color.limits=NULL, discrete.colors=NULL) {
  
  if (requireNamespace("rgl", quietly = TRUE)) {
  
    # Extract the view from the object
    if (!is.null(view)) {
      view <- object@plot.3d[[view]]
      reduction.use <- view$reduction.use
    }
    
    # Configure interpolation parameters
    if (is.null(view) & any(is.null(dim.1), is.null(dim.2), is.null(dim.3))) stop("Specify either a view or dimensions to plot.")
    if (is.null(dim.1)) dim.1 <- view$dim.1
    if (is.null(dim.2)) dim.2 <- view$dim.2
    if (is.null(dim.3)) dim.3 <- view$dim.3
    if (is.null(w.1)) {
      if (!is.null(view)) { w.1 <- view$w.1 }
      else { w.1 <- rep(1/length(dim.1),length(dim.1)) }
    }
    if (is.null(w.2)) {
      if (!is.null(view)) { w.2 <- view$w.2 }
      else { w.2 <- rep(1/length(dim.2),length(dim.2)) }
    }
    if (is.null(w.3)) {
      if (!is.null(view)) { w.3 <- view$w.3 }
      else { w.3 <- rep(1/length(dim.3),length(dim.3)) }
    }
    
    # Get the position data
    if (length(reduction.use) > 1) reduction.use <- reduction.use[1]
    if (tolower(reduction.use)=="pca") {
      data.plot <- object@pca.scores
      if (any(c(dim.1, dim.2, dim.3) > ncol(object@pca.scores))) stop("Dimensions requested were not previously calculated.")
      dim.1 <- paste0("PC", dim.1)
      dim.2 <- paste0("PC", dim.2)
      dim.3 <- paste0("PC", dim.3)
      data.plot <- object@pca.scores[,c(dim.1, dim.2, dim.3)]
    } else if (tolower(reduction.use)=="dm") {
      data.plot <- object@dm@eigenvectors
      if (any(c(dim.1, dim.2, dim.3) > ncol(object@dm@eigenvectors))) stop("Dimensions requested were not previously calculated.")
      dim.1 <- paste0("DC", dim.1)
      dim.2 <- paste0("DC", dim.2)
      dim.3 <- paste0("DC", dim.3)
      data.plot <- as.data.frame(object@dm@eigenvectors[,c(dim.1, dim.2, dim.3)])
    } else {
      stop("The reduction.use provided is invalid.")
    }
    
    # Determine cells to plot
    if(is.null(cells)) cells <- rownames(data.plot)
    
    # Interpolate the position data if desired
    data.plot <- interpolate.points(raw.coordinates = data.plot, coord.prefix = "", dim.1 = dim.1, dim.2 = dim.2, dim.3 = dim.3, w.1 = w.1, w.2 = w.2, w.3 = w.3, rows.subset = cells)
    rownames(data.plot) <- cells
    
    # Get colors, if no default is provided.
    if (is.null(continuous.colors)) continuous.colors <- defaultURDContinuousColors()
    
    # Get label data
    if (!is.null(label)) {
      color.data <- data.for.plot(object = object, label = label, label.type = label.type, as.color = T, as.discrete.list=T, cells.use = cells, continuous.colors=continuous.colors, colors.use = discrete.colors, continuous.color.limits = continuous.color.limits)
      data.plot$color <- color.data$data
    } else {
      data.plot$color <- "#000000"
    }
    
    # Figure out axis labels
    if (is.null(xlab)) xlab <- paste(dim.1, sep="", collapse="+")
    if (is.null(ylab)) ylab <- paste(dim.2, sep="", collapse="+")
    if (is.null(zlab)) zlab <- paste(dim.3, sep="", collapse="+")
    
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
    
    # Add the plot
    if (bounding.box) {
      rgl::plot3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], col=data.plot[,"color"], xlab=xlab, ylab=ylab, zlab=zlab, alpha=alpha, size=size, axes=T)
    } else {
      rgl::plot3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], col=data.plot[,"color"], xlab="", ylab="", zlab="", alpha=alpha, size=size, axes=F)
    }
    
    # Add title if desired
    if (!is.null(title)) {
      # Brief pause to make sure window has resized 
      Sys.sleep(0.2)
      rgl::bgplot3d({plot.new(); title(main=title, line=title.line, cex.main=title.cex)})
    }
    # Brief pause to give function time to finish before moving on
    Sys.sleep(0.2)
  } else {
    stop("Package rgl is required for this function. To install: install.packages('rgl')\n")
  }
}

#' Interpolate dimensionality reduction dimensions
#' 
#' Allows for fitting more than three dimensions into a dimensionality reduction plot.
#' 
#' @param raw.coordinates (data.frame) Rows are cells, and columns are dimensions to use for interpolation.
#' @param coord.prefix (Character) Append this prefix to the dimensions
#' @param dim.1 (Character vector) Dimensions to interpolate as dimension 1
#' @param dim.2 (Character vector) Dimensions to interpolate as dimension 2
#' @param dim.3 (Character vector) Dimensions to interpolate as dimension 3
#' @param w.1 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param w.2 (Numeric vector) Weight to use for the members of \code{dim.2}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param w.3 (Numeric vector) Weight to use for the members of \code{dim.3}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param rows.subset (Character vector) Which rows to subset before interpolating. (Default \code{NULL} uses all rows.)
#' 
#' @return nrows x 3 data.frame, where > 3 dimensions have been interpolated to produce a 3-dimensional embedding.
#' 
#' @keywords internal
interpolate.points <- function(raw.coordinates, coord.prefix="DC", dim.1, dim.2, dim.3, w.1=rep(1/length(dim.1),length(dim.1)), w.2=rep(1/length(dim.2),length(dim.2)), w.3=rep(1/length(dim.3),length(dim.3)), rows.subset=NULL) {
  # If subset of rows to do is not provided, then do all rows
  if (is.null(rows.subset)) rows.subset <- rownames(raw.coordinates)
  # Reformat dimensions
  dim.1 <- paste0(coord.prefix, dim.1)
  dim.2 <- paste0(coord.prefix, dim.2)
  dim.3 <- paste0(coord.prefix, dim.3)
  # Interpolate data into x, y, z data.frame
  interpolated.data <- data.frame(
    d1=apply(as.data.frame(raw.coordinates[rows.subset,dim.1]), 1, function(x) sum(x*w.1)),
    d2=apply(as.data.frame(raw.coordinates[rows.subset,dim.2]), 1, function(x) sum(x*w.2)),
    d3=apply(as.data.frame(raw.coordinates[rows.subset,dim.3]), 1, function(x) sum(x*w.3))
  )
  return(interpolated.data)
}

#' Store a 3D view for plotDim3D
#' 
#' This allows storing "views" which contain the information about which dimensionality reduction,
#' dimensions, weights, and 3D orientation to use to make plots. This is useful for
#' generating figures with many similar panels.
#' 
#' @param object An URD object
#' @param view.name (Character) What the name of the view should be.
#' @param reduction.use (Character) Dimensionality reduction to use (Diffusion Map ("dm") or PCA ("pca"))
#' @param dim.1 (Numeric vector) Component(s) to use on x-axis. If > 1, they will be interpolated according to weights in \code{w.1}.
#' @param dim.2 (Numeric vector) Component(s) to use on y-axis. If > 1, they will be interpolated according to weights in \code{w.2}.
#' @param dim.3 (Numeric vector) Component(s) to use on z-axis. If > 1, they will be interpolated according to weights in \code{w.3}.
#' @param w.1 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param w.2 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' @param w.3 (Numeric vector) Weight to use for the members of \code{dim.1}. Defaults to equal weights for all members. Important to note that these can be negative to reverse the contribution of particular dimensions.
#' 
#' @return An URD object with a new 3D view stored in \code{@@plot.3d}.
#' 
#' @export
plotDim3DStoreView <- function(object, view.name, reduction.use=c("dm", "pca"), dim.1=1, dim.2=2, dim.3=3, w.1=NULL, w.2=NULL, w.3=NULL) {
  if (requireNamespace("rgl", quietly = TRUE)) {
  
    # Check reduction.use parameter
    if (length(reduction.use) > 1) reduction.use <- reduction.use[1]
    reduction.use <- tolower(reduction.use)
    if (!(reduction.use %in% c("dm", "pca"))) stop("reduction.use must be 'dm' or 'pca'")
    
    # Generate weights if they're not provided
    if (is.null(w.1)) w.1 <- rep(1/length(dim.1),length(dim.1))
    if (is.null(w.2)) w.2 <- rep(1/length(dim.2),length(dim.2))
    if (is.null(w.3)) w.3 <- rep(1/length(dim.3),length(dim.3))
    
    # Grab rgl plot settings, if desired
    rgl.settings <- rgl::par3d(no.readonly=F)
  
    # Make a list of everything
    this.view <- list(reduction.use=reduction.use, dim.1=dim.1, dim.2=dim.2, dim.3=dim.3, w.1=w.1, w.2=w.2, w.3=w.3, rgl.settings=rgl.settings)
    
    # Store if in the views list and return object
    object@plot.3d[[view.name]] <- this.view
    return(object)
  } else {
    stop("Package rgl is required for this function. To install: install.packages('rgl')\n")
  }
}