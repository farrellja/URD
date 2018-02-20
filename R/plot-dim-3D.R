
plotDim3D <- function(object, reduction.use=c("dm", "pca"), view.use=NULL, dim.1=NULL, dim.2=NULL, dim.3=NULL, w.1=NULL, w.2=NULL, w.3=NULL, color.by=NULL, discrete=F, cells.plot=NULL, alpha=0.1, title=NULL, axis.labels=T, add=FALSE, colors.to.use=NULL, ...) {
  
  # Get data to plot
  if (length(reduction.use) > 1) reduction.use <- reduction.use[1]
  if (tolower(reduction.use)=="pca") {
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
  
  # Configure parameters
  if (is.null(view.use) & any(is.null(dim.1), is.null(dim.2), is.null(dim.3))) stop("Specify either a view (view.use) or dimensions (dim.1, dim.2, dim.3)")
  if (is.null(dim.1)) dim.1 <- view.use$dim.1
  if (is.null(dim.2)) dim.2 <- view.use$dim.2
  if (is.null(dim.3)) dim.3 <- view.use$dim.3
  if (is.null(w.1)) {
    if (!is.null(view.use)) { w.1 <- view.use$w.1 }
    else { w.1 <- rep(1/length(dim.1),length(dim.1)) }
  }
  if (is.null(w.2)) {
    if (!is.null(view.use)) { w.2 <- view.use$w.2 }
    else { w.2 <- rep(1/length(dim.2),length(dim.2)) }
  }
  if (is.null(w.3)) {
    if (!is.null(view.use)) { w.3 <- view.use$w.3 }
    else { w.3 <- rep(1/length(dim.3),length(dim.3)) }
  }
  
  # Plot all cells if not specified
  if(is.null(cells.plot)) cells.plot <- rownames(diff.data)
  
  # Get the data out, interpolation of dimensions if desired
  data.plot <- interpolate.points(diff.data, "DC", dim.1, dim.2, dim.3, w.1, w.2, w.3, cells.plot)
  
  # Make a color scale
  if (is.null(colors.to.use)) {
    if (is.null(color.by)) {
      colors.to.use <- 'black'
    } else {
      if (discrete) {
        score.to.use <- as.factor(diff.data[cells.plot,color.by])
        palette.to.use <- rainbow(length(levels(score.to.use)))
        colors.to.use <- palette.to.use[score.to.use]
      } else {
        colors.to.use <- cscale(diff.data[cells.plot,color.by], gradient_n_pal(colrs))
      }
    }
  }
  
  # Figure out axis labels
  if (axis.labels) {
    xlab <- paste("DC", dim.1, sep="", collapse="+")
    ylab <- paste("DC", dim.2, sep="", collapse="+")
    zlab <- paste("DC", dim.3, sep="", collapse="+")
  } else {
    xlab <- ""
    ylab <- ""
    zlab <- ""
  }
  
  # Make a base plot
  if (!add) {
    if (!is.null(view.use$rgl.settings)) {
      open3d(zoom=view.use$rgl.settings$zoom, scale=view.use$rgl.settings$scale, userMatrix=view.use$rgl.settings$userMatrix, windowRect=view.use$rgl.settings$windowRect)
      plot3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], xlab=xlab, ylab=ylab, zlab=zlab, alpha=alpha, col=colors.to.use, xlim=view.use$rgl.settings$bbox[1:2], ylim=view.use$rgl.settings$bbox[3:4], zlim=view.use$rgl.settings$bbox[5:6], ...)
      #open3d(params=view.use$rgl.settings)
    } else {
      open3d()
      plot3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], xlab=xlab, ylab=ylab, zlab=zlab, alpha=alpha, col=colors.to.use, ...)
    }
    
  } else {
    points3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], alpha=alpha, col=colors.to.use)
  }
  
  # Add title if desired
  if (!is.null(title)) {
    bgplot3d({
      plot.new()
      title(main=title, line=3)
    })
  }
}

# Function to normalize vectors
vector.norm <- function(x) sqrt(sum(x^2))

# Function to interpolate a set of points for plotting
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

# Store a 3D view for diffusion map renderings
store.3d.view <- function(views.list=NULL, dim.1=1, dim.2=2, dim.3=3, w.1=rep(1/length(dim.1),length(dim.1)), w.2=rep(1/length(dim.2),length(dim.2)), w.3=rep(1/length(dim.3),length(dim.3)), view.name=NULL, store.rgl.rotation=T) {
  # If not appending to a views list, make one
  if (is.null(views.list)) views.list <- list()
  # Make a default name if one isn't provided
  if(is.null(view.name)) view.name <- as.character(length(views.list) + 1)
  # Grab rgl plot settings, if desired
  if (store.rgl.rotation) {
    rgl.settings <- par3d(no.readonly=F)
  } else {
    rgl.settings <- NULL
  }
  # Make a list of everything
  this.view <- list(dim.1=dim.1, dim.2=dim.2, dim.3=dim.3, w.1=w.1, w.2=w.2, w.3=w.3, rgl.settings=rgl.settings)
  # Store if in the views list
  views.list[[view.name]] <- this.view
  # Return the final views list
  return(views.list)
}

# Plot cells in 3D according to DCs
plot.dm.3d <- function(object, label=NULL, view.use=NULL, dim.1=NULL, dim.2=NULL, dim.3=NULL, w.1=NULL, w.2=NULL, w.3=NULL, label.type="search", cells.plot=NULL, alpha=0.3, title=NULL, axis.labels=T, legend=T, add=FALSE, label.2=NULL, label.3=NULL, colors.use=NULL, ...) {
  
  if (is.null(label) && !is.null(label.2)) stop("If providing label.2, must provide label.")
  
  # Configure parameters
  if (!is.null(view.use) && (is.null(object@tree$view.list) || !(view.use %in% names(object@tree$view.list)))) stop(view.use, " is not in view.list. Use store.view\n")
  if (is.null(view.use) & any(is.null(dim.1), is.null(dim.2), is.null(dim.3))) stop("Specify either a view (view.use) or dimensions (dim.1, dim.2, dim.3)")
  if (is.null(dim.1)) dim.1 <- object@tree$view.list[[view.use]]$dim.1
  if (is.null(dim.2)) dim.2 <- object@tree$view.list[[view.use]]$dim.2
  if (is.null(dim.3)) dim.3 <- object@tree$view.list[[view.use]]$dim.3
  if (is.null(w.1)) {
    if (!is.null(view.use)) { w.1 <- object@tree$view.list[[view.use]]$w.1 }
    else { w.1 <- rep(1/length(dim.1),length(dim.1)) }
  }
  if (is.null(w.2)) {
    if (!is.null(view.use)) { w.2 <- object@tree$view.list[[view.use]]$w.2 }
    else { w.2 <- rep(1/length(dim.2),length(dim.2)) }
  }
  if (is.null(w.3)) {
    if (!is.null(view.use)) { w.3 <- object@tree$view.list[[view.use]]$w.3 }
    else { w.3 <- rep(1/length(dim.3),length(dim.3)) }
  }
  
  # Plot all cells if not specified
  if(is.null(cells.plot)) cells.plot <- rownames(object@dm@eigenvectors)
  
  # Get the data out, interpolation of dimensions if desired
  data.plot <- interpolate.points(object@dm@eigenvectors, "DC", dim.1, dim.2, dim.3, w.1, w.2, w.3, cells.plot)
  
  # Get them colors right
  if (!is.null(label)) {
    if (is.null(label.2)) {
      color.data <- data.for.plot(object = object, name = label, type = label.type, cells.use = cells.plot, as.color = T, as.discrete.list = T, colors.use=colors.use)
      data.plot$color <- color.data$data
    } else {
      color.data <- data.for.plot(object = object, name = label, type = label.type, cells.use = cells.plot, as.single.color = T, as.discrete.list = T)
      color.data.2 <- data.for.plot(object = object, name = label.2, type = label.type, cells.use = cells.plot, as.single.color = T, as.discrete.list = T)
      if (!is.null(label.3)) {
        color.data.3 <- data.for.plot(object = object, name = label.3, type = label.type, cells.use = cells.plot, as.single.color = T, as.discrete.list = T)
      } else {
        color.data.3 <- list(data=0)
      }
      data.plot$color <- rgb(color.data$data, color.data.2$data, color.data.3$data)
    }
  } else {
    data.plot$color <- 'black'
  }
  
  # Figure out axis labels
  if (axis.labels) {
    xlab <- paste("DC", dim.1, sep="", collapse="+")
    ylab <- paste("DC", dim.2, sep="", collapse="+")
    zlab <- paste("DC", dim.3, sep="", collapse="+")
    bbox=T
  } else {
    xlab <- ""
    ylab <- ""
    zlab <- ""
    bbox=F
  }
  
  # Make a base plot
  if (!add) {
    if (!is.null(view.use) && !is.null(object@tree$view.list[[view.use]]$rgl.settings)) {
      open3d(
        zoom=object@tree$view.list[[view.use]]$rgl.setting$zoom, 
        scale=object@tree$view.list[[view.use]]$rgl.setting$scale, 
        userMatrix=object@tree$view.list[[view.use]]$rgl.setting$userMatrix, 
        windowRect=object@tree$view.list[[view.use]]$rgl.setting$windowRect
      )
      plot3d(
        x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], 
        xlab=xlab, ylab=ylab, zlab=zlab, 
        alpha=alpha, col=data.plot[,"color"], 
        xlim=object@tree$view.list[[view.use]]$rgl.setting$bbox[1:2], 
        ylim=object@tree$view.list[[view.use]]$rgl.setting$bbox[3:4], 
        zlim=object@tree$view.list[[view.use]]$rgl.setting$bbox[5:6],
        axes=bbox, ...
      )
    } else {
      open3d()
      plot3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], col=data.plot[,"color"], xlab=xlab, ylab=ylab, zlab=zlab, alpha=alpha, axes=bbox, ...)
    }
  } else {
    points3d(x=data.plot[,"d1"], y=data.plot[,"d2"], z=data.plot[,"d3"], col=data.plot[,"color"], alpha=alpha)
  }
  
  # Add title if desired
  if (!is.null(title)) {
    bgplot3d({
      plot.new()
      title(main=title, line=3)
    })
  }
  
  # Add legend is desired
  if (legend && !is.null(label) && color.data$discrete) {
    legend3d("topright", legend = names(color.data$legend), pch=16, col=color.data$legend, cex=1, inset=c(0.02))
  }
}

# Store a 3D view for diffusion map renderings
store.3d.view <- function(object, view.name, dim.1=1, dim.2=2, dim.3=3, w.1=NULL, w.2=NULL, w.3=NULL, store.rgl.rotation=T) {
  # If not appending to a views list, make one
  if (is.null(object@tree$view.list)) object@tree$view.list <- list()
  # Generate weights if they're not provided
  if (is.null(w.1)) w.1 <- rep(1/length(dim.1),length(dim.1))
  if (is.null(w.2)) w.2 <- rep(1/length(dim.2),length(dim.2))
  if (is.null(w.3)) w.3 <- rep(1/length(dim.3),length(dim.3))
  # Grab rgl plot settings, if desired
  if (store.rgl.rotation) {
    rgl.settings <- par3d(no.readonly=F)
  } else {
    rgl.settings <- NULL
  }
  # Make a list of everything
  this.view <- list(dim.1=dim.1, dim.2=dim.2, dim.3=dim.3, w.1=w.1, w.2=w.2, w.3=w.3, rgl.settings=rgl.settings)
  # Store if in the views list
  object@tree$view.list[[view.name]] <- this.view
  # Return the final views list
  return(object)
}