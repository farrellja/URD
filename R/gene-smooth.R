
#' Gene Smooth: Fit expression of genes
#' 
#' This takes a group of genes and cells, averages gene expression in groups of cells (determined using a moving window through pseudotime), and then uses smoothing algorithms (e.g. LOESS or spline fitting) to describe the expression of each gene. The results are returned for use in \code{\link{plotSmoothFit}} to visualize the fits or several sets of results (e.g. from different trajectories) can be combined with \code{\link{cropSmoothFit}} and plotted with \code{\link{plotSmoothFitMultiCascade}}.
#'
#' @importFrom stats lowess smooth.spline
#'
#' @param object An URD object
#' @param pseudotime (Character) Name of pseudotime (i.e. a column name of \code{@@pseudotime})
#' @param cells (Character vector) Cells to include
#' @param genes (Character vector) Genes to include
#' @param method (Character) Which smoothing method to use ("lowess" calls \code{\link{stats::lowess}} and "spline" calls \code{\link{stats::smooth.spline}})
#' @param moving.window (Numeric) Number of bins to use per window
#' @param cells.per.window (Numeric or \code{NULL}) Size of bins (number of cells)
#' @param pseudotime.per.window (Numeric or \code{NULL}) Size of bins (pseudotime)
#' @param verbose (Logical) Display progress reports?
#' @param verbose.genes (Logical) Display a message when starting each gene?
#' 
#' @return (Named List):
#' \itemize{ 
#'   \item{\strong{\code{pt.windows}}: (List) Character vectors of cell IDs of cells in each window, named by either the mean, min, or max pseudotime of those cells (depending on \code{name.by})}
#'   \item{\strong{\code{mean.expression}}: (data.frame) Mean expression of each gene in each pseudotime window}
#'   \item{\strong{\code{scaled.expression}}: (data.frame) Mean expression of each gene in each pseudotime window, scaled to the maximum observed expression}
#'   \item{\strong{\code{mean.smooth}}: (data.frame) Value of fit curve for each pseudotime window (i.e. to predict value of mean.expression)}
#'   \item{\strong{\code{scaled.smooth}}: (data.frame) Value of fit curve for each pseudotime window, scaled to max expression (i.e. to predict value of scaled.expression)}
#'   \item{\strong{\code{mean.expression.red}}: (data.frame) Mean expression reduced to the same dimensions of \code{mean.smooth}. (When method is 'spline' if two windows have the same pseudotime, only a single prediction is made for that pseudotime, so this will average the two input windows.)}
#'   \item{\strong{\code{scaled.expression.red}}: (data.frame) Scaled expression reduced to the same dimensions of \code{scaled.smooth}.}
#'   \item{\strong{\code{method}}: (Character) Identify fit method ("impulse")}
#' }
#' 
#' @examples 
#' # Fit a gene expression curve using "spline" method for
#' # genes "endo.genes" and cells in segments 47 and 46 in the tree.
#' tentacle.spline <- geneSmoothFit(urd.object, method = "spline", 
#'   pseudotime = "pseudotime", 
#'   cells = cellsInCluster(urd.object, "segment", c("47", "46")), 
#'   genes = endo.genes, moving.window = 1, cells.per.window = 5, spar = 0.9
#' )
#' 
#' @export
## TO DO: Special case for considering single-cells only that skips apply step, and will run faster.
geneSmoothFit <- function(object, pseudotime, cells, genes, method=c("lowess", "spline"), moving.window=3, cells.per.window=NULL, pseudotime.per.window=NULL, verbose=T, verbose.genes=F, ...) {
  
  # Check inputs 
  method <- tolower(method[1])
  if (!(method %in% c("lowess", "spline"))) stop("Method must be 'lowess' or 'spline'")
  
  # Get moving window of cells by pseudotime
  if (verbose) print(paste0(Sys.time(), ": Calculating moving window expression."))
  pt.windows <- pseudotimeMovingWindow(object, pseudotime=pseudotime, cells=cells, moving.window=moving.window, cells.per.window=cells.per.window, pseudotime.per.window = pseudotime.per.window)
  # Calculate pseudotime parameters for each window
  pt.info <- as.data.frame(t(as.data.frame(lapply(pt.windows, function(cells) {
    pt <- object@pseudotime[cells,pseudotime]
    return(c(mean(pt), min(pt), max(pt), diff(range(pt))))
  }))))
  names(pt.info) <- c("mean","min","max","width")
  rownames(pt.info) <- 1:length(pt.windows)
  
  # Make aggregated expression data and scale it
  # (Note special case: if doing single cells can just copy expression data)
  if (moving.window == 1 && !is.null(cells.per.window) && cells.per.window == 1) {
    mean.expression <- as.data.frame(as.matrix(object@logupx.data[genes,cells]))
  } else {
    mean.expression <- as.data.frame(lapply(pt.windows, function(window.cells) apply(object@logupx.data[genes,window.cells,drop=F], 1, mean.of.logs)))
  }
  names(mean.expression) <- names(pt.windows)
  scaled.expression <- sweep(mean.expression, 1, apply(mean.expression, 1, max), "/")
  scaled.expression[which(apply(mean.expression, 1, max) == 0),] <- 0
  
  # Lowess/spline fit for unscaled data
  if (verbose) print(paste0(Sys.time(), ": Generating un-scaled fits."))
  x <- as.numeric(colnames(mean.expression))
  if (method == "lowess") {
    mean.smooth <- as.data.frame(t(as.data.frame(lapply(1:nrow(mean.expression), function(i) {
      if (verbose.genes) print(paste0(i, ": ", rownames(mean.expression)[i]))
      lf <- stats::lowess(x=x, y=mean.expression[i,], ...)
      return(lf$y)
    }))))
  } else if (method == "spline") {
    mean.smooth <- as.data.frame(t(as.data.frame(lapply(1:nrow(mean.expression), function(i) {
      if (verbose.genes) print(paste0(i, ": ", rownames(mean.expression)[i]))
      lf <- stats::smooth.spline(x=x, y=mean.expression[i,], ...)
      return(lf$y)
    }))))
  }
  rownames(mean.smooth) <- rownames(mean.expression)
  if (method == "lowess") colnames(mean.smooth) <- x
  if (method == "spline") colnames(mean.smooth) <- unique(x)
  
  # Lowess/spline fit for scaled data
  if (verbose) print(paste0(Sys.time(), ": Generating scaled fits."))
  if (method=="lowess") {
    scaled.smooth <- as.data.frame(t(as.data.frame(lapply(1:nrow(scaled.expression), function(i) {
      if (verbose.genes) print(paste0(i, ": ", rownames(scaled.expression)[i]))
      lf <- stats::lowess(x=x, y=scaled.expression[i,], ...)
      return(lf$y)
    }))))
  } else if (method=="spline") {
    scaled.smooth <- as.data.frame(t(as.data.frame(lapply(1:nrow(scaled.expression), function(i) {
      if (verbose.genes) print(paste0(i, ": ", rownames(scaled.expression)[i]))  
      lf <- stats::smooth.spline(x=x, y=scaled.expression[i,], ...)
      return(lf$y)
    }))))
  }
  rownames(scaled.smooth) <- rownames(scaled.expression)
  if (method == "lowess") colnames(scaled.smooth) <- x
  if (method == "spline") colnames(scaled.smooth) <- unique(x)
  
  # Generate list
  smoothed.fit <- list(
    pt.windows=pt.windows,
    mean.expression=mean.expression,
    scaled.expression=scaled.expression,
    mean.smooth=mean.smooth,
    scaled.smooth=scaled.smooth,
    method=method
  )
  
  # Dimension-matched expression data to go with spline fits
  if (method=="spline") {
    if (verbose) print(paste0(Sys.time(), ": Reducing mean expression data to same dimensions as spline fits."))
    smoothed.fit <- geneSmoothReduce(smoothed.fit)
  } else {
    smoothed.fit$mean.expression.red <- smoothed.fit$mean.expression
    smoothed.fit$scaled.expression.red <- smoothed.fit$scaled.expression
  }
  
  return(smoothed.fit)
}

#' Gene Smooth: Plot a single fit
#' 
#' This plots gene expression in groups of cells and a curve representing its mean expression (generated using a smoothing algorithm). It takes output from either \code{\link{geneSmoothFit}} or \code{\link{geneCascadeProcess}}. (Output from \code{\link{geneCascadeProcess}} may also be plotted with \code{\link{geneCascadeImpulsePlots}}.) If multiple genes are provided, they are plotted on the same plot in different colors if \code{multiplot=F} or as a grid of plots per gene if \code{multiplot=T}.
#' 
#' @param smoothed.fit (List) Output from either \code{\link{geneSmoothFit}} or \code{\link{geneCascadeProcess}}
#' @param genes (Character vector) Genes to include in the plot
#' @param scaled (Logical) Plot actual expression values (\code{FALSE}) or expression scaled to its maximum value (\code{TRUE})
#' @param multiplot (Logical) If multiple genes are provided, they are plotted on the same plot in different colors if \code{multiplot=F} or as a grid of plots per gene if \code{multiplot=T}
#' @param plot.data (Logical) Plot data points?
#' @param alpha.data (Numeric: 0-1) Transparency of data points
#' @param alpha.smooth (Numeric: 0-1) Transparency of curve
#' @param lwd.smooth (Numeric) Line width of curve
#' @param plot.title (Character) Title of plot (Default \code{NULL} is no title)
#' 
#' @return A ggplot2 object
#' 
#' @examples
#' # Calculate spline fit
#' expressed.spline <- geneSmoothFit(hydra.spumous, method = "spline", 
#'   pseudotime = "pseudotime", cells = colnames(hydra.spumous@logupx.data),
#'   genes = expressed.genes, moving.window = 1, cells.per.window = 5, 
#'   spar = 0.875
#' )
#' # Plot a few genes on the fit curves
#' genes.spline.plot <- c("t14194aep|WNT3_MOUSE", "t15597aep|WNT1_DANRE",
#'   "t20768aep|BRAC_CANLF", "t29725aep|BRAC_CHICK", "t22116aep|ETV1_MOUSE"
#' )
#' plotSmoothFit(expressed.spline, genes = genes.spline.plot, scaled = T)
#' 
#' @export
plotSmoothFit <- function(smoothed.fit, genes, scaled=T, multiplot=F, plot.data=T, alpha.data=0.2, alpha.smooth=1, lwd.smooth=1, plot.title=NULL) {
  
  # Validate parameters
  if (class(smoothed.fit) != "list") {
    if (class(smoothed.fit) == "URD") stop("plotSmoothFit requires a spline object; you have provided an URD object.") else stop("plotSmoothFit should be provided output from geneSmoothFit, which is class 'list'.")
  }
  if (!("mean.expression" %in% names(smoothed.fit))) {
    if ("mean.expression" %in% names(smoothed.fit[[1]])) {
      stop("You have provided a multi-segment spline, which should be plotted with the function plotSmoothFitMultiCascade.")
    } else {
      stop("The list provided does not seem to be output from the geneSmoothFit function.")
    }
  }
  
  # Grab scaled or un-scaled data
  if (scaled) {
    expression <- smoothed.fit$scaled.expression[genes,,drop=F]
    smooth <- smoothed.fit$scaled.smooth[genes,,drop=F]
  } else {
    expression <- smoothed.fit$mean.expression[genes,,drop=F]
    smooth <- smoothed.fit$mean.smooth[genes,,drop=F]
  }
  
  # Melt for ggplot
  expression$Gene <- rownames(expression)
  smooth$Gene <- rownames(smooth)
  expression.melt <- reshape2::melt(expression, id.vars="Gene", variable.name="Pseudotime", value.name="Expression")
  smooth.melt <- reshape2::melt(smooth, id.vars="Gene", variable.name="Pseudotime", value.name="Expression")
  expression.melt$Pseudotime <- as.numeric(as.character(expression.melt$Pseudotime))
  smooth.melt$Pseudotime <- as.numeric(as.character(smooth.melt$Pseudotime))
  
  if (multiplot && length(genes) > 1) color.by="Gene" else color.by="black"
  
  the.plot <- ggplot(data=smooth.melt, aes_string(x="Pseudotime", y="Expression", group="Gene", color="Gene")) + theme_bw()
  if (plot.data) the.plot <- the.plot + geom_point(data=expression.melt, alpha=alpha.data)
  the.plot <- the.plot + geom_line(alpha=alpha.smooth, lwd=lwd.smooth)
  
  
  # Add title
  if (is.null(plot.title) && length(genes) == 1) plot.title <- genes
  if (!is.null(plot.title)) the.plot <- the.plot + ggtitle(plot.title)
  
  # If multiplotting
  if (multiplot && length(genes) > 1) {
    # Turn on faceting
    the.plot <- the.plot + facet_wrap("~Gene")
  }
  if (multiplot || length(genes) == 1) {
    # Kill the legend
    the.plot <- the.plot + guides(color=F)
    # Turn colors off
    the.plot <- the.plot + scale_color_manual(values = rep('black', length(genes)))
  }
  
  return(the.plot)
}

#' Gene Smooth: Crop to particular pseudotimes
#' 
#' This crops a result from \code{\link{geneSmoothFit}} to a particular pseudotime range. Can be used to split a smoothed result into segments in order to label them separately or combine them for use with \code{\link{plotSmoothFitMultiCascade}}.
#' 
#' @param smoothed.fit (List) Result from \code{\link{geneSmoothFit}}
#' @param pt.min (Numeric) Minimum pseudotime to retain
#' @param pt.max (Numeric) Maximum pseudotime to retain
#' 
#' @return (List) Result from \code{\link{geneSmoothFit}}, but cropped to the range of \code{pt.min} and \code{pt.max}.
#' 
#' @examples 
#' # Identify the pseudotime of the branchpoint
#' pt.crop <- as.numeric(unlist(hydra.en.tree@tree$segment.pseudotime.limits)[1])
#' # Crop according to the pseudotime of the branchpoint
#' foot.only.spline <- cropSmoothFit(tentacle.spline, pt.max = pt.crop)
#' hypostome.only.spline <- cropSmoothFit(hypo.spline, pt.min = pt.crop)
#' tentacle.only.spline <- cropSmoothFit(tentacle.spline, pt.min = pt.crop)
#' # Combine into a list
#' # Names in the plots are determined by names of the smooth objects in the list
#' splines <- list(foot.only.spline, tentacle.only.spline, hypostome.only.spline)
#' names(splines) <- c("Foot/Body", "Tentacle", "Hypostome")
#' # Plot expression of genes on the curves
#' endoderm.genes.plot <- c("t14194aep|WNT3_MOUSE", "t20768aep|BRAC_CANLF", 
#'   "t25396aep|NKX26_HUMAN", "t18735aep|FOXA2_ORYLA")
#' plotSmoothFitMultiCascade(smoothed.fits = splines, 
#'   genes = endoderm.genes.plot, scaled = F, 
#'   colors = c(`Foot/Body` = "#FF8C00", Hypostome = "#1E90FF", 
#'     Tentacle = "#32CD32"
#'   )
#' )
#' 
#' @export
cropSmoothFit <- function(smoothed.fit, pt.min=-Inf, pt.max=Inf) {
  # Figure out which windows fit the pseudotime limits in the raw data
  pt <- as.numeric(names(smoothed.fit$pt.windows))
  pt.in <- which(pt >= pt.min & pt <= pt.max)
  
  # Figure out which windows fit the pseudotime limits in the smoothed data
  smoothed.pt <- as.numeric(colnames(smoothed.fit$scaled.smooth))
  smoothed.pt.in <- which(smoothed.pt >= pt.min & smoothed.pt <= pt.max)
  
  # Crop things to the dimensions that you want
  smoothed.fit$pt.windows <- smoothed.fit$pt.windows[pt.in]
  smoothed.fit$mean.expression <- smoothed.fit$mean.expression[,pt.in]
  smoothed.fit$scaled.expression <- smoothed.fit$scaled.expression[,pt.in]
  smoothed.fit$mean.smooth <- smoothed.fit$mean.smooth[,smoothed.pt.in]
  smoothed.fit$scaled.smooth <- smoothed.fit$scaled.smooth[,smoothed.pt.in]
  smoothed.fit$scaled.expression.red <- smoothed.fit$scaled.expression.red[,smoothed.pt.in]
  smoothed.fit$mean.expression.red <- smoothed.fit$mean.expression.red[,smoothed.pt.in]
  
  # Rename columns since R silently renames them!!!
  colnames(smoothed.fit$mean.expression) <- as.character(pt)[pt.in]
  colnames(smoothed.fit$scaled.expression) <- as.character(pt)[pt.in]
  colnames(smoothed.fit$mean.smooth) <- as.character(smoothed.pt)[smoothed.pt.in]
  colnames(smoothed.fit$scaled.smooth) <- as.character(smoothed.pt)[smoothed.pt.in]
  colnames(smoothed.fit$scaled.expression.red) <- as.character(smoothed.pt)[smoothed.pt.in]
  colnames(smoothed.fit$mean.expression.red) <- as.character(smoothed.pt)[smoothed.pt.in]
  
  return(smoothed.fit)
}

#' Gene Smooth: Combine smoothed fits
#' 
#' This function takes a list of smoothed.fits (i.e. output from \code{\link{geneSmoothFit}}) and combines them in the order provided into a single smoothed fit. This is useful, for instance, for combining multiple pieces of fits together for plotting on a heatmap with several panels.
#' 
#' @param fit.list (List) A list of smoothed.fits to combine
#' @return (List) Akin to a single output from \code{\link{geneSmoothFit}}.
#' @export
#' 
combineSmoothFit <- function(fit.list) {
  # Grab pseudotime names
  pt.names.full <- unlist(lapply(fit.list, function(x) colnames(x$mean.expression)))
  pt.names.red <- unlist(lapply(fit.list, function(x) colnames(x$mean.smooth)))
  
  # Combine the many elements of each entry in the list
  out.mean.expression <- do.call("cbind", lapply(fit.list, function(x) x$mean.expression))
  colnames(out.mean.expression) <- pt.names.full
  out.mean.smooth <- do.call("cbind", lapply(fit.list, function(x) x$mean.smooth))
  colnames(out.mean.smooth) <- pt.names.red
  out.mean.expression.red <- do.call("cbind", lapply(fit.list, function(x) x$mean.expression.red))
  colnames(out.mean.expression.red) <- pt.names.red
  out.method <- unique(unlist(lapply(fit.list, function(x) x$method)))
  if (length(out.method) > 1) warning("Smoothing method of all inputs was not the same.")
  out.pt.windows <- unlist(lapply(fit.list, function(x) x$pt.windows), recursive = F)
  names(out.pt.windows) <- pt.names.full
  
  # Re-scale scaled values according to the max of any present 
  out.scaled.expression <- sweep(out.mean.expression, 1, apply(out.mean.expression, 1, max), "/")
  out.scaled.smooth <- sweep(out.mean.smooth, 1, apply(out.mean.smooth, 1, max), "/")
  out.scaled.expression.red <- sweep(out.mean.expression.red, 1, apply(out.mean.expression.red, 1, max), "/")
  
  # Return a new list
  return(list(
    pt.windows=out.pt.windows,
    mean.expression=out.mean.expression,
    scaled.expression=out.scaled.expression,
    mean.smooth=out.mean.smooth,
    scaled.smooth=out.scaled.smooth,
    method=out.method,
    scaled.expression.red=out.scaled.expression.red,
    mean.expressoin.red=out.mean.expression.red
  ))
}

#' Gene Smooth: Plot multiple fits
#' 
#' This plots gene expression in groups of cells and a curve representing its mean expression (generated using a smoothing algorithm). It takes output from multiple runs of \code{\link{geneSmoothFit}} (i.e. for different trajectories or different segments of the same trajectory). The individual results from \code{\link{geneSmoothFit}} should be combined into a list, and the names of the elements will be used to define the names of the segments in the plot. If multiple genes are provided, a grid of plots is produced.
#' 
#' @param smoothed.fit (List of lists) A list containing several lists that are the output of \code{\link{geneSmoothFit}}. Each one will generate a separate segment in the plot, labeled according to its name in the list.
#' @param genes (Character vector) Genes to include in the plot(s)
#' @param colors (Character vector) Colors to use for each segment (names are matched against the names of elements in \code{smoothed.fit}).
#' @param scaled (Logical) Plot actual expression values (\code{FALSE}) or expression scaled to its maximum value (\code{TRUE})
#' @param plot.data (Logical) Plot data points?
#' @param alpha.data (Numeric: 0-1) Transparency of data points
#' @param alpha.smooth (Numeric: 0-1) Transparency of curve
#' @param lwd.smooth (Numeric) Line width of curve
#' @param ncol (Numeric) Number of columns to use if several plots are produced (i.e. if several \code{genes} are provided). (Default \code{NULL} will attempt to produce a square aspect ratio.)
#' 
#' @return A ggplot2 object
#' 
#' @examples 
#' # Identify the pseudotime of the branchpoint
#' pt.crop <- as.numeric(unlist(hydra.en.tree@tree$segment.pseudotime.limits)[1])
#' # Crop according to the pseudotime of the branchpoint
#' foot.only.spline <- cropSmoothFit(tentacle.spline, pt.max = pt.crop)
#' hypostome.only.spline <- cropSmoothFit(hypo.spline, pt.min = pt.crop)
#' tentacle.only.spline <- cropSmoothFit(tentacle.spline, pt.min = pt.crop)
#' # Combine into a list
#' # Names in the plots are determined by names of the smooth objects in the list
#' splines <- list(foot.only.spline, tentacle.only.spline, hypostome.only.spline)
#' names(splines) <- c("Foot/Body", "Tentacle", "Hypostome")
#' # Plot expression of genes on the curves
#' endoderm.genes.plot <- c("t14194aep|WNT3_MOUSE", "t20768aep|BRAC_CANLF", 
#'   "t25396aep|NKX26_HUMAN", "t18735aep|FOXA2_ORYLA")
#' plotSmoothFitMultiCascade(smoothed.fits = splines, 
#'   genes = endoderm.genes.plot, scaled = F, 
#'   colors = c(`Foot/Body` = "#FF8C00", Hypostome = "#1E90FF", 
#'     Tentacle = "#32CD32"
#'   )
#' )
#' 
#' @export
plotSmoothFitMultiCascade <- function(smoothed.fits, genes, colors=NULL, scaled=T, plot.data=T, alpha.data=0.2, alpha.smooth=1, lwd.smooth=1, ncol=NULL) {
  
  # Validate parameters
  if (class(smoothed.fits) != "list") {
    if (class(smoothed.fits) == "URD") stop("plotSmoothFitMultiCascade requires a list of spline objects; you have provided an URD object.") else stop("plotSmoothFitMultiCascade should be provided a list of outputs from geneSmoothFit, with list names set to the names desired for fit segments.")
  }
  if (!("mean.expression" %in% names(smoothed.fits[[1]]))) {
    if ("mean.expression" %in% names(smoothed.fits)) {
      stop("You have provided the output of geneSmoothFit, which should be plotted with the function plotSmoothFit. plotSmoothFitMultiCascade should be provided a list of outputs from geneSmoothFit, with list names set to the names desired for fit segments.")
    } else {
      stop("smoothed.fits does not seem a list of splines output from the geneSmoothFit function.")
    }
  }
  
  # Grab scaled or un-scaled data
  if (scaled) {
    expression <- lapply(1:length(smoothed.fits), function(i) {
      x <- smoothed.fits[[i]]$scaled.expression[genes,,drop=F]
      x$Gene <- rownames(x)
      x$Type <- names(smoothed.fits)[i]
      return(x)
    })
    smooth <- lapply(1:length(smoothed.fits), function(i) {
      x <- smoothed.fits[[i]]$scaled.smooth[genes,,drop=F]
      x$Gene <- rownames(x)
      x$Type <- names(smoothed.fits)[i]
      return(x)
    })
  } else {
    expression <- lapply(1:length(smoothed.fits), function(i) {
      x <- smoothed.fits[[i]]$mean.expression[genes,,drop=F]
      x$Gene <- rownames(x)
      x$Type <- names(smoothed.fits)[i]
      return(x)
    })
    smooth <- lapply(1:length(smoothed.fits), function(i) {
      x <- smoothed.fits[[i]]$mean.smooth[genes,,drop=F]
      x$Gene <- rownames(x)
      x$Type <- names(smoothed.fits)[i]
      return(x)
    })
  }
  
  # Melt for ggplot
  expression.melt <- do.call("rbind", lapply(expression, reshape2::melt, id.vars=c("Gene", "Type"), variable.name="Pseudotime", value.name="Expression"))
  smooth.melt <- do.call("rbind", lapply(smooth, reshape2::melt, id.vars=c("Gene", "Type"), variable.name="Pseudotime", value.name="Expression"))
  expression.melt$Pseudotime <- as.numeric(as.character(expression.melt$Pseudotime))
  smooth.melt$Pseudotime <- as.numeric(as.character(smooth.melt$Pseudotime))
  
  the.plot <- ggplot(data=smooth.melt, aes(x=Pseudotime, y=Expression, group=Type, color=Type)) + theme_bw() + facet_wrap(~Gene, scales="free_y", ncol=ncol)
  if (scaled) the.plot <- the.plot + ylim(-0.1, 1.1)
  if (plot.data) the.plot <- the.plot + geom_point(data=expression.melt, alpha=alpha.data)
  the.plot <- the.plot + geom_line(alpha=alpha.smooth, lwd=lwd.smooth)
  if (!is.null(colors)) the.plot <- the.plot + scale_color_manual(values=colors)
  return(the.plot)
}

#' Gene Smooth: Reduce dimensions of Smooth Fit
#' 
#' When multiple windows have the same pseudotime, and \code{\link{geneSmoothFit}}
#' is run with \code{method="spline"}, a single output point will be created,
#' creating smoothed expression data with different dimensions from the original
#' expression data. This function produces 'reduced' data frames of the same
#' dimensions and loads them into the smoothed.fit object by aggregating windows
#' with the same pseudotime.
#' 
#' @param smoothed.fit (List) Output from \code{\link{geneSmoothFit}}
#' @param agg.func (function) Function to use to aggregate columns with same
#' pseudotime. Default is mean.
#' @return (List) \code{smoothed.fit} input with \code{$mean.expression.red} and \code{$scaled.expression.red} added.
#' 
#' @export
#' @keywords internal
geneSmoothReduce <- function(smoothed.fit, agg.func=mean) {
  # Get unique column names
  cols <- unique(colnames(smoothed.fit$scaled.expression))
  # Find how many times each column name occurs
  cols.n <- table(colnames(smoothed.fit$scaled.expression))
  # Generate a new data frame that's reduced for scaled data
  scaled.expression.reduced <- as.data.frame(lapply(cols, function(col) {
    # Find columns with that column name
    col.i <- which(colnames(smoothed.fit$scaled.expression) == col)
    # If more than one column with that name, aggregate them.
    # Else, just return the single column
    if (cols.n[col] > 1) {
      return(apply(smoothed.fit$scaled.expression[,col.i], 1, agg.func))
    } else {
      return(smoothed.fit$scaled.expression[,col.i])
    }
  }))
  colnames(scaled.expression.reduced) <- cols
  # Generate a new data frame that's reduced for mean data
  mean.expression.reduced <- as.data.frame(lapply(cols, function(col) {
    # Find columns with that column name
    col.i <- which(colnames(smoothed.fit$mean.expression) == col)
    # If more than one column with that name, aggregate them.
    # Else, just return the single column
    if (cols.n[col] > 1) {
      return(apply(smoothed.fit$mean.expression[,col.i], 1, agg.func))
    } else {
      return(smoothed.fit$mean.expression[,col.i])
    }
  }))
  colnames(mean.expression.reduced) <- cols
  # Add back into the smoothed.fit object and return.
  smoothed.fit$scaled.expression.red <- scaled.expression.reduced
  smoothed.fit$mean.expression.red <- mean.expression.reduced
  return(smoothed.fit)
}

#' Reduce matrix
#' 
#' This takes a matrix/data.frame and combines all columns with the same column name
#' (for instance, the same pseudotime) using \code{agg.fun}, which defaults
#' to \code{mean}.
#' 
#' @param m (data.frame or Matrix) Input
#' @param agg.func (function) Function to use to aggregate columns with same
#' colname. Default is \code{mean}.
#' @return (data.frame) Input with columns with same colname aggregated.
#' 
#' @export
#' @keywords internal
matrixReduce <- function(m, agg.func=mean) {
  # Get unique column names
  cols <- unique(colnames(m))
  # Find how many times each column name occurs
  cols.n <- table(colnames(m))
  # Generate a new data frame that's reduced for scaled data
  m.reduced <- as.data.frame(lapply(cols, function(col) {
    # Find columns with that column name
    col.i <- which(colnames(m) == col)
    # If more than one column with that name, aggregate them.
    # Else, just return the single column
    if (cols.n[col] > 1) {
      return(apply(m[,col.i], 1, agg.func))
    } else {
      return(m[,col.i])
    }
  }))
  colnames(m.reduced) <- cols
  return(m.reduced)
}