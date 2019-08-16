#' Moving window through pseudotime
#' 
#' Generates lists of cells using a moving window. Size of windows is either determined by pseudotime (if \code{pseudotime.per.window} is set) or by number of cells (if \code{cells.per.window} is set). If both are set, windows are determined by number of cells, but then windows whose pseudotime differ by less than \code{pseudotime.per.window} are collapsed.
#' 
#' @param object A URD object
#' @param pseudotime (Character) Name of column in \code{@@pseudotime} to use for pseudotime
#' @param cells (Character vector) Names of cells to include
#' @param moving.window (Numeric) Number of bins to use per window
#' @param cells.per.window (Numeric or \code{NULL}) Size of bins (number of cells)
#' @param pseudotime.per.window (Numeric or \code{NULL}) Size of bins (pseudotime)
#' @param name.by (Character: "mean", "min", or "max") 
#' 
#' @return List of windows of cells, named by either the mean, min, or max pseudotime of those cells (depending on \code{name.by}).
#' 
#' @keywords internal
#' @export

pseudotimeMovingWindow <- function(object, pseudotime, cells, moving.window, cells.per.window=NULL, pseudotime.per.window=NULL, name.by=c("mean","min","max")) {
  # Configure function to work based on mode
  # This is a hacky patch...
  if (!is.null(pseudotime.per.window) && is.null(cells.per.window)) {
    moving.pt.window <- T
  } else if (is.null(pseudotime.per.window) && !is.null(cells.per.window)) {
    moving.pt.window <- F
  } else stop("Define either cells.per.window or pseudotime.per.window (but not both).")
  
  # Figure out pseudotime of cells
  pt <- object@pseudotime[cells,pseudotime]
  pt.order <- order(pt)
  
  # Figure out window size
  if (moving.pt.window) {
    n.windows <- round((max(pt,na.rm=T)-min(pt,na.rm=T)) / pseudotime.per.window)
    pt.cut <- cut(pt, n.windows, labels = FALSE)
    c.windows <- lapply(1:n.windows, function(x) return(cells[which(pt.cut==x)]))
    
  } else {
    n.windows <- round(length(cells) / cells.per.window)
    c.windows.end <- round((1:n.windows) * (length(cells) / n.windows))
    c.windows.start <- c(0,head(c.windows.end, -1))+1
    c.windows <- lapply(1:n.windows, function(window) return(cells[pt.order[c.windows.start[window]:c.windows.end[window]]]))
  }
  
  # Assign cells to windows
  i.windows <- embed(x=1:n.windows, dimension=moving.window)
  pt.windows <- lapply(1:dim(i.windows)[1], function(window) unlist(c.windows[i.windows[window,]]))
  
  # Remove any empty window (can occur in pseudotime.per.window mode)
  if (moving.pt.window) {
    pt.windows <- pt.windows[which(unlist(lapply(pt.windows, length))>0)]
  }
  
  # Figure out min/mean/max pseudotime of each window and name list.
  if (length(name.by) > 1) name.by <- name.by[1]
  namefunc <- get(name.by)
  names(pt.windows) <- round(unlist(lapply(pt.windows, function(window.cells) namefunc(object@pseudotime[window.cells,pseudotime]))), digits=3)
  
  return(pt.windows)
}

pseudotimeMovingWindow <- function(object, pseudotime, cells, moving.window, cells.per.window=NULL, pseudotime.per.window=NULL, name.by=c("mean","min","max")) {
  # Configure function to work based on mode
  # This is a hacky patch...
  if (!is.null(pseudotime.per.window) && is.null(cells.per.window)) {
    moving.pt.window <- T
  } else if (!is.null(cells.per.window)) {
    moving.pt.window <- F
  } else stop("Define either cells.per.window or pseudotime.per.window.")
  
  # Figure out pseudotime of cells
  pt <- object@pseudotime[cells,pseudotime]
  pt.order <- order(pt)
  
  # Figure out window size
  if (moving.pt.window) {
    n.windows <- round((max(pt,na.rm=T)-min(pt,na.rm=T)) / pseudotime.per.window)
    pt.cut <- cut(pt, n.windows, labels = FALSE)
    c.windows <- lapply(1:n.windows, function(x) return(cells[which(pt.cut==x)]))
  } else {
    n.windows <- round(length(cells) / cells.per.window)
    c.windows.end <- round((1:n.windows) * (length(cells) / n.windows))
    c.windows.start <- c(0,head(c.windows.end, -1))+1
    c.windows <- lapply(1:n.windows, function(window) return(cells[pt.order[c.windows.start[window]:c.windows.end[window]]]))
  }
  
  # Assign cells to windows
  i.windows <- embed(x=1:n.windows, dimension=moving.window)
  pt.windows <- lapply(1:dim(i.windows)[1], function(window) unlist(c.windows[i.windows[window,]]))
  
  # Remove any empty window (can occur in pseudotime.per.window mode)
  if (moving.pt.window) {
    pt.windows <- pt.windows[which(unlist(lapply(pt.windows, length))>0)]
  }
  
  # Figure out min/mean/max pseudotime of each window and name list.
  if (length(name.by) > 1) name.by <- name.by[1]
  namefunc <- get(name.by)
  names(pt.windows) <- round(unlist(lapply(pt.windows, function(window.cells) namefunc(object@pseudotime[window.cells,pseudotime]))), digits=3)
  
  # If both cells.per.window & pseudotime.per.window, collapse windows with not enough pseudotime diff
  if (!is.null(cells.per.window) && !is.null(pseudotime.per.window)) {
    window.pt.int <- floor(as.numeric(names(pt.windows)) / pseudotime.per.window)
    pt.windows.new <- lapply(unique(window.pt.int), function(i) {
      unique(unlist(pt.windows[which(window.pt.int == i)]))
    })
    names(pt.windows.new) <- round(unlist(lapply(pt.windows.new, function(window.cells) namefunc(object@pseudotime[window.cells,pseudotime]))), digits=3)
    pt.windows <- pt.windows.new
  }
  
  return(pt.windows)
}

#' Gene Cascade: Fit expression of genes
#' 
#' This takes a group of genes and cells, averages gene expression in groups of cells (determined using a moving window through pseudotime), and then fits either a linear, single sigmoid, or double sigmoid ("impulse") model to describe the expression of each gene. The results are returned for use in \code{\link{geneCascadeImpulsePlots}} to visualize the fits or \code{\link{geneCascadeHeatmap}} to plot gene expression cascades in heatmap format.
#' 
#' Thanks to Yiqun Wang for considerable improvements to the impulse fitting functions.
#' 
#' @param object An URD object
#' @param pseudotime (Character) Name of pseudotime (i.e. a column name of \code{@@pseudotime})
#' @param cells (Character vector) Cells to include
#' @param genes (Character vector) Genes to include
#' @param moving.window (Numeric) Number of bins to use per window
#' @param cells.per.window (Numeric or \code{NULL}) Size of bins (number of cells)
#' @param pseudotime.per.window (Numeric or \code{NULL}) Size of bins (pseudotime)
#' @param scale.data (Logical) If \code{TRUE}, fit impulse model to gene expression that has been scaled to its maximum expression in the chosen cells, otherwise fit to the original data.
#' @param k (Numeric) Number of sets of initial conditions to try for fitting expression of each gene
#' @param pt.windows (List or \code{NULL}) Cells that belong to each pseudotime window -- if \code{NULL} will determine automatically, or provide a list of character vectors of cell ids.
#' @param interpolate (Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points for choosing potential starting conditions. Default (\code{NULL}) does not do any interpolation. Interpolation is NOT used during actual function fitting.
#' @param pulse.only (Logical) If \code{TRUE}, filters out double sigmoid functions that are monotonous and prefers either single sigmoid or linear fit for them instead.
#' @param verbose (Logical) Print major status updates?
#' @param verbose.genes (Logical) Report fitting of each gene?
#' 
#' @return Named list: "pt.windows" ; "mean.expression" Mean expression of each gene in each pseudotime window; "scaled.expression" Scaled mean expression of each gene in each pseudotime window, if \code{scale.data=T}, otherwise same as \code{mean.expression}; "impulse.fits" List of fit parameters for each gene; "timing" (data.frame) of earliest onset time and latest offset time for each gene.
#' 
#' @return (Named List):
#' \itemize{ 
#'   \item{\strong{\code{pt.windows}}: (List) Character vectors of cell IDs of cells in each window, named by either the mean, min, or max pseudotime of those cells (depending on \code{name.by})}
#'   \item{\strong{\code{pt.info}}: (data.frame) Mean, min, max, and range of pseudotime in each window of cells.}
#'   \item{\strong{\code{mean.expression}}: (data.frame) Mean expression of each gene in each pseudotime window}
#'   \item{\strong{\code{scaled.expression}}: (data.frame) Mean expression of each gene in each pseudotime window, scaled to the maximum observed expression}
#'   \item{\strong{\code{mean.smooth}}: (data.frame) Value of fit curve for each pseudotime window (i.e. to predict value of mean.expression)}
#'   \item{\strong{\code{scaled.smooth}}: (data.frame) Value of fit curve for each pseudotime window, scaled to max expression (i.e. to predict value of scaled.expression)}
#'   \item{\strong{\code{impulse.fits}}: (List) Parameters of fit for each gene. See \code{\link{impulseFit}}.}
#'   \item{\strong{\code{timing}}: (data.frame) Earliest onset and latest offset pseudotime determined for each gene.}
#'   \item{\strong{\code{method}}: (Character) Identify fit method ("impulse")}
#'   \item{\strong{\code{fit.scaled}}: (Logical) Was data scaled before fitting?}
#' }
#' 
#' @export
geneCascadeProcess <- function(object, pseudotime, cells, genes, moving.window=3, cells.per.window=NULL, pseudotime.per.window=NULL, scale.data=T, k=50, pt.windows=NULL, interpolate=NULL, pulse.only=T, verbose=T, verbose.genes=F) {
  # Input checking
  genes.dontexist <- setdiff(genes, rownames(object@logupx.data))
  cells.dontexist <- setdiff(cells, colnames(object@logupx.data))
  genes <- intersect(genes, rownames(object@logupx.data))
  cells <- intersect(cells, colnames(object@logupx.data))
  if (length(genes.dontexist) > 0) warning("These genes IDs were not found: ", paste(genes.dontexist))
  if (length(cells.dontexist) > 0) warning("These cell IDs were not found: ", paste(cells.dontexist))
  if (length(genes) == 0) stop("No valid gene IDs were provided.")
  if (length(cells) == 0) stop("No valid cell IDs were provided.")
  # Get moving window of cells by pseudotime
  if (verbose) message(paste0(Sys.time(), ": Calculating moving window expression."))
  if(is.null(pt.windows)){
    pt.windows <- pseudotimeMovingWindow(object, pseudotime=pseudotime, cells=cells, moving.window=moving.window, cells.per.window=cells.per.window, pseudotime.per.window = pseudotime.per.window)
  }
  # Calculate pseudotime parameters for each window
  pt.info <- as.data.frame(t(as.data.frame(lapply(pt.windows, function(cells) {
    pt <- object@pseudotime[cells,pseudotime]
    return(c(mean(pt), min(pt), max(pt), diff(range(pt))))
  }))))
  names(pt.info) <- c("mean","min","max","width")
  rownames(pt.info) <- 1:length(pt.windows)
  
  motifs <- genes[which(genes%in%colnames(object@gene.sig.z))]
  genes_exp <- genes[which(genes%in%rownames(object@logupx.data))]
  mean.motif.expression <- c()
  mean.expression <- c()
  if(length(motifs)>0){
    mean.motif.expression <- as.data.frame(lapply(pt.windows, function(window.cells) apply(t(object@gene.sig.z[window.cells,motifs,drop=F]), 1, mean.of.logs)))
    colnames(mean.motif.expression) <- names(pt.windows)
  }
  if(length(genes_exp)>0){
    mean.expression <- as.data.frame(lapply(pt.windows, function(window.cells) apply(object@logupx.data[genes_exp,window.cells,drop=F], 1, mean.of.logs)))
    colnames(mean.expression) <- names(pt.windows)
  }
  # Make aggregated expression data and scale it
  mean.expression <- rbind(mean.motif.expression, mean.expression)
  names(mean.expression) <- names(pt.windows)
  scale.param <- apply(mean.expression, 1, max)
  scaled.expression <- sweep(mean.expression, 1, scale.param, "/")
  
  if (scale.data) fit.expression <- scaled.expression else fit.expression <- mean.expression
  
  # Do impulse model fitting for all genes
  if (verbose) message(paste0(Sys.time(), ": Fitting impulse model for all genes."))
  impulse.fits <- lapply(genes, function(g) {
    if (verbose.genes) message(paste0(Sys.time(), ":    ", g))
    impulseFit(x=as.numeric(names(fit.expression)), y=as.numeric(fit.expression[g,]), k = k, interpolate=interpolate)
  })
  names(impulse.fits) <- genes
  
  # Get out onset/offset times
  timing <- data.frame(
    time.on=unlist(lapply(impulse.fits, function(x) if(is.list(x)){return(min(x[['time.on']]))}else{return(x['time.on'])})),
    time.off=unlist(lapply(impulse.fits, function(x) if(is.list(x)){return(max(x[['time.off']]))}else{return(x['time.off'])})),
    row.names=genes, stringsAsFactors=F
  )
  
  # Produce impulse fitted parameters for compatibility with plotSmoothFit
  x <- as.numeric(colnames(fit.expression))
  fit.smooth <- as.data.frame(t(as.data.frame(lapply(impulse.fits, function(i) {
    m <- i$model
    if (i$type == 0) {
      m['a']*x + m['b']
    } else if (i$type == 1) {
      impulse.single(x=x, b1 = m['b1'], h0 = m['h0'], h1 = m['h1'], t1 = m['t1'])
    } else if (i$type == 2) {
      impulse.double(x=x, b1 = m['b1'], b2 = m['b2'], h0 = m['h0'], h1 = m['h1'], h2 = m['h2'], t1 = m['t1'], t2 = m['t2'])
    }
  }))))
  colnames(fit.smooth) <- x
  if (scale.data) {
    scaled.smooth <- fit.smooth
    mean.smooth <- sweep(fit.smooth, 1, scale.param, "*")
  } else {
    mean.smooth <- fit.smooth
    scaled.smooth <- sweep(fit.smooth, 1, scale.param, "/")
  }
  
  return(list(
    pt.windows=pt.windows,
    pt.info=pt.info,
    mean.expression=mean.expression,
    scaled.expression=scaled.expression,
    mean.smooth=mean.smooth,
    scaled.smooth=scaled.smooth,
    impulse.fits=impulse.fits,
    timing=timing,
    method="impulse",
    fit.scaled=scale.data
  ))
}

#' Plot Impulse Fits For a Gene Cascade
#' 
#' This plots the results of fitting gene expression data with impulse models.
#' (i.e. from \code{\link{geneCascadeProcess}}.)
#' The windowed data is plotted as black points for each gene, and the fit curve
#' is plotted. Its color indicates the type of fit that
#' was chosen: blue (linear), green (single sigmoid), red (double sigmoid).
#' Diamonds indicate determined gene onset (orange) and gene offset (blue) times.
#' 
#' @param cascade (List) A gene cascade, such as generated by \code{\link{geneCascadeProcess}}
#' @param genes (Character vector) Genes to plot (default \code{NULL} plots all genes that were fit by \code{\link{geneCascadeProcess}})
#' @param ymin0 (Logical) Should the minimum value of the y-axis be 0 (\code{TRUE}) or the minimum value of the fit data (\code{FALSE})
#' @param file (Path) A path to save the plot to as PDF (if NULL, display)
#' @export
geneCascadeImpulsePlots <- function(cascade, genes=NULL, ymin0=T, file=NULL) {
  if(is.null(genes)){
    genes <- names(cascade$impulse.fits)
  }
  ncol <- ceiling(sqrt(length(genes)))
  nrow <- ceiling(length(genes)/ncol)
  if (!is.null(file)) {
    pdf(file=file, width=ncol*4, height=nrow*4)
  }
  par(mfrow=c(nrow,ncol))
  x <- as.numeric(names(cascade$scaled.expression))
  for (g in genes) {
    # Plot gene expression
    if (cascade$fit.scaled) {
      if (ymin0) yl <- c(0, max(cascade$scaled.expression[g,])) else yl <- NULL
      plot(x, cascade$scaled.expression[g,], pch=16, main=g, xlab="Pseudotime", ylab="Expression (scaled)", ylim=yl)
    } else {
      if (ymin0) yl <- c(0, max(cascade$mean.expression[g,])) else yl <- NULL
      plot(x, cascade$mean.expression[g,], pch=16, main=g, xlab="Pseudotime", ylab="Expression", ylim=yl)
    }
    # Get impulse fit parameters
    i <- cascade$impulse.fits[[g]]
    # Determine x-values to use for plotting
    if(length(x)<100){
      x_ <- seq(min(x),max(x),length.out = 500)
    } else {
      x_ <- x
    }
    # Unpack impulse fit parameters
    if (is.list(i)){
      time.on <- i[["time.on"]]
      time.off <- i[["time.off"]]
      i <- c(type=i[["type"]], unlist(i[["model"]]))
    }
    # 
    if (!is.na(i['type'])) {
      if (i['type'] == 0) {
        abline(b=i["a"], a=i["b"], col=rgb(0,0,1,0.7), lwd=5)
      } else if (i['type'] == 1) {
        y_ <- impulse.single(x_, b1=i['b1'], h0=i['h0'], h1=i['h1'], t1=i['t1'])
        lines(x_, y_, col=rgb(0, 1, 0, 0.7), lwd=5)
        if(!is.infinite(time.on)){
          points(time.on, impulse.single(time.on,b1=i['b1'], h0=i['h0'], h1=i['h1'], t1=i['t1']), pch=18, cex=4, col='orange')
        }
        if(!is.infinite(time.off)){
          points(time.off, impulse.single(time.off,b1=i['b1'], h0=i['h0'], h1=i['h1'], t1=i['t1']), pch=18, cex=4, col='blue')
        }
      } else {
        y_ <- impulse.double(x_, b1=i['b1'], b2=i['b2'], h0=i['h0'], h1=i['h1'], h2=i['h2'], t1=i['t1'], t2=i['t2'])
        lines(x_, y_, col=rgb(1, 0, 0, 0.7), lwd=5)
        if(!is.null(time.on)){
          points(time.on, impulse.double(time.on, b1=i['b1'], b2=i['b2'], h0=i['h0'], h1=i['h1'], h2=i['h2'], t1=i['t1'], t2=i['t2']), pch=18, cex=4, col='orange')
        }
        if(!is.null(time.off)){
          points(time.off, impulse.double(time.off, b1=i['b1'], b2=i['b2'], h0=i['h0'], h1=i['h1'], h2=i['h2'], t1=i['t1'], t2=i['t2']), pch=18, cex=4, col='blue')
        }
      }
    }
  }  
  if (!is.null(file)) {
    dev.off()
  }
}

# New version that calculates actual timepoints
#' Plot Gene Cascade Heatmap
#' 
#' @importFrom RColorBrewer brewer.pal
#' @importFrom gplots heatmap.2
#' @importFrom scales gradient_n_pal
#' 
#' @param cascade (list) A processed gene cascade from \code{\link{geneCascadeProcess}}
#' @param color.scale (Character vector) List of color values to use for the heatmap
#' @param add.time (Character or NULL) Either the name of a column of \code{@@meta} that contains actual time information to label the x-axis of the heatmap or NULL to instead use pseudotime
#' @param times.annotate (Numeric vector)
#' @param annotation.list Color bar
#' @param row.font.size (Numeric) The font size of rows (gene names) should scale automatically, but this allows manual adjustment if needed.
#' @param max.font.size (Numeric) This is used to prevent the font from getting big enough that it runs off the page.
#' 
#' @export
geneCascadeHeatmap <- function(cascade, color.scale=RColorBrewer::brewer.pal(9, "YlOrRd"), add.time=NULL, times.annotate=seq(0,1,0.1), title="", annotation.list=NULL, row.font.size=1, max.font.size=0.9) {
  # Correct for NA timings
  timing <- cascade$timing
  timing[intersect(which(is.na(timing$time.on)), which(is.infinite(timing$time.off))), "time.on"] <- Inf
  gene.order <- order(timing$time.on, timing$time.off, na.last=F)
  cols <- scales::gradient_n_pal(color.scale)(seq(0,1,length.out = 50))
  if (!is.null(add.time)) {
    time <- unlist(lapply(cascade$pt.windows, function(cells) mean(object@meta[cells, add.time])))
  } else {
    time <- as.numeric(names(cascade$pt.windows))
  }
  time.lab <- rep("", length(time))
  for (annotate in times.annotate) {
    gt <- which(time >= annotate)
    if (length(gt)>0) time.lab[min(gt)] <- as.character(annotate)
  }
  if (!is.null(annotation.list)) {
    annot <- data.frame(
      gene=rownames(cascade$timing)[gene.order],
      type=factor(rep(NA, length(gene.order)), levels=unique(names(annotation.list))),
      row.names=rownames(cascade$timing)[gene.order],
      stringsAsFactors=F
    )
    for (l in names(annotation.list)) {
      annot[annotation.list[[l]], "type"] <- l
    }
    gplots::heatmap.2(as.matrix(cascade$scaled.expression[gene.order,]), Rowv=F, Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", key=F, labCol=time.lab, RowSideColors=as.character(annot$type), cexCol=0.8, cexRow=min((3.4-0.56*log(length(gene.order)))*row.font.size, max.font.size), margins = c(5,8), lwid=c(0.2,4), lhei=c(0.4, 4))
  } else {
    gplots::heatmap.2(as.matrix(cascade$scaled.expression[gene.order,]), Rowv=F, Colv=F, dendrogram="none", col=cols, trace="none", density.info="none", key=F, labCol=time.lab, cexCol=0.8, cexRow=min((3.4-0.56*log(length(gene.order)))*row.font.size, max.font.size), margins = c(5,8), lwid=c(0.3,4), lhei=c(0.4, 4))
  }
  title(title, line=1, adj=0.4)
}


