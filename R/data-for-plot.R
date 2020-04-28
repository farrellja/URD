#' Retrieve data for plotting
#' 
#' This retrieves data to annotate points on various plots and returns either
#' the data, or the data transformed into colors for plotting. This function is
#' called by most of URD's plotting routines.
#' 
#' By default, it searches for a label in the following order:
#' metadata (\code{"meta"}), group IDs (\code{"group"}), gene signatures (\code{"sig"}), 
#' gene expression (\code{"gene"}), unnormalized gene expression counts (\code{"counts"}),
#' pseudotime (\code{"pseudotime"}), NMF modules (\code{"nmf"}), principal components (\code{"pca"}), 
#' simulated diffusion results (\code{"diff.data"}).
#' However, data can be requested from a particular source by setting \code{label.type}. 
#' This can be used to get log-transformed but unnormalized expression data from
#' \code{count.data} by setting \code{label.type="counts"}.
#' 
#' @importFrom scales cscale gradient_n_pal hue_pal squish
#' 
#' @param object An URD object
#' @param label (Character) The label of the data to search for
#' @param label.type (Character) Where to look for the data. Default is "search" which looks in order: "meta", "group", "sig", "gene", "counts", "pseudotime", "pca", "diff.data"
#' @param cells.use (Character vector) Which cells to return information for (default is NULL, which returns all cells)
#' @param as.color (Logical) Return hex color values instead of the raw data
#' @param as.single.color (Logical) Return on a scale from 0-1 as a "single color" to feed into \code{rgb()} for more complex color generation.
#' @param as.discrete.list (Logical) If TRUE, returns a list (see below).
#' @param continuous.colors (Character vector) Colors to use to produce a continuous color scale (used if data is continuous) and \code{as.color=T}.
#' @param continuous.color.limits (Numeric vector, length 2) Limits to use for scaling continuous color data. Data outside the range is squished into range. If \code{NULL}, uses the range of the data.
#' @param colors.use (Character vector) Colors to use for discrete data (default is NULL, which will use palette) 
#' @param palette (Function) A palette (see \code{\link[grDevices]{Palettes}}) to use to generate colors for discrete data
#' @return The value returned depends on the flags set. A vector is returned if \code{as.discrete.list=FALSE} or a list is returned if \code{as.discrete.list=TRUE} (see below). The data is returned directly if \code{as.color=F} and \code{as.single.color=F}, or it is returned as color values to plot directly if either are \code{TRUE}.
#' @section Discrete List Elements (if \code{as.discrete.list=T}):
#' \describe{
#'   \item{$discrete}{(Logical) Is the data discrete or continuous?}
#'   \item{$data}{(Numeric or Character Vector) Either the data itself, or colors if \code{as.color=T} or \code{as.single.color=T}.}
#'   \item{$legend}{(Named Character Vector) Colors used, named by the data, returned if \code{as.color=T} and data is discrete}
#'   \item{$range}{(Numeric Vector) If data is continuous, the range of the data, otherwise NULL.}
#' }
#' @export
data.for.plot <- function(object, label, label.type=c("search", "meta", "group", "sig", "gene", "counts", "pseudotime", "nmf", "pca", "diff.data"), cells.use=NULL, as.color=F, as.single.color=F, as.discrete.list=F, continuous.colors=NULL, continuous.color.limits=NULL, colors.use=NULL) {
  
  # Default URD colors
  if (is.null(continuous.colors)) continuous.colors <- defaultURDContinuousColors()
  
  # Search order:
  # meta, group.ids, signatures, genes (logupx), genes (counts), pseudotime, PCA, diff.data
  
  # Case insensitive matching
  label.type <- tolower(label.type[1])
  
  # Check metadata
  if (label.type=="meta" | (label.type=="search" & label %in% colnames(object@meta))) {
    if (is.null(cells.use)) cells.use <- rownames(object@meta)
    data <- object@meta[cells.use,label]
    discrete <- "check"
  }
  
  # Check group.ids
  else if (label.type=="group" | (label.type=="search" & label %in% colnames(object@group.ids))) {
    if (is.null(cells.use)) cells.use <- rownames(object@group.ids)
    data <- object@group.ids[cells.use,label]
    discrete <- T
  }
  
  # Check signatures
  else if (label.type=="sig" | (label.type=="search" & label %in% colnames(object@gene.sig.z))) {
    if (is.null(cells.use)) cells.use <- rownames(object@gene.sig.z)
    data <- object@gene.sig.z[cells.use,label]
    discrete <- F
  }
  
  # Check genes
  else if (label.type=="gene" | (label.type=="search" & label %in% rownames(object@logupx.data))) {
    if (is.null(cells.use)) cells.use <- colnames(object@logupx.data)
    data <- object@logupx.data[label,cells.use]
    discrete <- F
  }
  
  # Check counts
  else if (label.type=="counts" | (label.type=="search" & label %in% rownames(object@count.data))) {
    if (is.null(cells.use)) cells.use <- colnames(object@count.data)
    data <- log2(object@count.data[label,cells.use]+1) 
    discrete <- F
  }
  
  # Check pseudotime
  else if (label.type=="pseudotime" | (label.type=="search" & label %in% colnames(object@pseudotime))) {
    if (is.null(cells.use)) cells.use <- rownames(object@pseudotime)
    data <- object@pseudotime[cells.use, label]
    discrete <- F
  }
  
  # Check NMF modules
  else if (label.type=="nmf.c1" | (label.type=="search" & label %in% colnames(object@nmf.c1))) {
    if (is.null(cells.use)) cells.use <- rownames(object@nmf.c1)
    data <- object@nmf.c1[cells.use, label]
    discrete <- F
  }
  
  # Check PC
  else if (label.type=="pca" | (label.type=="search" & label %in% colnames(object@pca.scores))) {
    if (is.null(cells.use)) cells.use <- rownames(object@pca.scores)
    data <- object@pca.scores[cells.use, label]
    discrete <- F
  }
  
  # Check diff.data
  else if (label.type=="diff.data" | (label.type=="search" & label %in% colnames(object@diff.data))) {
    if (is.null(cells.use)) cells.use <- rownames(object@diff.data)
    data <- object@diff.data[cells.use, label]
    discrete <- F
  }
  
  # Uh oh
  else { stop(paste("Cannot find", label, "in metadata, group.ids, signatures, genes, NMF modules, PCA, or pseudotime."))}
  
  data.range <- NULL
  
  # Is it discrete?
  if (discrete == "check") {
    if (class(data) == "factor") {
      discrete <- T
    } else if (class(data) == "numeric") {
      discrete <- F
    } else {
      discrete <- any(is.na(suppressWarnings(as.numeric(data))))  
    }
  }
  
  # If data is logical, convert to character
  if (class(data) == "logical") data <- as.character(data)
  
  # Convert to color values if desired
  if (as.color || as.single.color) {
    if (discrete) {
      data.to.color <- sort(unique(data))
      if (is.null(colors.use)) colors <- scales::hue_pal()(length(data.to.color)) else colors <- colors.use
      if (is.null(names(colors))) names(colors) <- data.to.color
      legend <- colors
      data <- colors[data]
    } else {
      data.range <- range(data)
      if (!as.single.color) {
        if (is.null(continuous.color.limits)) continuous.color.limits <- range(as.numeric(data))
        data <- cscale(squish(as.numeric(data), range=continuous.color.limits), gradient_n_pal(continuous.colors))
      } else {
        data <- (data - data.range[1]) / data.range[2]
      }
      legend <- NULL
    }
  } else {
    legend <- NULL
  }
  
  names(data) <- cells.use
  
  if (as.discrete.list) {
    return(list(discrete=discrete, data=data, legend=legend, range=data.range))
  } else {
    return(data)
  }
}

#' Default URD continuous colors
#' 
#' @return Vector of default colors
#' 
#' @export
#' 
#' @keywords internal
defaultURDContinuousColors <- function(with.grey=F, symmetric=F, evenly.spaced=F) {
  if (symmetric) {
    return(c("#B3006B", "#C04981", "#CA7197", "#D396AE", "#D9B9C5", "#DDDDDD", 
             "#BFD6B7", "#A0CE91", "#7FC66B", "#57BD43", "#0CB300"))
  } else if (with.grey) {
    return(c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
             "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
             "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
             "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
             "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"))
  } else if (evenly.spaced) { 
    return(c("#0000FF", "#0023FF", "#0046FF", "#0069FF", "#008CFF", "#00AFFF", 
             "#00D2FF", "#00F6FF", "#1AFFE4", "#3DFFC1", "#60FF9D", "#83FF7A", 
             "#A6FF57", "#CAFF34", "#EDFF11", "#FFF800", "#FFEC00", "#FFDF00", 
             "#FFD300", "#FFC700", "#FFBA00", "#FFAE00", "#FF9F00", "#FF8800", 
             "#FF7100", "#FF5A00", "#FF4300", "#FF2D00", "#FF1600", "#FF0000"))
  } else {
    return(c("#0000FF", "#0013FF", "#0028FF", "#003CFF", "#0050FF", "#0065FF", 
             "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
             "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
             "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
             "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"))
  }
}