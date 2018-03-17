#' @import Matrix
#' @import ggplot2
#' @importClassesFrom destiny DiffusionMap
NULL

#' URD class
#' 
#' @importClassesFrom destiny DiffusionMap
#' @importClassesFrom Matrix dgCMatrix
#' @exportClass URD
URD <- methods::setClass("URD", slots = c(
  count.data=c("dgCMatrix", NULL), 
  logupx.data=c("dgCMatrix", NULL), 
  meta="data.frame", 
  group="vector", 
  group.ids="data.frame", 
  active.group="character",
  pca.obj="list", 
  pca.sdev="vector", 
  pca.load="data.frame", 
  pca.scores="data.frame", 
  pca.sig="vector", 
  tsne.y="data.frame", 
  cell.names="vector", 
  gene.sig.z="data.frame", 
  dm=c("DiffusionMap",NULL), 
  var.genes="vector", 
  knn="list",
  plot.3d="list",
  diff.data="data.frame", 
  pseudotime="data.frame",
  pseudotime.stability="list", 
  tree="list",
  nmf.g=c("dgCMatrix", "matrix", NULL),
  nmf.c=c("dgCMatrix", "matrix", NULL)
))

#' URD show method
#'
#' Prevents R from crashing by trying to print all slots of an URD object
#' if a returned object is not stored in a variable.
#' @param object An URD object
#' @name show
#' @docType methods
#'
setMethod(
  f = "show",
  signature = "URD",
  definition = function(object) {
    cat(
      "URD object:",
      nrow(object@logupx.data),
      "genes x",
      ncol(object@logupx.data),
      "cells.\n"
    )
    invisible(NULL)
  }
)

#' Create a new URD object
#' 
#' Creates a new URD object. Provide expression data as UMI counts and (optionally)
#' additional metdata. This function will perform basic filtering of the data
#' (though the user may want to perform more advanced filtering before creating
#' the URD object, using \code{ds.meta.trim}), will normalize and log2 transform the
#' data, and will return an URD object with the original data in slot \code{count.data},
#' the normalized log-transformed data in slot \code{logupx.data}, the metadata in slot
#' \code{meta}, and an initial grouping of the data drawn from the cell names, up to the
#' first dash (-) or underscore (_) in slot \code{group.ids}.
#' 
#' @importClassesFrom Matrix dgCMatrix
#' 
#' @param count.data (Matrix or dgCMatrix) UMI expression data, with rows as genes and columns as cells
#' @param meta (data.frame) Metadata, with rows as cells (row names should match column names of \code{count.data})
#' @param min.cells (Numeric) Minimum number of cells that must express a gene to retain it
#' @param min.genes (Numeric) Minimum number of genes a cell must express to retain it
#' @param min.counts (Numeric) Minimum number of UMIs detected for a gene across the entire data to retain it
#' @param gene.max.cut (Numeric) Maximum number of UMIs observed for a gene in a single cell to retain it
#' @param max.genes.in.ram (Numeric) Number of genes to normalize and log-transform at a time (to prevent running out of memory as matrices become non-sparse during the process)
#' 
#' @return An URD object
#' 
#' @export
createURD <- function(count.data, meta=NULL, min.cells=3, min.genes=500, min.counts=10, gene.max.cut=5000, max.genes.in.ram=5000, verbose=T) {

  # Filter cells by nGenes
  if (verbose) message(paste0(Sys.time(), ": Filtering cells by number of genes."))
  num.genes <- apply(count.data, 2, function(x) sum(x > 0))
  names(num.genes) <- colnames(count.data)
  cells.enough.genes <- names(num.genes[which(num.genes>min.genes)])
  shhhh <- gc()
  # Filter genes by nCells
  if (verbose) message(paste0(Sys.time(), ": Filtering genes by number of cells."))
  num.cells <- apply(count.data[,cells.enough.genes], 1, function(x) sum(x > 0))            
  genes.enough.cells <- names(num.cells[which(num.cells>min.cells)])
  shhhh <- gc()
  # Filter genes by total counts
  if (verbose) message(paste0(Sys.time(), ": Filtering genes by number of counts across entire data."))
  num.counts <- apply(count.data[,cells.enough.genes], 1, sum)
  genes.enough.counts <- names(num.counts[which(num.counts>min.counts)])
  shhhh <- gc()
  # Filter genes by maximum expression
  if (verbose) message(paste0(Sys.time(), ": Filtering genes by maximum observed expression."))
  maxgene <- apply(count.data[, cells.enough.genes], 1, max)
  genes.not.above.max <- rownames(count.data)[maxgene <= gene.max.cut]
  shhhh <- gc()
  # Figure out genes to retain
  genes.use <- intersect(intersect(genes.enough.cells, genes.enough.counts), genes.not.above.max)

  # Create an URD object
  if (verbose) message(paste0(Sys.time(), ": Creating URD object."))
  object <- new("URD", count.data=as(count.data[genes.use, cells.enough.genes], "dgCMatrix"))
  shhhh <- gc()
  
  # Determine normalization factor
  if (verbose) message(paste0(Sys.time(), ": Determining normalization factors."))
  cs <- apply(object@count.data, 2, sum)
  norm_factors <- (10**ceiling(log10(median(cs))))/cs
  #if (verbose) message(summary(norm_factors))
  shhhh <- gc()
  
  # Log-normalize the data
  if (verbose) message(paste0(Sys.time(), ": Normalizing and log-transforming the data."))
  n.chunks <- ceiling(ncol(object@count.data) / max.genes.in.ram)
  n.cells <- ncol(object@count.data)
  lognorm.chunks <- lapply(1:n.chunks, function(chunk) {
    shhhh <- gc()
    i <- (chunk-1) * max.genes.in.ram + 1
    j <- min((chunk * max.genes.in.ram), n.cells)
    as(round(log2(sweep(object@count.data[,i:j], 2, norm_factors[i:j], "*")+1), digits=2), "dgCMatrix")
  })
  
  # Generate @logupx.data
  shhhh <- gc()
  object@logupx.data <- do.call(what = "cbind", lognorm.chunks)
  
  # Set up the grouping IDs
  if (verbose) message(paste0(Sys.time(), ": Finishing setup of the URD object."))
  initial.group <- factor(unlist(lapply(colnames(object@logupx.data),function(x) strsplit(x,"_|-")[[1]][1] )))
  names(initial.group) <- colnames(object@logupx.data)
  object@group.ids <- data.frame(initial.group)
  names(object@group.ids) <- "init"
  object@active.group <- "init"
  
  # Set up the metadata
  object@meta <- data.frame(n.Genes=num.genes[colnames(object@count.data)])
  object@meta[,"n.Trans"] <- apply(object@count.data, 2, sum)
  if (!is.null(meta)) {
    object@meta <- cbind(object@meta, meta[rownames(object@meta),])
  }
  
  if (verbose) message(paste0(Sys.time(), ": All done."))
  return(object)
} 

#' Default URD continuous colors
#' 
#' @return Vector of default colors
defaultURDContinuousColors <- function(with.grey=F) {
  if (with.grey) {
    return(c("#B2B2B2", "#9BABC2", "#7D9FD1", "#5A90E0", "#307DF0", "#0065FF", 
             "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
             "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
             "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
             "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"))
  } else {
    return(c("#0000FF", "#0013FF", "#0028FF", "#003CFF", "#0050FF", "#0065FF", 
           "#0078FF", "#008DFF", "#00A1FF", "#00B5FF", "#00CAFF", "#00DEFF", 
           "#00F2FF", "#27FFD7", "#8CFF71", "#F1FF0D", "#FFEE00", "#FFDB00", 
           "#FFC900", "#FFB700", "#FFA500", "#FF9200", "#FF8000", "#FF6D00", 
           "#FF5B00", "#FF4800", "#FF3600", "#FF2400", "#FF1200", "#FF0000"))
  }
}
