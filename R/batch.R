#' Batch correct data using ComBat
#' 
#' @param object An URD object
#' @param batch.column (Character) Name of column in \code{@@meta} that encodes batch information.
#' @param adjustment.columns (Character vector) Name of columns in \code{@@meta} that encode information to make sure not to correct (e.g. developmental stage, etc.)
#' @param reference.batch (Character) Name of reference batch to hold constant (default \code{NULL} adjusts all batches toward each other)
#' @param genes.correct (Character vector) Which genes to correct (default \code{NULL} corrects all genes)
#' @param thresh.out (Numeric) To keep matrix sparse, change values less than this to \code{min.out}
#' @param min.out (Numeric) Minimum gene expression value to return
#' @param max.out (Numeric) Maximum gene expression value to return (Default: capped at maximum gene expression value of original data + 1.)
#' @param parametric.prior (Logical) ComBat option: Assume normal distributed gene expression (Warning: extremely slow if FALSE)
#' @param prior.plots (Logical) ComBat option: Plots to compare empirical and parametric batch (see \code{\link[sva]{ComBat}} for more info.)
#' @param verbose (Logical) Print progress
#'
#' @return An URD object with logupx.data batch-corrected using ComBat.
#'
#' @export
batchComBat <- function(object, batch.column, adjustment.columns=NULL, reference.batch=NULL, genes.correct=NULL, thresh.out=.05, min.out=0, max.out=(ceiling(max(object@logupx.data))+1), parametric.prior=T, prior.plots=F, verbose=F) {
  # Bioconductor package 'sva' is required for ComBat batch correction.
  if (requireNamespace("sva", quietly=T)) {
    
    # Build model matrix if adjustment columns exist
    if (verbose) print(paste0(Sys.time(), ": Building model matrix."))
    if (!is.null(adjustment.columns)) {
      matrix.call <- paste0("model.matrix(~", paste(adjustment.columns, collapse="+"), ", data=object@meta)")
      model.matrix <- eval(parse(text=matrix.call))
    } else {
      model.matrix <- NULL
    }
    
    # Get genes.correct
    if (is.null(genes.correct)) genes.correct <- rownames(object@logupx.data)
    
    # Filter genes to those with variance > 0
    if (verbose) print(paste0(Sys.time(), ": Removing genes with no variance."))
    if (is.null(reference.batch)) {
      var <- apply(object@logupx.data[genes.correct,], 1, var)
      genes.use <- names(which(var>0))
    } else {
      ref.cells <- rownames(object@meta)[which(object@meta[,batch.column] == reference.batch)]
      var <- apply(object@logupx.data[genes.correct,ref.cells], 1, var)
      genes.use <- names(which(var>0))
    }
    if (verbose) print(paste("Retaining", length(genes.use), "of", nrow(object@logupx.data), "genes."))
    
    # Run ComBat
    if (verbose) print(paste0(Sys.time(), ": Starting ComBat."))
    batch.data <- as.matrix(sva::ComBat(dat=as.matrix(object@logupx.data[genes.use,]), batch=object@meta[colnames(object@logupx.data), batch.column], mod=model.matrix, prior.plots=prior.plots, par.prior=parametric.prior, ref.batch=reference.batch))
    
    # Crop extreme values occasionally produced by Combat
    if (verbose) print(paste0(Sys.time(), ": Correcting extreme values and re-sparsifying matrix."))
    batch.data[batch.data < thresh.out] <- min.out
    batch.data[batch.data > max.out] <- max.out
    
    # Replace logupx data with batch-corrected data.
    object@logupx.data <- as(round(batch.data, digits=4), "dgCMatrix")
    return(object)
  } else {
    stop("Package sva from bioConductor is required for this function. To install:\nsource('https://bioconductor.org/biocLite.R')\nbiocLite('sva')\n")
  }
}