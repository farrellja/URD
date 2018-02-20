#' Get z-scored data
#' 
#' Returns z-scored data, calculated on the log-transformed, normalized data. Eliminates
#' any genes that produce NA values, which can occur if they have standard deviation 0.
#' @param object An URD object
#' @param genes (Character vector) Genes to return in z-scored matrix (default, all genes)
#' @return Matrix of z-scored data (Warning: memory hog as data will no longer be sparse.)
#' @export
get.z.data <- function(object, genes=NULL) {
  # Z-scored data
  if (is.null(genes)) genes <- rownames(object@logupx.data)
  data.mean <- apply(object@logupx.data[genes,], 1, mean)
  data.sd <- apply(object@logupx.data[genes,], 1, sd)
  z.data <- as.matrix(sweep(sweep(object@logupx.data[genes,], 1, data.mean, "-"), 1, data.sd, "/"))
  z.data <- z.data[complete.cases(z.data),]
  return(z.data)
}

#' Get binary data
#' 
#' Returns a sparse matrix of 0 and 1 showing whether any expression of each gene
#' was observed in each cell.
#' @param object An URD object
#' @param genes (Character vector) Genes to return in binarized matrix (default, all genes)
#' @return Sparse matrix (dgCMatrix)
#' @export
get.binary.data <- function(object, genes=NULL) {
  # Binarize data
  if (is.null(genes)) genes <- rownames(object@count.data)
  binary.data <- object@count.data[genes,]
  binary.data[binary.data > 1] <- 1
  return(binary.data)
}

#' Get UPX (un-logged) data
#' 
#' Returns a sparse matrix of expression values that are normalized for library
#' size (UPX normalization), but are not log transformed.
#' @param object An URD object
#' @param genes (Character vector) Genes to return in UPX atrix (default, all genes)
#' @return Sparse matrix (dgCMatrix)
#' @export
get.upx.data <- function(object, genes=NULL) {
  # Get non-log upx data
  if (is.null(genes)) genes <- rownames(object@logupx.data)
  return(as((2^(object@logupx.data) - 1), "dgCMatrix"))
}
