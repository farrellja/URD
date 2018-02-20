#' Calculate tSNE projection of data
#' 
#' 
#' @importFrom Rtsne Rtsne
#' @param object An URD object
#' @param dim.use (Character) Whether to calculate the tSNE projection from PCs (\code{"pca"}) or diffusion components (\code{"dm"}). Default is from PCs.
#' @param which.dims (Numeric vector) Which PCs (or diffusion components) to use. Defaults to the significant PCs. (The default will probably work with diffusion components, though it is non-sensical in that case.)
#' @param perlexity (Numeric) Perplexity parameter for the tSNE
#' @param theta (Numeric) Speed/accuracy trade-off for Barnes-Hut approximation of tSNE (0-1, 0 is exact tSNE, higher is less accurate, default is 0.5)
#' @param max_iter (Numeric) Number of iterations to perform
#' @param verbose (Logical) Should Rtsne print progress updates?
#' @return An URD object with tSNE coordinates stored in the \code{tsne.y} slot.
#' @export
calcTsne <- function(object, dim.use=c("pca", "dm"), which.dims=which(object@pca.sig), perplexity=30, theta=0.5, max_iter=1000, verbose=FALSE) {
  if (length(dim.use) > 1) dim.use <- dim.use[1]
  if (dim.use == "pca") {
    if (length(which.dims) == 0) {
      stop ("No PCs are marked as significant. Specify dimensions to use (which.dims=).")
    }
    cells.use <- rownames(object@pca.scores)
    tsne.result <- Rtsne(as.matrix(object@pca.scores[,which.dims]), dims=2, pca=FALSE, perplexity=perplexity, theta=theta, max_iter=max_iter, verbose=verbose)
  } else if (dim.use == "dm") {
    cells.use <- intersect(colnames(object@logupx.data), rownames(object@dm.obj[[1]]@eigenvectors))
    tsne.result <- Rtsne(as.matrix(object@dm.obj[[1]]@eigenvectors[cells.use,which.dims]), dims=2, pca=FALSE, perplexity=perplexity, theta=theta, max_iter=max_iter, verbose=verbose)
  } else {
    stop("dim.use must be either pca or dm")
  }
  object@tsne.y <- data.frame(tsne.result$Y)
  names(object@tsne.y) <- c("tSNE1", "tSNE2")
  rownames(object@tsne.y) <- cells.use
  return(object)
}
