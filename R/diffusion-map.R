#' Import a pre-calculated diffusion map
#' 
#' @param object An URD object
#' @param dm A DiffusionMap object produced by destiny's \code{\link[destiny]{DiffusionMap}}.
#' @return An URD object
#' @export
importDM <- function(object, dm, cell.names=NULL) {
            # Make sure things are named
            if (is.null(cell.names)) {
              if (!is.null(rownames(dm@eigenvectors))) cell.names <- rownames(dm@eigenvectors)
              else if (!is.null(rownames(dm@transitions))) cell.names <- rownames(dm@transitions)
              else if (nrow(dm@eigenvectors) == ncol(object@logupx.data)) {
                warning("Cell names not provided, but diffusion map matches dimensions of object data, so assuming cell names are the same.")
                cell.names <- colnames(object@logupx.data)
              } else {
                stop("Please provide names of cells used to generate the diffusion map.")
              }
            }  
            if (is.null(rownames(dm@transitions))) rownames(dm@transitions) <- cell.names
            if (is.null(colnames(dm@transitions))) colnames(dm@transitions) <- cell.names
            if (is.null(rownames(dm@eigenvectors))) rownames(dm@eigenvectors) <- cell.names
            # Load diffusion map into the Dropseq object
            object@dm <- dm
            # Initialize the diff.data data frame
            object@diff.data <- data.frame(
              row.names = rownames(dm@eigenvectors)
            )
            # Initialize the pseudotime data frame
            object@pseudotime <- data.frame(
              row.names = rownames(dm@eigenvectors)
            )
            return(object)
}

#' Calculate a diffusion map
#' 
#' @importFrom destiny DiffusionMap find_dm_k find_sigmas
#' @param object An URD object
#' @param genes.use (Character vector) Genes to use (default is variable genes, NULL uses all genes)
#' @param cells.use (Character vector) Cells to use (default is NULL, which uses all cells)
#' @param knn (Numeric) Number of nearest neighbors to use (NULL takes a guess) 
#' @param sigma.use (NULL, numeric, or "local") Kernel width to use for the diffusion map (NULL uses destiny's global auto-detection procedure)
#' @param n_local (Numeric) If \code{sigma.use="local"}, the distance to the n-local(th) nearest neighbors are used for determining the kernel width for each cell.
#' @param distance (Character) Distance metric to use for determining transition probabilities.
#' @param density.norm (Logical) If TRUE, use density normalization
#' @param dcs.store (Numeric) How many diffusion components to store
#' @param verbose (Logical) Report determined values for auto-detected parameters?
#' @export
calcDM <- function(object, genes.use=object@var.genes, cells.use=NULL, knn=NULL, sigma.use=NULL, n_local=5:7, distance=c("euclidean", "cosine", "rankcor"), density.norm=T, dcs.store=200, verbose=T) {
  
  # Subset the data and convert to a matrix without row or column names because they crash destiny.
  if (is.null(genes.use) || length(genes.use) == 0) genes.use <- rownames(object@logupx.data)
  if (is.null(cells.use)) cells.use <- colnames(object@logupx.data)
  data.use <- t(object@logupx.data[genes.use, cells.use])
  rownames(data.use) <- NULL
  colnames(data.use) <- NULL
  data.use <- as.matrix(data.use)
  
  # Figure out sigma
  if (is.null(sigma.use)) {
    sigma.use <- find_sigmas(data.use, steps=25, verbose=F)@optimal_sigma
    if (verbose) print(paste("destiny determined an optimal global sigma of", round(sigma.use, digits=3)))
  } else if (is.numeric(sigma.use)) {
    if (verbose) print(paste("Using provided global sigma", sigma.use))
  } else if (sigma.use == "local") {
    if (verbose) print(paste("Using local sigma."))
  } else {
    stop("sigma must either be NULL, 'local', or numeric.")
  }
  
  # Figure out k-nearest neighbors
  if (is.null(knn)) {
    knn <- find_dm_k(length(cells.use))
    if (verbose) print(paste("destiny will use", knn, "nearest neighbors."))
  }
  
  # Calculate the diffusion map
  dm <- DiffusionMap(data.use, sigma=sigma.use, k=knn, n_eigs = dcs.store, density_norm = density.norm, distance=distance[1])
  
  # Name things because you have to remove all of the names
  rownames(dm@eigenvectors) <- cells.use
  rownames(dm@transitions) <- cells.use
  colnames(dm@transitions) <- cells.use
  
  # Store the resulting diffusion map
  object <- importDM(object, dm)
  
  return(object)
}