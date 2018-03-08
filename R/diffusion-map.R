#' Import a pre-calculated diffusion map
#' 
#' If destiny is run separately (such as on a cluster to vary parameters), this
#' allows quick addition of destiny's output to an URD object.
#' 
#' @param object An URD object
#' @param dm A DiffusionMap object produced by destiny's \code{\link[destiny]{DiffusionMap}}.
#' 
#' @return An URD object with the diffusion map stored in \code{@@dm}.
#' 
#' @examples
#' object <- importDM(object, dm.8)
#' 
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
#' URD uses the R package destiny to calculate a diffusion map.
#' 
#' The diffusion map
#' is stored in the \code{@@dm} slot of the URD object, and the transition probabilities
#' (\code{@@dm@@transitions}) are used extensively downstream in URD's pipeline.
#' Important parameters to vary while choosing a diffusion map to use are the number
#' of nearest neighbors considered (the square root of the number of cells
#' is a good starting point, or \code{NULL} will trigger destiny to try to determine
#' a good starting point), and the sigma of the Gaussian used to transform cells'
#' distance.
#' 
#' @importFrom destiny DiffusionMap find_dm_k find_sigmas
#' 
#' @param object An URD object
#' @param genes.use (Character vector) Genes to use (default is variable genes, NULL uses all genes)
#' @param cells.use (Character vector) Cells to use (default is NULL, which uses all cells)
#' @param knn (Numeric) Number of nearest neighbors to use (NULL takes a guess) 
#' @param sigma.use (NULL, numeric, or "local") Kernel width to use for the diffusion map (NULL uses destiny's global auto-detection procedure, and "local" determines a sigma for each cell based on the distance to its nearest neighbors.)
#' @param n_local (Numeric) If \code{sigma.use="local"}, the distance to the n-local(th) nearest neighbors are used for determining the kernel width for each cell.
#' @param distance (Character) Distance metric to use for determining transition probabilities.
#' @param density.norm (Logical) If TRUE, use density normalization
#' @param dcs.store (Numeric) How many diffusion components to store
#' @param verbose (Logical) Report determined values for auto-detected parameters?
#' 
#' @seealso \link[destiny]{DiffusionMap}
#' 
#' @return An URD object with the diffusion map stored in \code{@@dm}.
#' 
#' @examples 
#' # Allow destiny to determine most parameters on its own
#' object <- calcDM(object, knn=NULL, sigma.use=NULL)
#' 
#' # Refine the diffusion map by adjusting parameters manually
#' object <- calcDM(object, knn = 200, sigma.use = 8)
#' 
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

#' Create an edge list from the transition matrix of the diffusion map
#'
#' Returns a list of edges with transition probability from the diffusion
#' map in an URD object. Used by some plotting functions to show the
#' connectivity of the network if desired.
#'
#' @importFrom igraph graph_from_adjacency_matrix as_data_frame V induced_subgraph
#' 
#' @param object An URD object
#' @param cells (Character vector) Cells to calculate edges between (default \code{NULL} uses all cells in diffusion map.)
#' @param edges.return (Numeric) Number of edges to return (default \code{NULL} returns all edges.)
#' 
#' @return A data.frame with columns from (cell name), to (cell name), and weight (transition probability)
#' 
#' @examples 
#' edges <- edgesFromDM(object, cells=NULL, edges.return=100000)
#' 
#' @export
edgesFromDM <- function(object, cells=NULL, edges.return=NULL) {
  # Create igraph representation from transition matrix
  ig <- graph_from_adjacency_matrix(object@dm@transitions, weighted=T, mode="undirected")
  
  # Subset if only some cells are provided
  if (!is.null(cells)) {
    v.keep <- which(names(V(ig)) %in% cells)
    ig <- induced_subgraph(ig, vids=v.keep)
  }
  
  # Convert igraph representation to data frame
  ig_edges <- as_data_frame(ig)
  
  # Sample number of edges at random, if desired
  if (!is.null(edges.return) && (edges.return < nrow(ig_edges))) ig_edges <- ig_edges[sample(1:nrow(ig_edges), edges.return),]
  
  return(ig_edges)
}