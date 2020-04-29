#' @import Matrix
#' @import ggplot2
#' @importClassesFrom destiny DiffusionMap
NULL

#' URD class
#' 
#' All information associated with a reconstruction project is stored in an
#' URD object. Most functions in the URD package take one as input, and many of
#' them return an URD object that has been operated on. A new URD object is
#' created using the \code{\link{createURD}} function.
#' 
#' @slot count.data (Sparse Matrix) The initially provided count data expression matrix
#' @slot logupx.data (Sparse Matrix) The normalized and log2-transformed expression matrix
#' @slot meta (data.frame) Continuous and categorical information about cells provided during object creation
#' @slot group.ids (data.frame) Categorial information about cells, such as cluster identities and assignment within the tree
#' @slot var.genes (Character vector) Genes that have been determined as highly variable
#' @slot pca.load (data.frame) Loading of genes (rows) into principal components (columns)
#' @slot pca.score (data.frame) Score for each principal component (columns) in each cell (rows)
#' @slot pca.sdev (Numeric vector) Standard deviation of each principal component
#' @slot pca.sig (Logical vector) Whether each PC has been determined to be significant
#' @slot tsne.y (data.frame) tSNE coordinates for each cell (rows)
#' @slot dm (DiffusionMap from destiny) Diffusion map calculated by destiny
#' @slot diff.data (data.frame) Visitation frequencies from 
#' @slot pseudotime (data.frame) Determined pseudotime for each cell (row)
#' @slot pseudotime.stability (List) Contains pseudotime and visitation for subsets of floods/walks for assaying whether enough simulations have been performed
#' \describe{
#'   \item{\code{pseudotime}}{Determined pseudotime for each cell}
#'   \item{\code{walks.per.cell}}{Number of walks that visited each cell}
#' }
#' @slot tree (List) Contains information from building the tree
#' \describe{
#'   \item{\code{tips}}{(Vector) The tips used as input into building the tree}
#'   \item{\code{cells.in.tip}}{(List of character vectors) The cells that belong to each tip}
#'   \item{\code{pseudotime}}{(Named vector) The pseudotime of each cell used in construction of the tree}
#'   \item{\code{segments}}{(Vector) All segments present in the final tree}
#'   \item{\code{segment.pseudotime.limits}}{(data.frame) Start and end pseudotimes for each tree segment}
#'   \item{\code{segment.joins}}{(data.frame) Parent and child relationships for all segments in the tree and pseudotime of their breakpoints}
#'   \item{\code{segment.joins.initial}}{(data.frame) Parent and child relationships during construction of the tree -- this data is prior to collapsing any multifurcating branchpoints or removing any segments that are too short or assigned too few cells}
#'   \item{\code{pseudotime.breakpoint.details}}{(List) Contains information used during the determination of putative pseudotime breakpoints between each pair of segments. Only retained if \code{save.all.breakpoint.info=T} during \code{buildTree}}
#'   \item{\code{segment.divergence}}{(data.frame) Used during tree building, stores potential breakpoint pseudotime between each pair of segments remaining}
#'   \item{\code{cells.in.segment}}{(List of character vectors) Cells assigned to each segment of the tree by URD}
#'   \item{\code{cells.in.nodes}}{(List of character vectors) Cells assigned to each node of the tree by URD}
#'   \item{\code{node.mean.pseudotime}}{(Named numeric vector) Mean pseudotime of cells in each node}
#'   \item{\code{node.max.pseudotime}}{(Named numeric vector) Max pseudotime of cells in each node}
#'   \item{\code{edge.list}}{(data.frame) All node-node edges in the dendrogram}
#'   \item{\code{tree.igraph}}{(igraph) igraph representation of all segment-segment relationships}
#'   \item{\code{segment.layout}}{(data.frame) }
#'   \item{\code{tree.layout}}{(data.frame) Placement of node-node edges on the dendrogram representation}
#'   \item{\code{cell.layout}}{(data.frame) Placement of cells on the dendrogram representation}
#'   \item{\code{segment.names}}{(Named vector) Names for each terminal population for use in the dendrogram layout}
#'   \item{\code{segment.names.short}}{(Named vector) Short names for each terminal population to use in the force-directed layout}
#'   \item{\code{walks.force.edges}}{(data.frame) k-NN edge list based on visitation frequency used for force-directed layout}
#'   \item{\code{walks.force.layout}}{(data.frame) 3D coordinates for the force-directed layout}
#'   \item{\code{force.view.list}}{(List) Stored orientations for displaying the force-directed layout}
#'   \item{\code{force.view.default}}{(Character) The force-directed view that should be used if none is specified}
#' }
#' @slot nmf.g (Sparse or full matrix) Non-negative matrix factorization Gene x Module matrix
#' @slot nmf.c (Sparse of full matrix) Non-negative matrix factorization Module x Cell matrix
#' 
#' @importClassesFrom destiny DiffusionMap
#' @importClassesFrom Matrix dgCMatrix
#' 
#' @exportClass URD
#' 
#' @aliases URDclass
#' @name URDclass
URD <- methods::setClass("URD", slots = c(
  count.data=c("dgCMatrix", NULL), 
  logupx.data=c("dgCMatrix", NULL), 
  meta="data.frame", 
  group.ids="data.frame", 
  var.genes="vector", 
  knn="list",
  pca.sdev="vector", 
  pca.load="data.frame", 
  pca.scores="data.frame", 
  pca.sig="vector", 
  tsne.y="data.frame", 
  plot.3d="list",
  gene.sig.z="data.frame", 
  dm=c("DiffusionMap",NULL), 
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
#' @keywords internal
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
#' @importFrom methods new
#' @importFrom Matrix rowSums colSums
#' 
#' @param count.data (Matrix or dgCMatrix) UMI expression data, with rows as genes and columns as cells
#' @param meta (data.frame) Metadata, with rows as cells (row names should match column names of \code{count.data})
#' @param min.cells (Numeric) Minimum number of cells that must express a gene to retain it
#' @param min.genes (Numeric) Minimum number of genes a cell must express to retain it
#' @param min.counts (Numeric) Minimum number of UMIs detected for a gene across the entire data to retain it
#' @param gene.max.cut (Numeric) Maximum number of UMIs observed for a gene in a single cell to retain it. If \code{NA}, this step will be skipped.
#' @param max.genes.in.ram (Numeric) Number of genes to normalize and log-transform at a time (to prevent running out of memory as matrices become non-sparse during the process)
#' 
#' @return An URD object
#' 
#' @export
createURD <- function(count.data, meta=NULL, min.cells=3, min.genes=500, min.counts=10, gene.max.cut=5000, max.genes.in.ram=5000, verbose=T) {
  # Ensuring matrices are in sparse format
  if (class(count.data) != "dgCMatrix") {
    if (verbose) message(paste0(Sys.time(), ": Converting counts to sparse matrix."))
    count.data <- as(as.matrix(count.data), "dgCMatrix")
  }
  if (verbose) message(paste0(Sys.time(), ": Creating binary count matrix."))
  count.data.binary <- count.data > 0
  # Filter cells by nGenes
  if (verbose) message(paste0(Sys.time(), ": Filtering cells by number of genes."))
  num.genes <- Matrix::colSums(count.data.binary)
  names(num.genes) <- colnames(count.data.binary)
  cells.enough.genes <- names(num.genes[which(num.genes>min.genes)])
  shhhh <- gc()
  # Filter genes by nCells
  if (verbose) message(paste0(Sys.time(), ": Filtering genes by number of cells."))
  count.data.binary <- count.data.binary[,cells.enough.genes]
  num.cells <- Matrix::rowSums(count.data.binary)
  genes.enough.cells <- names(num.cells[which(num.cells>min.cells)])
  shhhh <- gc()
  # Filter genes by total counts
  if (verbose) message(paste0(Sys.time(), ": Filtering genes by number of counts across entire data."))
  num.counts <- Matrix::rowSums(count.data[,cells.enough.genes])
  genes.enough.counts <- names(num.counts[which(num.counts>min.counts)])
  shhhh <- gc()
  # Filter genes by maximum expression
  if (!is.na(gene.max.cut)) {
    if (verbose) message(paste0(Sys.time(), ": Filtering genes by maximum observed expression."))
    count.data.binary <- count.data[, cells.enough.genes] > gene.max.cut
    gene.over.max <- names(which(Matrix::rowSums(count.data.binary) > 0))
    shhhh <- gc()
  } else {
    gene.over.max <- c()
  }
  # Figure out genes to retain
  genes.use <- setdiff(intersect(genes.enough.cells, genes.enough.counts), gene.over.max)
  
  # Create an URD object
  if (verbose) message(paste0(Sys.time(), ": Creating URD object."))
  object <- methods::new("URD", count.data=as(count.data[genes.use, cells.enough.genes], "dgCMatrix"))
  shhhh <- gc()
  
  # Determine normalization factor
  if (verbose) message(paste0(Sys.time(), ": Determining normalization factors."))
  cs <- Matrix::colSums(object@count.data)
  norm_factors <- (10**ceiling(log10(median(cs))))/cs
  shhhh <- gc()
  
  # Log-normalize the data
  if (verbose) message(paste0(Sys.time(), ": Normalizing and log-transforming the data."))
  n.chunks <- ceiling(ncol(object@count.data) / max.genes.in.ram)
  n.cells <- ncol(object@count.data)
  lognorm.chunks <- lapply(1:n.chunks, function(chunk) {
    if (verbose) message(paste0(Sys.time(), ":   Chunk ", chunk, " of ", n.chunks))
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
  
  # Set up the metadata
  object@meta <- data.frame(n.Genes=num.genes[colnames(object@count.data)])
  object@meta[,"n.Trans"] <- apply(object@count.data, 2, sum)
  if (!is.null(meta)) {
    object@meta <- cbind(object@meta, meta[rownames(object@meta),])
  }
  
  if (verbose) message(paste0(Sys.time(), ": All done."))
  return(object)
} 
