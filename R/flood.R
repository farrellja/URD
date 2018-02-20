#' Build pseudotime transition matrix
#' 
#' This scales the transition probabilities so that maximum sum of transition
#' probabilities for a cell is 1. 
#' This modified matrix is used for calculating pseudotime. It requires
#' either the user to provide an URD object (\code{object=}) with a diffusion map
#' stored in slot \code{@@dm}, or a diffusion map directly (\code{dm=}).
#' 
#' This function can be called directly and the transition matrix stored, or
#' \code{\link{floodPseudotime}} will call this function automatically if
#' one isn't provided.
#' 
#' @param object An URD object
#' @param dm A Diffusion Map
#' @return A sparse matrix (dgCMatrix) of normalized transition probabilities.
#' @export
floodBuildTM <- function(object=NULL, dm=NULL) {
  if (is.null(object) & is.null(dm)) stop("Must provide either object or diffusion map.")
  # Normalize transition matrix so that maximum sum of transition probabilities is 1.
  if (!is.null(object)) { tm.flood <- object@dm@transitions }
  else if (!is.null(dm)) { tm.flood <- dm@transitions }
  tm.sum.max <- max(apply(tm.flood, 1, sum))
  tm.flood <- tm.flood / tm.sum.max
  return(tm.flood)
}

#' Combine probabilities
#' 
#' Calculates cumulative probability
#' @param x (Numeric Vector) Probabilities
#' @return (Numeric) Cumulative probability
combine.probs <- function(x) {
  1-prod(1-x)
}

#' A single iteration of flood pseudotime calculation.
#' 
#' This runs a single iteration of flood pseudotime. This function is meant to 
#' be called by \code{\link{floodPseudotime}} to run many iterations.
#' 
#' @param object An URD object
#' @param start.cells (Character vector) Cells to initialize as visited to start pseudotime calculation
#' @param minimum.cells.flooded (Numeric) Stop pseudotime calculation when fewer cells are newly visited in a given iteration, and assign remaining unvisited cells pseudotime NA.
#' @param tm.flood (Sparse matrix) A sparse matrix of normalized transition probabilities. If unprovided (\code{NULL}), this will be calculated automatically.
#' @param verbose.freq (Numeric) Report progress after this many iterations. 0 is silent.
#' 
#' @return (Numeric vector) The iteration that newly visited each cell.
floodPseudotimeCalc <- function(object, start.cells, minimum.cells.flooded=2, tm.flood=NULL, verbose.freq=0) {
  # Build a normalized transition matrix for flooding, unless it's passed to the function.
  if (is.null(tm.flood)) {
    if (verbose) message(paste0(Sys.time(), ": tm.flood not specified, so calculating normalized transition matrix."))
    tm.flood <- floodBuildTM(object)
  }
  
  # Start vector to store which cells have been visited and their pseudotime
  pseudotime <- rep(NA, dim(tm.flood)[1])
  names(pseudotime) <- rownames(tm.flood)
  i <- 0
  
  # Initialize with starting cells
  cells.visited <- start.cells
  cells.not.visited <- setdiff(names(pseudotime), cells.visited)
  pseudotime[start.cells] <- 0
  total.cells <- dim(tm.flood)[1]
  cells.newly.visited <- rep(NA, minimum.cells.flooded)
  
  # Loop until all cells have a pseudotime
  while (length(cells.visited) < (total.cells-1) & length(cells.newly.visited) >= minimum.cells.flooded) {
    # Increment step counter
    i <- i + 1
    if (i %% verbose.freq == 0) print(paste("Flooding, step", i, "-", paste0(round(length(cells.visited) / total.cells * 100, digits=1), "%:"), length(cells.visited), "of", total.cells, "cells visited."))
    # Calculate visitation probability for each cell
    visitation.prob <- apply(tm.flood[cells.not.visited,cells.visited], 1, combine.probs)
    # Figure out which cells are newly visited
    cells.newly.visited <- names(visitation.prob)[which(rbinom(length(visitation.prob), 1, visitation.prob) == 1)]
    # Record the visited cells
    cells.visited <- c(cells.visited, cells.newly.visited)
    cells.not.visited <- setdiff(cells.not.visited, cells.newly.visited)
    pseudotime[cells.newly.visited] <- i
  }
  
  return(pseudotime)
}

#' Calculate pseudotime by 'flooding'
#' 
#' This calculates pseudotime by performing a probabilistic breadth-first search
#' of the k-nearest neighbor graph that was used to generate the diffusion map
#' calculated on the data. 
#' 
#' A group of 'root' cells are pre-initialized as visited, an an iterative search
#' of the graph is performed. The probability that a cell is newly visited is the
#' cumulative transition probability from all cells that have previously been
#' visited in the graph. The process continues until either the entire graph is
#' visited, or fewer than \code{minimum.cells.flooded} new cells are visited in
#' a given iteration. 
#' 
#' @param object An URD object
#' @param root.cells (Character vector) Names of cells that will be considered the root and assigned pseudotime 0.
#' @param n (Numeric) Number of simulations to perform and average
#' @param minimum.cells.flooded (Numeric) Stop pseudotime calculation when fewer cells are newly visited in a given iteration, and assign remaining unvisited cells pseudotime NA. This is designed to prevent a few poor data points from throwing off the entire pseudotime calculation.
#' @param tm.flood (Sparse matrix) A sparse matrix of normalized transition probabilities (as returned by \code{\link{floodBuildTM}}. If unprovided (\code{NULL}), this will be calculated automatically.
#' @param verbose (Logical) Report on progress?
#' @param verbose.freq (Numeric) Give a report every this many iterations.
#' 
#' @return (data.frame) Rows are cells, columns are flood simulations, values are the iteration that newly visited each cell. This can be processed into pseudotime and stored in an URD object with \code{\link{floodPseudotimeProcess}}.
#' @export
floodPseudotime <- function(object, root.cells, n=20, minimum.cells.flooded=2, tm.flood=NULL, verbose=F, verbose.freq=10) {
  if (is.null(tm.flood)) tm.flood <- floodBuildTM(object)
  if (!verbose) verbose.freq <- 0
  floods <- as.data.frame(lapply(1:n, function(x) {
    if (verbose) print(paste("Starting flood number", x))
    floodPseudotimeCalc(object, start.cells = root.cells, minimum.cells.flooded = minimum.cells.flooded, tm.flood = tm.flood, verbose.freq=verbose.freq)
  }))
  colnames(floods) <- 1:n
  return(floods)
}



#' Process Flood Pseudotime
#' 
#' This processes the values returned by \code{\link{floodPseudotime}} and stores
#' the result as a pseudotime in an URD object.
#' 
#' @param object An URD object
#' @param floods (data.frame) A data.frame of flood visitation information, returned by \code{\link{floodPseudotime}}
#' @param floods.name (Character) The name to use for this pseudotime.
#' @param max.frac.NA (Numeric) Maximum fraction of simulations that assigned NA for a cell. Cells with more NAs will not be assigned a pseudotime.
#' @param pseudotime.fun (Function) Function to combine pseudotimes resulting from each simulation. Default (and only validated function) is mean.
#' @param stability.div (Numeric) Number of simulation subsamplings to calculate.
#' @return An URD object with pseudotime stored in \code{@@pseudotime} under the column name given in \code{floods.name}, with subsampling information stored in \code{@@pseudotime.stability}. The information stored in \code{@@pseudotime.stability} will be overwritten by any future pseudotime processing.
#' @export
floodPseudotimeProcess <- function(object, floods, floods.name="flood", max.frac.NA=0.4, pseudotime.fun=mean, stability.div=10) {
  # Make sure that for tips this doesn't accidentally refer to a columns index.
  floods.name <- as.character(floods.name)
  # If floods is a list of data.frames, combine them into a single data.frame
  if (class(floods) == "list") {
    floods <- do.call("cbind", floods)
  } else if (class(floods) != "data.frame") {
    stop("floods should be a list or data.frame.")
  }
  # Remove any cells with too many NAs as outliers
  frac.na <- apply(floods, 1, function(x) sum(is.na(x))) / dim(floods)[2]
  cells.toss <- names(which(frac.na > max.frac.NA))
  floods <- floods[setdiff(rownames(floods), cells.toss),]
  # Turn flood positions into relative positions
  # (e.g. instead of ordinal positioning, fractional position from 0 to 1)
  floods <- sweep(floods, 2, apply(floods, 2, max, na.rm=T), "/")
  # Divide the floods into sections for stability calculations
  floods.in.division <- ceiling(1:stability.div / stability.div * dim(floods)[2])
  ## CALCULATE PSEUDOTIME AND VISIT FREQUENCY WITH INCREASING AMOUNTS OF DATA TO TEST WHEN IT STABILIZES
  walks.per.cell <- as.data.frame(lapply(floods.in.division, function(n.floods) {
    if (n.floods == 1) return(as.numeric(!is.na(floods[,1])))
    else return (apply(floods[,1:n.floods], 1, function(x) sum(as.numeric(!is.na(x)))))
  }))
  pseudotime.stability <- as.data.frame(lapply(floods.in.division, function(n.floods) {
    if (n.floods == 1) return(floods[,1])
    else return (apply(floods[,1:n.floods], 1, mean, na.rm=T))
  }))
  colnames(walks.per.cell) <- seq(1, ncol(floods), length.out=stability.div)
  colnames(pseudotime.stability) <- seq(1, ncol(floods), length.out=stability.div)
  # Store the results: matrices for stability plots
  object@pseudotime.stability$pseudotime <- pseudotime.stability
  object@pseudotime.stability$walks.per.cell <- walks.per.cell
  # Store the results: after calculating with all walks, visit frequencies in diff.data
  final.visit.freq <- walks.per.cell[rownames(object@diff.data),dim(walks.per.cell)[2]]
  #final.visit.freq[is.na(final.visit.freq)] <- 0
  object@diff.data[,paste0("visitfreq.raw.", floods.name)] <- final.visit.freq
  object@diff.data[,paste0("visitfreq.log.", floods.name)] <- log10(final.visit.freq + 1)
  # Store the results: pseudotime
  object@pseudotime[,floods.name] <- NA
  object@pseudotime[rownames(pseudotime.stability), floods.name] <- pseudotime.stability[,stability.div]
  # Return the object
  return(object)
}

#' Find tip potential
#' 
#' This finds the neighbors in gene expression space within a given range, and calculates the proportion
#' of those neighbors that have earlier pseudotimes.
#' 
#' In samples that were not a time-course, the length of branches in pseudotime won't be constant, but will
#' instead represent the relative differences in gene expression from the progentior stem cells to the
#' final differentiated fates. Thus, rather than using pseudotime directly, this calculates whether cells'
#' neighborhoods are asymmetric, such that most cells in their neighborhood have earlier pseudotimes,
#' which indicates that they have a high likelihood of being tips.
#' 
#' In large data sets, it can be slow to calculate the distance between cells, 
#' and may be worth calculating in advance with \link{\code{cellDistByExpression}}
#' so that multiple parameters for \code{nn.dist.range} can be explored without
#' needing to recalculate.
#' 
#' @param object An URD object
#' @param pseudotime (Character) The column name of a previously calculated pseudotime
#' @param range.function (Function) The function used to summarize the delta pseudotimes for each cell (default is the proportion <= 0, but mean and other options can also work here)
#' @param genes.use (Character vector) Which genes to use for calculating the distance between cells. (Default is variable genes)
#' @param nn.dist.range (Numeric vector, length 2) The range of nearest neighbor distances to consider.
#' @param pseudotime.name (Character) Output as saved in \code{@@pseudotime} under this column name.
#' @param cell.dist (Matrix) Distance between cells (see \code{\link{cellDistByExpression}}; if NULL, will be calculated automatically.)
#' @return An URD object with \code{@@pseudotime$pseudotime.name} set to the potential for cells to be part of a tip. 
#' @export
tipPotential <- function(object, pseudotime, range.function=prop.nonexp, genes.use=object@var.genes, nn.dist.range=c(20,22), pseudotime.name="tip.potential", cell.dist=NULL) {
  # Get distance between all cells
  if (is.null(cell.dist)) cell.dist <- cellDistByExpression(object, genes.use=genes.use, return.as.matrix=T)
  # Find the neighbors for each cell in a defined distance range
  neighbors.in.range <- lapply(rownames(cell.dist), function(cell) colnames(cell.dist)[which(cell.dist[cell,] >= nn.dist.range[1] & cell.dist[cell,] <= nn.dist.range[2])])
  names(neighbors.in.range) <- rownames(cell.dist)
  # Find the relative pseudotime for those neighbors
  neighbor.rel.pt <- lapply(rownames(cell.dist), function(cell) {
    object@pseudotime[neighbors.in.range[[cell]], pseudotime] - object@pseudotime[cell,pseudotime]
  })
  names(neighbor.rel.pt) <- rownames(cell.dist)
  # Summarize delta pseudotimes 
  neighbor.rel.pt.processed <- unlist(lapply(neighbor.rel.pt, range.function))
  object@pseudotime[,pseudotime.name] <- neighbor.rel.pt.processed[rownames(object@pseudotime)]
  return(object)
}

#' Calculates the distance matrix between cells given their gene expression
#' 
#' Calculates distance between cells in gene expression space. Can be passed to \code{\link{tipPotential}}.
#' 
#' @importFrom stats dist
#' 
#' @param object An URD object
#' @param genes.use (Character vector) Which genes to use for distance calculation (default: variable genes. NULL is all genes)
#' @param return.as.matrix (Logical) Return value as a matrix (memory hungry, but used for other URD functions) or as a dist object.
#' 
#' @return Returns the distance between cells in gene expression space, either as a matrix
#' (if \code{return.as.matrix=T}) or a dist object.
#' 
#' @export
cellDistByExpression <- function(object, genes.use=object@var.genes, return.as.matrix=T) {
  # Get all genes if genes.use not specified.
  if (is.null(genes.use)) genes.use <- rownames(object@logupx.data)
  # Get distance between all cells
  d <- dist(t(object@logupx.data[genes.use,]))
  if (return.as.matrix) return(as.matrix(d)) else return(d)
}
