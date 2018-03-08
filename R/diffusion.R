#' Determine logistic parameters
#' 
#' @param object An URD object
#' @param pseudotime (Character)
#' @param optimal.cells.forward (Numeric)
#' @param max.cells.back (Numeric)
#' @param pseudotime.direction (Character: "<" or ">")
#' @param asymptote (Numeric)
#' @param do.plot (Logical)
#' @param print.values (Logical)
#' @return Logistic parameters
#' @export

pseudotimeDetermineLogistic <- function(object, pseudotime, optimal.cells.forward, max.cells.back, pseudotime.direction=">", asymptote=0.01, do.plot=T, print.values=T) {
  if (pseudotime.direction == ">") {
    sort.dec <- FALSE
  } else if (pseudotime.direction == "<") {
    sort.dec <- TRUE
  } else {
    stop("pseudotime.direction must be either \">\" or \"<\"")
  }
  pseudotime.vec <- sort(object@pseudotime[,pseudotime], decreasing=sort.dec)
  mean.pseudotime.back <- mean(pseudotime.vec[1:(length(pseudotime.vec) - max.cells.back)] - pseudotime.vec[(max.cells.back+1):length(pseudotime.vec)])
  mean.pseudotime.forward <- mean(pseudotime.vec[(optimal.cells.forward+1):length(pseudotime.vec)] - pseudotime.vec[1:(length(pseudotime.vec)-optimal.cells.forward)])
  x0 <- mean(c(mean.pseudotime.back, mean.pseudotime.forward))  
  k <- log(asymptote) / (x0 - mean.pseudotime.forward)
  if (do.plot) {
    x <- seq(2*mean.pseudotime.back, 2*mean.pseudotime.forward, length.out = 100)
    plot(x=x, y=logistic(x, x0, k), pch=20, xlab="Delta pseudotime", ylab="Chance of accepted move")
    abline(v=0, col='red')
    abline(v=mean.pseudotime.back, col='blue', lty=2)
    abline(v=mean.pseudotime.forward, col='blue', lty=2)
  }
  if (print.values) {
    print(paste0("Mean pseudotime back (~", max.cells.back, " cells) ", mean.pseudotime.back))
    print(paste0("Chance of accepted move to equal pseudotime is ", logistic(0, x0, k)))
    print(paste0("Mean pseudotime forward (~", optimal.cells.forward, " cells) ", mean.pseudotime.forward))
  }
  return(list(x0=x0, k=k))
}

#' Weight transition matrix by pseudotime
#' 
#' @param object An URD object
#' @param pseudotime (Character)
#' @param x0 (Numeric)
#' @param k (Numeric)
#' @param logistic.params (List)
#' @param pseudotime.direction (Character: ">" or "<")
#' @return Sparse Matrix (dgCMatrix) of transition probabilities, weighted by pseudotime
#' @export
pseudotimeWeightTransitionMatrix <- function(object, pseudotime, x0=NULL, k=NULL, logistic.params=NULL, pseudotime.direction=">") {
  # Check that pseudotime.direction is valid
  if (!(pseudotime.direction %in% c(">", "<"))) stop ("pseudotime.direction parameter must be either \">\" or \"<\"")
  # Unpack logistic.params
  if (!is.null(logistic.params)) {
    if (is.null(x0)) x0 <- logistic.params$x0
    if (is.null(k)) k <- logistic.params$k
  }
  # Check that logistic parameters were specified
  if (is.null(x0) | is.null(k)) stop("Either specify x0 and k, or provide the results of determine.logistic.parameters to logistic.params.")
  # Get pseudotime vector
  cell.names <- rownames(object@dm@transitions)
  pseudotime.vec <- object@pseudotime[cell.names,pseudotime]
  # Figure out which cells have pseudotime and only use those.
  cells.w.pt <- which(!is.na(pseudotime.vec))
  cell.names <- cell.names[cells.w.pt]
  pseudotime.vec <- pseudotime.vec[cells.w.pt]
  # Calculate delta-pt matrix
  time.transition.matrix <- t(sapply(pseudotime.vec, function(y) {
    logistic(x=(pseudotime.vec - y), x0=x0, k=k) 
  }))
  time.transition.matrix <- time.transition.matrix * object@dm@transitions[cells.w.pt, cells.w.pt]
  rownames(time.transition.matrix) <- cell.names
  colnames(time.transition.matrix) <- cell.names
  return(time.transition.matrix)
}

#' Simulate random walks
#' 
#' 
#' 
#' @param start.cells (Character vector) Cells to use as a starting pool. One is chosen at random.
#' @param transition.matrix (Matrix or Sparse Matrix) Transition matrix 
#' @param end.cells (Character vector)
#' @param n (Numeric) Number of walks to simulate
#' @param end.visits (Numeric) Number of visits to end.cells to do before stopping
#' @param verbose.freq (Numeric) Print a progress update, every \code{verbose.freq} walks. If 0, no progress updates are reported.
#' @param max.steps (Numeric) Maximum number of steps to take before aborting the walk, making the assumption that the walk has somehow gotten stuck. Returns \code{NULL} for those walks.
#' 
#' @return (Character vector) Names of cells visited during the random walk.
#' 
#' @export
simulateRandomWalk <- function(start.cells, transition.matrix, end.cells, n=10000, end.visits=1, verbose.freq=0, max.steps=5000) {
  return(lapply(1:n, function(i) {
    # Initialize with a starting cell from start.cell.list
    current.cell <- sample(start.cells, 1)
    if (verbose.freq > 0 && i %% verbose.freq == 0) print(paste0("Starting walk ", i, " from cell ", current.cell))
    diffusion.path <- c(current.cell)
    stops.in.endzone <- 0
    n.steps <- 0
    # Continue moving to new cells until within root.threshold cells of youngest cell
    while(stops.in.endzone < end.visits) {
      # Grab a new cell based on the weights from the transition probability matrix
      current.cell <- sample(rownames(transition.matrix), size=1, prob=transition.matrix[current.cell,])
      # Update and store information about the new cell.
      diffusion.path <- c(diffusion.path, current.cell)
      # Count the number of cells in the root you've visited
      if (current.cell %in% end.cells) stops.in.endzone <- stops.in.endzone + 1
      # If walk is too long, then it's probably stuck, so abandon.
      n.steps <- n.steps + 1
      if (n.steps > max.steps) {
        warning("Walk ", i, " length greater than ", max.steps, " so returning NULL.\n")
        return(NULL)
      }
    }
    # Return the list of cells walked through
    return(diffusion.path)
  }))
}

#' Simulate Random Walks From All Tips
#' 
#' @param object An URD object
#' @param tip.group.id (Character) The name of the clustering that defines tips
#' @param root.cells (Character vector) Names of cells that constitute the root
#' @param transition.matrix (Matrix or dgCMatrix) Biased transition matrix 
#' @param n.per.tip (Numeric) Number of walks to do per tip
#' @param root.visits (Numeric) Number of steps to take that visit a root.cell before stopping
#' @param max.visits (Numeric) Abandon walks that take more steps than this, as it likely means that it has gotten stuck. (Default is number of cells in the data. On large data, may want to lower this value.)
#' @param verbose (Logical) Whether to report on progress
#' @return List (of tips) of lists (of walk paths) of character vectors
#' @export
simulateRandomWalksFromTips <- function(object, tip.group.id, root.cells, transition.matrix, n.per.tip=10000, root.visits=1, max.steps=ncol(object@logupx.data), verbose=T) {
  # Get all tips from that id
  all.tips <- sort(setdiff(as.character(unique(object@group.ids[,tip.group.id])), NA))
  n.tips <- length(all.tips)
  # Run through tips
  tip.walks <- lapply(all.tips, function(tip) {
    if (verbose) print(paste0(Sys.time(), " - Starting random walks from tip ", tip))
    # Get tip cells
    tip.cells <- rownames(object@group.ids)[which(as.character(object@group.ids[,tip.group.id]) == tip)]
    # Do random walks
    walks <- simulateRandomWalk(start.cells=tip.cells, transition.matrix=transition.matrix, end.cells=root.cells, n=n.per.tip, end.visits=root.visits, verbose.freq=round(n.per.tip/10), max.steps=max.steps)
    return(walks)
  })
  names(tip.walks) <- all.tips
  return(tip.walks)
}

#' Process random walks and store in object
#' 
#' 
#' @importFrom reshape2 dcast melt
#' @importFrom data.table rbindlist
#' 
#' @param object An URD object
#' @param walks (List) List of character vectors of cells visited during random walks.
#' @param walks.name (Character) Name to use for storing walks
#' @param aggregate.fun (Function) Function to aggregate pseudotime (default: mean)
#' @param n.subsample (Numeric) Number of subsamplings to perform for calculating 
#' @param verbose (Logical) Report on progress
#' @return An URD object with cell visitation frequency stored in \code{@@diff.data}, a calculated pseudotime stored in \code{@@pseudotime}, and
#' subsampled data in \code{@@pseudotime.stability}.
#' @export
processRandomWalks <- function(object, walks, walks.name, aggregate.fun=mean, n.subsample=10, verbose=T) {
  # Make sure that for tips this doesn't accidentally refer to a columns index.
  walks.name <- as.character(walks.name)
  # Remove any walks that are NULL (which means they were cut off because they had too many steps)
  walks <- walks[!unlist(lapply(walks, is.null))]
  # Calculate cells' relative position in the diffusion list
  # (e.g. instead of ordinal positioning, fractional position from 0 to 1)
  hops.melt <- as.data.frame(rbindlist(lapply(walks, function(i) {
    data.frame(hops=(1:length(i))/length(i), cell=i, stringsAsFactors=F)
  })), stringsAsFactors=F)
  # Divide the walks up into sections for stability calculations
  walks.in.division <- ceiling(1:n.subsample / n.subsample * length(walks))
  walk.lengths <- unlist(lapply(walks, length))
  ## CALCULATE PSEUDOTIME AND VISIT FREQUENCY WITH INCREASING AMOUNTS OF DATA TO TEST WHEN IT STABILIZES
  # Initialize matrices to hold the results
  pseudotime.stability <- matrix(rep(NA, dim(object@diff.data)[1] * n.subsample), nrow = dim(object@diff.data)[1], ncol=n.subsample, dimnames = list(rownames(object@diff.data), walks.in.division))
  walks.per.cell <- matrix(rep(0, dim(object@diff.data)[1] * n.subsample), nrow = dim(object@diff.data)[1], ncol=n.subsample, dimnames = list(rownames(object@diff.data), walks.in.division))
  # Loop through each section
  for (section in 1:n.subsample) {
    if (verbose) print(paste("Calculating pseudotime with", walks.in.division[section], "walks."))
    # Get number of times each cell has been visited so far
    visit.freq <- table(unlist(walks[1:walks.in.division[section]]))
    walks.per.cell[names(visit.freq),section] <- visit.freq
    # Calculate pseudotime per cell, given number of walks so far
    hops.melt.rows <- sum(walk.lengths[1:walks.in.division[section]])
    hops.relative <- dcast(data=hops.melt[1:hops.melt.rows,], formula=cell~., fun.aggregate=aggregate.fun, value.var="hops")
    pseudotime.stability[hops.relative$cell,section] <- hops.relative[,"."]
  }
  # Store the results: matrices for stability plots
  object@pseudotime.stability$pseudotime <- pseudotime.stability
  object@pseudotime.stability$walks.per.cell <- walks.per.cell
  # Store the results: after calculating with all walks, visit frequencies in diff.data
  final.visit.freq <- walks.per.cell[rownames(object@diff.data),dim(walks.per.cell)[2]]
  final.visit.freq[is.na(final.visit.freq)] <- 0
  object@diff.data[,paste0("visitfreq.raw.", walks.name)] <- final.visit.freq
  object@diff.data[,paste0("visitfreq.log.", walks.name)] <- log10(final.visit.freq + 1)
  # Store the results: pseudotime
  object@pseudotime[,walks.name] <- NA
  object@pseudotime[hops.relative[,"cell"], walks.name] <- hops.relative[,"."]
  # Return the object
  return(object)
}

#' Process random walks from all tips
#' 
#' 
#' @param object An URD object
#' @param walks.list (List of lists) List of character vectors of cells visited during random walks.
#' @param aggregate.fun (Function) Function to aggregate pseudotime (default: mean)
#' @param n.subsample (Numeric) Number of subsamplings to perform for calculating 
#' @param verbose (Logical) Report on progress
#' @return An URD object with cell visitation frequency stored in \code{@@diff.data}, a calculated pseudotime stored in \code{@@pseudotime}, and
#' subsampled data in \code{@@pseudotime.stability}.
#' @export
processRandomWalksFromTips <- function(object, walks.list, aggregate.fun=mean, n.subsample=10, verbose=T) {
  # Check that walks.list is a proper list
  if (class(walks.list) != "list" || class(walks.list[[1]]) != "list") stop("walks.list must be a list (tips) of lists (walks).")
  if (is.null(names(walks.list))) stop("walks.list must be named according to tips.")
  # Get names
  tip.names <- names(walks.list)
  # Loop through
  for (tip in tip.names) {
    if (verbose) print(paste0(Sys.time(), " - Processing walks from tip ", tip))
    # Process the walks
    object <- processRandomWalks(object, walks = walks.list[[tip]], walks.name = tip, aggregate.fun = aggregate.fun, n.subsample = n.subsample, verbose = verbose)
  }
  # Return the object
  return(object)
}