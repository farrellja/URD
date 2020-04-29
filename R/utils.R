#' Corner of data frame or matrix
#' 
#' Returns the upper-left corner of a data.frame or matrix
#' @param data.frame Input data.frame or matrix
#' @param n (Numeric) Number of rows and columns to return (default: 5)
#' @export
corner <- function(data.frame, n=5) {
  width=n
  height=n
  if (width > dim(data.frame)[2]) width <- dim(data.frame)[2]
  if (height > dim(data.frame)[1]) height <- dim(data.frame)[1]
  data.frame[1:height,1:width]
}

#' Cut string into fields
#'
#' Reworks \code{\link{strsplit}} to act like Unix cut. Instead of returning a list, returns a vector of the chosen field.
#' @param x (Character) A character vector
#' @param delimiter (Character) Delimiter character to use to split
#' @param field (Numeric) Which field to return; count starts at left with 1.
#' @return Character vector
#' @export
cutstring <- function(x, delimiter, field) {
  unlist(lapply(strsplit(x, split = delimiter), function(y) y[field]))
}

#' Arrange multiple ggplots into several gridded pages
#' 
#' Arranges several plots into a multi-page gridded graph
#' @param grobs (List) List of grobs (grid objects, e.g. returned by ggplot functions)
#' @param ncol (Numeric) Number of columns in each plot (default: 2)
#' @param nrow (Numeric) Number of rows in each plot (default: 2)
#' @importFrom gridExtra grid.arrange
#' @export
gridArrangeMulti <- function(grobs, ncol=2, nrow=2) {
  # Figure out how to do the layout
  elements.per.page <- ncol * nrow
  num.pages <- ceiling(length(grobs) / elements.per.page)
  # Loop through pages and run grid.arrange
  for (this.page in 1:num.pages) {
    min.element <- (this.page - 1) * elements.per.page + 1
    max.element <- min(c(this.page * elements.per.page, length(grobs)))
    gridExtra::grid.arrange(grobs=grobs[min.element:max.element], ncol=ncol, nrow=nrow)
  }
}

#' Check whether input is a whole number
#' 
#' @param x (Numeric vector) Input value(s)
#' @export
is.wholenumber <- function(x) { 
  return(ceiling(x) == floor(x))
}

#' Equivalent of apply function for sparse matrices from Matrix
#' 
#' Matrix package will 'automagically' and silently convert dgCMatrix sparse matrices to a
#' dense matrix when they are passed to the apply function. While rowMeans, colMeans, rowSums,
#' and colSums functions exist to reduce need for apply function, sometimes other functions
#' are needed. For loops are remarkably inefficient, so this will chunk the original matrix into
#' pieces, process them using apply, convert the results back to sparse if needed, then assemble
#' into a final result and return. It's a mediocre balance between memory usage and speed.
#' 
#' @param x (dgCMatrix) Input
#' @param margin (Numeric) \code{1} for rows, \code{2} for columns
#' @param fun (Function)
#' @param max.dim (Numeric) Number of rows or columns to process at a time
#' @param verbose (Logical) Report progress reports on chunks
#' @param ... Additional arguments to pass to \code{fun}
#' 
#' @keywords internal
#' @export
sparseApply <- function(x, margin, fun, max.dim, verbose=F, ...) {
  if (class(x) != "dgCMatrix") stop ("Designed to be used with sparse matrices.")
  # Calculate sizes of chunks
  n.chunks <- ceiling(dim(x)[margin] / max.dim)
  dim.chunk <- dim(x)[margin]
  # If only one chunk, then just use apply even though it's inefficient
  if (n.chunks == 1) return(apply(x, margin, fun, ...))
  if (margin == 1) {
    # Cut data into chunks and process individually
    done.chunks <- lapply(1:n.chunks, function(chunk) {
      if (verbose) message(paste0(Sys.time(), ":   Chunk ", chunk, " of ", n.chunks))
      shhh <- gc()
      i <- (chunk-1) * max.dim + 1
      j <- min((chunk * max.dim), dim.chunk)
      k <- apply(x[i:j,], MARGIN = margin, FUN = fun, ...)
      # Make sure matrices stay sparse
      if (class(k) == "matrix") k <- as(k, "dgCMatrix")
      return(k)
    })
    # Put results together and return
    if (class(done.chunks[[1]]) == "dgCMatrix") {
      return(do.call("cbind", done.chunks))
    } else {
      return(unlist(done.chunks))
    }
  }
  if (margin == 2) {
    # Cut data into chunks and process individually
    done.chunks <- lapply(1:n.chunks, function(chunk) {
      if (verbose) message(paste0(Sys.time(), ":   Chunk ", chunk, " of ", n.chunks))
      shhh <- gc()
      i <- (chunk-1) * max.dim + 1
      j <- min((chunk * max.dim), dim.chunk)
      k <- apply(x[,i:j], MARGIN = margin, FUN = fun, ...)
      # Make sure matrices stay sparse
      if (class(k) == "matrix") k <- as(k, "dgCMatrix")
      return(k)
    })
    # Put results together and return
    if (class(done.chunks[[1]]) == "dgCMatrix") {
      return(do.call("rbind", done.chunks))
    } else {
      return(unlist(done.chunks))
    }
  }
}

#' Equivalent of sweep function for sparse matrices from Matrix
#' 
#' Matrix package will 'automagically' and silently convert dgCMatrix sparse matrices to a
#' dense matrix when they are passed to the apply function. For loops are remarkably inefficient, 
#' so this will chunk the original matrix into
#' pieces, process them using sweep, convert the results back to sparse, then assemble
#' into a final result and return. It's a mediocre balance between memory usage and speed.
#' 
#' @param x (dgCMatrix) Input
#' @param margin (Numeric) \code{1} for rows, \code{2} for columns
#' @param stats (Vector) Summary statistic to be swept out
#' @param fun (Function)
#' @param max.dim (Numeric) Number of rows or columns to process at a time
#' @param verbose (Logical) Report progress reports on chunks
#' @param ... Additional arguments to pass to \code{fun}
#' 
#' @keywords internal
#' @export
sparseSweep <- function(x, margin, stats, fun, max.dim, verbose=F, ...) {
  if (class(x) != "dgCMatrix") stop ("Designed to be used with sparse matrices.")
  # Calculate sizes of chunks
  n.chunks <- ceiling(dim(x)[margin] / max.dim)
  dim.chunk <- dim(x)[margin]
  # If only one chunk, then just use sweep
  if (n.chunks == 1) return(sweep(x, margin, stats, fun, ...))
  if (margin == 1) {
    # Cut data into chunks and process individually
    done.chunks <- lapply(1:n.chunks, function(chunk) {
      if (verbose) message(paste0(Sys.time(), ":   Chunk ", chunk, " of ", n.chunks))
      shhh <- gc()
      i <- (chunk-1) * max.dim + 1
      j <- min((chunk * max.dim), dim.chunk)
      k <- sweep(x[i:j,], MARGIN = margin, STATS = stats[i:j], FUN = fun, ...)
      # Make sure matrices stay sparse
      k <- as(k, "dgCMatrix")
      return(k)
    })
    # Put results together and return
    return(do.call("rbind", done.chunks))
  }
  if (margin == 2) {
    # Cut data into chunks and process individually
    done.chunks <- lapply(1:n.chunks, function(chunk) {
      if (verbose) message(paste0(Sys.time(), ":   Chunk ", chunk, " of ", n.chunks))
      shhh <- gc()
      i <- (chunk-1) * max.dim + 1
      j <- min((chunk * max.dim), dim.chunk)
      k <- sweep(x[,i:j], MARGIN = margin, STATS = stats[i:j], FUN = fun, ...)
      # Make sure matrices stay sparse
      k <- as(k, "dgCMatrix")
      return(k)
    })
    # Put results together and return
    return(do.call("cbind", done.chunks))
  }
}