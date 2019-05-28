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

