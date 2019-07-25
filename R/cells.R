#' Which cells
#' 
#' Identifies cells that meet particular criteria and returns their cell names. This
#' can be used to select populations for subsetting, differential expression, plotting, etc.
#' 
#' @param object An URD object
#' @param label (Character) The label of the data to search for
#' @param value (Character or numeric vector) If \code{label} is discrete (e.g. a clustering, developmental stage, treatment), one or more character values to match against. If \code{label} is continuous (e.g. pseudotime), should be a numeric vector or length 2 that describes the inclusive range of acceptable values.
#' @param label.type (Character) Where to look for the data. Default is "search" which looks in order: "meta", "group", "sig", "gene", "counts", "pseudotime", "pca", "diff.data"
#' 
#' @return (Character vector) Cell  that match the criteria specified.
#' 
#' @examples
#' # Find cells from dome stage
#' whichCells(urd.object, "STAGE", "ZFDOME")
#' 
#' # Find cells with 500-1000 genes detected
#' whichCells(urd.object, "NUM_GENES", c(500,1000))
#' 
#' # Find cells with pseudotime between 0.54-0.72
#' whichCells(urd.object, "pseudotime", c(0.54,0.72))
#' 
#' @export
whichCells <- function(object, label, value, label.type="search") {
  data <- data.for.plot(object, label=label, label.type=label.type, as.discrete.list=T)
  if (data$discrete) {
    return(names(data$data)[which(data$data %in% value)])
  } else {
    if (length(value) != 2) {
      stop("Label ", label, " is not discrete, so value should be length 2 and describe a range of acceptable values.")
    }
    return(names(data$data)[which(data$data >= value[1] & data$data <= value[2])])
  }
}