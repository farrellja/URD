#' Create a group.id from a list of cells
#' 
#' This function will create a column in \code{group.ids} that identifies
#' whether cells are part of a provided list of cells. This can be used to
#' highlight a group of cells in any of URD's plotting functions.
#' 
#' @param object An URD object
#' @param group.id (Character) Name of the \code{group.id} to store.
#' @param cells (Character vector) Name of cells in the group.
#' 
#' @return Returns an URD object, with a column in \code{@@group.ids}
#' marked \code{TRUE} for cells in \code{cells} and \code{FALSE} for
#' all other cells.
#' 
#' @export
groupFromCells <- function(object, group.id, cells) {
  object@group.ids[,group.id] <- rownames(object@group.ids) %in% cells
  return(object)
}