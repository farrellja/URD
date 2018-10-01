#' Parent of segment
#' 
#' @export
segParent <- function(object, segment, return.self.if.na=F, original.joins=F) {
  # Get correct segment.joins object out
  if (original.joins) sj <- object@tree$segment.joins.initial else sj <- object@tree$segment.joins
  parent <- sj[which(sj$child == segment),"parent"]
  if (!return.self.if.na) return(parent)
  if (length(parent) > 0) return(parent) else return(segment)
}

#' All Upstream Segments
#' 
#' @export
segParentAll <- function(object, segment, stop.at=NULL) {
  parents <- c()
  new.parent <- segment
  while (length(new.parent) > 0) {
    new.parent <- object@tree$segment.joins[which(object@tree$segment.joins$child %in% new.parent),"parent"]
    if (!is.null(stop.at) && new.parent %in% stop.at) break
    parents <- c(parents, new.parent)
  }
  return(parents)
}

#' Direct children of Segment
#' 
#' @export
segChildren <- function(object, segment) {
  object@tree$segment.joins[which(object@tree$segment.joins$parent %in% segment),"child"]
}

#' All children of segment
#' 
#' @export
segChildrenAll <- function(object, segment, include.self=F, original.joins=F, point.at.segment.joins.initial=original.joins) {
  children <- c()
  new.children <- segment
  while(length(new.children) > 0) {
    if (!original.joins) {
      new.children <- object@tree$segment.joins[which(object@tree$segment.joins$parent %in% new.children), "child"]
    } else {
      if (point.at.segment.joins.initial) {
        new.children <- unlist(object@tree$segment.joins.initial[which(object@tree$segment.joins.initial$parent %in% new.children), c("child.1", "child.2")])
      } else {
        new.children <- unlist(object@tree$segment.joins[which(object@tree$segment.joins$parent %in% new.children), c("child.1", "child.2")])
      }
    }
    children <- c(children, new.children)
  }
  if (include.self) children <- c(children, segment)
  children <- unique(children)
  return(children)
}

#' All children of segment
#' 
#' Returns the segment ids of all segments that are children (and grandchildren and great-grandchildren...)
#' or a parent segment specified.
#' 
#' @param object An URD object
#' @param segment (Character) Name of segment to start from
#' @param include.self (Logical) Should \code{segment} be included in the list of children?
#' @param original.joins (Logical) Should children be drawn from \code{@@tree$segment.joins} (\code{FALSE}) or
#' \code{@@tree$segment.joins.initial} (\code{TRUE})? (\code{@@tree$segment.joins.initial} reflects the process
#' of constructing the tree and includes all initial starting tips. \code{@@tree$segment.joins} reflects the
#' final structure of the tree, after segments have been removed due to insufficient members or pseudotime length.)
#' @param format (Character: "unary" or "binary") Reflects the storage structure of \code{@@tree$segment.joins}. \code{binary}
#' is only used internally, and all user queries should be \code{unary}.
#'  
#' @export
segChildrenAll <- function(object, segment, include.self=F, original.joins=F, format=c("unary", "binary")) {
  # Get desired segment.joins or segment.joins.initial out
  if (original.joins) sj <- object@tree$segment.joins.initial else sj <- object@tree$segment.joins
  # Initialize
  children <- c()
  new.children <- segment
  # Keep going down the tree until everything has been searched
  while(length(new.children) > 0) {
    if (format=="unary") {
      new.children <- sj[which(sj$parent %in% new.children), "child"]
    } else if (format=="binary") {
      new.children <- unlist(sj[which(sj$parent %in% new.children), c("child.1", "child.2")])
    } else stop("format must be 'unary' or 'binary'.")
    children <- c(children, new.children)
  }
  if (include.self) children <- c(children, segment)
  children <- unique(children)
  return(children)
}


#' Siblings of segments 
#' 
#' @export
segSiblings <- function(object, segment, include.self=T) {
  parent <- segParent(object, segment)
  sibs <- object@tree$segment.joins[which(object@tree$segment.joins$parent == parent), "child"]
  if (!include.self) sibs <- setdiff(sibs, segment)
  return(sibs)
}

# If segments have been named, can provide names to get segment numbers
#' Translate Segment Names
#' 
#' Allows access to segments by name
#' 
#' @importFrom plyr mapvalues
#' 
#' @param object An URD object
#' @param segments (Character vector) A vector of segment names for which to return the segment numbers.
#' 
#' @return Character vector of segment numbers.
#' 
#' @export
translateSegmentNames <- function(object, segments) {
  segments <- as.character(segments)
  segments <- plyr::mapvalues(segments, from=object@tree$segment.names, to=names(object@tree$segment.names), warn_missing = F)
  if(!all(segments %in% object@tree$segments)) warning("Some segment names were not found.")
  return(segments)
}

#' Extract cells along a lineage pathway
#' 
#' @param object An URD object
#' @param segments (Character vector) Terminal segment numbers or names
#' @param remove.root (Logical) 
#' 
#' @return A character vector of cell names
#' 
#' @export
cellsAlongLineage <- function(object, segments, remove.root=T) {
  segments <- translateSegmentNames(object, segments)
  segments.to.grab <- segments
  while (length(segments > 0)) {
    segments <- unique(object@tree$segment.joins[which(object@tree$segment.joins$child %in% segments),"parent"])
    segments.to.grab <- c(segments.to.grab, segments)
  }
  if (remove.root) segments.to.grab <- head(segments.to.grab, -1)
  return(unique(unlist(object@tree$cells.in.segment[segments.to.grab])))
}

#' Get all terminal segments from a tree
#' 
#' Identify the segments that do not have any descendants in URD's determined
#' tree structure
#' 
#' @param object An URD object
#' 
#' @return A character vector of segment names
#' 
#' @export
segTerminal <- function(object) {
  sort(unique(setdiff(object@tree$segment.joins$child, object@tree$segment.joins$parent)))
}