#' Load Tip Cells
#' 
#' This loads the cells that belong to each tip into the URD object. This information
#' is required if \code{weighted.fusion=T} during \code{\link{buildTree}} so that
#' new segment visitation can properly be determined using a weighted average, based
#' on the number of cells that were in each tip.
#' 
#' @param object An URD object
#' @param tips (Character) A clustering that was used as the tips (i.e. a column name of \code{@@group.ids})
#' 
#' @return An URD object with the cells that were used to start each tip stored in \code{@@tree$cells.in.tip}
#' 
#' @examples 
#' # Load the cells used for each tip into the URD object
#' axial.tree <- loadTipCells(axial, "tip.clusters")
#' 
#' @export
loadTipCells <- function(object, tips) {
  all.tips <- setdiff(unique(object@group.ids[,tips]), NA)
  cells.in.tip <- lapply(all.tips, function(tip) rownames(object@group.ids)[which(object@group.ids[,tips] == tip)])
  names(cells.in.tip) <- all.tips
  object@tree$cells.in.tip <- cells.in.tip
  object@tree$tips <- all.tips
  return(object)
}

#' Combine Tip Visitation
#' 
#' If two tips were used to do random walks, but then seem to be the exact same population,
#' the random walks can be combined (biased by the number of cells in each tip), to combine
#' their visitation frequencies. Requires \code{\link{loadTipCells}} to be run prior.
#' 
#' @param object An URD object
#' @param tip.1 (Character) Name of a tip to combine
#' @param tip.2 (Character) Name of another tip to combine
#' @param new.tip (Character) Name to call the combined tip
#' 
#' @return An URD object
#' 
#' @export
combineTipVisitation <- function(object, tip.1, tip.2, new.tip) {
  # Get number of cells in each tip
  if (is.null(object@tree$cells.in.tip)) stop("Run load.tip.cells first to properly weight means.")
  tip.n <- unlist(lapply(object@tree$cells.in.tip, length))
  # Get columns to act on
  cols <- paste0("visitfreq.raw.", c(tip.1, tip.2))
  cols.log <- paste0("visitfreq.log.", c(tip.1, tip.2))
  # Take weighted mean of visit frequency based on number of cells in each tip.
  newcol <- apply(object@diff.data[,cols], 1, weighted.mean, w=tip.n[c(tip.1,tip.2)])
  # Remove old visit frequencies from walk data
  object@diff.data <- object@diff.data[,setdiff(names(object@diff.data), c(cols, cols.log))]
  # Put new data into diff.data
  object@diff.data[,paste0("visitfreq.raw.", new.tip)] <- newcol
  object@diff.data[,paste0("visitfreq.log.", new.tip)] <- log10(newcol+1)
  # Combine tip.cell recordings
  new.tip.cells <- unique(unlist(object@tree$cells.in.tip[c(tip.1, tip.2)]))
  object@tree$cells.in.tip <- object@tree$cells.in.tip[setdiff(names(object@tree$cells.in.tip), c(tip.1, tip.2))]
  object@tree$cells.in.tip[[new.tip]] <- new.tip.cells
  # Remove tips from object@tree$tips
  object@tree$tips <- c(setdiff(object@tree$tips, c(tip.1, tip.2)), new.tip)
  # Return object
  return(object)
}

#' Name Segments
#' 
#' This stores the names of tips for use in plotting. \code{segment.names} is used on dendrogram-style plots
#' (\code{\link{plotTree}}). They are also used on force-directed layout plots \code{\link{plotTreeForce}},
#' unless \code{short.names} is set, in which case those shorter names are used. (Labels are horizontal on
#' the force-directed layout, so it is our recommendation that a series of 2-4 letter abbreviations for
#' each cell population are used.) For terminal segments of the tree that are not given a name, this function
#' will determine their children and (if named) combine the names of the segment's children. (This is in
#' case, for instance, two populations immediately fuse into a single segment.) Segment numbers to use in
#' this function are easily determined using \code{\link{plotTree}} with \code{label.segments=T}.
#' 
#' @param object An URD object
#' @param segments (Character Vector) Segment numbers to name
#' @param segment.names (Character Vector) Names to use for segment
#' @param short.names (Character Vector, optional) Short names to use when labeling tips in 3D layout
#' @param sep (Character) Separator to use for combining tip names, if needed. 
#' @export

nameSegments <- function(object, segments, segment.names, short.names=NULL, sep="+") {
  # Make sure everything is treated as a character vector, not numeric or a factor
  segments <- as.character(segments)
  segment.names <- as.character(segment.names)
  if (!is.null(short.names)) short.names <- as.character(short.names)
  
  # Get all terminal segments of the tree and original tips
  terminal <- segTerminal(object)
  
  # Determine which terminal segments are named
  original <- intersect(terminal, segments)
  
  # Associate segment numbers with their segment.names
  birth.names <- data.frame(
    segment=segments,
    segment.names=segment.names,
    row.names = segments,
    stringsAsFactors=F
  )
  if (!is.null(short.names)) birth.names$short.names <- short.names
  
  # Make a data.frame to keep track of the new names.
  name.them <- data.frame(
    segment=original,
    segment.names=birth.names[original, "segment.names"],
    stringsAsFactors=F
  )
  if (!is.null(short.names)) name.them$short.names <- birth.names[original, "short.names"]
  
  # Figure out which terminal tips are a combination of ones that exist
  combined <- setdiff(terminal, segments)
  for (tt in combined) {
    og.tips <- intersect(segChildrenAll(object, segment = tt, original.joins = T), segments)
    if (!is.null(short.names)) {
      new.name <- c(
        tt, 
        paste0(sort(birth.names[og.tips, "segment.names"]), collapse=sep),
        paste0(sort(birth.names[og.tips, "short.names"]), collapse=sep)
      ) 
    } else {
      new.name <- c(
        tt, 
        paste0(sort(birth.names[og.tips, "segment.names"]), collapse=sep)
      )
    }
    name.them <- rbind(name.them, new.name)
  }
  
  # Turn into named vectors and add to object
  to.name.segments <- name.them$segment.names
  names(to.name.segments) <- name.them$segment
  object@tree$segment.names <- to.name.segments
  if (!is.null(short.names)) {
    to.name.segments.short <- name.them$short.names
    names(to.name.segments.short) <- name.them$segment
    object@tree$segment.names.short <- to.name.segments.short
  }
  
  # If a force-directed layout has been generated
  if (!is.null(object@tree$walks.force.layout)) {
    # If label positions have been stored, update their names
    if (!is.null(object@tree$walks.force.labels)) {
      object@tree$walks.force.labels$name <- object@tree$segment.names[object@tree$walks.force.labels$seg]
      object@tree$walks.force.labels$name[is.na(object@tree$walks.force.labels$name)] <- ""
      if (!is.null(short.names)) {
        object@tree$walks.force.labels$name.short <- object@tree$segment.names.short[object@tree$walks.force.labels$seg]
        object@tree$walks.force.labels$name.short[is.na(object@tree$walks.force.labels$name.short)] <- ""
      }
    } else {
    # If label positions were not stored because segments were unnamed, generate labels
      object@tree$walks.force.labels <- treeForcePositionLabels(object)
    }
  }
  
  
  return(object)
}