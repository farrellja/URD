#' Layout Tree Dendrogram
#' 
#' This function is called automatically by \code{\link{buildTree}}. 
#' 
#' @importFrom igraph graph_from_edgelist 
#' @importFrom ggraph create_layout
#' 
#' @param object An URD object
#' 
#' @return An URD object with a dendrogram layout of segments stored in \code{object@tree$segment.layout}.
#' 
#' @export
#' 
#' @keywords internal
treeLayoutDendrogram <- function(object) {
  # Convert edge list to igraph format
  object@tree$tree.igraph <- graph_from_edgelist(el=as.matrix(object@tree$segment.joins[,1:2]))
  # Generate tree layout using igraph layout function
  igraph.layout <- create_layout(object@tree$tree.igraph, layout="dendrogram")
  # Label segments in layout
  igraph.layout$segment <- as.character(igraph.layout$name)
  # Add it to the object
  object@tree$segment.layout <- igraph.layout
  
  return(object)
}

#' Elaborate Tree Dendrogram With Nodes
#' 
#' This function is called automatically by \code{\link{buildTree}}.
#' 
#' @param object An URD object
#' 
#' @return An URD object with a dendrogram that includes a short segment for
#' each node stored in \code{object@tree$tree.layout}.
#' 
#' @export
#' 
#' @keywords internal
treeLayoutElaborate <- function(object) {
  # Grab igraph layout
  igraph.layout <- object@tree$segment.layout
  
  # Generate tree edges layout
  tree.layout <- as.data.frame(object@tree$edge.list, stringsAsFactors=F)
  names(tree.layout) <- c("node.1", "node.2")
  rownames(tree.layout) <- NULL
  tree.layout$segment.1 <- unlist(lapply(strsplit(tree.layout$node.1, "-"), function(x) x[1]))
  tree.layout$segment.2 <- unlist(lapply(strsplit(tree.layout$node.2, "-"), function(x) x[1]))
  tree.layout[,"x1"] <- igraph.layout[tree.layout$segment.1,"x"]
  tree.layout[,"y1"] <- object@tree$node.max.pseudotime[tree.layout$node.1]
  tree.layout[,"x2"] <- igraph.layout[tree.layout$segment.2,"x"]
  tree.layout[,"y2"] <- object@tree$node.max.pseudotime[tree.layout$node.2]
  tree.layout$do.mean <- F
  
  # Add the root
  # Determine root - must be the first node of the last segment added
  #root.node <- tail(object@tree$segment.joins, 1)$parent
  root.segment <- tail(object@tree$segments, 1)
  root.layout <- dim(tree.layout)[1] + 1
  tree.layout[root.layout, "node.1"] <- paste0(root.segment, "-0")
  tree.layout[root.layout, "node.2"] <- paste0(root.segment, "-1")
  tree.layout[root.layout, c("segment.1", "segment.2")] <- rep(root.segment, 2)
  tree.layout[root.layout, c("x1", "x2")] <- igraph.layout[root.segment,"x"]
  tree.layout[root.layout, "y1"] <- 0
  tree.layout[root.layout, "y2"] <- object@tree$node.max.pseudotime[paste0(root.segment,"-1")]
  tree.layout[root.layout, "do.mean"] <- F
  
  # Get rid of any NAs because node.max.pseudotime is NA (if there were no cells there)
  #tree.layout[is.na(tree.layout)] <- 0
  tree.layout[is.na(tree.layout$y1),"y1"] <- as.numeric(object@tree$segment.pseudotime.limits[tree.layout[is.na(tree.layout$y1),"segment.1"], "start"])
  tree.layout[is.na(tree.layout$y2),"y2"] <- as.numeric(object@tree$segment.pseudotime.limits[tree.layout[is.na(tree.layout$y2),"segment.2"], "end"])
  
  # Adjust for edges with different x-values
  to.square <- which(tree.layout$x1 != tree.layout$x2)
  squared.layout <- do.call("rbind", lapply(to.square, function(line) {
    df.out <- tree.layout[c(line,line),]
    new.node <- gsub("-1", "-0", df.out[1,"node.2"])
    df.out[1,"node.2"] <- new.node
    df.out[2,"node.1"] <- new.node
    df.out[1,"do.mean"] <- T
    df.out[1,"y2"] <- df.out[1,"y1"]
    df.out[2,"x1"] <- df.out[2,"x2"]
    df.out[2,"segment.1"] <- df.out[2,"segment.2"]
    return(df.out)
  }))
  rownames(squared.layout) <- NULL
  rows.keep <- setdiff(1:dim(tree.layout)[1], to.square)
  tree.layout <- rbind(tree.layout[rows.keep,], squared.layout)
  
  # Add it to the object
  object@tree$tree.layout <- tree.layout
  return(object)
}

#' Add cells to tree layout
#' 
#' This function is called automatically by \code{\link{buildTree}}.
#' 
#' @importFrom stats runif
#' 
#' @param object An URD object
#' @param pseudotime (Character) Pseudotime to use (i.e. a column name of \code{object@pseudotime}).
#' @param jitter (Numeric) Distance to jitter cells before placing next to dendrogram. (Relative to 1, which is the distance between terminal segments).
#' @param jitter.push (Numeric) Distance next to lines of dendrogram that should not include any cells.
#' 
#' @return An URD object with cell locations added to \code{object@tree$cell.layout}.
#' 
#' @export
#' 
#' @keywords internal
treeLayoutCells <- function(object, pseudotime, jitter=0.15, jitter.push=0.05) {
  # Grab tree layout
  tl <- object@tree$tree.layout
  # Layout all cells
  cell.layout <- do.call("rbind",lapply(names(object@tree$cells.in.nodes), function(node) {
    # Get layout information for the graph segment these cells are associated with
    tls <- tl[which(tl$node.2==node),]
    # Build a data frame for these cells
    cell.layout <- data.frame(
      cell = object@tree$cells.in.nodes[[node]],
      y = object@pseudotime[object@tree$cells.in.nodes[[node]], pseudotime],
      stringsAsFactors=F
    )
    n <- dim(cell.layout)[1]
    # Figure out the magnitude of each point's jitter
    cell.layout$jitter <- runif(n, min=jitter.push, max=jitter+jitter.push) 
    # Figure out which way cells should be jittered
    cell.layout$jitter.dir <- sample(c(-1,1), n, replace = T)
    # Compute final x-coordinate
    cell.layout$x <- rep(tls$x1, n) + cell.layout$jitter * cell.layout$jitter.dir
    return(cell.layout)
  }))
  cell.layout <- cell.layout[,c("cell","x","y")]
  rownames(cell.layout) <- cell.layout$cell
  object@tree$cell.layout <- cell.layout
  return(object)
}