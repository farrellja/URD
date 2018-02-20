#' Plot 2D Dendrogram of URD Tree
#' 
#' @import ggplot2
#' @importFrom stats aggregate
#' @export
plotTree <- function(object, label=NULL, label.type="search", title=label, legend.title="", plot.tree=T, tree.alpha=1, tree.size=1, plot.cells=T, cell.alpha=0.25, cell.size=0.5, label.x=T, label.segments=F, segment.colors=NULL, discrete.ignore.na=F, color.tree=T, hide.y.ticks=T) {
  
  # Grab various layouts from the object
  segment.layout <- object@tree$segment.layout
  tree.layout <- object@tree$tree.layout
  if (plot.cells) cell.layout <- object@tree$cell.layout
  
  # Initialize ggplot and do basic formatting
  the.plot <- ggplot()
  if (hide.y.ticks) {
    the.plot <- the.plot + scale_y_reverse(c(1,0), name="Pseudotime", breaks=NULL)
  } else {
    the.plot <- the.plot + scale_y_reverse(c(1,0), name="Pseudotime", breaks=seq(0, 1, 0.1))
  }
  the.plot <- the.plot + theme_bw() + theme(axis.ticks=element_blank(), panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  the.plot <- the.plot + labs(x="", title=title, color=legend.title)
  
  # Extract expression information
  if (!is.null(label)) {
    # Grab data to color by
    if (length(label) > 1) stop("Cannot plot by multiple labels simultaneously.")
    color.data <- data.for.plot(object, label=label, label.type=label.type, as.color=F, as.discrete.list = T, cells.use=rownames(object@diff.data))
    color.discrete <- color.data$discrete
    if (color.discrete && color.tree) {
      warning("Coloring the tree with discrete variables often crashes for reasons unknown.")
      color.tree <- F
    }
    color.data <- data.frame(cell=names(color.data$data), value=color.data$data, node=object@diff.data[,"node"], stringsAsFactors=F)
  }
  
  # Summarize expression information if plotting tree
  if (plot.tree && !is.null(label)) {
    if (!color.discrete) {
      # Mean expression per node
      node.data <- aggregate(color.data$value, by=list(color.data$node), FUN="mean.of.logs")
      rownames(node.data) <- node.data$Group.1
      node.data$n <- unlist(lapply(object@tree$cells.in.nodes, length))[node.data$Group.1]
    } else {
      # If uniform expression, then give that output, otherwise give NA.
      node.data <- aggregate(color.data$value, by=list(color.data$node), FUN="output.uniform", na.rm=discrete.ignore.na)
      rownames(node.data) <- node.data$Group.1
      node.data$n <- unlist(lapply(object@tree$cells.in.nodes, length))[node.data$Group.1]
    }
    
    # Color segments according to their expression of their end node
    # (Replace -0 nodes with -1 for getting expression data.)
    tree.layout$node.1 <- gsub("-0","-1",tree.layout$node.1)
    tree.layout$node.2 <- gsub("-0","-1",tree.layout$node.2)
    tree.layout[,"expression"] <- node.data[tree.layout$node.2,"x"]
  }  
  
  # Add cells to graph
  if (plot.cells) {
    if (!is.null(label)) {
      # Add color info to cell.layout
      if (color.discrete) {
        cell.layout$expression <- as.factor(color.data[cell.layout$cell, "value"])
      } else {
        cell.layout$expression <- color.data[cell.layout$cell, "value"]
      }
      # With color
      the.plot <- the.plot + geom_point(data=cell.layout, aes(x=x,y=y,color=expression), alpha=cell.alpha, size=cell.size)
    } else {
      # Just plain black if no label
      the.plot <- the.plot + geom_point(data=cell.layout, aes(x=x,y=y), alpha=cell.alpha, size=cell.size)
    }
  }
  
  # Add tree to graph
  if (plot.tree) {
    if (!is.null(label) && color.tree) {
      # With color, if desired
      the.plot <- the.plot + geom_segment(data=tree.layout, aes(x=x1, y=y1, xend=x2, yend=y2, color=expression), alpha=tree.alpha, size=tree.size, lineend="square") 
    } else {
      # Just plain black if no label
      the.plot <- the.plot + geom_segment(data=tree.layout, aes(x=x1, y=y1, xend=x2, yend=y2), color='black', alpha=tree.alpha, size=tree.size, lineend="square")
    }
  }
  
  # Add color
  if (!is.null(label) && !color.discrete) the.plot <- the.plot + scale_color_gradientn(colors=defaultURDContinuousColors(with.grey=T))
  
  # Label segment names along the x-axis?
  if (label.x) {
    if ("segment.names" %in% names(object@tree)) {
      # Add segment names to segment.layout
      segment.layout$name <- object@tree$segment.names[segment.layout$segment]
      tip.layout <- segment.layout[complete.cases(segment.layout),]
    } else {
      # Find terminal tips
      tip.layout <- segment.layout[which(segment.layout$segment %in% object@tree$tips),]
      tip.layout$name <- as.character(tip.layout$segment)
    }
    the.plot <- the.plot + scale_x_continuous(breaks=as.numeric(tip.layout$x), labels=as.character(tip.layout$name))
    if (any(unlist(lapply(tip.layout$name, nchar)) > 2)) {
      the.plot <- the.plot + theme(axis.text.x = element_text(angle = 68, vjust = 1, hjust=1))
    }
  }
  
  # Label the segments with their number?
  if (label.segments) {
    segment.labels <- as.data.frame(segment.layout[,c("segment","x")])
    segment.labels$y <- apply(object@tree$segment.pseudotime.limits, 1, num.mean)[segment.labels$segment]
    the.plot <- the.plot + geom_label(data=segment.labels, aes(x=x, y=y, label=segment), alpha=0.5)
  }
  
  return(the.plot)
}

#' Is Output Uniform?
#' 
#' @export
output.uniform <- function(x, na.rm=F) {
  y <- unique(x)
  if (na.rm) y <- setdiff(y, NA)
  if (length(y) == 1) return(y) else return(NA)
}
