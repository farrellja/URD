#' Branchpoint Preference Layout
#' 
#' @param object An URD object
#' @param lineages.1 (Character) Segment(s) for the left side of the branchpoint
#' @param lineages.2 (Character) Segment(s) for the right side of the branchpoint
#' @param parent.of.lineages (Character) Segment that is upstream of the branchpoint
#' @param opposite.parent (Character vector) Siblings of the \code{parent.of.lineages}.
#' @param min.visit (Numeric) Minimum visitation frequency to appear in the layout
#' 
#' @return A branchpoint layout that can be used as input to \code{\link{plotBranchpointByGene}}
#' 
#' @export
branchpointPreferenceLayout <- function(object, pseudotime, lineages.1, lineages.2, parent.of.lineages, opposite.parent, min.visit=0) {
  
  # Grab cells along this lineage
  lineage.cells <- cellsAlongLineage(object, unique(c(lineages.1, lineages.2)), remove.root=F)
  lineage.df <- data.frame(
    cell=lineage.cells,
    pseudotime=object@pseudotime[lineage.cells,pseudotime],
    stringsAsFactors=F, row.names = lineage.cells
  )
  
  # Calculate their relative walk frequency
  lineage.df$walk.1 <- apply(object@diff.data[lineage.cells,paste0("visitfreq.raw.", lineages.1), drop=F], 1, max)
  lineage.df$walk.2 <- apply(object@diff.data[lineage.cells,paste0("visitfreq.raw.", lineages.2), drop=F], 1, max)
  
  lineage.df$max.1 <- lineages.1[apply(object@diff.data[lineage.cells,paste0("visitfreq.raw.", lineages.1), drop=F], 1, which.max)]
  lineage.df$max.2 <- lineages.2[apply(object@diff.data[lineage.cells,paste0("visitfreq.raw.", lineages.2), drop=F], 1, which.max)]
  
  # Calculate their preference at the branchpoint
  lineage.df$b.pref <- preference(x=lineage.df$walk.1, y=lineage.df$walk.2, signed=T)
  
  # Calculate preference for these lineages
  parent.vf <- apply(object@diff.data[lineage.cells,paste0("visitfreq.raw.", parent.of.lineages), drop=F], 1, mean)
  opposite.parent.vf <- apply(object@diff.data[lineage.cells,paste0("visitfreq.raw.", opposite.parent), drop=F], 1, mean)
  lineage.df$other.pref <- preference(x=parent.vf, y=opposite.parent.vf, signed=T)
  
  # How much were they visited?
  lineage.df$visited <- log10(apply(lineage.df[,c("walk.1","walk.2")], 1, sum)+1)
  
  # Get rid of cells that show no preference for NP lineage vs. rest of the embryo at all
  lineage.df <- lineage.df[lineage.df$other.pref >= 0.1,]
  lineage.df <- lineage.df[lineage.df$visited >= min.visit,]
  
  return(lineage.df)
}

#' Plot Gene Expression On Branchpoint Layout
#' 
#' @param object An URD object
#' @param lineage.df The output of \code{\link{branchpointPreferenceLayout}}
#' @param label (Character) Value to plot
#' @param label.2 (Character) If dual-color plotting, the second value to plot.
#' @param label.type (Character)
#' @param populations (Character vector, length 2) Labels for the two populations on the x-axis.
#' @param visited.size (Logical) Should the size of points reflect how frequently they were visited by the random walks (to emphasize the more confident data points)
#' @param point.alpha (Numeric) Transparency of points
#' @param pt.lim (Numeric vector) Pseudotime (y) axis limits. Default \code{NULL} lets ggplot determine them automatically.
#' @param color.scale (Character vector) Vector of colors to use as a color scale if \code{label} is continuous values. Does not apply if a second label is given, in which case the color scale will be red-green.
#' @param discrete.colors (Character vector) Vector of colors to use as a color scale if \code{label} is discrete values.
#' @param ylab (Character) Y-axis label
#' @param xlab (Character) X-axis label
#' @param title (Character) Plot title
#' @param axis.lines (Logical) Plot axis lines on the plot?
#' @param legend (Logical) Plot an expression legend?
#' @param fade.low (Numeric) The transparency of the lowest expression points are 
#' decreased by this factor in order to prevent them from blocking the visualization 
#' of expressing cells if they are intermixed.
#' 
#' @return A ggplot2 object
#' 
#' @export
plotBranchpointByGene <- function(object, lineage.df, label, label.2=NULL, label.type="search", populations=NULL, visited.size=T, point.alpha=0.2, pt.lim=NULL, color.scale=NULL, discrete.colors=NULL, ylab="Pseudotime", xlab="Preference", title=label, axis.lines=T, legend=T, fade.low=0.5) {
  if (is.null(color.scale)) color.scale <- defaultURDContinuousColors(with.grey=T)
  if (is.null(label.2)) {
    expression.data <- data.for.plot(object, label=label, label.type=label.type, cells.use=rownames(lineage.df), as.color=F, as.discrete.list = T)
    lineage.df$expression <- expression.data$data
    lineage.df$alpha <- point.alpha
    if (fade.low > 0 & !expression.data$discrete) {
      er <- range(lineage.df$expression)
      fade <- diff(er)*2/9 + er[1]
      fade.alpha <- (fade - lineage.df$expression) / fade * fade.low
      fade.alpha[fade.alpha < 0] <- 0
      lineage.df$alpha <- 1 - fade.alpha
    }
  } else {
    expression.data.1 <- data.for.plot(object, label=label, label.type=label.type, cells.use=rownames(lineage.df), as.single.color=T)
    expression.data.2 <- data.for.plot(object, label=label.2, label.type=label.type, cells.use=rownames(lineage.df), as.single.color=T)
    expression.data <- rgb(expression.data.1, 0.2*expression.data.1+0.2*expression.data.2, expression.data.2)
    lineage.df$expression <- expression.data
    lineage.df$alpha <- 1
  }
  if (visited.size) {
    the.plot <- ggplot(lineage.df, aes_string(y="pseudotime", x="b.pref", color="expression", size="visited", alpha="alpha"))
  } else {
    the.plot <- ggplot(lineage.df, aes_string(y="pseudotime", x="b.pref", color="expression", alpha="alpha")) 
  }
  the.plot <- the.plot + geom_point() + theme_bw() + scale_x_continuous(name=xlab, breaks=c(1,-1), labels=populations) + scale_size_continuous(range=c(0,2)) + scale_alpha_identity() + ggtitle(title)
  if (is.null(pt.lim)) the.plot <- the.plot + scale_y_reverse(name=ylab) else the.plot <- the.plot + scale_y_reverse(limits=pt.lim, name=ylab)
  if (is.null(label.2) && !expression.data$discrete) the.plot <- the.plot + scale_color_gradientn(colors=color.scale, name=label)
  if (!is.null(label.2)) the.plot <- the.plot + scale_color_identity()
  if (!is.null(discrete.colors) && expression.data$discrete) the.plot <- the.plot + scale_color_manual(values=discrete.colors)
  if (!legend) the.plot <- the.plot + guides(size=F, color=F)
  if (!axis.lines) the.plot <- the.plot + theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank())
  return(the.plot)
}

