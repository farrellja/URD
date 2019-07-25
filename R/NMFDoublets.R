#' NMF Doublets: Calculate module pair overlap
#' 
#' This considers NMF module expression within cells pairwise
#' to determine how many cells have overlapping expression of that pair of
#' modules. It generates a list that can be used as input into
#' \code{\link{NMFDoubletsPlotModuleThresholds}} to choose which module pairs to
#' use for selection and then \code{\link{NMFDoubletsDetermineCells}} to identify
#' cells that express multiple NMF modules that should not overlap and thus seem
#' to be doublets.
#' 
#' Two thresholds are used to determine whether a cell 'expresses'
#' a module -- this helps differentiate between totally distinct modules (in
#' which case the number of co-expressing cells is similar with either threshold)
#' from pairs of modules expressed in an overlapping gradient (in which case
#' the number of cells that co-express them increases dramatically as the threshold
#' is lowered).
#' 
#' Cells that express a pair of modules that are totally distinct are often
#' useful for identifying doublets -- where one detected cell contains cells that
#' are highly expressing two cell type programs that don't usually occur together.
#' Excluding module pairs that are expressed in an overlapping gradient is important,
#' however, as these often identify cells that are participating in a developmental
#' transition and should not be excluded.
#' 
#' This function requires that a matrix of NMF modules expression per cell is put
#' into slot \code{@@nmf.c1}. It should be a sparse matrix of class dgCMatrix
#' (\code{object@@nmf.c1 <- as(as.matrix(nmf.modules), 'dgCMatrix')}), where rows
#' are cells and columns are modules. In general, we generally normalize each
#' module from 0-1 across all cells.
#' 
#' @importFrom utils combn
#' 
#' @param object An URD object
#' @param modules.use (Character vector) Modules to consider (i.e. column names of \code{@@nmf.c1})
#' @param module.thresh.high (Numeric) High threshold to use for 'expression' of an NMF module
#' @param module.thresh.low (Numeric) Lower threshold to use for 'expression' of an NMF module
#' 
#' @return (List) For use in \code{\link{NMFDoubletsPlotModuleThresholds}} and \code{\link{NMFDoubletsDetermineCells}}
#' @export

NMFDoubletsDefineModules <- function(object, modules.use, module.thresh.high=0.3, module.thresh.low=0.15) {
  # Define a data.frame of module combinations
  module.combos <- as.data.frame(t(combn(modules.use, 2)), stringsAsFactors=F)
  names(module.combos) <- c("Mod1", "Mod2")
  rownames(module.combos) <- paste0(module.combos$Mod1, "-", module.combos$Mod2)
  
  # Determine cells that express each module
  cells.express.mod.high <- lapply(modules.use, function(mod) {
    rownames(object@nmf.c1)[which(object@nmf.c1[,mod] >= module.thresh.high)]
  })
  names(cells.express.mod.high) <- modules.use
  
  cells.express.mod.low <- lapply(modules.use, function(mod) {
    rownames(object@nmf.c1)[which(object@nmf.c1[,mod] >= module.thresh.low)]
  })
  names(cells.express.mod.low) <- modules.use
  
  # Figure out how many cells express each module
  n.cells.express.mod.high <- unlist(lapply(cells.express.mod.high, length))
  names(n.cells.express.mod.high) <- modules.use
  module.combos$n.Mod1.high <- n.cells.express.mod.high[module.combos$Mod1]
  module.combos$n.Mod2.high <- n.cells.express.mod.high[module.combos$Mod2]
  
  n.cells.express.mod.low <- unlist(lapply(cells.express.mod.low, length))
  names(n.cells.express.mod.low) <- modules.use
  module.combos$n.Mod1.low <- n.cells.express.mod.low[module.combos$Mod1]
  module.combos$n.Mod2.low <- n.cells.express.mod.low[module.combos$Mod2]
  
  # Figure out how many cells express each combination of modules
  module.combos$n.overlap.high <- apply(module.combos, 1, function(x){
    length(intersect(cells.express.mod.high[[x[1]]], cells.express.mod.high[[x[2]]]))
  })
  module.combos$frac.overlap.high <- round(module.combos$n.overlap.high / pmin(module.combos$n.Mod1.high, module.combos$n.Mod2.high),digits=4)
  
  module.combos$n.overlap.low <- apply(module.combos, 1, function(x){
    length(intersect(cells.express.mod.low[[x[1]]], cells.express.mod.low[[x[2]]]))
  })
  module.combos$frac.overlap.low <- round(module.combos$n.overlap.low / pmin(module.combos$n.Mod1.low, module.combos$n.Mod2.low),digits=4)
  
  # Difference in frac.overlap between the different thresholds
  module.combos$frac.overlap.diff <- module.combos$frac.overlap.low - module.combos$frac.overlap.high
  
  return(list(
    module.combos=module.combos,
    modules.use=modules.use,
    cells.express.mod.high=cells.express.mod.high,
    cells.express.mod.low=cells.express.mod.low
  ))
}

#' NMF Doublets: Plot module overlaps
#' 
#' This plots each module pair (calculated by \code{\link{NMFDoubletsDefineModules}})
#' in terms of the proportion of 'overlapping' cells (on the x-axis) and the change in
#' proportion of overlapping cells between the high and low thresholds (on the y-axis).
#' This can be used to select parameters for choosing module pairs that likely
#' identify doublets for use in \code{\link{NMFDoubletsDetermineCells}}. See documentation
#' of \code{\link{NMFDoubletsDefineModules}} for more information.
#' 
#' @import ggplot2
#' 
#' @param module.combos (List) Output from \code{\link{NMFDoubletsDefineModules}}
#' @param frac.overlap.max (Numeric) Value to highlight on x-axis as a potential parameter. (Maximum portion of cells that express both modules from a pair.)
#' @param frac.overlap.diff.max (Numeric) Value to highlight on y-axis as a potential parameter. (Maximum change in proportion of cells that express both modules from a pair as the threshold is lowered.)
#' @param modules.highlight (Character Vector) Any module pairs with both modules in this vector will be highlighted in green
#' 
#' @return A ggplot2 object
#' 
#' @export
NMFDoubletsPlotModuleThresholds <- function(module.combos, frac.overlap.max, frac.overlap.diff.max, modules.highlight=NULL) {
  if (class(module.combos)=="list") {
    module.combos <- module.combos$module.combos
  } else if (class(module.combos) != "data.frame") stop("module.combos must be either a list or data.frame")
  the.plot <- ggplot(data=module.combos, aes(x=frac.overlap.high, y=frac.overlap.diff)) + geom_point(alpha=0.1) + geom_hline(yintercept=frac.overlap.diff.max, color='red', lty=2) + geom_vline(xintercept=frac.overlap.max, color='blue', lty=2) + labs(x="Fraction cells expressing both modules at high threshold", y="Proportional increase at low threshold")
  if (!is.null(modules.highlight)) {
    module.combos.highlight <- module.combos[which(module.combos$Mod1 %in% modules.highlight & module.combos$Mod2 %in% modules.highlight),]
    if (nrow(module.combos.highlight) > 0) {
      the.plot <- the.plot + geom_point(data=module.combos.highlight, alpha=1, color='green')
    }
  }
  return(the.plot)
}

#' NMF Doublets: Identify cells that express non-overlapping NMF modules
#' 
#' This identifies which cells express multiple NMF modules that the user
#' has determined should probably not be co-expressed in single-cells, as
#' a way of identifying potential doublets in the data. 
#' 
#' Expression of NMF modules are considered pairwise in 
#' \code{\link{NMFDoubletsDefineModules}}. Thresholds for pairs of NMF modules that
#' probably define doublets (versus ones that might define legitimate developmental
#' transitions) can be determined by using the \code{\link{NMFDoubletsPlotModuleThresholds}}
#' and \code{\link{NMFDoubletsPlotModuleCombos}} functions. 
#' Two thresholds are used to determine whether a cell 'expresses'
#' a module -- this helps differentiate between totally distinct modules (in
#' which case the number of co-expressing cells is similar with either threshold)
#' from pairs of modules expressed in an overlapping gradient (in which case
#' the number of cells that co-express them increases dramatically as the threshold
#' is lowered).
#' 
#' Cells that express a pair of modules that are totally distinct are often
#' useful for identifying doublets -- where one detected cell contains cells that
#' are highly expressing two cell type programs that don't usually occur together.
#' Excluding module pairs that are expressed in an overlapping gradient is important,
#' however, as these often identify cells that are participating in a developmental
#' transition and should not be excluded.
#' 
#' This function requires that a matrix of NMF modules expression per cell is put
#' into slot \code{@@nmf.c1}. It should be a sparse matrix of class dgCMatrix
#' (\code{object@@nmf.c1 <- as(as.matrix(nmf.modules), 'dgCMatrix')}), where rows
#' are cells and columns are modules. In general, we generally normalize each
#' module from 0-1 across all cells.
#' 
#' @param object An URD object
#' @param module.combos (List) Output from \code{\link{NMFDoubletsDefineModules}}
#' @param module.expressed.thresh (Numeric) Threshold for considering a module 'expressed' in a cell
#' @param frac.overlap.max (Numeric) Only pairs of NMF modules where the maximum portion of cells that express both is less than \code{frac.overlap.max} are used.
#' @param frac.overlap.diff.max (Numeric) Only pairs of NMF modules where the portion of co-expressing cells increases less than \code{frac.overlap.diff.max} when the detection threshold is lowered are used. (See \code{\link{NMFDoubletsDefineModules}} for more information.)
#' 
#' @return (Character Vector) Cell Names that express multiple distinct NMF modules
#' 
#' @export
NMFDoubletsDetermineCells <- function(object, module.combos, module.expressed.thresh, frac.overlap.max, frac.overlap.diff.max) {
  
  # Define cells that express modules highly enough for cropping
  cells.express.mod.crop <- lapply(module.combos$modules.use, function(mod) {
    rownames(object@nmf.c1)[which(object@nmf.c1[,mod] >= module.expressed.thresh)]
  })
  names(cells.express.mod.crop) <- module.combos$modules.use
  
  # Module combos to use to define doublets
  module.combos.to.use.for.doublets <- which(module.combos$module.combos$frac.overlap.high <= frac.overlap.max & module.combos$module.combos$frac.overlap.diff <= frac.overlap.diff.max)
  
  # Determine cells that are doublets
  module.combos.doublets <- module.combos$module.combos[module.combos.to.use.for.doublets,]
  
  potential.doublets <- unique(unlist(lapply(1:nrow(module.combos.doublets), function(i) {
    intersect(cells.express.mod.crop[[module.combos.doublets[i,"Mod1"]]],cells.express.mod.crop[[module.combos.doublets[i,"Mod2"]]])
  })))
  
  return(potential.doublets)
  
}

#' Plot cells' expression of modules
#'
#' For inspecting why particular cells do or don't get called as doublets. Makes a
#' dot plot to represent NMF module expression within a particular group of cells.
#' Expression level is represented by color and whether the cell passes the threshold
#' for expression of that module is represented by shape.
#' 
#' @param object An URD object
#' @param modules (Character vector) Modules to plot (Column names of @nmf.c1)
#' @param cells (Character vector) Cells to plot on y-axis
#' @param module.expressed.thresh (Numeric) Threshold for module expression. Will change the shape of points to emphasize whether they made cut-offs.
#' 
#' @return A ggplot2 object
#' 
#' @export

NMFDoubletsPlotModulesInCell <- function(object, modules, cells, module.expressed.thresh=Inf) {
  data.plot <- as.data.frame(as.matrix(object@nmf.c1[cells, modules]))
  data.plot$CELL <- rownames(data.plot)
  data.plot.melt <- reshape2::melt(data.plot, id.vars="CELL", variable.name="MOD")
  data.plot.melt$above.thresh <- data.plot.melt$value >= module.expressed.thresh
  return(ggplot(data=data.plot.melt, aes(x=MOD, y=CELL, color=value, shape=above.thresh)) + theme_bw() + geom_point(size=3) + scale_color_gradientn(colours = defaultURDContinuousColors()) + theme(axis.text.x = element_text(angle = 60, hjust = 1)))
}

#' Plot NMF module combinations and dual-expressing cells
#' 
#' This function creates dimensionality reduction plots that show the expression
#' of pairs of NMF modules and the cells that are selected as doublets, given the
#' chosen set of parameters. It can used to inspect how cells are chosen and refine
#' thresholds to ensure that pairs of modules are not being chosen that might remove
#' legitimate transitions from the data. It produces a plot that shows the expression
#' of each pair of modules (one in red, one in green, overlap in yellow) and a second
#' plot that highlights which cells would be chosen as doublets, given the current
#' parameters. A series of these plots (chosen by \code{boundary}, \code{sort}, and
#' \code{n.plots} parameters) are saved for further inspection.
#' 
#' @param object An URD object
#' @param path (Character) Path to a directory to save plots
#' @param module.combos (List) Output from \code{\link{NMFDoubletsDefineModules}}
#' @param module.expressed.thresh (Numeric) Threshold for calling a cell as an expresser of a module
#' @param frac.overlap.max (Numeric) Threshold to determine which NMF module pairs to use for doublet detecction -- only pairs of NMF modules where the maximum portion of cells that express both is less than \code{frac.overlap.max} are used.
#' @param frac.overlap.diff.max (Numeric) Threshold to determine which NMF module pairs to use for doublet detecction -- only pairs of NMF modules where the portion of co-expressing cells increases less than \code{frac.overlap.diff.max} when the detection threshold is lowered are used.
#' @param boundary (Character) Should plots be made for module combinations that \code{"pass"} the thresholds or that would be \code{"discarded"}
#' @param sort (Character) Should plots be for those module combinations \code{"near"} the boundary or chosen at \code{"random"}
#' @param only.combos.with.doublets (Logical) Should plots be limited to those module combinations that identify doublet cells?
#' @param n.plots (Numeric) What is the maximum number of plots that should be generated?
#' @param width (Numeric) Width of each plot (in pixels)
#' @param height (Numeric) Height of each plot (in pixels)
#' @param ... Additional parameters to pass to \code{\link{plotDimDual}} and \code{\link{plotDimHighlight}} -- for instance to control which dimensionality reduction is plotted.
#' 
#' @return Nothing. A series of plots are generated using ggplot2 and saved to the directory specified in \code{path}.
#' 
#' @export

NMFDoubletsPlotModuleCombos <- function(object, path, module.combos, module.expressed.thresh, frac.overlap.max, frac.overlap.diff.max, boundary=c("pass", "discarded"), sort=c("near", "random"), only.combos.with.doublets=T, n.plots=50, width=1400, height=600, ...) {
  
  if (length(boundary) > 1) boundary <- boundary[1]
  if (length(sort) > 1) sort <- sort[1]
  
  # Grab the module combo information out of the module.combos result
  module.combos <- module.combos$module.combos
  
  # Apply threshold to either select modules that pass the threshold or don't
  if (tolower(boundary)=="pass") {
    module.combos.to.plot <- module.combos[which(module.combos$frac.overlap.high <= frac.overlap.max & module.combos$frac.overlap.diff <= frac.overlap.diff.max),]
    module.combos.to.plot$delta.frac.overlap.high <- frac.overlap.max - module.combos.to.plot$frac.overlap.high
    module.combos.to.plot$delta.frac.overlap.diff <- frac.overlap.diff.max - module.combos.to.plot$frac.overlap.diff
    module.combos.to.plot$delta.min <- pmin(module.combos.to.plot$delta.frac.overlap.high, module.combos.to.plot$delta.frac.overlap.diff)
    module.combos.to.plot$delta.max <- pmax(module.combos.to.plot$delta.frac.overlap.high, module.combos.to.plot$delta.frac.overlap.diff)
    module.combos.to.plot <- module.combos.to.plot[order(module.combos.to.plot$delta.min, module.combos.to.plot$delta.max, decreasing=F),]
  } else if (tolower(boundary)=="discarded") {
    module.combos.to.plot <- module.combos[which(module.combos$frac.overlap.high > frac.overlap.max | module.combos$frac.overlap.diff > frac.overlap.diff.max),]
    module.combos.to.plot$delta.frac.overlap.high <- module.combos.to.plot$frac.overlap.high - frac.overlap.max
    module.combos.to.plot$delta.frac.overlap.diff <- module.combos.to.plot$frac.overlap.diff - frac.overlap.diff.max
    module.combos.to.plot$delta.min <- pmin(module.combos.to.plot$delta.frac.overlap.high, module.combos.to.plot$delta.frac.overlap.diff)
    module.combos.to.plot$delta.max <- pmax(module.combos.to.plot$delta.frac.overlap.high, module.combos.to.plot$delta.frac.overlap.diff)
    module.combos.to.plot <- module.combos.to.plot[order(module.combos.to.plot$delta.max, module.combos.to.plot$delta.min, decreasing=F),]
  } else {
    stop('boundary must be "pass" or "discarded"')
  }
  
  # Determine cells that would be identified as doublets
  modules.use <- unique(c(module.combos.to.plot$Mod1, module.combos.to.plot$Mod2))
  cells.express.mod.crop <- lapply(modules.use, function(mod) {
    rownames(object@nmf.c1)[which(object@nmf.c1[,mod] >= module.expressed.thresh)]
  })
  names(cells.express.mod.crop) <- modules.use
  
  # If desired, crop to only those module pairs with doublets
  if (only.combos.with.doublets) {
    module.combos.to.plot$n.doublets <- unlist(lapply(1:nrow(module.combos.to.plot), function(i) {
      length(intersect(intersect(cells.express.mod.crop[[module.combos.to.plot[i,"Mod1"]]],cells.express.mod.crop[[module.combos.to.plot[i,"Mod2"]]]), colnames(object@logupx.data)))
    }))
    module.combos.to.plot <- module.combos.to.plot[which(module.combos.to.plot$n.doublets > 0),]
  }
  
  # Determine number of modules to plot
  n.plots <- min(n.plots, nrow(module.combos.to.plot))
  if (n.plots == 0) stop("No module combinations to plot are selected by the current parameters.")
  
  # Select module combinations near the boundary or at random
  if (tolower(sort)=="random") {
    module.combos.to.plot <- module.combos.to.plot[base::sort(sample(1:nrow(module.combos.to.plot), n.plots)),]
  } else if (tolower(sort)=="near") {
    module.combos.to.plot <- module.combos.to.plot[1:n.plots,]
  } else {
    stop('sort must be "random" or "near"')
  }
  
  # Actually loop through and generate the plots
  for (i in 1:n.plots) {
    p1.title <- paste0(module.combos.to.plot[i, "Mod1"], " + ", module.combos.to.plot[i, "Mod2"], " (Overlap: ", 100*module.combos.to.plot[i, "frac.overlap.high"], "%; Δ Overlap: ", 100*module.combos.to.plot[i, "frac.overlap.diff"], "%)")
    p1 <- plotDimDual(object, module.combos.to.plot[i,"Mod1"], module.combos.to.plot[i,"Mod2"], plot.title = p1.title, ...)
    overlap.cells <- intersect(intersect(cells.express.mod.crop[[module.combos.to.plot[i,"Mod1"]]],cells.express.mod.crop[[module.combos.to.plot[i,"Mod2"]]]), colnames(object@logupx.data))
    object <- groupFromCells(object, "overlap.cells", overlap.cells)
    p2.title <- paste0("(", length(overlap.cells), " doublets identified)")
    p2 <- plotDimHighlight(object, clustering = "overlap.cells", cluster="TRUE", plot.title=p2.title, highlight.color = "blue", legend.title="Doublet?", ...)
    file.path <- paste0(path, sprintf(fmt = "%04d", i), "-", module.combos.to.plot[i, "Mod1"], "-", module.combos.to.plot[i, "Mod2"], ".png")
    png(file=file.path, width=width, height=height)
    gridExtra::grid.arrange(grobs=list(p1,p2), ncol=2)
    dev.off()
  }
}