#' NMF Doublets: Calculate module pair overlap
#' 
#' This considers NMF modules (as loaded into slot \code{@@nmf.c1}) pairwise
#' to determine how many cells have overlapping expression of that pair of
#' modules. Two thresholds are used to determine whether a cell 'expresses'
#' a module -- this helps differentiate between totally distinct modules (in
#' which case the number of co-expressing cells is similar with either threshold)
#' from pairs of modules expressed in an overlapping gradient (in which case
#' the number of cells that co-express them increases dramatically as the threshold
#' is lowered).
#' 
#' @importFrom utils combn
#' 
#' @param object An URD object
#' @param modules.use (Character vector) Modules to consider (i.e. column names of \code{@@nmf.c1})
#' @param module.thresh.high (Numeric) High threshold to use for 'expression' of an NMF module
#' @param module.thresh.low (Numeric) Lower threshold to use for 'expression' of an NMF module
#' 
#' @return (List) For use in \code{\link{NMFDoubletsPlotModuleThresholds}}
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
#' in terms of the propotion of 'overlapping' cells (on the x-axis) and the change in
#' proportion of overlapping cells between the high and low thresholds (on the y-axis).
#' This can be used to select parameters for choosing module pairs that likely
#' identify doublets for use in \code{\link{NMFDoubletsDetermineCells}}.
#' 
#' @import ggplot2
#' 
#' @param module.combos (List) Output from \code{\link{NMFDoubletsDefineModules}}
#' @param frac.overlap.max (Numeric) Value to highlight on x-axis as a potential parameter.
#' @param frac.overlap.diff.max (Numeric) Value to highlight on y-axis as a potential parameter.
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
#' This uses NMF modules (as loaded into slot \code{@@nmf.c1}) to determine cells
#' that express pairs of non-overlapping NMF modules and are likely to be technical
#' doublets.
#' 
#' #' This considers NMF modules (as loaded into slot \code{@@nmf.c1}) pairwise
#' to determine how many cells have overlapping expression of that pair of
#' modules. Two thresholds are used to determine whether a cell 'expresses'
#' a module -- this helps differentiate between totally distinct modules (in
#' which case the number of co-expressing cells is similar with either threshold)
#' from pairs of modules expressed in an overlapping gradient (in which case
#' the number of cells that co-express them increases dramatically as the threshold
#' is lowered).
#' 
#' 
#' @param object An URD object
#' @param module.combos (List) Output from \code{\link{NMFDoubletsDefineModules}}
#' @param module.expressed.thresh (Numeric) Treshold for considering a module 'expressed' in a cell
#' @param frac.overlap.max (Numeric) 
#' @param frac.overlap.diff.max (Numeric)
#' @return (Character Vector) Cell Names that are potentially doublets
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
#' For inspecting why particular cells do or don't get called as doublets
#' Looks for module expression in object@gene.sig.z for now
#' @param object An URD object
#' @param modules (Character vector) Modules to plot (Currently column names of @gene.sig.z)
#' @param cells (Character vector) Cells to plot on y-axis
#' @param module.expressed.thresh (Numeric) Threshold for module expression. Will change the shape of points to emphasize whether they made cut-offs
#' @return A ggplot2 object
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
#' @param object An URD object
#' @param module.combos (List) Output from \code{\link{NMFDoubletsDefineModules}}
#' @param module.expressed.thresh (Numeric) Threshold for calling a cell as an expresser of a module
#' @param frac.overlap.max (Numeric) Maximum 
#' @param frac.overlap.diff.max (Numeric) Max
#' @param boundary (Character) Should plots be for module combinations that \code{"pass"} the thresholds or that would be \code{"discarded"}
#' @param sort (Character) Should plots be for those module combinations \code{"near"} the boundary or chosen at \code{"random"}
#' @param only.combos.with.doublets (Logical) Should plots be limited to those module combinations that identify doublet cells?
#' @param n.plots (Numeric) What is the maximum number of plots that should be generated?
#' @return A ggplot2 object
#' @export

NMFDoubletsPlotModuleCombos <- function(object, path, module.combos, module.expressed.thresh, frac.overlap.max, frac.overlap.diff.max, boundary=c("pass", "discarded"), sort=c("near", "random"), only.combos.with.doublets=T, n.plots=50, width=1400, height=600) {
  
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
    p1 <- plotDimDual(object, module.combos.to.plot[i,"Mod1"], module.combos.to.plot[i,"Mod2"], plot.title = p1.title)
    overlap.cells <- intersect(intersect(cells.express.mod.crop[[module.combos.to.plot[i,"Mod1"]]],cells.express.mod.crop[[module.combos.to.plot[i,"Mod2"]]]), colnames(object@logupx.data))
    object <- groupFromCells(object, "overlap.cells", overlap.cells)
    p2.title <- paste0("(", length(overlap.cells), " doublets identified)")
    p2 <- plotDimHighlight(object, clustering = "overlap.cells", cluster="TRUE", plot.title=p2.title, highlight.color = "blue", legend.title="Doublet?")
    file.path <- paste0(path, sprintf(fmt = "%04d", i), "-", module.combos.to.plot[i, "Mod1"], "-", module.combos.to.plot[i, "Mod2"], ".png")
    png(file=file.path, width=width, height=height)
    gridExtra::grid.arrange(grobs=list(p1,p2), ncol=2)
    dev.off()
  }
}