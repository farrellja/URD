#' Find variable genes
#' 
#' Single-cell RNAseq data is noisy, so we perform our analyses using only those genes that exhibit greater variability than those of similar expression levels. In theory, those genes have biological variability across cells in addition to their technical variability.
#' 
#' A null mathematical model is built to model the relationship between average UMI counts and  coefficient of variation (CV) across all genes, based on a negative binomial distribution that incorporates sampling noise and relative library size. Those genes that have a CV greater than the null model (threshold determined by \code{diffCV.cutoff}) are chosen as variable.
#' 
#' If \code{do.plot=T}, produces three plots: the first shows the relative library sizes and the gamma distribution fit to them. The second shows a histogram of each gene's CV ratio to the null for its mean expression level and the \code{diffCV.cutoff} threshold chosen. The third shows each gene's mean expression and CV, the determined null model (in pink), and whether the gene was selected as variable (green genes were variable).
#' 
#' @references Pandey S, Shekhar K, Regev A, and Schier AF. Comprehensive Identification and Spatial Mapping of Habenular Neuronal Types Using Single-Cell RNA-Seq. 2018. Current Biology 28(7):1052-1065. DOI: https://doi.org/10.1016/j.cub.2018.02.040
#' 
#' @importFrom MASS fitdistr
#' @importFrom Matrix rowMeans
#' 
#' @param object An URD object
#' @param cells.fit (Character Vector) Cells to use for finding variable genes (if \code{NULL}, uses all cells.)
#' @param set.object.var.genes (Logical) Return an object with \code{@@var.genes} set (if \code{TRUE}) or return a character vector of variable genes (if \code{FALSE})
#' @param diffCV.cutoff (Numeric) Difference in coefficient of variation (CV) between null distribution & genes that will be considered variable (Difference is in log space, so this amounts to a fold-change)
#' @param mean.min (Numeric) Genes must have this minimum mean expression to be selected (Use to eliminate noisy lowly expressed genes)
#' @param mean.max (Numeric) Genes must have less than this maximum mean expression to be selected (Use to eliminate the high end if the null distribution fits very poorly there)
#' @param main.use (Character) Title to display for the overall three-panel plot
#' @param do.plot (Logical) Whether or not to display plots
#' @param max.genes.in.ram (Numeric) Maximum number of genes to process at a time. This can be used to prevent memory errors on large data sets.
#' 
#' @return Either an URD object with \code{@@var.genes} set (if \code{set.object.var.genes=T}) or character vector of variable genes (if \code{set.object.var.genes=F})
#' 
#' @examples
#' # Find a list of cells from each stage.
#' stages <- unique(object@meta$STAGE)
#' cells.each.stage <- lapply(stages, function(stage) rownames(object@meta)[which(object@meta$STAGE ==                                                                                 stage)])
#' # Compute variable genes for each stage.
#' var.genes.by.stage <- lapply(1:length(stages), function(n) findVariableGenes(object, cells.fit = cells.each.stage[[n]], set.object.var.genes = F, diffCV.cutoff = 0.3, mean.min = 0.005, mean.max = 100, main.use = stages[[n]], do.plot = T))
#' 
#' @export

findVariableGenes <- function(object, cells.fit=NULL, set.object.var.genes=T, diffCV.cutoff=0.5, mean.min=0.005, mean.max=100, main.use="", do.plot=T, max.genes.in.ram=1.6e9/ncol(object@count.data)) {
  
  # Format for the plot
  if (do.plot) par(mfrow=c(1,3), oma=c(0,0,2,0))
  
  # If cells.fit not provided, use all cells.
  if (is.null(cells.fit)) cells.fit <- colnames(object@count.data)
  
  # Calculate empirical mean, var and CV
  mean_emp <- Matrix::rowMeans(object@count.data[,cells.fit])
  var_emp <- sparseApply(object@count.data[,cells.fit], 1, var, max.dim = max.genes.in.ram, verbose=T)
  cv_emp <- sqrt(var_emp) / mean_emp
  
  # Fit gamma distribution to 'size factors' to build negative binomial
  a <- Matrix::colSums(object@count.data[,cells.fit])
  size_factor <- a/mean(a)
  fit <- fitdistr(size_factor, "Gamma")
  if (do.plot) {
    hist(size_factor, 50, probability=TRUE, xlab="UMIs per cell / mean UMIs per cell", main=paste0("Size Factors & Gamma Fit (a=", round(fit$estimate[1], 2), ")"))
    curve(dgamma(x, shape=fit$estimate[1], rate=fit$estimate[2]),from=0, to=quantile(size_factor, 0.999), add=TRUE, col="red")
  }
  
  # Create null negative binomial distribution
  #   Gamma distributions of individual genes are just scaled versions. If X ~ Gamma(a,b)
  #   then cX ~ Gamma(a, b/c)
  a_i <- rep(fit$estimate[1], length(mean_emp))
  names(a_i) <- names(mean_emp)
  b_i <- fit$estimate[2] / mean_emp
  names(b_i) <- names(mean_emp)
  mean_NB <- a_i / b_i
  var_NB <- a_i*(1+b_i) / (b_i^2)
  cv_NB <- sqrt(var_NB)/mean_NB
  
  # Calculate difference in genes' CV and null CV
  diffCV = log(cv_emp) - log(cv_NB)
  if (do.plot) {
    hist(diffCV, 150, xlab="log(CVgene / CVnull)", main="Diff CV")
    abline(v=diffCV.cutoff, lty=2, col='red')
  }
  pass.cutoff <- names(diffCV)[which(diffCV > diffCV.cutoff & (mean_emp > mean.min & mean_emp < mean.max))]
  
  # Plot variable gene selection
  if (do.plot) {
    plot(mean_emp,cv_emp,pch=16,cex=0.5,col="black",xlab="Mean Counts",ylab="CV (counts)", log="xy", main = "Selection of Variable Genes")
    #curve(sqrt(1/x), add=TRUE, col="red", log="xy", lty=2, lwd=2)
    or = order(mean_NB)
    lines(mean_NB[or], cv_NB[or], col="magenta", lwd=2)
    points(mean_emp[pass.cutoff], cv_emp[pass.cutoff], pch=16, cex=0.5, col='green')
    title(main.use, outer=T)
  }
  
  # Return list of variable genes
  if (set.object.var.genes) {
    object@var.genes <- pass.cutoff
    return(object)
  } else {
    return(pass.cutoff)
  }
}

