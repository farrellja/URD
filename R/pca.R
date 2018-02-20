#' Calculate principal components
#' 
#' This performs principal component analysis on centered and scaled expression data,
#' which is stored in \code{pca.scores} and \code{pca.load}. Then, the principal components
#' that are significant are estimated using the Marchenko Pastur Law (by calling
#' \code{\link{pca.marchenko.pastur}}, which is stored in \code{pca.sig}. By default, 
#' twice as many PCs are saved as are determined as significant, which can be adjusted
#' with the \code{store.thresh} parameter.
#' 
#' @importFrom gmodels fast.prcomp
#' 
#' @param object URD object
#' @param genes.use (Character Vector) Genes to use for principal components analysis (default: stored variable genes. Set NULL to use all genes.)
#' @param pcs.store (Numeric) Number of PCs to retain (if NULL, will determine using \code{store.thresh})
#' @param store.thresh (Numeric) If \code{pcs.store} isn't specified, stores the number of significant PCs times this number
#' @param do.print (Logical) Report determined Marchenko-Pastur values for significant PCs.
#' @param verbose (Logical) Whether to report on progress
#' 
#' @return An URD object, with loading of genes into PCs in slot \code{pca.load}, PC scores for each cell in slot \code{pca.scores}, and the significance of each PC stored in slot \code{pca.sig}.
#' 
#' @export
calcPCA <- function(object, genes.use=object@var.genes, pcs.store=NULL, store.thresh=2, do.print=T, verbose=T) {
            # Check whether there ARE variable genes
            if (length(genes.use) == 0) {
              warning("Variable genes have not been stored. Using all genes instead.")
              genes.use <- rownames(object@logupx.data)
            }
            if (is.null(genes.use)) genes.use <- rownames(object@logupx.data)
  
            # Get z-scored data for using in PCA.
            if (verbose) print(paste0(Sys.time(), ": Centering and scaling data."))
            data.use <- get.z.data(object, genes=genes.use)
            pc.genes <- intersect(genes.use, rownames(data.use))
            
            # Remove genes with zero variation
            if (verbose) print(paste0(Sys.time(), ": Removing genes with no variation."))
            pc.genes.var <- apply(data.use[pc.genes,],1,function(x) var(x))
            pc.genes <- pc.genes[pc.genes.var>0]
            data.use <- data.use[pc.genes,]
            
            # Do PCA calculation
            if (verbose) print(paste0(Sys.time(), ": Calculating PCA."))
            pca.obj <- fast.prcomp(t(data.use),center=FALSE,scale=FALSE)
            
            # Estimate significant PCs
            if (verbose) print(paste0(Sys.time(), ": Estimating significant PCs."))
            pca.sig <- pcaMarchenkoPastur(M=dim(data.use)[1], N=dim(data.use)[2], pca.sdev = pca.obj$sdev)
            
            # Figure out how many PCs to store
            if (is.null(pcs.store)) {
              max.sig.pc <- max(which(pca.sig))
              pcs.store <- ceiling(max.sig.pc * store.thresh)
            }
            if (do.print) print(paste("Storing", pcs.store, "PCs."))
            
            # Store values
            object@pca.scores <- data.frame(pca.obj$x[,1:pcs.store])
            object@pca.load <- data.frame(pca.obj$rotation[,1:pcs.store])
            object@pca.sdev <- pca.obj$sdev[1:pcs.store]
            object@pca.sig <- pca.sig[1:pcs.store]
            return(object)
}

#' Marchenko-Pastur Significant PCs
#' 
#' The Marchenko Pastur Law (MP) predicts the theoretical upper and lower bounds
#' on the null distribution of eigenvalues for an MxN random matrix. We take
#' significant principal components (PCs) as those with eigenvalues greater than
#' the maximum eigenvalue predicted for random data. This function assumes that
#' the data has mean 0 and variance 1 (i.e. that the data has been centered and
#' scaled). This called automatically by \code{\link{calc.PCA}} and the results
#' are stored in slot \code{pca.sig}.
#' 
#' @param M (Numeric) Number of rows in input data
#' @param N (Numeric) Number of columns in input data
#' @param pca.sdev (Numeric vector) Standard deviations for each principal component
#' @param do.print (Logical) Whether to report the results
#' @return Logical vector of whether each PC is significant.
#' 
#' @export
pcaMarchenkoPastur <- function(M, N, pca.sdev, do.print=T) {
  pca.eigenvalue <- (pca.sdev)^2
  marchenko.pastur.max <- (1+sqrt(M/N))^2
  pca.sig <- pca.eigenvalue > marchenko.pastur.max
  if (do.print) {
    print(paste("Marchenko-Pastur eigenvalue null upper bound:", marchenko.pastur.max))
    print(paste(length(which(pca.eigenvalue > marchenko.pastur.max)), "PCs have larger eigenvalues."))
  }
  return(pca.sig)
}

#' PC Top Loaded Genes
#' 
#' @param object An URD object
#' @param pcs.print (Numeric vector) Which PCs to report? (Default: significant PCs)
#' @param genes.print (Numeric) How many genes to report
#' 
#' @export
pcTopGenes <- function(object, pcs.print=which(object@pca.sig), genes.print=20) {
  top.loadings <- lapply(pcs.print, function(pc) {
    gene.order <- order(abs(object@pca.load[,pc]), decreasing = T)
    pc.top <- data.frame(
      gene=rownames(object@pca.load)[gene.order[1:genes.print]],
      loading=object@pca.load[gene.order[1:genes.print],pc],
      stringsAsFactors = F
    )
    return(pc.top)
  })
  names(top.loadings) <- pcs.print
  return(top.loadings)
}