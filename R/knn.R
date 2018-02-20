
#' Calculate k-nearest neighbor graph
#' 
#' Calculates and stores a k-nearest neighbor graph based on Euclidean distance in 
#' gene expression space. This is useful for finding and removing outliers prior
#' to calculating a diffusion map.
#' @importFrom RANN nn2
#' @param object An URD object
#' @param genes.use Genes to use for distance (default is variable genes, NULL is all genes.)
#' @param nn (Numeric) Number of nearest neighbors
#' @return An URD object with k-nn graph in slot \code{knn}.
#' @export
calcKNN <- function(object, genes.use=object@var.genes, nn=100) {
  # genes.use is NULL, use all genes.
  if (is.null(genes.use)) genes.use <- rownames(object@logupx.data)
  # Calculate nearest neighbor network
  nearest <- RANN::nn2(t(object@logupx.data[genes.use,]),t(object@logupx.data[genes.use,]),k=max(nn)+1, treetype = "kd", searchtype="priority")
  # Throw away first neighbor because it is the query itself.
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1]
  # Change IDs to cell names
  nearest$nn.cell <- apply(nearest$nn.idx, 2, function(idx) colnames(object@logupx.data)[idx])
  # Set rownames
  rownames(nearest$nn.idx) <- colnames(object@logupx.data)
  rownames(nearest$nn.dists) <- colnames(object@logupx.data)
  rownames(nearest$nn.cell) <- colnames(object@logupx.data)
  # Store parameters used to calculate
  nearest$nn <- nn
  nearest$genes.use <- genes.use
  # Store in URD object and return
  object@knn <- nearest
  return(object)
}

#' Find outliers in terms of their distance in the k-nearest neighbor network
#' 
#' Plots density of cells, given their distance in gene expression to their
#' nearest and nth nearest neighbors. 
#' 
#' Since the diffusion map is calcated on a k-nearest neighbor graph in gene
#' expression space, cells that have unusual distances to their nearest
#' neighbors in a k-nearest neighbor graph often cause poor resulting diffusion
#' maps. Cropping cells based on their distance to their nearest neighbor, and
#' cropping cells that have unusually large distances to an nth nearest neighbor,
#' given the distance to their nearest neighbor, generally produces better results.
#' @param object An URD object
#' @param nn.1 (Numeric) Nearest neighbor to compare on x-axis
#' @param nn.2 (Numeric) Nearest neighbor to compare on y-axis
#' @param x.max (Numeric) Maximum distance to nearest neighbor \code{nn.1} (green line)
#' @param slope.r (Numeric) Slope of red line
#' @param int.r (Numeric) y-intercept of red line
#' @param slope.b (Numeric) Slope of blue line
#' @param int.b (Numeric) y-intercept of blue line
#' @param title (Character) Title of the plot
#' @return Character vector of outlier cells (or non-outliers, if \code{invert=T})
#' @export
knnOutliers <- function(object, nn.1=1, nn.2=20, x.max, slope.r, int.r, slope.b, int.b, title="", invert=F) {
  if (length(object@knn) == 0) stop ("Must calculate a nearest neighbor graph first with calcKNN.")
  gg.data <- as.data.frame(object@knn$nn.dists)
  plot(ggplot(data=gg.data, aes_string(x=paste0("V",nn.1), y=paste0("V",nn.2))) + geom_bin2d(bins=80) + scale_fill_gradientn(trans="log2", colours = defaultURDContinuousColors()) + labs(x=paste("Distance to neighbor", nn.1), y=paste("Distance to neighbor", nn.2), title=title) + geom_abline(slope=slope.r, intercept=int.r, color='red') + geom_abline(slope=slope.b, intercept=int.b, color='blue') + geom_vline(xintercept=x.max, color='green'))
  ks.trim.green <- which(gg.data$V1 > x.max)
  ks.trim.red <- which(gg.data$V20 > int.r + slope.r*gg.data$V1)
  ks.trim.blue <- which(gg.data$V20 > int.b + slope.b*gg.data$V1)
  ks.dm.outliers <- rownames(gg.data)[sort(unique(c(ks.trim.green, ks.trim.red, ks.trim.blue)))]
  if (invert) return(setdiff(rownames(gg.data), ks.dm.outliers)) else return(ks.dm.outliers)
}