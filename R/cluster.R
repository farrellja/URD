#' Spectral Graph Clustering
#' 
#' Computes a spectral graph-based clustering, either on principal components or
#' diffusion components. A k-nearest neighbor graph is computed, its edges are
#' optionally weighted by the Jaccard overlap of cells' neighborhoods, and then
#' a clustering is computed with either the Louvain or Infomap algorithms. The
#' resultant clustering is stored as a column in the slot \code{@@group.ids}. The
#' number of detected clusters is sensitive to the nearest neighbor parameter,
#' and several should be investigated.
#' 
#' @importFrom RANN nn2
#' @importFrom reshape2 melt
#' @importFrom igraph cluster_louvain cluster_infomap graph.data.frame
#' 
#' @param object An URD object
#' @param dim.use (Character) Calculate on principal components (\code{pca}) or diffusion components (\code{dm})
#' @param cells.use (Character vector) Which cells to include in the clustering (default is NULL, which uses all cells)
#' @param which.dims (Numeric vector) Which PCs (or diffusion components) to use. Defaults to the significant PCs. (The default will probably work with diffusion components, though it is non-sensical in that case.)
#' @param num.nn (Numeric or numeric vector) How many nearest-neighbors to use in the k-nn graph. (If multiple values are provided, multiple clusterings are calculated.)
#' @param do.jaccard (Logical) Weight edges in the k-nn graph according to their Jaccard overlap?
#' @param method (Character) Clustering method to use (\code{Louvain})
#' @param group.id (Character) Prefix to use for clustering name (Default is method). Number of nearest neighbors is appended.
#' 
#' @return An URD object with cluster identities saved in \code{@@group.ids} in the column named \code{group.id}.
#' 
#' @examples
#' # Try several different nearest neighbor parameters
#' # Output will be stored as Infomap-10, Infomap-15, ... in object.6s.mnn@group.ids
#' object.6s.mnn <- graphClustering(object.6s.mnn, num.nn = c(10,15,20,30,40), 
#' method="Infomap", do.jaccard = T)
#' 
#' # Cluster on the diffusion map instead of PCA
#' # Output will be stored as Louvain-DM-10, Louvain-DM-20, Louvain-DM-30, etc.
#' object <- graphClustering(object, dim.use="dm", num.nn = c(10,20,30,40), method="Louvain", do.jaccard=T, group.id="Louvain-DM")
#' 
#' @export
graphClustering <- function(object, dim.use=c("pca","dm"), cells.use=NULL, which.dims=which(object@pca.sig), num.nn=30, do.jaccard=TRUE, method=c("Louvain", "Infomap"), group.id=method) {
  if (length(dim.use) > 1) dim.use <- dim.use[1]
  if (length(method) > 1) method <- method[1]
  if (dim.use=="pca") {
    # Get proper PCA data to use
    if (is.null(cells.use)){
      data.use=object@pca.scores[,which.dims]
    } else {
      data.use=object@pca.scores[cells.use,which.dims]
    }
  } else if (dim.use=="dm") {
    # Get proper DM data to use
    if (is.null(cells.use)){
      data.use=object@dm@eigenvectors[,which.dims]
    } else {
      data.use=object@dm@eigenvectors[cells.use,which.dims]
    }
  } else {
    stop("dim.use must be either pca or dm")
  }
  # Get nearest neighbor graph for largest # of NNs requested
  nearest <- nn2(data.use,data.use,k=max(num.nn)+1, treetype = "bd", searchtype="priority")
  nearest$nn.idx <- nearest$nn.idx[,-1]
  nearest$nn.dists <- nearest$nn.dists[,-1]
  
  for (this.nn in num.nn) {
    # Get edges
    edges = melt(t(nearest$nn.idx[,1:this.nn]))
    colnames(edges) = c("B", "A", "C")
    edges = edges[,c("A","B","C")]
    edges$B = edges$C
    edges$C=1
    # Remove repetitions
    edges = unique(transform(edges, A = pmin(A,B), B=pmax(A,B)))
    
    if (do.jaccard){
      NN = nearest$nn.idx[,1:this.nn]
      jaccard_dist = apply(edges, 1, function(x) length(intersect(NN[x[1], ],NN[x[2], ]))/length(union(NN[x[1], ], NN[x[2], ])) )
      edges$C = jaccard_dist
      edges = subset(edges, C != 0)
      edges$C = edges$C/max(edges$C)
    }
    names(edges) <- c("V1", "V2", "weight")
    edges$V1 <- rownames(data.use)[edges$V1]
    edges$V2 <- rownames(data.use)[edges$V2]
    
    g <- graph.data.frame(edges, directed=F)
    if (tolower(method)=="louvain") {
      graph.out = cluster_louvain(g)
    } else if (tolower(method)=="infomap") { 
      graph.out = cluster_infomap(g)
    } else {
      stop("Method should be either Louvain or Infomap.")
    }
    
    clust.assign = factor(graph.out$membership, levels=sort(unique(graph.out$membership)))
    names(clust.assign) = graph.out$names
    k=order(table(clust.assign), decreasing = TRUE)
    new.levels = rep(1,length(unique(graph.out$membership)))
    new.levels[k] = 1:length(unique(graph.out$membership))
    levels(clust.assign) = new.levels
    clust.assign = factor(clust.assign, levels=1:length(unique(graph.out$membership)))
    object@meta$clust = NULL
    object@meta[names(clust.assign),"clust"]=clust.assign
    
    this.group.id <- paste0(group.id, "-", this.nn)
    object@group.ids[names(clust.assign),this.group.id] <- as.character(clust.assign)
  }
  return(object) 
}

#' Return names of cells in a given cluster
#' 
#' @param object An URD object
#' @param clustering (Character) Name of a column in \code{@@group.ids}
#' @param cluster (Character vector) Name of 1 or more clusters to return.
#' 
#' @return Character vector of cell names
#' 
#' @examples 
#' # Get members of cluster 15 from Louvain-10 clustering
#' cellsInCluster(object, clustering="Louvain=10", cluster="15")
#' 
#' # Get members of segments 29, 32, and 79 from reconstructed tree
#' cellsInCluster(object, clustering="segment", cluster=c("29","32","79"))
#' 
#' @export
cellsInCluster <- function(object, clustering, cluster) {
  rownames(object@group.ids)[which(object@group.ids[,clustering] %in% cluster)]
}

#' Cluster centroids in gene expression space
#' 
#' @param object An URD object
#' @param clustering (Character) Name of a clustering (i.e. a column in \code{@@group.ids})
#' @param genes.use (Character vector) Genes to use for calculating distance (default: variable genes, NULL is all genes)
#' 
#' @return A data.frame with mean gene expression per cluster in rows, and clusters as columns.
#' 
#' @examples 
#' cc <- clusterCentroids(object, "Louvain-15", genes.use=object@var.genes)
#' 
#' @export
clusterCentroids <- function(object, clustering, genes.use=object@var.genes) {
  # Get genes to use.
  if (is.null(genes.use)) genes.use <- rownames(object@logupx.data)
  
  # Take mean of gene expression per cluster
  cc <- aggregate(
    x=as.matrix(t(object@logupx.data[genes.use,])), 
    by=list(as.factor(object@group.ids[colnames(object@logupx.data), clustering])), 
    FUN=mean
  )
  
  # Reformat
  rownames(cc) <- cc$Group.1
  cc <- cc[,-1]
  cc <- as.data.frame(t(cc))
  
  return(cc)
}

#' Cluster Differential Expression
#' 
#' Calculates the differentially expressed genes between each cluster and its
#' nearest cluster (as determined by centroids in gene expression space). This
#' can be used to fuse nearby clusters that don't have significant gene expression
#' differences between them. This function uses the binomial differential expression
#' test (see \link{markersBinom}.)
#'
#' @param object An URD object
#' @param clustering (Character) Name of a clustering (i.e. a column in \code{@@group.ids})
#' @param genes.dist (Character vector) Genes to use for calculating distance (default: variable genes, NULL is all genes)
#' @param effect.size (Numeric) Minimum log fold-change for two genes to be considered differentially expressed
#' @param p.thresh (Numeric) Minimum significance value for teo genes to be considered expressed
#' @param frac.must.express (Numeric) Gene must be expressed in at least this fraction of cells in one of the two clusters to be considered.
#' @param genes.de (Character vector) Genes to consider for differential expression (default: NULL is all genes)
#' @param verbose (Logical) Report on progress
#' 
#' @return A list with entries:
#' \itemize{
#' \item{$n.de:} a data frame describing each pair of clusters tested and the number of differentially expressed genes between them
#' \item{$genes:} a list of data frames for each pair of clusters with the differential gene expression test results
#' }
#' 
#' @examples
#' # Test all genes for differential expression again each cluster's nearest 
#' # cluster (by centroid distance in variable gene expression space)
#' cde <- clusterDE(object, clustering="Louvain-15", genes.dist=object@var.genes, 
#' effect.size=log(2), p.thresh=0.01, frac.must.express=0.1, genes.de=NULL, verbose=T)
#' 
#' @export
clusterDE <- function(object, clustering, genes.dist=object@var.genes, effect.size=log(2), p.thresh=0.01, frac.must.express=0.1, genes.de=NULL, verbose=T) {
  # Get cluster centroids
  cc <- clusterCentroids(object, clustering=clustering, genes.use=genes.dist)
  
  # Get cluster distances
  cd <- dist(t(cc))
  
  # Get each cluster's closest neighbor
  cn <- data.frame(
    c1 = labels(cd),
    stringsAsFactors=F
  )
  cn$c2 <- apply(as.matrix(cd), 1, function(x) labels(cd)[order(x)[2]])
  
  # Eliminate duplicated pairs
  cn <- unique(t(apply(cn, 1, sort)))
  
  # Calculate differential expression
  cn.de <- lapply(1:nrow(cn), function(cr) {
    if (verbose) message(paste0(Sys.time(), ": Clusters ", cn[cr,1], " and ", cn[cr,2]))
    markersBinom(object, clust.1=cn[cr,1], clust.2=cn[cr,2], clustering=clustering, effect.size=effect.size, p.thresh=p.thresh, frac.must.express=frac.must.express, genes.use=genes.de)
  })
  names(cn.de) <- apply(cn, 1, paste0, collapse="-")
  
  # Figure out number of differential genes per pair
  cn <- cbind(as.data.frame(cn), unlist(lapply(cn.de, nrow)))
  names(cn) <- c("clust.1", "clust.2", "n.genes.de")
  
  # Put into a list
  return(list(
    n.de=cn,
    genes=cn.de
  ))
}