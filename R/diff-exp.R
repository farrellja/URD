#' Test for differential gene expression in two populations using a binomial test
#' 
#' @importFrom stats pbinom p.adjust
#' 
#' @param object An URD object
#' @param pseudotime (Character)
#' @param clust.1 (Character)
#' @param clust.2 (Character)
#' @param cells.1 (Character vector)
#' @param cells.2 (Character vector)
#' @param clustering (Character) Name of a clustering (i.e. a column in \code{@@group.ids})
#' @param effect.size (Numeric) Minimum log fold-change for two genes to be considered differentially expressed
#' @param p.thresh (Numeric) Minimum significance value for teo genes to be considered expressed
#' @param frac.must.express (Numeric) Gene must be expressed in at least this fraction of cells in one of the two clusters to be considered.
#' @param genes.use (Character vector) Genes to compare, default is NULL (all genes)
#' 
#' @return (data.frame)
#' 
#' @references 
#' Shekhar et al., 2016, Cell 166, 1308–1323
#' 
#' @export
markersBinom <- function(object, pseudotime, clust.1=NULL,clust.2=NULL,cells.1=NULL,cells.2=NULL,clustering=NULL,effect.size=log(2),p.thresh=.01,frac.must.express=0.1,genes.use=NULL) {
            if (is.null(genes.use)) genes.use <- rownames(object@logupx.data)
            
            if (is.null(clustering)) {
              # Must provide clustering if either cell list should be pulled from clustering
              if (!is.null(clust.1) | !is.null(clust.2)) stop("If clust.1 or clust.2 is set, clustering must also be provided.")
            } else {
              clust.use <- object@group.ids[,clustering]
              names(clust.use) <- rownames(object@group.ids)
            }
            
            if (is.null(cells.1)) {
              if (is.null(clust.1)) stop("Must provide either cells.1 or clust.1")
              cells.1 <- names(clust.use[which(clust.use%in%clust.1)])
            }
            
            if (is.null(cells.2)) {
              if (is.null(clust.2)) {
                clust.2 <- "rest"
                cells.2 <- setdiff(names(clust.use), cells.1)
              } else {
                cells.2 <- names(clust.use[which(clust.use%in%clust.2)])
              }
            }
            
            # Calculate proportion of expressing cells and regularize
            m = apply(object@logupx.data[genes.use, cells.2], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #2
            m1 = m; m1[m==0]=1; # Regularization. Add a pseudo count of 1 to unexpressed markers to avoid false positives
            n = apply(object@logupx.data[genes.use, cells.1], 1, function(x) sum(x>0)) #Number of cells expressing marker in cluster #1
            n1 = n; n1[n==0]=1; # Regularization.  Add a pseudo count of 1 to unexpressed markers to avoid false positives
            log_fold_express = log(n*length(cells.2)/(m*length(cells.1))) #log proportion of expressing cells
            
            #Test for enrichments in cluster #1
            #Find the probability of finding n or more cells +ve for marker in cluster 1 given the fraction in cluster 2
            pv1 <- pbinom(n, length(cells.1), m1/length(cells.2), lower.tail = FALSE) + dbinom(n, length(cells.1), m1/length(cells.2))
            d1 <- data.frame(log.effect=log_fold_express,pval=pv1)
            d1$pval.fdr <- p.adjust(d1$pval, method="fdr")
            d1 <- subset(d1, log.effect >= effect.size)
            
            #Enrichments in cells.2
            #Find the probability of finding m or more cells +ve for marker in cluster 2 given the fraction in cluster 1
            pv2 <- pbinom(m, length(cells.2), n1/length(cells.1), lower.tail=FALSE) + dbinom(m, length(cells.2), n1/length(cells.1))
            d2 <- data.frame(log.effect=log_fold_express,pval=pv2)
            d2$pval.fdr <- p.adjust(d2$pval, method="fdr")
            d2 <- subset(d2, log.effect <= -effect.size)
            
            # Combine enrichments
            result <- rbind(d1, d2);
            
            # Fraction expressing & mean transcripts
            result$posFrac.1 <- round(n/length(cells.1),3)[rownames(result)]
            result$posFrac.2 <- round(m/length(cells.2),3)[rownames(result)]
            result$nTrans.1 <- apply(object@count.data[rownames(result), cells.1], 1, function(x) round(mean(x),3))
            result$nTrans.2 <- apply(object@count.data[rownames(result), cells.2], 1, function(x) round(mean(x),3))
            
            # Trim to those results that are expressed in a sufficient proportion of the population.
            if (!is.null(clust.2) && clust.2=="rest"){
              genes.include = result$posFrac.1 >= frac.must.express & result$pval.fdr <= p.thresh
            } else{
              genes.include = (result$posFrac.1 >= frac.must.express | result$posFrac.2 >= frac.must.express) & result$pval.fdr <= p.thresh
            }
            result <- result[genes.include,]
            result <- result[order(abs(result$log.effect), decreasing=TRUE),]
            
            # Name columns & return result.
            colnames(result)[4:7] <- paste(c("posFrac", "posFrac", "nTrans", "nTrans"), c(clust.1, clust.2, clust.1, clust.2), sep="_")
            
            return(result)
      } 

#' Test for markers of a population using a precision-recall curve
#' 
#' @importFrom stats p.adjust
#' 
#' @param object An URD object
#' @param clust.1 (Character)
#' @param clust.2 (Character)
#' @param cells.1 (Character vector)
#' @param cells.2 (Character vector)
#' @param clustering (Character) Name of a clustering (i.e. a column in \code{@@group.ids})
#' @param effect.size (Numeric) Minimum log fold-change for two genes to be considered differentially expressed
#' @param frac.must.express (Numeric) Gene must be expressed in at least this fraction of cells in one of the two clusters to be considered.
#' @param frac.min.diff (Numeric) Fraction of cells expressing the gene must be at least this different between two populations to be considered.
#' @param genes.use (Character vector) Genes to compare, default is NULL (all genes)
#' 
#' @return (data.frame)
#' 
#' @export
markersAUCPR <- function(object, clust.1=NULL, clust.2=NULL, cells.1=NULL, cells.2=NULL, clustering=NULL, effect.size=0.25, frac.must.express=0.1, frac.min.diff=0, genes.use=NULL) {
  if (is.null(genes.use)) genes.use <- rownames(object@logupx.data)
  
  if (is.null(clustering)) {
    # Must provide clustering if either cell list should be pulled from clustering
    if (!is.null(clust.1) | !is.null(clust.2)) stop("If clust.1 or clust.2 is set, clustering must also be provided.")
  } else {
    clust.use <- object@group.ids[,clustering]
    names(clust.use) <- rownames(object@group.ids)
  }
  
  if (is.null(cells.1)) {
    if (is.null(clust.1)) stop("Must provide either cells.1 or clust.1")
    cells.1 <- names(clust.use[which(clust.use%in%clust.1)])
  }
  
  if (is.null(cells.2)) {
    if (is.null(clust.2)) {
      clust.2 <- "rest"
      cells.2 <- setdiff(names(clust.use), cells.1)
    } else {
      cells.2 <- names(clust.use[which(clust.use%in%clust.2)])
    }
  }
  
  # Figure out proportion expressing and mean expression
  genes.data <- data.frame(
    frac.1=round(apply(object@logupx.data[genes.use, cells.1, drop=F], 1, prop.exp), digits=3),
    frac.2=round(apply(object@logupx.data[genes.use, cells.2, drop=F], 1, prop.exp), digits=3),
    exp.1=round(apply(object@logupx.data[genes.use, cells.1, drop=F], 1, mean.of.logs), digits=3),
    exp.2=round(apply(object@logupx.data[genes.use, cells.2, drop=F], 1, mean.of.logs), digits=3)
  )
  genes.data$exp.fc <- genes.data$exp.1 - genes.data$exp.2
  
  # Throw out genes that don't mark either population or obviously change between the two
  # populations to reduce downstream computation.
  genes.use <- names(which(
    (apply(genes.data[,c("frac.1", "frac.2")], 1, max) > frac.must.express) &
    (apply(genes.data[,c("frac.1", "frac.2")], 1, function(x) abs(diff(x))) > frac.min.diff) & 
    (apply(genes.data[,c("exp.1", "exp.2")], 1, function(x) abs(diff(x))) > effect.size)
  ))
  genes.data <- genes.data[genes.use,]
  
  # Calculate area under precision-recall curve for each gene (and set to 0 if NA is returned)
  genes.data$AUCPR <- unlist(lapply(genes.use, function(gene) differentialAUCPR(object@logupx.data[gene,cells.1], object@logupx.data[gene,cells.2])))
  genes.data[is.na(genes.data$AUCPR),"AUCPR"] <- 0
  
  # Order by AUC
  genes.data <- genes.data[order(genes.data$AUCPR, decreasing=T),c("AUCPR", "exp.fc", "frac.1", "frac.2", "exp.1", "exp.2")]
  names(genes.data)[3:6] <- paste(c("posFrac", "posFrac", "nTrans", "nTrans"), c(clust.1, clust.2, clust.1, clust.2), sep="_")
  
  # Return
  return(genes.data)
} 

#' Area under precision-recall curve for determining marker specificity.
#' 
#' @importFrom caTools trapz
#' @importFrom ROCR prediction performance
#'
#' @param x (Numeric vector) Gene expression values from one population
#' @param y (Numeric vector) Gene expression values from the other population
#' 
#' @return (Numeric) Area under the precision-recall curve
differentialAUCPR = function (x, y) {
  prediction.use <- ROCR::prediction(predictions = c(x, y), labels = c(rep(x = 1, length(x = x)), rep(x = 0, length(x = y))), label.ordering = 0:1)
  perf1 <- ROCR::performance(prediction.use, "prec", "rec")
  is.nanvals = is.nan(perf1@x.values[[1]]) | is.nan(perf1@y.values[[1]])
  auc = caTools::trapz(perf1@x.values[[1]][!is.nanvals], 
              perf1@y.values[[1]][!is.nanvals])
  auc.use <- round(x = auc, digits = 3)
  return(auc.use)
}

# V5 Things looked pretty good, but now you changed what you do with NMF modules.
# So — now do comparisons segment-by-segment and only accept ones that beat all of their siblings and lose to no sibling.
# only.return.global (logical): will cull the results to those markers that pass a test across all cells.

#' Test for differential gene expression along tree using binomial test
#' 
#' This test performs differential expression testing along URD's recovered tree
#' structure using a binomial expression test (see: \code{\link{markersBinom}}).
#' 
#' It starts from the segment provided in \code{tips} and compares it to each of
#' its siblings (and their descendants) matched to the same pseudotime limits as
#' the segment under consideration. Markers that are upregulated in a segment compared
#' to either any of its siblings (if \code{must.beat.all.sibs=FALSE}) or compared
#' to each of its siblings (if \code{must.beat.all.sibs=TRUE}) become putative
#' markers.
#' 
#' 
#' @param object An URD object
#' @param pseudotime (Character) Name of pseudotime to use
#' @param tips (Character vector) The tip of the trajectory to find differential expression for (either as tip number or name)
#' @param log.effect.size (Numeric) Minimum fold-change difference in expression to be considered differential (default, \code{log(2)} is 2-fold change)
#' @param p.thresh (Numeric) Significance threshold for being considered differential
#' @param frac.must.express (Numeric) Fraction of cells of interest that must express the gene to consider it differential.
#' @param genes.use (Character vector) Which genes to test for differential expression (default \code{NULL} is all genes)
#' @param root (Character) Which segment of the tree (by number) to use as the root for differential express (default \code{NULL} will autodetect, but this can be useful if you do not want to traverse up the entire tree)
#' @param only.return.global (Logical) Genes must also pass an overall in-trajectory vs. out-of-trajectory test across all cells in the data to be considered differential (default, \code{FALSE})
#' @param must.beat.sibs (Numeric) At multi-furcated branchpoints, independent comparisons are made to each sibling. For a gene to be considered differential it must be a marker against this proportion of its sibling branches. (1 means it must be a marker against all other branches, 0 means it must be a marker against a single other branch.)
#' @param report.stats (Logical) Should information be returned about all of the comparisons? (nGene, nTrans, pseudotime, n.cells) If \code{TRUE}, this function returns a list instead of a data.frame
#' 
#' @return a data.frame of gene expression results if \code{report.stats=F} or a list with entries \code{diff.exp} and \code{stats} if \code{report.stats=T}.
#' 
#' @export
binomTestAlongTree <- function(object, pseudotime, tips, log.effect.size=log(2), p.thresh=.01, frac.must.express=0.1, genes.use=NULL, root=NULL, only.return.global=F, must.beat.sibs=0.5, report.stats=F) {
  # Translate provided tip names into numbers if needed
  tips <- translateSegmentNames(object, tips)
  
  # Initialize
  current.tip <- tips
  tip.list <- c()
  markers <- list()
  cells.1.all <- c()
  cells.2.all <- c()
  stats <- data.frame(stringsAsFactors=F)
  
  # Loop through tree, collecting two sets of cells -- one on the trajectory from
  # population of interest to tip, and one pseudotime matched on sibling segments.
  while(length(current.tip) > 0) {
    # Find parent
    parent.to.tip <- segParent(object, current.tip)
    # Find other side(s) of branchpoints
    opposing.tips <- segSiblings(object, current.tip, include.self = F)
    # Compare each opposing branch separately.
    for (opposing.tip in opposing.tips) {
      segs.to.consider <- c(opposing.tip, segChildrenAll(object, opposing.tip))
      if (length(segs.to.consider) > 0) {
        tip.list <- c(tip.list, current.tip)
        # Get from both branches
        cells.1 <- unique(unlist(object@tree$cells.in.segment[current.tip]))
        cells.2 <- unique(unlist(object@tree$cells.in.segment[segs.to.consider]))
        # Limit them to the same pseudotime range
        pt.limits <- object@tree$segment.pseudotime.limits[current.tip,]
        cells.2 <- cells.2[which(object@pseudotime[cells.2,pseudotime] >= as.numeric(pt.limits[1]) & object@pseudotime[cells.2,pseudotime] <= as.numeric(pt.limits[2]))]
        # Add cells to the list
        cells.1.all <- c(cells.1.all, cells.1)
        cells.2.all <- c(cells.2.all, cells.2)
        # Do a differential expression test
        if (report.stats) {
          n.1 <- length(cells.1)
          n.2 <- length(cells.2)
          pt.1.mean <- mean(object@pseudotime[cells.1, pseudotime])
          pt.2.mean <- mean(object@pseudotime[cells.2, pseudotime])
          genes.1.mean <- mean(object@meta[cells.1, "n.Genes"])
          genes.2.mean <- mean(object@meta[cells.2, "n.Genes"])
          trans.1.mean <- mean(object@meta[cells.1, "n.Trans"])
          trans.2.mean <- mean(object@meta[cells.2, "n.Trans"])
          pt.1.median <- median(object@pseudotime[cells.1, pseudotime])
          pt.2.median <- median(object@pseudotime[cells.2, pseudotime])
          genes.1.median <- median(object@meta[cells.1, "n.Genes"])
          genes.2.median <- median(object@meta[cells.2, "n.Genes"])
          trans.1.median <- median(object@meta[cells.1, "n.Trans"])
          trans.2.median <- median(object@meta[cells.2, "n.Trans"])
          stats <- rbind(stats, c(n.1, n.2, pt.1.mean, pt.2.mean, pt.1.median, pt.2.median, genes.1.mean, genes.2.mean, genes.1.median, genes.2.median, trans.1.mean, trans.2.mean, trans.1.median, trans.2.median))
        }
        markers[[(length(markers)+1)]] <- markersBinom(object=object, cells.1=cells.1, cells.2=cells.2, genes.use=genes.use, effect.size=log.effect.size, p.thresh=p.thresh, frac.must.express=frac.must.express)
        names(markers)[length(markers)] <- paste0(current.tip, "-", opposing.tip)
      }
    }
    # Move upwards one branch
    current.tip <- parent.to.tip
  }
  
  # If a marker is negatively enriched in a *downstream* segment in the tree, then remove it as a marker.
  depleted.by.segment <- lapply(unique(tip.list), function(tip) {
    unique(unlist(lapply(markers[0:(min(which(tip.list == tip))-1)], function(markers.downstream) {
      rownames(markers.downstream)[markers.downstream$log.effect < 0 ]
    })))
  })
  names(depleted.by.segment) <- unique(tip.list)
  
  markers.2 <- lapply(1:length(tip.list), function(x) {
    tm <- markers[[x]]
    tm[setdiff(rownames(tm), depleted.by.segment[[tip.list[x]]]),]
  })
  names(markers.2) <- names(markers)
  
  # If user provides a root, get rid of the markers from it and earlier segments.
  if (!is.null(root) && root %in% tip.list) {
    keep.markers.until <- which(tip.list == root)-1
    markers.2 <- markers.2[1:keep.markers.until]
    tip.list <- tip.list[1:keep.markers.until]
  }
  
  # # Now, determine the remaining markers
  # if (must.beat.all.sibs) {
  #   markers.remain <- unique(unlist(lapply(unique(tip.list), function(tip) {
  #     Reduce(intersect, lapply(markers.2[which(tip.list == tip)], function(m) rownames(m)[m$log.effect > 0]))
  #   })))
  # } else {
  #   markers.remain <- unique(unlist(lapply(markers.2, function(m) rownames(m)[m$log.effect > 0])))
  # }

  # Now, determine how many siblings each marker beats and curate them.
  markers.3 <- markers.2
  for (vs in unique(tip.list)) {
    # Figure out which comparisons were for a group of sibs
    sibs.use <- which(tip.list==vs)
    # How many times was each marker found?
    sib.table <- table(unlist(lapply(markers.2[sibs.use], rownames)))
    # Which markers were found enough times?
    markers.keep <- names(which(sib.table >= (length(sibs.use) * must.beat.sibs)))
    # Curate individual tables
    for (this.sib in sibs.use) {
      markers.3[[this.sib]] <- markers.3[[this.sib]][intersect(rownames(markers.3[[this.sib]]), markers.keep),]
    }
  }
  markers.remain <- unique(unlist(lapply(markers.3, rownames)))
    
  # Summarize data by calculating across _all_ cells considered, and also max (by log.effect) in any branch.
  if (only.return.global) {
    markers.summary <- markersBinom(object=object, cells.1=cells.1.all, cells.2=cells.2.all, genes.use=markers.remain, frac.must.express = frac.must.express, effect.size=log.effect.size, p.thresh=p.thresh)
  } else {
    markers.summary <- markersBinom(object=object, cells.1=cells.1.all, cells.2=cells.2.all, genes.use=markers.remain, frac.must.express = 0, effect.size=0, p.thresh=1)
  }
  names(markers.summary) <- c("log.effect.all", "pval.all", "pval.fdr.all", "posFrac_lineage", "posFrac_rest", "nTrans_lineage", "nTrans_rest")
  markers.max <- lapply(rownames(markers.summary), function(marker) {
    id <- which.max(unlist(lapply(markers.3, function(x) x[marker,"log.effect"])))
    return(markers.3[[id]][marker,])  
  })
  markers.max <- do.call(what = "rbind", markers.max)
  names(markers.max) <- c("log.effect.maxbranch", "pval.maxbranch", "pval.fdr.maxbranch", "posFrac_maxBranch", "posFrac_opposingBranch", "nTrans_maxBranch", "nTrans_opposingBranch")
  markers.all <- cbind(markers.summary, markers.max)
  
  if (report.stats) {
    names(stats) <- c("n.1", "n.2", "pt.1.mean", "pt.2.mean", "pt.1.median", "pt.2.median", "genes.1.mean", "genes.2.mean", "genes.1.median", "genes.2.median", "trans.1.mean", "trans.2.mean", "trans.1.median", "trans.2.median")
    return(list(
      diff.exp=markers.all,
      stats=stats
    ))
  } else {
    return(markers.all)
  }
}

#' Test for differential gene expression along tree using precision-recall test
#' 
#' This test performs differential expression testing along URD's recovered tree
#' structure using the area under a precision-recall curve (see: \code{\link{markersAUCPR}})
#' to determine how good individual genes are as markers of a lineage.
#' 
#' It starts from the segment provided in \code{tips} and compares it to each of
#' its siblings (and their descendants) matched to the same pseudotime limits as
#' the segment under consideration. Markers that are upregulated in a segment compared
#' to either any of its siblings (if \code{must.beat.all.sibs=FALSE}) or compared
#' to each of its siblings (if \code{must.beat.all.sibs=TRUE}) become putative
#' markers.
#' 
#' 
#' @param object An URD object
#' @param pseudotime (Character) Name of pseudotime to use
#' @param tips (Character vector) The tip of the trajectory to find differential expression for (either as tip number or name)
#' @param log.effect.size (Numeric) Minimum fold-change difference in mean expression to be considered differential (default, \code{log(2)} is 2-fold change)
#' @param auc.factor (Numeric) The precision-recall AUC is determined for a random classifier is determined
#' based on the size of populations. To be considered differential, genes must have an AUC this factor
#' multiplied by the expected AUC of a random classifier.
#' @param max.auc.threshold (Numeric) This acts as an upper bound for how high the AUC must be for a gene to be considered differential.
#' @param frac.must.express (Numeric) Fraction of cells of interest that must express the gene to consider it differential.
#' @param frac.min.diff (Numeric) Minimum difference in fraction of cells expressing a gene to consider it differential.
#' @param genes.use (Character vector) Which genes to test for differential expression (default \code{NULL} is all genes)
#' @param root (Character) Which segment of the tree (by number) to use as the root for differential express (default \code{NULL} will autodetect, but this can be useful if you do not want to traverse up the entire tree)
#' @param segs.to.skip (Character vector) Are there any segments in the tree that should not be used for testing differential expression? (Default, \code{NULL} will test all segments.)
#' @param only.return.global (Logical) Genes must also pass an overall in-trajectory vs. out-of-trajectory test across all cells in the data to be considered differential (default, \code{FALSE})
#' @param must.beat.sibs (Numeric) At multi-furcated branchpoints, independent comparisons are made to each sibling. For a gene to be considered differential it must be a marker against this proportion of its sibling branches. (1 means it must be a marker against all other branches, 0 means it must be a marker against a single other branch.)
#' @param report.debug (Logical) If \code{TRUE}, this function returns a list instead of a data.frame, with 
#' \code{$stats} containing information about all of the comparisons? (nGene, nTrans, pseudotime, n.cells) and
#' \code{$marker.chain} containing the markers from each branchpoint.
#'
#' @return a data.frame of gene expression results if \code{report.debug=F} or a list with entries \code{diff.exp} and \code{stats} if \code{report.debug=T}.
#' 
#' @export
aucprTestAlongTree <- function(object, pseudotime, tips, log.effect.size=0.25, auc.factor=1, max.auc.threshold=1, frac.must.express=0.1, frac.min.diff=0.1, genes.use=NULL, root=NULL, segs.to.skip=NULL, only.return.global=F, must.beat.sibs=0.5, report.debug=F) {
  # Translate provided tip names into numbers if needed
  tips <- translateSegmentNames(object, tips)
  
  # Initialize
  current.tip <- tips
  tip.list <- c()
  markers <- list()
  anti.markers <- list()
  cells.1.all <- c()
  cells.2.all <- c()
  stats <- data.frame(stringsAsFactors=F)
  
  # Loop through tree, collecting two sets of cells -- one on the trajectory from
  # population of interest to tip, and one pseudotime matched on sibling segments.
  while(length(current.tip) > 0) {
    # Find parent
    parent.to.tip <- segParent(object, current.tip)
    # If the segment is no good, skip it.
    if (!(current.tip %in% segs.to.skip)) {
      # Find other side(s) of branchpoints
      opposing.tips <- segSiblings(object, current.tip, include.self = F)
      # Compare each opposing branch separately.
      for (opposing.tip in opposing.tips) {
        segs.to.consider <- c(opposing.tip, segChildrenAll(object, opposing.tip))
        if (length(segs.to.consider) > 0) {
          tip.list <- c(tip.list, current.tip)
          # Get cells from both branches
          cells.1 <- unique(unlist(object@tree$cells.in.segment[current.tip]))
          cells.2 <- unique(unlist(object@tree$cells.in.segment[segs.to.consider]))
          # Limit cell populations to the same pseudotime range
          pt.limits <- object@tree$segment.pseudotime.limits[current.tip,]
          cells.2 <- cells.2[which(object@pseudotime[cells.2,pseudotime] >= as.numeric(pt.limits[1]) & object@pseudotime[cells.2,pseudotime] <= as.numeric(pt.limits[2]))]
          # Add cells to the list keeping track of the globally compared cells
          cells.1.all <- c(cells.1.all, cells.1)
          cells.2.all <- c(cells.2.all, cells.2)
          # If keeping track of stats, add information to the data.frame
          if (report.debug) {
            n.1 <- length(cells.1)
            n.2 <- length(cells.2)
            pt.1.mean <- mean(object@pseudotime[cells.1, pseudotime])
            pt.2.mean <- mean(object@pseudotime[cells.2, pseudotime])
            genes.1.mean <- mean(object@meta[cells.1, "n.Genes"])
            genes.2.mean <- mean(object@meta[cells.2, "n.Genes"])
            trans.1.mean <- mean(object@meta[cells.1, "n.Trans"])
            trans.2.mean <- mean(object@meta[cells.2, "n.Trans"])
            pt.1.median <- median(object@pseudotime[cells.1, pseudotime])
            pt.2.median <- median(object@pseudotime[cells.2, pseudotime])
            genes.1.median <- median(object@meta[cells.1, "n.Genes"])
            genes.2.median <- median(object@meta[cells.2, "n.Genes"])
            trans.1.median <- median(object@meta[cells.1, "n.Trans"])
            trans.2.median <- median(object@meta[cells.2, "n.Trans"])
            stats <- rbind(stats, c(n.1, n.2, pt.1.mean, pt.2.mean, pt.1.median, pt.2.median, genes.1.mean, genes.2.mean, genes.1.median, genes.2.median, trans.1.mean, trans.2.mean, trans.1.median, trans.2.median))
          }
          # Test for markers of these cells with AUCPR test
          these.markers <- markersAUCPR(object=object, cells.1=cells.1, cells.2=cells.2, genes.use=genes.use, effect.size=log.effect.size, frac.must.express=frac.must.express, frac.min.diff=frac.min.diff)
          these.markers$AUCPR.ratio <- these.markers$AUCPR / aucprThreshold(cells.1, cells.2, factor=1, max.auc=Inf)
          # Since AUCPR is not symmetric, must test explicitly for markers of the other branch ("anti-marker")
          these.anti.markers <- markersAUCPR(object=object, cells.1=cells.2, cells.2=cells.1, genes.use=genes.use, effect.size=log.effect.size, frac.must.express=frac.must.express, frac.min.diff=frac.min.diff)
          # Limit markers to those that pass the minimum AUC and are positive markers of their respective branch
          markers[[(length(markers)+1)]] <- these.markers[these.markers$AUCPR >= aucprThreshold(cells.1, cells.2, factor = auc.factor, max.auc = max.auc.threshold) & these.markers$exp.fc > 0,]
          anti.markers[[(length(anti.markers)+1)]] <- these.anti.markers[these.anti.markers$AUCPR >= aucprThreshold(cells.2, cells.1, factor = auc.factor, max.auc = max.auc.threshold) & these.anti.markers$exp.fc > 0,]
          names(markers)[length(markers)] <- paste0(current.tip, "-", opposing.tip)
          names(anti.markers)[length(anti.markers)] <- paste0(current.tip, "-", opposing.tip)
        }
      }
    }
    # Move upwards one branch
    current.tip <- parent.to.tip
    # But, if a root is specified and you've gotten there, then stop.
    if (!is.null(root) && length(parent.to.tip) > 0 && parent.to.tip == root) {
        current.tip <- NULL
    }
  }
  
  # If a marker is negatively enriched in a *downstream* segment in the tree (i.e. is an anti-marker), then remove it as a marker.
  depleted.by.segment <- lapply(unique(tip.list), function(tip) {
    unique(unlist(lapply(anti.markers[0:(min(which(tip.list == tip))-1)], function(markers.downstream) {
      rownames(markers.downstream)
    })))
  })
  names(depleted.by.segment) <- unique(tip.list)
  
  markers.2 <- lapply(1:length(tip.list), function(x) {
    tm <- markers[[x]]
    tm[setdiff(rownames(tm), depleted.by.segment[[tip.list[x]]]),]
  })
  names(markers.2) <- names(markers)
  
  # If user provides a root, get rid of the markers from it and earlier segments.
  if (!is.null(root) && root %in% tip.list) {
    keep.markers.until <- which(tip.list == root)-1
    markers.2 <- markers.2[1:keep.markers.until]
    tip.list <- tip.list[1:keep.markers.until]
  }

  # Now, determine how many siblings each marker beats and curate them.
  markers.3 <- markers.2
  for (vs in unique(tip.list)) {
    # Figure out which comparisons were for a group of sibs
    sibs.use <- which(tip.list==vs)
    # How many times was each marker found?
    sib.table <- table(unlist(lapply(markers.2[sibs.use], rownames)))
    # Which markers were found enough times?
    markers.keep <- names(which(sib.table >= (length(sibs.use) * must.beat.sibs)))
    # Curate individual tables
    for (this.sib in sibs.use) {
      markers.3[[this.sib]] <- markers.3[[this.sib]][intersect(rownames(markers.3[[this.sib]]), markers.keep),]
    }
  }
  markers.remain <- unique(unlist(lapply(markers.3, rownames)))
  
  # Summarize data by calculating across _all_ cells considered, and also max and min passing (by AUCPR.ratio) in any branch.
  if (only.return.global) {
    markers.summary <- markersAUCPR(object=object, cells.1=cells.1.all, cells.2=cells.2.all, genes.use=markers.remain, frac.must.express = frac.must.express, effect.size=log.effect.size, frac.min.diff=frac.min.diff)
    markers.summary <- markers.summary[markers.summary$AUCPR >= aucprThreshold(cells.1.all, cells.2.all, factor = auc.factor, max.auc = max.auc.threshold),]
  } else {
    markers.summary <- markersAUCPR(object=object, cells.1=cells.1.all, cells.2=cells.2.all, genes.use=markers.remain, frac.must.express = 0, effect.size=0, frac.min.diff=0)
  }
  markers.summary$AUCPR.ratio <- markers.summary$AUCPR / aucprThreshold(cells.1.all, cells.2.all, factor=auc.factor, max.auc=Inf)
  names(markers.summary) <- c("AUCPR.all", "expfc.all", "posFrac_lineage", "posFrac_rest", "nTrans_lineage", "nTrans_rest", "AUCPR.ratio.all")
  markers.max.segment <- c()
  markers.max <- lapply(rownames(markers.summary), function(marker) {
    id <- which.max(unlist(lapply(markers.3, function(x) x[marker,"AUCPR.ratio"])))
    to.return <- markers.3[[id]][marker,]
    markers.max.segment <<- c(markers.max.segment, names(markers.3)[id])
    names(to.return) <- c("AUCPR_maxBranch", "expfc.maxBranch", "posFrac_maxBranch", "posFrac_opposingMaxBranch", "nTrans_maxBranch", "nTrans_opposingMaxBranch", "AUCPR.ratio.maxBranch")
    return(to.return)  
  })
  markers.max <- do.call(what = "rbind", markers.max)
  markers.min.segment <- c()
  markers.min <- lapply(rownames(markers.summary), function(marker) {
    id <- which.min(unlist(lapply(markers.3, function(x) x[marker,"AUCPR.ratio"])))
    to.return <- markers.3[[id]][marker,]
    markers.min.segment <<- c(markers.min.segment, names(markers.3)[id])
    names(to.return) <- c("AUCPR_minBranch", "expfc.minBranch", "posFrac_minBranch", "posFrac_opposingMinBranch", "nTrans_minBranch", "nTrans_opposingMinBranch", "AUCPR.ratio.minBranch")
    return(to.return)  
  })
  markers.min <- do.call(what = "rbind", markers.min)
  markers.all <- cbind(markers.summary, markers.max, markers.min)
  markers.all$segment.maxBranch <- markers.max.segment
  markers.all$segment.minBranch <- markers.min.segment
  
  if (report.debug) {
    names(stats) <- c("n.1", "n.2", "pt.1.mean", "pt.2.mean", "pt.1.median", "pt.2.median", "genes.1.mean", "genes.2.mean", "genes.1.median", "genes.2.median", "trans.1.mean", "trans.2.mean", "trans.1.median", "trans.2.median")
    return(list(
      diff.exp=markers.all,
      stats=stats, 
      marker.chain=markers
    ))
  } else {
    return(markers.all)
  }
}

#' Test for differential gene expression along tree using binomial test
#' 
#' @param object An URD object
#' @param cells.1 (Character vector) Cells in population of interest
#' @param cells.2 (List) Cells to compare against. If a list, each entry is a population to compare against, and markers must beat all populations.
#' @param label (Character) Label to use to split cells. Must be a column of \code{group.ids}
#' @param groups (List) List of label values to include in each group.
#' @param log.effect.size (Numeric) Minimum fold-change difference in mean expression to be considered differential (default, \code{log(2)} is 2-fold change)
#' @param auc.factor (Numeric) The precision-recall AUC is determined for a random classifier is determined
#' based on the size of populations. To be considered differential, genes must have an AUC this factor
#' multiplied by the expected AUC of a random classifier.
#' @param min.auc.threshold (Numeric) This acts as a lower bound for the AUC threshold, no matter how unbalanced the populations are.
#' @param max.auc.threshold (Numeric) This acts as an upper bound for how high the AUC must be for a gene to be considered differential.
#' @param frac.must.express (Numeric) Fraction of cells of interest that must express the gene to consider it differential.
#' @param frac.min.diff (Numeric) Minimum difference in fraction of cells expressing a gene to consider it differential.
#' @param genes.use (Character vector) Which genes to test for differential expression (default \code{NULL} is all genes)
#' @param min.groups.to.mark (Numeric) For how many of the groups must a gene be a marker in order to be considered differential?

#' @param report.debug (Logical) If \code{TRUE}, this function returns a list instead of a data.frame, with 
#' \code{$stats} containing information about all of the comparisons? (nGene, nTrans, pseudotime, n.cells) and
#' \code{$marker.chain} containing the markers from each branchpoint.
#'
#' @return a data.frame of gene expression results if \code{report.debug=F} or a list with entries \code{diff.exp} and \code{stats} if \code{report.debug=T}.
#' 
#' @export 
aucprTestByFactor <- function(object, cells.1, cells.2, label, groups, log.effect.size=0.25, auc.factor=1, min.auc.thresh=0.1, max.auc.thresh=Inf, frac.must.express=0.1, frac.min.diff=0, genes.use=NULL, min.groups.to.mark=1, report.debug=F) {
  # Divide cells up by grouping
  cells.byfac.1 <- lapply(groups, function(g) intersect(cells.1, cellsInCluster(object, label, g)))
  cells.byfac.2 <- lapply(groups, function(g) {
    lapply(cells.2, function(c) {
      intersect(c, cellsInCluster(object, label, g))
  })})
  
  # Initialize lists to hold results
    # Non-curated results of AUCPR test
    markers.by.grouping.raw <- rep(list(rep(list(NULL), length(cells.2))), length(groups))
    # Markers that pass AUCPR thresholds
    markers.by.grouping.pass <- rep(list(rep(list(NULL), length(cells.2))), length(groups))
    # Markers within each group that beat all cells.2 populations
    markers.by.grouping.allcells2 <- rep(list(NULL), length(groups))
    # Raw AUCPR threshold for a random classifier
    thresh.by.grouping <- rep(list(rep(list(NULL), length(cells.2))), length(groups))
    # Modified threshold (multiplied by factor, bounded by min/max) used to determine markers
    thresh.mod.by.grouping <- rep(list(rep(list(NULL), length(cells.2))), length(groups))
  
  # Compare 
  for (n in 1:length(groups)) {
    for (c in 1:length(cells.2)) {
      # Determine AUCPR markers for each grouping + cells.2
      markers <- markersAUCPR(object = object, cells.1 = cells.byfac.1[[n]], cells.2=cells.byfac.2[[n]][[c]], effect.size = log.effect.size, frac.must.express = frac.must.express, frac.min.diff = frac.min.diff, genes.use = genes.use)
      markers.by.grouping.raw[[n]][[c]] <- markers
      # Determine AUCPR threshold for random classifier, and threshold to use for calling markers
      thresh <- aucprThreshold(cells.1 = cells.byfac.1[[n]], cells.2=cells.byfac.2[[n]][[c]], factor=1, max.auc = Inf)
      thresh.by.grouping[[n]][[c]] <- thresh
      thresh.mod <- min(max(min.auc.thresh, thresh*auc.factor), max.auc.thresh)
      thresh.mod.by.grouping[[n]][[c]] <- thresh.mod
      # Determine which markers pass threshold
      markers.pass <- markers[markers$AUCPR >= thresh.mod & markers$exp.fc > 0,]
      markers.by.grouping.pass[[n]][[c]] <- markers.pass
    }
    
    # Determine which markers beat all cells.2
    markers.by.grouping.allcells2[[n]] <- Reduce(intersect, lapply(markers.by.grouping.pass[[n]], rownames))
  }

    # Determine markers that mark the required number of groups
    marker.table <- table(unlist(markers.by.grouping.allcells2))
    markers.enough.groups <- names(which(marker.table >= min.groups.to.mark))
    
    # Assemble differential expression summary to return
    markers.for.report <- unlist(markers.by.grouping.raw, recursive = F)
    mr.names <- paste0(rep(paste0("G", 1:length(groups)), each=length(cells.2)), rep(paste0("-C", 1:length(cells.2)), length(groups)))
    best.comparison <- rep(NA, length(markers.enough.groups))
    names(best.comparison) <- markers.enough.groups
    marker.report <- lapply(markers.enough.groups, function(m) {
      best.aucpr <- which.max(unlist(lapply(markers.for.report, function(g) {
        if (m %in% rownames(g)) {
          g[rownames(g)==m,'AUCPR']
        } else {
          NA
        }
      })))
      best.comparison[m] <<- mr.names[best.aucpr]
      return(markers.for.report[[best.aucpr]][m,])
    })
    marker.report <- do.call("rbind", marker.report)
    names(marker.report) <- c("AUCPR_best", "exp.fc_best", "posFrac_1", "posFrac_2", "nTrans_1", "nTrans_2")
    marker.report$comparison_best <- best.comparison
    marker.report$groups_marked <- marker.table[rownames(marker.report)]
    marker.report <- marker.report[order(marker.report$groups_marked, marker.report$AUCPR_best, decreasing = T),]

    # Return values
    if (!report.debug) {
      return(marker.report)
    } else {
      return(list(
        diff.exp=marker.report,
        by.grouping.raw=markers.by.grouping.raw,
        by.grouping.passthresh=markers.by.grouping.pass,
        by.grouping.allcells2=markers.by.grouping.allcells2,
        by.grouping.aucrandom=thresh.by.grouping,
        by.grouping.aucthresh=thresh.mod.by.grouping,
        groups.marked=marker.table
      ))
    }
}

#' Precision-Recall AUC threshold
#' 
#' The AUC for a random classifier for a precision-recall curve is based on the 
#' number of cells in each population. This function determines the random AUC,
#' multiplies it by an arbitrary factor (how much better than random must a
#' classifier be?), and gives it an upper cap.
#' 
#' @param cells.1 (Character vector) List of cells
#' @param cells.2 (Character vector) List of cells
#' @param factor (Numeric) Factor to multiply AUC of random classifier by to
#' determine threshold
#' @param max.auc (Numeric) Upper limit of AUC threshold to return
#' 
#' @return (Numeric) AUC threshold
#' 
#' @export
aucprThreshold <- function(cells.1, cells.2, factor, max.auc) {
  # Random is P / (P+N)
  random.auc <- length(cells.1) / (length(cells.1) + length(cells.2))
  # Multiply be factor
  min.auc <- random.auc * factor
  # Don't let it go above max.auc
  min.auc <- min(min.auc, max.auc)
  return(min.auc)
}

# Perform a t-test to test differential expression between two groups of cells.
#' t-Test for Differential Expression
#' 
#' @export
deTTest <- function(data, cells.1, cells.2, genes.use=NULL, p.thresh, effect.size) {
  if (is.null(genes.use)) genes.use <- rownames(data)
  de <- as.data.frame(t(as.data.frame(lapply(genes.use, function(gene) {
    t <- t.test(x=as.numeric(data[gene,cells.1]), as.numeric(data[gene,cells.2]))
    return(c(t$p.value, t$estimate[1], t$estimate[2]))
  }))))
  names(de) <- c("p", "mean.1", "mean.2")
  rownames(de) <- genes.use
  de$logfc <- log(de$mean.1/de$mean.2)
  de <- de[which(de$p <= p.thresh & abs(de$logfc) >= effect.size),]
  return(de)
}

#' Differential NMF module expression testing along tree
#' 
#' @param object An URD object
#' @param tips (Character) Tip of lineage to do differential expression for
#' @param data (data.frame or matrix) Module expression by cell: rows are modules, columns are cells.
#' @param genelist (List) ?
#' @param pseudotime (Character) Pseudotime to use (i.e. a column name of \code{@@pseudotime}) 
#' @param exclude.upstream (Logical)
#' @param effect.size (Numeric)
#' @param p.thresh (Numeric)
#' @param min.expression (Numeric)
#' @param root (Character) 
#' @param must.beat.all.sibs (Logical)
#' 
#' @export
moduleTestAlongTree <- function(object, tips, data, genelist, pseudotime, exclude.upstream=F, effect.size=log(2), p.thresh=.01, min.expression=0.1, root=NULL, must.beat.all.sibs=T) {
  # Translate provided tip names into numbers if needed
  tips <- translateSegmentNames(object, tips)
  
  # Initialize
  current.tip <- tips
  tip.list <- c()
  markers <- list()
  cells.1.all <- c()
  cells.2.all <- c()
  
  # Loop through tree, collecting two sets of cells -- one on the trajectory from
  # population of interest to tip, and one pseudotime matched on sibling segments.
  while(length(current.tip) > 0) {
    # Find parent
    parent.to.tip <- segParent(object, current.tip)
    # Find other side(s) of branchpoints
    opposing.tips <- segSiblings(object, current.tip, include.self = F)
    # Compare each opposing branch separately.
    for (opposing.tip in opposing.tips) {
      segs.to.consider <- c(opposing.tip, segChildrenAll(object, opposing.tip))
      if (length(segs.to.consider) > 0) {
        tip.list <- c(tip.list, current.tip)
        # Get from both branches
        cells.1 <- unique(unlist(object@tree$cells.in.segment[current.tip]))
        cells.2 <- unique(unlist(object@tree$cells.in.segment[segs.to.consider]))
        # Limit them to the same pseudotime range
        pt.limits <- object@tree$segment.pseudotime.limits[current.tip,]
        cells.2 <- cells.2[which(object@pseudotime[cells.2,pseudotime] >= as.numeric(pt.limits[1]) & object@pseudotime[cells.2,pseudotime] <= as.numeric(pt.limits[2]))]
        # Double check that they are cells in the data matrix
        cells.1 <- cells.1[cells.1 %in% colnames(data)]
        cells.2 <- cells.2[cells.2 %in% colnames(data)]
        # Add cells to the list
        cells.1.all <- c(cells.1.all, cells.1)
        cells.2.all <- c(cells.2.all, cells.2)
        # Do a differential expression test
        markers[[(length(markers)+1)]] <- deTTest(data=data, cells.1=cells.1, cells.2=cells.2, effect.size=effect.size, p.thresh=p.thresh)
        names(markers)[length(markers)] <- paste0(current.tip, "-", opposing.tip)
      }
    }
    # Move upwards one branch
    current.tip <- parent.to.tip
  }
  
  # If a marker is negatively enriched in a *downstream* segment in the tree, them remove it as a marker.
  mods.depleted <- unique(unlist(lapply(markers, function(y) rownames(y)[y$logfc < 0& y$mean.2 >= min.expression])))
  
  # If user provides a root, get rid of the markers from it and earlier segments.
  if (!is.null(root)) {
    keep.markers.until <- min(which(unlist(lapply(strsplit(names(markers), "-"), function(x) x[1]==root))))-1
    markers <- markers[1:keep.markers.until]
    tip.list <- unique(tip.list[1:(which(tip.list == root)-1)])
  }
  
  # Get all modules that are upregulated
  if (must.beat.all.sibs) {
    marker.seg.considered <- unlist(lapply(strsplit(names(markers), "-"), function(x) x[1]))
    mods.upregulated <- unique(unlist(lapply(unique(tip.list), function(seg) {
      up <- lapply(markers[which(marker.seg.considered == seg)], function(y) rownames(y)[y$logfc > 0 & y$mean.1 >= min.expression])
      return(Reduce(intersect, up))
    })))
  } else {
    mods.upregulated <- unique(unlist(lapply(markers, function(y) rownames(y)[y$logfc > 0 & y$mean.1 >= min.expression])))
  }
  
  # What genes 
  genes <- setdiff(unlist(genelist[mods.upregulated]), unlist(genelist[mods.depleted]))
  
  return(list(genes=genes, mods.up=mods.upregulated, mods.down=mods.depleted))
}
