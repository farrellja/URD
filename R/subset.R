#' Subset an URD object
#' 
#' This subsets an URD object, given a list of cells. This subsets all results
#' without recalculating them, which could produce non-sensical values for many
#' things. (For instance, it would make more sense to recalculate the PCA on
#' the subset data, rather than keeping previous PCA scores for only a subset
#' of cells.)
#' 
#' @param object An URD object
#' @param cells.keep (Character vector) Names of cells to retain
#' 
#' @return An URD object
#' 
#' @examples 
#' object <- urdSubset(object, not.outlier.cells)
#' 
#' @export
urdSubset <- function(object, cells.keep) {
            # Make sure all cells are actually in the object
            cells.keep <- intersect(cells.keep, colnames(object@logupx.data))
            # Subset count.data if it's present
            if(!any(dim(object@count.data) == 0)) object@count.data <- object@count.data[,cells.keep]
            # Subset logupx data
            object@logupx.data <- object@logupx.data[,cells.keep]
            # Subset metadata
            object@meta <- object@meta[cells.keep,]
            # Subset groups
            object@group.ids <- object@group.ids[cells.keep, , drop=F]
            # Subset PCA if it's been calculated
            if(!any(dim(object@pca.scores) == 0)) object@pca.scores <- object@pca.scores[cells.keep,]
            # Subset tSNE
            if(!any(dim(object@tsne.y) == 0)) object@tsne.y <- object@tsne.y[cells.keep,]
            # Subset gene signatures
            if(!any(dim(object@gene.sig.z) == 0)) object@gene.sig.z <- object@gene.sig.z[cells.keep,]
            #if(!any(dim(object@gene.sig.p) == 0)) object@gene.sig.p <- object@gene.sig.p[cells.keep,]
            #if(!any(dim(object@gene.sig.raw) == 0)) object@gene.sig.raw <- object@gene.sig.raw[cells.keep,]
            # Subset diffusion map
            if (length(object@dm) > 0 && length(object@dm@eigenvectors) > 0) {
              ids.keep <- which(rownames(object@dm@eigenvectors) %in% cells.keep)
              object@dm@eigenvectors <- object@dm@eigenvectors[cells.keep,]
              object@dm@eigenvec0 <- object@dm@eigenvec0[ids.keep]
              object@dm@transitions <- object@dm@transitions[cells.keep,cells.keep]
              object@dm@d <- object@dm@d[ids.keep]
              object@dm@d_norm <- object@dm@d_norm[ids.keep]
            }
            # Subset pseudotime
            object@pseudotime <- object@pseudotime[cells.keep,,drop=F]
            if ("pseudotime" %in% names(object@pseudotime.stability)) object@pseudotime.stability$pseudotime <- object@pseudotime.stability$pseudotime[cells.keep,]
            if ("walks.per.cell" %in% names(object@pseudotime.stability)) object@pseudotime.stability$walks.per.cell <- object@pseudotime.stability$walks.per.cell[cells.keep,]
            # Subset diff.data
            if(!any(dim(object@diff.data) == 0)) object@diff.data <- object@diff.data[cells.keep,]
            # Tree????
            # Return cropped object
            return(object)
          }