#' Plot Cumulative Distribution of Reads per Cell Barcode
#' 
#' To determine the number of cells recovered in a Dropseq sample, the cumulative
#' distribution of the number of reads recovered per unique cell barcode is plotted.
#' The elbow in the curve suggests the boundary between "real cells" with many reads
#' and "ambient RNA" that has bound to the remainder of the beads (which have few reads).
#' 
#' @param gz.in (Character) Path to Dropseq cell read counts report (".report_cell_readcounts.txt.gz")
#' @param n.cells (Numeric) Draw a line at putative number of cells to use downstream
#' @param main (Character) Title for the plot
#' @param xmax (Numeric) Max value of plot x-axis
#' @export
dsCutoffPlot <- function(gz.in, n.cells=NULL, main="", xmax=5000) {
  a=read.table(gz.in, header=F, stringsAsFactors=F)
  x=cumsum(a$V1)
  x=x/max(x)
  plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", main=main, xlim=c(1,xmax))
  if (!is.null(n.cells)) abline(v=n.cells, col="red", lty=2)
}

#' Read 10X Reads per Cell barcode
#' 
#' This reads and aggregates the molecule info file output by the 10X pathway
#' into reads per cell barcode for use in \code{\link{txCutoffPlot}}. (This is
#' because it takes some time to read in this data, so many plots can be tried
#' with only having to read the data once.)
#' 
#' @importFrom rhdf5 h5read
#' @importFrom stats aggregate
#' 
#' @param molecule.info.path (Character) Path to the molecule_info.h5 file output by 10X.
#' @return Returns a data.frame with columns "barcode" (2-bit encoded cell barcode), 
#' "reads", "unmapped_reads", and "nonconf_mapped_reads" (as per 10X documentation),
#' and "total_reads", which is the sum of the 3 previous columns.
#' @export
txReadsPerCell <- function(molecule.info.path) {
  # Read in the relevant info from the stupid HDF5 file
  mi <- data.frame(
    barcode=h5read(file=molecule.info.path, name="barcode", bit64conversion='bit64'),
    reads=h5read(file=molecule.info.path, name="reads"),
    unmapped_reads=h5read(file=molecule.info.path, name="unmapped_reads"),
    nonconf_mapped_reads=h5read(file=molecule.info.path, name="nonconf_mapped_reads"),
    stringsAsFactors=F
  )
  # Aggregate by barcode
  mi.agg <- aggregate(mi[,2:4], by=list(mi[,1]), FUN=sum)
  mi.agg$total_reads <- apply(mi.agg[,2:4], 1, sum)
  # Sort by total reads
  mi.agg <- mi.agg[order(mi.agg$total_reads, decreasing = T),]
  rm(mi)
  shh <- gc(verbose = F)
  names(mi.agg) <- c("barcode", "reads", "unmapped_reads", "nonconf_mapped_reads", "total_reads")
  # Convert 10X cell barcodes to strings
  mi.agg$barcode <- txBarcodesIntegerToChar(mi.agg$barcode)
  # Return
  return(mi.agg)
}

#' Convert 10X Integer-stored barcodes to character ones for subsetting matrices
#' 
#' Are you kidding me with this? How much space does this save, guys?
#' 
#' @importFrom bit64 as.bitstring
#' 
#' @param x (Integer64 vector) 10X cell barcodes, stored as 64-bit integers, as is molecule_info.h5
#' 
#' @return (Character vector) 10X cell barcodes, as DNA base characters
txBarcodesIntegerToChar <- function (x) {
  binary.code <- c("A","C","G","T")
  names(binary.code) <- c("00", "01", "10", "11")
  binary.barcodes <- bit64::as.bitstring(x)
  character.barcodes <- unlist(lapply(binary.barcodes, function(bb) {
    bb.chunk <- substring(bb, seq(1,nchar(bb),2), seq(2,nchar(bb),2))[17:32]
    return(paste0(binary.code[bb.chunk], collapse=""))
  }))
  return(character.barcodes)
}

#' Plot Cumulative Distribution of Reads per Cell Barcode
#' 
#' To determine the number of cells recovered in a Dropseq sample, the cumulative
#' distribution of the number of reads recovered per unique cell barcode is plotted.
#' The elbow in the curve suggests the boundary between "real cells" with many reads
#' and "ambient RNA" that has bound to the remainder of the beads (which have few reads).
#' The 10X pipeline determines which barcodes constitute 'cells' in a different fashion,
#' but this produces an equivalent 
#' 
#' @param gz.in (Character) Path to Dropseq cell read counts report (".report_cell_readcounts.txt.gz")
#' @param n.cells (Numeric) Draw a line at putative number of cells to use downstream
#' @param main (Character) Title for the plot
#' @param xmax (Numeric) Max value of plot x-axis
#' @export
txCutoffPlot <- function(tx.reads.per.cell=NULL, molecule.info.path=NULL, n.cells=NULL, main="", xmax=5000) {
  # If reads per cell object is not provided, load one.
  if (is.null(tx.reads.per.cell)) {
    if (is.null(molecule.info.path)) stop("Must provide either tx.reads.per.cell (output of txReadsPerCell) or molecule.info.path")
    tx.reads.per.cell <- txReadsPerCell(molecule.info.path)
  }
  x=cumsum(tx.reads.per.cell$total_reads)
  x=x/max(x)
  plot(1:length(x), x, type='l', col="blue", xlab="cell barcodes sorted by number of reads [descending]", ylab="cumulative fraction of reads", main=main, xlim=c(1,xmax))
  if (!is.null(n.cells)) abline(v=n.cells, col="red", lty=2)
}

#' Read cell barcode statistics produced by Dropseq pipeline
#' 
#' The Dropseq pipeline produces output summary files that describe the number of
#' genes, UMIs, and reads per cell barcode. This reads them in from several files
#' and collates them into a single data.table.
#' 
#' @importFrom data.table rbindlist
#' 
#' @param dge.summary.reports (Character vector) Paths of DGE report files. (These files typically end with the suffix ".report_dge.summary.txt")
#' @param cell.readcount.reports (Character vector) Paths of cell readcount reports for getting number of reads per cell.  (These files typically end with the suffix ".report_cell_readcounts.txt.gz")
#' @param sample.names (Character vector) Names of samples.
#' @param name.cells.with.sample (Logical) Convert cell barcodes to cell names?
#' 
#' @return Returns a data.frame with cells as rows, and metadata as columns.
#' 
#' @export
dsReadStats <- function(dge.summary.reports, cell.readcount.reports=NULL, sample.names, name.cells.with.sample=T) {
  dge.info <- data.table::rbindlist(lapply(1:length(sample.names), function(i) {
    dge.in <- read.table(dge.summary.reports[i], header=T, stringsAsFactors=F)
    if (!is.null(cell.readcount.reports)) {
      raw.in <- read.table(cell.readcount.reports[i], header=F, stringsAsFactors=F)
      names(raw.in) <- c("NUM_READS", "CELL_BARCODE")
      dge.in <- merge(dge.in, raw.in, all.x=TRUE, by="CELL_BARCODE")
    }
    dge.in$SAMPLE <- sample.names[i]
    return(dge.in)
  }), use.names=TRUE)
  if (name.cells.with.sample && !is.null(sample.names)) dge.info$CELL <- paste0(dge.info$SAMPLE, "_", dge.info$CELL_BARCODE)
  return(as.data.frame(dge.info))
}

#' Read Dropseq Digital Gene Expression tables
#' 
#' Reads digital gene expression tables output by the Dropseq pipeline.
#' @param dge.file (Character) Path to digital gene expression table
#' @param dge.name (Character) Name to append to all cell barcodes if multiple tables will be combined (default: NULL)
#' @return DGE: a sparse Matrix (dgCMatrix) with rows of genes and columns of cells.
#' @export
dsReadDGE <- function(dge.file, dge.name=NULL) {
  # Input raw DGE table
  this.dge <- read.table(file=dge.file, header=T, stringsAsFactors = F, sep="\t", quote="")
  # Change rownames to genes.
  rownames(this.dge) <- this.dge[,1]
  this.dge <- this.dge[,-1]
  # If given a name, append it to all cell barcodes.
  if (!is.null(dge.name)) names(this.dge) <- paste0(dge.name, "_", names(this.dge))
  this.dge <- as(as.matrix(this.dge), "dgCMatrix")
  return (this.dge)
}

#' Read in Digital Gene Expression Matrix from 10X Cellranger
#' 
#' @importFrom Matrix readMM
#' 
#' @param output.path (Character) Path that contains barcodes.tsv, genes.tsv, and matrix.mtx from Cellranger
#' @param dge.name (Character) Name to append to all cell barcodes if multiple tables will be combined (default: NULL)
#' @param cells.read (Character vector) Cell names to trim the matrix to, if desired (will be matched after strip.num, before appending dge.name)
#' @param name.by (Character: "gene" or "id") Whether to name genes by gene ID, or gene name. If gene name, expression of gene IDs with the same name will be summed.
#' @param strip.num (Logical) Cellranger appends a batch number (e.g. "-1") to cell barcodes. If TRUE, this strips it away.
#' 
#' @return A dgCMatrix, where rows are genes and columns are cells.
#' 
#' @export
txReadDGE <- function(output.path, dge.name=NULL, cells.read=NULL, name.by=c("gene", "id"), strip.num=T) {
  # Check arguments
  name.by <- tolower(name.by[1])
  if (!(name.by %in% c("gene", "id"))) stop("name.by must be either 'gene' or 'id'.")
  # Read the expression data
  tx.mtx <- Matrix::readMM(paste0(output.path, "/matrix.mtx"))
  # Read in the barcodes
  tx.barcodes <- read.table(paste0(output.path, "/barcodes.tsv"), stringsAsFactors = F)
  # Strip sample number if desired
  if (strip.num) {
    tx.barcodes <- unlist(lapply(strsplit(tx.barcodes$V1, "-"), function(x) x[1]))
  } else {
    tx.barcodes <- tx.barcodes$V1
  }
  # If going to trim cells, figure out which to keep
  cells.keep <- which(tx.barcodes %in% cells.read)
  # Add DGE name if desired
  if (!is.null(dge.name)) tx.barcodes <- paste0(dge.name, tx.barcodes)
  # Add cell names to expression data and trim if desired
  colnames(tx.mtx) <- tx.barcodes
  if (!is.null(cells.read)) tx.mtx <- tx.mtx[,cells.keep]
  # Read in the genes
  tx.genes <- read.table(paste0(output.path, "/genes.tsv"), stringsAsFactors = F, sep="\t")
  names(tx.genes) <- c("id", "gene")
  # Add gene names
  if (name.by == "id") {
    rownames(tx.mtx) <- tx.genes$id
    tx.mtx <- as(tx.mtx, "dgCMatrix")
  } else {
    tx.mtx.ag <- aggregate(as.matrix(tx.mtx), by=list(tx.genes$gene), FUN=sum)
    rownames(tx.mtx.ag) <- tx.mtx.ag[,1]
    tx.mtx.ag <- tx.mtx.ag[,-1]
    tx.mtx <- as(as.matrix(tx.mtx.ag), "dgCMatrix")
    rm(tx.mtx.ag)
    shh <- gc(verbose=F)
  }
  # Return it
  return(tx.mtx)
}

#' Fuse Dropseq Digital Gene Expression tables
#' 
#' Combines multiple DGEs with different dimensions, keeping all genes present
#' in any DGE, and correctly setting values to 0 if that gene was not observed
#' in another sample.
#' @param list.dge (List) List of DGEs
#' @return Sparse matrix (dgCMatrix) with rows as genes and columns as cells.
#' @seealso \code{\link{ds.read.dge}} to read in DGEs to provide in \code{list.dge}.
#' @export
dsCombineDGE <- function(list.dge) {
  all.gene.names <- sort(unique(unlist(lapply(list.dge, function(x) rownames(x)))))
  # Expand each data frame to include a bunch of NAs where a gene was missing
  for (i in 1:length(list.dge)) {
    genes.not.in.this.matrix <- setdiff(all.gene.names, rownames(list.dge[[i]]))
    if (length(genes.not.in.this.matrix) > 0) {
      add.matrix <- Matrix(data=rep(0, length(genes.not.in.this.matrix) * dim(list.dge[[i]])[2]), nrow=length(genes.not.in.this.matrix), sparse=T)
      rownames(add.matrix) <- genes.not.in.this.matrix
      list.dge[[i]] <- rbind(list.dge[[i]], add.matrix)
      list.dge[[i]] <- list.dge[[i]][all.gene.names,]
    }
  }
  # Now cbind all the data frames
  combined.dge <- do.call(cbind,list.dge)
  return(combined.dge)
}


#' Plot and trim Dropseq DGE tables by metadata
#' 
#' @param dge.table (Matrix or data.frame) Gene expression table
#' @param ds.meta (data.frame) Metadata 
#' @param cell.name (Character) Column name of \code{ds.meta} that corresponds to column names of \code{dge.table} (default: "CELL") 
#' @param meta.name (Character) Column name of \code{ds.meta} to use for filtering cells
#' @param meta.min (Numeric) Minimum acceptable value of metadata (default: NULL)
#' @param meta.max (Numeric) Maximum acceptable value of metadata (default: NULL)
#' @param title (Character) Histogram title (default: "UMIs per cell")
#' @param xlim (Numeric) Histogram x-axis limits (default: NULL)
#' @param breaks (Numeric) Histogram bins (default: 150)
#' @export
dsMetaTrim <- function(dge.table, ds.meta, cell.name="CELL", meta.name, meta.min=NULL, meta.max=NULL, title="UMIs per cell", x.lim=NULL, breaks=150) {
  # Extract the relevant metadata
  this.meta <- ds.meta[,meta.name]
  names(this.meta) <- ds.meta[,cell.name]
  # Keep track of which cells to exclude
  cells.keep <- rep(TRUE, length(this.meta))
  names(cells.keep) <- ds.meta[,cell.name]
  # Append newline to title
  title <- paste0(title, "\n")
  # If low threshold, calculate cells to exclude and track
  if (!is.null(meta.min)) {
    low.meta <- length(which(this.meta < meta.min)) / length(this.meta) * 100
    cells.keep[which(this.meta < meta.min)] <- FALSE
    title <- paste0(title, "Low: ", round(low.meta, 2), "%")
    if (!is.null(meta.max)) title <- paste0(title, " / ")
  }
  # If high threshold, calculate cells to exclude and track
  if (!is.null(meta.max)) {
    high.meta <- length(which(this.meta > meta.max)) / length(this.meta) * 100
    cells.keep[which(this.meta > meta.max)] <- FALSE
    title <- paste0(title, "High: ", round(high.meta, 2), "%")
  }
  # Plot histogram
  if (!is.null(x.lim)) {
    hist(this.meta, breaks=breaks, main=title, xlab="", xlim=x.lim)
  } else {
    hist(this.meta, breaks=breaks, main=title, xlab="")
  }
  # Add trimming lines
  if (!is.null(meta.min)) abline(v=meta.min, col='red', lty=2)
  if (!is.null(meta.max)) abline(v=meta.max, col='red', lty=2)
  # Trim cells
  dge.table <- dge.table[,cells.keep[colnames(dge.table)]]
  # Throw out genes that are 0.
  dge.table <- dge.table[apply(dge.table, 1, function(x) any(x != 0)),]
  return(dge.table)
}