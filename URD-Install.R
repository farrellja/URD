# Install Bioconductor requirements: Biobase, S4Vectors, AnnotationDbi, destiny
# Needs to check R version because installation method changed between R 3.4 and 3.5

if (as.numeric(R.version$major) <= 3 && as.numeric(R.version$minor) < 5) {
  # R < 3.5
  message("Loading Bioconductor (biocLite) to install required packages.")
  source("https://bioconductor.org/biocLite.R")
  message("Installing required packages from Bioconductor.")
  biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny'), suppressUpdates=T)
} else {
  # R >= 3.5
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("Installing biocManager to install required Bioconductor packages.")
    install.packages("BiocManager")
  }
  message("Installing required packages from Bioconductor (BiocManager).")
  BiocManager::install(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny'), suppressUpdates=T)
}

# Check that Bioconductor installation went smoothly.
if (!requireNamespace("Biobase", quietly = TRUE)) {stop("Failed to install required package 'Biobase' from Bioconductor.")}
if (!requireNamespace("S4Vectors", quietly = TRUE)) {stop("Failed to install required package 'S4Vectors' from Bioconductor.")}
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {stop("Failed to install required package 'AnnotationDbi' from Bioconductor.")}
if (!requireNamespace("destiny", quietly = TRUE)) {stop("Failed to install required package 'destiny' from Bioconductor.")}

# Install URD (and devtools, if necessary)
if (requireNamespace("devtools", quietly = TRUE)) {
  message("Installing URD")
  devtools::install_github(repo="farrellja/URD")
} else {
  message("Installing devtools")
  install.packages("devtools")
  message("Installing URD")
  devtools::install_github(repo="farrellja/URD")
}

# Check that URD installed.
if (requireNamespace("URD", quietly = TRUE)) {
  message("URD installed successfully!")
  message('You can load it by typing: library("URD")')
  message('Try "?URD" for starting tips.')
} else {
  message("Something went wrong. It doesn't seem that URD installed correctly.")
}