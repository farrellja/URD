# Load Bioconductor
message("Loading Bioconductor to install required packages.")
source("https://bioconductor.org/biocLite.R")

# Install Bioconductor requirements: Biobase, S4Vectors, AnnotationDbi, destiny
message("Installing required packages from Bioconductor.")
biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny'), suppressUpdates=T)

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