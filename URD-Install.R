message("Loading Bioconductor to install required packages.")
source("https://bioconductor.org/biocLite.R")
message("Installing required packages from Bioconductor.")
biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny'))
if (requireNamespace("devtools", quietly = TRUE)) {
  message("Installing URD")
  devtools::install_github(repo="farrellja/URD")
} else {
  message("Installing devtools")
  install.packages("devtools")
  message("Installing URD")
  devtools::install_github(repo="farrellja/URD")
}
