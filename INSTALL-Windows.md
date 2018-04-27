# Installing URD (Windows):

## Quick installation

### 1. Install R and RStudio

*URD* is a package written in R and designed to be used in the RStudio interactive environment.

R can be obtained and installed from ([https://cran.rstudio.com](https://cran.rstudio.com)). 

Following installation of R, Rstudio can be obtained and installed from ([https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)). The free version of RStudio Desktop is sufficient.
        
### 2. Start RStudio

### 3. Install required R packages and URD

We wrote a script that will attempt to install all requirements for URD and then URD itself.

```source("https://raw.githubusercontent.com/farrellja/URD/master/URD-Install.R")```

## Troubleshooting: 

If something went wrong previously, you may have try installing some of URD's dependencies manually:

###### A. Install and attach the *devtools* package

*devtools* is required to compile and install URD, since it's distributed through GitHub.

```install.packages("devtools")```
     
###### B. Install required Bioconductor packages

Because these packages are deposited in Bioconductor and not CRAN, they must be installed manually before installing URD (otherwise its installation will fail.)

```source("https://bioconductor.org/biocLite.R")```
```biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny')```

Optionally, additional packages that are required for specific functions can be installed.

```biocLite(c('sva', 'rhdf5', 'scran'))```

If installing on a cluster, where you do not have write permissions to the main R libraries (for instance, if you receive errors like "ERROR: no permission to install to directory", you may need to specify a library location, such as:

```biocLite(c('sva', 'rhdf5', 'scran'), lib=[PATH TO DIRECTORY TO STORE R LIBRARIES WHERE YOU HAVE PERMISSION TO WRITE])```
     
###### C. Install URD

*URD* can be installed directly from the GitHub repository

```library(devtools)```
```install_github(repo="farrellja/URD")```

	