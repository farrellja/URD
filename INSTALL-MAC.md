# Installing URD

## For Mac OS X:

### 1. Install R and RStudio

*URD* is a package written in R and designed to be used in the RStudio interactive environment.

R can be obtained and installed from ([https://cran.rstudio.com](https://cran.rstudio.com)). 

Following installation of R, Rstudio can be obtained and installed from ([https://www.rstudio.com/products/rstudio/download/](https://www.rstudio.com/products/rstudio/download/)). The free version of RStudio Desktop is sufficient.

### 2. Install an X11 window client

3D plots produced by URD require an X11 windowing client to display. Our favorite for Mac is XQuartz ([https://www.xquartz.org](https://www.xquartz.org)). After installation, restart your computer.

### 3. Increase number of DLLs available to R

R has limit on the number of DLLs that can be loaded by linked packages. Currently, URD's dependencies require more DLLs than is allowed in a default installation of R. To increase the number of simultaneously allowed DLLs from 100 to 150, the *.Renviron* file must be modified.

From terminal, run:
```echo "R_MAX_NUM_DLLS=150" >> ~/.Renviron```
        
### 4. Start RStudio

### 5. Install required R packages and URD

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

	