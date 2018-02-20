---
output:
  pdf_document: default
  html_document: default
---
# Installing URD

## 1. Install an X11 window client

The 3D plots require an X11 windowing client to display. This might already be installed on your computer if you've previously installed the *rgl* package.

#### For Mac OS X:

Install XQuartz ([https://www.xquartz.org](https://www.xquartz.org)), and then restart your computer.
        
## 2. Increase number of DLLs available to R

R has a hard-coded limit on the number of DLLs that can be loaded by linked packages. Currently, URD imports more packages than are 'allowed' by default. Increasing the number in the environment slightly increases the base memory footprint of R.

To increase the number of allowed DLLs from 100 to 150, the line ```R_MAX_NUM_DLLS=150``` must be included in ~/.Renviron

#### For Mac OS X:

From terminal, run:
```echo "R_MAX_NUM_DLLS=150" >> ~/.Renviron```
        
## 3. Start R

## 4. Install and attach the *devtools* package

*devtools* is required to compile and install URD, since it's currently a local source package.

```install.packages("devtools")```
     
```library(devtools)```
     
## 5. Install required Bioconductor packages

Because these packages are deposited in Bioconductor and not CRAN, they must be installed manually before installing URD (otherwise its installation will fail.)

```source("https://bioconductor.org/biocLite.R")```

```biocLite(c('Biobase', 'S4Vectors', 'AnnotationDbi', 'destiny', 'sva', 'rhdf5', 'scran'))```

If installing on a cluster, where you do not have write permissions to the main R libraries (for instance, if you receive errors like "ERROR: no permission to install to directory", you may need to specify a library location, such as:

```biocLite(c('destiny', 'sva', 'rhdf5', 'scran'), lib=[PATH TO DIRECTORY TO STORE R LIBRARIES WHERE YOU HAVE PERMISSION TO WRITE])```
     
## 6. Install URD

You may have to update the path to reflect where URD lives on your computer.

```install("~/Dropbox/URD")```
     
## 7. Report inevitable problems to Jeff

	