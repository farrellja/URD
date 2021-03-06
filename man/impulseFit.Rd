% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impulse.R
\name{impulseFit}
\alias{impulseFit}
\title{Fit gene expression data with an impulse model}
\usage{
impulseFit(x, y, k = 50, interpolate = NULL, pulse.only = T, min.slope = 0.1)
}
\arguments{
\item{x}{(Numeric) Expression data to fit (e.g. pseudotime)}

\item{y}{(Numeric) Expression data to fit (e.g. expression level)}

\item{k}{(Numeric) Number of sets of starting conditions to try}

\item{interpolate}{(Numeric or NULL) If low number of data points, can interpolate them linearly to this number of points for choosing potential starting conditions. Default (\code{NULL}) does not do any interpolation. Interpolation is NOT used during actual function fitting.}

\item{pulse.only}{(Logical) If \code{TRUE}, filters out double sigmoid functions that are monotonous and prefers either single sigmoid or linear fit for them instead.}

\item{min.slope}{(Numeric) For \code{pulse.only} filtering, first derivative of double sigmoid function must include positive and negative values with absolute value > \code{min.slope}.}
}
\value{
(Named List, length 4) Parameters & information about best fitting model:
\itemize{ 
  \item{\strong{\code{type}}: 0/1/2 = linear/single sigmoid/double sigmoid}
  \item{\strong{\code{time.on}}: Gene onset pseudotimes}
  \item{\strong{\code{time.off}}: Gene offset pseudotimes}
  \item{\strong{\code{model}}: Parameters of the best model}
}
}
\description{
This fits either a linear, single sigmoid, or double sigmoid ("impulse") model to given x and y coordinates (e.g. expression data). This function is primarily used by \code{\link{geneCascadeProcess}}.
}
\details{
Thanks to Yiqun Wang for considerable improvements to the impulse fitting functions.
}
\references{
Chechik G, Oh E, Rando O, Weissman J, Regev A, Koller D. "Activity motifs reveal principles of timing in transcriptional control of the yeast metabolic network." Nat Biotechnol. 2008 Nov;26(11):1251-9. doi: 10.1038/nbt.1499
}
