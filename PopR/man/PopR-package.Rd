\name{PopR-package}
\alias{PopR-package}
\alias{PopR}
\docType{package}
\title{
Analysis of population structure in R using finite and infinite Bayesian mixture models
}
\description{
The package implements R wrappers for the analysis of population structure using Dirichlet process models coded in Julia. Julia is a recent open source, high performance language for technical computing. Please follow the instructions in the package url to install Julia and to set up necessary packages in Julia. Post processing is done by transforming the output with the exact linkage algorithm implemented in Julia as well, with plotting functions from the ape package providing tree plotting facilities. 
}
\details{
\tabular{ll}{
Package: \tab PopR\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-05-07\cr
License: \tab GPL (>= 2)\cr
}

}
\author{
Philipp Neubauer

Maintainer: Philipp Neubauer <neubauer.phil@gmail.com>

}
\references{
P. Neubauer, J.S. Shima and S.E. Swearer - Inferring dispersal and migrations from incomplete geochemical baselines: analysis of population structure using Bayesian infinite mixture models. DOI: 10.1111/2041-210X.12076

\url{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12076/abstract}
}

\keyword{ package }
\seealso{
\code{\link{elink.call}},\code{\link{as.phylogg}},\code{\link{DPM.call}},\code{\link{MM.call}}
}
