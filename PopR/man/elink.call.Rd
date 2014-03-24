\name{elink.call}
\alias{elink.call}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Call exact linkage algorithm on a matrix of class memberships
}
\description{
Calls Julia to use the exact linkage algorithm on a matrix of class memberships as produced by iterative estiamtion of mixtures
}
\usage{
elink.call(class.ids, path.to.julia = getwd(), elink_path = system.file("exec", package = "PopR"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class.ids}{
A matrix of class memberships as produced by iterative estiamtion of mixtures, where columns are MCMC samples and rows are mixed sample individuals. Can be the $class_ids field of the output produced by \code{\link{DPM.call}}
}
  \item{path.to.julia}{
provides the path to the julia executable - defaults to the working directory
}
  \item{elink_path}{
provides the path to the julia source script - defaults to the exec directory in the package root of PopR
}
}

\value{
A canonical representation of the tree, as produced used by matlabs tree plotting functions. Can be transformed into a phylo type object using \code{\link{as.phylogg}}
}
\references{
P. Neubauer, J.S. Shima and S.E. Swearer - Inferring dispersal and migrations from incomplete geochemical baselines: analysis of population structure using Bayesian infinite mixture models. DOI: 10.1111/2041-210X.12076

\url{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12076/abstract}
}
\author{
Philipp Neubauer
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{as.phylogg}},\code{\link{DPM.call}},\code{\link{MM.call}}
}
