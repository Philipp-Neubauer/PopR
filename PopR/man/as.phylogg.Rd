\name{as.phylogg}
\alias{as.phylogg}
\title{
converts output from exact linkage algorithm (Z) to ape 'phylo' format for tree visualization
}
\description{
Produce an ape 'phylo' format object for tree visualization
}
\usage{
as.phylogg(Z, N, tip.label = as.character(1:N))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Z}{
Output from exact linkage algorithm produced by \code{\link{elink.call}}
}
  \item{N}{
The number of tips (individual samples)
}
  \item{tip.label}{
tip.label are individual sample labels to be displayed on the tree
}
}

\value{
Phylogenetic tree object for use with the ape package.
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
\code{\link{elink.call}}
}
