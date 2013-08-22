\name{get_Pc}
\alias{get_Pc}

\title{
Get the co-assignment probability of two (sets) of individuals from a co-assignment tree
}
\description{
Get the co-assignment probability of two (sets) of individuals as given by their index in the dataset from a co-assignment tree
}
\usage{
get_Pc(inds1, inds2, phylo)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{inds1}{
The index in the mixed sample of an individual or set of individuals
}
  \item{inds2}{
The index in the mixed sample of an individual or set of individuals
}
  \item{phylo}{
A co-assignment tree as produced by \code{\link{as.phylogg}}
}
}
\details{
returns the co-assignment probability of ind1 and ind2
}
\value{
Pc, the co-assignment probability of ind1 and ind2
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
\code{\link{as.phylogg}}
}
