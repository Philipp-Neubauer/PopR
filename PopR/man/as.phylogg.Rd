\name{as.phylogg}
\alias{as.phylogg}
\title{
converts output from exact linkage algorithm (Z) to ape 'phylo' format for tree visualization
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Phylogenetic tree object for use with the ape package.
}
\references{
TBA
}
\author{
Philipp Neubauer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{elink.call}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
