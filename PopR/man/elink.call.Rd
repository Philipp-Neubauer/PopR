\name{elink.call}
\alias{elink.call}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calls Julia to use the exact linkage algorithm on a matrix of class memberships as produced by iterative estiamtion of mixtures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
elink.call(class.ids, path.to.julia = getwd(), elink_path = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{class.ids}{
A matrix of class memberships as produced by iterative estiamtion of mixtures, where columns are MCMC samples and rows are mixed sample individuals. Can be the $class_ids field of the output produced by \code{\link{call_DPM}}
}
  \item{path.to.julia}{
provides the path to the julia executable - defaults to the working directory
}
  \item{elink_path}{
provides the path to the julia source script - defaults to the exec directory in the package root of PopR
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A canonical representation of the tree, as produced used by matlabs tree plotting functions. Can be transformed into a phylo type object using \code{\link{as.phylogg}}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Philipp Neubauer
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{as.phylogg}},\code{\link{DPM.call}},\code{\link{MM.call}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
