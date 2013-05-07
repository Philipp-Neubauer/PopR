\name{MM.call}
\alias{MM.call}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calls Julia to conduct Bayesian inference in finite mixtures
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MM.call(datas = NULL, baseline = NULL, labels = NULL, conditional = FALSE, iters = 1000, thin = 10, np = 1, typeof = "N", path.to.julia = getwd(), call_MM_path = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datas}{
the mixed sample - an n by p matrix of n individuals and p variables
}
  \item{baseline}{
the baseline - an K by p matrix of K individuals and p variables
}
  \item{labels}{
the baseline group (source) memberships (optional) - a K integer vector
}
  \item{conditional}{
Should the analysis be conditional? TRUE or FALSE (default)
}
  \item{iters}{
number of MCMC iterations
}
  \item{thin}{
thinning interval of the Markov Chain
}
  \item{np}{
depricated
}
  \item{typeof}{
Normal model inference is supported for now with option 'N'.
}
  \item{path.to.julia}{
provides the path to the julia executable - defaults to the working directory. If elsewhere, and not in the system PATH, the path needs to be provided here.
}
  \item{call_MM_path}{
provides the path to the julia source script - defaults to the exec directory in the package root of PopR (where it should be)
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with fields:
 \item{class_ids}{
A matrix of class memberships where columns are MCMC samples and rows are mixed sample individuals. 
}
 \item{proportions}{
 Samples from the posterior distribution of class proportions at each MCMC iteration.
}
 \item{post_probas}{
 Samples from the posterior distribution of posterior probabilities of class memebership for each individual
}

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
\code{\link{elink.call}},\code{\link{as.phylogg}},\code{\link{DPM.call}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
