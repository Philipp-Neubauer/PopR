\name{DPM.call}
\alias{DPM.call}
\title{
Calls multiple Julia instances for MCMC inference in Dirichlet process models
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
DPM.call(datas = NULL, baseline = NULL, labels = NULL, learn = FALSE, iters = 1000, thin = 10, np = 2, typeof = "N", path.to.julia = getwd(), call_DPM_path = getwd())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{datas}{
the mixed sample - an n by p matrix of n individuals and p variables
}
  \item{baseline}{
the baseline (optional) - an K by p matrix of K individuals and p variables
}
  \item{labels}{
the baseline group (source) memberships (optional) - a K integer vector
}
  \item{learn}{
use a baseline? boolean TRUE or FALSE
}
  \item{iters}{
number of MCMC iterations
}
  \item{thin}{
thinning interval of the Markov Chain
}
  \item{np}{
Number of Processors - set to at least 2
}
  \item{typeof}{
Normal model inference is supported for now with option 'N'.
}
  \item{path.to.julia}{
provides the path to the julia executable - defaults to the working directory. If elsewhere, and not in the system PATH, the path needs to be provided here.
}
  \item{call_DPM_path}{
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
 \item{K_record}{
A vector giving the number of sources at each iteration.
}
 \item{alpha_record}{
 Samples from the posterior distribution of the Dirichlet process concentration parameter
}
 \item{k_0s}{
  Samples from the posterior distribution of the prior co-variance scale parameter k_0
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
\code{\link{elink.call}},\code{\link{as.phylogg}},\code{\link{MM.call}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
