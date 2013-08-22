\name{DPM.call}
\alias{DPM.call}
\title{
Calls multiple Julia instances for MCMC inference in Dirichlet process models
}
\description{
A wrapper function that outputs .csv files with parameter values and calls julia to read these in and perform MCMC for the DPM model
}
\usage{
DPM.call(datas = NULL, baseline = NULL, labels = NULL, learn = FALSE, iters = 1000, thin = 10, np = 1, typeof = "N", path.to.julia = getwd(), call_DPM_path = system.file("exec", package = "PopR"), a.0 = NULL, b.0
                 = NULL, ak.0 = NULL, bk.0 = NULL, v.0 =
                 NULL, lambda.0 = NULL)
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
\item{a.0}{
prior shape for the gamma prior for the Dirichlet concentration parameter, as produced by \code{\link{get_prior_ab}}
}
\item{b.0}{
prior scale for the gamma prior for the Dirichlet concentration parameter, as produced by \code{\link{get_prior_ab}}
}
\item{ak.0}{
prior shape for the gamma prior on kappa_0
}
\item{bk.0}{
prior scale for the gamma prior on kappa_0
}
\item{v.0}{
prior degrees of freedom v_0 (certainty) for lambda_0 - needs to be at least p+1, where p is the dimension of the data (this is the least informative choice)
}
\item{lambda.0}{
prior guess at the within source co-variance co-variance matrix. Best not fixed at extreme values but from sensible guesses from the data, or calculated from data.
}
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
P. Neubauer, J.S. Shima and S.E. Swearer - Inferring dispersal and migrations from incomplete geochemical baselines: analysis of population structure using Bayesian infinite mixture models. DOI: 10.1111/2041-210X.12076

\url{http://onlinelibrary.wiley.com/doi/10.1111/2041-210X.12076/abstract}
}
\author{
Philipp Neubauer
}

\seealso{
\code{\link{elink.call}},\code{\link{as.phylogg}},\code{\link{MM.call}}
}

