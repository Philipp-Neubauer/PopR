\name{MM.call}
\alias{MM.call}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calls Julia to conduct Bayesian inference in finite mixtures
}
\description{
A wrapper function that outputs .csv files with parameter values and calls julia to read these in and perform MCMC for the finite mixture model
}
\usage{
MM.call(datas = NULL, baseline = NULL, labels = NULL, conditional = FALSE, iters = 1000, thin = 10, np = 1, typeof = "N", path.to.julia = getwd(), call_MM_path = system.file("exec", package = "PopR"), v.0 =
                 NULL, lambda.0 = NULL)
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
\item{v.0}{
prior degrees of freedom nu_0 (certainty) for lambda_0 - needs to be at least p+1, where p is the dimension of the data (this is the least informative choice)
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
 \item{proportions}{
 Samples from the posterior distribution of class proportions at each MCMC iteration.
}
 \item{post_probas}{
 Samples from the posterior distribution of posterior probabilities of class memebership for each individual
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
\code{\link{elink.call}},\code{\link{as.phylogg}},\code{\link{DPM.call}}
}
