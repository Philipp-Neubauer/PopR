\name{get_prior_ab}
\alias{get_prior_ab}

\title{
produce an implicit prior over the number of instances (sources) k
}
\description{
matches the gamma parameters a and b to produce an implicit prior over the number of instances (sources) k that best resembles the explicit prior described by g
}
\usage{
get_prior_ab(n, g = "uniform", mu = NULL, var = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
number of individuals in the mixed sample
}
  \item{g}{
explicit prior. Can be a vector of length n with probability mass at each 1:n, or can be one of 'uniform' (default), 'Poisson','negbin'(negative binomial with mean mu and rate var), 'lnorm' (log-normal),'norm'
}
  \item{mu}{
mean for the specified explicit prior distributions
}
  \item{var}{
variance/rate for the specified explicit prior distributions
}
}
\details{
Produces a plot of the explicit/original prior and implicit/induced prior.
}
\value{
  \item{a}{
Gamma shape parameter to be used as prior in \code{\link{DPM.call}}
}
    \item{b}{
Gamma scale parameter to be used as prior in \code{\link{DPM.call}}
}
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
\code{\link{DPM.call}}
}
