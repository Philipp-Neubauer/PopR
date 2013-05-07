\name{get_prior_ab}
\alias{get_prior_ab}

\title{
matches the gamma parameters a and b to produce an implicit prior over the number of instances (sources) k that best resembles the explicit prior described by g 
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
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
Gamma shape parameter to be used as prior in \code{\link{call_DPM}}
}}
    \item{b}{
Gamma scale parameter to be used as prior in \code{\link{call_DPM}}
}}
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
\code{\link{DPM.call}}
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
