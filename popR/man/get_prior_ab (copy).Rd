\name{get_prior_ab}
\alias{get_prior_ab}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
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
%%     ~~Describe \code{n} here~~
}
  \item{g}{
%%     ~~Describe \code{g} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{var}{
%%     ~~Describe \code{var} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (n, g = "uniform", mu = NULL, var = NULL) 
{
    if (!is.numeric(g)) {
        if (g == "uniform") {
            g <- rep(1/n, n)
        }
        else if (g == "poisson") {
            g <- dpois(1:n, mu)
        }
        else if (g == "negbin") {
            g <- dnbinom(1:n, mu, var)
        }
        else if (g == "lnorm") {
            g <- dlnorm(1:n, mu, var)
        }
        else if (g == "norm") {
            g <- dnorm(1:n, mu, var)
        }
        else {
            stop("Distribution need to be one of the following: uniform (default),poisson,negbin,lnorm,norm")
        }
    }
    g = g/sum(g)
    s_recursive <- function(n, k, s_rec) {
        if (k == 1 & n == 1) {
            return(1)
        }
        else if (k == 1) {
            return(-(n - 1) * s_rec[n - 1, k])
        }
        else {
            return(s_rec[n - 1, k - 1] - (n - 1) * s_rec[n - 
                1, k])
        }
    }
    snk = matrix(, n, n)
    for (i in 1:n) {
        for (j in 1:(i)) {
            if (i == j) 
                snk[i, j] = 1
            else snk[i, j] = (s_recursive(i, j, snk))
        }
    }
    snk = abs(snk)
    KL_dist <- function(par, others) {
        a <- par[1]
        b <- par[2]
        n <- others$n
        g <- others$g
        snk <- others$snk
        pi_fun <- function(alpha, otros) {
            a <- otros$a
            b <- otros$b
            n <- otros$n
            k <- otros$k
            exp((k + a - 1) * log(alpha) - b * alpha + lgamma(alpha) - 
                lgamma(alpha + n))
        }
        KL = 0
        pi_k = vector(, n)
        for (k in 1:n) {
            nint <- integrate(pi_fun, lower = 0, upper = Inf, 
                otros = list(a = a, b = b, k = k, n = n), stop.on.error = F, 
                abs.tol = 0)$value
            pi_k[k] <- ((b^a * snk[n, k])/gamma(a)) * nint
        }
        KL <- sum(g * log(g/pi_k))
        return(KL)
    }
    res <- optim(par = c(0.5, 0.5), KL_dist, others = list(n = n, 
        g = g, snk = snk), lower = c(0.001, 1e-04), upper = c(100, 
        100), method = "L-BFGS-B")
    a <- res$par[1]
    b <- res$par[2]
    pi_fun <- function(alpha, otros) {
        a <- otros$a
        b <- otros$b
        n <- otros$n
        k <- otros$k
        exp((k + a - 1) * log(alpha) - b * alpha + lgamma(alpha) - 
            lgamma(alpha + n))
    }
    pi_k = vector(, n)
    for (k in 1:n) {
        nint <- integrate(pi_fun, lower = 0, upper = Inf, otros = list(a = a, 
            b = b, k = k, n = n), stop.on.error = T, abs.tol = 0)$value
        pi_k[k] <- (((b^a) * snk[n, k])/gamma(a)) * nint
    }
    pi_k = pi_k/sum(pi_k)
    plot(1:n, pi_k, pch = 16, col = 4, xlab = "number of sources", 
        ylab = "prior probability mass")
    points(1:n, g, pch = 17, col = 3)
    legend(n - n/3, (max(c(pi_k, g))), pch = 16:17, legend = c("induced prior", 
        "original prior"), col = 4:3)
    list(a = a, b = b)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
