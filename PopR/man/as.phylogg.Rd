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
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Z, N, tip.label = as.character(1:N)) 
{
    Z[Z > N] = Z[Z > N] + 1
    Z <- rbind(Z, c(N + 1, max(Z) + 1, 1))
    edgess <- (matrix(NA, max(Z), 2))
    edgel <- (matrix(NA, max(Z), 1))
    N = as.integer(N + 1)
    z = cbind(Z[, c(2, 1)], rep(NA, nrow(Z)))
    z[, 3] = (N + 1):(max(Z) + 1)
    l = Z[, 3]
    cnt = 1
    for (i in 1:nrow(z)) {
        if (all(z[i, 1:2] <= N)) {
            edgess[cnt, ] = as.integer(z[i, c(3, 1)])
            edgel[cnt, ] = l[i]
            edgess[cnt + 1, ] = as.integer(z[i, c(3, 2)])
            edgel[cnt + 1, ] = l[i]
            cnt = cnt + 2
        }
        else if (any(z[i, 1:2] <= N)) {
            edgess[cnt, ] = as.integer(z[i, c(3, 1)])
            edgel[cnt, ] = l[i] - l[z[, 3] == z[i, 1]]
            edgess[cnt + 1, ] = as.integer(z[i, c(3, 2)])
            edgel[cnt + 1, ] = l[i]
            cnt = cnt + 2
        }
        else {
            edgess[cnt, ] = as.integer(z[i, c(3, 1)])
            edgel[cnt, ] = l[i] - l[z[, 3] == z[i, 1]]
            edgess[cnt + 1, ] = as.integer(z[i, c(3, 2)])
            edgel[cnt + 1, ] = l[i] - l[z[, 3] == z[i, 2]]
            cnt = cnt + 2
        }
    }
    edgess[edgess[, 1] > N, 1] = -1L * (edgess[edgess[, 1] > 
        N, 1] - max(edgess)) + N + 1L
    edgess[edgess[, 2] > N, 2] = -1L * (edgess[edgess[, 2] > 
        N, 2] - max(edgess)) + N + 1L
    Zp <- list(edge = edgess, edge.length = edgel, tip.label = c(tip.label, 
        "OUT"), Nnode = N - 1L)
    class(Zp) <- "phylo"
    return(Zp)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
