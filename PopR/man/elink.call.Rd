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
%%     ~~Describe \code{class.ids} here~~
}
  \item{path.to.julia}{
%%     ~~Describe \code{path.to.julia} here~~
}
  \item{elink_path}{
%%     ~~Describe \code{elink_path} here~~
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
function (class.ids, path.to.julia = getwd(), elink_path = getwd()) 
{
    write.csv(file = "class_ids.csv", class.ids)
    if (.Platform$OS.type == "unix") {
        exec = file.path(path.to.julia, "./julia")
    }
    else {
        exec = file.path(path.to.julia, "julia/julia.bat")
    }
    command = paste(exec, file.path(elink_path, "elink.jl"), 
        getwd())
    system(command)
    link <- read.csv("linkages.csv", header = F)
    list(tree = link)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
