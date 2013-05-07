\name{MM.call}
\alias{MM.call}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calls Julia to conduct BAyesian inference in finite mixtures
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
%%     ~~Describe \code{datas} here~~
}
  \item{baseline}{
%%     ~~Describe \code{baseline} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{conditional}{
%%     ~~Describe \code{conditional} here~~
}
  \item{iters}{
%%     ~~Describe \code{iters} here~~
}
  \item{thin}{
%%     ~~Describe \code{thin} here~~
}
  \item{np}{
%%     ~~Describe \code{np} here~~
}
  \item{typeof}{
%%     ~~Describe \code{typeof} here~~
}
  \item{path.to.julia}{
%%     ~~Describe \code{path.to.julia} here~~
}
  \item{call_MM_path}{
%%     ~~Describe \code{call_MM_path} here~~
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
function (datas = NULL, baseline = NULL, labels = NULL, conditional = FALSE, 
    iters = 1000, thin = 10, np = 1, typeof = "N", path.to.julia = getwd(), 
    call_MM_path = getwd()) 
{
    write.csv(file = "single_priors.csv", c(a.0, b.0, k.0, ak.0, 
        bk.0, v.0, mu.0))
    write.csv(file = "matrix_priors.csv", lambda.0)
    write.csv(file = "datas.csv", datas)
    cond = ifelse(conditional == F, 0, 1)
    write.csv(file = "baseline.csv", baseline)
    write.csv(file = "labels.csv", labels)
    if (.Platform$OS.type == "unix") {
        exec = file.path(path.to.julia, "./julia")
        command = paste(exec, "-p", np, file.path(call_MM_path, 
            "call_MM.jl"), cond, iters, thin, typeof, getwd())
        system(command)
    }
    else {
        exec = file.path(path.to.julia, "/julia/julia.bat")
        command = c("-p", np, file.path(call_MM_path, "call_MM.jl"), 
            cond, iters, thin, typeof, getwd())
        system2(exec, command)
    }
    class_ids <- read.csv("source_ids.csv", header = F)
    proportions <- read.csv("proportions.csv", header = F)
    post_probas <- read.csv("post_probas.csv", header = F)
    list(class_ids = class_ids, proportions = proportions, post_probas = post_probas)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
