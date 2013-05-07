\name{DPM.call}
\alias{DPM.call}
Calls multiple Julia instances for MCMC inference in Dirichlet process models
\title{
%%  ~~function to do ... ~~
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
%%     ~~Describe \code{datas} here~~
}
  \item{baseline}{
%%     ~~Describe \code{baseline} here~~
}
  \item{labels}{
%%     ~~Describe \code{labels} here~~
}
  \item{learn}{
%%     ~~Describe \code{learn} here~~
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
  \item{call_DPM_path}{
%%     ~~Describe \code{call_DPM_path} here~~
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
function (datas = NULL, baseline = NULL, labels = NULL, learn = FALSE, 
    iters = 1000, thin = 10, np = 1, typeof = "N", path.to.julia = getwd(), 
    call_DPM_path = getwd()) 
{
    write.csv(file = "single_priors.csv", c(a.0, b.0, k.0, ak.0, 
        bk.0, v.0, mu.0))
    write.csv(file = "matrix_priors.csv", lambda.0)
    write.csv(file = "datas.csv", datas)
    bl = ifelse(learn == F, 0, 1)
    if (learn) {
        write.csv(file = "baseline.csv", baseline)
        write.csv(file = "labels.csv", labels)
    }
    if (.Platform$OS.type == "unix") {
        exec = file.path(path.to.julia, "./julia")
        command = paste(exec, "-p", np, file.path(call_DPM_path, 
            "call_DPM.jl"), bl, iters, thin, typeof, getwd())
        system(command)
    }
    else {
        exec = file.path(path.to.julia, "/julia/julia.bat")
        command = c("-p", np, file.path(call_DPM_path, "call_DPM.jl"), 
            bl, iters, thin, typeof, getwd())
        system2(exec, command)
    }
    class_ids <- read.csv("source_ids.csv", header = F)
    K_record <- read.csv("K_record.csv", header = F)
    alpha_record <- read.csv("gammas.csv", header = F)
    k_0s <- read.csv("k_0s.csv", header = F)
    list(class_ids = class_ids, K_record = K_record, alpha_record = alpha_record, 
        k_0s = k_0s)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
