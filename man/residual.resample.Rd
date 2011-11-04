\name{residual.resample}
\alias{residual.resample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
residual.resample(weights, num.samples = length(weights), resample.function = multinomial.resample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{num.samples}{
%%     ~~Describe \code{num.samples} here~~
}
  \item{resample.function}{
%%     ~~Describe \code{resample.function} here~~
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
function(weights,num.samples=length(weights),resample.function=multinomial.resample) {
  expected.num.samples = num.samples*weights
  deterministic.reps   = floor(expected.num.samples)
  sum.reps             = sum(deterministic.reps)
  residual             = num.samples-sum.reps

  deterministic.ids = rep2id(deterministic.reps)
  random.ids        = resample.function((expected.num.samples-deterministic.reps)/residual,residual)

  return(c(deterministic.ids,random.ids))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
