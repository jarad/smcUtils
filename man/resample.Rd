\name{resample}
\alias{resample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
resample(weights, num.samples = length(weights), method = c("multinomial", "residual", "stratified", "systematic", "branching"), nonuniformity = c("none", "ess", "cov", "entropy"), threshold = NULL, resample.function = multinomial.resample)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{num.samples}{
%%     ~~Describe \code{num.samples} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{nonuniformity}{
%%     ~~Describe \code{nonuniformity} here~~
}
  \item{threshold}{
%%     ~~Describe \code{threshold} here~~
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
function(weights, num.samples=length(weights), 
                    method = c("multinomial","residual","stratified","systematic","branching"),
                    nonuniformity = c("none","ess","cov","entropy"), threshold=NULL,
                    resample.function = multinomial.resample)
{
  method        <- match.arg(method)
  nonuniformity <- match.arg(nonuniformity)

  do.resample = FALSE
  switch(nonuniformity,
    "none"     = { do.resample = TRUE; },
    "ess"      = { if (is.null(threshold)) threshold=0.5*num.samples; 
                   if (ess.weights(weights)<threshold) do.resample = TRUE; },
    "cov"      = { if (is.null(threshold)) threshold=0.5*num.samples;
                   if (cov.weights(weights)>threshold) do.resample = TRUE },
    "entropy"  = { if (is.null(threshold)) threshold=0.5*log2(num.samples);
                   if (ent.weights(weights)<threshold) do.resample = TRUE }
  )

  if (do.resample) {
    switch(method,
      "multinomial" = { ids = multinomial.resample(weights, num.samples) },
      "residual"    = { ids = residual.resample(   weights, num.samples, resample.function) },
      "stratified"  = { ids = stratified.resample( weights, num.samples) },
      "systematic"  = { ids = systematic.resample( weights, num.samples) },
      "branching"   = { ids = branching.resample(  weights, num.samples) }
    )
    weights = rep(1/num.samples,num.samples)
  } else { 
    ids = 1:length(weights)
  }

  return(list(weights=weights,indices=ids))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
