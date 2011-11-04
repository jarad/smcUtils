\name{inverse.cdf.weights}
\alias{inverse.cdf.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
inverse.cdf.weights(weights, uniforms = runif(length(weights)))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{
%%     ~~Describe \code{weights} here~~
}
  \item{uniforms}{
%%     ~~Describe \code{uniforms} here~~
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
function(weights,uniforms=runif(length(weights))) {
  if (is.unsorted(uniforms)) uniforms = sort(uniforms)

  num.samples   = length(uniforms)
  ids           = integer(num.samples)
  weights.cusum = cumsum(weights)

  index = as.integer(1)
  one   = as.integer(1)
  for (i in 1:num.samples) {
    found = FALSE
    while (!found) {
      if (uniforms[i] > weights.cusum[index]) {
        index = index+one
      } else {
        found = TRUE
      }
    }
    ids[i] = index
  }
  return(ids)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
