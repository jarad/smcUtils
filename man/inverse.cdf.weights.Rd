\name{inverse.cdf.weights}
\alias{inverse.cdf.weights}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Compute weights by inverse cdf}

\description{Indices given by the inverse cdf method}

\usage{
inverse.cdf.weights(weights, uniforms = runif(length(weights)))
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{weights}{a vector of normalized probabilities}
  \item{uniforms}{a vector of uniform (0,1) random variables}
}

\details{
For each uniform u, returns the index j such that 
`cusum(weights)[j-1] < u < cusum(weights)[j]'.
}

\value{
  a vector of indices for the replicated particles
}

\author{Jarad Niemi}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{resampling}}
}
\examples{
ws = renormalize.weights(runif(6))
inverse.cdf.weights(ws)
inverse.cdf.weights(ws,runif(9,0,0.2))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
