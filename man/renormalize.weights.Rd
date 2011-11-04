\name{renormalize.weights}
\alias{renormalize.weights}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Renormalize weights}
\description{This function creates a discrete probability distribution from a 
set of unnormalized (log) weights.}
\usage{renormalize.weights(weights,log=FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{weights}{a vector of non-negative numbers}
  \item{log}{a logical indicating whether the weights are given as the logarithm
  of the weights}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{a vector of probabilities that sum to 1}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Jarad Niemi
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
\code{\link{resampling}}
}
\examples{
ws = runif(10)
renormalize.weights(ws)

log.ws = dnorm(rnorm(10,1,0),log=TRUE)
renormalize.weights(log.ws,log=TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
