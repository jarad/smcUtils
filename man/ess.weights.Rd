\name{ess.weights}
\alias{ess.weights}
\alias{cov.weights}
\alias{ent.weights}

\title{Measures of weight non-uniformity}

\description{These functions calculate heuristic measures of the effective
number of samples given a set of weights (probabilities).}

\usage{
ess.weights(weights, engine="R")
cov.weights(weights, engine="R")
ent.weights(weights, engine="R")
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{weights}{a vector of weights (probabilities)}
  \item{engine}{run using "R" or "C" code}
}

\details{
`ess.weights' calculates the effective sample size, namely `1/(sum(weights^2))'.
ESS has a minimum of 1 and a maximum equal to `length(weights)' when weights 
are uniform.

`cov.weights' calculates the coefficient of variation of the weights, namely
`var(weights)/mean(weights)^2'. CoV has a minimum of 0 when weights are uniform 
and a maximum equal to `length(weights)'.

`ent.weights' calculates the entropy of the weights, namely 
`-sum(weights * log2(weights))'. Entropy has a minimum of 0 and a 
maximum equal to `log2(length(weights))' when weights are uniform. (For numerical 
stability, the log term is actually calculated with `weights+.Machine$double.eps',
which can cause the observed minimum to be less than 0.)
}

\value{
a scalar indicating how uniform the weights are 
}

\references{
Liu, J. (2004) _Monte Carlo Strategies in Scientific Computing_

Doucet, A., dr Freitas, N., and Gordon, N. (2001) _Sequential Monte Carlo 
Methods in Practice_
}

\author{
Jarad Niemi
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{resample}}
}

\examples{
ws = renormalize(runif(10))
ess.weights(ws)
cov.weights(ws)
ent.weights(ws)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

