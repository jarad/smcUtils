\name{renormalize.weights-deprecated}
\alias{renormalize.weights-deprecated}
\alias{renormalize.weights}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Renormalize weights}

\description{This function creates a discrete probability distribution from a 
set of unnormalized (log) weights.}

\usage{renormalize.weights(weights, log=FALSE, engine="R")}

\arguments{
  \item{weights}{a vector of numbers (>0 if log=FALSE)}
  \item{log}{a logical indicating whether the weights are given as the logarithm
  of the weights}
  \item{engine}{run using "R" or "C" code}      
}
\value{a vector of probabilities that sum to 1}

\author{Jarad Niemi}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{renormalize}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

\keyword{ htest }
