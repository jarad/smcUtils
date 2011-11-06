\name{rep2id}
\alias{rep2id}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{rep2id}

\description{Turns a vector of replicate numbers into a vector }

\usage{rep2id(rep)}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rep}{a vector giving the number of replicates for each component
%%     ~~Describe \code{rep} here~~
}
}

\value{
a vector giving indices
}
\author{
Jarad Niemi
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{resampling}}
}
\examples{
# The goal here is to get functionality that might be expected from
# rep(1:4,each=1:4) # but the rep function only takes the first value in `each'
all.equal(rep(1:4,each=1:4),rep2id(rep(1,4)))
rep2id(1:4)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ internal }
