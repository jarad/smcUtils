\name{Unbiased resampling}
\alias{resampling}
\alias{multinomial.resample}  
\alias{residual.resample}
\alias{stratified.resample}
\alias{systematic.resample}
\alias{branching.resample}

\title{Resampling functions}

\description{A set of resampling functions with unbiased number of replicates.}

\usage{
multinomial.resample(weights, num.samples = length(weights), engine="R")
residual.resample(weights, num.samples = length(weights), engine="R", rrf = "stratified")
stratified.resample( weights, num.samples = length(weights), engine="R")
systematic.resample( weights, num.samples = length(weights), engine="R")
branching.resample(  weights, num.samples = length(weights), engine="R")
}

\arguments{
  \item{weights}{a vector of normalized weights}
  \item{num.samples}{the number of samples to return (for `branching.resample', 
  `num.samples' is the expected number of samples as the actual number is random)}   
  \item{rrf}{for residual resampling, the resampling function to use on the residual}
  \item{engine}{run using "R" or "C" code}      
}

\details{
'multinomial.resample' samples component i with probability `weights[i]', 
repeats this sampling `num.samples' times, and returns indices for the sampled
components.

'residual.resample' deterministically copies `floor(weights)' number of 
each component and then performs `rrf' on the remainder.

`stratified.resample' draws `num.samples' uniform random variables on the 
((i-1)/num.samples,i/num.samples) intervals of (0,1). It then uses the
inverse.cdf.weights function to determine which components to sample.

`systematic.resample' draws 1 uniform random variable on (0,1/num.samples),
builds a sequence of `num.samples' numbers by sequentially adding `1/num.samples',
and then uses `inverse.cdf.weights' to determine which components to 
sample. 

`branching.resample' deterministically copies `floor(weights)' number of components
and then draws another component i with probability equal to the residual for 
that component. Note: the actual number of components after resampling is random.
}

\value{Returns a vector of length `num.samples' with indices for sampled components.}

\references{Douc, R., Cappe, O.,  Moulines, E. (2005) Comparison of Resampling 
Schemes for Particle Filtering. _Proceedings of the 4th International Symposium 
on Image and Signal Processing and Analysis_

Carpenter, J., Clifford, P., Fearnhead, P. An improved particle filter for 
non-linear problems. _IEEE proceedings - Radar, Sonar and Navigation_ *146*, 2-7
}

\author{Jarad Niemi}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{resample}},\code{\link{renormalize}},
}

\examples{
ws = renormalize(runif(10))
multinomial.resample(ws)
residual.resample(ws,rrf="stratified")          
stratified.resample(ws,15)
systematic.resample(ws)
branching.resample(ws)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
