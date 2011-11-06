\name{resample}
\alias{resample}

\title{Resample}

\description{A set of resampling functions with unbiased number of replicates.}

\usage{
resample(weights, num.samples=length(weights),
         method = c("multinomial","residual","stratified","systematic","branching"),
         nonuniformity = c("none","ess","cov","entropy"), threshold = NULL,
         resample.function = multinomial.resample)
}

\arguments{
  \item{weights}{a vector of normalized weights}
  \item{num.samples}{a scalar indicating the number of samples to return (for `branching.resample',
  `num.samples' is the expected number of samples as the actual number is random)}
  \item{method}{a character string indicating the resampling method to use. One
  of `"multinomial"' (default), `"residual"', `"stratified"', `"systematic"', or
  `"branching"', can be abbreviated.}
  \item{nonuniformity}{a character string indicating which hueristic to use for
  measuring weight nonuniformity. One of `"none"' (default), `"ess"', `"cov"', or
  `"entropy"', can be abbreviated. If "none", resampling will always be performed.}
  \item{threshold}{a scalar indicating when to resample. If `nonuniformity="none"',
  resampling is always performed regardless of threshold. If `nonuniformity="ess"'
  or "entropy", resampling is performed when `nonuniformity<threshold'. If
  `nonuniformity="cov"', resampling is performed when `nonuniformity>
  threshold'. The default threshold is `0.5*num.samples' for "ess" and "cov", and
  `0.5*log2(num.samples)' for "entropy".}
  \item{resample.function}{the resampling function to use on the remainder
  (only used when method = "residual")}
}

\value{
  \item{weights}{the component weights}
  \item{indices}{an integer vector containing the indices of resampled components.
  If no resampling was performed, then `indices=1:length(weights)'.}
}

\author{Jarad Niemi}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{resampling}}}

\examples{
ws = renormalize.weights(runif(10))
resample(ws)
resample(ws,13)
resample(ws,method="residual")
resample(ws,method="residual",resample.function="stratified")
resample(ws,17,"stratified","ess",5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }
