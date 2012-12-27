\name{resample}
\alias{resample}

\title{Resample}

\description{A set of resampling functions with unbiased number of replicates.}

\usage{
resample(weights, num.samples=length(weights),
         method = c("stratified","residual","multinomial","systematic","branching"),
         nonuniformity = c("none","ess","cov","entropy"), threshold = 0.5,
         rrf = "stratified", engine="R", log=TRUE, normalized=FALSE)
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
  \item{threshold}{a scalar in [0,1] indicating when to resample. If `nonuniformity="none"',
  resampling is always performed regardless of threshold. If `nonuniformity="ess"'
  or "entropy", resampling is performed when `nonuniformity/num.samples<threshold'. If
  `nonuniformity="cov"', resampling is performed when `nonuniformity/log2(num.samples)>
  threshold'.}
  \item{rrf}{for residual resampling, the resampling function to use on the residual}
  \item{engine}{run using "R" or "C" code}      
  \item{log}{if TRUE, the weights are assumed to be logged}
  \item{normalized}{if TRUE, the weights are assumed to be normalized and no error checking is done}
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
ws = dnorm(0:5, log=TRUE)
resample(ws)
resample(ws,13)
resample(ws,method="residual")
resample(ws,method="residual",rrf="stratified")
resample(ws,17,"stratified","ess",0.5)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ htest }

