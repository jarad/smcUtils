\name{smcUtils-package}
\alias{smcUtils-package}
\alias{smcUtils}
\docType{package}
\title{Provides utility functions for sequential Monte Carlo.}
\description{Contains functions for resampling, evaluating non-uniformity, and 
normalizing particle weights in sequential Monte Carlo, a.k.a. particle 
filtering.}
\details{
\tabular{ll}{
Package: \tab smcUtils\cr
Type: \tab Package\cr
Version: \tab 0.1.2\cr
Date: \tab 2011-11-06\cr
License: \tab GPL-3\cr
LazyLoad: \tab yes\cr
}
The workhorse of the package is `resample' which will perform multinomial, 
residual, stratified, systematic, or branching resampling. It will perform 
the resampling when a non-uniformity heuristic crosses a threshold. 
}

\author{
Author: Jarad Niemi

Maintainer: Jarad Niemi <niemi@iastate.edu>

}
\references{
Douc, R., Cappe, O.,  Moulines, E. (2005) Comparison of Resampling 
Schemes for Particle Filtering. _Proceedings of the 4th International Symposium 
on Image and Signal Processing and Analysis_

Carpenter, J., Clifford, P., Fearnhead, P. An improved particle filter for 
non-linear problems. _IEE proceedings - Radar, Sonar and Navigation_ *146*, 2-7
}
% Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
\keyword{ internal }
