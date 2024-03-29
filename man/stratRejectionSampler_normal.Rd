% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratRejectionSampler_normal.R
\name{stratRejectionSampler_normal}
\alias{stratRejectionSampler_normal}
\title{Stratified rejection sampler for multivariate normal point process}
\usage{
stratRejectionSampler_normal(
  numPoints,
  lowerCoords,
  upperCoords,
  s,
  windowIntensities,
  sd
)
}
\arguments{
\item{numPoints}{Number of spatial points to generate.}

\item{lowerCoords, upperCoords}{Matrices of lower and upper x- and y-coordinates of a set of detection windows. One row for each window.}

\item{s}{Vector of x- and y-coordinates of of the isotropic multivariate normal distribution mean.}

\item{windowIntensities}{Vector of integrated intensities over all detection windows.}

\item{sd}{Standard deviation of the isotropic multivariate normal distribution.}
}
\value{
A matrix of x- and y-coordinates of the generated points. One row corresponds to one point.
}
\description{
Simulate data using a stratified rejection sampler from a point process with an isotropic multivariate normal decay kernel.
}
\examples{
numPoints <- 10
lowerObsCoords <- matrix(c(0, 0, 1, 0, 0, 1, 1, 1), nrow = 4, byrow = TRUE)
upperObsCoords <- matrix(c(1, 1, 2, 1, 1, 2, 2, 2), nrow = 4, byrow = TRUE)
s <- c(1, 1)
windowIntensities <- c(1:4)
sd <- 0.1
set.seed(0)
stratRejectionSampler_normal(numPoints, lowerObsCoords, upperObsCoords, s, windowIntensities, sd)

}
\author{
Joseph D. Chipperfield and Wei Zhang
}
