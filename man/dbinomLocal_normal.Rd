% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbinomLocal_normal.R
\name{dbinomLocal_normal}
\alias{dbinomLocal_normal}
\alias{rbinomLocal_normal}
\title{Local evaluation of a binomial SCR observation process}
\usage{
dbinomLocal_normal(
  x,
  detNums,
  detIndices,
  size,
  p0,
  sigma,
  s,
  trapCoords,
  localTrapsIndices,
  localTrapsNum,
  resizeFactor = 1,
  habitatGrid,
  indicator = 1,
  log = 0
)

rbinomLocal_normal(
  n = 1,
  detNums,
  detIndices,
  size,
  p0,
  sigma,
  s,
  trapCoords,
  localTrapsIndices,
  localTrapsNum,
  resizeFactor = 1,
  habitatGrid,
  indicator = 1
)
}
\arguments{
\item{x}{Vector of individual detection frequencies, as returned by the \code{getSparseY} function 
(padded with -1's to maintain the square structure of the observation data).}

\item{detNums}{Number of detections recorded in x, as returned by the \code{getSparseY} function.}

\item{detIndices}{Vector of the detector indices where the detections in x were recorded, as returned by the \code{getSparseY} function.}

\item{size}{Vector of the number of trials (zero or more) for each trap (trapCoords).}

\item{p0}{Baseline detection probability used in the half-normal detection function.}

\item{sigma}{Scale parameter of the half-normal detection function.}

\item{s}{Individual activity center x- and y-coordinates.}

\item{trapCoords}{Matrix of x- and y-coordinates of all traps.}

\item{localTrapsIndices}{Matrix of indices of local traps around each habitat grid cell, as returned by the \code{getLocalObjects} function.}

\item{localTrapsNum}{Vector of numbers of local traps around all habitat grid cells, as returned by the getLocalObjects function.}

\item{resizeFactor}{Aggregation factor used in the \code{getLocalObjects} function to reduce the number of habitat grid cells to retrieve local traps for.}

\item{habitatGrid}{Matrix of habitat grid cells indices.}

\item{indicator}{Logical argument, specifying whether the individual is available for detection.}

\item{log}{Logical argument, specifying whether to return the log-probability of the distribution.}

\item{n}{Integer specifying the number of realisations to generate.  Only n = 1 is supported.}
}
\value{
The log-likelihood value associated with the vector of detections, given the location of the activity center (s),
 and the half-normal detection function : \eqn{p = p0 * exp(-d^2 / \sigma^2)}.
}
\description{
The \code{dbinomLocal_normal} distribution is a NIMBLE custom distribution which can be used to model 
the binomial observations (x) of a single individual over a set of detectors defined by their 
coordinates (trapCoords). The distribution assumes that the detection probability at any detector 
follows a half-normal function of the distance between the individual's activity center (s) and the detector location.
}
\details{
The \code{dbinomLocal_normal} distribution incorporates three features to increase computation efficiency:
\enumerate{
\item A local evaluation of the detection probability calculation (see Milleret et al. (2019) <doi:10.1002/ece3.4751> for more details).
\item It uses a sparse matrix representation (x, detIndices, detNums) of the observation data to reduce the size of objects to be processed.
\item It uses an indicator (indicator) to shortcut calculations for individuals unavailable for detection.
}

The \code{dbinomLocal_normal} distribution requires that the x- and y- coordinates should be scaled to the habitat (\code{\link{scaleCoordsToHabitatGrid}})
}
\examples{
## define model code
code <- nimbleCode({
    psi ~ dunif(0,1)
    p0 ~ dunif(0,1)
    sigma ~ dunif(0,100)
    N <- sum(z[1:M])
    for(i in 1:M) {
        s[i, 1] ~ dunif(0, 100)
        s[i, 2] ~ dunif(0, 100)
        z[i] ~ dbern(psi)
        y[i,1:maxDetNum] ~ dbinomLocal_normal(detNums,
                                                 detIndices,
                                                 size,
                                                 p0,
                                                 sigma,
                                                 s[i,1:2],
                                                 trapCoords,
                                                 localTrapsIndices,
                                                 localTrapsNum,
                                                 resizeFactor,
                                                 habitatGrid,
                                                 z[i])
    }
})
 
## create NIMBLE model object
## Rmodel <- nimbleModel(code, ...)
 
## use model object for MCMC, etc.

}
\author{
Cyril Milleret
}
