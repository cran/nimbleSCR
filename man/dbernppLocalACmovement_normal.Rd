% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dbernppLocalACmovement_normal.R
\name{dbernppLocalACmovement_normal}
\alias{dbernppLocalACmovement_normal}
\alias{rbernppLocalACmovement_normal}
\title{Local evaluation of a Bernoulli point process for activity center movement (normal kernel)}
\usage{
dbernppLocalACmovement_normal(
  x,
  lowerCoords,
  upperCoords,
  s,
  sd,
  baseIntensities,
  habitatGrid,
  habitatGridLocal,
  resizeFactor = 1,
  localHabWindowIndices,
  numLocalHabWindows,
  numGridRows,
  numGridCols,
  numWindows,
  log = 0
)

rbernppLocalACmovement_normal(
  n,
  lowerCoords,
  upperCoords,
  s,
  sd,
  baseIntensities,
  habitatGrid,
  habitatGridLocal,
  resizeFactor = 1,
  localHabWindowIndices,
  numLocalHabWindows,
  numGridRows,
  numGridCols,
  numWindows
)
}
\arguments{
\item{x}{Vector of x- and y-coordinates of a single spatial point (typically AC location at time t+1) scaled to the habitat (see (\code{\link{scaleCoordsToHabitatGrid}}).}

\item{lowerCoords, upperCoords}{Matrices of lower and upper x- and y-coordinates of all habitat windows. One row for each window.
Each window should be of size 1x1 (after rescaling if necessary).}

\item{s}{Vector of x- and y-coordinates of the isotropic bivariate normal distribution mean (i.e. the AC location).}

\item{sd}{Standard deviation of the isotropic bivariate normal distribution.}

\item{baseIntensities}{Vector of baseline habitat intensities for all habitat windows.}

\item{habitatGrid}{Matrix of habitat window indices. Cell values should correspond to the order of habitat windows in  \code{lowerCoords} and \code{upperCoords}. 
When the habitat grid only consists of a single row or column of windows, an additional row or column of dummy indices has to be added because the \code{nimble} model code requires a matrix.}

\item{habitatGridLocal}{Matrix of rescaled habitat grid cells indices, from localIndices returned by the  \code{getLocalObjects} function (}

\item{resizeFactor}{Aggregation factor used in the  \code{getLocalObjects} function to reduce the number of habitat grid cells.}

\item{localHabWindowIndices}{Matrix of indices of local habitat windows around each local habitat grid cell (\code{habitatGridLocal}), from localIndices returned by the \code{getLocalObjects} function.}

\item{numLocalHabWindows}{Vector of numbers of local habitat windows around all habitat grid cells, as returned by the getLocalObjects function (object named \code{numLocalIndices}). 
The ith number gives the number of local (original) habitat windows for the ith (rescaled) habitat window.}

\item{numGridRows, numGridCols}{Numbers of rows and columns of the \code{habitatGrid}.}

\item{numWindows}{Number of habitat windows. This value (positive integer) is used to truncate \code{lowerCoords} and \code{upperCoords} 
so that extra rows beyond \code{numWindows} are ignored.}

\item{log}{Logical argument, specifying whether to return the log-probability of the distribution.}

\item{n}{Integer specifying the number of realisations to generate.  Only n = 1 is supported.}
}
\value{
The (log) probability density of the observation vector \code{x}.
}
\description{
Density and random generation functions of the Bernoulli point process for activity center movement between occasions based on a bivariate normal distribution and local evaluation.
}
\details{
The \code{dbernppLocalACmovement_normal} distribution is a NIMBLE custom distribution which can be used to model and simulate
movement of activity centers between consecutive occasions in open population models.
The distribution assumes that the new individual activity center location (\emph{x})
follows an isotropic multivariate normal centered on the previous activity center (\emph{s}) with standard deviation (\emph{sd}).
The local evaluation technique is implemented.
}
\examples{

# Creat habitat grid
habitatGrid <- matrix(c(1:(4^2)), nrow = 4, ncol=4, byrow = TRUE)
coordsHabitatGridCenter <- matrix(c(0.5, 3.5,
                                    1.5, 3.5,
                                    2.5, 3.5,
                                    3.5, 3.5,
                                    0.5, 2.5,
                                    1.5, 2.5,
                                    2.5, 2.5,
                                    3.5, 2.5,
                                    0.5, 1.5,
                                    1.5, 1.5,
                                    2.5, 1.5,
                                    3.5, 1.5,
                                    0.5, 0.5,
                                    1.5, 0.5,
                                    2.5, 0.5,
                                    3.5, 0.5), ncol = 2,byrow = TRUE)
colnames(coordsHabitatGridCenter) <- c("x","y")
# Create habitat windows
lowerCoords <- coordsHabitatGridCenter-0.5
upperCoords <- coordsHabitatGridCenter+0.5
colnames(lowerCoords) <- colnames(upperCoords) <- c("x","y")
# Plot check
plot(lowerCoords[,"y"]~lowerCoords[,"x"],pch=16, xlim=c(0,4), ylim=c(0,4),col="red") 
points(upperCoords[,"y"]~upperCoords[,"x"],col="red",pch=16) 
points(coordsHabitatGridCenter[,"y"]~coordsHabitatGridCenter[,"x"],pch=16) 

# Rescale coordinates 
ScaledLowerCoords <- scaleCoordsToHabitatGrid(coordsData =  lowerCoords,
                                              coordsHabitatGridCenter = coordsHabitatGridCenter)
ScaledUpperCoords <- scaleCoordsToHabitatGrid(coordsData =  upperCoords,
                                              coordsHabitatGridCenter = coordsHabitatGridCenter)
ScaledUpperCoords$coordsDataScaled[,2] <- ScaledUpperCoords$coordsDataScaled[,2] + 1
ScaledLowerCoords$coordsDataScaled[,2] <- ScaledLowerCoords$coordsDataScaled[,2] - 1
habitatMask <- matrix(1, nrow = 4, ncol=4, byrow = TRUE)
# Create local objects 
HabWindowsLocal <- getLocalObjects(habitatMask = habitatMask,
                                   coords = coordsHabitatGridCenter,
                                   dmax=4,
                                   resizeFactor = 1,
                                   plot.check = TRUE
)

s <- c(1, 1) # Currrent activity center location
sd <- 0.1
numWindows <- nrow(coordsHabitatGridCenter)
baseIntensities <- rep(1,numWindows)
numRows <- nrow(habitatGrid)
numCols <- ncol(habitatGrid)

# The log probability density of moving from (1,1) to (1.2, 0.8) 
dbernppLocalACmovement_normal(x = c(1.2, 0.8), lowerCoords, upperCoords, s,
                              sd, baseIntensities, habitatGrid, 
                              HabWindowsLocal$habitatGrid, HabWindowsLocal$resizeFactor,
                              HabWindowsLocal$localIndices, HabWindowsLocal$numLocalIndices,
                              numRows, numCols, numWindows, log = TRUE)


}
\references{
W. Zhang, J. D. Chipperfield, J. B. Illian, P. Dupont, C. Milleret, P. de Valpine and R. Bischof. 2020. 
A hierarchical point process model for spatial capture-recapture data. bioRxiv. DOI 10.1101/2020.10.06.325035 

C. Milleret, P. Dupont, C. Bonenfant, H. Brøseth, Ø. Flagstad, C. Sutherland and R. Bischof. 2019. 
A local evaluation of the individual state-space to scale up Bayesian spatial capture-recapture. Ecology and Evolution 9:352-363
}
\author{
Wei Zhang and Cyril Milleret
}
