## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 8, fig.height = 10) 

## ----echo = FALSE-------------------------------------------------------------
##load("C:/Users/cymi/Downloads/WolverineData.RData")
#load("~/Downloads/WolverineData.RData")
## So it works for everyone
if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else baseDir <- NULL

if(is.null(baseDir)) {
    load(system.file("extdata", "WolverineData.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/WolverineData.RData"))
}

## ---- warning = FALSE, message = FALSE----------------------------------------
## Load packages
library(nimble)
library(nimbleSCR)
library(basicMCMCplots)

## ---- warning = FALSE, message = FALSE----------------------------------------
## Create habitat grid
coordsHabitatGridCenter <- cbind(rep(seq(75, 5, by = -10), 10),
                                 sort(rep(seq(5, 100, by = 10), 8)))
colnames(coordsHabitatGridCenter) <- c("x","y")

## Create trap grid
coordsObsCenter <- cbind(rep(seq(15, 65, by = 10), 8),
                         sort(rep(seq(15, 85, by = 10), 6)))
colnames(coordsObsCenter) <- c("x","y")

## Plot check
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"],
     xlim = c(0,80), ylim = c(0,100),
     pch = 1, cex = 1.5) 
points(coordsObsCenter[,"y"] ~ coordsObsCenter[,"x"], col="red", pch=16 ) 
par(xpd=TRUE)
legend(x = 7, y = 13,
       legend=c("Habitat window centers", "Observation window centers"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1,16),
       col=c("black", "red"),
       bty = 'n')


## ---- warning = FALSE, message = FALSE----------------------------------------
## Rescale coordinates
scaledObjects <- scaleCoordsToHabitatGrid(
  coordsData = coordsObsCenter,            
  coordsHabitatGridCenter = coordsHabitatGridCenter)

## Get lower and upper cell coordinates
lowerAndUpperCoords <- getWindowCoords(
  scaledHabGridCenter = scaledObjects$coordsHabitatGridCenterScaled,
  scaledObsGridCenter = scaledObjects$coordsDataScaled,
  plot.check = F)


## ---- warning = FALSE, message = FALSE----------------------------------------
modelCode <- nimbleCode({
  ##---- SPATIAL PROCESS 
  ## Prior for AC distribution parameter
  habCoeffSlope ~ dnorm(0, sd = 10)
  
  ## Intensity of the AC distribution point process
  habIntensity[1:numHabWindows] <- exp(habCoeffSlope * habCovs[1:numHabWindows])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  
  ## AC distribution
  for(i in 1:M){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
      numGridRows =  numGridRows,
      numGridCols = numGridCols
    )
  }
  
  ##---- DEMOGRAPHIC PROCESS
  ## Prior for data augmentation
  psi ~ dunif(0,1)
  
  ## Data augmentation
  for (i in 1:M){
    z[i] ~ dbern(psi)
  }
  
  ##---- DETECTION PROCESS
  ## Priors for detection parameters
  sigma ~ dunif(0, 50)
  detCoeffInt ~ dnorm(0, sd = 10)
  detCoeffSlope ~ dnorm(0, sd = 10)
  
  ## Intensity of the detection point process
  detIntensity[1:numObsWindows] <- exp(detCoeffInt + detCoeffSlope * detCovs[1:numObsWindows]) 
  
  ## Detection process
  for (i in 1:M){
    y[i, 1:numMaxPoints, 1:3] ~ dpoisppDetection_normal(
      lowerCoords = obsLoCoords[1:numObsWindows, 1:2],
      upperCoords = obsUpCoords[1:numObsWindows, 1:2],
      s = sxy[i, 1:2],
      sd = sigma,
      baseIntensities = detIntensity[1:numObsWindows],
      numMaxPoints = numMaxPoints,
      numWindows = numObsWindows,
      indicator = z[i]
    )
  }
  
  ##---- DERIVED QUANTITIES
  ## Number of individuals in the population
  N <- sum(z[1:M])
})

## ---- warning = FALSE, message = FALSE----------------------------------------
sigma <- 1
psi <- 0.6
detCoeffInt <- 0.1
detCoeffSlope <- 0.5
habCoeffSlope <- -1.5

## ---- warning = FALSE, message = FALSE----------------------------------------
M <- 150

## ---- warning = FALSE, message = FALSE----------------------------------------
numMaxPoints <- 19 + 1

## ---- warning = FALSE, message = FALSE----------------------------------------
detCovs <- runif(dim(lowerAndUpperCoords$lowerObsCoords)[1],-1,1)
habCovs <- runif(dim(lowerAndUpperCoords$lowerHabCoords)[1],-1,1)

## ---- warning = FALSE, message = FALSE----------------------------------------
nimConstants <- list( M = M,
                      numObsWindows = dim(lowerAndUpperCoords$lowerObsCoords)[1],
                      numMaxPoints = numMaxPoints,
                      numHabWindows = dim(lowerAndUpperCoords$upperHabCoords)[1],
                      habitatGrid = lowerAndUpperCoords$habitatGrid,
                      numGridRows = dim(lowerAndUpperCoords$habitatGrid)[1],
                      numGridCols = dim(lowerAndUpperCoords$habitatGrid)[2])

nimData <- list( obsLoCoords = lowerAndUpperCoords$lowerObsCoords,
                 obsUpCoords = lowerAndUpperCoords$upperObsCoords,
                 lowerHabCoords = lowerAndUpperCoords$lowerHabCoords,
                 upperHabCoords = lowerAndUpperCoords$upperHabCoords,
                 detCovs = detCovs,
                 habCovs = habCovs)

## ---- warning = FALSE, message = FALSE----------------------------------------
nimInits <- list( psi = psi,
                  sigma = sigma,
                  detCoeffInt = detCoeffInt,
                  detCoeffSlope = detCoeffSlope,
                  habCoeffSlope = habCoeffSlope)

## ---- warning = FALSE, message = FALSE----------------------------------------
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,
                      calculate = F)

## ---- warning = FALSE, message = FALSE----------------------------------------
nodesToSim <- model$getDependencies(names(nimInits), self = F)
set.seed(1)
model$simulate(nodesToSim, includeData = FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
N <- sum(model$z)

## ---- warning = FALSE, message = FALSE----------------------------------------
i = 7
## Number of detections for individual i
model$y[i,1,1]

## Plot of the habitat and trap grids 
plot( scaledObjects$coordsHabitatGridCenterScaled[,"y"] ~ scaledObjects$coordsHabitatGridCenterScaled[,"x"],
      pch = 1, cex = 0.5)
rect( xleft = lowerAndUpperCoords$lowerHabCoords[,1] ,
     ybottom = lowerAndUpperCoords$lowerHabCoords[,2] ,
     xright = lowerAndUpperCoords$upperHabCoords[,1],
     ytop = lowerAndUpperCoords$upperHabCoords[,2],
     col = adjustcolor("red", alpha.f = 0.4),
     border = "red")

rect( xleft = lowerAndUpperCoords$lowerObsCoords[,1] ,
     ybottom = lowerAndUpperCoords$lowerObsCoords[,2] ,
     xright = lowerAndUpperCoords$upperObsCoords[,1],
     ytop = lowerAndUpperCoords$upperObsCoords[,2],
     col = adjustcolor("blue",alpha.f = 0.4),
     border = "blue")

## Plot the activity center of individual i
points( model$sxy[i, 2] ~ model$sxy[i, 1],
        col = "orange", pch = 16)

## Plot detections of individual i
dets <- model$y[i,2:model$y[i,1,1], ]
points( dets[,2] ~ dets[,1],
        col = "green", pch = 16)

par(xpd = TRUE)
legend(x = -1, y = 13,
       legend = c("Habitat windows",
                  "Observation windows",
                  "Simulated AC",
                  "Detections"),
       pt.cex = c(1,1),
       horiz = T,
       pch = c(16, 16, 16, 16),
       col = c("red", "blue", "orange", "green"),
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
nimData1 <- nimData
nimData1$y <- model$y
nimInits1 <- nimInits
nimInits1$z <- model$z
nimInits1$sxy <- model$sxy

## Create and compile the NIMBLE model
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData1,
                      inits = nimInits1,
                      check = F,
                      calculate = F)
cmodel <- compileNimble(model)
## Check the initial log-likelihood 
cmodel$calculate()

## ---- warning = FALSE, message = FALSE----------------------------------------
MCMCconf <- configureMCMC(model = model,
                          monitors  = c("N","sigma","psi","detCoeffInt",
                                        "detCoeffSlope","habCoeffSlope"),
                          control = list(reflective = TRUE),
                          thin = 1)
MCMC <- buildMCMC(MCMCconf)

## ----eval = F-----------------------------------------------------------------
#  
#  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
#  
#  ## Run MCMC
#  MCMCRuntime <- system.time(samples <- runMCMC( mcmc = cMCMC,
#                                                        nburnin = 500,
#                                                        niter = 10000,
#                                                        nchains = 3,
#                                                        samplesAsCodaMCMC = TRUE))
#  

## ----eval = F, echo = FALSE---------------------------------------------------
#  save(samples, MCMCRuntime, file = file.path(baseDir,"nimbleSCR/inst/extdata/pointProcess_samples.RData"))
#  #file.path("C:/Personal_Cloud/OneDrive/Work/nimbleSCR/nimbleSCR/inst/extdata","pointProcess_samples.RData"))

## ----echo = FALSE-------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "pointProcess_samples.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/pointProcess_samples.RData"))
}

## ---- warning = FALSE, message = FALSE----------------------------------------

## Print runtime
MCMCRuntime

## Traceplots and density plots for the tracked parameters
chainsPlot(samples)

## ---- warning = FALSE, message = FALSE----------------------------------------
modelCodeSemiCompleteLikelihood <- nimbleCode({
  #----- SPATIAL PROCESS
  ## Priors
  habCoeffInt ~ dnorm(0, sd = 10)
  habCoeffSlope ~ dnorm(0, sd = 10)

  ## Intensity of the AC distribution point process
  habIntensity[1:numHabWindows] <- exp(habCoeffInt + habCoeffSlope * habCovs[1:numHabWindows])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sum(habIntensity[1:numHabWindows] ))

  ## AC distribution
  for(i in 1:nDetected){
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
      numGridRows =  numGridRows,
      numGridCols = numGridCols
    )
  }

  ##---- DEMOGRAPHIC PROCESS
  ## Number of individuals in the population
  N ~ dpois(sumHabIntensity)
  ## Number of detected individuals
  nDetectedIndiv ~ dbin(probDetection, N)

  ##---- DETECTION PROCESS
  ## Probability that an individual in the population is detected at least once
  ## i.e. 1 - void probability over all detection windows
  probDetection <- 1 - marginalVoidProbNumIntegration(
    quadNodes = quadNodes[1:nNodes, 1:2, 1:numHabWindows],
    quadWeights = quadWeights[1:numHabWindows],
    numNodes = numNodes[1:numHabWindows],
    lowerCoords = obsLoCoords[1:numObsWindows, 1:2],
    upperCoords = obsUpCoords[1:numObsWindows, 1:2],
    sd = sigma,
    baseIntensities = detIntensity[1:numObsWindows],
    habIntensities = habIntensity[1:numHabWindows],
    sumHabIntensity = sumHabIntensity,
    numObsWindows = numObsWindows,
    numHabWindows = numHabWindows
  )

  ## Priors for detection parameters
  sigma ~ dunif(0, 50)
  detCoeffInt ~ dnorm(0, sd = 10)
  detCoeffSlope ~ dnorm(0, sd = 10)

  ## Intensity of the detection point process
  detIntensity[1:numObsWindows] <- exp(detCoeffInt + detCoeffSlope * detCovs[1:numObsWindows])
  ## Detection process
  ## Note that this conditions on the fact that individuals are detected (at least once)
  ## So, at the bottom of this model code we deduct log(probDetection) from the log-likelihood
  ## function for each individual
  for (i in 1:nDetected){
    y[i, 1:numMaxPoints, 1:3] ~ dpoisppDetection_normal(
      lowerCoords = obsLoCoords[1:numObsWindows, 1:2],
      upperCoords = obsUpCoords[1:numObsWindows, 1:2],
      s = sxy[i, 1:2],
      sd = sigma,
      baseIntensities = detIntensity[1:numObsWindows],
      numMaxPoints = numMaxPoints,
      numWindows = numObsWindows,
      indicator = 1
    )
  }
  ## Normalization: normData can be any scalar in the data provided when building the model
  ## The dnormalizer is a custom distribution defined for efficiency, where the input data
  ## does not matter. It makes it possible to use the general dpoippDetection_normal function
  ## when either data augmentation or the SCDL is employed
  logDetProb <- log(probDetection)
  normData ~ dnormalizer(logNormConstant = -M * logDetProb)
})

## ---- warning = FALSE, message = FALSE----------------------------------------
idDetected <- which(nimData1$y[,1,1] > 0)
## Subset data to detected individuals only
nimData1$y <- nimData1$y[idDetected,,]
## Provide the number of detected individuals as constant
nimConstants$nDetected <- length(idDetected)
## With this model, we also need to provide the number of detected individuals as data for the estimation of population size.
nimData1$nDetectedIndiv <- length(idDetected)
## As mentioned above, "normData" can take any value.
nimData1$normData <- 1

## We also provide initial values for the new parameters that need to be estimated
nimInits1$N <- 100
nimInits1$habCoeffInt <- 0.5
nimInits1$sxy <- nimInits1$sxy[idDetected,]

## ---- warning = FALSE, message = FALSE----------------------------------------
## Number of equal subintervals for each dimension of a grid cell
nPtsPerDim <- 2
## Number of points to use for the numerical integration for each grid cell
nNodes <- nPtsPerDim^2
## Generate midpoint nodes coordinates for numerical integration using the "getMidPointNodes" function
nodesRes <- getMidPointNodes( nimData1$lowerHabCoords,
                           nimData1$upperHabCoords,
                           nPtsPerDim)
## Add this info to the data and constant objects
nimData1$quadNodes <- nodesRes$quadNodes
nimData1$quadWeights <- nodesRes$quadWeights
nimData1$numNodes <- rep(nNodes,dim(nimData1$lowerHabCoords)[1])
nimConstants$nNodes <- dim(nodesRes$quadNodes)[1]

## ---- eval = T, warning = FALSE, message = FALSE------------------------------
model <- nimbleModel(code = modelCodeSemiCompleteLikelihood,
                     constants = nimConstants,
                     data = nimData1,
                     inits = nimInits1,
                     check = F,
                     calculate = F)

model$calculate()
cmodel <- compileNimble(model)
cmodel$calculate()


MCMCconf <- configureMCMC(model = model,
                          monitors =  c("N","sigma","probDetection","habCoeffInt", "detCoeffInt","detCoeffSlope","habCoeffSlope"),
                          control = list(reflective = TRUE),
                          thin = 1)

MCMC <- buildMCMC(MCMCconf)

## ---- eval = F, warning = FALSE, message = FALSE------------------------------
#  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
#  
#  ## Run MCMC
#  MCMCRuntime1 <- system.time(samples1 <- runMCMC( mcmc = cMCMC,
#                                                       nburnin = 500,
#                                                       niter = 10000,
#                                                       nchains = 3,
#                                                       samplesAsCodaMCMC = TRUE))
#  

## ----eval = F, echo = FALSE---------------------------------------------------
#  save(samples1, MCMCRuntime1, file = file.path(baseDir,"nimbleSCR/inst/extdata/pointProcess_samples1.RData"))

## ----echo = F-----------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "pointProcess_samples1.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/pointProcess_samples1.RData"))
}

## ---- warning = FALSE, message = FALSE----------------------------------------

MCMCRuntime1

## Plot check
chainsPlot(samples1)

