## ---- warning = FALSE, message = FALSE----------------------------------------
library(nimble)
library(basicMCMCplots)
library(coda)

## ---- warning = FALSE, message = FALSE----------------------------------------
library(nimbleSCR)

## -----------------------------------------------------------------------------
code <- nimbleCode({
  ## priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  p0 ~ dunif(0, 1)
  ## loop over individuals
  for(i in 1:n.individuals) {
    ## AC coordinates
    sxy[i,1] ~ dunif(0, x.max)
    sxy[i,2] ~ dunif(0, y.max)
    ## habitat constraint 
    ones[i] ~ dHabitatMask( s = sxy[i,1:2],
                            xmin = lowerCoords[1],
                            xmax = upperCoords[1],
                            ymin = lowerCoords[2],
                            ymax = upperCoords[2],
                            habitat = habitat.mx[1:y.max,1:x.max])
    ## latent dead/alive indicators
    z[i] ~ dbern(psi)
    ## likelihood
    y[i, 1:nMaxDetectors] ~ dbinomLocal_normal( detNums = nbDetections[i],
                                                   detIndices = yDets[i,1:nMaxDetectors],
                                                   size = trials[1:n.detectors],
                                                   p0 = p0,
                                                   s = sxy[i,1:2],
                                                   sigma = sigma,
                                                   trapCoords = detector.xy[1:n.detectors,1:2],
                                                   localTrapsIndices = detectorIndex[1:n.cells,1:maxNBDets],
                                                   localTrapsNum = nDetectors[1:n.cells],
                                                   resizeFactor = ResizeFactor,
                                                   habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                   indicator = z[i])
  }
  ## derived quantity: total population size
  N <- sum(z[1:n.individuals])
})

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

## ----eval = FALSE-------------------------------------------------------------
#  load("WolverineData.RData")

## -----------------------------------------------------------------------------
data <- list(y = my.jags.input$y,
             z = my.jags.input$z,
             detector.xy = my.jags.input$detector.xy,
             habitat.mx = my.jags.input$habitat.mx,
             ones = my.jags.input$OK,
             lowerCoords = c(0,0),
             upperCoords = c(
               dim(my.jags.input$habitat.mx)[2],
               dim(my.jags.input$habitat.mx)[1]),
             trials = rep(1, dim(my.jags.input$detector.xy)[1]))
constants <- list(n.individuals = my.jags.input$n.individuals,
                  n.detectors = dim(my.jags.input$detector.xy)[1],
                  y.max = dim(my.jags.input$habitat.mx)[1],
                  x.max = dim(my.jags.input$habitat.mx)[2])
inits <- list(sxy = inits.1$sxy,
              z = inits.1$z,
              p0 = 0.05,
              psi = 0.5,
              sigma = 6)

## ---- fig.width = 6, fig.height = 7-------------------------------------------
set.seed(2)

DetectorIndex <- getLocalObjects(habitatMask = data$habitat.mx,
                               coords = data$detector.xy,
                               dmax = 38,
                               resizeFactor = 24)

constants$y.maxDet <- dim(DetectorIndex$habitatGrid)[1]
constants$x.maxDet <- dim(DetectorIndex$habitatGrid)[2]
constants$ResizeFactor <- DetectorIndex$resizeFactor
constants$n.cells <- dim(DetectorIndex$localIndices)[1]
constants$maxNBDets <- DetectorIndex$numLocalIndicesMax
data$detectorIndex <- DetectorIndex$localIndices
data$nDetectors <- DetectorIndex$numLocalIndices
data$habitatIDDet <- DetectorIndex$habitatGrid

## -----------------------------------------------------------------------------
ySparse <- getSparseY(x = my.jags.input$y)
data$y <- ySparse$y[,,1]  
data$yDets <- ySparse$detIndices[,,1]
data$nbDetections <- ySparse$detNums[,1]
constants$nMaxDetectors <- ySparse$maxDetNums

## ----eval = FALSE-------------------------------------------------------------
#  Rmodel <- nimbleModel(code, constants, data, inits)

## ----message = FALSE, eval = FALSE--------------------------------------------
#  conf <- configureMCMC(Rmodel, monitors = c("N", "sigma", "p0"), print = FALSE)
#  conf$removeSamplers("sxy")
#  ACnodes <- paste0("sxy[", 1:constants$n.individuals, ", 1:2]")
#  for(node in ACnodes) {
#    conf$addSampler(target = node,
#                    type = "RW_block",
#                    control = list(adaptScaleOnly = TRUE),
#                    silent = TRUE)
#  }
#  Rmcmc <- buildMCMC(conf)

## ----eval = FALSE-------------------------------------------------------------
#  Cmodel <- compileNimble(Rmodel)
#  Cmcmc <- compileNimble(Rmcmc, project = Rmodel)
#  MCMC_runtime <- system.time(
#    samples <- runMCMC(Cmcmc, niter = 10000)
#  )

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  save(samples, MCMC_runtime, file = file.path(baseDir,"nimbleSCR/tests/testthat/wolverine_samples.RData"))

## ----echo = FALSE-------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "wolverine_samples.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/wolverine_samples.RData"))
}

## -----------------------------------------------------------------------------
round(MCMC_runtime[3] / 60, 1)

## -----------------------------------------------------------------------------
round(effectiveSize(samples),2) 

## -----------------------------------------------------------------------------
round(effectiveSize(samples)/MCMC_runtime[3],2)  

## -----------------------------------------------------------------------------
round(samplesSummary(samples), 2)

## -----------------------------------------------------------------------------
chainsPlot(samples)

