## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 8, fig.height = 10) 

## ----echo = FALSE-------------------------------------------------------------
## So it works for everyone
if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else if(Sys.info()['user'] == 'admin') {                  ## Soumen
  baseDir <- '~/GitHubSD/nimbleSCR/' 
} else baseDir <- NULL

## ---- warning = FALSE, message = FALSE----------------------------------------
# LOAD PACKAGES
library(nimble)
library(nimbleSCR)
library(basicMCMCplots)

## ---- warning = FALSE, message = FALSE----------------------------------------
# CREATE HABITAT GRID 
coordsHabitatGridCenter <- cbind(rep(seq(11.5, 0.5, by=-1), 12),
                                 sort(rep(seq(0.5, 11.5, by=1), 12)))
colnames(coordsHabitatGridCenter) <- c("x","y")

# CREATE TRAP GRID
trapCoords <- cbind(rep(seq(2.5, 9.5,by=1),8),
               sort(rep(seq(2.5, 9.5,by=1),8)))
colnames(trapCoords) <- c("x","y")

# PLOT CHECK
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch = 1, cex = 1.5) #pch=16) 
points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16 ) 
par(xpd=TRUE)
legend(x = 7, y = 13,
       legend=c("Habitat grid centers", "Traps"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1,16),
       col=c("black", "red"),
       bty = 'n')


## ---- warning = FALSE, message = FALSE----------------------------------------
ScaledtrapCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
                                             coordsHabitatGridCenter = coordsHabitatGridCenter)$coordsDataScaled
habitatMask <- matrix(1, nrow = 12, ncol= 12, byrow = TRUE)

## ---- warning = FALSE, message = FALSE----------------------------------------
trapLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = ScaledtrapCoords,
                             dmax = 7,
                             plot.check = FALSE
)

## ---- warning = FALSE, message = FALSE----------------------------------------
modelCode <- nimbleCode({
  ## priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  p0 ~ dunif(0, 1)
  ## loop over individuals
  for(i in 1:M) {
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
    y[i, 1:lengthYCombined] ~ dbinomLocal_normal(size = trials[1:n.traps],
                                                 p0 = p0,
                                                 s = sxy[i,1:2],
                                                 sigma = sigma,
                                                 trapCoords = trapCoords[1:n.traps,1:2],
                                                 localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                                 localTrapsNum = nTraps[1:n.cells],
                                                 habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                 indicator = z[i],
                                                 lengthYCombined = lengthYCombined)
   }
  ## derived quantity: total population size
  N <- sum(z[1:M])
})

## ---- warning = FALSE, message = FALSE----------------------------------------
# PARAMETERS
p0 <- 0.2
sigma <- 2
psi <- 0.5

## ---- warning = FALSE, message = FALSE----------------------------------------
M <- 200

## ---- warning = FALSE, message = FALSE----------------------------------------
lengthYCombined <- 1 + trapLocal$numLocalIndicesMax*2

## ---- warning = FALSE, message = FALSE----------------------------------------
nimConstants <- list(M = M,
                     n.traps = dim(ScaledtrapCoords)[1],
                     y.max = dim(habitatMask)[1],
                     x.max = dim(habitatMask)[2],
                     y.maxDet = dim(trapLocal$habitatGrid)[1],
                     x.maxDet = dim(trapLocal$habitatGrid)[2],
                     n.cells = dim(trapLocal$localIndices)[1],
                     maxNBDets = trapLocal$numLocalIndicesMax,
                     trapIndex = trapLocal$localIndices,
                     nTraps = trapLocal$numLocalIndices,
                     habitatIDDet = trapLocal$habitatGrid,
                     lengthYCombined = lengthYCombined)


nimData <- list(trapCoords = ScaledtrapCoords,
                habitat.mx = habitatMask,
                ones = rep(1, nimConstants$M),
                lowerCoords = c(min(coordsHabitatGridCenter[,1]) - 0.5, min(coordsHabitatGridCenter[,2]) - 0.5),
                upperCoords = c(max(coordsHabitatGridCenter[,1]) + 0.5, max(coordsHabitatGridCenter[,2]) + 0.5),
                trials = rep(1, dim(ScaledtrapCoords)[1])
)
# We set the parameter values as inits
nimInits <- list(p0 = p0,
                 psi = psi,
                 sigma = sigma)

## ---- warning = FALSE, message = FALSE----------------------------------------
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = F)  

## ---- warning = FALSE, message = FALSE----------------------------------------
# FIRST WE GET THE NODES TO SIMULATE
nodesToSim <- model$getDependencies(c("sxy", "z"), self=T)
# THEN WE SIMULATE THOSE NODES 
set.seed(100)
model$simulate(nodesToSim, includeData = FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
N <- sum(model$z)

## ---- warning = FALSE, message = FALSE----------------------------------------
i <- 5
#NUMBER OF DETECTIONS SIMULATED FOR INDIVIDUAL i
model$y[i,1]

#LET'S PLOT THEM
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch = 1, cex = 1.5) #pch=16) 
points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16 )  
# PLOT ITS ACTIVITY CETENR 
points(model$sxy[i, 2] ~ model$sxy[i, 1],col="orange", pch=16) 

whichdets <- model$y[i, (trapLocal$numLocalIndicesMax + 2):(trapLocal$numLocalIndicesMax+model$y[i, 1] + 1)]
points(ScaledtrapCoords[whichdets, "y"] ~
         ScaledtrapCoords[whichdets, "x"], col="blue", pch=16) 

par(xpd=TRUE)
legend(x = -1, y = 13,
       legend=c("Habitat grid centers", "Traps", "Simulated AC", "Detections"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1, 16, 16, 16),
       col=c("black", "red", "orange", "blue"),
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
nimData1 <- nimData
nimData1$y <- model$y
nimInits1 <- nimInits
nimInits1$z <- model$z
nimInits1$sxy <- model$sxy

# CREATE AND COMPILE THE NIMBLE MODEL
model <- nimbleModel( code = modelCode,
                      constants = nimConstants,
                      data = nimData1,
                      inits = nimInits1,
                      check = F,
                      calculate = F)

model$calculate()

## ----eval = FALSE-------------------------------------------------------------
#  cmodel <- compileNimble(model)
#  cmodel$calculate()

## -----------------------------------------------------------------------------
MCMCconf <- configureMCMC(model = model,
                          monitors = c("N", "sigma", "p0","psi"),
                          control = list(reflective = TRUE),
                          thin = 1)

## Add block sampling of sxy coordinates see Turek et al 2021 for further details 
MCMCconf$removeSamplers("sxy")
ACnodes <- paste0("sxy[", 1:nimConstants$M, ", 1:2]")
for(node in ACnodes) {
  MCMCconf$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

MCMC <- buildMCMC(MCMCconf)

## ----eval = FALSE-------------------------------------------------------------
#  cMCMC <- compileNimble(MCMC, project = model)
#  
#  niter <- 5000
#  nburnin <- 1000
#  nchains <- 3
#  
#  # RUN THE MCMC
#  MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
#                                                               niter = niter,
#                                                               nburnin = nburnin,
#                                                               nchains = nchains,
#                                                               samplesAsCodaMCMC = TRUE))

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  save(myNimbleOutput, niter, nburnin, nchains, MCMCRuntime,
#       file = file.path(baseDir,"nimbleSCR/inst/extdata/simulate_and_fit_SCR_models_samples.RData"))

## ----echo = FALSE-------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "simulate_and_fit_SCR_models_samples.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/simulate_and_fit_SCR_models_samples.RData"))
}

## -----------------------------------------------------------------------------
chainsPlot(myNimbleOutput, line = c(N, p0, N/M, sigma))

## ---- warning = FALSE, message = FALSE----------------------------------------
set.seed(1)
trapCovs <- cbind( runif(nimConstants$n.traps,-1, 1),
                   runif(nimConstants$n.traps,-1, 1))

## ---- warning = FALSE, message = FALSE----------------------------------------
nimData$trapCovs <-  trapCovs

## ---- warning = FALSE, message = FALSE----------------------------------------
nimInits$betaTraps <- c(-2, 2)

## ---- warning = FALSE, message = FALSE----------------------------------------
modelCodeTrap <- nimbleCode({
  ## priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  p0 ~ dunif(0, 1)
  
  #trap covariates 
  betaTraps[1] ~ dunif(-5,5)
  betaTraps[2] ~ dunif(-5,5)
  p0Traps[1:n.traps] <- ilogit(logit(p0) + betaTraps[1] * trapCovs[1:n.traps, 1] +
                                           betaTraps[2] * trapCovs[1:n.traps, 2])
    
  ## loop over individuals
  for(i in 1:M) {
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
    y[i, 1:lengthYCombined] ~ dbinomLocal_normal( size = trials[1:n.traps],
                                                  p0Traps = p0Traps[1:n.traps],
                                                  s = sxy[i,1:2],
                                                  sigma = sigma,
                                                  trapCoords = trapCoords[1:n.traps,1:2],
                                                  localTrapsIndices = trapIndex[1:n.cells,1:maxNBDets],
                                                  localTrapsNum = nTraps[1:n.cells],
                                                  habitatGrid = habitatIDDet[1:y.maxDet,1:x.maxDet],
                                                  indicator = z[i],
                                                  lengthYCombined = lengthYCombined)
   }
  ## derived quantity: total population size
  N <- sum(z[1:M])
})

## ---- warning = FALSE, message = FALSE----------------------------------------
model <- nimbleModel( code = modelCodeTrap,
                      constants = nimConstants,
                      data = nimData,
                      inits = nimInits,
                      check = F,       
                      calculate = T)  

## ---- warning = FALSE, message = FALSE----------------------------------------
# FIRST WE GET THE NODES TO SIMULATE
nodesToSim <- model$getDependencies(c("sxy", "z"), self=T)
# THEN WE SIMULATE THOSE NODES 
set.seed(100)
model$simulate(nodesToSim, includeData = FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
N <- sum(model$z)

## ---- warning = FALSE, message = FALSE----------------------------------------
i <- 5
#NUMBER OF DETECTIONS SIMULATED FOR INDIVIDUAL i
model$y[i,1]

#LET'S PLOT THEM
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch = 1, cex = 1.5) #pch=16) 
points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16 )  
# PLOT ITS ACTIVITY CETENR 
points(model$sxy[i, 2] ~ model$sxy[i, 1],col="orange", pch=16) 

whichdets <- model$y[i, (trapLocal$numLocalIndicesMax + 2):(trapLocal$numLocalIndicesMax+model$y[i, 1] + 1)]
points(ScaledtrapCoords[whichdets, "y"] ~
         ScaledtrapCoords[whichdets, "x"], col="blue", pch=16) 

par(xpd=TRUE)
legend(x = -1, y = 13,
       legend=c("Habitat grid centers", "Traps", "Activity center", "Detections"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1, 16, 16, 16),
       col=c("black", "red", "orange", "blue"),
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
nimData1 <- nimData
nimData1$y <- model$y
nimInits1 <- nimInits
nimInits1$z <- model$z
nimInits1$sxy <- model$sxy

# CREATE AND COMPILE THE NIMBLE MODEL
model <- nimbleModel( code = modelCodeTrap,
                      constants = nimConstants,
                      data = nimData1,
                      inits = nimInits1,
                      check = F,
                      calculate = F)
model$calculate()

## ----eval = FALSE-------------------------------------------------------------
#  cmodel <- compileNimble(model)
#  cmodel$calculate()

## -----------------------------------------------------------------------------
MCMCconf <- configureMCMC(model = model,
                          monitors = c("N", "sigma", "p0","psi", "betaTraps"),
                          control = list(reflective = TRUE),
                          thin = 1)

## Add block sampling of sxy coordinates see Turek et al 2021 for further details 
MCMCconf$removeSamplers("sxy")
ACnodes <- paste0("sxy[", 1:nimConstants$M, ", 1:2]")
for(node in ACnodes) {
  MCMCconf$addSampler(target = node,
                  type = "RW_block",
                  control = list(adaptScaleOnly = TRUE),
                  silent = TRUE)
}

MCMC <- buildMCMC(MCMCconf)

## ----eval = FALSE-------------------------------------------------------------
#  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
#  
#  # RUN THE MCMC
#  MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
#                                                               niter = niter,
#                                                               nburnin = nburnin,
#                                                               nchains = nchains,
#                                                               samplesAsCodaMCMC = TRUE))

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  save(myNimbleOutput, MCMCRuntime,
#       file = file.path(baseDir,"nimbleSCR/inst/extdata/simulate_and_fit_SCR_models_samples2.RData"))

## ----echo = FALSE-------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "simulate_and_fit_SCR_models_samples2.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/simulate_and_fit_SCR_models_samples2.RData"))
}

## -----------------------------------------------------------------------------
chainsPlot(myNimbleOutput, line = c(N, nimInits$betaTraps, p0, N/M, sigma))

