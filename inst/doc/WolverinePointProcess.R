## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 5, fig.height = 5) 

## ---- echo = FALSE------------------------------------------------------------
if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else if(Sys.info()['user'] == 'weizhang'){               
    baseDir <- '/Users/weizhang/Documents/GitHub/nimbleSCR/' ## Wei 
} else baseDir <- NULL

if(is.null(baseDir)) {
    load(system.file("extdata", "WolverinePointProcess_data.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/WolverinePointProcess_data.RData"))
}

## ---- warning = FALSE, message = FALSE----------------------------------------
library(nimble)
library(nimbleSCR)
library(basicMCMCplots)
library(coda)

## -----------------------------------------------------------------------------
modelCode1 <- nimbleCode({  
  ##------ SPATIAL PROCESS 
  ## Intercept and slope for the log-linear model for habitat selection intensity
  habCoeffInt ~ dnorm(0, sd = 10)
  habCoeffSlope ~ dnorm(0, sd = 10)
  ## Habitat intensity for each habitat window
  habIntensity[1:numHabWindows] <- exp(habCoeffInt + habCoeffSlope * habCovs[1:numHabWindows])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  ## Activity centres of the observed individuals: a Bernoulli point process
  for(i in 1:numIdDetected){
    sxy[i,1:2] ~ dbernppAC(
      lowerCoords = habLoCoords[1:numHabWindows, 1:2],
      upperCoords = habUpCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i
  
  ##----- DEMOGRAPHIC PROCESS
  ## Number of individuals in the population
  N ~ dpois(sumHabIntensity)
  ## Number of detected individuals
  nDetectedIndiv ~ dbin(probDetection, N)

  ##----- DETECTION PROCESS 
  ## Scale for the multivariate normal detection function
  sigma ~ dunif(0,10)
  ## Intercept and slope parameters for the log-linear model for detection intensity
  for(c in 1:numCounties){
    detCoeffInt[c] ~ dnorm(0, sd = 10)
  }#c
  for(cc in 1:numDetCovs){
    detCoeffSlope[cc] ~ dnorm(0, sd = 10)
  }#cc
  ## Baseline detection intensity for each detection window
  for(j in 1:numDetWindows){
    detIntensity[j] <- exp( detCoeffInt[detCounties[j]] +
                              detCoeffSlope[1] * detCovs[j,1] +
                              detCoeffSlope[2] * detCovs[j,2] +
                              detCoeffSlope[3] * detCovs[j,3])
  }#j
  ## Detections of the observed individuals conditional on their activity centers
  for(i in 1:numIdDetected) {
    y[i,1:(maxDetections+1),1:3] ~ dpoisppDetection_normal(
      lowerCoords = detLoCoords[1:numDetWindows,1:2],
      upperCoords = detUpCoords[1:numDetWindows,1:2],
      s = sxy[i,1:2],
      sd = sigma,
      baseIntensities = detIntensity[1:numDetWindows],
      numMaxPoints = maxDetections,
      numWindows = numDetWindows,
      indicator = 1)
  }#i

  ## The probability that an individual in the population is detected at least once
  ## i.e. one minus the void probability over all detection windows
  probDetection <- 1 - marginalVoidProbNumIntegration(
    quadNodes = quadNodes[1:maxNumNodes,1:2,1:numHabWindows],
    quadWeights = quadWeights[1:numHabWindows],
    numNodes = numNodes[1:numHabWindows],
    lowerCoords = detLoCoords[1:numDetWindows,1:2],
    upperCoords = detUpCoords[1:numDetWindows,1:2],
    sd = sigma,
    baseIntensities = detIntensity[1:numDetWindows],
    habIntensities = habIntensity[1:numHabWindows],
    sumHabIntensity = sumHabIntensity,
    numObsWindows = numDetWindows,
    numHabWindows = numHabWindows
  )
  logDetProb <- log(probDetection)
  normData ~ dnormalizer(logNormConstant = -numIdDetected*logDetProb)

})

## -----------------------------------------------------------------------------
modelCode2 <- nimbleCode({
  ##------ SPATIAL PROCESS
  habCoeffSlope ~ dnorm(0, sd = 10)
  ## Habitat intensity for each habitat window
  habIntensity[1:numHabWindows] <- exp(habCoeffSlope * habCovs[1:numHabWindows])
  for(i in 1:M){
    sID[i] ~ dcat(habIntensity[1:numHabWindows])
  }#i

  ##----- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  for(i in 1:M){
    z[i] ~ dbern(psi)
  }#i
  ## Number of individuals in the population
  N <- sum(z[1:M])

  ##----- DETECTION PROCESS
  ## Scale for the multivariate normal detection function
  sigma ~ dunif(0,10)
  ## Intercept and slope parameters for the log-linear model for detection intensity
  for(c in 1:numCounties){
    detCoeffInt[c] ~ dnorm(0, sd = 10)
  }#c
  for(cc in 1:numDetCovs){
    detCoeffSlope[cc] ~ dnorm(0, sd = 10)
  }#cc
  ## Baseline detection intensity for each detection window
  for(j in 1:numDetWindows){
    lambdaTraps[j] <- exp( detCoeffInt[detCounties[j]] +
                             detCoeffSlope[1] * detCovs[j,1] +
                             detCoeffSlope[2] * detCovs[j,2] +
                             detCoeffSlope[3] * detCovs[j,3])
  }#j

  ## Detections of the observed individuals conditional on their activity centers
  for(i in 1:M) {
    y[i,1:lengthYCombined] ~ dpoisLocal_normal(
      lambdaTraps = lambdaTraps[1:numDetWindows],
      trapCoords = detCoords[1:numDetWindows,1:2],
      s = habCoords[sID[i],1:2],
      sigma = sigma,
      localTrapsIndices = localTrapsIndices[1:numHabWindows,1:numDetWindows],
      localTrapsNum = localTrapsNum[1:numHabWindows],
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      lengthYCombined = lengthYCombined,
      indicator = z[i])
  }#i
})

## -----------------------------------------------------------------------------
modelCode3 <- nimbleCode({
  ##------ SPATIAL PROCESS
  habCoeffSlope ~ dnorm(0, sd = 10)
  ## Habitat intensity for each habitat window
  habIntensity[1:numHabWindows] <- exp(habCoeffSlope * habCovs[1:numHabWindows])
  sumHabIntensity <- sum(habIntensity[1:numHabWindows])
  logHabIntensity[1:numHabWindows] <- log(habIntensity[1:numHabWindows])
  logSumHabIntensity <- log(sumHabIntensity)
  ## Activity centres of the observed individuals: a bernoulli point process
  for(i in 1:M){
    sxy[i,1:2] ~ dbernppAC(
      lowerCoords = habLoCoords[1:numHabWindows,1:2],
      upperCoords = habUpCoords[1:numHabWindows,1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:y.max,1:x.max],
      numGridRows = y.max,
      numGridCols = x.max)
  }#i

  ##----- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  for(i in 1:M){
    z[i] ~ dbern(psi)
  }#i
  ## Number of individuals in the population
  N <- sum(z[1:M])

  ##----- DETECTION PROCESS
  ## Scale for the half-normal detection function
  sigma ~ dunif(0,10)
  ## Intercept and slope parameters for the log-linear model for detection intensity
  for(c in 1:numCounties){
    detCoeffInt[c] ~ dnorm(0, sd = 10)
  }#c
  for(cc in 1:numDetCovs){
    detCoeffSlope[cc] ~ dnorm(0, sd = 10)
  }#cc
  ## Baseline detection intensity for each detection window
  for(j in 1:numDetWindows){
    detIntensity[j] <- exp( detCoeffInt[detCounties[j]] +
                              detCoeffSlope[1] * detCovs[j,1] +
                              detCoeffSlope[2] * detCovs[j,2] +
                              detCoeffSlope[3] * detCovs[j,3])
  }#j

  ## Detections of the observed individuals conditional on their activity centers
  for(i in 1:M){
    y[i,1:(maxDetections+1),1:3] ~ dpoisppDetection_normal(
      lowerCoords = detLoCoords[1:numDetWindows,1:2],
      upperCoords = detUpCoords[1:numDetWindows,1:2],
      s = sxy[i,1:2],
      sd = sigma,
      baseIntensities = detIntensity[1:numDetWindows],
      numMaxPoints = maxDetections,
      numWindows = numDetWindows,
      indicator = z[i])
  }#i
})

## -----------------------------------------------------------------------------
modelCode4 <- nimbleCode({
  ##------ SPATIAL PROCESS
  habCoeffSlope ~ dnorm(0, sd = 10)
  ## Habitat intensity for each habitat window
  habIntensity[1:numHabWindows] <- exp(habCoeffSlope * habCovs[1:numHabWindows])
  for(i in 1:M){
    sID[i] ~ dcat(habIntensity[1:numHabWindows])
  }#i

  ##----- DEMOGRAPHIC PROCESS
  psi ~ dunif(0,1)
  for(i in 1:M){
    z[i] ~ dbern(psi)
  }#i
  ## Number of individuals in the population
  N <- sum(z[1:M])

  ##----- DETECTION PROCESS
  ## Scale for the half-normal detection function
  sigma ~ dunif(0,10)
  ## Intercept and slope parameters for the log-linear model for detection intensity
  for(c in 1:numCounties){
    detCoeffInt[c] ~ dnorm(0, sd = 10)
  }#c
  for(cc in 1:numDetCovs){
    detCoeffSlope[cc] ~ dnorm(0, sd = 10)
  }#cc
  ## Baseline detection intensity for each detection window
  for(j in 1:numDetWindows){
    detIntensity[j] <- exp( detCoeffInt[detCounties[j]] +
                              detCoeffSlope[1] * detCovs[j,1] +
                              detCoeffSlope[2] * detCovs[j,2] +
                              detCoeffSlope[3] * detCovs[j,3])
  }#j

  ## Detections of the observed individuals conditional on their activity centers
  for(i in 1:M){
    y[i,1:(maxDetections+1),1:3] ~ dpoisppDetection_normal(
      lowerCoords = detLoCoords[1:numDetWindows,1:2],
      upperCoords = detUpCoords[1:numDetWindows,1:2],
      s = habCoords[sID[i],1:2],
      sd = sigma,
      baseIntensities = detIntensity[1:numDetWindows],
      numMaxPoints = maxDetections,
      numWindows = numDetWindows,
      indicator = z[i])
  }#i
})

## ---- eval = FALSE------------------------------------------------------------
#  load("WolverinePointProcess_Data.RData")

## ----eval = FALSE-------------------------------------------------------------
#  Rmodel1 <- nimbleModel(modelCode1, constants1, data1, inits1)
#  Rmodel2 <- nimbleModel(modelCode2, constants2, data2, inits2)
#  Rmodel3 <- nimbleModel(modelCode3, constants3, data3, inits3)
#  Rmodel4 <- nimbleModel(modelCode4, constants4, data4, inits4)

## ----message = FALSE, eval = FALSE--------------------------------------------
#  conf1 <- configureMCMC( Rmodel1,
#                         monitors = params1,
#                         print = FALSE)
#  conf2 <- configureMCMC( Rmodel2,
#                         monitors = params2,
#                         print = FALSE)
#  conf3 <- configureMCMC( Rmodel3,
#                         monitors = params3,
#                         print = FALSE)
#  conf4 <- configureMCMC( Rmodel4,
#                         monitors = params4,
#                         print = FALSE)
#  Rmcmc1 <- buildMCMC(conf1)
#  Rmcmc2 <- buildMCMC(conf2)
#  Rmcmc3 <- buildMCMC(conf3)
#  Rmcmc4 <- buildMCMC(conf4)

## ----eval = FALSE-------------------------------------------------------------
#  Cmodel1 <- compileNimble(Rmodel1)
#  Cmcmc1 <- compileNimble(Rmcmc1, project = Rmodel1)
#  MCMC_runtime1 <- system.time(
#    MCMC_samples1 <- runMCMC(Cmcmc1,
#                             niter = 11000,
#                             nburnin = 1000,
#                             nchains = 4))
#  
#  Cmodel2 <- compileNimble(Rmodel2)
#  Cmcmc2 <- compileNimble(Rmcmc2, project = Rmodel2)
#  MCMC_runtime2 <- system.time(
#    MCMC_samples2 <- runMCMC(Cmcmc2,
#                             niter = 11000,
#                             nburnin = 1000,
#                             nchains = 4))
#  
#  Cmodel3 <- compileNimble(Rmodel3)
#  Cmcmc3 <- compileNimble(Rmcmc3, project = Rmodel3)
#  MCMC_runtime3 <- system.time(
#    MCMC_samples3 <- runMCMC(Cmcmc3,
#                             niter = 11000,
#                             nburnin = 1000,
#                             nchains = 4))
#  
#  Cmodel4 <- compileNimble(Rmodel4)
#  Cmcmc4 <- compileNimble(Rmcmc4, project = Rmodel4)
#  MCMC_runtime4 <- system.time(
#    MCMC_samples4 <- runMCMC(Cmcmc4,
#                             niter = 11000,
#                             nburnin = 1000,
#                             nchains = 4))
#  
#  MCMC_samples1_summary <- samplesSummary(do.call(rbind,MCMC_samples1))
#  MCMC_samples2_summary <- samplesSummary(do.call(rbind,MCMC_samples2))
#  MCMC_samples3_summary <- samplesSummary(do.call(rbind,MCMC_samples3))
#  MCMC_samples4_summary <- samplesSummary(do.call(rbind,MCMC_samples4))
#  
#  MCMC_samples1_ess <- effectiveSize(MCMC_samples1)
#  MCMC_samples2_ess <- effectiveSize(MCMC_samples2)
#  MCMC_samples3_ess <- effectiveSize(MCMC_samples3)
#  MCMC_samples4_ess <- effectiveSize(MCMC_samples4)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  ##save(MCMC_samples1, MCMC_samples2, MCMC_samples3, MCMC_samples4,
#  ##     MCMC_runtime1, MCMC_runtime2, MCMC_runtime3, MCMC_runtime4,
#  ##     file = file.path(baseDir,
#  ##                      "nimbleSCR/inst/extdata/WolverinePointProcess_samples.RData"))
#  
#  save(MCMC_samples1_summary, MCMC_samples2_summary, MCMC_samples3_summary, MCMC_samples4_summary,
#       MCMC_samples1_ess, MCMC_samples2_ess, MCMC_samples3_ess, MCMC_samples4_ess,
#       MCMC_runtime1, MCMC_runtime2, MCMC_runtime3, MCMC_runtime4,
#       file = file.path(baseDir,
#                        "nimbleSCR/inst/extdata/WolverinePointProcess_summary.RData"))

## ----echo = FALSE-------------------------------------------------------------
##if(is.null(baseDir)) {
##    load(system.file("extdata",
##                     "WolverinePointProcess_samples.RData",
##                     package = "nimbleSCR"))
##} else {
##    load(file.path(baseDir,"nimbleSCR/inst/extdata/WolverinePointProcess_samples.RData"))
##}

if(is.null(baseDir)) {
    load(system.file("extdata",
                     "WolverinePointProcess_summary.RData",
                     package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/WolverinePointProcess_summary.RData"))
}

## -----------------------------------------------------------------------------
round(MCMC_samples1_summary, 2)
round(MCMC_samples2_summary, 2)
round(MCMC_samples3_summary, 2)
round(MCMC_samples4_summary, 2)

## ----echo = FALSE, eval = FALSE-----------------------------------------------
#  chainsPlot(MCMC_samples1, var = c("N","sigma"))
#  chainsPlot(MCMC_samples2, var = c("N","sigma"))
#  chainsPlot(MCMC_samples3, var = c("N","sigma"))
#  chainsPlot(MCMC_samples4, var = c("N","sigma"))

## -----------------------------------------------------------------------------
round(MCMC_samples1_ess,2)[c("N","sigma")] 
round(MCMC_samples2_ess,2)[c("N","sigma")] 
round(MCMC_samples3_ess,2)[c("N","sigma")] 
round(MCMC_samples4_ess,2)[c("N","sigma")] 

## -----------------------------------------------------------------------------
round(MCMC_samples1_ess,2)["N"]/round(MCMC_runtime1[3],1)
round(MCMC_samples2_ess,2)["N"]/round(MCMC_runtime2[3],1)
round(MCMC_samples3_ess,2)["N"]/round(MCMC_runtime3[3],1)
round(MCMC_samples4_ess,2)["N"]/round(MCMC_runtime4[3],1)

