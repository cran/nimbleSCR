## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width = 8, fig.height = 8) 

## ---- warning = FALSE, message = FALSE----------------------------------------
# LOAD PACKAGES
library(coda)
library(nimble)
library(nimbleSCR)
library(basicMCMCplots)

## ---- warning = FALSE, message = FALSE, include=FALSE-------------------------
# FunD <- "E:/RovQuantSD_Codes"
# source(file.path(FunD, "nimbleSCR_HRfunctions/HRA_nimble_SD.R"))
# source(file.path(FunD, "nimbleSCR_HRfunctions/dbinomLocal_normalPlateau.R"))
if(Sys.info()['user'] == 'dturek') {
  baseDir <- '~/github/nimble/nimbleSCR/'                   ## Daniel
} else if(Sys.info()['user'] == 'pidu') {
  baseDir <- 'C:/Users/pidu/PACKAGES/nimbleSCR/'            ## Pierre
} else if(Sys.info()['user'] == 'cymi') {
  baseDir <- 'C:/Personal_Cloud/OneDrive/Work/nimbleSCR/'   ## Cyril
} else if(Sys.info()['user'] == 'arichbi') {
  baseDir <- 'C:/PROJECTS/nimbleSCR/'                       ## Richard
} else if(Sys.info()['user'] == 'admin') {                    ## Soumen
    baseDir <- '~/GitHubSD/nimbleSCR/'
} else baseDir <- NULL


## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 6---------
# CREATE HABITAT GRID 
coordsHabitatGridCenter <- cbind(rep(seq(29.5, 0.5, by=-1), 30),
                                 sort(rep(seq(0.5, 29.5, by=1), 30)))
colnames(coordsHabitatGridCenter) <- c("x","y")

# CREATE TRAP GRID
trapCoords <- cbind(rep(seq(5.5, 24.5,by=1),20),
               sort(rep(seq(5.5, 24.5,by=1),20)))
colnames(trapCoords) <- c("x","y")

# PLOT CHECK
plot(coordsHabitatGridCenter[,"y"] ~ coordsHabitatGridCenter[,"x"], pch = 1, cex = 1.5) #pch=16) 
points(trapCoords[,"y"] ~ trapCoords[,"x"], col="red", pch=16 ) 
par(xpd=TRUE)
legend(x = 7, y = 33.3,
       legend=c("Habitat grid centers", "Traps"),
       pt.cex = c(1.5,1),
       horiz = T,
       pch=c(1,16),
       col=c("black", "red"),
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
## RESCALE COORDINATES
ScaledCoords <- scaleCoordsToHabitatGrid(coordsData =  trapCoords,
                                         coordsHabitatGridCenter = coordsHabitatGridCenter)

lowerAndUpperCoords <- getWindowCoords(
  scaledHabGridCenter = ScaledCoords$coordsHabitatGridCenterScaled,
  scaledObsGridCenter = ScaledCoords$coordsDataScaled,
  plot.check = F)

## LOCAL EVALUATION
habitatMask <- matrix(1, nrow = 30, ncol= 30, byrow = TRUE)

trapLocal <- getLocalObjects(habitatMask = habitatMask,
                             coords = ScaledCoords$coordsDataScaled,
                             dmax = 10,
                             plot.check = FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
lengthYCombined <- 1 + 2*trapLocal$numLocalIndicesMax 

## ---- warning = FALSE, message = FALSE----------------------------------------
modelCode <- nimbleCode({
  ## priors
  psi ~ dunif(0, 1)
  sigma ~ dunif(0, 50)
  w ~ dunif(0, 20)
  p0 ~ dunif(0, 1)
  
  ## loop over individuals
  for(i in 1:M) {
    ## AC coordinates
    sxy[i, 1:2] ~ dbernppAC(
      lowerCoords = lowerHabCoords[1:numHabWindows, 1:2],
      upperCoords = upperHabCoords[1:numHabWindows, 1:2],
      logIntensities = logHabIntensity[1:numHabWindows],
      logSumIntensity = logSumHabIntensity,
      habitatGrid = habitatGrid[1:numGridRows,1:numGridCols],
      numGridRows =  numGridRows,
      numGridCols = numGridCols
    )
    
    ## latent dead/alive indicators
    z[i] ~ dbern(psi)
    ## likelihood
    y[i, 1:lengthYCombined] ~ dbinomLocal_normalPlateau(
      size = trials[1:n.traps],
      p0 = p0,
      s = sxy[i,1:2],
      sigma = sigma,
      w = w,
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
p0 <- 0.25
sigma <- 1
w <- 1.5
psi <- 0.5

## ---- warning = FALSE, message = FALSE----------------------------------------
M <- 500

## ---- warning = FALSE, message = FALSE----------------------------------------
nimConstants <- list(M = M,
                     n.traps = dim(ScaledCoords$coordsDataScaled)[1],
                     y.maxDet = dim(trapLocal$habitatGrid)[1],
                     x.maxDet = dim(trapLocal$habitatGrid)[2],
                     n.cells = dim(trapLocal$localIndices)[1],
                     maxNBDets = trapLocal$numLocalIndicesMax,
                     nTraps = trapLocal$numLocalIndices,
                     lengthYCombined = lengthYCombined,
                     numHabWindows = dim(lowerAndUpperCoords$upperHabCoords)[1],
                     numGridRows = dim(lowerAndUpperCoords$habitatGrid)[1],
                     numGridCols = dim(lowerAndUpperCoords$habitatGrid)[2]
                     )

habIntensity = rep(1, dim(lowerAndUpperCoords$upperHabCoords)[1])
logSumHabIntensity <- log(sum(habIntensity))
logHabIntensity <- log(habIntensity)
               
nimData <- list(trapCoords = ScaledCoords$coordsDataScaled,
                trapIndex = trapLocal$localIndices,
                lowerHabCoords = lowerAndUpperCoords$lowerHabCoords,
                upperHabCoords = lowerAndUpperCoords$upperHabCoords,
                habitatGrid = lowerAndUpperCoords$habitatGrid,
                logHabIntensity = logHabIntensity,
                logSumHabIntensity = logSumHabIntensity,
                habitatIDDet = trapLocal$habitatGrid,
                trials = rep(1, dim(ScaledCoords$coordsDataScaled)[1])
)

## ---- warning = FALSE, message = FALSE----------------------------------------
nimInits <- list(p0 = p0,
                 psi = psi,
                 sigma = sigma,
                 w = w)

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
set.seed(1234)
model$simulate(nodesToSim, includeData = FALSE)

## ---- warning = FALSE, message = FALSE----------------------------------------
N <- sum(model$z)

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

MCMCconf <- configureMCMC(model = model,
                          monitors = c("N", "sigma", "p0","w","psi"),
                          control = list(reflective = TRUE),
                          thin = 1)

samplerConfList <- MCMCconf$getSamplers()
print(samplerConfList[1:5])

## ----eval = FALSE-------------------------------------------------------------
#  cmodel <- compileNimble(model)

## ---- warning = FALSE, message = FALSE----------------------------------------
# sigma
control <- samplerConfList[[2]]$control
control$log <- F
control$reflective <- T
control$adaptive <- T
control$scale <- 1
samplerConfList[[2]]$setControl(control)

# w
control <- samplerConfList[[3]]$control
control$log <- F
control$reflective <- T
control$adaptive <- T
control$scale <- 1
samplerConfList[[3]]$setControl(control)

# Use this modified list of samplerConf objects in the MCMC configuration
MCMCconf$setSamplers(samplerConfList)

# Retrieve the current ordering of sampler execution
ordering <- MCMCconf$getSamplerExecutionOrder()
len <- length(ordering)
MCMCconf$setSamplerExecutionOrder(c(ordering[1], 
                                    rep(ordering[2], 10), # sigma
                                    rep(ordering[3], 10), # w
                                    ordering[4:len]))

## ---- eval = FALSE, warning = FALSE, message = FALSE--------------------------
#  MCMC <- buildMCMC(MCMCconf)

## ---- eval = FALSE, warning = FALSE, message = FALSE--------------------------
#  cMCMC <- compileNimble(MCMC, project = model, resetFunctions = TRUE)
#  
#  # RUN THE MCMC
#  niter <- 10000
#  burnin <- 2000
#  nchains <- 3
#  MCMCRuntime <- system.time(myNimbleOutput <- runMCMC( mcmc = cMCMC,
#                                                               nburnin = burnin,
#                                                               niter = niter,
#                                                               nchains = nchains,
#                                                               samplesAsCodaMCMC = TRUE))
#  

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  save(myNimbleOutput, niter, burnin, nchains, MCMCRuntime,
#       file = file.path(baseDir,"nimbleSCR/inst/extdata/dbinomLocal_normalPlateau_samples.RData"))

## ----echo = FALSE-------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "dbinomLocal_normalPlateau_samples.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/dbinomLocal_normalPlateau_samples.RData"))
}

## ---- warning = FALSE, message = FALSE, fig.width = 7, fig.height = 10--------
print(MCMCRuntime)
#plot check 
chainsPlot(myNimbleOutput, line = c(N, p0, N/M, sigma, w))

## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4---------
## HALF-NORMAL PLATEAU FUNCTION
detfun.HNP <- function(D, p0, sigma, w){
  bool <- D <= w
  which_not_bool <- which(!bool)
  dens <- numeric(length = length(D))  
  dens[bool] <- 1
  adjustedD2 <- (D[which_not_bool]-w)^2
  dens[which_not_bool] <- exp(-(adjustedD2)/(2.0*sigma*sigma))
  return(p0*dens)
}
par(mfrow=c(1,1))
x <- seq(0, 15, len=10001)
plot(x, detfun.HNP(D=x, p0=0.2, sigma = 1, w = 1), type = "l", lwd=2, col = "orange",
     main = "Half-normal plateau detection function", xlab = 'd', ylab = 'Detection probability') 
lines(x, detfun.HNP(D=x, p0=0.1, sigma = 1.5, w = 3), lwd=2, lty = "solid", col = "darkcyan")
focal.vars <- c('p0', 'sigma', 'w') 
focal.values <- cbind(c(0.2, 0.1), # p0
                      c(1, 1.5), # sigma
                      c(1, 3)) # w 
legend.items <- apply(focal.values, 1, 
                      function(s){paste(paste(focal.vars,"=",s,sep=" "),collapse=", ")})
legend("topright",
       legend = legend.items,
       lty = c('solid', 'solid'), 
       bg="transparent",
       col=c('orange', "darkcyan"),
       cex = 0.75,
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
prob <- 0.95
paramnames.hr <- c("HRradius", "HRarea")
params <- c(sigma, w)
names(params) <- c("sigma", "w")
HRAnim <- getHomeRangeArea( x = params, detFun = 1, prob = prob, d = 6, 
                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
                      nBreaks = 800, tol = 1E-5, nIter = 2000)

# Different values of argument "detFun"
# 0 = Half-normal, 1 = Half-normal plateau, 2 = Exponential,
# 3 = Aysmmetric logistic, 4 = Bimodal, 5 = Donut.
HR.hnp <- c(HRAnim$run())
names(HR.hnp) <- paramnames.hr
print(HR.hnp)


## ---- eval = FALSE, warning = FALSE, message = FALSE--------------------------
#  myNimbleOutput <- as.mcmc.list(myNimbleOutput)
#  if(is.list(myNimbleOutput)){posteriorSamples <- do.call(rbind, myNimbleOutput)}
#  if(!is.list(myNimbleOutput)){posteriorSamples <- as.matrix(myNimbleOutput)}
#  theseIter <- round(seq(1, nrow(posteriorSamples), by = 10))
#  posteriorSamples <- posteriorSamples[theseIter,names(params)]
#  HRAnim.mat <- getHomeRangeArea(x = posteriorSamples, detFun = 1, prob = prob, d = 6,
#                           xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
#                           nBreaks = 800, tol = 1E-5, nIter = 2000)
#  
#  cHRAnim.arr <- compileNimble(HRAnim.mat, resetFunctions = TRUE)
#  
#  HRA.Runtime <- system.time(
#    HR.chain <- cHRAnim.arr$run()
#  )

## ----eval = FALSE, echo = FALSE-----------------------------------------------
#  save(HR.chain, HRA.Runtime,
#       file = file.path(baseDir,"nimbleSCR/inst/extdata/HRA_chain_runtime.RData"))

## ----echo = FALSE-------------------------------------------------------------
if(is.null(baseDir)) {
    load(system.file("extdata", "HRA_chain_runtime.RData", package = "nimbleSCR"))
} else {
    load(file.path(baseDir,"nimbleSCR/inst/extdata/HRA_chain_runtime.RData"))
}

## -----------------------------------------------------------------------------
print(HRA.Runtime)
dimnames(HR.chain)[[2]] <- paramnames.hr

## ---- warning = FALSE, message = FALSE----------------------------------------
HRest <- do.call(rbind, lapply(c(1:2), function(j){
  c(mean(HR.chain[,j], na.rm = T), sd(HR.chain[,j], na.rm = T), quantile(HR.chain[,j], probs = c(0.025, 0.975), na.rm = T))
}))
dimnames(HRest) <- list(paramnames.hr, c("MEAN", "SD", "2.5%", "97.5%"))

cat("Numerical estimates using MCMC samples: \n", sep = "")
print(HRest)

## ---- warning = FALSE, message = FALSE, fig.width = 7, fig.height = 5---------
HR.mcmc.sample <- as.mcmc.list(lapply(1:nchains, function(x) {
  niterPerChain <- nrow(HR.chain) / nchains
  theseRows <- ((x - 1) * niterPerChain + 1):(x * niterPerChain)
  this.MCMC.chain <- mcmc(HR.chain[theseRows, ])
  return(this.MCMC.chain)
})) 
names(HR.mcmc.sample) <- paste0("chain", 1:nchains)
chainsPlot(HR.mcmc.sample, line = c(HR.hnp))

## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4---------
## HALF-NORMAL FUNCTION
detfun.HN <- function(D, p0, sigma){ p0*exp(-D*D/(2*sigma*sigma))}
par(mfrow=c(1,1))
x <- seq(0, 15, len=10001)
plot(x, detfun.HN(D=x, p0=0.25, sigma = 2), type = "l", lwd=2, col = "orange", main = "Half-normal detection function", xlab = 'd', ylab = 'Detection probability') 
lines(x, detfun.HN(D=x, p0=0.15, sigma = 3), lwd=2, lty = "solid", col = "darkcyan")
focal.vars <- c('p0', 'sigma') 
focal.values <- cbind(c(0.25, 0.15), # p0
                      c(2, 3)) # sigma
                      
legend.items <- apply(focal.values, 1, 
                      function(s){paste(paste(focal.vars,"=",s,sep=" "),collapse=", ")})
legend("topright",
       legend = legend.items,
       lty = c('solid', 'solid'), 
       bg="transparent",
       col=c('orange', "darkcyan"),
       cex = 0.75,
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
# Half-normal 
sigma = 2
params <- c(sigma)
names(params) <- c("sigma")

HRAnim <- getHomeRangeArea(x = params, detFun = 0, prob = prob, d = 6, 
                      xlim = c(0, dim(habitatMask)[2]), ylim = c(0, dim(habitatMask)[1]),
                      nBreaks = 800, tol = 1E-5, nIter = 2000)

HR.hn <- c(HRAnim$run())
names(HR.hn) <- paramnames.hr
print(HR.hn)

## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4---------
## EXPONENTIAL FUNCTION
detfun.EXP <- function(D, p0, rate){ p0*exp(-rate*D)}
par(mfrow=c(1,1))
x <- seq(0, 15, len=10001)
plot(x, detfun.EXP(D=x, p0=0.25, rate = 1/2), type = "l", lwd=2, col = "orange", main = "Exponential detection function", xlab = 'd', ylab = 'Detection probability') 
lines(x, detfun.EXP(D=x, p0=0.15, rate = 1/3), lwd=2, lty = "solid", col = "darkcyan")
focal.vars <- c('p0', 'rate') 
focal.values <- cbind(c(0.25, 0.15), # p0
                      c(1/2, round(1/3,2))) # rate
                     
legend.items <- apply(focal.values, 1, 
                      function(s){paste(paste(focal.vars,"=",s,sep=" "),collapse=", ")})
legend("topright",
       legend = legend.items,
       lty = c('solid', 'solid'), 
       bg="transparent",
       col=c('orange', "darkcyan"),
       cex = 0.75,
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
# Exponential
rate = 1/2
params <- c(rate)
names(params) <- c("rate")
HRAnim <- getHomeRangeArea(x = params, detFun = 2, prob = prob, d = 6, 
                      xlim = c(0, nimConstants$x.max), ylim = c(0, nimConstants$y.max),
                      nBreaks = 800, tol = 1E-5, nIter = 2000)
HR.exp <- c(HRAnim$run())
names(HR.exp) <- paramnames.hr
print(HR.exp)

## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4---------
## ASYMMETRIC LOGISTICS FUNCTION
detfun.AL <- function(D, p0, sigma, alpha.a, alpha.b){
  Anuslope <- (2*abs(alpha.a)*alpha.b)/(1+alpha.b)
  fx <- 1/ (1+(D/sigma)^alpha.b)
  den <- 1+fx*((D/sigma)^(alpha.a))+(1-fx)*((D/sigma)^(alpha.a*alpha.b))
  detp <- p0/den
  return(detp)
}

par(mfrow=c(1,1))
x <- seq(0, 15, len=10001)
plot(x, detfun.AL(D=x, p0=0.3, sigma = 2, alpha.a = 5, alpha.b = 1), type = "l", lwd=2, col = "orange", main = "Asymmetric logistic detection function", xlab = 'd', ylab = 'Detection probability') 
lines(x, detfun.AL(D=x, p0=0.15, sigma = 3, alpha.a = 10, alpha.b = 1), lwd=2, lty = "solid", col = "darkcyan")
focal.vars <- c('p0', 'sigma', 'alpha.a', 'alpha.b') 
focal.values <- cbind(c(0.3, 0.15), # p0
                      c(2, 3), # sigma
                      c(5, 10), # alpha.a
                      c(1, 1)) # alpha.b
                     
legend.items <- apply(focal.values, 1, 
                      function(s){paste(paste(focal.vars,"=",s,sep=" "),collapse=", ")})
legend("topright",
       legend = legend.items,
       lty = c('solid', 'solid'), 
       bg="transparent",
       col=c('orange', "darkcyan"),
       cex = 0.75,
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
# Asymmetric logistic
sigma = 2
alpha.a = 5 
alpha.b = 1
params <- c(sigma, alpha.a, alpha.b)
names(params) <- c("sigma", "alpha.a", "alpha.b")
HRAnim <- getHomeRangeArea(x = params, detFun = 3, prob = prob, d = 6, 
                      xlim = c(0, nimConstants$x.max), ylim = c(0, nimConstants$y.max),
                      nBreaks = 800, tol = 1E-5, nIter = 2000)
HR.al <- c(HRAnim$run())
names(HR.al) <- paramnames.hr
print(HR.al)

## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4---------
## BIMODAL FUNCTION
detfun.BI <- function(D, p0.a, sigma.a, p0.b, sigma.b,  w){
  densa <- p0.a *  exp(-D*D/(2.0*sigma.a*sigma.a)) 
  densb <-  p0.b * exp(-(D-w)*(D-w)/(2.0*sigma.b*sigma.b))
  dens <- densa + densb
  return(dens)
}

par(mfrow=c(1,1))
x <- seq(0, 15, len=10001)
plot(x, detfun.BI(D=x, p0.a=0.25, sigma.a = 0.5, p0.b = 0.15, sigma.b = 1, w = 2), type = "l", lwd=2, col = "orange", main = "Bimodal detection function", xlab = 'd', ylab = 'Detection probability') 
lines(x, detfun.BI(D=x, p0.a=0.05, sigma.a = 1.5, p0.b = 0.1, sigma.b = 0.5, w = 3), lwd=2, lty = "solid", col = "darkcyan")
focal.vars <- c('p0.a', 'sigma.a', 'p0.b', 'sigma.b', 'w') 
focal.values <- cbind(c(0.25, 0.05), # p0.a
                      c(0.5, 1.5), # sigma.a
                      c(0.15, 0.1), # p0.a
                      c(1, 0.5), # sigma.b
                      c(2, 3)) # w 
legend.items <- apply(focal.values, 1, 
                      function(s){paste(paste(focal.vars,"=",s,sep=" "),collapse=", ")})
legend("topright",
       legend = legend.items,
       lty = c('solid', 'solid'), 
       bg="transparent",
       col=c('orange', "darkcyan"),
       cex = 0.75,
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
# Bimodal 
p0.a = 0.25
sigma.a = 0.5
p0.b = 0.15 
sigma.b = 1
w = 2
params <- c(sigma.a, sigma.b, p0.a, p0.b, w)
names(params) <- c("sigma.a", "sigma.b", "p0.a", "p0.b", "w")
HRAnim <- getHomeRangeArea(x = params, detFun = 4, prob = prob, d = 6, 
                      xlim = c(0, nimConstants$x.max), ylim = c(0, nimConstants$y.max),
                      nBreaks = 800, tol = 1E-5, nIter = 2000)
HR.bi <- c(HRAnim$run())
names(HR.bi) <- paramnames.hr
print(HR.bi)

## ---- warning = FALSE, message = FALSE, fig.width = 6, fig.height = 4---------
## DONUT FUNCTION
detfun.DN <- function(D, p0, sigma.a, sigma.b, w){
  bool <- D <= w
  which_not_bool <- which(!bool)
  dens <- numeric(length = length(D)) 
  adjustedD2.bool <- (D[bool]-w)^2
  adjustedD2.notbool <- (D[which_not_bool]-w)^2
  dens[bool] <- exp(-adjustedD2.bool/(2.0*sigma.a*sigma.a))
  dens[which_not_bool] <- exp(-adjustedD2.notbool/(2.0*sigma.b*sigma.b)) 
  return(p0*dens)
}

par(mfrow=c(1,1))
x <- seq(0, 15, len=10001)
plot(x, detfun.DN(D=x, p0=0.2, sigma.a = 1.5, sigma.b = 1, w = 1), type = "l", lwd=2, col = "orange",
     main = "Donut detection function", xlab = 'd', ylab = 'Detection probability') 
lines(x, detfun.DN(D=x, p0=0.1, sigma.a = 2, sigma.b = 1.5, w = 3), lwd=2, lty = "solid", col = "darkcyan")
focal.vars <- c('p0', 'sigma.a', 'sigma.b', 'w') 
focal.values <- cbind(c(0.2, 0.1), # p0
                      c(1.5, 2), # sigma.a
                      c(1, 1.5), # sigma.b
                      c(1, 3)) # w 
legend.items <- apply(focal.values, 1, 
                      function(s){paste(paste(focal.vars,"=",s,sep=" "),collapse=", ")})
legend("topright",
       legend = legend.items,
       lty = c('solid', 'solid'), 
       bg="transparent",
       col=c('orange', "darkcyan"),
       cex = 0.75,
       bty = 'n')

## ---- warning = FALSE, message = FALSE----------------------------------------
# Donut 
sigma.a = 1.5 
sigma.b = 1 
w = 1 
params <- c(sigma.a, sigma.b, w)
names(params) <- c("sigma.a", "sigma.b", "w")
HRAnim <- getHomeRangeArea(x = params, detFun = 5, prob = prob, d = 6, 
                      xlim = c(0, nimConstants$x.max), ylim = c(0, nimConstants$y.max),
                      nBreaks = 800, tol = 1E-5, nIter = 2000)
HR.dn <- c(HRAnim$run())
names(HR.dn) <- paramnames.hr
print(HR.dn)

## ---- warning = FALSE, message = FALSE----------------------------------------
sigma <- 2
q<-qchisq(prob,2)
radius<-sigma*sqrt(q)
area<-pi*(radius^2)
HR.true<-c(radius,area)
names(HR.true) <- paramnames.hr
cat("Analytical estimates: \n", sep = "")
print(HR.true)

