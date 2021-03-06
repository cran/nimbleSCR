.onAttach <- function(libname, pkgname) {
    
    packageStartupMessage("Warning message:\n'getLocalTraps' and 'dbinom_sparseLocalSCR' are deprecated.\nUse 'getLocalObjects' and 'dbinomLocal_normal' instead.")
    
    suppressMessages({

        registerDistributions(
            list(
                dDispersal_exp = list(
                    BUGSdist = 'dDispersal_exp(s, rate)',
                    types = c('value = double(1)', 's = double(1)', 'rate = double()'),
                    discrete = FALSE,
                    mixedSizes = FALSE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)


        registerDistributions(
            list(
                dHabitatMask = list(
                    BUGSdist = 'dHabitatMask(s, xmax, xmin, ymax, ymin, habitatMask)',
                    types = c('s = double(1)', 'habitatMask = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)

        registerDistributions(
            list(
                dbinom_sparseLocalSCR = list(
                    BUGSdist = 'dbinom_sparseLocalSCR(detNums, detIndices, size, p0, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor, habitatGrid, indicator)',
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)
        
        registerDistributions(
            list(
                dbinomLocal_normal = list(
                    BUGSdist = 'dbinomLocal_normal(detNums, detIndices, size, p0, sigma, s, trapCoords, localTrapsIndices, localTrapsNum, resizeFactor, habitatGrid, indicator)',
                    types = c('value = double(1)', 'detIndices = double(1)', 'size = double(1)', 's = double(1)', 'trapCoords = double(2)', 'localTrapsIndices = double(2)', 'localTrapsNum = double(1)', 'habitatGrid = double(2)'),
                    discrete = TRUE,
                    mixedSizes = TRUE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)
        

        registerDistributions(
            list(
                dbinom_vector = list(
                    BUGSdist = 'dbinom_vector(size, prob)',
                    types = c('value = double(1)', 'size = double(1)', 'prob = double(1)'),
                    discrete = TRUE,
                    mixedSizes = FALSE,
                    pqAvail = FALSE
                )
            ),
            verbose = FALSE)

    })
}
