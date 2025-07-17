library(survival)
library(EvidenceSynthesis)

# row = grid[1,]
row = list(treatedFraction = 0.9, hazardRatio = 2, nSites = 5, maxN = 2000, nStrata = 5, tau = 0, seed = 0.898389684967697)
analysisType = "cohort method"
doRandomEffects = TRUE
iterations = 1000
set.seed(row$seed)
doSimulation <- function(row, iterations, simulationFolder, analysisType = "sccs", doRandomEffects = TRUE) {
  
  eval <- function(estimate, trueRr, trueTau, time, type) {
    if (is.null(estimate$seLogRr)) {
      estimate$seLogRr <- (log(estimate$ci95Lb) - log(estimate$ci95Ub)) / (2*qnorm(0.975))
    }
    if (is.na(estimate$logRr) | is.infinite(estimate$logRr) | is.na(estimate$seLogRr) | is.infinite(estimate$seLogRr)) {
      estimate$logRr <- 0
      estimate$seLogRr <- 999
      estimate$ci95Lb <- 0
      estimate$ci95Ub <- 999
    }
    if (is.null(estimate$tau)) {
      estimate$tau <- 0
      estimate$tauCi95Lb <- 0
      estimate$tauCi95Ub <- 999
    }
    if (is.null(estimate$seTau)) {
      estimate$seTau <- (estimate$tauCi95Ub - estimate$tauCi95Lb) / (2*qnorm(0.975))
    }
    estimate$seTau
    result <- data.frame(coverageMu = estimate$ci95Lb <= trueRr & estimate$ci95Ub >= trueRr,
                         errorMu = estimate$logRr - log(trueRr),
                         precisionMu = 1/(estimate$seLogRr ^ 2),
                         coverageTau = estimate$tauCi95Lb <= trueTau & estimate$tauCi95Ub >= trueTau,
                         errorTau = estimate$tau - trueTau,
                         precisionTau = 1/(estimate$seTau ^ 2))
    result$type <- type
    result$time <- time["user.self"] + time["sys.self"]
    return(result)
  }
  
  generateFileName <- function(i, row) {
    return(file.path(simulationFolder, 
                     paste0(gsub("\\.", "", paste(paste0(names(row), row), collapse = "_")), "_", i, ".rds")))
  }
  
  createApproximations <- function(populations, approximationType) {
    getSccsFormula <- function(population) {
      xColnames <- colnames(population)[grep("x[0-9]+", colnames(population))]
      formula <- as.formula(paste("y ~ a +", paste0(xColnames, collapse = " + "), " + strata(stratumId) + offset(log(time))"))
      return(formula)
    }
    # population = populations[[1]]
    # approximationType = "grid"
    fitModelInDatabase <- function(population, approximationType) {
      if (nrow(population) == 0) {
        return(NULL)
      }
      if (analysisType == "sccs") {
        cyclopsData <- Cyclops::createCyclopsData(getSccsFormula(population), data = population, modelType = "cpr")
        treatmentVariable <- "a"
      } else {
        cyclopsData <- Cyclops::createCyclopsData(Surv(time, y) ~ x + strata(stratumId), data = population, modelType = "cox")
        treatmentVariable <- "x"
      }
      cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
      approximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, treatmentVariable, approximation = approximationType)
      return(approximation)
    }
    data <- lapply(populations, fitModelInDatabase, approximationType = approximationType)
    if (approximationType %in% c("grid with gradients", "adaptive grid")) {
      data <- data[!sapply(data, is.null)]
    } else {
      data <- do.call("rbind", data)
    }
    return(data)
  }
  
  estimateUsingStandardMetaAnalysis <- function(approximations) {
    
    doMa <- function() {
      meta <- meta::metagen(TE = approximations$logRr, 
                            seTE = approximations$seLogRr, 
                            studlab = rep("", nrow(approximations)), 
                            byvar = NULL,
                            sm = "RR",
                            iterations = 3000)
      s <- summary(meta)
      rnd <- s$random
      return(data.frame(rr = exp(rnd$TE),
                        ci95Lb = exp(rnd$lower),
                        ci95Ub = exp(rnd$upper),
                        logRr = rnd$TE,
                        seLogRr = rnd$seTE,
                        tau = s$tau,
                        tauCi95Lb = s$lower.tau,
                        tauCi95Ub = s$upper.tau))
    }
    result <- tryCatch({
      doMa()
    }, error = function(e) {
      data.frame(rr = 1,
                 ci95lb = -Inf,
                 ci95ub = Inf,
                 logRr = 0,
                 seLogRr = Inf)
    }
    )
    return(result)
  }
  
  doIteration <- function(i, row, simulationFolder) {
    fileName <- generateFileName(i, row)
    if (file.exists(fileName)) {
      return(readRDS(fileName))
    } else {
      if (analysisType == "sccs") {
        settings <- EvidenceSynthesis::createSccsSimulationSettings(
          nSites = row$nSites,
          n = round(runif(row$nSites, 1000, row$maxN)),
          minBackgroundRate = 0.001,
          maxBackgroundRate = 0.01,
          atRiskTimeFraction = row$atRiskTimeFraction,
          rateRatio = row$rateRatio,
          timeCovariates = 5,
          randomEffectSd = row$tau
        )
        trueEffectSize <- row$rateRatio
      } else {
        settings <- EvidenceSynthesis::createSimulationSettings(
          nSites = row$nSites,
          n = round(runif(row$nSites, 1000, row$maxN)),
          minBackgroundHazard = 0.00001,
          maxBackgroundHazard = 0.0001,
          treatedFraction = row$treatedFraction,
          hazardRatio = row$hazardRatio,
          randomEffectSd = row$tau,
          nStrata = row$nStrata
        )
        trueEffectSize <- row$hazardRatio
      }
      populations <<- EvidenceSynthesis::simulatePopulations(settings)
      
      # Normal
      timeNormal <- system.time(
        approxs <- createApproximations(populations, "normal")
      )
      estimate <- tryCatch({
        EvidenceSynthesis::computeFixedEffectMetaAnalysis(approxs)
      }, error = function(e) {
        data.frame(rr = 1,
                   lb = -Inf,
                   ub = Inf,
                   logRr = 0,
                   seLogRr = Inf)
      })
      maFixedFxEstimate <- data.frame(rr = estimate$rr,
                                      ci95Lb = estimate$lb,
                                      ci95Ub = estimate$ub,
                                      logRr = estimate$logRr,
                                      seLogRr = estimate$seLogRr)
      if (doRandomEffects) {
        maRandomFxEstimate <- estimateUsingStandardMetaAnalysis(approxs)
        estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(approxs)
        normalRandomFxEstimate <- data.frame(rr = exp(estimate$mu),
                                             ci95Lb = exp(estimate$mu95Lb),
                                             ci95Ub = exp(estimate$mu95Ub),
                                             logRr = estimate$mu,
                                             tau = estimate$tau,
                                             tauCi95Lb = estimate$tau95Lb,
                                             tauCi95Ub = estimate$tau95Ub)
      }
      
      # Grid
      timeGrid <- system.time(
        approxs <- createApproximations(populations, "grid")
      )
      estimate <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(approxs)
      gridFixedFxEstimate <- data.frame(rr = estimate$rr,
                                        ci95Lb = estimate$lb,
                                        ci95Ub = estimate$ub,
                                        logRr = estimate$logRr,
                                        seLogRr = estimate$seLogRr)
      if (doRandomEffects) {
        estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(approxs)
        gridRandomFxEstimate <- data.frame(rr = exp(estimate$mu),
                                           ci95Lb = exp(estimate$mu95Lb),
                                           ci95Ub = exp(estimate$mu95Ub),
                                           logRr = estimate$mu,
                                           tau = estimate$tau,
                                           tauCi95Lb = estimate$tau95Lb,
                                           tauCi95Ub = estimate$tau95Ub)
      }
      
      # Custom parametric function
      timeCustom <- system.time(
        approxs <- createApproximations(populations, "custom")
      )
      estimate <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(approxs)
      customFixedFxEstimate <- data.frame(rr = estimate$rr,
                                          ci95Lb = estimate$lb,
                                          ci95Ub = estimate$ub,
                                          logRr = estimate$logRr,
                                          seLogRr = estimate$seLogRr)
      if (doRandomEffects) {
        estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(approxs)
        customRandomFxEstimate <- data.frame(rr = exp(estimate$mu),
                                             ci95Lb = exp(estimate$mu95Lb),
                                             ci95Ub = exp(estimate$mu95Ub),
                                             logRr = estimate$mu,
                                             tau = estimate$tau,
                                             tauCi95Lb = estimate$tau95Lb,
                                             tauCi95Ub = estimate$tau95Ub)
      }
      
      # Hermite
      timeHermite<- system.time(
        approxs <- createApproximations(populations, "grid with gradients")
      )
      estimate <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(approxs)
      hermiteFixedFxEstimate <- data.frame(rr = estimate$rr,
                                           ci95Lb = estimate$lb,
                                           ci95Ub = estimate$ub,
                                           logRr = estimate$logRr,
                                           seLogRr = estimate$seLogRr)
      if (doRandomEffects) {
        estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(approxs)
        hermiteRandomFxEstimate <- data.frame(rr = exp(estimate$mu),
                                              ci95Lb = exp(estimate$mu95Lb),
                                              ci95Ub = exp(estimate$mu95Ub),
                                              logRr = estimate$mu,
                                              tau = estimate$tau,
                                              tauCi95Lb = estimate$tau95Lb,
                                              tauCi95Ub = estimate$tau95Ub)
      }
      
      result <- rbind(
        eval(maFixedFxEstimate, trueEffectSize, settings$randomEffectSd, timeNormal, "Traditional fixed-effects"),
        eval(gridFixedFxEstimate, trueEffectSize, settings$randomEffectSd, timeGrid, "Grid fixed-effects"),
        eval(customFixedFxEstimate, trueEffectSize, settings$randomEffectSd, timeCustom, "Custom fixed-effects"),
        eval(hermiteFixedFxEstimate, trueEffectSize, settings$randomEffectSd, timeHermite, "Hermite interpolation fixed-effects")
      )
      if (doRandomEffects) {
        result <- rbind(
          result,
          eval(maRandomFxEstimate, trueEffectSize, settings$randomEffectSd, timeNormal, "Traditional random-effects"),
          eval(normalRandomFxEstimate, trueEffectSize, settings$randomEffectSd, timeNormal, "Normal random-effects"),
          eval(gridRandomFxEstimate, trueEffectSize, settings$randomEffectSd, timeGrid, "Grid random-effects"),
          eval(customRandomFxEstimate, trueEffectSize, settings$randomEffectSd, timeCustom, "Custom random-effects"),
          eval(hermiteRandomFxEstimate, trueEffectSize, settings$randomEffectSd, timeHermite, "Hermite interpolation random-effects")
        )
      }
      saveRDS(result, fileName)
      return(result)
    }
  }
  set.seed(row$seed)
  row$seed <- NULL
  fileName <- file.path(simulationFolder, 
                        paste0(gsub("\\.", "", paste(paste0(names(row), row), collapse = "_")), ".rds"))
  if (file.exists(fileName)) {
    return(readRDS(fileName))
  } else {  
    ParallelLogger::logTrace("Simulation using ", paste(paste(colnames(row), row, sep = " = "), collapse = ", "))
    results <- plyr::llply(1:iterations, doIteration, row = row, simulationFolder = simulationFolder)
    results <- do.call("rbind", results)
    rownames(row) <- NULL
    results <- cbind(results, row)
    saveRDS(results, fileName)
    
    toDelete <- fileName <- generateFileName(1:iterations, row)
    unlink(toDelete)
    
    ParallelLogger::logTrace("Done")
  }
  return(results)
}

runSimulationStudy <- function(simulationFolder, grid, analysisType = "sccs", doRandomEffects = TRUE, maxCores = 16) {
  dir.create(simulationFolder, recursive = TRUE)
  
  unlink(file.path(simulationFolder, "log.txt"))
  ParallelLogger::addDefaultFileLogger(file.path(simulationFolder, "log.txt"), name = "MY_FILE_LOGGER")
  ParallelLogger::addDefaultErrorReportLogger(file.path(simulationFolder, "errorReport.txt"), name = "MY_ERROR_LOGGER")
  on.exit(ParallelLogger::unregisterLogger("MY_FILE_LOGGER"))
  on.exit(ParallelLogger::unregisterLogger("MY_ERROR_LOGGER"), add = TRUE)
  iterations <- 1000
  
  set.seed(1)
  grid$seed <- runif(nrow(grid))
  
  cluster <- ParallelLogger::makeCluster(maxCores)
  on.exit(ParallelLogger::stopCluster(cluster), add = TRUE)
  ParallelLogger::clusterRequire(cluster, "survival")
  results <- ParallelLogger::clusterApply(
    cluster = cluster, 
    x = split(grid, 1:nrow(grid)), 
    fun = doSimulation, 
    iterations = iterations, 
    simulationFolder = simulationFolder,
    analysisType = analysisType,
    doRandomEffects = doRandomEffects
  )
  results <- do.call(rbind, results)
  return(results)
}

exportResultsForShiny <- function(results, fileName) {
  # Find Non-estimables
  results$nonEstimable <- (results$precisionMu == 1/(999^2))
  writeLines(paste(sum(results$nonEstimable) , "of", nrow(results), " where unestimable"))
  results$coverageMu[results$nonEstimable] <- NA
  results$errorMu[results$nonEstimable] <- NA
  results$precisionMu[results$nonEstimable] <- NA
  results$coverageTau[results$nonEstimable] <- NA
  results$errorTau[results$nonEstimable] <- NA
  results$precisionTau[results$nonEstimable] <- NA
  
  aggregateMetric <- function(metric) {
    form <- paste(metric, "~ type + atRiskTimeFraction + rateRatio + nSites + maxN + tau")
    agg <- aggregate(as.formula(form), data = results, mean)
    colnames(agg)[colnames(agg) == metric] <- "value"
    agg$metric <- metric
    return(agg)
  }
  results$precisionMu <- log(results$precisionMu)
  results$biasMu <- results$errorMu
  results$mseMu <- results$errorMu ^ 2
  results$precisionTau <- log(results$precisionTau)
  results$biasTau <- results$errorTau
  results$mseTau <- results$errorTau ^ 2
  metrics <- c("coverageMu", "biasMu", "mseMu", "precisionMu", "coverageTau", "biasTau", "mseTau", "precisionTau", "nonEstimable")
  
  outputTable <- lapply(metrics, aggregateMetric)
  outputTable <- do.call(rbind, outputTable)
  idx <- outputTable$metric %in% c("precisionMu", "precisionTau")
  outputTable$value[idx] <- exp(outputTable$value[idx])
  
  outputTable$metric <- gsub("Mu", "", outputTable$metric)
  
  # Pretty names for metrics
  outputTable$metric[outputTable$metric == "coverage"] <- "Coverage"
  outputTable$metric[outputTable$metric == "bias"] <- "Bias"
  outputTable$metric[outputTable$metric == "mse"] <- "MSE"
  outputTable$metric[outputTable$metric == "precision"] <- "Precision"
  outputTable$metric[outputTable$metric == "nonEstimable"] <- "Non-Estimable"
  outputTable$metric[outputTable$metric == "coverageTau"] <- "CoverageTau"
  outputTable$metric[outputTable$metric == "biasTau"] <- "BiasTau"
  outputTable$metric[outputTable$metric == "mseTau"] <- "MSETau"
  outputTable$metric[outputTable$metric == "precisionTau"] <- "PrecisionTau"
  
  saveRDS(outputTable, file.path("inst/shinyApps/ResultsExplorer/data", fileName))
}
