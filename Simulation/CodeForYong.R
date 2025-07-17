# Example code demonstrating evidence synthesis across multiple sites using likelihood profiles:

# Install a special version of EvidenceSynthesis that contains code to simulate SCCS data:
remotes::install_github("ohdsi/EvidenceSynthesis", ref = "pooled_sccs")

# Simulate SCCS data:
settings <- EvidenceSynthesis::createSccsSimulationSettings(nSites = 10,
                                                            n = 2000,
                                                            minBackgroundRate = 0.001,
                                                            maxBackgroundRate = 0.01,
                                                            atRiskTimeFraction = 0.05,
                                                            rateRatio = 2,
                                                            timeCovariates = 5,
                                                            randomEffectSd = 0.5)
populations <- EvidenceSynthesis::simulatePopulations(settings)

# Fit models for each population, and create likelihood profile:
library(survival)
profiles <- list()

i <- 1
for (i in seq_along(populations)) {
  population <- populations[[i]]
  xColnames <- colnames(population)[grep("x[0-9]+", colnames(population))]
  formula <- as.formula(paste("y ~ a +", paste0(xColnames, collapse = " + "), " + strata(stratumId) + offset(log(time))"))
  cyclopsData <- Cyclops::createCyclopsData(formula, data = population, modelType = "cpr")
  cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
  
  # Example: Using Cyclops' getCyclopsProfileLogLikelihood() function, which uses the adaptive grid approach:
  profiles[[i]] <- Cyclops::getCyclopsProfileLogLikelihood(cyclopsFit, "a", bounds = c(log(0.1), log(10)))
  # The source code of this function is here: https://github.com/OHDSI/Cyclops/blob/main/R/ModelFit.R#L892
  # Essentially, this function calls Cyclops:::fixedGridProfileLogLikelihood() (source code: https://github.com/OHDSI/Cyclops/blob/main/R/ModelFit.R#L987)
  # to evaluate the log-likelihood at a set of points. The idea of the adaptive grid is to get a
  # good approximation with as few likelihood evaluations.
}

# Perform random-effects meta-analysis:
maEstimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(profiles)
sprintf("Meta-analytic estimate (mu): %0.2f (%0.2f-%0.2f). True value: %0.2f", exp(maEstimate$mu), exp(maEstimate$mu95Lb), exp(maEstimate$mu95Ub), settings$rateRatio)

sprintf("Estimate of tau: %0.2f (%0.2f-%0.2f). True value: %0.2f", maEstimate$tau, maEstimate$tau95Lb, maEstimate$tau95Ub, settings$randomEffectSd)
