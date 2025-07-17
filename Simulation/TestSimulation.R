# Some code to see if the SCCS simulation works
source("Simulation/SimulateSccs.R")
library(Cyclops) 
library(survival)
# Single site simulation
settings <- createSccsSimulationSettings(nSites = 1,
                                         n = 10000,
                                         atRiskTimeFraction = 0.1,
                                         timePartitions = 10,
                                         timeCovariates = 5,
                                         timeEffectSize = 1,
                                         minBackgroundRate = 0.1,
                                         maxBackgroundRate = 1,
                                         rateRatio = 2,
                                         randomEffectSd = 0)
population <- simulateSccsPopulations(settings)[[1]]
cyclopsData <- createCyclopsData(getFormula(population), data = population, modelType = "cpr")
fit <- fitCyclopsModel(cyclopsData)

exp(coef(fit)[1])
exp(confint(fit, parm = "a"))[2:3]

sum(population$y)
settings$n * mean(c(settings$minBackgroundRate, settings$maxBackgroundRate))

# Multiple site simulation
settings <- createSccsSimulationSettings(nSites = 10,
                                         n = 10000,
                                         atRiskTimeFraction = 0.1,
                                         timePartitions = 10,
                                         timeCovariates = 5,
                                         timeEffectSize = 1,
                                         minBackgroundRate = 0.1,
                                         maxBackgroundRate = 1,
                                         rateRatio = 2,
                                         randomEffectSd = 0)
populations <- simulateSccsPopulations(settings)
profiles <- list()
for (i in seq_along(populations)) {
  population <- populations[[i]]
  xColnames <- colnames(population)[grep("x[0-9]+", colnames(population))]
  formula <- as.formula(paste("y ~ a +", paste0(xColnames, collapse = " + "), " + strata(stratumId) + offset(log(time))"))
  cyclopsData <- createCyclopsData(formula, data = population, modelType = "cpr")
  fit <- fitCyclopsModel(cyclopsData)
  profile <- getCyclopsProfileLogLikelihood(fit, "a", bounds = c(log(0.1), log(10)))
  profiles[[i]] <- profile
}
results <- EvidenceSynthesis::computeBayesianMetaAnalysis(profiles)
sprintf("%0.2f (%0.2f-%0.2f)", exp(results$mu), exp(results$mu95Lb), exp(results$mu95Ub))
# [1] "1.96 (1.90-2.02)"
