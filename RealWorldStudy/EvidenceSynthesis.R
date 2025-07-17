library(dplyr)
library(survival)
source("RealWorldStudy/DatabaseDetails.R")

likelihoodData <- readRDS(file.path(rootFolder, "likelihoodData.rds"))

approximations <- sort(names(likelihoodData))

tableData <- list()
# approximation = "hermite"
# outcome = "outcome_2"
for (approximation in approximations) {
  message("Synthesizing evidence using ", approximation, " approximation")
  for (outcome in names(likelihoodData[[approximation]])) {
    message("- For ", outcome)
    as <- likelihoodData[[approximation]][[outcome]]
    if (approximation %in% c("custom", "grid", "normal")) {
      as <- bind_rows(as)
    }
    estimate <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(as)
    data <- tibble(
      type = "Fixed-effects",
      approximation = approximation,
      outcome = outcome,
      rr = estimate$rr, 
      lb = estimate$lb, 
      ub = estimate$ub
    )
    tableData[[length(tableData) + 1]] <- data
    
    if (approximation == "normal") {
      ma <- meta::metagen(
        TE = as$logRr,
        seTE = as$seLogRr,
        studlab = as.character(seq_along(as$logRr))
      )
      rfx <- summary(ma)$random
      data <- tibble(
        type = "Traditional random-effects",
        approximation = approximation,
        outcome = outcome,
        rr = exp(rfx$TE), 
        lb = exp(rfx$lower),
        ub = exp(rfx$upper)
      )
      tableData[[length(tableData) + 1]] <- data
    }
    estimate <- EvidenceSynthesis::computeBayesianMetaAnalysis(as)
    data <- tibble(
      type = "Random-effects",
      approximation = approximation,
      outcome = outcome,
      rr = exp(estimate$mu),
      lb = exp(estimate$mu95Lb),
      ub = exp(estimate$mu95Ub)
    )
    tableData[[length(tableData) + 1]] <- data
  }
}
tableData <- bind_rows(tableData)
saveRDS(tableData, "RealWorldStudy/tableData.rds")














hoi <- 1080
nc <- 2
normalApprox <- list()
adaptiveGridApproxHoi <- list()
adaptiveGridApproxNc <- list()
source("RealWorldStudy/DatabaseDetails.R")
for (i in 1:nrow(databases)) {
  database <- databases[i, ]
  message(sprintf("Fetching likelihood functions for %s", database$databaseId))
  dbEstimates <- getResultsSummary(database$outputFolder)
  dbEstimates$databaseId <- database$databaseId
  normalApprox[[i]] <- dbEstimates
  fileRef <- getFileReference(database$outputFolder)
  model <- readRDS(file.path(database$outputFolder, fileRef$sccsModelFile[fileRef$outcomeId == hoi]))
  if (!is.null(model$logLikelihoodProfiles) && length(model$logLikelihoodProfiles) > 0) {
    adaptiveGridApproxHoi[[length(adaptiveGridApproxHoi) + 1]] <- model$logLikelihoodProfiles[[1]]
    names(adaptiveGridApproxHoi)[length(adaptiveGridApproxHoi)] <- database$databaseId
  } 
  model <- readRDS(file.path(database$outputFolder, fileRef$sccsModelFile[fileRef$outcomeId == nc]))
  if (!is.null(model$logLikelihoodProfiles) && length(model$logLikelihoodProfiles) > 0) {
    adaptiveGridApproxNc[[length(adaptiveGridApproxNc) + 1]] <- model$logLikelihoodProfiles[[1]]
    names(adaptiveGridApproxNc)[length(adaptiveGridApproxNc)] <- database$databaseId
  } 
  
}
normalApprox <- do.call(rbind, normalApprox)

estimateHoi <- EvidenceSynthesis::computeBayesianMetaAnalysis(adaptiveGridApproxHoi)
estimateNc <- EvidenceSynthesis::computeBayesianMetaAnalysis(adaptiveGridApproxNc)

EvidenceSynthesis::plotMetaAnalysisForest(
  data = adaptiveGridApproxHoi,
  labels = names(adaptiveGridApproxHoi),
  estimate = estimateHoi,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maHoiAdaptiveGridRfx.png"
)

EvidenceSynthesis::plotMetaAnalysisForest(
  data = adaptiveGridApproxNc,
  labels = names(adaptiveGridApproxNc),
  estimate = estimateNc,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maNcAdaptiveGridRfx.png"
)

estimateHoi <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(adaptiveGridApproxHoi)
estimateNc <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(adaptiveGridApproxNc)

EvidenceSynthesis::plotMetaAnalysisForest(
  data = adaptiveGridApproxHoi,
  labels = names(adaptiveGridApproxHoi),
  estimate = estimateHoi,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maHoiAdaptiveGridFfx.png"
)

EvidenceSynthesis::plotMetaAnalysisForest(
  data = adaptiveGridApproxNc,
  labels = names(adaptiveGridApproxNc),
  estimate = estimateNc,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maNcAdaptiveGridFfx.png"
)

normalApproxHoi <- normalApprox[normalApprox$outcomeId == hoi & normalApprox$databaseId %in% names(adaptiveGridApproxHoi), ]
normalApproxNc <- normalApprox[normalApprox$outcomeId == nc & normalApprox$databaseId %in% names(adaptiveGridApproxNc), ]
estimateHoi <- EvidenceSynthesis::computeBayesianMetaAnalysis(as.data.frame(normalApproxHoi))
estimateNc <- EvidenceSynthesis::computeBayesianMetaAnalysis(as.data.frame(normalApproxNc))

EvidenceSynthesis::plotMetaAnalysisForest(
  data = normalApproxHoi,
  labels = normalApproxHoi$databaseId,
  estimate = estimateHoi,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maHoiNormalRfx.png"
)

EvidenceSynthesis::plotMetaAnalysisForest(
  data = normalApproxNc,
  labels = normalApproxNc$databaseId,
  estimate = estimateNc,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maNcNormalRfx.png"
)

estimateHoi <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(as.data.frame(normalApproxHoi))
estimateNc <- EvidenceSynthesis::computeFixedEffectMetaAnalysis(as.data.frame(normalApproxNc))

EvidenceSynthesis::plotMetaAnalysisForest(
  data = normalApproxHoi,
  labels = normalApproxHoi$databaseId,
  estimate = estimateHoi,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maHoiNormalFfx.png"
)

EvidenceSynthesis::plotMetaAnalysisForest(
  data = normalApproxNc,
  labels = normalApproxNc$databaseId,
  estimate = estimateNc,
  xLabel = "Incidence rate ratio",
  fileName = "RealWorldStudy/maNcNormalFfx.png"
)

            