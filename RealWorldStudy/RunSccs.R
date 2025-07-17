library(SelfControlledCaseSeries)

# Create analysis settings -----------------------------------------------------
getDbSccsDataArgs <- createGetDbSccsDataArgs(
  studyStartDates = "20090101",
  studyEndDates = "20101231"
)
createStudyPopulationArgs <- createCreateStudyPopulationArgs(
  naivePeriod = 365,
  firstOutcomeOnly = FALSE
)
covarExposureOfInt <- createEraCovariateSettings(
  label = "Exposure of interest",
  includeEraIds = "exposureId",
  start = 1,
  startAnchor = "era start",
  end = 28,
  endAnchor = "era start",
  profileLikelihood = TRUE,
  exposureOfInterest = TRUE
)
covarPreExposure <- createEraCovariateSettings(
  label = "Pre-exposure",
  includeEraIds = "exposureId",
  start = -30,
  end = -1,
  endAnchor = "era start"
)
seasonalityCovariateSettings <- createSeasonalityCovariateSettings(
  seasonKnots = 5
)
createSccsIntervalDataArgs <- createCreateSccsIntervalDataArgs(
  eraCovariateSettings = list(covarExposureOfInt, covarPreExposure),
  seasonalityCovariateSettings = seasonalityCovariateSettings
)
fitSccsModelArgs <- createFitSccsModelArgs()
sccsAnalysisList <- list(
  createSccsAnalysis(
    analysisId = 1,
    description = "SCCS, tar 1-28 days",
    getDbSccsDataArgs = getDbSccsDataArgs,
    createStudyPopulationArgs = createStudyPopulationArgs,
    createIntervalDataArgs = createSccsIntervalDataArgs,
    fitSccsModelArgs = fitSccsModelArgs
  )
)
exposures <- list(createExposure(1))
exposuresOutcomeList <- list(
  createExposuresOutcome(outcomeId = 2, exposures = exposures),
  createExposuresOutcome(outcomeId = 1080, exposures = exposures)
)

# Run analyses -----------------------------------------------------------------
source("RealWorldStudy/DatabaseDetails.R")
dir.create(rootFolder)
for (i in 1:nrow(databases)) {
  database <- databases[i, ]
  message(sprintf("Running SCCS for %s", database$name))
  sccsMultiThreadingSettings = createDefaultSccsMultiThreadingSettings(parallel::detectCores())
  dir.create(database$folder, recursive = TRUE)
  runSccsAnalyses(
    connectionDetails = connectionDetails,
    cdmDatabaseSchema = database$cdmDatabaseSchema,
    exposureDatabaseSchema = database$cohortDatabaseSchema,
    exposureTable = database$cohortTable,
    outcomeDatabaseSchema = database$cohortDatabaseSchema,
    outcomeTable = database$cohortTable,
    exposuresOutcomeList = exposuresOutcomeList,
    sccsAnalysisList = sccsAnalysisList,
    outputFolder = database$folder,
    sccsMultiThreadingSettings = sccsMultiThreadingSettings
  )
}

# Collect estimates ------------------------------------------------------------
estimates <- list()
source("RealWorldStudy/DatabaseDetails.R")
for (i in 1:nrow(databases)) {
  database <- databases[i, ]
  message(sprintf("Running SCCS for %s", database$name))
  dbEstimates <- getResultsSummary(database$folder)
  dbEstimates$databaseId <- database$name
  estimates[[i]] <- dbEstimates
}
estimates <- do.call(rbind, estimates)

# Create likelihood approximations -----------------------------------------------------------------
source("RealWorldStudy/DatabaseDetails.R")
likelihoodData <- list()
for (i in 1:nrow(databases)) {
  database <- databases[i, ]
  message(sprintf("Building approximations for %s", database$name))
  fileRef <- SelfControlledCaseSeries::getFileReference(database$folder)
  for (j in seq_len(nrow(fileRef))) {
    sccsIntervalData <- SelfControlledCaseSeries::loadSccsIntervalData(file.path(database$folder, fileRef$sccsIntervalDataFile[j]))
    if (sccsIntervalData$outcomes %>% count() %>% pull() != 0 && 
        sccsIntervalData$covariates %>% filter(covariateId == 1000) %>% count() %>% pull() != 0) {
      
      # Person-level data for pooling:
      population <- EvidenceSynthesis::prepareSccsIntervalData(sccsIntervalData, 1000)
      # Make x covariates unique to each database (no common seasonality)
      xNames <- names(population)[grep("^x", names(population))] 
      names(population)[grep("^x", names(population))] <- paste(xNames, i, sep = "")
      xNamesZeroes <- paste(xNames, rep(seq_len(nrow(databases))[-i], each = length(xNames)), sep = "")
      zeroes <- array(0, dim = c(nrow(population), length(xNamesZeroes)))
      colnames(zeroes) = xNamesZeroes
      population <- cbind(population, zeroes)
      likelihoodData[["population"]][[paste("outcome", fileRef$outcomeId[j], sep = "_")]][[database$name]] <- population
      
      cyclopsData <- Cyclops::convertToCyclopsData(sccsIntervalData$outcomes, sccsIntervalData$covariates, modelType = "cpr", addIntercept = FALSE)
      cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData, control = Cyclops::createControl(threads = 8))
      
      # Normal:
      approximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, parameter = "1000", approximation = "normal")
      likelihoodData[["normal"]][[paste("outcome", fileRef$outcomeId[j], sep = "_")]][[database$name]] <- approximation
      
      # Grid:
      approximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, parameter = "1000", approximation = "grid")
      likelihoodData[["grid"]][[paste("outcome", fileRef$outcomeId[j], sep = "_")]][[database$name]] <- approximation

      # Grid with gradients:
      approximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, parameter = "1000", approximation = "grid with gradients")
      likelihoodData[["hermite"]][[paste("outcome", fileRef$outcomeId[j], sep = "_")]][[database$name]] <- approximation
      
      # Custom:
      approximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, parameter = "1000", approximation = "custom")
      likelihoodData[["custom"]][[paste("outcome", fileRef$outcomeId[j], sep = "_")]][[database$name]] <- approximation
    }
  }
}
saveRDS(likelihoodData, file.path(rootFolder, "likelihoodData.rds"))
