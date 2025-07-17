# Run study ----------------------------------------------------------------------
source("Simulation/SimulationFunctions.R")

simulationFolder <- "e:/EsSimulations/randomFxCm"

grid <- expand.grid(treatedFraction = c(0.1, 0.5, 0.9),
                    hazardRatio = c(1, 2, 4),
                    nSites = c(5, 10, 20),
                    maxN = c(2000, 4000),
                    nStrata = c(5),
                    tau = c(0, 0.25, 0.5, 1))

results <- runSimulationStudy(
  simulationFolder = simulationFolder,
  grid = grid,
  analysisType = "Cohort method",
  doRandomEffects = TRUE, 
  maxCores = 10
)
saveRDS(results, "Simulation/cmRandomFxSimulationResults.rds")
