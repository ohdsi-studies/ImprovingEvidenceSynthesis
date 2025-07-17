# Run study ------------------------
source("Simulation/SimulationFunctions.R")

simulationFolder <- "e:/EsSimulations/fixedFxSccs"

grid <- expand.grid(atRiskTimeFraction = c(0.0625, 0.125, 0.25, 0.5),
                    rateRatio = c(1, 2, 4),
                    nSites = c(5, 10, 20),
                    maxN = c(2000, 3000, 4000, 5000),
                    tau = 0)

results <- runSimulationStudy(
  simulationFolder = simulationFolder,
  grid = grid,
  analysisType = "sccs",
  doRandomEffects = FALSE, 
  maxCores = 12
)
saveRDS(results, "Simulation/fixedFxSimulationResults.rds")
exportResultsForShiny(results, "fixedFx.rds")

library(dplyr)
results |>
  group_by(type, maxN) |>
  summarise(medianTime = median(time))
# type                                 maxN medianTime
# <chr>                               <dbl>      <dbl>
# 1 Custom fixed-effects                 2000      0.800
# 2 Custom fixed-effects                 3000      0.750
# 3 Custom fixed-effects                 4000      0.790
# 4 Custom fixed-effects                 5000      0.840
# 5 Grid fixed-effects                   2000      3.37 
# 6 Grid fixed-effects                   3000      3.35 
# 7 Grid fixed-effects                   4000      3.63 
# 8 Grid fixed-effects                   5000      3.94 
# 9 Hermite interpolation fixed-effects  2000      0.280
# 10 Hermite interpolation fixed-effects  3000      0.230
# 11 Hermite interpolation fixed-effects  4000      0.240
# 12 Hermite interpolation fixed-effects  5000      0.250
# 13 Traditional fixed-effects            2000      0.200
# 14 Traditional fixed-effects            3000      0.210
# 15 Traditional fixed-effects            4000      0.230
# 16 Traditional fixed-effects            5000      0.240