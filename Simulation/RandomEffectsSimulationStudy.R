# Run study ----------------------------------------------------------------------
source("Simulation/SimulationFunctions.R")

simulationFolder <- "e:/EsSimulations/randomFxSccs"
grid <- expand.grid(atRiskTimeFraction = c(0.0625, 0.125, 0.25, 0.5),
                    rateRatio = c(1, 2, 4),
                    nSites = c(5, 10, 20),
                    maxN = c(2000, 4000),
                    tau = c(0, 0.25, 0.5, 1))

results <- runSimulationStudy(
  simulationFolder = simulationFolder,
  grid = grid,
  analysisType = "sccs",
  doRandomEffects = TRUE, 
  maxCores = 10
)
saveRDS(results, "Simulation/randomFxSimulationResults.rds")
exportResultsForShiny(results, "randomFx.rds")

library(dplyr)
results |>
  group_by(type, maxN) |>
  summarise(medianTime = median(time),
            p90Time = quantile(time, 0.9))

# type                                  maxN medianTime p90Time
# 1 Custom fixed-effects                  2000      0.730   4.08 
# 2 Custom fixed-effects                  4000      0.730   2.63 
# 3 Custom random-effects                 2000      0.730   4.08 
# 4 Custom random-effects                 4000      0.730   2.63 
# 5 Grid fixed-effects                    2000      3.11   10.7  
# 6 Grid fixed-effects                    4000      3.38    8.53 
# 7 Grid random-effects                   2000      3.11   10.7  
# 8 Grid random-effects                   4000      3.38    8.53 
# 9 Hermite interpolation fixed-effects   2000      0.25    0.710
# 10 Hermite interpolation fixed-effects   4000      0.220   0.520
# 11 Hermite interpolation random-effects  2000      0.25    0.710
# 12 Hermite interpolation random-effects  4000      0.220   0.520
# 13 Normal random-effects                 2000      0.190   0.450
# 14 Normal random-effects                 4000      0.200   0.480
# 15 Traditional fixed-effects             2000      0.190   0.450
# 16 Traditional fixed-effects             4000      0.200   0.480
# 17 Traditional random-effects            2000      0.190   0.450
# 18 Traditional random-effects            4000      0.200   0.480

