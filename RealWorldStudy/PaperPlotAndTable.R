library(survival)
library(ggplot2)
library(dplyr)
source("RealWorldStudy/DatabaseDetails.R")

# Plot real-world example likelihood data ----------------------------------------------------------
# likelihoodData <- readRDS( file.path(rootFolder, "likelihoodData.rds"))
likelihoodData <- readRDS("/Users/schuemie/Data/likelihoodData.rds")

approximations <- names(likelihoodData)
vizData <- list()
countData <- list()
x <- seq(log(0.1), log(10), length.out = 200)
for (approximation in approximations) {
  for (outcome in names(likelihoodData[[approximation]])) {
    if (outcome == "outcome_2") {
      outcomeName <- "GBS"
    } else {
      outcomeName <- "Hypothermia"
    }
    for (database in names(likelihoodData[[approximation]][[outcome]])) {
      a <- likelihoodData[[approximation]][[outcome]][[database]]
      if (approximation == "population") {
        exposedCaseIds <- a %>%
          filter(a == 1) %>%
          distinct(stratumId) %>%
          pull()
        exposedCaseData <- a %>%
          filter(stratumId %in% exposedCaseIds)
        countData[[length(countData) + 1]] <- tibble(
          caseCount = a %>%
            distinct(stratumId) %>%
            count() %>%
            pull(),
          exposedCaseCount = length(exposedCaseIds),
          fractionExposed = exposedCaseData %>%
            filter(a == 1) %>%
            summarize(sum(time)) %>%
            pull() / exposedCaseData %>%
            summarize(sum(time)) %>%
            pull() ,
          database = database,
          outcome = outcomeName
        )
      } else {
        if (approximation == "normal") {
          points <- tibble(
            x = x,
            y = dnorm(x, mean = a$logRr, sd = a$seLogRr, log = TRUE),
            type = "Normal"
          )
        } else if (approximation == "grid") {
          points <- tibble(
            x = as.numeric(names(a)),
            y = as.vector(a),
            type = "Likelihood"
          )
        } else if (approximation == "hermite") {
          points <- tibble(
            x = x,
            y = EvidenceSynthesis::hermiteInterpolation(x, a),
            type = "Hermite interpolation"
          )
        } else if (approximation == "custom") {
          points <- tibble(
            x = x,
            y = EvidenceSynthesis::customFunction(x, a$mu, a$sigma, a$gamma),
            type = "Custom"
          )
        }
        points <- points %>%
          mutate(outcome = outcomeName,
                 database = database,
                 y = y - max(y))
        vizData[[length(vizData) + 1]] <- points
      }
    }
  }
}
vizData <- bind_rows(vizData)
vizData$type <- factor(vizData$type, levels = c("Likelihood","Custom", "Hermite interpolation", "Normal"))
countData <- bind_rows(countData)
countLabels <- countData %>%
  mutate(x = 1,
         y = 1,
         label = sprintf("Cases: %s\nCases with exposure: %s\nFraction at risk: %.2f", 
                         formatC(caseCount, big.mark = ","), 
                         formatC(exposedCaseCount, big.mark = ","), 
                         fractionExposed))
breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)
plot <- ggplot(vizData, aes(x = x, y = y)) +
  geom_line(aes(group = type, color = type, linetype = type, size = type), alpha = 0.9) +
  geom_rect(xmin = log(0.25) + 0.1, xmax = log(4) - 0.1, ymin = -6, ymax = -3.2, fill = rgb(1, 1, 1, 0.8), color = "black", alpha = 0.8, data = countLabels) +
  geom_text(x = log(0.25) + 0.2, y = -3.4, aes(label = label), hjust = "left", vjust = "top", size = 3, lineheight = 0.9, data = countLabels) +
  # scale_size_manual(values = c(1, 2,4, 2,4, 2,4) * 0.5) +
  scale_size_manual(values = c(0.5, 1.25, 1.25, 1.25)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  scale_color_manual(values = c("#000000", "#FBC511", "#EB6622", "#11A08A")) +
  scale_y_continuous("Log-likelihood", limits = c(-6.2, 0.2), expand = c(0,0)) +
  scale_x_continuous("Incidence Rate Ratio", limits = c(log(0.1), log(10)), breaks = log(breaks), labels = breaks, expand = c(0,0)) +
  facet_grid(database~outcome) +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(color = "#DDDDDD"), 
        legend.title = element_blank(), 
        legend.key = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.border = element_rect(fill = NA), 
        strip.background = element_blank(),
        legend.position = "top")
plot
ggsave("RealWorldStudy/LikelihoodCurves.png", plot, width = 6, height = 6.6)



# Create real-world example table ------------------------------------------------------------------
library(dplyr)

tableData <- readRDS("RealWorldStudy/tableData.rds")

table <- tableData |>
  mutate(irr = sprintf("%0.2f (%0.2f-%0.2f)", rr, lb, ub),
         name = if_else(outcome == "outcome_2", "GBS", "Hypothermia"),
         approximation = case_when(
           approximation == "population" ~  "Pooled",
           approximation == "normal" ~  "Normal",
           approximation == "grid" ~  "Grid",
           approximation == "hermite" ~  "Hermite interpolation",
           approximation == "custom" ~  "Custom"
         )) |>
  select(name, type, approximation, irr) |>
  arrange(name, type, approximation)


reshapedTable <- table |>
  tidyr::pivot_wider(id_cols = c("type", "approximation"), values_from = irr, names_from = name) |>
  mutate(label = paste(approximation, tolower(type))) |>
  arrange(type) |>
  select(-approximation, -type) |>
  relocate(label, 1)

readr::write_csv(reshapedTable, "RealWorldStudy/realWorldTable.csv")

# Plot MA estimates ----------------
breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 6, 8, 10)
subset <- tableData[grepl("random", tableData$type, ignore.case = TRUE), ]
plot <- ggplot(subset, aes(x = rr, y = approximation , xmin = lb, xmax = ub)) +
  geom_vline(xintercept = breaks, colour = "#AAAAAA", lty = 1, size = 0.2) +
  geom_vline(xintercept = 1, size = 0.5) +
  geom_errorbarh(height = 0.15) +
  geom_point(size = 3, shape = 23, fill = "#FFFFFF") +
  scale_x_continuous("Hazard Ratio", trans = "log10", breaks = breaks, labels = breaks) +
  coord_cartesian(xlim = c(0.2, 5)) +
  facet_grid(~outcome) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        legend.position = "none",
        panel.border = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks = element_blank(),
        plot.margin = grid::unit(c(0,0,0.1,0), "lines"))
plot


