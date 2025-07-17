library(ggplot2)
library(dplyr)

# Plot real-world example likelihood data ----------------------------------------------------------
likelihoodData <- readRDS( file.path(rootFolder, "likelihoodData.rds"))

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
        } else if (approximation == "adaptiveGrid") {
          points <- tibble(
            x = a$point,
            y = a$value,
            type = "Adaptive Grid"
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
vizData$type <- factor(vizData$type, levels = c("Likelihood", "Adaptive Grid", "Custom", "Normal"))
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
  geom_rect(xmin = log(0.25) + 0.1, xmax = log(4) - 0.1, ymin = -6, ymax = -3.2, fill = rgb(0.1, 0.1, 0.1, 0.6), color = "white", alpha = 0.8, data = countLabels) +
  geom_text(x = log(0.25) + 0.2, y = -3.4, aes(label = label), hjust = "left", vjust = "top", size = 3, lineheight = 0.9, color = "white", data = countLabels) +
  scale_size_manual(values = c(1, 2, 2, 2, 2, 2, 2) * 0.5) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  scale_color_manual(values = c("#DDDDDD", "#FBC511", "#EB6622", "#11A08A")) +
  scale_y_continuous("Partial log-likelihood", limits = c(-6.2, 0.2), expand = c(0,0)) +
  scale_x_continuous("Incidence Rate Ratio", limits = c(log(0.1), log(10)), breaks = log(breaks), labels = breaks, expand = c(0,0)) +
  facet_grid(database~outcome) +
  theme(panel.grid.minor = element_blank(), 
        panel.background = element_blank(), 
        panel.grid.major = element_line(color = "#444444"), 
        legend.title = element_blank(), 
        legend.key = element_blank(),
        axis.ticks.x = element_blank(), 
        panel.border = element_rect(fill = NA), 
        strip.background = element_blank(),
        plot.background = ggplot2::element_rect(fill = "black", color = NA),
        legend.background = ggplot2::element_rect(fill = "black", color = NA),
        text = ggplot2::element_text(color = "white"),
        axis.text = ggplot2::element_text(color = "white"),
        strip.text = ggplot2::element_text(color = "white"),
        legend.position = "right")
plot
ggsave("RealWorldStudy/DarkLikelihoodCurves.png", plot, width = 8, height = 6)

     
