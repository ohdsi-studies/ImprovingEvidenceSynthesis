getReferenceValues <- function(metrics) {
  ref <- data.frame()
  if ("Bias" %in% metrics) {
    ref <- rbind(ref, data.frame(value = 0,
                                 metric = "Bias")) 
  }
  if ("Coverage" %in% metrics) {
    ref <- rbind(ref, data.frame(value = 0.95,
                                 metric = "Coverage")) 
  }
  if ("MSE" %in% metrics) {
    ref <- rbind(ref, data.frame(value = 0,
                                 metric = "MSE")) 
  }
  if ("Non-Estimable" %in% metrics) {
    ref <- rbind(ref, data.frame(value = 0,
                                 metric = "Non-Estimable")) 
  }
  return(ref)
}

plotViolin <- function(group = "Fixed-effects", 
                       types = c("Traditional", "Custom", "Normal", "Hermite interpolation", "Grid"),
                       metrics = c("Bias", "Coverage", "MSE", "Non-Estimable", "Precision"),
                       xParameter = "atRiskTimeFraction",
                       xLabel = "At-risk time fraction",
                       theme = "dark") {
  if (group == "Fixed-effects") {
    data <- readRDS("inst/shinyApps/ResultsExplorer/data/fixedFx.rds")
    data$type[data$type == "Traditional fixed-effects"] <- "Normal fixed-effects"
  } else if (group == "Random-effects") {
    data <- readRDS("inst/shinyApps/ResultsExplorer/data/randomFx.rds")
    data <- data[!grepl("Tau", data$metric), ]
  } else if (group == "Cohort method random-effects") {
    data <- readRDS("inst/shinyApps/ResultsExplorer/data/cmRandomFx.rds")
    data <- data[!grepl("Tau", data$metric), ]
    group <- "Random-effects"
  } else {
    stop("Unknown group: ", group)
  }
  data <- data[grepl(group, data$type, ignore.case = TRUE), ]
  data <- data[grepl(paste(types, collapse = "|"), data$type), ]
  data$type <- gsub(paste0(" ", group), "", data$type, ignore.case = TRUE)
  
  # Pivot data:
  pivotData <- function(simParam, data) {
    data$parameterValue <- data[, simParam]
    data$simParam <- simParam
    data[simParam] <- NULL
    return(data)
  }
  vizData <- pivotData(xParameter, data)
  ref <- getReferenceValues(vizData$metric)
  vizData <- vizData[vizData$metric %in% metrics,]
  ref <- ref[ref$metric %in% metrics,]
  
  if (theme == "dark") {
    bgColor <- "black"
    textColor <- "white"
    gridColor <- "#444444"
    refColor <- "#AAAAAA"
    alpha <- 0.6
    plot <- ggplot2::ggplot(vizData, ggplot2::aes(x = factor(parameterValue), y = value, fill = type, color = type))
  } else {
    bgColor <- "white"
    textColor <- "black"
    gridColor <- "#DDDDDD"
    refColor <- "black"
    alpha <- 0.4
    plot <- ggplot2::ggplot(vizData, ggplot2::aes(x = factor(parameterValue), y = value, fill = type))
  }
  
  plot <- plot + ggplot2::geom_violin(position = ggplot2::position_dodge(0.9), scale = "width", alpha = alpha) 
  if (nrow(ref) > 0) {
    plot <- plot + ggplot2::geom_hline(ggplot2::aes(yintercept = value), data = ref, linetype = "dashed", color = refColor)
  }
  plot <- plot +
    ggplot2::scale_x_discrete(xLabel) +
    ggplot2::scale_fill_manual(values = c("#FBC511", "#69AED5", "#EB6622", "#11A08A", "#336B91")) +
    ggplot2::scale_color_manual(values = c("#FBC511", "#69AED5", "#EB6622", "#11A08A", "#336B91")) +
    ggplot2::facet_grid(metric~., scales = "free", switch = "both") +
    ggplot2::theme(legend.position = "top",
                   legend.title = ggplot2::element_blank(),
                   legend.key = ggplot2::element_blank(),
                   legend.background = ggplot2::element_blank(),
                   text = ggplot2::element_text(color = textColor),
                   axis.text = ggplot2::element_text(color = textColor),
                   strip.text = ggplot2::element_text(color = textColor),
                   panel.grid.major.y = ggplot2::element_line(color = gridColor), 
                   panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_rect(fill = bgColor, color = NA),
                   axis.title.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_rect(fill = NA), 
                   strip.placement = "outside",
                   strip.background = ggplot2::element_blank())
  #plot
  return(plot)
}
# For presentation
plot <- plotViolin(group = "Fixed-effects", 
                   theme = "dark",
                   metrics = c("Bias", "Coverage"))
ggplot2::ggsave("Simulation/FixedFxDark.png", plot, width = 6, height = 3, dpi = 300)
                   

plot <- plotViolin(group = "Random-effects", 
                   theme = "dark",
                   metrics = c("Bias", "Coverage"))
ggplot2::ggsave("Simulation/RandomFxDark.png", plot, width = 6, height = 3, dpi = 300)

plot <- plotViolin(group = "Random-effects", 
                   theme = "light",
                   metrics = c("Bias", "Coverage"))
ggplot2::ggsave("Simulation/RandomFxLight.png", plot, width = 6, height = 3, dpi = 300)

plot <- plotViolin(group = "Random-effects", 
                   theme = "dark",
                   xParameter = "nSites",
                   xLabel = "Number of sites",
                   metrics = c("Bias", "Coverage"))
ggplot2::ggsave("Simulation/RandomFxDarkNsites.png", plot, width = 6, height = 3, dpi = 300)

# For paper
plot <- plotViolin(group = "Fixed-effects", 
                   theme = "light")
ggplot2::ggsave("Simulation/FixedFxLight.png", plot, width = 6, height = 6, dpi = 300)

plot <- plotViolin(group = "Random-effects", 
                   theme = "light",
                   xParameter = "nSites",
                   xLabel = "Number of sites")
ggplot2::ggsave("Simulation/RandomFxLight.png", plot, width = 6, height = 6, dpi = 300)

plot <- plotViolin(group = "Cohort method random-effects", 
                   theme = "light",
                   xParameter = "nSites",
                   xLabel = "Number of sites")
ggplot2::ggsave("Simulation/CmRandomFxLight.png", plot, width = 6, height = 6, dpi = 300)


# Plot example of Hermite interpolation, grid with gradients ---------------------------------------
library(survival)
library(dplyr)
library(ggplot2)
set.seed(5)
settings <- EvidenceSynthesis::createSccsSimulationSettings(nSites = 1,
                                                            n = 2000,
                                                            minBackgroundRate = 0.001,
                                                            maxBackgroundRate = 0.01,
                                                            atRiskTimeFraction = 0.05,
                                                            rateRatio = 2,
                                                            timeCovariates = 5,
                                                            randomEffectSd = 0)
population <- EvidenceSynthesis::simulatePopulations(settings)[[1]]
cyclopsData <- Cyclops::createCyclopsData(y ~ a + x1 + x2 + x3 + x4 + x5 + strata(stratumId) + offset(log(time)), data = population, modelType = "cpr")
cyclopsFit <- Cyclops::fitCyclopsModel(cyclopsData)
approximation <- EvidenceSynthesis::approximateLikelihood(cyclopsFit, "a", approximation = "grid with gradients")


vizData <- EvidenceSynthesis:::getLikelihoodCoordinates(approximation, c(0.1, 10))
vizData$trueLl <- EvidenceSynthesis:::getLikelihoodProfile(cyclopsFit, parameter = "a", vizData$x)
normalizationFactor <- max(vizData$trueLl)
vizData <- vizData |>
  mutate(trueLl = trueLl - normalizationFactor)
approximation <- approximation |>
  mutate(value = value - normalizationFactor) 

# delta <- 0.1
# gradientLineData <- approximation |>
#   transmute(xMin = point - delta,
#             xMax = point + delta,
#             yMin = value - derivative * delta,
#             yMax = value + derivative * delta)

segmentLength <- 0.4
gradientLineData <- approximation |>
  transmute(
    # Calculate delta based on desired segment length and derivative
    delta = segmentLength / (2 * sqrt(1 + derivative^2)),
    xMin = point - delta,
    xMax = point + delta,
    yMin = value - derivative * delta,
    yMax = value + derivative * delta
  )

breaks <- c(0.1, 0.25, 0.5, 1, 2, 4, 8)
xLimits <- c(min(c(gradientLineData$xMin, gradientLineData$xMax)), max(c(gradientLineData$xMin, gradientLineData$xMax)))
yLimits <- c(min(c(gradientLineData$yMin, gradientLineData$yMax)), max(c(gradientLineData$yMin, gradientLineData$yMax)))

# For presentation: separate plots -----------------------------------------------------------------
ggplot(vizData, aes(x = x, y = trueLl)) +
  geom_line(alpha = 0.5, size = 1) +
  scale_x_continuous("Incidence Rate Ratio", breaks = log(breaks), labels = breaks, limits = xLimits) +
  scale_y_continuous("Log likelihood", limits = yLimits) +
  theme(legend.position = "top",
                 legend.title = element_blank(),
                 legend.key = element_blank(),
                 legend.background = element_blank(),
                 panel.grid.major = element_line(color = "#DDDDDD"), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 panel.border = element_rect(fill = NA), 
                 strip.placement = "outside",
                 strip.background = element_blank())
ggsave("Simulation/gridWithGradients1.png", width = 6, height = 3, dpi = 300)

ggplot(vizData, aes(x = x, y = trueLl)) +
  geom_line(alpha = 0.5, size = 1) +
  geom_point(aes(x = point, y = value), shape = 16, size = 3, color = "#EB6622", data = approximation)+
  geom_segment(aes(x = xMin, y = yMin, xend = xMax, yend = yMax), size = 2, color = "#EB6622", alpha = 0.7, data = gradientLineData) +
  scale_x_continuous("Incidence Rate Ratio", breaks = log(breaks), labels = breaks, limits = xLimits) +
  scale_y_continuous("Log likelihood", limits = yLimits) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "#DDDDDD"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA), 
        strip.placement = "outside",
        strip.background = element_blank())
ggsave("Simulation/gridWithGradients2.png", width = 6, height = 3, dpi = 300)


ggplot(vizData, aes(x = x, y = trueLl)) +
  geom_line(alpha = 0.5, size = 1) +
  geom_point(aes(x = point, y = value), shape = 16, size = 3, color = "#EB6622", data = approximation)+
  geom_line(aes(y = y), color = "#EB6622", size = 2, alpha = 0.6) +
  scale_x_continuous("Incidence Rate Ratio", breaks = log(breaks), labels = breaks, limits = xLimits) +
  scale_y_continuous("Log likelihood", limits = yLimits) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "#DDDDDD"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA), 
        strip.placement = "outside",
        strip.background = element_blank())

ggsave("Simulation/gridWithGradients3.png", width = 6, height = 3, dpi = 300)

# For paper: one plot ------------------------------------------------------------------------------
trueLlAll <- bind_rows(
  vizData |>
    mutate(y = trueLl, 
           label = "True likelihood"),
  vizData |>
    mutate(y = trueLl, 
           label = "Grid with gradients"),
  vizData |>
    mutate(y = trueLl, 
           label = "Hermite interpolation"),
)

gradientLineDataAll <- gradientLineData |>
    mutate(label = "Grid with gradients")

pointsAll <- approximation |>
  mutate(label = "Grid with gradients")

hermiteInterpolationAll <- vizData |>
  mutate(y = y, 
         label = "Hermite interpolation")

setFactor <- function(data) {
  data$label <- factor(data$label, levels = c("True likelihood", "Grid with gradients", "Hermite interpolation"))
  return(data)
}
trueLlAll <- setFactor(trueLlAll)
gradientLineDataAll <- setFactor(gradientLineDataAll)
pointsAll <- setFactor(pointsAll)
hermiteInterpolationAll <- setFactor(hermiteInterpolationAll)

ggplot(trueLlAll, aes(x = x, y = trueLl)) +
  geom_line(alpha = 0.5, size = 0.5) +
  geom_point(aes(x = point, y = value), shape = 16, size = 1.75, color = "#EB6622", data = pointsAll)+
  geom_segment(aes(x = xMin, y = yMin, xend = xMax, yend = yMax), size = 1.25, color = "#EB6622", alpha = 0.7, data = gradientLineDataAll) +
  geom_line(, color = "#EB6622", size = 1.25, alpha = 0.6, data = hermiteInterpolationAll) +
  scale_x_continuous("Incidence Rate Ratio", breaks = log(breaks), labels = breaks, limits = xLimits) +
  scale_y_continuous("Log likelihood", limits = yLimits) +
  facet_grid(~label) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.key = element_blank(),
        legend.background = element_blank(),
        panel.grid.major = element_line(color = "#DDDDDD"), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA), 
        strip.placement = "outside",
        strip.background = element_blank())
ggsave("Simulation/gridWithGradientsAll.png", width = 6, height = 2, dpi = 300)


