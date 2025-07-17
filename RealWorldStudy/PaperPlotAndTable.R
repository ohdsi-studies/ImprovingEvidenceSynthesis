library(survival)
library(ggplot2)
library(dplyr)
source("RealWorldStudy/DatabaseDetails.R")

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
vizData$type <- factor(vizData$type, levels = c("Likelihood", "Normal", "Hermite interpolation", "Custom"))
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
  scale_size_manual(values = c(1, 2, 2, 2, 2, 2, 2) * 0.5) +
  scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) +
  scale_color_manual(values = c("#000000", RColorBrewer::brewer.pal(4, "Set2"))) +
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

# Plot Hermite interpolation with gradients ------------------------
likelihoodData <- readRDS( file.path(rootFolder, "likelihoodData.rds"))






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


#' @export
plotViolin <- function(group = "Fixed-effects", 
                       types = c("Traditional", "Custom", "Normal", "Skew normal", "Grid"),
                       xParameter = NULL,
                       xLabel = NULL) {
  if (group == "Fixed-effects") {
    data <- readRDS("../inst/shinyApps/ResultsExplorer/data/fixedFx.rds")
    data$type[data$type == "Traditional fixed-effects"] <- "Normal fixed-effects"
    if (is.null(xParameter)) {
      xParameter <- "treatedFraction"
    } 
    if (is.null(xLabel)) {
      xLabel <- "Treated Fraction"
    }
  } else {
    data <- readRDS("../inst/shinyApps/ResultsExplorer/data/randomFx.rds")
    data <- data[!grepl("Tau", data$metric), ]
    data <- data[grepl(paste(group, collapse = "|"), data$type, ignore.case = TRUE), ]
    if (is.null(xParameter)) {
      xParameter <- "nSites"
    } 
    if (is.null(xLabel)) {
      xLabel <- "Number of Sites"
    }    
  }
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
  # vizData$type <- gsub(" ", "\n", vizData$type)
  
  plot <- ggplot(vizData, aes(x = factor(parameterValue), y = value, fill = type)) 
  if (nrow(ref) > 0) {
    plot <- plot + geom_hline(aes(yintercept = value), data = ref, linetype = "dashed")
  }
  plot <- plot + geom_violin(position = position_dodge(0.9), scale = "width", alpha = 0.4) +
    scale_x_discrete(xLabel) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(5, "Set2")) +
    facet_grid(metric~., scales = "free", switch = "both") +
    theme(legend.position = "top",
          legend.title = element_blank(),
          legend.key = element_blank(),
          panel.grid.major.y = element_line(color = "#DDDDDD"), 
          panel.grid.major.x = element_blank(),
          panel.background = element_blank(),
          axis.title.y = element_blank(),
          panel.border = element_rect(fill = NA), 
          strip.placement = "outside",
          strip.background = element_blank())
  # plot
  return(plot)
}

#' @export
createCustomVsGridTable <- function() {
  
  doTest <- function(group) {
    metric <- group$metric[1]
    group$metric <- NULL
    group <- as.data.frame(tidyr::pivot_wider(group, names_from = "type", values_from = "value"))
    customCol <- which(grepl("Custom", colnames(group)))
    gridCol <- which(grepl("Grid", colnames(group)))
    test <- t.test(group[, customCol], group[, gridCol], paired = TRUE, alternative = "two.sided")
    return(data.frame(metric = metric,
                      meanCustom = mean(group[, customCol]),
                      meanGrid =  mean(group[, gridCol]),
                      meanDifference = test$estimate,
                      p = test$p.value,
                      row.names = NULL))
  }
  
  formatTable <- function(results) {
    results$meanCustom <- formatC(results$meanCustom, digits = 3, format = "f")
    results$meanGrid <- formatC(results$meanGrid, digits = 3, format = "f")
    results$meanDifference <- formatC(results$meanDifference, digits = 3, format = "f")
    results$p <- formatC(results$p, digits = 3, format = "f", )
    results$p[grepl("NaN", results$p)] <- ""
    results$p[grepl("0.000", results$p)] <- "$<$0.001"
    results <- results[order(results$metric), ]
    colnames(results) <- c("Metric", "Mean (Custom)", "Mean (Grid)", "Mean diff.", "p")
    rownames(results) <- NULL
    return(results)
  }
  
  data <- readRDS("../inst/shinyApps/ResultsExplorer/data/fixedFx.rds")
  # data <- readRDS("inst/shinyApps/ResultsExplorer/data/fixedFx.rds")
  groups <- split(data, data$metric)
  resultsFixed <- lapply(groups, doTest)
  resultsFixed <- do.call("rbind", resultsFixed) 
  
  data <- readRDS("../inst/shinyApps/ResultsExplorer/data/randomFx.rds")
  # data <- readRDS("inst/shinyApps/ResultsExplorer/data/randomFx.rds")
  data <- data[!grepl("Tau", data$metric), ]
  data <- data[grepl("Custom|Grid", data$type), ]
  data <- data[grepl("random-effects", data$type), ]
  groups <- split(data, data$metric)
  resultsRandom <- lapply(groups, doTest)
  resultsRandom <- do.call("rbind", resultsRandom) 
  
  resultsFixed <- formatTable(resultsFixed)
  resultsRandom <- formatTable(resultsRandom)
  
  
  resultsFixed$Type <- "Fixed-effects"
  resultsRandom$Type <- "Random-effects"
  
  results <- rbind(resultsFixed, resultsRandom)
  results <- cbind(data.frame(Type = results[, ncol(results), ]), results[, -ncol(results)])
  return(results)
  # results <- cbind(resultsFixed, resultsRandom[, -1])
  # 
  # return(kableExtra::add_header_above(
  #     knitr::kable(results, 
  #                  caption = caption, 
  #                  col.names = kableExtra::linebreak(colnames(results)),
  #                  align = c("l", rep("r", 8))),
  #     c("", "Fixed-effects" = 4, "Random-effects" = 4)))
}
