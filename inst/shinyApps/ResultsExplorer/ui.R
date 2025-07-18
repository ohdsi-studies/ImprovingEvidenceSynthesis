library(shiny)
library(DT)
source("widgets.R")

shinyUI(
  fluidPage(style = "width:1500px;",
            titlePanel("Simulation Results"),
            tabsetPanel(id = "mainTabsetPanel",
                        tabPanel("About",
                                 HTML("<p>TODO</p>"),
                                 HTML("<p>For review use only")
                        ),
                        tabPanel("Fixed-effects simulations",       
                                 fluidRow(
                                   column(3,
                                          checkboxGroupInput("typeFixed", "Meta-analysis algorithm", choices = typesFixed, selected = typesFixed[!grepl("Pooled", typesFixed)]),
                                          lapply(simParamsFixed, createSimParamWidget, results = resultsFixed, suffix = "Fixed"),
                                          checkboxGroupInput("metricFixed", "Metric", choices = metricsFixed, selected = metricsFixed)
                                          ),  
                                   column(9,
                                          tabsetPanel(id = "fixedEffectsTabsetPanel",
                                                      type = "pills",
                                                      tabPanel("Violin plots",
                                                               plotOutput("mainViolinPlotFixed", height = 800),
                                                               div(align = "center",
                                                                   radioButtons("simParamFixedRadioButton", "", choices = simParamsFixed, selected = simParamsFixed[1], inline = TRUE)
                                                               ),
                                                               uiOutput("mainViolinCaptionFixed")
                                                      ),tabPanel("Scatter plots",
                                                               uiOutput("hoverInfoPlotFixed"),
                                                               plotOutput("mainPlotFixed", height = 800, hover = hoverOpts("plotHoverMainPlotFixed", delay = 100, delayType = "debounce")),
                                                               uiOutput("mainCaptionFixed")
                                                      ),
                                                      tabPanel("Rankings",
                                                               plotOutput("rankPlotFixed", height = 800),
                                                               uiOutput("rankCaptionFixed")
                                                      )
                                          )
                                   )
                                 )
                        ),
                        tabPanel("Random-effects simulations",       
                                 fluidRow(
                                   column(3,
                                          checkboxGroupInput("typeRandom", "Meta-analysis algorithm", choices = typesRandom, selected = typesRandom[grepl("random", typesRandom)]),
                                          lapply(simParamsRandom, createSimParamWidget, results = resultsRandom, suffix = "Random"),
                                          checkboxGroupInput("metricRandom", "Metric", choices = metricsRandom, selected = metricsRandom[!grepl("Tau", metricsRandom)])
                                   ), 
                                   column(9,
                                          tabsetPanel(id = "randomEffectsTabsetPanel",
                                                      type = "pills",
                                                      tabPanel("Violin plots",
                                                               plotOutput("mainViolinPlotRandom", height = 800),
                                                               div(align = "center",
                                                                   radioButtons("simParamRandomRadioButton", "", choices = simParamsRandom, selected = simParamsRandom[1], inline = TRUE)
                                                               ),
                                                               uiOutput("mainViolinCaptionRandom")
                                                               ),
                                                      tabPanel("Scatter plots",
                                                               uiOutput("hoverInfoPlotRandom"),
                                                               plotOutput("mainPlotRandom", height = 800, hover = hoverOpts("plotHoverMainPlotRandom", delay = 100, delayType = "debounce")),
                                                               uiOutput("mainCaptionRandom")
                                                      ),
                                                      tabPanel("Rankings",
                                                               plotOutput("rankPlotRandom", height = 800),
                                                               uiOutput("rankCaptionRandom")
                                                      )
                                          )
                                   )
                                 )
                        ),
                        tabPanel("Cohort method random-effects simulations",       
                                 fluidRow(
                                   column(3,
                                          checkboxGroupInput("typeCmRandom", "Meta-analysis algorithm", choices = typesCmRandom, selected = typesRandom[grepl("random", typesCmRandom)]),
                                          lapply(simParamsCmRandom, createSimParamWidget, results = resultsCmRandom, suffix = "CmRandom"),
                                          checkboxGroupInput("metricCmRandom", "Metric", choices = metricsCmRandom, selected = metricsCmRandom[!grepl("Tau", metricsCmRandom)])
                                   ), 
                                   column(9,
                                          tabsetPanel(id = "cmRandomEffectsTabsetPanel",
                                                      type = "pills",
                                                      tabPanel("Violin plots",
                                                               plotOutput("mainViolinPlotCmRandom", height = 800),
                                                               div(align = "center",
                                                                   radioButtons("simParamCmRandomRadioButton", "", choices = simParamsCmRandom, selected = simParamsCmRandom[1], inline = TRUE)
                                                               ),
                                                               uiOutput("mainViolinCaptionCmRandom")
                                                      ),
                                                      tabPanel("Scatter plots",
                                                               uiOutput("hoverInfoPlotCmRandom"),
                                                               plotOutput("mainPlotCmRandom", height = 800, hover = hoverOpts("plotHoverMainPlotCmRandom", delay = 100, delayType = "debounce")),
                                                               uiOutput("mainCaptionCmRandom")
                                                      ),
                                                      tabPanel("Rankings",
                                                               plotOutput("rankPlotCmRandom", height = 800),
                                                               uiOutput("rankCaptionCmRandom")
                                                      )
                                          )
                                   )
                                 )
                        )
            )
  )
)
