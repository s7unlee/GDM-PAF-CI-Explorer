library(shiny)
library(shinyjs)
library(dplyr)
library(purrr)
library(openxlsx)
library(bslib)

ui <- fluidPage(
  theme = bs_theme(
    bg = "#1a1a1a", fg = "#eaeaea",
    primary = "#f39c12", secondary = "#8e44ad",
    base_font = font_google("Roboto")
  ),
  shinyjs::useShinyjs(),
  navbarPage(
    "GDM-PAF CI Explorer",
    
    # Tab 1: Single Input Calculator
    tabPanel("Single Input Calculator",
             sidebarLayout(
               sidebarPanel(
                 div(id = "infoInputs",
                     helpText("Please enter your name, affiliation, and email for the analysis to proceed. Your information will not be stored."),
                     textInput("inputName", "Name", width = "180px"),
                     textInput("inputAffiliation", "Affiliation", width = "180px"),
                     textInput("inputEmail", "Email", width = "180px"),
                     actionButton("confirmInfo", "Confirm", class = "btn-primary")
                 ),
                 tags$br(), tags$br(),
                 tags$p("Please enter your inputs:"),
                 div(id = "inputFields",
                     fluidRow(
                       column(2, numericInput("inputTp1", "Tp", value = 1, min = 0, width = "120px")),
                       column(2, numericInput("inputPe1", "Pe", value = 0.1, min = 0, max = 1, width = "120px")),
                       column(2, numericInput("inputRR1", "RR", value = 1, min = 0, width = "120px")),
                       column(2, numericInput("inputLower1", "Lower", value = 1, min = 0, width = "120px")),
                       column(2, numericInput("inputUpper1", "Upper", value = 1, min = 0, width = "120px"))
                     )
                 ),
                 actionButton("addRow", "+", class = "btn-secondary", style = "margin-bottom: 10px; font-size: 12px; padding: 4px 8px;"),
                 tags$br(), tags$br(),
                 tags$p("Or upload an Excel file with your inputs:"),
                 fileInput("inputFile", "Upload Excel File", accept = c(".xlsx")),
                 downloadButton("downloadData", "Download Results", class = "btn-primary")
               ),
               mainPanel(
                 tableOutput("outputTable"),
                 verbatimTextOutput("variableDescription"),
                 tags$div(
                   style = "display: inline; color: #eaeaea;",
                   tags$p("The format of the Excel file should be as follows:", style = "display: inline;"),
                   tags$a(href = "https://doi.org/10.3961/jpmph.24.272", "doi.org/10.3961/jpmph.24.272", style = "display: inline; color: #f39c12;"),
                   tags$br(), tags$br(),
                   tags$p("Github:", style = "display: inline;"),
                   tags$a(href = "https://github.com/s7unlee/GDM-PAF-CI-Explorer", "https://github.com/s7unlee/GDM-PAF-CI-Explorer", style = "display: inline; color: #f39c12;")
                 ),
                 tableOutput("excelFormat")
               )
             )
    ),
    
    # Tab 2: Polytomous Exposure Calculator
    tabPanel("Polytomous Exposure Calculator",
             sidebarLayout(
               sidebarPanel(
                 tags$p("Please enter your inputs:"),
                 div(id = "polyInputFields",
                     fluidRow(
                       column(2, numericInput("polyInputTp1", "Tp", value = 1, min = 0, width = "120px")),
                       column(2, numericInput("polyInputPe1", "Pe", value = 0.1, min = 0, max = 1, width = "120px")),
                       column(2, numericInput("polyInputRR1", "RR", value = 1, min = 0, width = "120px")),
                       column(2, numericInput("polyInputLower1", "Lower", value = 1, min = 0, width = "120px")),
                       column(2, numericInput("polyInputUpper1", "Upper", value = 1, min = 0, width = "120px"))
                     )
                 ),
                 actionButton("addPolyRow", "+", class = "btn-secondary", style = "margin-bottom: 10px; font-size: 12px; padding: 4px 8px;"),
                 tags$br(), tags$br(),
                 tags$p("Or upload an Excel file with your inputs:"),
                 fileInput("polyFile", "Upload Excel File", accept = c(".xlsx")),
                 downloadButton("downloadPolyResults", "Download Results", class = "btn-primary")
               ),
               mainPanel(
                 tableOutput("polyOutputTable"),
                 verbatimTextOutput("polyVariableDescription"),
                 tableOutput("polyExcelFormat")
               )
             )
    ),
    
    # Tab 3: PAF Summation Calculator
    tabPanel("PAF Summation Calculator",
             sidebarLayout(
               sidebarPanel(
                 tags$p("Please enter your PAF values:"),
                 div(id = "pafInputFields",
                     fluidRow(
                       column(3, numericInput("pafInput1", "PAF", value = 0.1, min = 0, max = 1, width = "120px"))
                     )
                 ),
                 actionButton("addPafRow", "+", class = "btn-secondary", style = "margin-bottom: 10px; font-size: 12px; padding: 4px 8px;"),
                 tags$br(), tags$br(),
                 tags$p("Or upload an Excel file with your PAF values:"),
                 fileInput("pafFile", "Upload Excel File", accept = c(".xlsx")),
                 downloadButton("downloadPafResults", "Download Results", class = "btn-primary")
               ),
               mainPanel(
                 tableOutput("pafOutputTable"),
                 verbatimTextOutput("pafVariableDescription"),
                 tags$div(
                   style = "display: inline; color: #eaeaea;",
                   tags$p("Example Table Format:", style = "display: inline;"),
                   tags$br(), tags$br(),
                   tableOutput("pafExampleTable"),
                   tags$br(),
                   tags$p("The PAF summation formula from the following study:", style = "display: inline; color: #f39c12;"),
                   tags$a(href = "https://doi.org/10.1093/oxfordjournals.aje.a121617", "10.1093/oxfordjournals.aje.a121617", style = "display: inline; color: #f39c12;")
                 )
               )
             )
    )
  )
)

server <- function(input, output, session) {
  rowCount <- reactiveVal(1)
  polyRowCount <- reactiveVal(1)
  pafRowCount <- reactiveVal(1)
  
  # Add new row of inputs for Single Input Calculator
  observeEvent(input$addRow, {
    rowCount(rowCount() + 1)
    insertUI(
      selector = "#inputFields",
      ui = fluidRow(
        column(2, numericInput(paste0("inputTp", rowCount()), "Tp", value = 1, min = 0, width = "120px")),
        column(2, numericInput(paste0("inputPe", rowCount()), "Pe", value = 0.1, min = 0, max = 1, width = "120px")),
        column(2, numericInput(paste0("inputRR", rowCount()), "RR", value = 1, min = 0, width = "120px")),
        column(2, numericInput(paste0("inputLower", rowCount()), "Lower", value = 1, min = 0, width = "120px")),
        column(2, numericInput(paste0("inputUpper", rowCount()), "Upper", value = 1, min = 0, width = "120px"))
      )
    )
  })
  
  # Add new row of inputs for Polytomous Exposure Calculator
  observeEvent(input$addPolyRow, {
    polyRowCount(polyRowCount() + 1)
    insertUI(
      selector = "#polyInputFields",
      ui = fluidRow(
        column(2, numericInput(paste0("polyInputTp", polyRowCount()), "Tp", value = 1, min = 0, width = "120px")),
        column(2, numericInput(paste0("polyInputPe", polyRowCount()), "Pe", value = 0.1, min = 0, max = 1, width = "120px")),
        column(2, numericInput(paste0("polyInputRR", polyRowCount()), "RR", value = 1, min = 0, width = "120px")),
        column(2, numericInput(paste0("polyInputLower", polyRowCount()), "Lower", value = 1, min = 0, width = "120px")),
        column(2, numericInput(paste0("polyInputUpper", polyRowCount()), "Upper", value = 1, min = 0, width = "120px"))
      )
    )
  })
  
  # Add new row of inputs for PAF Summation Calculator
  observeEvent(input$addPafRow, {
    pafRowCount(pafRowCount() + 1)
    insertUI(
      selector = "#pafInputFields",
      ui = fluidRow(
        column(3, numericInput(paste0("pafInput", pafRowCount()), "PAF", value = 0.1, min = 0, max = 1, width = "120px"))
      )
    )
  })
  
  
  # Tab 1: Single Input Calculator
  reactiveData <- reactive({
    req(input$inputName, input$inputAffiliation, input$inputEmail)
    
    if (!is.null(input$inputFile)) {
      data <- read.xlsx(input$inputFile$datapath)
      data$CIR <- data$Upper / data$Lower
    } else {
      Tp <- sapply(1:rowCount(), function(i) input[[paste0("inputTp", i)]])
      RR <- sapply(1:rowCount(), function(i) input[[paste0("inputRR", i)]])
      Lower <- sapply(1:rowCount(), function(i) input[[paste0("inputLower", i)]])
      Upper <- sapply(1:rowCount(), function(i) input[[paste0("inputUpper", i)]])
      Pe <- sapply(1:rowCount(), function(i) input[[paste0("inputPe", i)]])
      CIR <- Upper / Lower
      data <- data.frame(Tp = Tp, RR = RR, CIR = CIR, Pe = Pe)
    }
    
    data1 <- dplyr::rename(data, Tp = Tp, Pe = Pe, RR = RR, CIR = CIR) %>%
      dplyr::mutate(
        beta = log(RR),
        Var.Pe = (Pe * (1 - Pe)) / Tp,
        O = Pe / (1 - Pe),
        logse = log(CIR) / (2 * qnorm(0.975)),
        Z = abs(beta / logse),
        Pval = exp(-0.717 * Z - 0.416 * Z * Z),
        Var.beta = (logse)^2,
        AF = (((Pe * (RR - 1)) / (Pe * (RR - 1) + 1))),
        Delta.Var.AF = (((RR - 1)^2) * Var.Pe + ((Pe * RR)^2) * Var.beta) / ((1 + Pe * (RR - 1))^4),
        Delta.low = AF - (qnorm(1 - 0.025)) * sqrt(Delta.Var.AF),
        Delta.up = AF + (qnorm(1 - 0.025)) * sqrt(Delta.Var.AF),
        Green.Var.AF = ((O / Tp) * ((RR - 1)^2 + Var.beta * RR^2) + Var.beta * O^2 * RR^2) * (1 - AF)^2 / (1 + O)^2,
        Green.low = (1 - (1 - (AF)) * exp((qnorm(1 - 0.025)) * sqrt(Green.Var.AF))),
        Green.up = (1 - (1 - (AF)) * exp(-(qnorm(1 - 0.025)) * sqrt(Green.Var.AF)))
      )
    
    Monte <- dplyr::select(data1, beta, Var.beta, Tp, Pe)
    
    Monte.low <- c()
    Monte.up <- c()
    Monte.RR <- c()
    Monte.AF <- c()
    Monte.Pe <- c()
    
    for (i in 1:nrow(Monte)) {
      if (i %% 1000 == 0) {
        cat(paste0("Monte CI Processed ", i, " rows...", "\n\n"))
      }
      
      beta <- data1$beta[i]
      Var.beta <- data1$Var.beta[i]
      Tp <- data1$Tp[i]
      Pe <- data1$Pe[i]
      
      set.seed(1234)
      MonteRR <- rlnorm(10000, beta, sqrt(Var.beta))
      MontePe <- rbinom(10000, Tp, Pe) / Tp
      MonteAF <- (MontePe * (MonteRR - 1.0)) / (1.0 + MontePe * (MonteRR - 1.0))
      
      Monte.AF[i] <- quantile(MonteAF, 0.5)
      Monte.RR[i] <- quantile(MonteRR, 0.5)
      Monte.Pe[i] <- quantile(MontePe, 0.5)
      Monte.low[i] <- quantile(MonteAF, 0.5 - 0.475)
      Monte.up[i] <- quantile(MonteAF, 0.5 + 0.475)
    }
    
    monte.results <- cbind(Monte.RR, Monte.Pe, Monte.AF, Monte.low, Monte.up) %>% data.frame
    colnames(monte.results) <- c("Monte.RR", "Monte.Pe", "Monte.AF", "Monte.low", "Monte.up")
    
    data2 <- cbind(data1, monte.results)
    
    delta <- 0.01
    
    compute_delta_var_AF <- function(Tp, Pe, RR, Var_beta) {
      Var.Pe = (Pe * (1 - Pe)) / Tp
      (((RR - 1)^2) * Var.Pe + ((Pe * RR)^2) * Var_beta) / ((1 + Pe * (RR - 1))^4)
    }
    
    compute_green_var_AF <- function(Tp, Pe, RR, Var_beta, AF) {
      O = Pe / (1 - Pe)
      ((O / Tp) * ((RR - 1)^2 + Var_beta * RR^2) + Var_beta * O^2 * RR^2) * (1 - AF)^2 / (1 + O)^2
    }
    
    data3 <- data2 %>%
      mutate(
        Delta.Var.AF.Tp = map_dbl(1:n(), ~compute_delta_var_AF(Tp[.x] * (1 + delta), Pe[.x], RR[.x], Var.beta[.x])),
        Delta.Var.AF.Pe = map_dbl(1:n(), ~compute_delta_var_AF(Tp[.x], Pe[.x] * (1 + delta), RR[.x], Var.beta[.x])),
        Delta.Var.AF.RR = map_dbl(1:n(), ~compute_delta_var_AF(Tp[.x], Pe[.x], RR[.x] * (1 + delta), Var.beta[.x])),
        Delta.Var.AF.Var.beta = map_dbl(1:n(), ~compute_delta_var_AF(Tp[.x], Pe[.x], RR[.x], Var.beta[.x] * (1 + delta)))
      ) %>%
      mutate(
        Delta.Sen.Tp = (Delta.Var.AF.Tp - Delta.Var.AF) / (Tp * delta),
        Delta.Sen.Pe = (Delta.Var.AF.Pe - Delta.Var.AF) / (Pe * delta),
        Delta.Sen.RR = (Delta.Var.AF.RR - Delta.Var.AF) / (RR * delta),
        Delta.Sen.Var.beta = (Delta.Var.AF.Var.beta - Delta.Var.AF) / (Var.beta * delta)
      ) %>%
      mutate(
        Green.Var.AF.Tp = map_dbl(1:n(), ~compute_green_var_AF(Tp[.x] * (1 + delta), Pe[.x], RR[.x], Var.beta[.x], AF[.x])),
        Green.Var.AF.Pe = map_dbl(1:n(), ~compute_green_var_AF(Tp[.x], Pe[.x] * (1 + delta), RR[.x], Var.beta[.x], AF[.x])),
        Green.Var.AF.RR = map_dbl(1:n(), ~compute_green_var_AF(Tp[.x], Pe[.x], RR[.x] * (1 + delta), Var.beta[.x], AF[.x])),
        Green.Var.AF.Var.beta = map_dbl(1:n(), ~compute_green_var_AF(Tp[.x], Pe[.x], RR[.x], Var.beta[.x] * (1 + delta), AF[.x]))
      ) %>%
      mutate(
        Green.Sen.Tp = (Green.Var.AF.Tp - Green.Var.AF) / (Tp * delta),
        Green.Sen.Pe = (Green.Var.AF.Pe - Green.Var.AF) / (Pe * delta),
        Green.Sen.RR = (Green.Var.AF.RR - Green.Var.AF) / (RR * delta),
        Green.Sen.Var.beta = (Green.Var.AF.Var.beta - Green.Var.AF) / (Var.beta * delta)
      ) %>% data.frame
    
    data3
  })
  
  output$outputTable <- renderTable({
    reactiveData()
  })
  
  output$excelFormat <- renderTable({
    data.frame(
      Tp = c(1e+03, 1e+03, 1e+04, 1e+05, 1e+03, 1e+04, 1e+05, 1e+03, 1e+05, 1e+05, 1e+04, 1e+04),
      Pe = c(0.01, 0.90, 0.50, 0.10, 0.10, 0.01, 0.90, 0.10, 0.90, 0.10, 0.01, 0.50),
      RR = c(1.2, 1.2, 2.0, 5.0, 10.0, 1.2, 1.2, 1.2, 2.0, 1.2, 5.0, 1.2),
      Lower = c(1.0524696, 1.0524696, 1.7541160, 4.3852901, 8.7705802, 0.9486833, 0.9486833, 0.8485281, 1.4142136, 0.6928203, 2.8867513, 0.3794733),
      Upper = c(1.368211, 1.368211, 2.280351, 5.700877, 11.401754, 1.517893, 1.517893, 1.697056, 2.828427, 2.078461, 8.660254, 3.794733)
    )
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(reactiveData(), file)
    }
  )
  
  output$variableDescription <- renderText({
    paste(" TP: Total population", "\n",
          "Pe: prevalence of the risk factor in the populaiton", "\n",
          "RR: relative risk", "\n",
          "Lower: Lower bound of 95% CI for RR", "\n",
          "Upper: Upper bound of 95% CI for RR", "\n",
          "Var.Pe: variance of Pe", "\n",
          "O: Odds of Pe", "\n",
          "logse: the standard error of log(RR)", "\n",
          "Z: Absolute value of beta divided by logse",
          "Pval: P-value", "\n",
          "Var.beta: square of logse", "\n",          
          "AF: Attributable fraction", "\n",          
          "Delta.Var.AF, Delta.low, Delta.up: Variance, lower and upper 95% CI of AF using Delta method", "\n",                    
          "Green.Var.AF, Green.low, Green.up: Variance, lower and upper 95% CI of AF using Greenland method", "\n",                  
          "Monte.RR, Monte.Pe, Monte.AF, Monte.low, Monte.up: Median and 95% CI for the RR, Pe, and AF from the Monte Carlo method", "\n",                  
          "Sensitivity results: Derivatives of AF, Delta.Var.AF, and Green.Var.AF with respect to each of the inputs (beta, Var.beta, Tp, Pe) are calculated.", "\n"                  
    )
  })
  
  observeEvent(input$confirmInfo, {
    shinyjs::hide("infoInputs")
  })
  
  # Tab 2: Polytomous Exposure Calculator
  reactivePolyData <- reactive({
    if (!is.null(input$polyFile)) {
      data <- read.xlsx(input$polyFile$datapath)
      validate(
        need(all(c("Category", "Tp", "Pe", "RR", "Lower", "Upper") %in% colnames(data)), 
             "The uploaded file must contain the columns: Category, Tp, Pe, RR, Lower, Upper.")
      )
    } else {
      Tp <- sapply(1:polyRowCount(), function(i) input[[paste0("polyInputTp", i)]])
      RR <- sapply(1:polyRowCount(), function(i) input[[paste0("polyInputRR", i)]])
      Lower <- sapply(1:polyRowCount(), function(i) input[[paste0("polyInputLower", i)]])
      Upper <- sapply(1:polyRowCount(), function(i) input[[paste0("polyInputUpper", i)]])
      Pe <- sapply(1:polyRowCount(), function(i) input[[paste0("polyInputPe", i)]])
      data <- data.frame(Tp = Tp, RR = RR, Lower = Lower, Upper = Upper, Pe = Pe)
    }
    
    # Polytomous Exposure Equation
    data <- data %>%
      mutate(
        Pi_RR_minus1 = Pe * (RR - 1)  # p_i * (RR_i - 1)
      ) 
    
    # Summation of Pi_RR_minus1 across all categories
    sum_Pi_RR_minus1 <- sum(data$Pi_RR_minus1)
    PAF <- sum_Pi_RR_minus1 / (sum_Pi_RR_minus1 + 1)
    
    # Monte Carlo simulation
    set.seed(1234)
    Monte.PAF <- replicate(10000, {
      MonteRR <- rlnorm(nrow(data), log(data$RR), log(data$Upper / data$Lower) / (2 * qnorm(0.975)))
      MontePe <- rbinom(nrow(data), data$Tp, data$Pe) / data$Tp
      MontePi_RR_minus1 <- MontePe * (MonteRR - 1)
      MonteSum_Pi_RR_minus1 <- sum(MontePi_RR_minus1)
      MonteSum_Pi_RR_minus1 / (MonteSum_Pi_RR_minus1 + 1)
    })
    
    # Monte Carlo CI
    Monte.PAF.low <- quantile(Monte.PAF, 0.025)
    Monte.PAF.up <- quantile(Monte.PAF, 0.975)
    Monte.PAF.median <- median(Monte.PAF)
    
    # Add results to the dataframe
    data <- data %>%
      mutate(Sum_Pi_RR_minus1 = sum_Pi_RR_minus1,
             PAF = PAF,
             Monte.PAF.median = Monte.PAF.median,
             Monte.PAF.low = Monte.PAF.low,
             Monte.PAF.up = Monte.PAF.up)
    
  })
  
  output$polyOutputTable <- renderTable({
    reactivePolyData()
  })
  
  output$polyExcelFormat <- renderTable({
    data.frame(
      Category = c("<18.5", "25.0-26.9", "27.0-29.9", "\u226530"),
      Tp = c(1000.00, 2000.00, 1500.00, 500.00),
      Pe = c(0.20, 0.30, 0.25, 0.10),
      RR = c(1.50, 2.00, 2.50, 3.00),
      Lower = c(1.20, 1.80, 2.20, 2.70),
      Upper = c(1.80, 2.20, 2.80, 3.30)
    )
  })
  
  output$downloadPolyResults <- downloadHandler(
    filename = function() {
      paste("polytomous_results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(reactivePolyData(), file)
    }
  )
  
  output$polyVariableDescription <- renderText({
    paste(
      "The Polytomous Exposure Equation is used to calculate Population Attributable Fractions (PAF) as follows:",
      "Sum_Pi_RR_minus1 = Σ (Pe_i * (RR_i - 1))",
      "PAF = Sum_Pi_RR_minus1 / (Sum_Pi_RR_minus1 + 1)",
      "Input columns in the Excel file must include: Category, Tp, Pe, RR, Lower, and Upper.",
      sep = "\n"
    )
  })
  
  
  # Tab 3: Reactive data for PAF Summation Calculator
  reactivePafData <- reactive({
    if (!is.null(input$pafFile)) {
      data <- read.xlsx(input$pafFile$datapath)
      validate(
        need("PAF" %in% colnames(data), "The uploaded file must contain a column named 'PAF'.")
      )
      pafValues <- data$PAF
    } else {
      pafValues <- sapply(1:pafRowCount(), function(i) input[[paste0("pafInput", i)]])
    }
    pafSum <- 1 - prod(1 - pafValues, na.rm = TRUE)
    data.frame(PAF_Values = pafValues, Summed_PAF = pafSum)
  })
  
  output$pafOutputTable <- renderTable({
    reactivePafData()
  })
  
  output$pafExampleTable <- renderTable({
    data.frame(
      PAF = c(0.1, 0.2, 0.15)
    )
  })
  
  output$pafVariableDescription <- renderText({
    "The formula for PAF summation is 1 - (1 - PAFa) * (1 - PAFb) * ... . Provide PAF values either manually or via Excel upload."
  })
  
  output$downloadPafResults <- downloadHandler(
    filename = function() {
      paste("paf_results-", Sys.Date(), ".xlsx", sep = "")
    },
    content = function(file) {
      write.xlsx(reactivePafData(), file)
    }
  )
}

shinyApp(ui = ui, server = server)
