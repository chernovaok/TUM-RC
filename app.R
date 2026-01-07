library(shiny)
library(dplyr)
library(ggplot2)
library(DT) # For the data table output

# ----------------------------------------------------------------------
# 1. RISK FUNCTIONS ----
# ----------------------------------------------------------------------

# Risk function for TUM-RC  risk calculator
# Outputs risk (0 to 1) of any-grade significant cancer ( ISUP >= 2)

coef <- matrix(c(-2.07771, 0.04981, 0.48253, 1.28603, -1.20072, 0.88053, -1.11740,
                 -1.78074, 0.05017, 0.41618, 1.28363, -1.25050, 0.91657, NA,
                 -2.53834, 0.05433, 0.46774, 1.40723, -1.21223, NA, -1.11284,
                 -6.47440, 0.02571, 0.21028, 1.27058, NA, 1.03662, -1.27395,
                 -2.27881, 0.05417, 0.40401, 1.41712, -1.25813, NA, NA,
                 -6.34199, 0.02497, 0.11358, 1.27386, NA, 1.05844, NA,
                 -7.11820, 0.03122, 0.19274, 1.41818, NA, NA, -1.27009,
                 -7.03464, 0.03036, 0.10207, 1.43306, NA, NA, NA),
               byrow = T, ncol = 7)
colnames(coef) <- c("(Intercept)","age","log2_psa","pirads", "log2_vol", "dre","priorbiopsy")

# age, psa and pirads are mandatory
# age= 55
# psa= 4
# volume = NA
# dre = "abnormal"
# pirads = "3"
# priorbiopsy = "yes"
risk_tum <- function(
    age,           # Numeric: Patient's age at biopsy.
    psa,           # Numeric: Prostate-Specific Antigen level.
    dre,           # DRE Abnormality (0=No, 1=Yes, Unknown).
    priorbiopsy,   # Previous negative biopsy (0=No, 1=Yes, Unknown).
    volume, # Numeric: Prostate Volume.
    pirads  # PIRADS (2-5)
) {
  
  a = is.na(volume)
  b = is.na(dre)
  c = is.na(priorbiopsy)
  
  log2_vol <- ifelse(!a, log(volume, 2), NA)
  dre_abnormal <- ifelse(!b, ifelse(dre == "abnormal", 1, 0), NA)
  priorbiopsy_yes <- ifelse(!c, ifelse(priorbiopsy == "yes", 1, 0), NA)
  log2_psa <- log(psa, 2)
  pirads <- as.numeric(pirads)
  
  vec <- c(age, log2_psa, pirads, log2_vol, dre_abnormal, priorbiopsy_yes)
  
  # Calculate Log-Odds
  log_odds <- sum(coef[is.na(coef[,"log2_vol"]) == a &
                         is.na(coef[,"dre"]) == b &
                         is.na(coef[,"priorbiopsy"]) == c, -1] * vec,
                  na.rm=T) + 
    coef[is.na(coef[,"log2_vol"]) == a &
           is.na(coef[,"dre"]) == b &
           is.na(coef[,"priorbiopsy"]) == c, 1] # intercept
  
  # Calculate risk
  risk.outcome <- exp(log_odds) / (1 + exp(log_odds))
  
  return(risk.outcome)  
}  


# Risk function for PCRC-MRI 
# Kinnaird A, Brisbane W, Kwan L, et al. A prostate cancer risk calculator: Use of clinical and
# magnetic resonance imaging data to predict biopsy outcome in North American men. Can Urol Assoc
# J 2022;16(3):E161-6. http://dx.doi.org/10.5489/cuaj.7380

risk_ucla <- function(
    age,           # Patient's age at biopsy.
    psa,           # Prostate-Specific Antigen level.
    dre,           # DRE Abnormality (0=No/Unknown, 1=Yes).
    priorbiopsy,
    race,   
    volume,        # Prostate Volume 
    pirads 
) {
  
  a = is.na(volume)
  c = is.na(priorbiopsy)
  
  dre_abnormal <- if (dre == "abnormal") 1 else 0
  priorbiopsy_yes <- ifelse(!c, ifelse(priorbiopsy == "yes", 1, 0), NA)
  pirads <- as.numeric(pirads)
  
  # Define Model Parameters
  coef_intercept    <- -4.62990
  coef_age          <- 0.05470
  
  # Race coefficients (0:African American, 1:Asian, 2:Caucasian(Ref), 3:Other, 4:Unknown)
  race_coeff_map <- c(
    "African American" = 0.23780,
    "Asian" = -0.77190,
    "Caucasian" = 0,
    "Other" = 0.18840,
    "Unknown" = -0.08560
  )
  
  coef_PSA          <- 0.01920
  coef_DRE          <- 0.86460
  coef_PrevBiopsy   <- -0.61630
  coef_ProstateVolume <- -0.01790
  coef_PSADensity   <- 0.88020 # Coefficient for binary PSAD variable
  
  # PIRADS score coefficients
  pirads_coeff_map <- c(
    "2" = 0,
    "3" = 0.51020,
    "4" = 1.37510,  
    "5" = 2.76270
  )
  
  coef_pirads <- pirads_coeff_map[as.character(pirads)]
  
  if (!a & !c){
    # Calculate binary PSAD density
    psa_density_value <- psa / volume
    # If PSAD is >= 0.15, the binary variable is 1, otherwise 0
    binary_psad <- ifelse(psa_density_value >= 0.15, 1, 0)
    
    # Calculate MRI Model Log-Odds
    log_odds <- coef_intercept + 
      coef_age * age + 
      race_coeff_map[race] + 
      coef_PSA * psa + 
      coef_DRE * dre_abnormal + 
      coef_PrevBiopsy * priorbiopsy_yes + 
      coef_ProstateVolume * volume + 
      coef_PSADensity * binary_psad + 
      coef_pirads 
    
    # Calculate risk 
    risk.outcome = exp(log_odds) / (1 + exp(log_odds))
  }else{
    risk.outcome <- NA
  }
  return(risk.outcome)
}

# SPCC
risk_stanford <- function(
    age,           # Numeric: Patient's age (ptAge)
    psa,           # Numeric: Prostate-Specific Antigen level (psa)
    volume,        # Numeric: Prostate Volume (prostateVolume)
    history,       # Character: Clinical History ("Biopsy Naive", "Active Surveillance", "Prior Negative Biopsy")
    race,          # Character: Patient Race ("Asian", "Black", "Hispanic", "White", "Other or Unknown")
    pirads         # Character: PIRADS Score ("PIRADS 3", "PIRADS 4", "PIRADS 5")
) {
  
  pirads <- as.numeric(pirads)
  # Define Model coefficients
  
  # Fixed Intercept
  const_Intercept <- -4.326
  
  # PSA Density (psaD) is calculated as psa / volume
  psaD <- psa / volume
  
  # Coefficients and Variances for fixed terms (psaDCoeff and bAge are used in the equation)
  const_psaDCoeff     <- 7.245 # Coefficient for PSA Density
  const_psaDCoeffVar  <- 0.859 # Variance for PSA Density
  const_bAge          <- 3.099 # Coefficient for Age (applied to age/100)
  const_bAgeVar       <- 1.405 # Variance for Age
  
  
  # --- 2. Define Lookups for Categorical Variables (Coefficients and Variances) ---
  
  # Clinical History (aHx and aHxVar)
  history_coeff_map <- c(
    "Biopsy Naive" = 0,
    "Active Surveillance" = 0.139,
    "Prior Negative Biopsy" = -0.6112
  )
  
  # PIRADS Score
  pirads_coeff_map <- c(
    "2" = 0,
    "3" = 0,
    "4" = 1.051,
    "5" = 1.717
  )
  
  # Race (bRace and bRaceVar)
  race_coeff_map <- c(
    "Asian" = -0.1757,
    "Black" = -0.145,
    "Hispanic" = -0.288,
    "White" = 0.4736,
    "Other or Unknown" = 0.0
  )
  
  # --- 3. Extract Coefficients and Variances for the specific patient ---
  
  aHx      <- history_coeff_map[history]
  bPIRADS  <- pirads_coeff_map[as.character(pirads)]
  bRace    <- race_coeff_map[race]
  
  # Check for NA (unmatched input strings)
  if (is.na(aHx) || is.na(bPIRADS) || is.na(bRace)) {
    stop("Invalid input for history, PIRADS, or race. Check spelling.")
  }
  
  # Calculate Log-Odds ---
  
  log_odds <- const_Intercept +
    aHx +
    (psaD * const_psaDCoeff) +
    bRace +
    (const_bAge * age / 100) +
    bPIRADS
  
  # Calculate risk
  risk.outcome <- exp(log_odds) / (1 + exp(log_odds))
  
  return(risk.outcome)
}

# ----------------------------------------------------------------------
# 2. UI DEFINITION -----
# ----------------------------------------------------------------------

ui <- fluidPage(
  tags$head(
    tags$style(HTML("
      .main-header {
        color: #004d40;
        font-weight: 600;
        border-bottom: 2px solid #004d40;
        padding-bottom: 10px;
        margin-bottom: 20px;
      }
      .sidebar {
        background-color: #f7f7f7;
        padding: 15px;
        border-right: 1px solid #ddd;
      }
      .data-table-container {
        margin-top: 20px;
        font-size: 1.1em;
      }

      /* Style for the Calculate Button */
      #calculate_button {
        width: 100%;
        margin-top: 20px;
        background-color: #004d40; 
        color: white; 
        font-weight: bold;
        padding: 10px;
        border-radius: 5px;
        border: none;
      }
      #calculate_button:hover {
        background-color: #00897b;
      }
    "))
  ),
  
  # Sidebar Layout
  sidebarLayout(
    # --- Sidebar Panel for Inputs ---
    
    sidebarPanel(class = "sidebar",
                 width = 4,
                 h3("Patient Characteristics"),
                 # Fixed Inputs
                 numericInput("age", "Age", 55, min = 55.0, max = 90.0),
                 numericInput("psa", "PSA [ng/ml]", 4, min = 0.1, max = 100), 
                 numericInput("volume", "Prostate Volume [ml] (leave blank if unknown)", 40), # Default volume set
                 
                 radioButtons(
                   "pirads",
                   "PIRADS",
                   choices = c("2" = 2,"3" = 3, "4" = 4, "5" = 5),
                   selected = 3,
                   inline = TRUE
                 ),
                 
                 div(
                   #style = "margin-top: 20px; padding-top: 10px; border-top: 1px solid #ddd;",
                   
                   radioButtons(
                     "priornegbiopsy",
                     "Prior negative biopsy:",
                     choices = c("No" = "no", "Yes" = "yes", "Unknown" = NA),
                     selected = "no",
                     inline = TRUE
                   ),
                   radioButtons(
                     "dre",
                     "Digital rectal examination",
                     choices = c("Abnormal" = "abnormal", "Normal" = "normal", "Unknown" = NA),
                     selected = "abnormal",
                     inline = TRUE
                   )
                 ),
                 
                 # --- Action Button for Calculation ---
                 actionButton("calculate_button", "Calculate Risk")
    ),
    
    # --- Main Panel for Outputs ---
    mainPanel(
      width = 8,
      h3(textOutput("result_title")),
      
      # Risk Output Table Container
      div(class = "data-table-container",
          dataTableOutput("risk_table_output")
      ),
      
      p(em(" ")),
      
      p("PCRC-MRI: ",
        a("www.uclahealth.org/departments/urology/iuo/research/prostate-cancer/risk-calculator-mri-guided-biopsy-pcrc-mri", 
          href = "https://www.uclahealth.org/departments/urology/iuo/research/prostate-cancer/risk-calculator-mri-guided-biopsy-pcrc-mri", 
          target = "_blank")
      ),
      p("SPCC: ",
        a("www.med.stanford.edu/ucil/nomogram.html", 
          href = "https://www.med.stanford.edu/ucil/nomogram.html", 
          target = "_blank")
      )
    )
  )
)

# ----------------------------------------------------------------------
# 3. SERVER DEFINITION -----
# ----------------------------------------------------------------------

server <- function(input, output, session) {
  
  # 1. Event Reactive Prediction Calculation (Runs ONLY on button click)
  risk_prediction_event <- eventReactive(input$calculate_button, {
    
    # --- Input Validation ---
    # Ensure all required inputs are present. Volume must be a number > 0.
    validate(
      #need(!is.na(input$volume) && input$volume > 0, "Prostate Volume [ml] must be a valid number greater than 0."),
      need(input$psa > 0, "PSA [ng/ml] must be a valid number greater than 0."),
      need(input$age >= 55 && input$age <= 90, "Age must be between 55 and 90."),
      need(input$volume >= 0 || is.na(input$volume), "Prostate Volume [ml] must be a valid number greater than 0 or missing.")
    )
    
    # --- Calculation ---
    tryCatch({
      # TUM calculation
      prob_tum_val <- risk_tum(
        psa = input$psa, 
        volume = input$volume, 
        age = input$age, 
        pirads = input$pirads, 
        dre = input$dre, 
        priorbiopsy = input$priornegbiopsy
      )
      
      # UCLA calculation 
      prob_ucla_val <- risk_ucla(
        psa = input$psa, 
        volume = input$volume, 
        age = input$age, 
        pirads = input$pirads, 
        dre = input$dre, 
        race = "Caucasian", 
        priorbiopsy = input$priornegbiopsy
      )
      prob_stanf_val <- risk_stanford(
        psa = input$psa,
        volume = input$volume,
        age = input$age,
        pirads = input$pirads,
        race = "White",
        history = ifelse(input$priornegbiopsy == "yes", "Prior Negative Biopsy", "Biopsy Naive")
      )
      # Create the data frame for output
      prob_df = data.frame(
        Risk_Category = c("TUM-RC", "PCRC-MRI", "SPCC"),
        Probability = c(prob_tum_val, prob_ucla_val, prob_stanf_val),
        stringsAsFactors = FALSE
      )
      return(prob_df)
      
    }, error = function(e) {
      # Print error for debugging and return a useful message
      print(paste("Error in risk calculation:", e$message)) 
      showNotification(paste("Calculation Error:", e$message), type = "error")
      return(NULL) # Return NULL on error
    })
  })
  
  # 2. Output Rendering (Title)
  output$result_title <- renderText({
    # Use the event reactive to trigger the title update
    req(risk_prediction_event()) 
    "Risk of clinically significant prostate cancer (ISUP >=2)"
  })
  
  # 3. Output Table (Displays the risks)
  output$risk_table_output <- renderDataTable({
    prob_data <- risk_prediction_event() # Triggered by the action button
    
    if (is.null(prob_data)) {
      return(NULL)
    }
    
    # Format probabilities as percentages
    prob_data$Probability_Pct <- ifelse(!is.na(prob_data$Probability), sprintf("%.1f%%", prob_data$Probability * 100), NA)
    
    # Select and rename columns for display
    display_data <- prob_data %>%
      select(Risk_Category, Probability_Pct) %>%
      #rename(" " = Risk_Category, "" = Probability_Pct)
      rename("Risk Model" = Risk_Category, "Probability" = Probability_Pct)
    
    # Render the data table with nice formatting
    datatable(display_data, 
              options = list(dom = 't', paging = FALSE, ordering = FALSE,
                             headerCallback = JS("function(thead, data, start, end, display){ $(thead).remove(); }"),
                             columnDefs = list(
                               list(className = 'dt-right', targets = c(1))
                             )
              ), 
              rownames = FALSE,
              caption = " ") 
  })
  
}

# Run the application
shinyApp(ui = ui, server = server) # Uncomment this line to run the app
