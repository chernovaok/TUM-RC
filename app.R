# Oksana Chernova, TUM, 2026

library(shiny)
library(bslib)
library(dplyr)

# ----------------------------------------------------------------------
# 1. RISK FUNCTIONS ----
# ----------------------------------------------------------------------

# Risk function for TUM-RC  risk calculator
# Outputs risk of clinically significant prostate cancer (ISUP >= 2)

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
  theme = bs_theme(version = 5, primary = "#004d40"),
  
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
      .custom-card {
        background-color: #ffffff;    /* Card background color */
        border: 1px solid #004d40;    /* Teal border to match your header */
        border-radius: 8px;           /* Rounded corners */
        padding: 20px;                /* Space inside the box */
        margin-top: 20px;             /* Space above the box */
        box-shadow: 0 4px 6px rgba(0,0,0,0.1); /* Subtle shadow for depth */
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
                 numericInput("age", "Age", NA, min = 30.0, max = 90.0),
                 numericInput("psa", "PSA [ng/ml]", NA, min = 0.1, max = 100), 
                 numericInput("volume", "Prostate Volume [ml] (leave blank if unknown)", NA), 
                 
                 radioButtons(
                   "pirads",
                   "PI-RADS",
                   choices = c("2" = 2,"3" = 3, "4" = 4, "5" = 5),
                   selected = NA,
                   inline = TRUE
                 ),
                 radioButtons(
                   "priornegbiopsy",
                   "Prior negative biopsy",
                   choices = c("No" = "no", "Yes" = "yes", "Unknown" = NA),
                   selected = character(0),
                   inline = TRUE
                 ),
                 radioButtons(
                   "dre",
                   "Digital rectal examination",
                   choices = c("Abnormal" = "abnormal", "Normal" = "normal", "Unknown" = NA),
                   selected = character(0),
                   inline = TRUE
                 ),
                 
                 # --- Action Button for Calculation ---
                 actionButton("calculate_button", "Calculate Risk")
    ),
    
    # --- Main Panel for Outputs ---
    mainPanel(
      width = 8,
      h3(textOutput("result_title")),
      p(" "),
      
      
      card(
        #class = "border-info",
        #class = "bg-info text-white", # Applies the "secondary" theme color
        #card_header("TUM-RC", class = "fw-bold"),
        card_body(
          div(
            style = "font-size: 1.1em; font-weight: bold;",
            textOutput("tum_risk")
          )
        ),
        markdown("
          **Description:** 67% csPCa, 1228 biopsies, 74% PI-RADS 4,5, Single institution in Munich, Germany <br>
          **Includes:** 11.4% PI-RADS 2 <br>
          **Excludes:** PI-RADS 1, prior PCa diagnosis, PSA ≥ 100 ng/mL <br>
          **Not in tool:** family history, race/ethnicity"
                 #Online risk tool: Not available, only code for local
        )
      ),
      card(
        #class = "bg-secondary text-white", # Applies the "secondary" theme color
        #card_header("SPCC", class = "fw-bold"),
        card_body(
          div(
            style = "font-size: 1.1em; font-weight: bold;",
            textOutput("spcc_risk")
          ),
          markdown("
         **Description:** 54% csPCa, 1922 biopsies, 75%  PI-RADS 4,5, Multiple institutions in Stanford, Yale, and Alabama, USA <br>
         **Includes:** active surveillance patients with prior prostate cancer diagnosis <br>
         **Excludes:** PI-RADS 1 or 2 <br>
         **Not in tool:** DRE, family history <br>
         **Online risk tool:**  https://www.med.stanford.edu/ucil/nomogram.html 
                     "),
        )
      ),
      
      card(
        #class = "bg-secondary text-white", # Applies the "secondary" theme color
        #card_header("PCRC-MRI", class = "fw-bold"),
        card_body(
          div(
            style = "font-size: 1.1em; font-weight: bold;",
            textOutput("ucla_risk")
          )
        ),
        markdown("
         **Description:** 40% csPCa, 2354 biopsies, 52%  PI-RADS 4,5, Multiple institutions Cornell, NY and Los Angeles, CA, USA <br>
         **Includes:** 19% PI-RADS 1, 2<br>
         **Excludes:** prior PCa diagnosis<br>
         **Not in tool:** family history <br>
         **Online risk tool:**  https://www.uclahealth.org/departments/urology/iuo/research/prostate-cancer/risk-calculator-mri-guided-biopsy-pcrc-mri 
        ")
      )
    )
  )
)

# ----------------------------------------------------------------------
# 3. SERVER DEFINITION -----
# ----------------------------------------------------------------------

server <- function(input, output, session) {
  
  # Output Rendering (Title)
  output$result_title <- renderText({
    "Risk Calculators for ISUP ≥ 2 prostate cancer"
  })
  
  # Render text for TUM-RC
  output$tum_risk <- renderText({
    validate(
      need(input$psa > 0 && input$psa <= 100, "PSA [ng/ml] must be a valid number between 0 and 100."),
      need(input$age >= 30 && input$age <= 90, "Age must be between 30 and 90."),
      need((input$volume >= 0 && input$volume <= 300) || is.na(input$volume), "Prostate Volume [ml] must be a valid number between 0 and 300 or missing.")
    )
    
    prob_tum_val <- risk_tum(
      psa = input$psa, 
      volume = input$volume, 
      age = input$age, 
      pirads = input$pirads, 
      dre = input$dre, 
      priorbiopsy = input$priornegbiopsy
    )
    
    out <- ifelse(!is.na(prob_tum_val), sprintf("%.0f%%", prob_tum_val * 100), NA)
    HTML(paste0( "TUM Risk Calculator risk is ", out ))
    # out |>
    #   format(big.mark = ",")
  }) %>% bindEvent(input$calculate_button)
  
  # Render text for SPCC
  output$spcc_risk <- renderText({
    validate(
      need(input$psa > 0 && input$psa <= 100, "PSA [ng/ml] must be a valid number between 0 and 100."),
      need(input$age >= 30 && input$age <= 90, "Age must be between 30 and 90."),
      need((input$volume >= 0 && input$volume <= 300), "Prostate Volume [ml] must be a valid number between 0 and 300."),
      need(input$pirads %in% c("3", "4", "5"), "PI-RADS must be ≥ 3."),
      need(input$priornegbiopsy %in% c("yes", "no"), "Prior negative biopsy can not be unknown.")
    )
    
    prob_stanf_val <- risk_stanford(
      psa = input$psa,
      volume = input$volume,
      age = input$age,
      pirads = input$pirads,
      race = "White",
      history = ifelse(input$priornegbiopsy == "yes", "Prior Negative Biopsy", "Biopsy Naive")
    )
    
    out <- ifelse(!is.na(prob_stanf_val), sprintf("%.0f%%", prob_stanf_val * 100), NA)
    HTML(paste0( "Stanford Prostate Cancer Calculator risk is ", out ))
  }) %>% bindEvent(input$calculate_button)
  
  # Render text for PCRC-MRI
  output$ucla_risk <- renderText({
    validate(
      need(input$psa > 0 && input$psa <= 100, "PSA [ng/ml] must be a valid number between 0 and 100."),
      need(input$age >= 30 && input$age <= 90, "Age must be between 30 and 90."),
      need((input$volume >= 0 && input$volume <= 300), "Prostate Volume [ml] must be a valid number between 0 and 300.")
    )
    
    prob_ucla_val <- risk_ucla(
      psa = input$psa, 
      volume = input$volume, 
      age = input$age, 
      pirads = input$pirads, 
      dre = input$dre, 
      race = "Caucasian", 
      priorbiopsy = input$priornegbiopsy
    )
    out <- ifelse(!is.na(prob_ucla_val), sprintf("%.0f%%", prob_ucla_val * 100), NA)
    
    HTML(paste0(
      "Prostate Cancer Risk Calculator - MRI risk is ",  out))
  }) %>% bindEvent(input$calculate_button)
  
}

# Run the application
shinyApp(ui = ui, server = server) 

