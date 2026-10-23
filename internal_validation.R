# Oksana Chernova, TUM, 2026
# Internal validation

library(openxlsx)
library(dplyr)
library(pROC)
library(dplyr)

# Read data ----
data <- read.csv(paste0(getwd(),"/MRIbiopsy21.01.2026.csv"))                 

dim(data)

# Create new variables, define factor levels ----
data <- data %>% mutate(
  psa = as.numeric(psa, na.strings = "NA"),
  PIRADS = as.numeric(data$PIRADS, na.strings = "NA"),
  PIRADScat = case_when(
    PIRADS == 0 ~ "PI-RADS=0",
    PIRADS == 1 ~ "PI-RADS=1",
    PIRADS == 2 ~ "PI-RADS=2",
    PIRADS == 3 ~ "PI-RADS=3",
    PIRADS == 4 ~ "PI-RADS=4",
    PIRADS == 5 ~ "PI-RADS=5",
    TRUE ~ NA
  ),
  csPCa = factor(ifelse(ISUP >= 2, "csPCa", "non-csPCa"), levels = c("non-csPCa", "csPCa")),
  csPCa2 = ifelse(ISUP >= 2, 1, 0)
)


# Select some columns
datavar = data[,c("age", "psa", "volume", "PIRADS", "DRE", "priornegbiopsy", "csPCa", "csPCa2")]

datavar <- datavar %>%
  rename(pirads = PIRADS,
         dre = DRE)

data_filt <- datavar %>%
  filter(!is.na(dre) & !is.na(volume) & !is.na(pirads) & !is.na(priornegbiopsy))

data_filt <- data_filt %>%
  mutate(id = row_number())

# ----------------------------------------------------------------------
# RISK FUNCTIONS ----
# ----------------------------------------------------------------------

# Risk function for TUM-RC  risk calculator
# Outputs risk of clinically significant prostate cancer (ISUP >= 2)

coef <- matrix(c(-2.27918, 0.05139, 0.46219, 1.31431, -1.19057, 0.91551, -1.09070,
                 -1.99885, 0.05169, 0.4002,  1.31446, -1.23796, 0.93983, NA,
                 -2.70773, 0.0555,  0.44646, 1.43709, -1.20305, NA, -1.08037,
                 -6.59522, 0.02665, 0.19985, 1.29563, NA, 1.07673, -1.25034,
                 -2.45922, 0.05533, 0.38753, 1.44691, -1.24684, NA, NA,
                 -6.47392, 0.026,   0.10883, 1.30009, NA, 1.08885, NA,
                 -7.22844, 0.03191, 0.18077, 1.44682, NA, NA, -1.23983, 
                 -7.15298, 0.03116, 0.09626, 1.46109, NA, NA, NA),
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
  coef_PSADensity   <- 0.88020 
  
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

# Risk function for SPCC ----
# Wang NN, Zhou SR, Chen L, et al. The stanford prostate cancer calculator: 
# Development and external validation of online nomograms incorporating PI-RADS scores to predict clinically significant prostate cancer. Urol Oncol. 2021;39(12):831.e19-831.e27. doi: 10.1016/j.urolonc.2021.06.004
# https://www.med.stanford.edu/ucil/nomogram.html
risk_stanford <- function(
    age,           # Numeric: Patient's age (ptAge)
    psa,           # Numeric: Prostate-Specific Antigen level (psa)
    volume,        # Numeric: Prostate Volume (prostateVolume)
    history,       # Character: Clinical History ("Biopsy Naive", "Active Surveillance", "Prior Negative Biopsy")
    race,          # Character: Patient Race ("Asian", "Black", "Hispanic", "White", "Other or Unknown")
    pirads         # Character: PIRADS Score ("3", "4", "5")
) {
  
  # Define Model coefficients
  # Fixed Intercept
  const_Intercept <- -4.326
  
  # PSA Density (psaD) is calculated as psa / volume
  psaD <- psa / volume
  
  # Coefficients for psaD and age
  const_psaDCoeff     <- 7.245 # Coefficient for PSA Density
  const_bAge          <- 3.099 # Coefficient for Age (applied to age/100)
  
  # Define Coefficients for Categorical Variables
  
  # Clinical History (aHx)
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
  
  # Race
  race_coeff_map <- c(
    "Asian" = -0.1757,
    "Black" = -0.145,
    "Hispanic" = -0.288,
    "White" = 0.4736,
    "Other or Unknown" = 0.0
  )
  
  # Extract Coefficients
  aHx      <- history_coeff_map[history]
  bPIRADS  <- pirads_coeff_map[as.character(pirads)]
  bRace    <- race_coeff_map[race]
  
  # Check for NA (unmatched input strings)
  if (is.na(aHx) || is.na(bRace)) {
    stop("Invalid input for history, PIRADS, or race. Check spelling.")
  }
  
  pirads <- as.numeric(pirads)
  
  # Calculate risk
  if (pirads %in% 3:5 ){
    log_odds <- const_Intercept +
      aHx +
      (psaD * const_psaDCoeff) +
      bRace +
      (const_bAge * age / 100) +
      bPIRADS
    
    # Calculate risk
    risk.outcome <- exp(log_odds) / (1 + exp(log_odds))
  } else{
    risk.outcome <- NA
  }
  
  return(risk.outcome)
}


# Calculate risk for TUM-RC, PCRC-MRI, SPCC
est_risk <- function(dsin) {
  # Create history vector once
  dsout <- dsin
  dsout$history <- ifelse(dsout$priornegbiopsy == "yes", "Prior Negative Biopsy", "Biopsy Naive")
  
  # Pre-allocate columns with NAs
  n <- nrow(dsout)
  dsout$`TUM-RC`   <- numeric(n)
  dsout$`PCRC-MRI` <- numeric(n)
  dsout$SPCC       <- numeric(n)
  
  for (i in 1:n) {
    # Standardize inputs to variables to make the function calls cleaner
    row <- dsout[i, ]
    
    dsout$`TUM-RC`[i] <- risk_tum(
      psa = row$psa, 
      volume = row$volume, 
      age = row$age, 
      pirads = row$pirads, 
      dre = row$dre, 
      priorbiopsy = row$priornegbiopsy
    )
    
    dsout$`PCRC-MRI`[i] <- risk_ucla(
      psa = row$psa, 
      volume = row$volume, 
      age = row$age, 
      pirads = row$pirads, 
      dre = row$dre, 
      priorbiopsy = row$priornegbiopsy,
      race = "Caucasian")
    
    dsout$SPCC[i] <- risk_stanford(
      psa = row$psa,
      volume = row$volume,
      age = row$age,
      pirads = row$pirads,
      race = "White",
      history = row$history)
  }
  
  return(dsout)
}

risk_compl <- est_risk(data_filt)

# Read web-scraping outputs ----
risk_ny <- read.csv(paste0(getwd(),"/msprc.csv"))
risk_bc <- read.csv(paste0(getwd(),"/BCN2.csv"))
risk_er <- read.csv(paste0(getwd(),"/ERSPC34.csv"))

risk_compl$`MSP-RC` <- risk_ny$cspca_risk
risk_compl$`BCN2-RC` <- risk_bc$BCN2.RC/100
risk_compl$`ERSPC34-RC` <- risk_er$ERSPC34.RC/100

# data set for model
data_comp <- risk_compl |> mutate(
  log2_vol = log(volume, 2),
  log2_psa = log(psa, 2),
  dre = relevel(factor(dre), ref="normal"),
  priornegbiopsy = relevel(factor(priornegbiopsy), ref="no")
)

# Internal validation ---------------------------------------------------------
set.seed(123) # Set seed for reproducibility
# ------------------------------------------------------------------------------
# STEP 1: Assign Stratified 5-Fold Split
# ------------------------------------------------------------------------------
# Stratify based on csPCa2 outcome to maintain balanced cases across folds
data_comp$cv_group <- createFolds(data_comp$csPCa2, k = 5, list = FALSE)

# ------------------------------------------------------------------------------
# STEP 2: Define Risk Prediction Function for Logistic Regression
# ------------------------------------------------------------------------------
risk_cv <- function(test, train) {
  # Fit logistic regression model on training fold
  rg <- glm(
    formula = csPCa2 ~ age + log2_psa + pirads + log2_vol + dre + priornegbiopsy,
    data    = train,
    family  = binomial(link = 'logit')
  )
  
  # Predict absolute risk (predicted probability) on test fold
  test$pred_risk_cv <- predict(rg, newdata = test, type = "response")
  
  return(test)
}

# ------------------------------------------------------------------------------
# STEP 3: Execute Cross-Validation Across 5 Folds
# ------------------------------------------------------------------------------
# Clean, vectorized loop over folds 1 through 5
valid_list <- lapply(1:5, function(k) {
  test_data  <- data_comp[data_comp$cv_group == k, ]
  train_data <- data_comp[data_comp$cv_group != k, ]
  risk_cv(test = test_data, train = train_data)
})

# Recombine validated subsets into a single dataframe
valid <- do.call(rbind, valid_list)

# ------------------------------------------------------------------------------
# Compute Calibration (O/E Ratio) and AUC
# ------------------------------------------------------------------------------

# Calculate Overall Observed-to-Expected (O/E) Ratio
c_oe <- sum(valid$csPCa2) / sum(valid$pred_risk_cv)

# Calculate AUC
roc_obj <- pROC::roc(valid$csPCa2, valid$pred_risk_cv, quiet = TRUE)
auc_tum <- as.numeric(pROC::auc(roc_obj))

# AUC and O/E per fold
fold_summary <- valid %>%
  group_by(cv_group) %>%
  summarise(
    Observed = sum(csPCa2),
    Expected = sum(pred_risk_cv),
    OE_Ratio = round(sum(csPCa2) / sum(pred_risk_cv), 4),
    AUC      = round(as.numeric(auc(roc(csPCa2, pred_risk_cv, quiet = TRUE))), 4)
  )

print(fold_summary)
apply(fold_summary, 2, mean)

# Internal validation using stratified bootstrap -----
B <- 500 # Number of bootstrap resamples
N <- nrow(valid)

boot_auc <- numeric(B)
boot_c <- numeric(B)

for (b in 1:B) {
  set.seed(b)
  # Sample with replacement
  boot_indices <- sample(1:N, size = N, replace = TRUE)
  boot_data    <- valid[boot_indices, ]
  
  # Performance on bootstrap sample
  boot_auc[b]  <- as.numeric(auc(roc(boot_data$csPCa2, boot_data$pred_risk_cv, quiet = TRUE)))
  boot_c[b]    <- sum(boot_data$csPCa2) / sum(boot_data$pred_risk_cv)
}

# ------------------------------------------------------------------------------
# STEP 5: Report Calibration Measures 
# ------------------------------------------------------------------------------
results <- data.frame(
  Metric   = c("AUC TUM", "O/E TUM"),
  value    = c(auc_tum, c_oe),
  CI_Lower = c(quantile(boot_auc, 0.025), quantile(boot_c, 0.025)),
  CI_Upper = c(quantile(boot_auc, 0.975), quantile(boot_c, 0.975))
)

library(flextable)
results |> 
  flextable()|> 
  colformat_double(digits = 2)


# Internal validation using 50:50 bootstrap (not used) -----
# model1 <- glm(formula = csPCa2 ~ age + log2_psa + pirads + log2_vol + dre + priornegbiopsy,
#               data = data_comp,
#               family=binomial(link='logit'))
# summary(model1)
# round(exp(model1$coefficients), 2)
# 
# set.seed(123) # Set seed for reproducibility
# B <- 1000     # Number of bootstrap iterations
# 
# # 1. Separate case and control indices
# case_indices    <- which(data_comp$csPCa2 == 1)
# control_indices <- which(data_comp$csPCa2 == 0)
# 
# n_cases    <- length(case_indices)
# n_controls <- length(control_indices)
# 
# # Vectors to store metrics across iterations
# auc_vec_tum <- numeric(B)
# auc_vec_avg <- numeric(B)
# oe_vec_tum  <- numeric(B)
# oe_vec_avg  <- numeric(B) 
# 
# # Formula for re-fitting
# mod_formula <- csPCa2 ~ age + log2_psa + pirads + log2_vol + dre + priornegbiopsy
# 
# for (i in 1:B) {
#   # 2. Stratified random sampling: 50% cases, 50% controls for training
#   train_cases    <- sample(case_indices, size = floor(n_cases / 2), replace = FALSE)
#   train_controls <- sample(control_indices, size = floor(n_controls / 2), replace = FALSE)
#   
#   train_idx <- c(train_cases, train_controls)
#   
#   # Training and Test sets
#   train_data <- data_comp[train_idx, ]
#   test_data  <- data_comp[-train_idx, ]
#   
#   # 3. Fit new model on training data
#   fit_boot <- glm(mod_formula, data = train_data, family = binomial(link = 'logit'))
#   
#   # 4. Predict probabilities on test data
#   pred_tum <- predict(fit_boot, newdata = test_data, type = "response")
#   actuals  <- test_data$csPCa2
#   
#   # Predict based on SPCC and PCRC-MRI
#   pred_spcc <- test_data$SPCC
#   pred_pcrc <- test_data$`PCRC-MRI`
#   
#   # Calculate average of TUM-RC and SPCC and PCRC-MRI
#   pred_avg <- rowMeans(cbind(pred_tum, pred_spcc, pred_pcrc), na.rm = TRUE)
#   
#   # 5. Evaluate Metrics
#   # AUC TUM
#   roc_obj_tum    <- pROC::roc(actuals, pred_tum, quiet = TRUE)
#   auc_vec_tum[i] <- as.numeric(pROC::auc(roc_obj_tum))
#   # AUC average
#   roc_obj_avg    <- pROC::roc(actuals, pred_avg, quiet = TRUE)
#   auc_vec_avg[i] <- as.numeric(pROC::auc(roc_obj_avg))
#   
#   # Accuracy (Observed / Expected Ratio) 
#   obs_cases <- sum(actuals)
#   # O/E TUM
#   exp_cases_tum <- sum(pred_tum)
#   oe_vec_tum[i]  <- obs_cases / exp_cases_tum
#   # O/E Average
#   exp_cases_avg <- sum(pred_avg)
#   oe_vec_avg[i] <- obs_cases / exp_cases_avg
# }
# 
# # 6. Summarize results with 95% Percentile Confidence Intervals
# results <- data.frame(
#   Metric   = c("AUC TUM", "O/E TUM", "AUC Avg", "O/E Avg"),
#   Mean     = c(mean(auc_vec_tum), mean(oe_vec_tum), 
#                mean(auc_vec_avg), mean(oe_vec_avg)),
#   CI_Lower = c(quantile(auc_vec_tum, 0.025), quantile(oe_vec_tum, 0.025), 
#                quantile(auc_vec_avg, 0.025), quantile(oe_vec_avg, 0.025)),
#   CI_Upper = c(quantile(auc_vec_tum, 0.975), quantile(oe_vec_tum, 0.975),
#                quantile(auc_vec_avg, 0.975), quantile(oe_vec_avg, 0.975))
# )
# 
# library(flextable)
# results |> 
#   flextable()|> 
#   colformat_double(digits = 3)
  
