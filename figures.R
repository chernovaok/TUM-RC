<<<<<<< HEAD
# Oksana Chernova, TUM, 2026
# Generate figures for the manuscript

library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)

# Set working directory ----
setwd("//nas.ads.mwn.de/tuma/m12/Projects/TUM MRI biopsy")

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

# Summary by PIRADScat
datapi_summary <- data %>%
  group_by(PIRADScat) %>%
  summarise(count = n())

# Merge counts back to original data
datapirads <- data %>%
  left_join(datapi_summary, by = "PIRADScat")

# Create the scatter plot with jitter 
set.seed(25)
f1 <- ggplot(datapirads, aes(x = factor(PIRADS), y = psa, fill = csPCa)) + 
  geom_jitter(width = 0.3,
              shape = 21,          # Shape 21 allows for a border
              color = "black",     # The boundary color
              size = 2.4,          # Adjust size as needed
              stroke = 0.5) +      # Thickness of the black boundary)
  labs(y = "PSA [ng/ml]", 
       x=""
       ) +
  scale_y_continuous(trans='log10',
                     breaks= c(0.3, 1, 2, 4, 10, 20, 50, 100, 400) 
                     ,labels = function(x) ifelse(x < 1, as.character(x), sprintf("%.0f", x))
  )+
  scale_x_discrete(labels = paste0(datapi_summary$PIRADScat, "\n(n = ", datapi_summary$count, ")"))+
  scale_fill_manual(values = c("csPCa" = "red3", "non-csPCa" = "greenyellow"))+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black", fill=NA)
  )

# Figure 1 ----
png(paste0(getwd(),"/Figures/Fig11.png"), width = 11, height = 8.25, units = "in", res = 300) #, pointsize = 4
f1
dev.off()

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
risk_ny <- read.csv("C:/Users/ochernova/mydocs/calculator/web/msprc.csv")
risk_bc <- read.csv("C:/Users/ochernova/mydocs/calculator/web/BCN2.csv")
risk_er <- read.csv("C:/Users/ochernova/mydocs/calculator/web/ERSPC34.csv")

risk_compl$`MSP-RC` <- risk_ny$cspca_risk
risk_compl$`BCN2-RC` <- risk_bc$BCN2.RC/100
risk_compl$`ERSPC34-RC` <- risk_er$ERSPC34.RC/100

# sample size 
table(is.na(risk_compl$`TUM-RC`), useNA = "always" ) # 1233
table(is.na(risk_compl$`PCRC-MRI`), useNA = "always" ) # 1233
table(is.na(risk_compl$`SPCC`), useNA = "always" ) # 1175
table(is.na(risk_compl$`MSP-RC`), useNA = "always" ) # 1233
table(is.na(risk_compl$`BCN2-RC`), useNA = "always" ) # 1233
table(is.na(risk_compl$`ERSPC34-RC`), useNA = "always" ) # 937


# Calibration ----
# smooth = "gam"
library(pmcalibration)

get_pcc <- function(data, y_var, p_var, smooth = "gam", ci = "boot") {
  
  # 1. Extract the columns from the dataframe
  # We use deparse(substitute()) to handle column names with hyphens like TUM-RC
  y_vec <- data[[deparse(substitute(y_var))]]
  p_vec <- data[[deparse(substitute(p_var))]]
  
  # 2. Run the calibration analysis
  cal_obj <- pmcalibration(
    y = y_vec, 
    p = p_vec, 
    smooth = smooth, 
    ci = ci
  )
  
  # 3. Print the summary if needed
  #cat("\n==============================================\n")
  #cat("Calibration Summary for:", deparse(substitute(p_var)), "\n")
  #print(summary(cal_obj))
  #cat("==============================================\n")
  
  # 4. Extract and return the curve data for plotting
  return(pmcalibration::get_curve(cal_obj, conf_level = 0.95))
}

# 1. TUM-RC
pcctum <- get_pcc(risk_compl, csPCa2, `TUM-RC`)

# 2. PCRC-MRI (UCLA)
pccucla <- get_pcc(risk_compl, csPCa2, `PCRC-MRI`)

# 3. SPCC (Stanford)
risk_spcc <- risk_compl %>% filter( !is.na(`SPCC`))
pccstanford <- get_pcc(risk_spcc, csPCa2, SPCC)

# 4. MSP-RC (NY)
pccny <- get_pcc(risk_compl, csPCa2, `MSP-RC`)

# 5. BCN2-RC
pccba <- get_pcc(risk_compl, csPCa2, `BCN2-RC`)

# 6. ERSPC34-RC
risk_erspc <- risk_compl %>% filter( !is.na(`ERSPC34-RC`))
pccer <- get_pcc(risk_erspc, csPCa2, `ERSPC34-RC`)


# Calibration plots ----
plot_calibration_hist <- function(cal_data, hist_data = risk_compl, x_hist, 
                                  title = "TUM-RC", titley = "Observed (%)") {
  
  ggplot(data = cal_data, aes(x = p, y = p_c, ymin = lower, ymax = upper)) +
    # Reference identity line
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    # Calibration line and confidence ribbon
    geom_line() +
    geom_ribbon(alpha = 1/2, fill = "lightblue") +
    # Force axes to 0-1
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    #scale_y_continuous(labels = function(x) paste0(x * 100, "%"))+
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    scale_x_continuous(labels = function(x) paste0(x * 100))+
    
    # Labels
    labs(
      title = title,
      x = "Predicted (%)",
      y = titley
    ) +
    theme_bw(base_size = 14) +
    # Histogram layer for risk distribution
    geom_histogram(
      data = hist_data, 
      aes(x = {{ x_hist }}, y = after_stat(density) * 0.1), 
      binwidth = 0.05, 
      inherit.aes = FALSE, 
      alpha = 1/2, 
      fill = "grey"
    )+
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(colour = "black", fill=NA))
}

gba <- plot_calibration_hist(
  cal_data = pccba, 
  hist_data = risk_compl, 
  x_hist = `BCN2-RC`, 
  title = "BCN2-RC")

ger <- plot_calibration_hist(
  cal_data = pccer, 
  hist_data = risk_erspc, 
  x_hist = `ERSPC34-RC`, 
  titley = "",
  title = "ERSPC34-RC")

gny <- plot_calibration_hist(
  cal_data = pccny, 
  hist_data = risk_compl, 
  x_hist = `MSP-RC`,
  titley = "",
  title = "MSP-RC")

gucla <- plot_calibration_hist(
  cal_data = pccucla, 
  hist_data = risk_compl, 
  x_hist = `PCRC-MRI`, 
  title = "PCRC-MRI")

gstanford <- plot_calibration_hist(
  cal_data = pccstanford, 
  hist_data = risk_spcc, 
  x_hist = SPCC, 
  titley = "",
  title = "SPCC")

gtum <- plot_calibration_hist(
  cal_data = pcctum, 
  hist_data = risk_compl, 
  x_hist = `TUM-RC`, 
  titley = "",
  title = "TUM-RC")

# Figure 2 ----
png(paste0(getwd(),"/Figures/Fig2.png"), width = 11, height = 8.25, units = "in", res = 300)
gba + ger + gny + gucla + gstanford 
dev.off()

# Risk comparison TUM vs others ----
plot_risk_comparison <- function(data, x_var, y_var, color_var = csPCa,
                                 titley = "TUM-RC (%)") {
  
  # 1. Create the base scatter plot
  p1 <- ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }}, fill = {{ color_var }})) +
    # Add the scatter points
    geom_point(alpha = 1,
               shape = 21,          # Shape 21 allows for a border
               color = "black",     # The boundary color
               size = 1.4,          # Adjust size as needed
               stroke = 0.5) + 
    
    # Add the identity line (y = x)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    
    # Ensure the axes have the same scale and range
    # Force axes to 0-1
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    #scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    #scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    scale_x_continuous(labels = function(x) paste0(x * 100))+
    scale_fill_manual(values = c("csPCa" = "red3", "non-csPCa" = "greenyellow"))+
    # Labels and Theme
    labs(title = "", 
         x = paste0(deparse(substitute(x_var)), " (%)"), 
         y = titley
         ) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16))+
    theme(
      legend.position = "none",
      legend.justification = c("left", "top"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(colour = "black", fill = NA),
      plot.title = element_text(hjust = 0.5)
    )
  return(p1)
}

p_tum_ba <- plot_risk_comparison(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = `BCN2-RC`, 
  color_var = csPCa
)

p_tum_er <- plot_risk_comparison(
  data = risk_erspc, 
  y_var = `TUM-RC`, 
  x_var = `ERSPC34-RC`,
  titley = "",
  color_var = csPCa
)

p_tum_ny <- plot_risk_comparison(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = `MSP-RC`, 
  titley = "",
  color_var = csPCa
)

p_tum_pcrc <- plot_risk_comparison(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = `PCRC-MRI`, 
  color_var = csPCa
)

p_tum_spcc <- plot_risk_comparison(
  data = risk_spcc, 
  y_var = `TUM-RC`, 
  x_var = `SPCC`, 
  titley = "",
  color_var = csPCa
)

# Figure 3 ----
png(paste0(getwd(),"/Figures/Fig3.png"), width = 11, height = 8.25, units = "in", res = 300) #, pointsize = 4
p_tum_ba + p_tum_er + p_tum_ny + p_tum_pcrc + p_tum_spcc 
dev.off()


# Appendix ----
# Panel plot, Spearman correlation
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { 
  usr <- par("usr"); 
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y, method = "spearman", use = "complete.obs") 
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste0(prefix, txt) 
  if(missing(cex.cor)) 
    cex.cor <- 0.8/strwidth(txt) 
  text(0.5, 0.5, txt, cex = cex.cor * r) 
} 

panel.d <- function(x,...) { 
  usr <- par("usr") 
  on.exit(par(usr = usr))
  par(usr=c(usr[1:2],0,5)) 
  #hist(x)
  lines(density(x, na.rm = T), col="palevioletred",lwd=2) 
} 

my_panel_smooth <- function(x, y, ...) {
  panel.smooth(x, y, 
               col.smooth = "black", # Line color
               #lwd = 2,             # Line thickness
               lty = 1,              # Line type (1=solid, 2=dashed)
               ...)
}
r = c(0,1)

rsk <- risk_compl[,c("BCN2-RC", "ERSPC34-RC","MSP-RC", "PCRC-MRI", "SPCC", "TUM-RC", "csPCa")]

# Figure Appendix 1 ----
png(paste0(getwd(),"/Figures/FigAp1.png"), width = 11, height = 8.25, units = "in", res = 300)
pairs(rsk[,c("BCN2-RC", "ERSPC34-RC", "MSP-RC", "PCRC-MRI", "SPCC", "TUM-RC")],
      lower.panel = my_panel_smooth, 
      upper.panel = panel.cor, 
      gap=0, 
      row1attop=FALSE,
      pch=21,
      bg = ifelse(rsk$csPCa == "non-csPCa", "greenyellow", "red3"),
      diag.panel = panel.d
)
dev.off()

# Violin plots ----
# R function for violin-boxplot combination
plot_risk_violin <- function(data, y_var, x_var = csPCa, title = "",
                             titley= "Estimated risk of csPCa (%)") {
  
  ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }}, fill = {{ x_var }})) +
    
    # A. Violin Plot for density shape
    geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
    
    # B. Box Plots for Median and IQR
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") + 
    
    # C. Color and Labels
    scale_fill_manual(values = c("non-csPCa" = "greenyellow", "csPCa" = "red3")) +
    labs(
      title = title,
      y = titley,
      x = "",
      fill = "Group"
    ) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16))+
    # D. Theme adjustments
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(colour = "black", fill=NA),
          legend.position = "none"
    )
}


rba <- plot_risk_violin(
  data = risk_compl, 
  y_var = `BCN2-RC`, 
  x_var = csPCa, 
  title = "BCN2-RC"
)

rer <- plot_risk_violin(
  data = risk_erspc, 
  y_var = `ERSPC34-RC`, 
  x_var = csPCa, 
  titley = "",
  title = "ERSPC34-RC"
)

rny <- plot_risk_violin(
  data = risk_compl, 
  y_var = `MSP-RC`, 
  x_var = csPCa, 
  titley = "",
  title = "MSP-RC"
)

rucla <- plot_risk_violin(
  data = risk_compl, 
  y_var = `PCRC-MRI`, 
  x_var = csPCa, 
  title = "PCRC-MRI"
)

rstanford <- plot_risk_violin(
  data = risk_spcc, 
  y_var = SPCC, 
  x_var = csPCa, 
  titley = "",
  title = "SPCC"
)

rtum <- plot_risk_violin(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = csPCa, 
  titley = "",
  title = "TUM-RC"
)

# Figure Appendix 2 ----
png(paste0(getwd(),"/Figures/FigAp2.png"), width = 11, height = 8.25, units = "in", res = 300)
rba + rer + rny  + rucla + rstanford + rtum 
dev.off()
=======
# Oksana Chernova, TUM, 2026
# Generate figures for the manuscript

library(openxlsx)
library(dplyr)
library(ggplot2)
library(ggExtra)
library(patchwork)

# Set working directory ----
setwd("//nas.ads.mwn.de/tuma/m12/Projects/TUM MRI biopsy")

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

# Summary by PIRADScat
datapi_summary <- data %>%
  group_by(PIRADScat) %>%
  summarise(count = n())

# Merge counts back to original data
datapirads <- data %>%
  left_join(datapi_summary, by = "PIRADScat")

# Create the scatter plot with jitter 
set.seed(25)
f1 <- ggplot(datapirads, aes(x = factor(PIRADS), y = psa, fill = csPCa)) + 
  geom_jitter(width = 0.3,
              shape = 21,          # Shape 21 allows for a border
              color = "black",     # The boundary color
              size = 2.4,          # Adjust size as needed
              stroke = 0.5) +      # Thickness of the black boundary)
  labs(y = "PSA [ng/ml]", 
       x=""
       ) +
  scale_y_continuous(trans='log10',
                     breaks= c(0.3, 1, 2, 4, 10, 20, 50, 100, 400) 
                     ,labels = function(x) ifelse(x < 1, as.character(x), sprintf("%.0f", x))
  )+
  scale_x_discrete(labels = paste0(datapi_summary$PIRADScat, "\n(n = ", datapi_summary$count, ")"))+
  scale_fill_manual(values = c("csPCa" = "red3", "non-csPCa" = "greenyellow"))+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 16))+
  theme(legend.position = "none",
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_rect(colour = "black", fill=NA)
  )

# Figure 1 ----
png(paste0(getwd(),"/Figures/Fig11.png"), width = 11, height = 8.25, units = "in", res = 300) #, pointsize = 4
f1
dev.off()

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
risk_ny <- read.csv("C:/Users/ochernova/mydocs/calculator/web/msprc.csv")
risk_bc <- read.csv("C:/Users/ochernova/mydocs/calculator/web/BCN2.csv")
risk_er <- read.csv("C:/Users/ochernova/mydocs/calculator/web/ERSPC34.csv")

risk_compl$`MSP-RC` <- risk_ny$cspca_risk
risk_compl$`BCN2-RC` <- risk_bc$BCN2.RC/100
risk_compl$`ERSPC34-RC` <- risk_er$ERSPC34.RC/100

# sample size 
table(is.na(risk_compl$`TUM-RC`), useNA = "always" ) # 1233
table(is.na(risk_compl$`PCRC-MRI`), useNA = "always" ) # 1233
table(is.na(risk_compl$`SPCC`), useNA = "always" ) # 1175
table(is.na(risk_compl$`MSP-RC`), useNA = "always" ) # 1233
table(is.na(risk_compl$`BCN2-RC`), useNA = "always" ) # 1233
table(is.na(risk_compl$`ERSPC34-RC`), useNA = "always" ) # 937


# Calibration ----
# smooth = "gam"
library(pmcalibration)

get_pcc <- function(data, y_var, p_var, smooth = "gam", ci = "boot") {
  
  # 1. Extract the columns from the dataframe
  # We use deparse(substitute()) to handle column names with hyphens like TUM-RC
  y_vec <- data[[deparse(substitute(y_var))]]
  p_vec <- data[[deparse(substitute(p_var))]]
  
  # 2. Run the calibration analysis
  cal_obj <- pmcalibration(
    y = y_vec, 
    p = p_vec, 
    smooth = smooth, 
    ci = ci
  )
  
  # 3. Print the summary if needed
  #cat("\n==============================================\n")
  #cat("Calibration Summary for:", deparse(substitute(p_var)), "\n")
  #print(summary(cal_obj))
  #cat("==============================================\n")
  
  # 4. Extract and return the curve data for plotting
  return(pmcalibration::get_curve(cal_obj, conf_level = 0.95))
}

# 1. TUM-RC
pcctum <- get_pcc(risk_compl, csPCa2, `TUM-RC`)

# 2. PCRC-MRI (UCLA)
pccucla <- get_pcc(risk_compl, csPCa2, `PCRC-MRI`)

# 3. SPCC (Stanford)
risk_spcc <- risk_compl %>% filter( !is.na(`SPCC`))
pccstanford <- get_pcc(risk_spcc, csPCa2, SPCC)

# 4. MSP-RC (NY)
pccny <- get_pcc(risk_compl, csPCa2, `MSP-RC`)

# 5. BCN2-RC
pccba <- get_pcc(risk_compl, csPCa2, `BCN2-RC`)

# 6. ERSPC34-RC
risk_erspc <- risk_compl %>% filter( !is.na(`ERSPC34-RC`))
pccer <- get_pcc(risk_erspc, csPCa2, `ERSPC34-RC`)


# Calibration plots ----
plot_calibration_hist <- function(cal_data, hist_data = risk_compl, x_hist, 
                                  title = "TUM-RC", titley = "Observed (%)") {
  
  ggplot(data = cal_data, aes(x = p, y = p_c, ymin = lower, ymax = upper)) +
    # Reference identity line
    geom_abline(intercept = 0, slope = 1, lty = 2) +
    # Calibration line and confidence ribbon
    geom_line() +
    geom_ribbon(alpha = 1/2, fill = "lightblue") +
    # Force axes to 0-1
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    #scale_y_continuous(labels = function(x) paste0(x * 100, "%"))+
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    scale_x_continuous(labels = function(x) paste0(x * 100))+
    
    # Labels
    labs(
      title = title,
      x = "Predicted (%)",
      y = titley
    ) +
    theme_bw(base_size = 14) +
    # Histogram layer for risk distribution
    geom_histogram(
      data = hist_data, 
      aes(x = {{ x_hist }}, y = after_stat(density) * 0.1), 
      binwidth = 0.05, 
      inherit.aes = FALSE, 
      alpha = 1/2, 
      fill = "grey"
    )+
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16))+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(colour = "black", fill=NA))
}

gba <- plot_calibration_hist(
  cal_data = pccba, 
  hist_data = risk_compl, 
  x_hist = `BCN2-RC`, 
  title = "BCN2-RC")

ger <- plot_calibration_hist(
  cal_data = pccer, 
  hist_data = risk_erspc, 
  x_hist = `ERSPC34-RC`, 
  titley = "",
  title = "ERSPC34-RC")

gny <- plot_calibration_hist(
  cal_data = pccny, 
  hist_data = risk_compl, 
  x_hist = `MSP-RC`,
  titley = "",
  title = "MSP-RC")

gucla <- plot_calibration_hist(
  cal_data = pccucla, 
  hist_data = risk_compl, 
  x_hist = `PCRC-MRI`, 
  title = "PCRC-MRI")

gstanford <- plot_calibration_hist(
  cal_data = pccstanford, 
  hist_data = risk_spcc, 
  x_hist = SPCC, 
  titley = "",
  title = "SPCC")

gtum <- plot_calibration_hist(
  cal_data = pcctum, 
  hist_data = risk_compl, 
  x_hist = `TUM-RC`, 
  titley = "",
  title = "TUM-RC")

# Figure 2 ----
png(paste0(getwd(),"/Figures/Fig2.png"), width = 11, height = 8.25, units = "in", res = 300)
gba + ger + gny + gucla + gstanford 
dev.off()

# Risk comparison TUM vs others ----
plot_risk_comparison <- function(data, x_var, y_var, color_var = csPCa,
                                 titley = "TUM-RC (%)") {
  
  # 1. Create the base scatter plot
  p1 <- ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }}, fill = {{ color_var }})) +
    # Add the scatter points
    geom_point(alpha = 1,
               shape = 21,          # Shape 21 allows for a border
               color = "black",     # The boundary color
               size = 1.4,          # Adjust size as needed
               stroke = 0.5) + 
    
    # Add the identity line (y = x)
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    
    # Ensure the axes have the same scale and range
    # Force axes to 0-1
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    #scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    #scale_x_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    scale_x_continuous(labels = function(x) paste0(x * 100))+
    scale_fill_manual(values = c("csPCa" = "red3", "non-csPCa" = "greenyellow"))+
    # Labels and Theme
    labs(title = "", 
         x = paste0(deparse(substitute(x_var)), " (%)"), 
         y = titley
         ) +
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16))+
    theme(
      legend.position = "none",
      legend.justification = c("left", "top"),
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      panel.background = element_rect(colour = "black", fill = NA),
      plot.title = element_text(hjust = 0.5)
    )
  return(p1)
}

p_tum_ba <- plot_risk_comparison(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = `BCN2-RC`, 
  color_var = csPCa
)

p_tum_er <- plot_risk_comparison(
  data = risk_erspc, 
  y_var = `TUM-RC`, 
  x_var = `ERSPC34-RC`,
  titley = "",
  color_var = csPCa
)

p_tum_ny <- plot_risk_comparison(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = `MSP-RC`, 
  titley = "",
  color_var = csPCa
)

p_tum_pcrc <- plot_risk_comparison(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = `PCRC-MRI`, 
  color_var = csPCa
)

p_tum_spcc <- plot_risk_comparison(
  data = risk_spcc, 
  y_var = `TUM-RC`, 
  x_var = `SPCC`, 
  titley = "",
  color_var = csPCa
)

# Figure 3 ----
png(paste0(getwd(),"/Figures/Fig3.png"), width = 11, height = 8.25, units = "in", res = 300) #, pointsize = 4
p_tum_ba + p_tum_er + p_tum_ny + p_tum_pcrc + p_tum_spcc 
dev.off()


# Appendix ----
# Panel plot, Spearman correlation
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) { 
  usr <- par("usr"); 
  on.exit(par(usr = usr))
  par(usr = c(0, 1, 0, 1)) 
  r <- cor(x, y, method = "spearman", use = "complete.obs") 
  txt <- format(c(r, 0.123456789), digits = digits)[1] 
  txt <- paste0(prefix, txt) 
  if(missing(cex.cor)) 
    cex.cor <- 0.8/strwidth(txt) 
  text(0.5, 0.5, txt, cex = cex.cor * r) 
} 

panel.d <- function(x,...) { 
  usr <- par("usr") 
  on.exit(par(usr = usr))
  par(usr=c(usr[1:2],0,5)) 
  #hist(x)
  lines(density(x, na.rm = T), col="palevioletred",lwd=2) 
} 

my_panel_smooth <- function(x, y, ...) {
  panel.smooth(x, y, 
               col.smooth = "black", # Line color
               #lwd = 2,             # Line thickness
               lty = 1,              # Line type (1=solid, 2=dashed)
               ...)
}
r = c(0,1)

rsk <- risk_compl[,c("BCN2-RC", "ERSPC34-RC","MSP-RC", "PCRC-MRI", "SPCC", "TUM-RC", "csPCa")]

# Figure Appendix 1 ----
png(paste0(getwd(),"/Figures/FigAp1.png"), width = 11, height = 8.25, units = "in", res = 300)
pairs(rsk[,c("BCN2-RC", "ERSPC34-RC", "MSP-RC", "PCRC-MRI", "SPCC", "TUM-RC")],
      lower.panel = my_panel_smooth, 
      upper.panel = panel.cor, 
      gap=0, 
      row1attop=FALSE,
      pch=21,
      bg = ifelse(rsk$csPCa == "non-csPCa", "greenyellow", "red3"),
      diag.panel = panel.d
)
dev.off()

# Violin plots ----
# R function for violin-boxplot combination
plot_risk_violin <- function(data, y_var, x_var = csPCa, title = "",
                             titley= "Estimated risk of csPCa (%)") {
  
  ggplot(data, aes(x = {{ x_var }}, y = {{ y_var }}, fill = {{ x_var }})) +
    
    # A. Violin Plot for density shape
    geom_violin(trim = TRUE, alpha = 0.6, scale = "width") +
    
    # B. Box Plots for Median and IQR
    geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") + 
    
    # C. Color and Labels
    scale_fill_manual(values = c("non-csPCa" = "greenyellow", "csPCa" = "red3")) +
    labs(
      title = title,
      y = titley,
      x = "",
      fill = "Group"
    ) +
    scale_y_continuous(labels = function(x) paste0(x * 100))+
    theme(
      axis.title.x = element_text(size = 16),
      axis.text.x = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      axis.title.y = element_text(size = 16))+
    # D. Theme adjustments
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_rect(colour = "black", fill=NA),
          legend.position = "none"
    )
}


rba <- plot_risk_violin(
  data = risk_compl, 
  y_var = `BCN2-RC`, 
  x_var = csPCa, 
  title = "BCN2-RC"
)

rer <- plot_risk_violin(
  data = risk_erspc, 
  y_var = `ERSPC34-RC`, 
  x_var = csPCa, 
  titley = "",
  title = "ERSPC34-RC"
)

rny <- plot_risk_violin(
  data = risk_compl, 
  y_var = `MSP-RC`, 
  x_var = csPCa, 
  titley = "",
  title = "MSP-RC"
)

rucla <- plot_risk_violin(
  data = risk_compl, 
  y_var = `PCRC-MRI`, 
  x_var = csPCa, 
  title = "PCRC-MRI"
)

rstanford <- plot_risk_violin(
  data = risk_spcc, 
  y_var = SPCC, 
  x_var = csPCa, 
  titley = "",
  title = "SPCC"
)

rtum <- plot_risk_violin(
  data = risk_compl, 
  y_var = `TUM-RC`, 
  x_var = csPCa, 
  titley = "",
  title = "TUM-RC"
)

# Figure Appendix 2 ----
png(paste0(getwd(),"/Figures/FigAp2.png"), width = 11, height = 8.25, units = "in", res = 300)
rba + rer + rny  + rucla + rstanford + rtum 
dev.off()
>>>>>>> 57ef5b96a9f867d262e658be1202c59ddc849995
