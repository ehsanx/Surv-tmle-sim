n <- 30
# Ensure the surv.TMLE function definition (with formula fixes) is loaded
# setwd("E:/GitHub/Surv-tmle") # Set your working directory if needed
source("../source/SurvTMLE2.R") # Make sure this points to the corrected function

# Load required packages
# library(SuperLearner) # Not strictly needed for parametric models
library(dplyr)
library(tidyr)
# library(lubridate) # Not used in this version after date conversion
library(data.table)
library(mlogit)       # Dependency for surv.TMLE (parametric multinomial A)
library(polspline)    # Dependency for surv.TMLE (polyclass for SL multinomial A)
library(SuperLearner) # Dependency for surv.TMLE (if using SL)
library(rstudioapi)
library(ggplot2)

setwd(dirname(getActiveDocumentContext()$path))

# Load RHC dataset
# rhc_data <- read.csv("rhc_data.csv")
rhc_data <- read.csv("https://hbiostat.org/data/repo/rhc.csv", row.names = 1) # Use first column as row names if it's just an index

# --- Data Preparation ---

# Create binary exposure variable
rhc_data$RHC.use <- ifelse(rhc_data$swang1 == "RHC", 1, 0)

# Define event indicator: death (already 0/1 based on source description)
rhc_data$death <- ifelse(rhc_data$death == "Yes", 1, 0) # Already done in source? Check data if needed. Assuming it's 0/1

# Define time-to-event: discharge or death
# Need robust date conversion
# Study Admission Date
# rhc_data$sadmdte <- as.Date(as.character(rhc_data$sadmdte), format = "%m/%d/%Y")
# Hospital Discharge Date
# rhc_data$dschdte <- as.Date(as.character(rhc_data$dschdte), format = "%m/%d/%Y")
# Date of Death
# rhc_data$dthdte <- as.Date(as.character(rhc_data$dthdte), format = "%m/%d/%Y")

rhc_data$time <- as.numeric(rhc_data$dschdte - rhc_data$sadmdte)
# Impute missing discharge date with death date if death occurred
death_indices <- which(rhc_data$death == 1 & is.na(rhc_data$time))
rhc_data$time[death_indices] <- as.numeric(rhc_data$dthdte[death_indices] - rhc_data$sadmdte[death_indices])

# Handle cases where time might still be NA (e.g., alive but no discharge date?)
# Decide on a maximum follow-up or how to handle these otherwise.
# For now, let's assume time is calculated for all or censoring logic handles it.
# If time is still NA after imputation, na.omit will remove them anyway if 'time' is checked,
# but better to understand why they are NA.

# Cap follow-up at 30 days
max_follow_up <- 30
rhc_data$time_orig <- rhc_data$time # Keep original time if needed
rhc_data$time <- pmin(rhc_data$time, max_follow_up, na.rm = TRUE) # Cap time
rhc_data$death[rhc_data$time >= max_follow_up & !is.na(rhc_data$time)] <- 0 # Censor administratively at 30 days

# Create discrete time Y_t and C_t
n_time <- max_follow_up
wide_data <- rhc_data

for (t in 1:n_time) {
  # Event at or before time t
  wide_data[[paste0("Y", t)]] <- ifelse(!is.na(wide_data$death) & wide_data$death == 1 & !is.na(wide_data$time) & wide_data$time <= t, 1, 0)
  # Censored before time t (alive and discharged/administratively censored before t)
  wide_data[[paste0("C", t)]] <- ifelse(!is.na(wide_data$death) & wide_data$death == 0 & !is.na(wide_data$time) & wide_data$time < t, 1, 0)
}


# Force categorical variables to be factor (example)
factor_vars <- c("sex", "dnr1", "race", "income", "ninsclas",
                 "cardiohx", "chfhx", "dementhx", "psychhx",
                 "chrpulhx", "renalhx", "liverhx", "gibledhx",
                 "malighx", "immunhx", "transhx", "amihx")
# Ensure factor_vars exist in wide_data
factor_vars <- factor_vars[factor_vars %in% names(wide_data)]
wide_data[factor_vars] <- lapply(wide_data[factor_vars], factor)


# Define baseline covariates (ensure these columns exist)
L0var <- c("age", "sex", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1", "pafi1",
           "alb1", "hema1", "bili1", "crea1", "sod1", "pot1", "paco21", "ph1", "dnr1", "race",
           "income", "ninsclas", "cardiohx", "chfhx", "dementhx", "psychhx", "chrpulhx", "renalhx",
           "liverhx", "gibledhx", "malighx", "immunhx", "transhx", "amihx")
L0var <- L0var[L0var %in% names(wide_data)] # Keep only covariates present in data

# Ensure numeric types where expected
numeric_vars <- c("age", "aps1", "scoma1", "meanbp1", "wblc1", "hrt1", "resp1", "temp1", "pafi1",
                  "alb1", "hema1", "bili1", "crea1", "sod1", "pot1", "paco21", "ph1")
numeric_vars <- numeric_vars[numeric_vars %in% names(wide_data)]
for(var in numeric_vars) {
  if (is.factor(wide_data[[var]])) {
    warning(paste("Converting factor", var, "to numeric. Check levels carefully."))
    # This conversion might be problematic if factor levels aren't numeric strings
    wide_data[[var]] <- as.numeric(as.character(wide_data[[var]]))
  } else {
    wide_data[[var]] <- as.numeric(wide_data[[var]])
  }
}


# Repeat same covariates across time (only baseline needed for this TMLE version)
Lvar <- replicate(n_time, L0var, simplify = FALSE) # Lvar is list of names

# Define variable name vectors
Yvar <- paste0("Y", 1:n_time)
Cvar <- paste0("C", 1:n_time)
Avar <- "RHC.use"

# --- Handle Missing Data using na.omit ---

# Define columns essential for the analysis (predictors and outcomes/censoring)
# We need Avar, L0var, and all Y/C variables because a missing Y/C at any point invalidates the sequence.
cols_to_check <- c(Avar, L0var, Yvar, Cvar)
cols_to_check <- unique(cols_to_check) # Ensure no duplicates in the check list
cols_to_check <- cols_to_check[cols_to_check %in% names(wide_data)] # Ensure all check columns exist

original_rows <- nrow(wide_data)
print(paste("Original number of observations:", original_rows))

# Keep only rows with complete cases for the essential columns
wide_data_complete <- wide_data[complete.cases(wide_data[, cols_to_check]), ]

complete_rows <- nrow(wide_data_complete)
print(paste("Number of observations after na.omit:", complete_rows))
print(paste("Number of observations removed:", original_rows - complete_rows))

if (complete_rows == 0) {
  stop("No complete cases remaining after na.omit. Check input data and variables in cols_to_check.")
}
if ((original_rows - complete_rows) / original_rows > 0.1) {
  warning("More than 10% of observations removed by na.omit. Consider imputation instead.")
}
dim(wide_data_complete)

cat.drop <- c()
for (col in L0var) {
  x <- wide_data_complete[[col]]
  if ((is.factor(x) || is.character(x)) && length(unique(x)) > 2) {
    cat.drop <- c(cat.drop, col)
  }
}

num.scale <- c()
for (col in L0var) {
  x <- wide_data_complete[[col]]
  if (is.numeric(x) && max(x) > 1) {
    num.scale <- c(num.scale, col)
  }
}

L0var_01 <- setdiff(L0var, cat.drop)

for (col in L0var_01) {
  if (col %in% names(wide_data_complete)) {
    x <- wide_data_complete[[col]]
    
    # If character, convert to factor first
    if (is.character(wide_data_complete[[col]])) {
      wide_data_complete[[col]] <- as.factor(wide_data_complete[[col]])
    }
    
    # Scaling numeric variables if needed
    if (col %in% num.scale && is.numeric(wide_data_complete[[col]])) {
      wide_data_complete[[col]] <- as.numeric(scale(wide_data_complete[[col]]))
    }
    
    # Check levels and convert accordingly
    if (is.factor(wide_data_complete[[col]])) {
      lvls <- levels(wide_data_complete[[col]])
      
      if (all(sort(lvls) == c("No", "Yes"))) {
        wide_data_complete[[col]] <- as.numeric(wide_data_complete[[col]] == "Yes")
      } else if (all(sort(lvls) == c("Female", "Male"))) {
        wide_data_complete[[col]] <- as.numeric(wide_data_complete[[col]] == "Male")
        names(wide_data_complete)[names(wide_data_complete) == col] <- "sex.male"
      } else {
        # General case: convert to numeric
        wide_data_complete[[col]] <- as.numeric(wide_data_complete[[col]])
      }
    }
  }
}

L0var_01[L0var_01 == "sex"] <- "sex.male"
Lvar_01 <- replicate(n_time, L0var_01, simplify = FALSE)
# --- Run TMLE ---

# Ensure the corrected surv.TMLE function is loaded in your environment
start_time <- Sys.time()
result.para <- surv.TMLE(dat = wide_data_complete, # USE THE COMPLETE DATA
                         Yvar = Yvar[1:n],
                         Cvar = Cvar[1:n],
                         Avar = Avar,
                         Lvar = Lvar_01[1:n],  # Pass the list of covariate names
                         L0var = NULL, lookback = 1,
                         Ymod = "parametric",
                         Cmod = "parametric", 
                         Amod = "parametric",
                         MSM.form = ~ RHC.use + time, # Example MSM
                         gbound = 0.025,
                         Print = TRUE)
end_time <- Sys.time()
elapsed.time.para <- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time <- Sys.time()
result.sl2 <- surv.TMLE(dat = wide_data_complete, # USE THE COMPLETE DATA
                        Yvar = Yvar[1:n],
                        Cvar = Cvar[1:n],
                        Avar = Avar,
                        Lvar = Lvar_01[1:n],  # Pass the list of covariate names
                        L0var = NULL, lookback = 1,
                        Ymod = "SL",
                        Cmod = "SL",
                        Amod = "SL",
                        MSM.form = ~ RHC.use + time, # Example MSM
                        SL.library = c("SL.glm", "SL.glm.interaction"),
                        gbound = 0.025,
                        Print = TRUE)
end_time <- Sys.time()
elapsed.time.sl2 <- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time <- Sys.time()
result.sl4 <- surv.TMLE(dat = wide_data_complete, # USE THE COMPLETE DATA
                        Yvar = Yvar[1:n],
                        Cvar = Cvar[1:n],
                        Avar = Avar,
                        Lvar = Lvar_01[1:n],  # Pass the list of covariate names
                        L0var = NULL, lookback = 1,
                        Ymod = "SL",
                        Cmod = "SL",
                        Amod = "SL",
                        MSM.form = ~ RHC.use + time, # Example MSM
                        SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth"),
                        gbound = 0.025,
                        Print = TRUE)
end_time <- Sys.time()
elapsed.time.sl4 <- as.numeric(difftime(end_time, start_time, units = "secs"))

# Show ATE estimates
print("ATE estimates:")
print(result.para$ATE)
print(result.sl2$ATE)
print(result.sl4$ATE)

# Show MSM results if calculated
if (!is.null(result.para$MSM)) {
  print("MSM estimates: (parametric)")
  print(result.para$MSM)
}
if (!is.null(result.para$MSM)) {
  print("MSM estimates: (SL[2])")
  print(result.sl2$MSM)
}
if (!is.null(result.para$MSM)) {
  print("MSM estimates: (SL[4])")
  print(result.sl4$MSM)
}

elapsed_times <- data.frame(elapsed_time = numeric(3))
elapsed_times[1,1] <- elapsed.time.para
elapsed_times[2,1] <- elapsed.time.sl2
elapsed_times[3,1] <- elapsed.time.sl4
rownames(elapsed_times) <- c("parametric", "sl2", "sl4")
elapsed_times.rhc <- elapsed_times

ATE_rhc <- tibble(timepoint = c(rep(1:30, 3)),
                  est = -c(result.para$ATE[2:31],
                           result.sl2$ATE[2:31],
                           result.sl4$ATE[2:31]), 
                  se = c(result.para$SE.ATE[2:31],
                         result.sl2$SE.ATE[2:31],
                         result.sl4$SE.ATE[2:31]),
                  lower_ci = -c(result.para$UL.ATE[2:31],
                                result.sl2$UL.ATE[2:31],
                                result.sl4$UL.ATE[2:31]),
                  upper_ci = -c(result.para$LL.ATE[2:31],
                                result.sl2$LL.ATE[2:31],
                                result.sl4$LL.ATE[2:31]),
                  method = c(rep("parametric", 30),
                             rep("SL_2", 30),
                             rep("SL_4", 30)))

dodge_width <- 0.6
plot.ATE_rhc <- ggplot(ATE_rhc, aes(x = timepoint, y = est, color = method)) +
  geom_line(position = position_dodge(width = dodge_width)) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, color = method), 
                alpha = 0.4,
                position = position_dodge(width = dodge_width)) +
  labs(xlab = "Timepoint", ylab = "ATE", color = "Method",
       title = "RHC ATE over Timepoint by Method")

plot.time_rhc <- ggplot(elapsed_times.rhc, 
                        aes(x = factor(rownames(elapsed_times.rhc)),
                            y = elapsed_time,
                            color = factor(rownames(elapsed_times.rhc)))) +
  geom_col(width = 0.5, fill = "grey") +
  geom_text(aes(label = round(elapsed_time)), 
            vjust = -0.5, size = 3, fontface = "bold") +
  labs(x = "Method", y = "Elapsed Time (s)", 
       title = "RHC Elapsed Time (in Seconds) by Method") +
  theme(legend.position = "none")

ggsave("RHC.ATE_over_Timepoint.png", plot = plot.ATE_rhc, width = 8, height = 6)
ggsave("RHC.Elapsed_Time_by_Method.png", plot = plot.time_rhc, width = 6, height = 6)
