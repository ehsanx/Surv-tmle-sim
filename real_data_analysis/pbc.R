n <- 30
# Load necessary packages
# install.packages("survival") # If not already installed
# install.packages("SuperLearner") # If not already installed
# install.packages("dplyr") # If not already installed
# install.packages("mlogit") # If not already installed
# install.packages("polspline") # If not already installed

library(survival)
library(SuperLearner)
library(dplyr)
library(mlogit)       # Dependency for surv.TMLE (if used)
library(polspline)    # Dependency for surv.TMLE (if used)

# --- Load Corrected surv.TMLE Function ---
# Ensure this path points to the script with the corrected function definition
print("Sourcing corrected surv.TMLE function...")
source("../source/SurvTMLE2.R", local = TRUE) # MODIFY PATH AS NEEDED

# --- Load and Prepare PBC Data ---
print("Loading PBC data...")
data(pbc, package = "survival")

# Optional: Subset to randomized trial participants with non-missing treatment
pbc_subset <- pbc[1:312, ]
pbc_subset <- pbc_subset[!is.na(pbc_subset$trt), ]
print(paste("Using", nrow(pbc_subset), "observations from PBC trial."))

# Define Exposure (Avar) - Treatment variable 'trt'
# trt: 1 = D-penicillamine, 2 = placebo
# Let's code A=1 for D-penicillamine, A=0 for Placebo
pbc_subset$A <- ifelse(pbc_subset$trt == 1, 1, 0)
Avar <- "A"

# Define Outcome Event
# status: 0=censored, 1=transplant, 2=dead
# Define event as death
pbc_subset$event <- ifelse(pbc_subset$status == 2, 1, 0)

# Define Time variable (in days)
# time: number of days between registration and the earlier of death,
#       transplantation, or study analysis time
# Use 'time' directly

# Define Discretization Parameters
K <- 30 # Number of intervals (e.g., roughly quarterly for ~7 years)
interval_length <- 85 # days per interval (2550 days / 30 intervals)
max_time <- K * interval_length # Maximum follow-up time in days (~7 years)

print(paste("Discretizing time into K=", K, "intervals of length", interval_length, "days. Max follow-up:", max_time, "days."))

# Cap follow-up time and define event status at max_time
pbc_subset$time_orig <- pbc_subset$time
pbc_subset$time_capped <- pmin(pbc_subset$time, max_time)
# If censored/transplanted after max_time, event is 0 at max_time
# If died after max_time, event is 0 (administratively censored at max_time)
pbc_subset$event[pbc_subset$time > max_time] <- 0
# Keep event=1 if died exactly at max_time? Let's assume event occurs within interval ending at max_time
# event = 1 if status==2 and time <= max_time
# event = 0 if status==0 or status==1 or (status==2 and time > max_time)
pbc_subset$event <- ifelse(pbc_subset$status == 2 & pbc_subset$time <= max_time, 1, 0)


# Create Wide Format Data (Yt, Ct)
wide_data_pbc <- pbc_subset
Yvar <- paste0("Y", 1:K)
Cvar <- paste0("C", 1:K)

for (t in 1:K) {
  time_horizon <- t * interval_length
  # Yt: Died at or before time_horizon
  wide_data_pbc[[Yvar[t]]] <- ifelse(wide_data_pbc$event == 1 & wide_data_pbc$time_capped <= time_horizon, 1, 0)
  # Ct: Censored (status 0 or 1) before time_horizon (and did not die before)
  wide_data_pbc[[Cvar[t]]] <- ifelse(wide_data_pbc$event == 0 & wide_data_pbc$time_capped < time_horizon, 1, 0)
}

# --- Define Covariates (L0var) ---
# Choose baseline covariates expected to be confounders
L0var <- c("age", "sex", "bili", "albumin", "protime", "stage", "edema", "ascites", "hepato", "spiders", "chol", "alk.phos", "ast", "platelet", "trig", "copper")
L0var <- L0var[L0var %in% names(wide_data_pbc)] # Keep only those present

print("Selected baseline covariates (L0var):")
print(L0var)

# Ensure correct variable types
# Factors: sex, ascites, hepato, spiders, stage, edema (coded 0, 0.5, 1)
factor_vars_pbc <- c("sex", "ascites", "hepato", "spiders", "stage", "edema")
factor_vars_pbc <- factor_vars_pbc[factor_vars_pbc %in% L0var] # Check if selected
wide_data_pbc[factor_vars_pbc] <- lapply(wide_data_pbc[factor_vars_pbc], factor)

# Numerics: Ensure others are numeric
numeric_vars_pbc <- L0var[!L0var %in% factor_vars_pbc]
for(var in numeric_vars_pbc) {
  if (!is.numeric(wide_data_pbc[[var]])) {
    print(paste("Converting", var, "to numeric"))
    wide_data_pbc[[var]] <- as.numeric(wide_data_pbc[[var]])
  }
}

# Define Lvar (replicating L0var for time-fixed baseline covariates)
Lvar <- replicate(K, L0var, simplify = FALSE)

# --- Handle Missing Data using na.omit ---
cols_to_check <- c(Avar, L0var, Yvar, Cvar)
cols_to_check <- unique(cols_to_check)
cols_to_check <- cols_to_check[cols_to_check %in% names(wide_data_pbc)]

original_rows <- nrow(wide_data_pbc)
print(paste("Original number of observations:", original_rows))

pbc_complete <- wide_data_pbc[complete.cases(wide_data_pbc[, cols_to_check]), ]

complete_rows <- nrow(pbc_complete)
print(paste("Number of observations after na.omit:", complete_rows))
print(paste("Number of observations removed:", original_rows - complete_rows))

if (complete_rows == 0) {
  stop("No complete cases remaining after na.omit for PBC data. Check input data and variables in cols_to_check.")
}
if ((original_rows - complete_rows) / original_rows > 0.5) {
  warning("More than 50% of PBC observations removed by na.omit. Consider imputation or fewer covariates.")
}


# --- Run TMLE ---
print("Running surv.TMLE on prepared PBC data...")

# Ensure A.SL.library, C.SL.library, Y.SL.library are defined or use SL.library
# Example using SL for all and a simple library
# sl_lib <- c("SL.mean") # Keep it simple for testing

num.scale <- c()
for (col in L0var) {
  x <- pbc_complete[[col]]
  if (is.numeric(x) && max(x) > 1) {
    num.scale <- c(num.scale, col)
  }
}

for (col in L0var) {
  if (col %in% names(pbc_complete)) {
    # If character, convert to factor first
    if (is.character(pbc_complete[[col]])) {
      pbc_complete[[col]] <- as.factor(pbc_complete[[col]])
    }
    
    # Scaling numeric variables if needed
    if (col %in% num.scale && is.numeric(pbc_complete[[col]])) {
      pbc_complete[[col]] <- as.numeric(scale(pbc_complete[[col]]))
    }
    
    pbc_complete[[col]] <- as.numeric(pbc_complete[[col]])
    
    # Check levels and convert accordingly
    if (col == "sex") {
      pbc_complete[[col]] <- as.numeric(pbc_complete[[col]] == "m")
      names(pbc_complete)[names(pbc_complete) == col] <- "sex.male"
    } else {
      pbc_complete[[col]] <- as.numeric(pbc_complete[[col]])
    }
  }
}

L0var[L0var == "sex"] <- "sex.male"
Lvar <- replicate(K, L0var, simplify = FALSE)

start_time <- Sys.time()
result_pbc.para <- surv.TMLE(dat = pbc_complete,
                             Yvar = Yvar[1:n],
                             Cvar = Cvar[1:n],
                             Avar = Avar,
                             Lvar = Lvar[1:n],
                             L0var = NULL, lookback = 1,
                             Ymod = "parametric", 
                             Cmod = "parametric", 
                             Amod = "parametric", 
                             MSM.form = ~ A + time,
                             gbound = 0.025,
                             Print = TRUE)
end_time <- Sys.time()
elapsed.time.para <- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time <- Sys.time()
result_pbc.sl2 <- surv.TMLE(dat = pbc_complete,
                             Yvar = Yvar[1:n],
                             Cvar = Cvar[1:n],
                             Avar = Avar,
                             Lvar = Lvar[1:n],
                             L0var = NULL, lookback = 1,
                             Ymod = "SL", 
                             Cmod = "SL", 
                             Amod = "SL", 
                             SL.library = c("SL.glm", "SL.glm.interaction"), 
                             MSM.form = ~ A + time,
                             gbound = 0.025,
                             Print = TRUE)
end_time <- Sys.time()
elapsed.time.sl2 <- as.numeric(difftime(end_time, start_time, units = "secs"))

start_time <- Sys.time()
result_pbc.sl4 <- surv.TMLE(dat = pbc_complete,
                             Yvar = Yvar[1:n],
                             Cvar = Cvar[1:n],
                             Avar = Avar,
                             Lvar = Lvar[1:n],
                             L0var = NULL, lookback = 1,
                             Ymod = "SL", 
                             Cmod = "SL", 
                             Amod = "SL", 
                             SL.library = c("SL.mean", "SL.glm", "SL.glm.interaction", "SL.earth"), 
                             MSM.form = ~ A + time,
                             gbound = 0.025,
                             Print = TRUE)
end_time <- Sys.time()
elapsed.time.sl4 <- as.numeric(difftime(end_time, start_time, units = "secs"))

# --- Show Results ---
if (!is.null(result_pbc.para)) {
  print("--- PBC Analysis Results ---")
  print("Survival Estimates (St):")
  print(result_pbc.para$St)
  
  if (!is.null(result_pbc.para$MSM)) {
    print("MSM estimates:")
    print(result_pbc.para$MSM)
  }
  
  print("ATE estimates:")
  print(result_pbc.para$ATE)
} else {
  print("surv.TMLE failed to complete.")
}

# --- Show Results ---
if (!is.null(result_pbc.sl4)) {
  print("--- PBC Analysis Results ---")
  print("Survival Estimates (St):")
  print(result_pbc.sl4$St)
  
  if (!is.null(result_pbc.sl4$MSM)) {
    print("MSM estimates:")
    print(result_pbc.sl4$MSM)
  }
  
  print("ATE estimates:")
  print(result_pbc.sl4$ATE)
} else {
  print("surv.TMLE failed to complete.")
}

elapsed_times <- data.frame(elapsed_time = numeric(3))
elapsed_times[1,1] <- elapsed.time.para
elapsed_times[2,1] <- elapsed.time.sl2
elapsed_times[3,1] <- elapsed.time.sl4
rownames(elapsed_times) <- c("parametric", "sl2", "sl4")
elapsed_times.pbc <- elapsed_times
elapsed_times.pbc

ATE_pbc <- tibble(timepoint = c(rep(1:30, 3)),
                  est = -c(result_pbc.para$ATE[2:31],
                           result_pbc.sl2$ATE[2:31],
                           result_pbc.sl4$ATE[2:31]), 
                  se = c(result_pbc.para$SE.ATE[2:31],
                         result_pbc.sl2$SE.ATE[2:31],
                         result_pbc.sl4$SE.ATE[2:31]),
                  lower_ci = -c(result_pbc.para$UL.ATE[2:31],
                                result_pbc.sl2$UL.ATE[2:31],
                                result_pbc.sl4$UL.ATE[2:31]),
                  upper_ci = -c(result_pbc.para$LL.ATE[2:31],
                                result_pbc.sl2$LL.ATE[2:31],
                                result_pbc.sl4$LL.ATE[2:31]),
                  method = c(rep("parametric", 30),
                             rep("SL_2", 30),
                             rep("SL_4", 30)))

dodge_width <- 0.6
plot.ATE_pbc <- ggplot(ATE_pbc, aes(x = timepoint, y = est, color = method)) +
  geom_line(position = position_dodge(width = dodge_width)) + 
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci, color = method), 
                alpha = 0.4,
                position = position_dodge(width = dodge_width)) +
  labs(xlab = "Timepoint", ylab = "ATE", color = "Method",
       title = "PBC ATE over Timepoint by Method")

plot.time_pbc <- ggplot(elapsed_times.pbc, 
                        aes(x = factor(rownames(elapsed_times.pbc)),
                            y = elapsed_time,
                            color = factor(rownames(elapsed_times.pbc)))) +
  geom_col(width = 0.5, fill = "grey") +
  geom_text(aes(label = round(elapsed_time)), 
            vjust = -0.5, size = 3, fontface = "bold") +
  labs(x = "Method", y = "Elapsed Time (s)", 
       title = "PBC Elapsed Time (in Seconds) by Method") +
  theme(legend.position = "none")

ggsave("PBC.ATE_over_Timepoint.png", plot = plot.ATE_pbc, width = 8, height = 6)
ggsave("PBC.Elapsed_Time_by_Method.png", plot = plot.time_pbc, width = 6, height = 6)

