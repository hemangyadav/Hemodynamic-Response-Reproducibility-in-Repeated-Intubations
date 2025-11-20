library(pacman)
p_load(tidyverse, data.table, lme4, performance, broom.mixed, lmerTest)

# Load the analysis results ------------
icc_results <- readRDS("icc_analysis_results.rds")
paired_cohort <- readRDS("cleaned_datasets.rds")$paired_cohort

# 1. Fit revised model and extract coefficients -------------
cat("1. Fitting revised model (without CHF and IntubationNumber)...\n")

# Fit model with lmerTest for p-values
model_revised <- lmerTest::lmer(EffectiveMAP_bin_0.5 ~ 
                                  Baseline_EffMAP +
                                  SOFA + AGE_AT_ADMISSION + PATIENT_SEX + 
                                  HOSPITAL_ADMISSION_WEIGHT_KG +
                                  EJECTION_FRACTION_THREE_MONTHS + ef_missing +
                                  HTN + CKD + 
                                  (1 | StudyID),
                                data = paired_cohort)

# Extract fixed effects with confidence intervals and p-values
fixed_effects <- broom.mixed::tidy(model_revised, effects = "fixed", conf.int = TRUE, conf.level = 0.95)
setDT(fixed_effects)

# Clean up variable names 
fixed_effects[, Variable := c(
  "Intercept",
  "Baseline Effective MAP (per mmHg)",
  "SOFA Score (per point)",
  "Age (per year)",
  "Male Sex",
  "Weight (per kg)",
  "Ejection Fraction (per %)",
  "Missing Ejection Fraction",
  "Hypertension",
  "Chronic Kidney Disease"
)]

# Format table
model_table <- fixed_effects[, .(
  Variable,
  Estimate = round(estimate, 3),
  SE = round(std.error, 3),
  `95% CI Lower` = round(conf.low, 3),
  `95% CI Upper` = round(conf.high, 3),
  `t-value` = round(statistic, 3),
  `P-value` = ifelse(p.value < 0.001, "<0.001", sprintf("%.3f", p.value))
)]

print(model_table)
fwrite(model_table, "Table3_Adjusted_Model_Coefficients.csv")

# 2. Extract random effects ------------
cat("\n2. Extracting random effects variance components...\n")

vc <- as.data.frame(VarCorr(model_revised))
var_between <- vc$vcov[vc$grp == "StudyID"]
var_within <- vc$vcov[vc$grp == "Residual"]
sd_between <- vc$sdcor[vc$grp == "StudyID"]
sd_within <- vc$sdcor[vc$grp == "Residual"]
icc_value <- var_between / (var_between + var_within)

random_effects <- data.table(
  Component = c(
    "Between-patient variance (patient-specific intercept)",
    "Within-patient residual variance",
    "ICC (proportion of variance between patients)"
  ),
  Variance = c(round(var_between, 2), round(var_within, 2), NA),
  `Standard Deviation` = c(round(sd_between, 3), round(sd_within, 3), NA),
  ICC = c(NA, NA, round(icc_value, 4))
)

print(random_effects)
fwrite(random_effects, "Table3_Random_Effects.csv")

# 3. ICCs ---------
# Get complete cases for all three timepoints
complete_data <- paired_cohort[complete.cases(paired_cohort[, .(
  EffectiveMAP_bin_0.5, EffectiveMAP_bin_1, EffectiveMAP_bin_1.5,
  Baseline_EffMAP, SOFA, AGE_AT_ADMISSION, PATIENT_SEX, 
  HOSPITAL_ADMISSION_WEIGHT_KG, EJECTION_FRACTION_THREE_MONTHS,
  ef_missing, HTN, CKD, StudyID
)])]

# Refactor StudyID to drop unused levels
complete_data$StudyID <- factor(complete_data$StudyID)

results <- data.table()

for(tp in c("0.5", "1", "1.5")) {
  outcome_var <- paste0("EffectiveMAP_bin_", tp)
  
  formula_str <- paste0(outcome_var, " ~ Baseline_EffMAP + SOFA + AGE_AT_ADMISSION + ",
                        "PATIENT_SEX + HOSPITAL_ADMISSION_WEIGHT_KG + ",
                        "EJECTION_FRACTION_THREE_MONTHS + ef_missing + HTN + CKD + (1 | StudyID)")
  
  model <- lmer(as.formula(formula_str), data = complete_data, REML = TRUE)
  
  # Get variance components
  vc <- VarCorr(model)
  var_between <- as.numeric(attr(vc$StudyID, "stddev"))^2
  var_within <- attr(vc, "sc")^2
  icc <- var_between / (var_between + var_within)
  
  # Sample sizes
  n_patients <- nlevels(complete_data$StudyID)
  n_obs <- nrow(complete_data)
  
  # Standard error using delta method
  se_var_between <- sqrt(2 * var_between^2 / n_patients)
  se_var_within <- sqrt(2 * var_within^2 / n_obs)
  
  var_total <- var_between + var_within
  d_icc_d_vb <- var_within / (var_total^2)
  d_icc_d_vw <- -var_between / (var_total^2)
  
  se_icc <- sqrt((d_icc_d_vb^2 * se_var_between^2) + (d_icc_d_vw^2 * se_var_within^2))
  
  # 95% CI
  icc_lower <- max(0, icc - 1.96 * se_icc)
  icc_upper <- min(1, icc + 1.96 * se_icc)
  
  results <- rbind(results, data.table(
    Timepoint = paste(tp, "hours"),
    N = n_obs,
    ICC = round(icc, 4),
    `95% CI Lower` = round(icc_lower, 4),
    `95% CI Upper` = round(icc_upper, 4)
  ))
}

print(results)
fwrite(results, "ICC_Results.csv")

