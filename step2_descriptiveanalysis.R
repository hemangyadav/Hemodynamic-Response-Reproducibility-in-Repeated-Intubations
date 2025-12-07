setwd("_____")

library(pacman)
p_load(tidyverse, data.table, lme4, performance)

#1 Load Cleaned Data -----------
cat("\n1. Loading cleaned data...\n")
datasets <- readRDS("cleaned_datasets.rds")
all_intubations <- datasets$all_intubations
paired_cohort <- datasets$paired_cohort

#--- Keep only first reintubation for patients with 3+ intubations ---

# Count intubations per patient
n_intubs <- paired_cohort[, .N, by = .(StudyID, IntubationNumber)]
cat(sprintf("Patients with 3+ intubations: %d\n", 
            n_intubs[N > 1, uniqueN(StudyID)]))

# Convert datetime columns
paired_cohort[, IntubationDatetime := as.POSIXct(IntubationDatetime)]
paired_cohort[, FIRST_INTUBATION := as.POSIXct(FIRST_INTUBATION)]

# Identify first reintubation for each patient (by actual intubation datetime)
first_reintub <- paired_cohort[IntubationNumber == 2, 
                               .SD[which.min(IntubationDatetime)],  
                               by = StudyID]

# Keep all IntubationNumber=1 plus only first reintubation
paired_cohort_clean <- rbind(
  paired_cohort[IntubationNumber == 1],
  first_reintub
)

cat(sprintf("Original events: %d\n", nrow(paired_cohort)))
cat(sprintf("Cleaned events: %d (kept first reintubation only)\n", nrow(paired_cohort_clean)))
cat(sprintf("Unique patients: %d\n", uniqueN(paired_cohort_clean$StudyID)))

# Replace paired_cohort with cleaned version
paired_cohort <- paired_cohort_clean

#2 Table 1: Patient Characteristics -----------
baseline_chars <- paired_cohort[IntubationNumber == 1, .(
  StudyID, AGE_AT_ADMISSION, PATIENT_SEX, HOSPITAL_ADMISSION_WEIGHT_KG,
  SOFA, EJECTION_FRACTION_THREE_MONTHS, ef_missing,
  HTN, CKD, CHF, MI
)]

table1 <- data.table(
  Variable = c("Age (years)", "Female sex", "Weight (kg)", "SOFA score", 
               "Ejection fraction (%)", "EF missing", 
               "Hypertension", "Chronic kidney disease", 
               "Congestive heart failure", "Myocardial infarction"),
  
  Value = c(
    sprintf("%.0f [%.0f-%.0f]", 
            median(baseline_chars$AGE_AT_ADMISSION, na.rm = TRUE),
            quantile(baseline_chars$AGE_AT_ADMISSION, 0.25, na.rm = TRUE),
            quantile(baseline_chars$AGE_AT_ADMISSION, 0.75, na.rm = TRUE)),
    
    sprintf("%d (%.1f%%)", 
            sum(baseline_chars$PATIENT_SEX == "Female", na.rm = TRUE),
            100 * mean(baseline_chars$PATIENT_SEX == "Female", na.rm = TRUE)),
    
    sprintf("%.1f [%.1f-%.1f]", 
            median(baseline_chars$HOSPITAL_ADMISSION_WEIGHT_KG, na.rm = TRUE),
            quantile(baseline_chars$HOSPITAL_ADMISSION_WEIGHT_KG, 0.25, na.rm = TRUE),
            quantile(baseline_chars$HOSPITAL_ADMISSION_WEIGHT_KG, 0.75, na.rm = TRUE)),
    
    sprintf("%.0f [%.0f-%.0f]", 
            median(baseline_chars$SOFA, na.rm = TRUE),
            quantile(baseline_chars$SOFA, 0.25, na.rm = TRUE),
            quantile(baseline_chars$SOFA, 0.75, na.rm = TRUE)),
    
    sprintf("%.0f [%.0f-%.0f]", 
            median(baseline_chars$EJECTION_FRACTION_THREE_MONTHS, na.rm = TRUE),
            quantile(baseline_chars$EJECTION_FRACTION_THREE_MONTHS, 0.25, na.rm = TRUE),
            quantile(baseline_chars$EJECTION_FRACTION_THREE_MONTHS, 0.75, na.rm = TRUE)),
    
    sprintf("%d (%.1f%%)", 
            sum(baseline_chars$ef_missing, na.rm = TRUE),
            100 * mean(baseline_chars$ef_missing, na.rm = TRUE)),
    
    sprintf("%d (%.1f%%)", 
            sum(baseline_chars$HTN == TRUE, na.rm = TRUE),
            100 * mean(baseline_chars$HTN == TRUE, na.rm = TRUE)),
    
    sprintf("%d (%.1f%%)", 
            sum(baseline_chars$CKD == TRUE, na.rm = TRUE),
            100 * mean(baseline_chars$CKD == TRUE, na.rm = TRUE)),
    
    sprintf("%d (%.1f%%)", 
            sum(baseline_chars$CHF == TRUE, na.rm = TRUE),
            100 * mean(baseline_chars$CHF == TRUE, na.rm = TRUE)),
    
    sprintf("%d (%.1f%%)", 
            sum(baseline_chars$MI == TRUE, na.rm = TRUE),
            100 * mean(baseline_chars$MI == TRUE, na.rm = TRUE))
  )
)

print(table1)
fwrite(table1, "Table1_PatientCharacteristics.csv")

#3 MAP Trajectory -----------
map_trajectory <- paired_cohort[, .(
  N = .N,
  Baseline_MAP = median(MAP_bin_0, na.rm = TRUE),
  MAP_0to05h = median(MAP_bin_0.5, na.rm = TRUE),
  MAP_05to1h = median(MAP_bin_1, na.rm = TRUE),
  MAP_1to15h = median(MAP_bin_1.5, na.rm = TRUE),
  Baseline_EffMAP = median(EffectiveMAP_bin_0, na.rm = TRUE),
  EffMAP_0to05h = median(EffectiveMAP_bin_0.5, na.rm = TRUE),
  EffMAP_05to1h = median(EffectiveMAP_bin_1, na.rm = TRUE),
  EffMAP_1to15h = median(EffectiveMAP_bin_1.5, na.rm = TRUE)
), by = IntubationNumber]

print(map_trajectory)
fwrite(map_trajectory, "MAP_Trajectory_ByIntubation.csv")

#4 Within-Patient Correlation -----------

# Create wide format
wide_data <- dcast(paired_cohort, 
                   StudyID + AGE_AT_ADMISSION + PATIENT_SEX + SOFA + 
                     EJECTION_FRACTION_THREE_MONTHS + ef_missing + 
                     HTN + CKD + CHF ~ IntubationNumber,
                   value.var = c("MAP_bin_0", "MAP_bin_0.5", "MAP_bin_1", "MAP_bin_1.5",
                                 "EffectiveMAP_bin_0", "EffectiveMAP_bin_0.5", 
                                 "EffectiveMAP_bin_1", "EffectiveMAP_bin_1.5",
                                 "delta_EffMAP_0to05", "delta_EffMAP_0to10", "delta_EffMAP_0to15"))

correlations <- data.table(
  TimePoint = c("Baseline", "0-0.5h", "0.5-1h", "1-1.5h", 
                "Delta 0-0.5h", "Delta 0-1h", "Delta 0-1.5h"),
  
  Raw_MAP_Correlation = c(
    cor(wide_data$MAP_bin_0_1, wide_data$MAP_bin_0_2, use = "complete.obs"),
    cor(wide_data$MAP_bin_0.5_1, wide_data$MAP_bin_0.5_2, use = "complete.obs"),
    cor(wide_data$MAP_bin_1_1, wide_data$MAP_bin_1_2, use = "complete.obs"),
    cor(wide_data$MAP_bin_1.5_1, wide_data$MAP_bin_1.5_2, use = "complete.obs"),
    NA, NA, NA
  ),
  
  Effective_MAP_Correlation = c(
    cor(wide_data$EffectiveMAP_bin_0_1, wide_data$EffectiveMAP_bin_0_2, use = "complete.obs"),
    cor(wide_data$EffectiveMAP_bin_0.5_1, wide_data$EffectiveMAP_bin_0.5_2, use = "complete.obs"),
    cor(wide_data$EffectiveMAP_bin_1_1, wide_data$EffectiveMAP_bin_1_2, use = "complete.obs"),
    cor(wide_data$EffectiveMAP_bin_1.5_1, wide_data$EffectiveMAP_bin_1.5_2, use = "complete.obs"),
    cor(wide_data$delta_EffMAP_0to05_1, wide_data$delta_EffMAP_0to05_2, use = "complete.obs"),
    cor(wide_data$delta_EffMAP_0to10_1, wide_data$delta_EffMAP_0to10_2, use = "complete.obs"),
    cor(wide_data$delta_EffMAP_0to15_1, wide_data$delta_EffMAP_0to15_2, use = "complete.obs")
  )
)

print(correlations)
fwrite(correlations, "Within_Patient_Correlations.csv")

#5 ICC Analysis - Unadjusted -----------
# Raw MAP
icc_raw_unadj <- lmer(MAP_bin_0.5 ~ 1 + (1 | StudyID), data = paired_cohort)
icc_raw_value <- performance::icc(icc_raw_unadj)

# Effective MAP
icc_eff_unadj <- lmer(EffectiveMAP_bin_0.5 ~ 1 + (1 | StudyID), data = paired_cohort)
icc_eff_value <- performance::icc(icc_eff_unadj)

cat(sprintf("\n   Raw MAP ICC: %.3f\n", icc_raw_value$ICC_adjusted))
cat(sprintf("   Effective MAP ICC: %.3f\n", icc_eff_value$ICC_adjusted))

#6 ICC Analysis - Adjusted -----------
# Raw MAP - adjusted
icc_raw_adj <- lmer(MAP_bin_0.5 ~ 
                      Baseline_MAP +
                      SOFA + AGE_AT_ADMISSION + PATIENT_SEX + 
                      HOSPITAL_ADMISSION_WEIGHT_KG +
                      HTN + CKD + CHF +
                      (1 | StudyID),
                    data = paired_cohort)

icc_raw_adj_value <- performance::icc(icc_raw_adj)

# Effective MAP - adjusted
icc_eff_adj <- lmer(EffectiveMAP_bin_0.5 ~ 
                      Baseline_EffMAP +
                      SOFA + AGE_AT_ADMISSION + PATIENT_SEX + 
                      HOSPITAL_ADMISSION_WEIGHT_KG +
                      HTN + CKD + CHF +
                      (1 | StudyID),
                    data = paired_cohort)

icc_eff_adj_value <- performance::icc(icc_eff_adj)

cat(sprintf("\n   Adjusted Raw MAP ICC: %.3f\n", icc_raw_adj_value$ICC_adjusted))
cat(sprintf("   Adjusted Effective MAP ICC: %.3f\n", icc_eff_adj_value$ICC_adjusted))

#7 ICC Comparison Table -----------
icc_comparison <- data.table(
  Model = c("Unadjusted - Raw MAP", 
            "Unadjusted - Effective MAP",
            "Adjusted - Raw MAP",
            "Adjusted - Effective MAP"),
  
  ICC = c(
    icc_raw_value$ICC_adjusted,
    icc_eff_value$ICC_adjusted,
    icc_raw_adj_value$ICC_adjusted,
    icc_eff_adj_value$ICC_adjusted
  ),
  
  Interpretation = c(
    ifelse(icc_raw_value$ICC_adjusted > 0.5, "High", 
           ifelse(icc_raw_value$ICC_adjusted > 0.3, "Moderate", "Low")),
    ifelse(icc_eff_value$ICC_adjusted > 0.5, "High", 
           ifelse(icc_eff_value$ICC_adjusted > 0.3, "Moderate", "Low")),
    ifelse(icc_raw_adj_value$ICC_adjusted > 0.5, "High", 
           ifelse(icc_raw_adj_value$ICC_adjusted > 0.3, "Moderate", "Low")),
    ifelse(icc_eff_adj_value$ICC_adjusted > 0.5, "High", 
           ifelse(icc_eff_adj_value$ICC_adjusted > 0.3, "Moderate", "Low"))
  )
)

print(icc_comparison)
fwrite(icc_comparison, "ICC_Comparison.csv")

#8 ICC by Time Stratum -----------

if("HoursBetweenIntubations" %in% names(paired_cohort)) {
  paired_cohort[, TimeBetween_Strata := cut(HoursBetweenIntubations,
                                            breaks = c(0, 24, 72, Inf),
                                            labels = c("<24h", "24-72h", ">72h"))]
  
  time_strata_icc <- data.table()
  
  for(strat in c("<24h", "24-72h", ">72h")) {
    subset_data <- paired_cohort[TimeBetween_Strata == strat]
    
    if(nrow(subset_data) >= 20 && uniqueN(subset_data$StudyID) >= 10) {
      
      tryCatch({
        model <- lmer(EffectiveMAP_bin_0.5 ~ 
                        Baseline_EffMAP +
                        SOFA + AGE_AT_ADMISSION + PATIENT_SEX + 
                        HOSPITAL_ADMISSION_WEIGHT_KG +
                        HTN + CKD + CHF +
                        (1 | StudyID),
                      data = subset_data)
        
        icc_val <- performance::icc(model)
        
        time_strata_icc <- rbind(time_strata_icc, data.table(
          Stratum = strat,
          N_Patients = uniqueN(subset_data$StudyID),
          N_Events = nrow(subset_data),
          ICC = icc_val$ICC_adjusted
        ))
      }, error = function(e) {
        cat(sprintf("   Warning: Could not fit model for %s stratum\n", strat))
      })
    } else {
      cat(sprintf("   Skipping %s (insufficient data)\n", strat))
    }
  }
  
  if(nrow(time_strata_icc) > 0) {
    print(time_strata_icc)
    fwrite(time_strata_icc, "ICC_By_Time_Stratum.csv")
    cat("   Saved: ICC_By_Time_Stratum.csv\n")
  }
}

#--- SENSITIVITY: ICC across all timepoints ------
timepoint_icc <- data.table()

for(tp in c("0.5", "1", "1.5")) {
  outcome_col <- paste0("EffectiveMAP_bin_", tp)
  
  # Skip if too much missing data
  pct_complete <- 100 * mean(!is.na(paired_cohort[[outcome_col]]))
  if(pct_complete < 50) {
    cat(sprintf("Skipping %s hours (%.1f%% complete)\n", tp, pct_complete))
    next
  }
  
  cat(sprintf("\nAnalyzing timepoint: %s hours (%.1f%% complete)\n", tp, pct_complete))
  
  # Unadjusted ICC
  tryCatch({
    formula_unadj <- as.formula(paste0(outcome_col, " ~ Baseline_EffMAP + (1 | StudyID)"))
    model_unadj <- lmer(formula_unadj, data = paired_cohort, REML = TRUE)
    icc_unadj <- performance::icc(model_unadj)
    
    # Adjusted ICC
    formula_adj <- as.formula(paste0(outcome_col, " ~ Baseline_EffMAP + SOFA + AGE_AT_ADMISSION + ",
                                     "PATIENT_SEX + HOSPITAL_ADMISSION_WEIGHT_KG + ",
                                     "HTN + CKD + CHF + (1 | StudyID)"))
    model_adj <- lmer(formula_adj, data = paired_cohort, REML = TRUE)
    icc_adj <- performance::icc(model_adj)
    
    timepoint_icc <- rbind(timepoint_icc, data.table(
      Timepoint = paste0(tp, " hours"),
      N = sum(complete.cases(paired_cohort[, c(outcome_col, "Baseline_EffMAP", "SOFA"), with = FALSE])),
      Unadjusted_ICC = icc_unadj$ICC_adjusted,
      Adjusted_ICC = icc_adj$ICC_adjusted
    ))
  }, error = function(e) {
    cat(sprintf("  Error fitting model for %s hours: %s\n", tp, e$message))
  })
}

if(nrow(timepoint_icc) > 0) {
  print(timepoint_icc)
  fwrite(timepoint_icc, "ICC_By_Timepoint.csv")
}

#9 Model Summaries -----------
print(summary(icc_eff_adj))

#10 Save All Results -----------
icc_results <- list(
  unadjusted = list(
    raw_map = icc_raw_unadj,
    effective_map = icc_eff_unadj,
    raw_icc = icc_raw_value,
    eff_icc = icc_eff_value
  ),
  adjusted = list(
    raw_map = icc_raw_adj,
    effective_map = icc_eff_adj,
    raw_icc = icc_raw_adj_value,
    eff_icc = icc_eff_adj_value
  ),
  comparison = icc_comparison,
  wide_data = wide_data,
  correlations = correlations
)

saveRDS(icc_results, "icc_analysis_results.rds")

cat(sprintf("  - Unadjusted Raw MAP ICC: %.3f\n", icc_raw_value$ICC_adjusted))
cat(sprintf("  - Unadjusted Effective MAP ICC: %.3f\n", icc_eff_value$ICC_adjusted))
cat(sprintf("  - Adjusted Raw MAP ICC: %.3f\n", icc_raw_adj_value$ICC_adjusted))
cat(sprintf("  - Adjusted Effective MAP ICC: %.3f\n", icc_eff_adj_value$ICC_adjusted))
