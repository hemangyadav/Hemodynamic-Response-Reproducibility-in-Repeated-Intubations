setwd("_____")

library(pacman)
p_load(tidyverse, data.table, lme4)

#1 Load Data -----------
reintub_all <- fread('reintubation_all_anonymized.csv')
reintub_diff <- fread('reintubation_diff_anonymized.csv')
reintub_same <- fread('reintubation_same_anonymized.csv')

#2 Extract Key Variables -----------
# Core vars that should be in ALL datasets
core_vars_common <- c(
  "StudyID", "SOFA", "AGE_AT_ADMISSION", "PATIENT_SEX", 
  "HOSPITAL_ADMISSION_WEIGHT_KG", "EJECTION_FRACTION_THREE_MONTHS"
)

# Variables specific to first intubation dataset
first_specific <- c("ICU_UNIT", "ICU_START_DATETIME", "FIRST_INTUBATION", 
                    "FIRST_PROP_TIME", "IS_REINTUBATION", "SAME_ICU")

# Comorbidity vars ONLY in reintub_all (patient-level characteristics)
comorbid_vars <- c("HTN", "CKD", "CHF", "MI")

# Get MAP bins
map_cols_all <- c("MAP_firsticu_bin_0", "MAP_firsticu_bin_0.5", 
                  "MAP_firsticu_bin_1", "MAP_firsticu_bin_1.5")
map_cols_diff <- c("MAP_difficu_bin_0", "MAP_difficu_bin_0.5", 
                   "MAP_difficu_bin_1", "MAP_difficu_bin_1.5")

# Get NE equivalents (per kg per min versions)
ne_cols <- c("NE_EQUIVALENT_-0.5_perkg_permin", "NE_EQUIVALENT_0.0_perkg_permin",
             "NE_EQUIVALENT_0.5_perkg_permin", "NE_EQUIVALENT_1.0_perkg_permin",
             "NE_EQUIVALENT_1.5_perkg_permin", "NE_EQUIVALENT_2.0_perkg_permin")

# Extract from first intubation (has all variables including comorbidities)
first_intub <- reintub_all[, c(core_vars_common, first_specific, comorbid_vars, 
                               map_cols_all, ne_cols), with = FALSE]
first_intub[, IntubationNumber := 1]
first_intub[, IntubationDatetime := as.POSIXct(FIRST_INTUBATION)]

# Extract patient-level info to merge with reintubations
patient_info <- first_intub[, c("StudyID", "ICU_UNIT", "ICU_START_DATETIME", 
                                "FIRST_INTUBATION", "FIRST_PROP_TIME", 
                                "IS_REINTUBATION", "SAME_ICU", comorbid_vars), with = FALSE]

# Extract from diff ICU reintubations
if(nrow(reintub_diff) > 0) {
  cols_to_extract_diff <- intersect(c(core_vars_common, map_cols_diff, ne_cols, 
                                      "VENT2_DIFFICU_START_DATETIME"), 
                                    names(reintub_diff))
  second_intub_diff <- reintub_diff[, cols_to_extract_diff, with = FALSE]
  second_intub_diff[, IntubationNumber := 2]
  second_intub_diff[, IntubationDatetime := as.POSIXct(VENT2_DIFFICU_START_DATETIME)]
  second_intub_diff[, VENT2_DIFFICU_START_DATETIME := NULL]
  # Merge patient-level info from first intubation
  second_intub_diff <- merge(second_intub_diff, patient_info, by = "StudyID", all.x = TRUE)
} else {
  second_intub_diff <- data.table()
}

# Extract from same ICU reintubations
if(nrow(reintub_same) > 0) {
  # Check which MAP columns exist
  map_cols_same_check <- intersect(c("MAP_sameicu_bin_0", "MAP_sameicu_bin_0.5", 
                                     "MAP_sameicu_bin_1", "MAP_sameicu_bin_1.5",
                                     "MAP_firsticu_bin_0", "MAP_firsticu_bin_0.5", 
                                     "MAP_firsticu_bin_1", "MAP_firsticu_bin_1.5"),
                                   names(reintub_same))
  
  cols_to_extract_same <- intersect(c(core_vars_common, map_cols_same_check, ne_cols,
                                      "SAMEREINTUB_VENT_START_DATETIME"), 
                                    names(reintub_same))
  second_intub_same <- reintub_same[, cols_to_extract_same, with = FALSE]
  second_intub_same[, IntubationNumber := 2]
  second_intub_same[, IntubationDatetime := as.POSIXct(SAMEREINTUB_VENT_START_DATETIME)]
  second_intub_same[, SAMEREINTUB_VENT_START_DATETIME := NULL]
  # Merge patient-level info from first intubation
  second_intub_same <- merge(second_intub_same, patient_info, by = "StudyID", all.x = TRUE)
} else {
  second_intub_same <- data.table()
}

#3 Standardize MAP Column Names -----------
setnames(first_intub, 
         old = c("MAP_firsticu_bin_0", "MAP_firsticu_bin_0.5", 
                 "MAP_firsticu_bin_1", "MAP_firsticu_bin_1.5"),
         new = c("MAP_bin_0", "MAP_bin_0.5", "MAP_bin_1", "MAP_bin_1.5"))

if(nrow(second_intub_diff) > 0) {
  setnames(second_intub_diff,
           old = c("MAP_difficu_bin_0", "MAP_difficu_bin_0.5", 
                   "MAP_difficu_bin_1", "MAP_difficu_bin_1.5"),
           new = c("MAP_bin_0", "MAP_bin_0.5", "MAP_bin_1", "MAP_bin_1.5"))
}

if(nrow(second_intub_same) > 0) {
  # Rename MAP columns
  old_map_names <- grep("MAP.*bin", names(second_intub_same), value = TRUE)
  new_map_names <- gsub("MAP_.*icu_bin_", "MAP_bin_", old_map_names)
  setnames(second_intub_same, old = old_map_names, new = new_map_names)
}

# Rename NE columns to shorter names
ne_old <- c("NE_EQUIVALENT_-0.5_perkg_permin", "NE_EQUIVALENT_0.0_perkg_permin",
            "NE_EQUIVALENT_0.5_perkg_permin", "NE_EQUIVALENT_1.0_perkg_permin",
            "NE_EQUIVALENT_1.5_perkg_permin", "NE_EQUIVALENT_2.0_perkg_permin")
ne_new <- c("NE_m05", "NE_0", "NE_05", "NE_10", "NE_15", "NE_20")

for(dt in list(first_intub, second_intub_diff, second_intub_same)) {
  if(nrow(dt) > 0) {
    setnames(dt, old = ne_old, new = ne_new, skip_absent = TRUE)
  }
}

#4 Calculate Effective MAP -----------
calculate_effective_MAP <- function(raw_MAP, NE_mcg_kg_min) {
  # If raw MAP is missing, return NA
  if_else(is.na(raw_MAP), 
          NA_real_,
          # If NE is missing or zero, no adjustment needed - use raw MAP
          if_else(is.na(NE_mcg_kg_min) | NE_mcg_kg_min == 0,
                  raw_MAP,
                  # Otherwise apply NE adjustment
                  if_else(NE_mcg_kg_min <= 0.20,
                          raw_MAP - (NE_mcg_kg_min / 0.05) * 5,
                          raw_MAP - ((0.20 / 0.05) * 5 + ((NE_mcg_kg_min - 0.20) / 0.05) * 3))))
}

# Map time bins to NE columns
time_map <- list(
  "0" = "NE_0",
  "0.5" = "NE_05",
  "1" = "NE_10",
  "1.5" = "NE_15"
)

for(dt_name in c("first_intub", "second_intub_diff", "second_intub_same")) {
  dt <- get(dt_name)
  if(nrow(dt) == 0) next
  
  for(time_bin in names(time_map)) {
    map_col <- paste0("MAP_bin_", time_bin)
    ne_col <- time_map[[time_bin]]
    eff_col <- paste0("EffectiveMAP_bin_", time_bin)
    
    if(map_col %in% names(dt) && ne_col %in% names(dt)) {
      dt[, (eff_col) := calculate_effective_MAP(get(map_col), get(ne_col))]
    }
  }
  
  assign(dt_name, dt)
}

#5 Handle Missing EF -----------
for(dt_name in c("first_intub", "second_intub_diff", "second_intub_same")) {
  dt <- get(dt_name)
  if(nrow(dt) == 0) next
  
  dt[, ef_missing := is.na(EJECTION_FRACTION_THREE_MONTHS)]
  assign(dt_name, dt)
}

median_ef <- median(first_intub$EJECTION_FRACTION_THREE_MONTHS, na.rm = TRUE)
cat(sprintf("   Median EF: %.1f%%\n", median_ef))

first_intub[is.na(EJECTION_FRACTION_THREE_MONTHS), 
            EJECTION_FRACTION_THREE_MONTHS := median_ef]
if(nrow(second_intub_diff) > 0) {
  second_intub_diff[is.na(EJECTION_FRACTION_THREE_MONTHS), 
                    EJECTION_FRACTION_THREE_MONTHS := median_ef]
}
if(nrow(second_intub_same) > 0) {
  second_intub_same[is.na(EJECTION_FRACTION_THREE_MONTHS), 
                    EJECTION_FRACTION_THREE_MONTHS := median_ef]
}

#6 Combine All Intubations -----------
all_intubations <- rbindlist(list(first_intub, second_intub_diff, second_intub_same), 
                             fill = TRUE)
setorder(all_intubations, StudyID, IntubationNumber)

cat(sprintf("   Total events: %d from %d patients\n", 
            nrow(all_intubations), uniqueN(all_intubations$StudyID)))

#7 Create Paired Cohort -----------
paired_patients <- all_intubations[, .N, by = StudyID][N >= 2, StudyID]
paired_cohort <- all_intubations[StudyID %in% paired_patients]

cat(sprintf("   Paired cohort: %d events from %d patients\n", 
            nrow(paired_cohort), length(paired_patients)))

#8 Calculate MAP Changes -----------
for(dt_name in c("all_intubations", "paired_cohort")) {
  dt <- get(dt_name)
  
  # Raw MAP changes
  dt[, `:=`(
    delta_MAP_0to05 = MAP_bin_0.5 - MAP_bin_0,
    delta_MAP_0to10 = MAP_bin_1 - MAP_bin_0,
    delta_MAP_0to15 = MAP_bin_1.5 - MAP_bin_0,
    pct_MAP_0to05 = 100 * (MAP_bin_0.5 - MAP_bin_0) / MAP_bin_0,
    pct_MAP_0to10 = 100 * (MAP_bin_1 - MAP_bin_0) / MAP_bin_0,
    pct_MAP_0to15 = 100 * (MAP_bin_1.5 - MAP_bin_0) / MAP_bin_0
  )]
  
  # Effective MAP changes
  dt[, `:=`(
    delta_EffMAP_0to05 = EffectiveMAP_bin_0.5 - EffectiveMAP_bin_0,
    delta_EffMAP_0to10 = EffectiveMAP_bin_1 - EffectiveMAP_bin_0,
    delta_EffMAP_0to15 = EffectiveMAP_bin_1.5 - EffectiveMAP_bin_0,
    pct_EffMAP_0to05 = 100 * (EffectiveMAP_bin_0.5 - EffectiveMAP_bin_0) / EffectiveMAP_bin_0,
    pct_EffMAP_0to10 = 100 * (EffectiveMAP_bin_1 - EffectiveMAP_bin_0) / EffectiveMAP_bin_0,
    pct_EffMAP_0to15 = 100 * (EffectiveMAP_bin_1.5 - EffectiveMAP_bin_0) / EffectiveMAP_bin_0
  )]
  
  # Add baseline variables
  dt[, Baseline_MAP := MAP_bin_0]
  dt[, Baseline_EffMAP := EffectiveMAP_bin_0]
  
  assign(dt_name, dt)
}

#9 Calculate Time Between Intubations -----------
if(nrow(paired_cohort) > 0) {

  intub_times <- paired_cohort[, .(
    FirstIntubTime = IntubationDatetime[IntubationNumber == 1][1],
    SecondIntubTime = IntubationDatetime[IntubationNumber == 2][1]
  ), by = StudyID]
  
  intub_times[, HoursBetweenIntubations := 
                as.numeric(difftime(SecondIntubTime, FirstIntubTime, units = "hours"))]
  
  # Remove if exists and merge back
  if("HoursBetweenIntubations" %in% names(paired_cohort)) {
    paired_cohort[, HoursBetweenIntubations := NULL]
  }
  
  paired_cohort <- merge(paired_cohort, 
                         intub_times[, .(StudyID, HoursBetweenIntubations)], 
                         by = "StudyID", all.x = TRUE)
  
  cat(sprintf("   Median time between: %.1f hours (%.1f days)\n", 
              median(intub_times$HoursBetweenIntubations, na.rm = TRUE),
              median(intub_times$HoursBetweenIntubations, na.rm = TRUE) / 24))
}

#10 Save Cleaned Datasets -----------
fwrite(all_intubations, "cleaned_all_intubations.csv")
fwrite(paired_cohort, "cleaned_paired_cohort.csv")

saveRDS(list(
  all_intubations = all_intubations,
  paired_cohort = paired_cohort
), "cleaned_datasets.rds")

