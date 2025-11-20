# Hemodynamic-Response-Reproducibility-in-Repeated-Intubations

## Overview
This repository contains the complete analysis pipeline for a study examining whether hemodynamic responses to intubation represent reproducible, patient-specific traits in critically ill patients. Using a large cohort of patients who underwent multiple intubation events, we evaluated the consistency of blood pressure responses within individuals across time.

**Key Finding:** Hemodynamic responses to intubation demonstrate minimal patient-specific consistency (adjusted ICC = 0.090 at 30 minutes), with consistency rapidly decaying over time. This suggests that hemodynamic responses are not stable individual traits.

## Repository Structure

### Analysis Scripts

The analysis pipeline consists of five sequential R scripts:

1. **`step0_anonymize_data.R`** - Data anonymization
   - Generates unique study IDs for patient protection
   - Creates anonymized versions of all datasets
   - Outputs: `reintubation_all_anonymized.csv`, `reintubation_diff_anonymized.csv`, `reintubation_same_anonymized.csv`

2. **`step1_clean_data.R`** - Data cleaning and preparation
   - Extracts key variables from multiple datasets
   - Standardizes column names across datasets
   - **Calculates Effective MAP** (adjusts raw MAP for vasopressor support)
   - Handles missing ejection fraction data
   - Creates paired cohort with first and second intubations
   - Calculates MAP changes and time intervals between intubations
   - Outputs: `cleaned_all_intubations.csv`, `cleaned_paired_cohort.csv`, `cleaned_datasets.rds`

3. **`step2_descriptiveanalysis.R`** - Statistical analysis
   - Generates patient characteristics table (Table 1)
   - Calculates MAP trajectories by intubation number
   - Computes within-patient correlations
   - **Performs ICC analysis** (both unadjusted and adjusted models)
   - Conducts time-stratified sensitivity analysis
   - Outputs: Multiple CSV files with results and `icc_analysis_results.rds`

4. **`step3_exportdata.R`** - Model coefficient extraction  
   - Extracts fixed effects coefficients with confidence intervals
   - Calculates variance components and ICC
   - Computes ICC confidence intervals using delta method
   - Analyzes ICC across multiple timepoints
   - Outputs: `Table3_Adjusted_Model_Coefficients.csv`, `Table3_Random_Effects.csv`, `ICC_Results.csv`

5. **`step4_visualizations.R`**  
   - Creates figures
   - Outputs: Figure1_Scatterplot, Figure2_BlandAltman, Figure3_Spaghetti, Figure4_ICC_ByTime (PDF and PNG formats)

## Methods Summary

### Study Population
The analysis includes patients who underwent multiple intubation events. Each patient contributed one first intubation and their earliest subsequent reintubation event.

### Effective MAP Calculation
A critical methodological innovation is the calculation of "Effective MAP," which adjusts raw mean arterial pressure for vasopressor support. This allows more accurate comparison of hemodynamic responses independent of pharmacologic blood pressure augmentation. See 'Hemodynamic Effects of Guideline-Based Sedation in Mechanically Ventilated Adults: A Multicenter Observational Analysis': https://pubmed.ncbi.nlm.nih.gov/40911767/. 

**Formula:**
- For patients NOT on vasopressors: Effective MAP = Raw MAP
- For patients on vasopressors (NEE in mcg/kg/min):
  - If NEE ≤ 0.20: Effective MAP = Raw MAP - (NEE / 0.05) × 5
  - If NEE > 0.20: Effective MAP = Raw MAP - ((0.20 / 0.05) × 5 + ((NEE - 0.20) / 0.05) × 3)

This adjustment assumes a 5 mmHg MAP increase per 0.05 mcg/kg/min norepinephrine equivalent up to 0.20 mcg/kg/min, then 3 mmHg per 0.05 mcg/kg/min thereafter.

### Statistical Analysis

**Primary Analysis:** Intraclass correlation coefficient (ICC) estimated using linear mixed-effects models with patient-specific random intercepts. The ICC quantifies the proportion of variance in hemodynamic response attributable to stable patient characteristics.

**Adjusted Model:** Controlled for baseline effective MAP, SOFA score, age, sex, weight, ejection fraction (with missing indicator), hypertension, chronic kidney disease, and congestive heart failure.

**Sensitivity Analyses:**
1. ICC across multiple post-intubation timepoints (30 min, 1 hour, 1.5 hours)
2. Time-stratified analysis by interval between intubations (<24h, 24-72h, >72h)

## Software Requirements

### R Packages
```r
pacman::p_load(
  tidyverse,      # Data manipulation and visualization
  data.table,     # Efficient data handling
  lme4,           # Mixed-effects models
  lmerTest,       # P-values for mixed models
  performance,    # ICC calculation
  broom.mixed,    # Model tidying
  ggplot2,        # Plotting
  patchwork       # Figure composition
)
```

### R Version
- R 4.4.2 or later recommended
