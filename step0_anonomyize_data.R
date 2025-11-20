setwd("_____")

library(pacman)
p_load(tidyverse,
       data.table,
       lmer)

#1 Data Loading -----------
reintub_all <- read_csv('Reintubation_ALL_Jul17rev_ComorbiditiesDoseAdjust.csv')
reintub_diff <- read_csv('Reintubation_ALLDifferent_Jul17_DoseAdjust.csv')
reintub_same <- read_csv('Reintubation_ALLSame_Jul17_DoseAdjust.csv')


#2 Anonymization -----------
unique_patients <- unique(reintub_all$CLINIC_NUMBER)
anonymization_key <- data.table(
  CLINIC_NUMBER = unique_patients,
  StudyID = paste0("PT", sprintf("%05d", seq_along(unique_patients)))
)

reintub_all_anon <- reintub_all %>%
  left_join(anonymization_key, by = "CLINIC_NUMBER") %>%
  select(-CLINIC_NUMBER) %>%
  select(StudyID, everything())

reintub_diff_anon <- reintub_diff %>%
  left_join(anonymization_key, by = "CLINIC_NUMBER") %>%
  select(-CLINIC_NUMBER) %>%
  select(StudyID, everything())

reintub_same_anon <- reintub_same %>%
  left_join(anonymization_key, by = "CLINIC_NUMBER") %>%
  select(-CLINIC_NUMBER) %>%
  select(StudyID, everything())

#3 Export Anonymized Data -----------
write_csv(anonymization_key, 'anonymization_key.csv')
write_csv(reintub_all_anon, 'reintubation_all_anonymized.csv')
write_csv(reintub_diff_anon, 'reintubation_diff_anonymized.csv')
write_csv(reintub_same_anon, 'reintubation_same_anonymized.csv')

#4 Display Structure -----------
cat(sprintf("Total unique patients: %d\n", length(unique_patients)))
cat(sprintf("All reintubations: %d rows\n", nrow(reintub_all_anon)))
cat(sprintf("Different ICU reintubations: %d rows\n", nrow(reintub_diff_anon)))
cat(sprintf("Same ICU reintubations: %d rows\n\n", nrow(reintub_same_anon)))

cat("--- Reintubation All ---\n")
str(reintub_all_anon)
cat("\n--- Reintubation Different ICU ---\n")
str(reintub_diff_anon)
cat("\n--- Reintubation Same ICU ---\n")
str(reintub_same_anon)


