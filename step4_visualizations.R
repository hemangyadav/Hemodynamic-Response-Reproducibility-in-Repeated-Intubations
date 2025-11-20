library(pacman)
p_load(tidyverse, data.table, ggplot2, patchwork)

#1 Load Data -----------
cat("\n1. Loading cleaned data...\n")
data <- readRDS("cleaned_datasets.rds")
paired_cohort <- data$paired_cohort

#2 Prepare Wide Format -----------
paired_cohort[, FIRST_INTUBATION := as.POSIXct(FIRST_INTUBATION)]
first_reintub <- paired_cohort[IntubationNumber == 2, 
                               .SD[which.min(FIRST_INTUBATION)], 
                               by = StudyID]
paired_clean <- rbind(
  paired_cohort[IntubationNumber == 1],
  first_reintub
)

wide_data <- dcast(paired_clean, 
                   StudyID ~ IntubationNumber,
                   value.var = c("EffectiveMAP_bin_0", "EffectiveMAP_bin_0.5", 
                                 "EffectiveMAP_bin_1", "EffectiveMAP_bin_1.5",
                                 "delta_EffMAP_0to05"))

# Remove patients with missing data
wide_data <- wide_data[complete.cases(wide_data[, .(EffectiveMAP_bin_0.5_1, EffectiveMAP_bin_0.5_2)])]

#3 Figure 1: Scatterplot -----------
fig1 <- ggplot(wide_data, aes(x = EffectiveMAP_bin_0.5_1, y = EffectiveMAP_bin_0.5_2)) +
  geom_point(alpha = 0.3, size = 1.5, color = "steelblue") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "darkblue", linewidth = 1) +
  labs(
    title = "Hemodynamic Response Reproducibility Across Intubation Events",
    subtitle = sprintf("N = %d patients | r = %.3f", 
                       nrow(wide_data),
                       cor(wide_data$EffectiveMAP_bin_0.5_1, wide_data$EffectiveMAP_bin_0.5_2, use = "complete.obs")),
    x = "Effective MAP at First Intubation (0-30 min, mmHg)",
    y = "Effective MAP at Second Intubation (0-30 min, mmHg)"
  ) +
  coord_fixed(xlim = c(40, 110), ylim = c(40, 110)) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray30")
  )

ggsave("Figure1_Scatterplot.pdf", fig1, width = 8, height = 7)
ggsave("Figure1_Scatterplot.png", fig1, width = 8, height = 7, dpi = 300)

#4 Figure 2: Bland-Altman Plot -----------

wide_data[, `:=`(
  mean_MAP = (EffectiveMAP_bin_0.5_1 + EffectiveMAP_bin_0.5_2) / 2,
  diff_MAP = EffectiveMAP_bin_0.5_2 - EffectiveMAP_bin_0.5_1
)]

mean_diff <- mean(wide_data$diff_MAP, na.rm = TRUE)
sd_diff <- sd(wide_data$diff_MAP, na.rm = TRUE)
loa_upper <- mean_diff + 1.96 * sd_diff
loa_lower <- mean_diff - 1.96 * sd_diff

fig2 <- ggplot(wide_data, aes(x = mean_MAP, y = diff_MAP)) +
  geom_jitter(alpha = 0.3, size = 1.5, color = "steelblue") +
  geom_hline(yintercept = mean_diff, linetype = "solid", color = "darkblue", linewidth = 1) +
  geom_hline(yintercept = loa_upper, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = loa_lower, linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "gray50") +
  annotate("text", x = 130, y = mean_diff + 2, 
           label = sprintf("Mean: %.1f mmHg", mean_diff), hjust = 1) +
  annotate("text", x = 130, y = loa_upper + 2, 
           label = sprintf("+1.96 SD: %.1f mmHg", loa_upper), hjust = 1) +
  annotate("text", x = 130, y = loa_lower - 2, 
           label = sprintf("-1.96 SD: %.1f mmHg", loa_lower), hjust = 1) +
  labs(
    title = "Bland-Altman Plot: Agreement Between Intubation Events",
    subtitle = sprintf("N = %d patients", nrow(wide_data)),
    x = "Mean Effective MAP (mmHg)",
    y = "Difference in Effective MAP (Second - First, mmHg)"
  ) +
  theme_classic(base_size = 12) +
  ylim(-60, 60) + 
  xlim(40,140) + 
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray30")
  )

ggsave("Figure2_BlandAltman.pdf", fig2, width = 8, height = 8)
ggsave("Figure2_BlandAltman.png", fig2, width = 8, height = 8, dpi = 300)

#5 Figure 3: Spaghetti Plot -----------

# Sample 50 random patients for visualization
set.seed(123)
sample_ids <- sample(unique(paired_clean$StudyID), min(50, uniqueN(paired_clean$StudyID)))
sample_data <- paired_clean[StudyID %in% sample_ids]

# Reshape for plotting
spaghetti_data <- melt(sample_data, 
                       id.vars = c("StudyID", "IntubationNumber"),
                       measure.vars = c("EffectiveMAP_bin_0", "EffectiveMAP_bin_0.5", 
                                        "EffectiveMAP_bin_1", "EffectiveMAP_bin_1.5"),
                       variable.name = "Timepoint",
                       value.name = "EffectiveMAP")

spaghetti_data[, Timepoint := factor(Timepoint, 
                                     levels = c("EffectiveMAP_bin_0", "EffectiveMAP_bin_0.5",
                                                "EffectiveMAP_bin_1", "EffectiveMAP_bin_1.5"),
                                     labels = c("0", "0.5", "1", "1.5"))]
spaghetti_data[, Timepoint := as.numeric(as.character(Timepoint))]

fig3 <- ggplot(spaghetti_data, aes(x = Timepoint, y = EffectiveMAP, group = StudyID)) +
  geom_line(aes(color = factor(IntubationNumber)), alpha = 0.3, linewidth = 0.5) +
  stat_summary(aes(group = IntubationNumber, color = factor(IntubationNumber)),
               fun = median, geom = "line", linewidth = 2) +
  scale_color_manual(
    name = "Intubation",
    values = c("1" = "steelblue", "2" = "coral"),
    labels = c("1" = "First", "2" = "Second")
  ) +
  labs(
    title = "Individual Patient Trajectories Across Intubation Events",
    subtitle = sprintf("Sample of %d patients (thin lines) with median trajectory (thick lines)", 
                       length(sample_ids)),
    x = "Time Post-Intubation (hours)",
    y = "Effective MAP (mmHg)"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray30"),
    legend.position = "bottom"
  )

ggsave("Figure3_Spaghetti.pdf", fig3, width = 10, height = 7)
ggsave("Figure3_Spaghetti.png", fig3, width = 10, height = 7, dpi = 300)
cat("   Saved: Figure3_Spaghetti.pdf/png\n")

#6 Figure 4: ICC by Time Stratum -----------
cat("\n6. Creating Figure 4: ICC by Time Between Intubations...\n")

# Load ICC results
time_icc <- fread("ICC_By_Time_Stratum.csv")

# Define the correct chronological order
level_order <- c("<24h", "24-72h", ">72h")

# Re-order the 'Stratum' column  
time_icc$Stratum <- factor(time_icc$Stratum, levels = level_order)

fig4 <- ggplot(time_icc, aes(x = Stratum, y = ICC)) +
  geom_col(fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f\n(n=%d)", ICC, N_Patients)), 
            vjust = -0.5, size = 4) +
  # geom_hline(yintercept = 0.5, linetype = "dashed", color = "red", linewidth = 0.8) +
  # annotate("text", x = 2.5, y = 0.52, label = "High consistency threshold", color = "red", hjust = 0.5) +
  labs(
    title = "ICC and Time Between Intubations",
    x = "Time Between Intubations",
    y = "Intraclass Correlation Coefficient (ICC)"
  ) +
  ylim(0, 0.4) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(color = "gray30")
  )

ggsave("Figure4_ICC_ByTime.pdf", fig4, width = 8, height = 6)
ggsave("Figure4_ICC_ByTime.png", fig4, width = 8, height = 6, dpi = 300)

