library(tidyverse)
library(irr)        # for ICC — install.packages("irr") if needed
library(patchwork)  # for combining plots — install.packages("patchwork") if needed
library(lme4)

# =============================================================================
# DETECTOR AGREEMENT ANALYSIS: ICC, RMSE, and Bland-Altman
# =============================================================================
#
# This script takes a different angle than the GLMM ratio approach.
# Instead of asking "does Chorus detect more?", we ask:
#   1. How well do the two detectors AGREE? (ICC)
#   2. By how many calls do they typically DISAGREE? (RMSE)
#   3. Does the disagreement grow with activity level? (Bland-Altman)
#
# We do this at each temporal scale (10-min, 1-hour, nightly) to see
# whether agreement improves or degrades as you aggregate.
# =============================================================================

# counts_10, counts_60, and counts_night already exist in your R environment from the previous script.


# =============================================================================
# STEP 1: PIVOT TO PAIRED FORMAT
# =============================================================================
# To compute ICC and RMSE, we need paired observations: for each
# site × time-bin × species, one column for Chorus count and one for Other_FS.
# Bins where only one detector was active get dropped (no pair to compare).

make_paired <- function(count_data, time_col) {
  # time_col is the name of the time-bin column (string)
  # For nightly data, use "night"
  
  count_data %>%
    select(site_id, detector_type, all_of(time_col), spp, n_calls) %>%
    pivot_wider(
      names_from  = detector_type,
      values_from = n_calls,
      values_fill = list(n_calls = 0)
    ) %>%
    # Drop rows where either detector has no data (shouldn't happen if
    # you built the grid correctly, but just in case)
    filter(!is.na(Other_FS) & !is.na(Anabat_Chorus))
}

paired_10    <- make_paired(counts_10, "bin_10min")
paired_60    <- make_paired(counts_60, "hour_bin")
paired_night <- make_paired(counts_night, "night")

cat("Paired observations:\n")
cat("  10-min:", nrow(paired_10), "\n")
cat("  1-hour:", nrow(paired_60), "\n")
cat("  Night: ", nrow(paired_night), "\n\n")


# =============================================================================
# STEP 2: COMPUTE ICC AND RMSE PER SPECIES PER SCALE
# =============================================================================

compute_agreement <- function(paired_data, scale_label) {
  
  paired_data %>%
    group_by(spp) %>%
    summarise(
      scale = scale_label,
      n_pairs = n(),
      
      # --- ICC (two-way, agreement, single measures) ---
      # "agreement" flavour penalises systematic bias (Chorus always higher),
      # not just correlation. This is what you want.
      # If both columns are all zeros, ICC is undefined — handle gracefully.
      icc = tryCatch({
        icc_result <- icc(
          cbind(Other_FS, Anabat_Chorus),
          model  = "twoway",
          type   = "agreement",
          unit   = "single"
        )
        icc_result$value
      }, error = function(e) NA_real_),
      
      # --- RMSE (in raw call counts) ---
      rmse = sqrt(mean((Anabat_Chorus - Other_FS)^2)),
      
      # --- Mean absolute difference ---
      mae = mean(abs(Anabat_Chorus - Other_FS)),
      
      # --- Mean difference (bias direction) ---
      # Positive = Chorus detects more on average
      mean_diff = mean(Anabat_Chorus - Other_FS),
      
      # --- Percentage of bins where both agree on zero ---
      pct_both_zero = round(100 * mean(Other_FS == 0 & Anabat_Chorus == 0), 1),
      
      # --- Percentage of bins where they agree exactly ---
      pct_exact_match = round(100 * mean(Other_FS == Anabat_Chorus), 1),
      
      .groups = "drop"
    )
}

agreement_10    <- compute_agreement(paired_10, "10min")
agreement_60    <- compute_agreement(paired_60, "1hour")
agreement_night <- compute_agreement(paired_night, "night")

agreement_all <- bind_rows(agreement_10, agreement_60, agreement_night)
agreement_all$scale <- factor(agreement_all$scale, levels = c("10min", "1hour", "night"))

cat("=== AGREEMENT METRICS ===\n\n")
print(as.data.frame(agreement_all %>% arrange(spp, scale)), row.names = FALSE, digits = 3)

write.csv(agreement_all, "results_agreement_metrics.csv", row.names = FALSE)


# =============================================================================
# STEP 3: ICC AND RMSE COMPARISON PLOT
# =============================================================================
# Two-panel figure: ICC by scale (left) and RMSE by scale (right).
# This replaces the hard-to-read triple forest plot.

# Detector agreement (ICC)
p_icc <- agreement_all %>%
  filter(!is.na(icc)) %>%
  ggplot(aes(x = icc, y = reorder(spp, icc), colour = scale, shape = scale)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(
    "Temporal scale",
    values = c("10min" = "#E66100", "1hour" = "#5D3A9B", "night" = "#1A85FF")
  ) +
  scale_shape_manual(
    "Temporal scale",
    values = c("10min" = 16, "1hour" = 17, "night" = 15)
  ) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  geom_vline(xintercept = c(0.5, 0.75), linetype = "dotted", colour = "grey60") +
  labs(x = "ICC (agreement)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())
p_icc
ggsave("plot_icc.png", p_icc, width = 7, height = 7, dpi = 300)


# Typical disagreement (RMSE)
p_rmse <- agreement_all %>%
  ggplot(aes(x = rmse, y = reorder(spp, rmse), colour = scale, shape = scale)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  scale_colour_manual(
    "Temporal scale",
    values = c("10min" = "#E66100", "1hour" = "#5D3A9B", "night" = "#1A85FF")
  ) +
  scale_shape_manual(
    "Temporal scale",
    values = c("10min" = 16, "1hour" = 17, "night" = 15)
  ) +
  labs(x = "RMSE (calls)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom", panel.grid.minor = element_blank())

p_rmse

ggsave("plot_rmse.png", p_rmse, width = 7, height = 7, dpi = 300)



# =============================================================================
# STEP 4: BLAND-ALTMAN PLOTS
# =============================================================================
# A Bland-Altman plot shows:
#   x-axis: average of both detectors  (Chorus + Other_FS) / 2
#   y-axis: difference                  (Chorus - Other_FS)
#
# What to look for:
#   - Points centered on y = 0 → no systematic bias
#   - Points above 0 → Chorus detects more
#   - Fan shape (wider spread at higher x) → bias grows with activity
#   - Flat band → bias is constant regardless of activity level
#
# We make one Bland-Altman per species, faceted, at the NIGHTLY scale
# (clearest signal, fewest zeros). We add the mean difference (bias line)
# and ±1.96 SD "limits of agreement".

# --- BLAND-ALTMAN (nightly, all species) -------------------------------------
# Small tweak: increase point size slightly, ensure LOESS doesn't crash
# on species with very few points

ba_night <- paired_night %>%
  mutate(
    average    = (Anabat_Chorus + Other_FS) / 2,
    difference = Anabat_Chorus - Other_FS
  )

ba_lines <- ba_night %>%
  group_by(spp) %>%
  summarise(
    mean_diff = mean(difference),
    upper_loa = mean(difference) + 1.96 * sd(difference),
    lower_loa = mean(difference) - 1.96 * sd(difference),
    .groups   = "drop"
  )

p_ba <- ggplot(ba_night, aes(x = average, y = difference)) +
  geom_hline(yintercept = 0, colour = "grey40", linewidth = 0.4) +
  geom_point(alpha = 0.3, size = 1.5, colour = "#5D3A9B") +
  geom_hline(data = ba_lines, aes(yintercept = mean_diff),
             linetype = "dashed", colour = "#E66100", linewidth = 0.6) +
  geom_hline(data = ba_lines, aes(yintercept = upper_loa),
             linetype = "dotted", colour = "#E66100", linewidth = 0.5) +
  geom_hline(data = ba_lines, aes(yintercept = lower_loa),
             linetype = "dotted", colour = "#E66100", linewidth = 0.5) +
  # Only add LOESS if there are enough non-zero points
  geom_smooth(method = "loess", se = FALSE, colour = "#1A85FF",
              linewidth = 0.7, span = 0.8) +
  facet_wrap(~ spp, scales = "free", ncol = 4) +
  labs(
    x = "Average nightly count  (Chorus + Other_FS) / 2",
    y = "Difference  (Chorus \u2212 Other_FS)"
  ) +
  theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey95"),
    panel.grid.minor = element_blank()
  )
p_ba
ggsave("plot_bland_altman_nightly.png", p_ba, width = 14, height = 10, dpi = 300)



# =============================================================================
# STEP 5: SUMMARY TABLE 
# =============================================================================

summary_table <- agreement_all %>%
  filter(scale == "night") %>%
  select(spp, icc, rmse, mean_diff) %>%
  left_join(
    agreement_all %>%
      select(spp, scale, icc) %>%
      pivot_wider(names_from = scale, values_from = icc, names_prefix = "icc_"),
    by = "spp"
  ) %>%
  mutate(
    # Does ICC change much across scales?
    icc_range = abs(icc_10min - icc_night),
    # Interpretation
    bias_direction = case_when(
      mean_diff > 1  ~ "Chorus detects more",
      mean_diff < -1 ~ "Chorus detects less",
      TRUE           ~ "Similar"
    )
  ) %>%
  arrange(desc(abs(mean_diff)))



