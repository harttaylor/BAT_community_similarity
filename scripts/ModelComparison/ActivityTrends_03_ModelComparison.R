# Compare trend model results between Chorus and SM2 detector outputs
# Checks, for each species: trend direction, significance, and size of difference

library(dplyr)
library(readr)
library(knitr)
library(gt)
library(tidyr)
library(ggplot2)

# ---- Settings ----
chorus_file <- "Data/Analyzed/ActivityModels/bc.cov.estimates.singlet.Chorus.csv"
sm2_file    <- "Data/Analyzed/ActivityModels/bc.cov.estimates.singlet.SM2.csv"
alpha       <- 0.05

# ---- Read and tidy ----
read_est <- function(file, suffix) {
  read_csv(file, show_col_types = FALSE) %>%
    select(SpeciesGroup, est = Estimate, se = SE, p = p.value, n) %>%
    rename_with(~ paste0(.x, "_", suffix), -SpeciesGroup)
}

chorus <- read_est(chorus_file, "C")
sm2    <- read_est(sm2_file,    "S")

# Label a single trend as increase / decrease / n.s.
classify_trend <- function(est, p) {
  case_when(
    p < alpha & est > 0 ~ "Sig. increase",
    p < alpha & est < 0 ~ "Sig. decrease",
    TRUE                ~ "n.s."
  )
}

# ---- Compare ----
comparison <- full_join(chorus, sm2, by = "SpeciesGroup") %>%
  mutate(
    trend_C        = classify_trend(est_C, p_C),
    trend_S        = classify_trend(est_S, p_S),
    same_direction = sign(est_C) == sign(est_S),
    same_sig       = (p_C < alpha) == (p_S < alpha),
    est_diff       = est_C - est_S,
    # Difference expressed in units of the average SE of the two models
    # (the two models share data, so this is a rough guide, not a formal test)
    diff_in_SE     = abs(est_diff) / ((se_C + se_S) / 2),
    n_match        = n_C == n_S,
    agreement = case_when(
      same_direction & same_sig                 ~ "Agree",
      !same_direction & p_C >= alpha & p_S >= alpha ~ "Sign flip, both n.s.",
      same_direction & !same_sig                ~ "Significance differs",
      TRUE                                      ~ "Disagree"
    )
  ) %>%
  arrange(SpeciesGroup)

# ---- Display table ----
display_table <- comparison %>%
  transmute(
    Species        = SpeciesGroup,
    `Chorus est (p)` = sprintf("%+.3f (%.3f)", est_C, p_C),
    `SM2 est (p)`    = sprintf("%+.3f (%.3f)", est_S, p_S),
    `Chorus trend`   = trend_C,
    `SM2 trend`      = trend_S,
    `Diff (in SE)`   = sprintf("%.2f", diff_in_SE),
    `n (C/S)`        = ifelse(n_match, as.character(n_C), paste0(n_C, "/", n_S)),
    Agreement        = agreement
  )

print(kable(display_table, align = "lllllrrl"))

# ---- Summary ----
cat("\nSame direction:   ", sum(comparison$same_direction), "of", nrow(comparison), "\n")
cat("Same significance:", sum(comparison$same_sig),       "of", nrow(comparison), "\n")
cat("Sample size mismatch:",
    paste(comparison$SpeciesGroup[!comparison$n_match], collapse = ", "), "\n")
print(count(comparison, agreement))

# ---- Save ----
write_csv(comparison, "Data/Analyzed/ActivityModels/activitytrend_comparison_Chorus_vs_SM2.csv")

# Compare species-specific activity trends (GLMMs): Chorus vs SM2 data
# Input:  activitytrend_comparison_Chorus_vs_SM2.csv
# Output: fig_activity_trend_comparison.png
#         Trend estimate +/- 95% CI (est +/- 1.96 SE) for each data source

act <- comparison

# order species by mean estimate across both data sources
sp_order <- act |>
  mutate(m = (est_C + est_S) / 2) |>
  arrange(m) |>
  pull(SpeciesGroup)

# long format: one row per species x data source
act_long <- act |>
  select(SpeciesGroup, est_C, se_C, p_C, est_S, se_S, p_S) |>
  pivot_longer(-SpeciesGroup,
               names_to = c(".value", "detector"),
               names_sep = "_") |>
  mutate(detector = recode(detector, C = "Chorus", S = "SM2"),
         lo95 = est - 1.96 * se,
         hi95 = est + 1.96 * se,
         sig  = p < 0.05,
         SpeciesGroup = factor(SpeciesGroup, levels = sp_order))

pda <- position_dodge(width = 0.6)

pa <- ggplot(act_long, aes(x = est, y = SpeciesGroup,
                          colour = detector, group = detector)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95),
                 height = 0, linewidth = 0.8, position = pd) +
  geom_point(aes(shape = sig), size = 2.2, position = pda,
             fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE` = "p < 0.05", `FALSE` = "n.s."),
                     breaks = c("TRUE", "FALSE"), name = NULL) +
  scale_colour_manual(values = det_cols, name = "Data source") +
  labs(x = "Activity trend estimate (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"),
        legend.position = "bottom")
pa

ggsave("Figures/fig_activity_trend_comparison.png", pa, width = 6.5, height = 6.5, dpi = 600)
