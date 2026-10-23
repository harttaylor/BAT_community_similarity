# Compare trend model results between Chorus and SM2 detector outputs
# Checks, for each species: trend direction, significance, and size of difference

library(dplyr)
library(readr)
library(knitr)
library(gt)

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
write_csv(comparison, "activitytrend_comparison_Chorus_vs_SM2.csv")

# ---- Optional: formatted table with highlighting ----
 display_table %>%
   gt() %>%
   tab_header(title = "Trend estimates: Chorus vs SM2") %>%
   tab_style(style = cell_fill(color = "#FDE9C8"),
             locations = cells_body(rows = Agreement != "Agree"))