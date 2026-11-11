# =============================================================================
# QUESTION 2: How much recording time per night is needed to see the
# Chorus/SM2 difference?
#
# Kaleidoscope IDs only.
#
# OUTPUTS (results/ folder)
#   fig2_sampling_time.png       Fig. 2: Chorus/SM2 ratio (95% CI) per species
#                                for each amount of recording time per night
#   table2_sampling_time.csv     ratio, 95% CI and p for every species and window
#   table2b_shortest_window.csv  shortest window in which each species differed
#
# DATA
#   Each night starts at the first identified bat pass at that site that night
#   (recorded by either unit, so the start time is the same for both units).
#   Windows: the first 10 min, 30 min, 1 h and 2 h of each night, and the whole
#   night. For each window, each unit's passes inside the window are summed over
#   all paired nights at the site: two totals per site per species, exactly as
#   in Question 1. "Whole night" is the Question 1 data.
#
# MODEL: the Question 1 model, unchanged
#   n_calls ~ detector_type + (1 | site_id), negative binomial (glmmTMB)
#   The "Whole night" results must match Table 1 exactly.
#   A species is fitted in a window only if it has >= 20 passes in that window
#   (both units combined).
#
# WHAT TO LOOK FOR
#   Shorter windows = less data, so confidence intervals get wider.
#   - Ratio about the same in every window, CI only excludes 1 in longer
#     windows: the difference is real but needs that much recording to detect.
#   - Ratio changes with window length: the difference depends on the time of
#     night (e.g. different species or activity levels early in the night).
# =============================================================================

library(tidyverse)
library(glmmTMB)

setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity")
dir.create("results", showWarnings = FALSE)

bats <- readRDS("data/processed/bat_data_kaleidoscope.rds")


# -----------------------------------------------------------------------------
# MINUTES SINCE THE START OF EACH NIGHT
# -----------------------------------------------------------------------------

# Start of each night = first identified pass at the site (either unit, any species)
night_start <- bats$dat %>%
  group_by(site_id, night) %>%
  summarise(start = min(Timestamp), .groups = "drop")

passes <- bats$dat %>%
  filter(spp %in% bats$spp_keep) %>%
  left_join(night_start, by = c("site_id", "night")) %>%
  mutate(minutes = as.numeric(difftime(Timestamp, start, units = "mins")))

windows <- c("10 min" = 10, "30 min" = 30, "1 h" = 60, "2 h" = 120, "Whole night" = Inf)

# Every site x unit x species, so a unit with no passes in a window gets a 0
all_rows <- crossing(site_id       = unique(bats$night_counts$site_id),
                     detector_type = c("Legacy", "Chorus"),
                     spp           = bats$spp_keep)


# -----------------------------------------------------------------------------
# MODEL AND RESULT EXTRACTION (same as 02_abundance.R)
# -----------------------------------------------------------------------------

model_formula <- n_calls ~ detector_type + (1 | site_id)

fit_nb <- function(d) {
  m <- glmmTMB(model_formula, family = nbinom2, data = d)
  if (!isTRUE(m$sdr$pdHess)) {   # not converged: same model, different optimizer
    m <- glmmTMB(model_formula, family = nbinom2, data = d,
                 control = glmmTMBControl(optimizer = optim, optArgs = list(method = "BFGS")))
  }
  m
}

chorus_effect <- function(m) {
  ct <- summary(m)$coefficients$cond["detector_typeChorus", ]
  tibble(ratio     = exp(ct[["Estimate"]]),
         ci_lo     = exp(ct[["Estimate"]] - 1.96 * ct[["Std. Error"]]),
         ci_hi     = exp(ct[["Estimate"]] + 1.96 * ct[["Std. Error"]]),
         p_value   = ct[["Pr(>|z|)"]],
         converged = isTRUE(m$sdr$pdHess))
}


# -----------------------------------------------------------------------------
# FIT EVERY SPECIES IN EVERY WINDOW
# -----------------------------------------------------------------------------

results <- list()

for (w in names(windows)) {
  
  # Each unit's passes inside the window, summed per site (zeros filled in)
  totals <- passes %>%
    filter(minutes < windows[[w]]) %>%
    count(site_id, detector_type, spp, name = "n_calls")
  totals <- all_rows %>%
    left_join(totals, by = c("site_id", "detector_type", "spp")) %>%
    mutate(n_calls = replace_na(n_calls, 0),
           detector_type = factor(detector_type, levels = c("Legacy", "Chorus")))
  
  for (sp in bats$spp_keep) {
    d <- filter(totals, spp == sp)
    if (sum(d$n_calls) < 20) next
    results[[paste(w, sp)]] <- bind_cols(
      tibble(window = w, species = sp,
             passes_chorus = sum(d$n_calls[d$detector_type == "Chorus"]),
             passes_sm2    = sum(d$n_calls[d$detector_type == "Legacy"])),
      chorus_effect(fit_nb(d)))
  }
}

results <- bind_rows(results)
results$window <- factor(results$window, levels = names(windows))
results <- arrange(results, species, window)
print(results, n = Inf, width = Inf)
write.csv(results, "results/table2_sampling_time.csv", row.names = FALSE)

# Shortest window in which each species showed a difference (p < 0.05)
shortest <- results %>%
  filter(p_value < 0.05) %>%
  arrange(species, window) %>%
  group_by(species) %>%
  slice(1) %>%
  ungroup() %>%
  select(species, shortest_window = window, ratio, ci_lo, ci_hi)
print(shortest, n = Inf)
write.csv(shortest, "results/table2b_shortest_window.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# FIGURE 2: ratio for each recording window, one panel per species
# -----------------------------------------------------------------------------
# Same colours as Fig. 1: purple/orange = Chorus higher/lower (p < 0.05).
# All panels share the same y-axis so species can be compared. Intervals
# wider than the axis (rare species) are cut off at the panel edge.
# Panels ordered as in Fig. 1 (whole-night ratio, highest first).

results$effect <- case_when(
  results$p_value < 0.05 & results$ratio > 1 ~ "Chorus higher",
  results$p_value < 0.05 & results$ratio < 1 ~ "Chorus lower",
  TRUE                                       ~ "No difference"
)
panel_order <- results %>% filter(window == "Whole night") %>% arrange(desc(ratio)) %>% pull(species)
results$species <- factor(results$species, levels = panel_order)

fig2 <- ggplot(results, aes(x = window, y = ratio)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_line(aes(group = species), colour = "grey80", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi, colour = effect), width = 0.2, linewidth = 0.5) +
  geom_point(aes(colour = effect), size = 1.8) +
  scale_y_log10(breaks = c(0.25, 0.5, 1, 2, 4), labels = c("0.25", "0.5", "1", "2", "4")) +
  scale_x_discrete(labels = c("10", "30", "60", "120", "All")) +
  coord_cartesian(ylim = c(0.25, 5)) +
  scale_colour_manual(values = c("Chorus higher" = "#5D3A9B",
                                 "Chorus lower"  = "#E66100",
                                 "No difference" = "grey60"), name = NULL) +
  facet_wrap(~ species, ncol = 5) +
  labs(x = "Recording time per night (minutes)", y = "Detection ratio (Chorus / SM2)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", strip.background = element_blank())

ggsave("results/fig4_sampling_time.png", fig2, width = 9, height = 6.5, dpi = 600)
ggsave("results/fig4_sampling_time.pdf", fig2, width = 9, height = 6.5)
