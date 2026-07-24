# =============================================================================
# TEMPORAL SCALE ANALYSIS: Does detector bias depend on aggregation grain?
# =============================================================================

library(tidyverse)
library(lme4)
library(glmmTMB)
library(DHARMa)

# --- Load cleaned data (from data-prep script) ---------------------
# This is the pass-level data BEFORE aggregation to night counts.
# It has one row per detected bat pass with a Timestamp.

#dat <- read.csv("bat_detections.csv", stringsAsFactors = FALSE)
dat$Timestamp <- as.POSIXct(dat$Timestamp, tz = "America/Vancouver")
dat$night     <- as.Date(dat$night)
dat$detector_type <- factor(dat$detector_type, levels = c("Other_FS", "Anabat_Chorus"))

# If your bin columns don't exist yet, create them:
dat$hour_bin  <- floor_date(dat$Timestamp, unit = "1 hour")
dat$bin_10min <- floor_date(dat$Timestamp, unit = "10 minutes")


# --- Species to model (same filter as your main analysis) --------------------
spp_keep <- dat %>%
  group_by(spp) %>%
  summarise(total = n(), site_nights = n_distinct(paste(site_id, night))) %>%
  filter(total >= 50, site_nights >= 10) %>%
  pull(spp)

paste(spp_keep, collapse = ", ")


# =============================================================================
# STEP 1: BUILD COUNT DATASETS AT EACH TEMPORAL SCALE
# =============================================================================
# The critical step: for each scale, we need:
#   (a) All active detector × time-bin combinations (the "grid")
#   (b) Observed counts per species per bin
#   (c) Fill missing combos with zeros
#
# The grid defines "when was each detector listening?" — any bin where a
# detector recorded ANYTHING (any species, any call) was active. Bins with
# no detections of a focal species get a zero, not NA.

# ---- 10-minute scale --------------------------------------------------------

grid_10<- dat %>% distinct(site_id, detector_type, night, bin_10min)

obs_10 <- dat %>%
  filter(spp %in% spp_keep) %>%
  group_by(site_id, detector_type, night, bin_10min, spp) %>%
  summarise(n_calls = n(), .groups = "drop")

counts_10 <- crossing(grid_10, spp = spp_keep) %>%
  left_join(obs_10, by = c("site_id", "detector_type", "night", "bin_10min", "spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

nrow(counts_10)
round(100 * mean(counts_10$n_calls == 0), 1)


# ---- 1-hour scale -----------------------------------------------------------

grid_60 <- dat %>% distinct(site_id, detector_type, night, hour_bin)

obs_60 <- dat %>%
  filter(spp %in% spp_keep) %>%
  group_by(site_id, detector_type, night, hour_bin, spp) %>%
  summarise(n_calls = n(), .groups = "drop")

counts_60 <- crossing(grid_60, spp = spp_keep) %>%
  left_join(obs_60, by = c("site_id", "detector_type", "night", "hour_bin", "spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

nrow(counts_60)
round(100 * mean(counts_60$n_calls == 0), 1)


# ---- Nightly scale ( already have this — rebuild for consistency) --------

grid_night <- dat %>% distinct(site_id, detector_type, night)

obs_night <- dat %>%
  filter(spp %in% spp_keep) %>%
  group_by(site_id, detector_type, night, spp) %>%
  summarise(n_calls = n(), .groups = "drop")

counts_night <- crossing(grid_night, spp = spp_keep) %>%
  left_join(obs_night, by = c("site_id", "detector_type", "night", "spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

nrow(counts_night)
print(round(100 * mean(counts_night$n_calls == 0), 1))


# =============================================================================
# STEP 2: DIAGNOSTIC — sample sizes and zero proportions by scale & species
# =============================================================================
# this helps decide if finer scales are  viable for rare species.
# If a species has >95% zeros at 30-min, the model will struggle or need
# zero-inflation. may want to drop finest scale for rare species.

diag_table <- bind_rows(
  counts_10 %>% group_by(spp) %>%
    summarise(scale = "10min", n_rows = n(), n_nonzero = sum(n_calls > 0),
              pct_zero = round(100 * mean(n_calls == 0), 1), .groups = "drop"),
  counts_60    %>% group_by(spp) %>%
    summarise(scale = "1hour", n_rows = n(), n_nonzero = sum(n_calls > 0),
              pct_zero = round(100 * mean(n_calls == 0), 1), .groups = "drop"),
  counts_night %>% group_by(spp) %>%
    summarise(scale = "night", n_rows = n(), n_nonzero = sum(n_calls > 0),
              pct_zero = round(100 * mean(n_calls == 0), 1), .groups = "drop")
) %>%
  arrange(spp, scale)


write.csv(diag_table, "temporal_scale_diagnostics.csv", row.names = FALSE)


# =============================================================================
# STEP 3: FIT MODELS AT EACH SCALE
# =============================================================================
# 
# Model structure:
#   30-min & 1-hour:  n_calls ~ detector_type + (1|site_id/night)
#     - The nested random effect accounts for:
#       (a) site-level variation (some sites are busier)
#       (b) night-level variation within sites (some nights are busier)
#       (c) the fact that bins within a night are NOT independent
#
#   Nightly:          n_calls ~ detector_type + (1|site_id)
#     - Only site-level random effect (each row is already a full night)
#
# We try glmer.nb first. If zero-inflation is extreme (DHARMa flags it),
# we fall back to glmmTMB with ziformula = ~1.
#
# We wrap everything in a function so we don't repeat code 48 times.

fit_one_species <- function(species_code, count_data, scale_label) {
  
  d <- filter(count_data, spp == species_code)
  
  # Skip if too few non-zero observations
  if (sum(d$n_calls > 0) < 10) {
    cat("  ", species_code, "@", scale_label, ": skipped (< 10 non-zero obs)\n")
    return(data.frame(
      spp = species_code, scale = scale_label,
      ratio = NA, ci_lo = NA, ci_hi = NA, p_value = NA,
      model_type = "skipped", converged = NA,
      n_rows = nrow(d), n_nonzero = sum(d$n_calls > 0),
      pct_zero = round(100 * mean(d$n_calls == 0), 1)
    ))
  }
  
  # Choose random effect structure based on scale
  if (scale_label == "night") {
    formula_nb <- n_calls ~ detector_type + (1 | site_id)
    formula_zi <- n_calls ~ detector_type + (1 | site_id)
  } else {
    formula_nb <- n_calls ~ detector_type + (1 | site_id / night)
    formula_zi <- n_calls ~ detector_type + (1 | site_id / night)
  }
  
  # --- Attempt 1: glmer.nb (standard NB GLMM) ---
  model_type <- "glmer.nb"
  converged  <- TRUE
  
  m <- tryCatch(
    glmer.nb(formula_nb, data = d,
             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000))),
    error   = function(e) NULL,
    warning = function(w) {
      # Catch convergence warnings but still return the model
      m_inner <- suppressWarnings(
        glmer.nb(formula_nb, data = d,
                 control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 50000)))
      )
      attr(m_inner, "had_warning") <- TRUE
      m_inner
    }
  )
  
  if (!is.null(m) && !is.null(attr(m, "had_warning"))) converged <- FALSE
  
  # --- Attempt 2: if glmer.nb failed, try glmmTMB with zero-inflation ---
  if (is.null(m)) {
    model_type <- "glmmTMB_zi"
    m <- tryCatch(
      glmmTMB(formula_zi, ziformula = ~1, family = nbinom2, data = d),
      error = function(e) NULL
    )
    if (is.null(m)) {
      cat("  ", species_code, "@", scale_label, ": both models failed\n")
      return(data.frame(
        spp = species_code, scale = scale_label,
        ratio = NA, ci_lo = NA, ci_hi = NA, p_value = NA,
        model_type = "failed", converged = FALSE,
        n_rows = nrow(d), n_nonzero = sum(d$n_calls > 0),
        pct_zero = round(100 * mean(d$n_calls == 0), 1)
      ))
    }
  }
  
  # --- Extract results ---
  if (model_type == "glmer.nb") {
    coef_val <- fixef(m)["detector_typeAnabat_Chorus"]
    se_val   <- sqrt(diag(vcov(m)))["detector_typeAnabat_Chorus"]
    p_val    <- summary(m)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]
  } else {
    # glmmTMB stores conditional (non-ZI) fixed effects in $cond
    coef_val <- fixef(m)$cond["detector_typeAnabat_Chorus"]
    se_val   <- sqrt(diag(vcov(m)$cond))["detector_typeAnabat_Chorus"]
    p_val    <- summary(m)$coefficients$cond["detector_typeAnabat_Chorus", "Pr(>|z|)"]
  }
  
  ratio <- exp(coef_val)
  ci_lo <- exp(coef_val - 1.96 * se_val)
  ci_hi <- exp(coef_val + 1.96 * se_val)
  
  cat("  ", species_code, "@", scale_label, ":", model_type,
      " ratio =", round(ratio, 3),
      " CI [", round(ci_lo, 3), "-", round(ci_hi, 3), "]",
      " p =", round(p_val, 4), "\n")
  
  data.frame(
    spp = species_code, scale = scale_label,
    ratio = ratio, ci_lo = ci_lo, ci_hi = ci_hi, p_value = p_val,
    model_type = model_type, converged = converged,
    n_rows = nrow(d), n_nonzero = sum(d$n_calls > 0),
    pct_zero = round(100 * mean(d$n_calls == 0), 1)
  )
}


# --- Run all species × all scales -------------------------------------------
cat("\n=== FITTING MODELS ===\n\n")

all_results <- bind_rows(
  # 10-min scale
  map_dfr(spp_keep, ~fit_one_species(.x, counts_10, "10min")),
  # 1-hour scale
  map_dfr(spp_keep, ~fit_one_species(.x, counts_60, "1hour")),
  # Nightly scale
  map_dfr(spp_keep, ~fit_one_species(.x, counts_night, "night"))
)

all_results$scale <- factor(all_results$scale, levels = c("10min", "1hour", "night"))


print(as.data.frame(all_results), row.names = FALSE)

write.csv(all_results, "results_temporal_scale_comparison.csv", row.names = FALSE)

library(ggstance)
install.packages("ggstance")
# =============================================================================
# STEP 4: COMPARISON PLOT — ratio by scale for each species
# =============================================================================
# This is the key figure. For each species, we show the detection ratio at
# all three scales side by side. If the ratios are consistent, the dots
# line up vertically. If they diverge, temporal scale matters.

plot_dat <- all_results %>% filter(!is.na(ratio))
plot_dat$spp <- factor(plot_dat$spp, levels = rev(sort(unique(plot_dat$spp))))

ggplot(plot_dat, aes(x = ratio, y = spp, colour = scale, shape = scale)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbar(
    aes(xmin = ci_lo, xmax = ci_hi),
    width = 0.3, linewidth = 0.5,
    position = position_dodge(width = 0.6),
    orientation = "y"
  ) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25", "0.5", "1", "2", "4")
  ) +
  scale_colour_manual(values = c("10min" = "#E66100", "1hour" = "#5D3A9B", "night" = "#1A85FF")) +
  scale_shape_manual(values = c("10min" = 16, "1hour" = 17, "night" = 15)) +
  labs(
    x = "Detection ratio: Chorus / Other_FS  (log scale)",
    y = NULL,
    colour = "Temporal scale",
    shape  = "Temporal scale",
    title  = "Detector bias across temporal aggregation scales"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "bottom",
    panel.grid.major.y = element_line(colour = "grey90"),
    panel.grid.minor   = element_blank()
  )

ggsave("plot_temporal_scale_comparison.png", width = 10, height = 8, dpi = 150)





    