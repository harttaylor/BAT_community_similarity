# =============================================================================
# Abundance analysis: Do the Chorus and the SM2 record different amounts of bat activity?
#
#   Data:   night_counts (from 01), SM2 sites only: nightly passes per species
#           per detector, zeros included.
#   Model:  one negative binomial GLMM per species (glmmTMB, nbinom2)
#             n_calls ~ detector_type + (1 | site_id/night)
#           (1 | site_id/night) = random effects for site and for site-night.
#           The site-night term pairs the two detectors' counts from the same
#           night (shared weather and bat activity), like a paired test.
#   Zero-inflation: a zero-inflated version is kept only if it lowers AIC by > 2.
#
# =============================================================================
#
# Runs the SAME steps and the SAME model on both classifier datasets:
#   kaleidoscope (main analysis) and sonobat (classifier comparison).
#
#   Paired check: each site's total passes per detector compared with a
#           paired Wilcoxon signed-rank test (sites as pairs); a simple,
#           assumption-light confirmation of the GLMM.
#   Total activity: the same model on all species summed per night.
#   Diagnostics: DHARMa simulated residuals, every species.
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(DHARMa)

setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity")
dir.create("results", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# MODEL FUNCTION (identical for every species and both classifiers)
# -----------------------------------------------------------------------------

fit_model <- function(d) {
  formulas <- list(
    "site/night + unit" = n_calls ~ detector_type + (1 | site_id/night) + (1 | site_id:detector_type),
    "site + unit"       = n_calls ~ detector_type + (1 | site_id) + (1 | site_id:detector_type),
    "site/night"        = n_calls ~ detector_type + (1 | site_id/night),
    "site"              = n_calls ~ detector_type + (1 | site_id)
  )
  for (re in names(formulas)) {
    m_nb <- glmmTMB(formulas[[re]], family = nbinom2, data = d)
    m_zi <- glmmTMB(formulas[[re]], family = nbinom2, ziformula = ~1, data = d)
    zi   <- isTRUE(AIC(m_zi) < AIC(m_nb) - 2)
    m    <- if (zi) m_zi else m_nb
    if (isTRUE(m$sdr$pdHess)) break
  }
  ct <- summary(m)$coefficients$cond["detector_typeChorus", ]
  b  <- ct[["Estimate"]]
  se <- ct[["Std. Error"]]
  
  # SD of the unit random effect (NA if the model used has no unit term)
  vc <- VarCorr(m)$cond[["site_id:detector_type"]]
  unit_sd <- if (is.null(vc)) NA else sqrt(vc[1, 1])
  
  list(model = m,
       result = tibble(ratio = exp(b), ci_lo = exp(b - 1.96 * se),
                       ci_hi = exp(b + 1.96 * se), p_value = ct[["Pr(>|z|)"]],
                       unit_sd = unit_sd,
                       site_ratio_lo = exp(b - 1.96 * sqrt(2) * unit_sd),
                       site_ratio_hi = exp(b + 1.96 * sqrt(2) * unit_sd),
                       random_effects = re, zero_inflated = zi))
}

common_names <- c(
  ANPA = "Pallid bat",             COTO = "Townsend's big-eared bat",
  EPFU = "Big brown bat",          EUMA = "Spotted bat",
  LABO = "Eastern red bat",        LACI = "Hoary bat",
  LANO = "Silver-haired bat",      MYCA = "California myotis",
  MYCI = "W. small-footed myotis", MYEV = "Long-eared myotis",
  MYLU = "Little brown myotis",    MYSE = "N. long-eared myotis",
  MYTH = "Fringed myotis",         MYVO = "Long-legged myotis",
  MYYU = "Yuma myotis",            PAHE = "Canyon bat"
)

all_results <- list()


# =============================================================================
# SAME ANALYSIS FOR EACH CLASSIFIER
# =============================================================================

for (classifier in c("kaleidoscope", "sonobat")) {
  
  print(classifier)
  bats <- readRDS(paste0("data/processed/bat_data_", classifier, ".rds"))
  sm2  <- filter(bats$night_counts, legacy_model == "SM2")
  
  # Sites, site-nights and passes used
  print(c(sites = n_distinct(sm2$site_id),
          site_nights = n_distinct(sm2$site_id, sm2$night),
          passes = sum(sm2$n_calls)))
  
  # Overdispersion check (ratio well above 1 justifies negative binomial)
  m_pois <- glm(n_calls ~ detector_type + spp + site_id, data = sm2, family = poisson)
  print(deviance(m_pois) / df.residual(m_pois))
  
  # ---- Total activity (all species combined) --------------------------------
  total <- sm2 %>%
    group_by(site_id, night, detector_type) %>%
    summarise(n_calls = sum(n_calls), .groups = "drop")
  print(fit_model(total)$result)
  
  # ---- One model per species ------------------------------------------------
  results <- list()
  models  <- list()
  for (sp in bats$spp_keep) {
    d   <- filter(sm2, spp == sp)
    fit <- fit_model(d)
    results[[sp]] <- bind_cols(tibble(species = sp, passes = sum(d$n_calls)), fit$result)
    models[[sp]]  <- fit$model
  }
  results <- bind_rows(results)
  results$p_adj       <- p.adjust(results$p_value, method = "BH")
  results$common_name <- common_names[results$species]
  results$classifier  <- classifier
  results <- arrange(results, ratio)
  
  # ---- Simple paired check: one Chorus-vs-SM2 pair per site -----------------
  # Each site's total passes per detector, then a paired Wilcoxon signed-rank
  # test across sites (the classic paired test, with sites as the pairs).
  # Its p-values should broadly agree with the GLMM's.
  paired_check <- sm2 %>%
    group_by(spp, site_id, detector_type) %>%
    summarise(passes = sum(n_calls), .groups = "drop") %>%
    pivot_wider(names_from = detector_type, values_from = passes) %>%
    group_by(spp) %>%
    summarise(sites = n(),
              median_site_ratio = median((Chorus + 0.5) / (Legacy + 0.5)),
              p_paired_wilcoxon = suppressWarnings(
                wilcox.test(Chorus, Legacy, paired = TRUE)$p.value))
  results <- left_join(results, paired_check, by = c("species" = "spp"))
  
  print(select(results, species, passes, ratio, ci_lo, ci_hi, p_value,
               site_ratio_lo, site_ratio_hi, median_site_ratio, p_paired_wilcoxon,
               random_effects, zero_inflated), n = Inf)
  write.csv(results, paste0("results/NB_models_", classifier, ".csv"), row.names = FALSE)
  all_results[[classifier]] <- results
  
  # ---- Per-site ratios (how consistent is the ratio between sites?) ---------
  # Raw ratio at each site: (Chorus passes + 0.5) / (SM2 passes + 0.5).
  # The 0.5 avoids dividing by zero. Sites with < 5 passes on both detectors
  # are left out (ratio too noisy to read).
  site_ratios <- sm2 %>%
    group_by(spp, site_id, detector_type) %>%
    summarise(passes = sum(n_calls), .groups = "drop") %>%
    pivot_wider(names_from = detector_type, values_from = passes) %>%
    filter(Chorus + Legacy >= 5) %>%
    mutate(ratio = (Chorus + 0.5) / (Legacy + 0.5))
  write.csv(site_ratios, paste0("results/site_ratios_", classifier, ".csv"), row.names = FALSE)
  
  # Sites whose total-activity ratio is furthest from 1 (check for hardware problems)
  print(sm2 %>%
          group_by(site_id, detector_type) %>%
          summarise(passes = sum(n_calls), .groups = "drop") %>%
          pivot_wider(names_from = detector_type, values_from = passes) %>%
          mutate(ratio = (Chorus + 0.5) / (Legacy + 0.5)) %>%
          arrange(desc(abs(log(ratio)))), n = Inf)
  
  p_sites <- ggplot(site_ratios, aes(x = ratio, y = spp)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_jitter(height = 0.15, width = 0, alpha = 0.6, size = 1.8, colour = "#5D3A9B") +
    scale_x_log10() +
    labs(x = "Chorus / SM2 ratio at each site", y = NULL) +
    theme_classic(base_size = 12)
  ggsave(paste0("results/plot_site_ratios_", classifier, ".png"), p_sites,
         width = 6.5, height = 6, dpi = 600)
  
  # ---- Diagnostics (DHARMa) -------------------------------------------------
  # p < 0.05 = possible problem; check that species' plot.
  dharma_dir <- paste0("results/dharma_", classifier)
  dir.create(dharma_dir, showWarnings = FALSE)
  diagnostics <- list()
  for (sp in names(models)) {
    sim <- simulateResiduals(models[[sp]], n = 1000)
    png(file.path(dharma_dir, paste0(sp, ".png")), width = 1000, height = 500)
    plot(sim, title = sp)
    dev.off()
    diagnostics[[sp]] <- tibble(
      species        = sp,
      p_uniformity   = testUniformity(sim, plot = FALSE)$p.value,
      p_dispersion   = testDispersion(sim, plot = FALSE)$p.value,
      p_zero_inflate = testZeroInflation(sim, plot = FALSE)$p.value
    )
  }
  diagnostics <- bind_rows(diagnostics)
  print(diagnostics, n = Inf)
  write.csv(diagnostics, paste0("results/NB_diagnostics_", classifier, ".csv"), row.names = FALSE)
  
  # ---- Forest plot (coloured by unadjusted p < 0.05) ------------------------
  pd <- results
  pd$effect <- case_when(
    pd$p_value < 0.05 & pd$ratio > 1 ~ "Chorus higher",
    pd$p_value < 0.05 & pd$ratio < 1 ~ "Chorus lower",
    TRUE                             ~ "No difference"
  )
  pd$species <- factor(pd$species, levels = pd$species)
  
  p_forest <- ggplot(pd, aes(x = ratio, y = species, colour = effect)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.6) +
    geom_point(size = 2.6) +
    scale_x_log10(breaks = c(0.125, 0.25, 0.5, 1, 2, 4),
                  labels = c("0.125", "0.25", "0.5", "1", "2", "4")) +
    scale_colour_manual(values = c("Chorus higher" = "#5D3A9B",
                                   "Chorus lower"  = "#E66100",
                                   "No difference" = "grey60"), name = NULL) +
    labs(x = "Detection ratio (Chorus / SM2)", y = NULL) +
    theme_classic(base_size = 12) +
    theme(legend.position = "bottom", axis.text = element_text(colour = "grey20"))
  print(p_forest)
  ggsave(paste0("results/plot_NB_", classifier, ".png"), p_forest, width = 6.5, height = 6, dpi = 600)
  ggsave(paste0("results/plot_NB_", classifier, ".pdf"), p_forest, width = 6.5, height = 6)
}


# =============================================================================
# KALEIDOSCOPE vs SONOBAT, side by side
# =============================================================================

compare <- bind_rows(all_results)
compare$species <- factor(compare$species, levels = sort(unique(compare$species), decreasing = TRUE))

print(compare %>%
        select(species, classifier, passes, ratio, ci_lo, ci_hi, p_value) %>%
        arrange(species, classifier), n = Inf)

p_compare <- ggplot(compare, aes(x = ratio, y = species, colour = classifier)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.3, linewidth = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(size = 2.4, position = position_dodge(width = 0.6)) +
  scale_x_log10() +
  scale_colour_manual(values = c("kaleidoscope" = "#5D3A9B", "sonobat" = "#E66100"),
                      labels = c("Kaleidoscope", "SonoBat"), name = NULL) +
  labs(x = "Detection ratio (Chorus / SM2)", y = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "bottom")
print(p_compare)
ggsave("results/plot_classifier_comparison.png", p_compare, width = 6.5, height = 7, dpi = 600)

