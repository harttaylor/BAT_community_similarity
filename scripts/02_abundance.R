# =============================================================================
# QUESTION 1: For each species, does the Chorus record a different
# number of bat passes than the SM2 recording beside it?
#
# Kaleidoscope IDs only. The SonoBat check is in 05_supp_classifier.R and
# does not change anything here.
#
# OUTPUTS (results/ folder)
#   fig1_species_ratios.png      Fig. 1: Chorus/SM2 ratio per species
#   table1_species_ratios.csv    Table 1: passes per detector, ratio, 95% CI, p
#   total_activity.csv           all species combined
#   tableS1_model_checks.csv     residual checks for every species (supplement)
#   model_checks/<species>.png   residual plot for every species
#
# DATA
#   For each site, each unit's passes are summed over that site's paired nights
#   (nights when both units were recording). This gives two numbers per site
#   per species: the Chorus total and the SM2 total. Both totals cover exactly
#   the same nights, so effort is equal within each site (no offset needed).
#
# THE MODEL (the same for every species, no model selection)
#   n_calls ~ detector_type + (1 | site_id)
#   negative binomial GLMM (glmmTMB, family = nbinom2)
#
#   n_calls         total passes of one species on one unit at one site
#   detector_type   Chorus vs SM2  <- the effect we are testing
#   (1 | site_id)   pairs the two units at each site (some sites have more bats)
#   This is the count-data version of a paired t-test with sites as the pairs:
#   each site gives one Chorus vs SM2 comparison, so the unit pairs (one per
#   site) are the replicates.
#
#   Result: exp(detector coefficient) = Chorus / SM2 ratio of passes.
#           1.5 = Chorus recorded 50% more; 0.5 = Chorus recorded half as many.
#   P-values: unadjusted (BH-adjusted values kept in the CSV for reviewers).
#   Negative binomial: counts are far more variable than a Poisson model allows
#           (overdispersion ratio printed below).
#   Convergence: if a model does not converge with the default optimizer it is
#           refitted with another optimizer (BFGS). Same model, same terms.
#           The 'converged' column records the outcome.
#
# MODEL CHECKS (DHARMa simulated residuals, every species)
#   uniformity      do the residuals follow the distribution the model expects?
#   dispersion      is the spread of the data what the model expects?
#   zero_inflation  are there more zeros than the model expects?
#   p < 0.05 flags a possible problem; check the residual plot.
# =============================================================================

library(tidyverse)
library(glmmTMB)
library(DHARMa)

setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity")
dir.create("results/model_checks", recursive = TRUE, showWarnings = FALSE)

bats <- readRDS("data/processed/bat_data_kaleidoscope.rds")

# Data used (numbers for the methods)
print(c(sites       = n_distinct(bats$night_counts$site_id),
        site_nights = n_distinct(bats$night_counts$site_id, bats$night_counts$night),
        species     = length(bats$spp_keep),
        passes      = sum(bats$night_counts$n_calls)))

# Sum each unit's passes over the paired nights: two rows per site per species
site_totals <- bats$night_counts %>%
  group_by(spp, site_id, detector_type) %>%
  summarise(n_calls = sum(n_calls), .groups = "drop")

# Overdispersion: ratio well above 1 = use negative binomial, not Poisson
m_pois <- glm(n_calls ~ detector_type + spp + site_id, data = site_totals, family = poisson)
print(deviance(m_pois) / df.residual(m_pois))


# -----------------------------------------------------------------------------
# MODEL AND RESULT EXTRACTION
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
# TOTAL ACTIVITY (all species combined)
# -----------------------------------------------------------------------------

total <- site_totals %>%
  group_by(site_id, detector_type) %>%
  summarise(n_calls = sum(n_calls), .groups = "drop")

total_result <- chorus_effect(fit_nb(total))
print(total_result)
write.csv(total_result, "results/total_activity.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# ONE MODEL PER SPECIES
# -----------------------------------------------------------------------------

results <- list()
models  <- list()

for (sp in bats$spp_keep) {
  d <- filter(site_totals, spp == sp)
  m <- fit_nb(d)
  results[[sp]] <- bind_cols(
    tibble(species       = sp,
           passes_chorus = sum(d$n_calls[d$detector_type == "Chorus"]),
           passes_sm2    = sum(d$n_calls[d$detector_type == "Legacy"]),
           sites_with_species = sum(tapply(d$n_calls, d$site_id, sum) > 0)),
    chorus_effect(m))
  models[[sp]] <- m
}

results <- bind_rows(results)
results$p_adj <- p.adjust(results$p_value, method = "BH")   # kept for reviewers

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
results$common_name <- common_names[results$species]
results <- arrange(results, ratio)

print(results, n = Inf, width = Inf)
write.csv(results, "results/table1_species_ratios.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# MODEL CHECKS (DHARMa)
# -----------------------------------------------------------------------------

checks <- list()
for (sp in names(models)) {
  sim <- simulateResiduals(models[[sp]], n = 1000)
  
  png(paste0("results/model_checks/", sp, ".png"), width = 1000, height = 500)
  plot(sim, title = sp)
  dev.off()
  
  checks[[sp]] <- tibble(
    species          = sp,
    p_uniformity     = testUniformity(sim, plot = FALSE)$p.value,
    p_dispersion     = testDispersion(sim, plot = FALSE)$p.value,
    p_zero_inflation = testZeroInflation(sim, plot = FALSE)$p.value
  )
}
checks <- bind_rows(checks)
print(checks, n = Inf)
write.csv(checks, "results/tableS1_model_checks.csv", row.names = FALSE)


# -----------------------------------------------------------------------------
# FIGURE 1: Chorus/SM2 ratio per species
# -----------------------------------------------------------------------------
# Purple/orange = direction of a significant difference (p < 0.05)

pd <- results
pd$effect <- case_when(
  pd$p_value < 0.05 & pd$ratio > 1 ~ "Chorus higher",
  pd$p_value < 0.05 & pd$ratio < 1 ~ "Chorus lower",
  TRUE                             ~ "No difference"
)
pd$species <- factor(pd$species, levels = pd$species)   # sorted by ratio

fig1 <- ggplot(pd, aes(x = ratio, y = species, colour = effect)) +
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
print(fig1)
ggsave("results/fig1_species_ratios.png", fig1, width = 6.5, height = 6, dpi = 600)
ggsave("results/fig1_species_ratios.pdf", fig1, width = 6.5, height = 6)
