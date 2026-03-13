
library(DHARMa)
library(glmmTMB)
# =============================================================================
# CHECK OVERDISPERSION
# =============================================================================
# Before choosing our model, we check whether counts are overdispersed.
# We fit a simple Poisson GLM and check if residual deviance >> residual df.
# A ratio much larger than 1 means overdispersion — use negative binomial.
# We do this on all species pooled together just to get a quick read.

m_pois_check <- glm(n_calls ~ detector_type + spp + site_id,
                    data   = night_counts,
                    family = poisson)

dispersion_ratio <- deviance(m_pois_check) / df.residual(m_pois_check)
cat("\nOverdispersion ratio:", round(dispersion_ratio, 1), "\n")
cat("(If >> 1, data are overdispersed and we should use negative binomial)\n")


# =============================================================================
# FIT ONE NEGATIVE BINOMIAL MODEL PER SPECIES
# =============================================================================
# The model is:
#   n_calls ~ detector_type + (1 | site_id)
#
# The random effect (1 | site_id) is important because:
#   1. The same site is measured many nights — those nights are not independent.
#   2. Sites differ hugely in overall bat activity (busy site vs quiet site).
#   The random effect accounts for both of these things.
#
# Because Other_FS is the reference level, the detector_type coefficient tells
# us how Anabat_Chorus compares. A positive coefficient = Chorus detects more.
# We exponentiate (exp()) the coefficient to get a ratio: Chorus / Other_FS.
# A ratio of 1.5 means Chorus detects 50% more passes than Other_FS.
# A ratio of 0.7 means Chorus detects 30% fewer passes than Other_FS.
#
# use glmer.nb() which fits a negative binomial model with random effects.

# ---- MYLU: Little brown myotis ----------------------------------------------
dat_MYLU <- filter(night_counts, spp == "MYLU")
m_MYLU <- glmer.nb(n_calls ~ detector_type + (1|site_id), data = dat_MYLU)
summary(m_MYLU)
res_MYLU <- simulateResiduals(fittedModel = m_MYLU)
plot(res_MYLU)


# Pull out the key numbers
coef_MYLU   <- fixef(m_MYLU)["detector_typeAnabat_Chorus"]
se_MYLU     <- sqrt(diag(vcov(m_MYLU)))["detector_typeAnabat_Chorus"]
ratio_MYLU  <- exp(coef_MYLU)
ci_lo_MYLU  <- exp(coef_MYLU - 1.96 * se_MYLU)
ci_hi_MYLU  <- exp(coef_MYLU + 1.96 * se_MYLU)
p_MYLU      <- summary(m_MYLU)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYLU — detection ratio (Chorus / Other_FS):", round(ratio_MYLU, 3),
    "  95% CI:", round(ci_lo_MYLU, 3), "-", round(ci_hi_MYLU, 3),
    "  p =", round(p_MYLU, 4), "\n")

sim <- simulateResiduals(m_MYLU, plot = TRUE)
# Left panel: QQ plot — points should fall on diagonal
# Right panel: residuals vs predicted — should be flat, no pattern

testZeroInflation(sim)
# p < 0.05 here means more zeros than the NB model predicts — 
# you'd want to switch to zero-inflated NB in glmmTMB:
# glmmTMB(n_calls ~ detector_type + (1|site_id), 
#         ziformula = ~1, family = nbinom2)

testDispersion(sim)
# Should be non-significant if NB handled overdispersion adequately


# ---- LANO: Silver-haired bat ------------------------------------------------

dat_LANO <- filter(night_counts, spp == "LANO")

m_LANO <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_LANO)
summary(m_LANO)

coef_LANO   <- fixef(m_LANO)["detector_typeAnabat_Chorus"]
se_LANO     <- sqrt(diag(vcov(m_LANO)))["detector_typeAnabat_Chorus"]
ratio_LANO  <- exp(coef_LANO)
ci_lo_LANO  <- exp(coef_LANO - 1.96 * se_LANO)
ci_hi_LANO  <- exp(coef_LANO + 1.96 * se_LANO)
p_LANO      <- summary(m_LANO)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nLANO — detection ratio (Chorus / Other_FS):", round(ratio_LANO, 3),
    "  95% CI:", round(ci_lo_LANO, 3), "-", round(ci_hi_LANO, 3),
    "  p =", round(p_LANO, 4), "\n")

# Residual checks 
sim <- simulateResiduals(m_LANO, plot = TRUE)

testZeroInflation(sim)

testDispersion(sim)

# ---- MYCI: Western small-footed myotis --------------------------------------

dat_MYCI <- filter(night_counts, spp == "MYCI")

m_MYCI <- glmmTMB(n_calls ~ detector_type + (1 | site_id), ziformula = ~1, family = nbinom2, data = dat_MYCI)
summary(m_MYCI)

coef_MYCI   <- fixef(m_MYCI)["detector_typeAnabat_Chorus"]
se_MYCI     <- sqrt(diag(vcov(m_MYCI)))["detector_typeAnabat_Chorus"]

ratio_MYCI  <- exp(coef_MYCI)
ci_lo_MYCI  <- exp(coef_MYCI - 1.96 * se_MYCI)
ci_hi_MYCI  <- exp(coef_MYCI + 1.96 * se_MYCI)
p_MYCI      <- summary(m_MYCI)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYCI — detection ratio (Chorus / Other_FS):", round(ratio_MYCI, 3),
    "  95% CI:", round(ci_lo_MYCI, 3), "-", round(ci_hi_MYCI, 3),
    "  p =", round(p_MYCI, 4), "\n")

# Residual checks 
sim_zi <- simulateResiduals(m_MYCI, plot = TRUE)

testZeroInflation(sim_zi)

testDispersion(sim)

# ---- EPFU: Big brown bat ----------------------------------------------------

dat_EPFU <- filter(night_counts, spp == "EPFU")

m_EPFU <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_EPFU)
summary(m_EPFU)

coef_EPFU   <- fixef(m_EPFU)["detector_typeAnabat_Chorus"]
se_EPFU     <- sqrt(diag(vcov(m_EPFU)))["detector_typeAnabat_Chorus"]
ratio_EPFU  <- exp(coef_EPFU)
ci_lo_EPFU  <- exp(coef_EPFU - 1.96 * se_EPFU)
ci_hi_EPFU  <- exp(coef_EPFU + 1.96 * se_EPFU)
p_EPFU      <- summary(m_EPFU)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nEPFU — detection ratio (Chorus / Other_FS):", round(ratio_EPFU, 3),
    "  95% CI:", round(ci_lo_EPFU, 3), "-", round(ci_hi_EPFU, 3),
    "  p =", round(p_EPFU, 4), "\n")


# ---- MYCA: California myotis ------------------------------------------------

dat_MYCA <- filter(night_counts, spp == "MYCA")

m_MYCA <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYCA)
summary(m_MYCA)

coef_MYCA   <- fixef(m_MYCA)["detector_typeAnabat_Chorus"]
se_MYCA     <- sqrt(diag(vcov(m_MYCA)))["detector_typeAnabat_Chorus"]
ratio_MYCA  <- exp(coef_MYCA)
ci_lo_MYCA  <- exp(coef_MYCA - 1.96 * se_MYCA)
ci_hi_MYCA  <- exp(coef_MYCA + 1.96 * se_MYCA)
p_MYCA      <- summary(m_MYCA)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYCA — detection ratio (Chorus / Other_FS):", round(ratio_MYCA, 3),
    "  95% CI:", round(ci_lo_MYCA, 3), "-", round(ci_hi_MYCA, 3),
    "  p =", round(p_MYCA, 4), "\n")


# ---- LACI: Hoary bat --------------------------------------------------------

dat_LACI <- filter(night_counts, spp == "LACI")

m_LACI <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_LACI)
summary(m_LACI)

coef_LACI   <- fixef(m_LACI)["detector_typeAnabat_Chorus"]
se_LACI     <- sqrt(diag(vcov(m_LACI)))["detector_typeAnabat_Chorus"]
ratio_LACI  <- exp(coef_LACI)
ci_lo_LACI  <- exp(coef_LACI - 1.96 * se_LACI)
ci_hi_LACI  <- exp(coef_LACI + 1.96 * se_LACI)
p_LACI      <- summary(m_LACI)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nLACI — detection ratio (Chorus / Other_FS):", round(ratio_LACI, 3),
    "  95% CI:", round(ci_lo_LACI, 3), "-", round(ci_hi_LACI, 3),
    "  p =", round(p_LACI, 4), "\n")


# ---- MYYU: Yuma myotis ------------------------------------------------------

dat_MYYU <- filter(night_counts, spp == "MYYU")

m_MYYU <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYYU)
summary(m_MYYU)

coef_MYYU   <- fixef(m_MYYU)["detector_typeAnabat_Chorus"]
se_MYYU     <- sqrt(diag(vcov(m_MYYU)))["detector_typeAnabat_Chorus"]
ratio_MYYU  <- exp(coef_MYYU)
ci_lo_MYYU  <- exp(coef_MYYU - 1.96 * se_MYYU)
ci_hi_MYYU  <- exp(coef_MYYU + 1.96 * se_MYYU)
p_MYYU      <- summary(m_MYYU)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYYU — detection ratio (Chorus / Other_FS):", round(ratio_MYYU, 3),
    "  95% CI:", round(ci_lo_MYYU, 3), "-", round(ci_hi_MYYU, 3),
    "  p =", round(p_MYYU, 4), "\n")


# ---- MYVO: Long-legged myotis -----------------------------------------------

dat_MYVO <- filter(night_counts, spp == "MYVO")

m_MYVO <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYVO)
summary(m_MYVO)

coef_MYVO   <- fixef(m_MYVO)["detector_typeAnabat_Chorus"]
se_MYVO     <- sqrt(diag(vcov(m_MYVO)))["detector_typeAnabat_Chorus"]
ratio_MYVO  <- exp(coef_MYVO)
ci_lo_MYVO  <- exp(coef_MYVO - 1.96 * se_MYVO)
ci_hi_MYVO  <- exp(coef_MYVO + 1.96 * se_MYVO)
p_MYVO      <- summary(m_MYVO)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYVO — detection ratio (Chorus / Other_FS):", round(ratio_MYVO, 3),
    "  95% CI:", round(ci_lo_MYVO, 3), "-", round(ci_hi_MYVO, 3),
    "  p =", round(p_MYVO, 4), "\n")


# ---- MYEV -------------------------------------------------------
# Note: rarer species — model may give warnings about convergence.
# If glmer.nb fails, switch to glm.nb (see note at bottom of script).

dat_MYEV <- filter(night_counts, spp == "MYEV")

m_MYEV <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYEV)
summary(m_MYEV)

coef_MYEV   <- fixef(m_MYEV)["detector_typeAnabat_Chorus"]
se_MYEV     <- sqrt(diag(vcov(m_MYEV)))["detector_typeAnabat_Chorus"]
ratio_MYEV  <- exp(coef_MYEV)
ci_lo_MYEV  <- exp(coef_MYEV - 1.96 * se_MYEV)
ci_hi_MYEV  <- exp(coef_MYEV + 1.96 * se_MYEV)
p_MYEV      <- summary(m_MYEV)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYEV — detection ratio (Chorus / Other_FS):", round(ratio_MYEV, 3),
    "  95% CI:", round(ci_lo_MYEV, 3), "-", round(ci_hi_MYEV, 3),
    "  p =", round(p_MYEV, 4), "\n")

# ---- MYSE --------------------------------------

dat_MYSE <- filter(night_counts, spp == "MYSE")

m_MYSE <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYSE)
summary(m_MYSE)

coef_MYSE   <- fixef(m_MYSE)["detector_typeAnabat_Chorus"]
se_MYSE     <- sqrt(diag(vcov(m_MYSE)))["detector_typeAnabat_Chorus"]
ratio_MYSE  <- exp(coef_MYSE)
ci_lo_MYSE  <- exp(coef_MYSE - 1.96 * se_MYSE)
ci_hi_MYSE  <- exp(coef_MYSE + 1.96 * se_MYSE)
p_MYSE      <- summary(m_MYSE)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYSE — detection ratio (Chorus / Other_FS):", round(ratio_MYSE, 3),
    "  95% CI:", round(ci_lo_MYSE, 3), "-", round(ci_hi_MYSE, 3),
    "  p =", round(p_MYSE, 4), "\n")


# ---- MYTH ----------------------------------------------------

dat_MYTH <- filter(night_counts, spp == "MYTH")

m_MYTH <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYTH)
summary(m_MYTH)

coef_MYTH   <- fixef(m_MYTH)["detector_typeAnabat_Chorus"]
se_MYTH     <- sqrt(diag(vcov(m_MYTH)))["detector_typeAnabat_Chorus"]
ratio_MYTH  <- exp(coef_MYTH)
ci_lo_MYTH  <- exp(coef_MYTH - 1.96 * se_MYTH)
ci_hi_MYTH  <- exp(coef_MYTH + 1.96 * se_MYTH)
p_MYTH      <- summary(m_MYTH)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYTH — detection ratio (Chorus / Other_FS):", round(ratio_MYTH, 3),
    "  95% CI:", round(ci_lo_MYTH, 3), "-", round(ci_hi_MYTH, 3),
    "  p =", round(p_MYTH, 4), "\n")


# ---- LABO ------------------------------------------------

dat_LABO <- filter(night_counts, spp == "LABO")

m_LABO <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_LABO)
summary(m_LABO)

coef_LABO   <- fixef(m_LABO)["detector_typeAnabat_Chorus"]
se_LABO     <- sqrt(diag(vcov(m_LABO)))["detector_typeAnabat_Chorus"]
ratio_LABO  <- exp(coef_LABO)
ci_lo_LABO  <- exp(coef_LABO - 1.96 * se_LABO)
ci_hi_LABO  <- exp(coef_LABO + 1.96 * se_LABO)
p_LABO      <- summary(m_LABO)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nLABO — detection ratio (Chorus / Other_FS):", round(ratio_LABO, 3),
    "  95% CI:", round(ci_lo_LABO, 3), "-", round(ci_hi_LABO, 3),
    "  p =", round(p_LABO, 4), "\n")


# ---- ANPA --------------------------------------------------------

dat_ANPA <- filter(night_counts, spp == "ANPA")

m_ANPA <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_ANPA)
summary(m_ANPA)

coef_ANPA   <- fixef(m_ANPA)["detector_typeAnabat_Chorus"]
se_ANPA     <- sqrt(diag(vcov(m_ANPA)))["detector_typeAnabat_Chorus"]
ratio_ANPA  <- exp(coef_ANPA)
ci_lo_ANPA  <- exp(coef_ANPA - 1.96 * se_ANPA)
ci_hi_ANPA  <- exp(coef_ANPA + 1.96 * se_ANPA)
p_ANPA      <- summary(m_ANPA)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nANPA — detection ratio (Chorus / Other_FS):", round(ratio_ANPA, 3),
    "  95% CI:", round(ci_lo_ANPA, 3), "-", round(ci_hi_ANPA, 3),
    "  p =", round(p_ANPA, 4), "\n")


# ---- MYYU ------------------------------------------------------

dat_MYYU <- filter(night_counts, spp == "MYYU")

m_MYYU <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_MYYU)
summary(m_MYYU)

coef_MYYU   <- fixef(m_MYYU)["detector_typeAnabat_Chorus"]
se_MYYU     <- sqrt(diag(vcov(m_MYYU)))["detector_typeAnabat_Chorus"]
ratio_MYYU  <- exp(coef_MYYU)
ci_lo_MYYU  <- exp(coef_MYYU - 1.96 * se_MYYU)
ci_hi_MYYU  <- exp(coef_MYYU + 1.96 * se_MYYU)
p_MYYU      <- summary(m_MYYU)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nMYYU — detection ratio (Chorus / Other_FS):", round(ratio_MYYU, 3),
    "  95% CI:", round(ci_lo_MYYU, 3), "-", round(ci_hi_MYYU, 3),
    "  p =", round(p_MYYU, 4), "\n")


# ---- PAHE -----------------------------------------------

dat_PAHE <- filter(night_counts, spp == "PAHE")

m_PAHE <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_PAHE)
summary(m_PAHE)

coef_PAHE   <- fixef(m_PAHE)["detector_typeAnabat_Chorus"]
se_PAHE     <- sqrt(diag(vcov(m_PAHE)))["detector_typeAnabat_Chorus"]
ratio_PAHE  <- exp(coef_PAHE)
ci_lo_PAHE  <- exp(coef_PAHE - 1.96 * se_PAHE)
ci_hi_PAHE  <- exp(coef_PAHE + 1.96 * se_PAHE)
p_PAHE      <- summary(m_PAHE)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nPAHE — detection ratio (Chorus / Other_FS):", round(ratio_PAHE, 3),
    "  95% CI:", round(ci_lo_PAHE, 3), "-", round(ci_hi_PAHE, 3),
    "  p =", round(p_PAHE, 4), "\n")


# ---- COTO -------------------------------------------------------
# Note: rarer species — model may give warnings about convergence.
# If glmer.nb fails, switch to glm.nb (see note at bottom of script).

dat_COTO <- filter(night_counts, spp == "COTO")

m_COTO <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_COTO)
summary(m_COTO)

coef_COTO   <- fixef(m_COTO)["detector_typeAnabat_Chorus"]
se_COTO     <- sqrt(diag(vcov(m_COTO)))["detector_typeAnabat_Chorus"]
ratio_COTO  <- exp(coef_COTO)
ci_lo_COTO  <- exp(coef_COTO - 1.96 * se_COTO)
ci_hi_COTO  <- exp(coef_COTO + 1.96 * se_COTO)
p_COTO      <- summary(m_COTO)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nCOTO — detection ratio (Chorus / Other_FS):", round(ratio_COTO, 3),
    "  95% CI:", round(ci_lo_COTO, 3), "-", round(ci_hi_COTO, 3),
    "  p =", round(p_COTO, 4), "\n")

# ---- EUMA -------------------------------------------------------
# Note: rarer species — model may give warnings about convergence.
# If glmer.nb fails, switch to glm.nb (see note at bottom of script).

dat_EUMA <- filter(night_counts, spp == "EUMA")

m_EUMA <- glmer.nb(n_calls ~ detector_type + (1 | site_id), data = dat_EUMA)
summary(m_EUMA)

coef_EUMA   <- fixef(m_EUMA)["detector_typeAnabat_Chorus"]
se_EUMA     <- sqrt(diag(vcov(m_EUMA)))["detector_typeAnabat_Chorus"]
ratio_EUMA  <- exp(coef_EUMA)
ci_lo_EUMA  <- exp(coef_EUMA - 1.96 * se_EUMA)
ci_hi_EUMA  <- exp(coef_EUMA + 1.96 * se_EUMA)
p_EUMA      <- summary(m_EUMA)$coefficients["detector_typeAnabat_Chorus", "Pr(>|z|)"]

cat("\nEUMA — detection ratio (Chorus / Other_FS):", round(ratio_EUMA, 3),
    "  95% CI:", round(ci_lo_EUMA, 3), "-", round(ci_hi_EUMA, 3),
    "  p =", round(p_EUMA, 4), "\n")


# =============================================================================
# COMPILE RESULTS — all 16 species with actual values
# =============================================================================

results <- data.frame(
  species = c("MYLU","LANO","MYCI","EPFU","MYCA","LACI","MYYU","MYVO",
              "MYEV","MYSE","MYTH","LABO","ANPA","PAHE","COTO","EUMA"),
  
  ratio   = c(1.077, 1.027, 1.513, 1.102, 1.410, 1.013, 1.092, 1.911,
              1.305, 1.944, 1.006, 0.594, 1.712, 0.309, 1.073, 0.782),
  
  ci_lo   = c(0.873, 0.815, 1.099, 0.854, 1.120, 0.774, 0.864, 1.427,
              0.958, 1.344, 0.688, 0.418, 0.876, 0.128, 0.736, 0.438),
  
  ci_hi   = c(1.328, 1.294, 2.083, 1.422, 1.776, 1.325, 1.382, 2.561,
              1.778, 2.810, 1.473, 0.843, 3.348, 0.751, 1.566, 1.398),
  
  p_value = c(0.4913, 0.8209, 0.0112, 0.4534, 0.0035, 0.9271, 0.4616, 0.0000,
              0.0917, 0.0004, 0.9737, 0.0036, 0.1159, 0.0095, 0.7131, 0.4072),
  
  # Note on EUMA: model failed to converge — exclude from inference
  # Fill in EUMA values above once model is run or leave as NA
  
  # Convergence warnings on EPFU and MYEV — results treated with caution (flagged below)
  converged = c(TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE,
                FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, NA)
)

# Add common names as a reference column
results$common_name <- c(
  "Little brown myotis",       # MYLU
  "Silver-haired bat",         # LANO
  "W. small-footed myotis",    # MYCI
  "Big brown bat",             # EPFU
  "California myotis",         # MYCA
  "Hoary bat",                 # LACI
  "Yuma myotis",               # MYYU
  "Long-legged myotis",        # MYVO
  "Long-eared myotis",         # MYEV
  "N. long-eared myotis",      # MYSE
  "Fringed myotis",            # MYTH
  "Silver-haired bat",         # LABO
  "Pallid bat",                # ANPA
  "Canyon bat",                # PAHE
  "Townsend's big-eared bat",  # COTO
  "Spotted bat"                # EUMA
)

# Add approximate peak call frequency (kHz) — useful for interpreting patterns
# These are approximate characteristic frequencies from the literature
results$peak_freq_kHz <- c(
  45,   # MYLU  — high freq Myotis
  27,   # LANO  — low-mid freq
  55,   # MYCI  — high freq Myotis
  28,   # EPFU  — low-mid freq
  55,   # MYCA  — high freq Myotis
  20,   # LACI  — very low freq
  45,   # MYYU  — high freq Myotis
  40,   # MYVO  — high freq Myotis
  35,   # MYEV  — mid freq Myotis
  45,   # MYSE  — high freq Myotis
  25,   # MYTH  — mid freq Myotis (lower than many Myotis)
  27,   # LABO  — low-mid freq
  35,   # ANPA  — mid freq
  50,   # PAHE  — high freq
  30,   # COTO  — mid freq
  11    # EUMA  — very low freq (spotted bat)
)

# p_adj: not correcting across species (independent datasets, separate inferential units)
# but we add it here anyway for transparency — students can compare corrected vs uncorrected
results$p_adj <- p.adjust(results$p_value, method = "BH")

# Significance based on RAW p-value (our chosen approach — see methods justification)
results$sig <- "ns"
results$sig[results$p_value < 0.1]   <- "."
results$sig[results$p_value < 0.05]  <- "*"
results$sig[results$p_value < 0.01]  <- "**"
results$sig[results$p_value < 0.001] <- "***"

# Flag convergence warnings
results$sig[!results$converged & !is.na(results$converged)] <-
  paste0(results$sig[!results$converged & !is.na(results$converged)], "†")

# Sort by ratio
results <- results[order(results$ratio, na.last = TRUE), ]

cat("\n=== RESULTS SUMMARY ===\n")
cat("ratio > 1: Chorus detects MORE than Other_FS\n")
cat("ratio < 1: Chorus detects LESS than Other_FS\n")
cat("† = convergence warning, treat with caution\n\n")

print(
  results[, c("species","common_name","ratio","ci_lo","ci_hi",
              "p_value","p_adj","sig","peak_freq_kHz")],
  row.names = FALSE, digits = 3
)

write.csv(results, "results_NB_models.csv", row.names = FALSE)
cat("\nSaved: results_NB_models.csv\n")


# =============================================================================
# FOREST PLOT — coloured by direction of effect (CI-based, not p-value)
# =============================================================================
# Colouring logic: we colour by whether the 95% CI lies entirely above 1
# (consistently more), entirely below 1 (consistently less), or overlaps 1
# (uncertain direction). This is more honest than colouring by p-value
# because it communicates the consistency and direction of the effect directly.
# Both "Chorus detects more" and "Chorus detects less" use blue tones
# (colourblind-safe) — darker for consistent effects, grey for uncertain.

# Remove species with failed models
plot_dat <- results[!is.na(results$ratio), ]
plot_dat$species <- factor(plot_dat$species, levels = plot_dat$species)

# Classify by CI direction — no reference to p-values
plot_dat$effect <- "CI overlaps 1 (uncertain)"
plot_dat$effect[plot_dat$ci_lo > 1] <- "CI entirely above 1 (Chorus detects more)"
plot_dat$effect[plot_dat$ci_hi < 1] <- "CI entirely below 1 (Chorus detects less)"

ggplot(plot_dat, aes(x = ratio, y = species, colour = effect)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.35, linewidth = 0.7) +
  geom_point(size = 3.5) +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4),
    labels = c("0.25\n(4x less)", "0.5\n(2x less)", "1\n(equal)",
               "2\n(2x more)", "4\n(4x more)")
  ) +
  scale_colour_manual(
    values = c(
      "CI entirely above 1 (Chorus detects more)" = "#5D3A9B",   # dark blue
      "CI entirely below 1 (Chorus detects less)" = "#E66100",   # light blue
      "CI overlaps 1 (uncertain)"                 = "grey65"
    )
  ) +
  # Annotate peak frequency on right side
  geom_text(aes(x = 4.5, label = paste0(peak_freq_kHz, " kHz")),
            colour = "grey40", size = 3, hjust = 0) +
  coord_cartesian(clip = "off") +
  labs(
    x      = "Detection ratio: Chorus / SM2-Swift  (log scale)",
    y      = NULL,
    colour = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position      = "bottom",
    legend.text          = element_text(size = 10),
    panel.grid.major     = element_blank(),
    panel.grid.minor     = element_blank(),
    plot.margin          = margin(5, 65, 5, 5)   # extra right margin for kHz labels
  )

ggsave("plot_NB_results.png", width = 9, height = 7, dpi = 150)



