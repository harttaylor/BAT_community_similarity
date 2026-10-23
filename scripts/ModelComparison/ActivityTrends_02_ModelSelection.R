library(glmmTMB)
library(ggplot2)
library(GGally)
library(ggmap)
library(insight)
library(kableExtra)
library(lmerTest)
library(tidyverse)
library(ggforce)
library(huxtable)
library(AICcmodavg)
library(glmpathcr)
library(glmmTMB)
library(dplyr)
library(interactions)

#Need to run the TrendDataSetUp.R code first with the raw data frame. 
#That code will do data checks, and will give you the data frame that can be used to run this code

#Identify which species should be excluded from the analysis due to low sample numbers.
#excluded species in stationary
SpeciesExclude = c('PAHE')
#excluded species in transects
SpeciesExcludeTrans = c('PAHE')

##Model Selection for stationary models in report
# ── 1.  Define candidate variable sets per species ──────────────────────────
# Each species gets a named list of candidate formulas.
# Add as many candidates as you like; names will appear in the AIC table.

candidate_vars <- list(
  ANPA = list(
    m1 = "Nightly_Min_Temp",
    m2 = "Nightly_Mean_Temp"
  ),
  COTO = list(
    m1 = "dd,Nightly_Mean_Temp",
    m2 = "Water_Nearby,Nightly_Min_Temp, dd",
    m3 = "dd,Nightly_Min_Temp"
  ),
  EPFU = list(
    m1 = "moon*Distance_to_Clutter__m_,Nightly_Mean_Temp",
    m2 = "Nightly_Mean_Windsp,moon*Distance_to_Clutter__m_,Nightly_Min_Temp",
    m3 = "Nightly_Mean_Temp"
  ),
  EUMA = list(
    m1 = "Nightly_Min_Temp",
    m2 = "Nightly_Mean_Temp"
  ),
  LABO = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Mean_Temp",
    m2 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Min_Temp",
    m3 = "Water_Nearby,Nightly_Mean_Temp"
  ),
  LACI = list(
    m1 = "moon*Distance_to_Clutter__m_,dd,Nightly_Mean_Temp,I(YearS^2)",
    m2 = "moon*Distance_to_Clutter__m_,ddclim,Nightly_Min_Temp",
    m3 = "dd,Nightly_Mean_Temp"
  ),
  LANO = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Mean_Temp",
    m2 = "Nightly_Mean_RH,Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Min_Temp",
    m3 = "Nightly_Min_Temp"
  ),
  MYCA = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Mean_Temp",
    m2 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Min_Temp",
    m3 = "Water_Nearby,Nightly_Mean_Temp"
  ),
  MYCI = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Mean_Temp",
    m2 = "Nightly_Min_Temp,Nightly_Mean_Windsp*Distance_to_Clutter__m_",
    m3 = "Nightly_Min_Temp"
  ),
  MYEV = list(
    m1 = "dd,Nightly_Mean_Temp",
    m2 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,dd,Nightly_Min_Temp",
    m3 = "dd,Nightly_Mean_Temp"
  ),
  MYLU = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Mean_Temp",
    m2 = "Water_Nearby,Nightly_Mean_Temp",
    m3 = "Nightly_Mean_Temp"
  ),
  MYSE = list(
    m1 = "Nightly_Mean_RH,moon*Distance_to_Clutter__m_,Water_Nearby",
    m2 = "Nightly_Mean_RH,Water_Nearby",
    m3 = ""
  ),
  MYTH = list(
    m1 = "dd"
  ),
  MYVO = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Min_Temp",
    m2 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,Nightly_Mean_Temp",
    m3 = "Nightly_Mean_Temp"
  ),
  MYYU = list(
    m1 = "Nightly_Mean_Windsp*Distance_to_Clutter__m_,Water_Nearby,dd",
    m2 = "dd"
  )
)

# ── 2.  Helper: fit one glmmTMB model given a var string ────────────────────

fit_one_model <- function(cell.summary, var_string) {
  
  # Parse variables and drop rows with NAs in any covariate
  vars_vec <- if (nchar(trimws(var_string)) == 0) {
    character(0)
  } else {
    strsplit(var_string, ",")[[1]]
  }
  
  check <- unlist(sapply(strsplit(vars_vec, split = "*", fixed = TRUE), `[`, simplify = "TRUE"))
  check <- check[!is.na(check) & nchar(trimws(check)) > 0]
  
  # Strip any I(...) wrappers before trying to use as column names
  check_cols <- gsub("I\\((.*)\\)", "\\1", check)
  check_cols <- gsub("\\^.*", "", check_cols)          # remove ^2 etc.
  check_cols <- trimws(check_cols)
  check_cols <- check_cols[check_cols %in% names(cell.summary)]
  
  data_clean <- tidyr::drop_na(cell.summary, any_of(check_cols))
  
  if (nrow(data_clean) < 5) return(NULL)   # not enough data to fit
  
  # Build formula
  rhs_covs <- if (length(vars_vec) > 0) paste("+", paste(vars_vec, collapse = " + ")) else ""
  form <- paste0(
    "total.detect ~ YearS + log(quad.nights) + (1|GRTS_Cell_ID/Quadrant)",
    rhs_covs
  )
  
  tryCatch(
    glmmTMB::glmmTMB(
      formula = formula(form),
      family  = nbinom2(link = "log"),
      data    = data_clean,
      se      = TRUE,
      verbose = FALSE
    ),
    error   = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("  WARNING:", conditionMessage(w), "\n"); NULL }
  )
}


# ── 3.  Main loop: fit all candidates per species, rank by AIC ───────────────

sploc <- read.csv("Data/Raw/ActivityModels/sploc.csv",check.names = FALSE, row.names = 1)
bc.fits.cov.singlet <- plyr::dlply(bat.data.long, "SpeciesGroup", function(x) {
  
  sp <- x$SpeciesGroup[1]
  cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Species:", sp, "\n")
  
  # ── Retrieve candidate var sets for this species
  sp_candidates <- candidate_vars[[sp]]
  if (is.null(sp_candidates)) {
    cat("  No candidate variables defined — skipping\n")
    return(NULL)
  }
  
  # ── Filter to surveyed GRTS cells
  spgr_row  <- sploc[sp, ]
  grts_list <- names(spgr_row)[spgr_row == TRUE]
  x         <- x[x$GRTS_Cell_ID %in% grts_list, ]
  cat("  GRTS cells:", n_distinct(x$GRTS_Cell_ID), "\n")
  
  # ── Summarise to cell × year × quadrant
  cell.summary <- plyr::ddply(
    x, c("Year", "GRTS_Cell_ID", "Quadrant"), plyr::summarize,
    YearS              = Year[1] - 2016,
    total.detect       = sum(SpeciesSingleton),
    quad.nights        = length(SpeciesSingleton),
    Nightly_Mean_Temp  = mean(Nightly_Mean_Temp,  na.rm = TRUE),
    Nightly_Min_Temp   = mean(Nightly_Min_Temp,   na.rm = TRUE),
    Nightly_Max_Temp   = mean(Nightly_Max_Temp,   na.rm = TRUE),
    Nightly_Mean_RH    = mean(Nightly_Mean_RH,    na.rm = TRUE),
    Nightly_Mean_Windsp= mean(Nightly_Mean_Windsp,na.rm = TRUE),
    dd                 = mean(dd,                  na.rm = TRUE),
    jNight             = mean(jNight),
    moon_fraction_above_horizon   = mean(moon_fraction_above_horizon,na.rm = TRUE),
    moon               = mean(moon_illumination, na.rm = TRUE),
    Feature_Sampled    = unique(Feature_Sampled,   na.rm = TRUE)[1],
    Habitat_Type       = unique(Habitat_Type,      na.rm = TRUE)[1],
    Distance_to_Clutter__m_        = median(`Distance_to_Clutter_(m)` , na.rm = TRUE),
    Percent_Clutter    = mean(Percent_Clutter,     na.rm = TRUE),
    Water_Nearby       = median(Water_Nearby,      na.rm = TRUE)
  )
  cell.summary$countperquad <- cell.summary$total.detect / cell.summary$quad.nights
  
  # LACI outlier trim
  if (sp == "LACI") {
    cell.summary <- cell.summary[cell.summary$countperquad < 75, ]
  }
  
  # ── Fit every candidate model
  cat("  Fitting", length(sp_candidates), "candidate models...\n")
  
  fits <- lapply(names(sp_candidates), function(model_name) {
    cat("  →", model_name, ":", sp_candidates[[model_name]], "\n")
    fit_one_model(cell.summary, sp_candidates[[model_name]])
  })
  names(fits) <- names(sp_candidates)
  
  # ── Build AIC comparison table
  aic_table <- do.call(rbind, lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    if (is.null(fit)) {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = NA, converged = FALSE, stringsAsFactors = FALSE)
    } else {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = AIC(fit), converged = TRUE, stringsAsFactors = FALSE)
    }
  }))
  
  aic_table <- aic_table[order(aic_table$AIC), ]
  aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC, na.rm = TRUE)
  
  cat("\n  AIC comparison for", sp, ":\n")
  print(aic_table, row.names = FALSE)
  
  best_model_name <- aic_table$model[1]
  cat("\n  ✔ Best model:", best_model_name, "(AIC =", round(aic_table$AIC[1], 2), ")\n")
  
  # ── Return everything
  list(
    fits      = fits,
    aic_table = aic_table,
    best_fit  = fits[[best_model_name]],
    best_vars = sp_candidates[[best_model_name]]
  )
})


# ── 4.  Summary AIC table across all species ─────────────────────────────────

all_aic <- do.call(rbind, lapply(names(bc.fits.cov.singlet), function(sp) {
  res <- bc.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  cbind(species = sp, res$aic_table)
}))

cat("\n\n═══════════ FULL AIC SUMMARY ═══════════\n")
print(all_aic, row.names = FALSE)


# ── 5.  Extract best fits only (same structure as your original output) ───────

best_fits <- lapply(bc.fits.cov.singlet, `[[`, "best_fit")

bc.cov.slope.summary.singlet <- plyr::ldply(best_fits, function(x){
  if(is.null(x$fit))return(NULL)
  summary.table <- summary(x)
  #browser()
  slope   <- as.data.frame(summary.table$coefficients$cond[2,, drop=FALSE])
  slope$intercept <- summary.table$coefficients$cond[1,,drop=FALSE]
  slope$n <- nobs(x)
  slope
})


bc.cov.slope.summary.singlet <- plyr::rename(bc.cov.slope.summary.singlet, c("Pr(>|z|)"="p.value","Std. Error"="SE",".id"="SpeciesGroup"))
# bc.cov.slope.summary.singlet$model <- paste(bc.cov.slope.summary.singlet$Var1, bc.cov.slope.summary.singlet$Var2, bc.cov.slope.summary.singlet$Var3, bc.cov.slope.summary.singlet$Var4, bc.cov.slope.summary.singlet$Var5, bc.cov.slope.summary.singlet$Var6, bc.cov.slope.summary.singlet$Var7, sep = " + ")
temptable <- bc.cov.slope.summary.singlet[,c("SpeciesGroup","Estimate","SE","p.value")]

kable(temptable, row.names=FALSE, 
      caption="Estimated trends on logarithmic scale in bat detections for semipool, fullpool and singlet models-SM2",
      col.names=c("Species","Estimate","SE","P-value"),
      digits=c(0,  2,2,4))  %>% 
  add_header_above(c(" "=1, "AutoID BC W Covariates"=3)) %>%
  column_spec(column=c(1),       width="12cm") %>%
  column_spec(column=c(2:4),       width="2cm") %>%
  kable_styling("bordered",position = "center", full_width=FALSE, latex_options = "HOLD_position")  ####%$%$


write.csv(bc.cov.slope.summary.singlet, file="Data/ActivityModels/bc.cov.estimates.singlet.SM2.csv", row.names=FALSE)


##-------------------------Make the vars lists best on the best fit models---------------------
varlists <- do.call(rbind, lapply(names(bc.fits.cov.singlet), function(sp) {
  res <- bc.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  data.frame(
    vars = res$best_vars,
    sp   = sp,
    stringsAsFactors = FALSE
  )
}))

write.csv(varlists,"Data/ActivityModels/ModelVariables-SM2.csv")

##Model Selection for transect models in report
# ── 1. Define candidate variable sets ──────────────────────────────────────

candidate_vars_transect <- list(
  AUTO = list(
    # Single covariates
    m1  = "Nightly_Mean_Temp",
    m2  = "Nightly_Mean_RH",
    m3  = "Nightly_Mean_Windsp",
    m4  = "moon",
    # Two-variable combinations
    m5  = "Nightly_Mean_Temp,Nightly_Mean_RH",
    m6  = "Nightly_Mean_Temp,Nightly_Mean_Windsp",
    m7  = "Nightly_Mean_Temp,moon",
    m8  = "Nightly_Mean_RH,Nightly_Mean_Windsp",
    m9  = "Nightly_Mean_RH,moon",
    m10 = "Nightly_Mean_Windsp,moon",
    # Three-variable combinations
    m11 = "Nightly_Mean_Temp,Nightly_Mean_RH,Nightly_Mean_Windsp",
    m12 = "Nightly_Mean_Temp,Nightly_Mean_RH,moon",
    m13 = "Nightly_Mean_Temp,Nightly_Mean_Windsp,moon",
    m14 = "Nightly_Mean_RH,Nightly_Mean_Windsp,moon",
    # Full model
    m15 = "Nightly_Mean_Temp,Nightly_Mean_RH,Nightly_Mean_Windsp,moon"
  )
)


# ── 2. Helper: fit one transect glmmTMB model given a var string ────────────

fit_one_transect_model <- function(cell.summary, var_string) {
  
  # Parse variables and drop rows with NAs in any covariate
  vars_vec <- if (nchar(trimws(var_string)) == 0) {
    character(0)
  } else {
    strsplit(var_string, ",")[[1]]
  }
  
  check <- unlist(sapply(strsplit(vars_vec, split = "*", fixed = TRUE), `[`, simplify = "TRUE"))
  check <- check[!is.na(check) & nchar(trimws(check)) > 0]
  
  # Strip any I(...) wrappers before trying to use as column names
  check_cols <- gsub("I\\((.*)\\)", "\\1", check)
  check_cols <- gsub("\\^.*", "", check_cols)
  check_cols <- trimws(check_cols)
  check_cols <- check_cols[check_cols %in% names(cell.summary)]
  
  data_clean <- tidyr::drop_na(cell.summary, any_of(check_cols))

  if (nrow(data_clean) < 5) return(NULL)
  
  # Build formula — keeps transect.length and single random effect
  rhs_covs <- if (length(vars_vec) > 0) paste("+", paste(vars_vec, collapse = " + ")) else ""
  form <- paste0(
    "total.detect ~ YearS + log(transect.length) + (1|GRTS_Cell_ID)",
    rhs_covs
  )
  
  tryCatch(
    glmmTMB::glmmTMB(
      formula = formula(form),
      family  = nbinom1(link = "log"),
      data    = data_clean,
      se      = TRUE,
      verbose = FALSE
    ),
    error   = function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL },
    warning = function(w) { cat("  WARNING:", conditionMessage(w), "\n"); NULL }
  )
}


# ── 3. Main loop: fit all candidates per species, rank by AIC ───────────────

bc.auto.transect.fits.cov.singlet <- plyr::dlply(bat.transect.data.long, "SpeciesGroup", function(x) {
  
  sp <- x$SpeciesGroup[1]
  cat("\n\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
  cat("Species:", sp, "\n")
  
  # ── Retrieve candidate var sets (same for all species here)
  sp_candidates <- candidate_vars_transect[["AUTO"]]
  
  # ── Summarise to cell × year
  cell.summary <- plyr::ddply(x, c("Year", "GRTS_Cell_ID"), plyr::summarize,
                              YearS               = Year[1] - 2015,
                              total.detect        = sum(SpeciesSingleton),
                              transect.nights     = length(SpeciesSingleton),
                              transect.length     = mean(Tlength),
                              Nightly_Mean_Temp   = mean(Nightly_Mean_Temp,   na.rm = TRUE),
                              Nightly_Mean_RH     = mean(Nightly_Mean_RH,     na.rm = TRUE),
                              Nightly_Mean_Windsp = mean(Nightly_Mean_Windsp, na.rm = TRUE),
                              moon_fraction_above_horizon            = mean(moon_fraction_above_horizon,na.rm = TRUE),
                              moon      = mean(moon_illumination,      na.rm = TRUE))
  
  cell.summary$GRTS_Cell_ID <- as.factor(cell.summary$GRTS_Cell_ID)
  
  # ── Fit every candidate model
  cat("  Fitting", length(sp_candidates), "candidate models...\n")
  
  fits <- lapply(names(sp_candidates), function(model_name) {
    cat("  →", model_name, ":", sp_candidates[[model_name]], "\n")
    fit_one_transect_model(cell.summary, sp_candidates[[model_name]])
  })
  names(fits) <- names(sp_candidates)
  
  # ── Build AIC comparison table
  aic_table <- do.call(rbind, lapply(names(fits), function(model_name) {
    fit <- fits[[model_name]]
    if (is.null(fit)) {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = NA, converged = FALSE, stringsAsFactors = FALSE)
    } else {
      data.frame(model = model_name, vars = sp_candidates[[model_name]],
                 AIC = AIC(fit), converged = TRUE, stringsAsFactors = FALSE)
    }
  }))
  
  aic_table <- aic_table[order(aic_table$AIC), ]
  aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC, na.rm = TRUE)
  
  cat("\n  AIC comparison for", sp, ":\n")
  print(aic_table, row.names = FALSE)
  
  best_model_name <- aic_table$model[1]
  cat("\n  ✔ Best model:", best_model_name, "(AIC =", round(aic_table$AIC[1], 2), ")\n")
  
  # ── Return everything
  list(
    fits      = fits,
    aic_table = aic_table,
    best_fit  = fits[[best_model_name]],
    best_vars = sp_candidates[[best_model_name]]
  )
})


# ── 4. Summary AIC table across all species ──────────────────────────────────

all_aic_transect <- do.call(rbind, lapply(names(bc.auto.transect.fits.cov.singlet), function(sp) {
  res <- bc.auto.transect.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  cbind(species = sp, res$aic_table)
}))

cat("\n\n═══════════ FULL AIC SUMMARY ═══════════\n")
print(all_aic_transect, row.names = FALSE)


# ── 5. Extract best fits only ─────────────────────────────────────────────────

best_fits_transect <- lapply(bc.auto.transect.fits.cov.singlet, `[[`, "best_fit")

bc.transect.cov.slope.summary.singlet <- plyr::ldply(best_fits_transect, function(x){
  if (is.null(x)) return(NULL)                     # no model fitted
  if (!inherits(x, "glmmTMB")) return(NULL)        # try-error, NA, etc.
  cf <- tryCatch(summary(x)$coefficients$cond, error = function(e) NULL)
  if (is.null(cf) || nrow(cf) < 2) return(NULL)    # summary failed or intercept-only
  slope <- as.data.frame(cf[2, , drop = FALSE])
  slope$term         <- rownames(cf)[2]
  slope$intercept    <- cf[1, "Estimate"]
  slope$intercept.se <- cf[1, "Std. Error"]
  slope
})
#bc.slope.summary
bc.transect.cov.slope.summary.singlet <- plyr::rename(bc.transect.cov.slope.summary.singlet,
                                                      c("Pr(>|z|)"="p.value","Std. Error"="SE",".id"="SpeciesGroup"))


bc.transect.cov.slope.summary.table.singlet <- bc.transect.cov.slope.summary.singlet[,c("SpeciesGroup","Estimate","SE","p.value")]
bc.transect.cov.slope.summary.table.singlet[,"p.value"] <- insight::format_p(bc.transect.cov.slope.summary.table.singlet[,"p.value"])

write.csv(bc.transect.cov.slope.summary.table.singlet, file="Data/Analyzed/ActivityTrends/bc.transect.cov.estimates.csv", row.names=FALSE)


# ── 6. Save variable list for best models ────────────────────────────────────

varlists_transect <- do.call(rbind, lapply(names(bc.auto.transect.fits.cov.singlet), function(sp) {
  res <- bc.auto.transect.fits.cov.singlet[[sp]]
  if (is.null(res)) return(NULL)
  data.frame(
    vars = res$best_vars,
    sp   = sp,
    stringsAsFactors = FALSE
  )
}))

write.csv(varlists_transect, "Data/Analyzed/ActivityTrends/ModelVariables_Transect.csv")


##Running regional models
###-------------------------------Full Regional Models Updated 2025---------------------------------------
# Get all regions (excluding NA)
regions <- unique(bat.data.long$Region[!is.na(bat.data.long$Region)])

# Run models across all regions
all.region.fits <- plyr::dlply(
  bat.data.long[!bat.data.long$Year %in% c(2015) &
                  !is.na(bat.data.long$SpeciesGroup), ],
  c("Region"),
  function(region_data) {
    
    current_region <- region_data$Region[1]
    cat("\n\n===== Fitting models for region:", current_region, "=====\n")
    
    region_data$GRTS_Cell_ID <- as.character(region_data$GRTS_Cell_ID)
    
    plyr::dlply(region_data, c("SpeciesGroup"), function(x) {
      
      cat("\n\nFit cov singleton - Start to fit species groups", x$SpeciesGroup[1], "\n")
      
      # Identify variables of interest
      print(strsplit(varlists[varlists$sp == x$SpeciesGroup[1], "vars"], ",")[[1]])
      vars <- strsplit(varlists[varlists$sp == x$SpeciesGroup[1], "vars"], ",")[[1]]
      
      check <- unlist(sapply(strsplit(vars, split = '*', fixed = TRUE), `[`, simplify = "TRUE"))
      print(check)
      x <- x %>% drop_na(any_of(check))
      
      vars <- append(c(""), vars)
      
      # Filter to GRTS cells where species was surveyed
      spgr_row <- sploc[x$SpeciesGroup[1], ]
      grts_list <- names(spgr_row)[spgr_row == TRUE]
      grts_in_both <- intersect(unique(x$GRTS_Cell_ID), grts_list)
      
      if (length(grts_in_both) > 0) {
        x <- x[x$GRTS_Cell_ID %in% grts_in_both, ]
        cat("Filtered to", length(grts_in_both), "GRTS cells present in both data and sploc\n")
      } else {
        cat("No overlapping GRTS cells found - using all data for this species\n")
      }
      
      if (nrow(x) == 0) {
        cat("No data remaining after filtering for species group\n")
        return(NULL)
      }
      
      print(n_distinct(x$GRTS_Cell_ID))
      
      cell.summary <- plyr::ddply(
        x, c("Year", "GRTS_Cell_ID", "Quadrant"), plyr::summarize,
        YearS              = Year[1] - 2016,
        total.detect       = sum(SpeciesSingleton),
        quad.nights        = length(SpeciesSingleton),
        Nightly_Mean_Temp  = mean(Nightly_Mean_Temp,  na.rm = TRUE),
        Nightly_Min_Temp   = mean(Nightly_Min_Temp,   na.rm = TRUE),
        Nightly_Max_Temp   = mean(Nightly_Max_Temp,   na.rm = TRUE),
        Nightly_Mean_RH    = mean(Nightly_Mean_RH,    na.rm = TRUE),
        Nightly_Mean_Windsp= mean(Nightly_Mean_Windsp,na.rm = TRUE),
        dd                 = mean(dd,                  na.rm = TRUE),
        jNight             = mean(jNight),
        moon_fraction_above_horizon   = mean(moon_fraction_above_horizon,na.rm = TRUE),
        moon               = mean(moon_illumination, na.rm = TRUE),
        Feature_Sampled    = unique(Feature_Sampled,   na.rm = TRUE)[1],
        Habitat_Type       = unique(Habitat_Type,      na.rm = TRUE)[1],
        Distance_to_Clutter__m_        = median(`Distance_to_Clutter_(m)` , na.rm = TRUE),
        Percent_Clutter    = mean(Percent_Clutter,     na.rm = TRUE),
        Water_Nearby       = median(Water_Nearby,      na.rm = TRUE)
      )
      cell.summary$countperquad <- cell.summary$total.detect / cell.summary$quad.nights
      
      if (!is.na(x$SpeciesGroup[1]) && x$SpeciesGroup[1] == "LACI") {
        cell.summary <- cell.summary[(cell.summary$countperquad < 75), ]
      }
      
      if (nrow(cell.summary) == 0) {
        cat("No data remaining after summary/filtering\n")
        return(NULL)
      }
      
      form <- paste("total.detect ~ YearS + log(quad.nights) + (1|GRTS_Cell_ID/Quadrant)",
                    paste(vars, collapse = " + "), collapse = "")
      cat("Formula:", form, "\n")
      
      fit <- tryCatch({
        glmmTMB::glmmTMB(formula = formula(form), family = nbinom2(link = "log"),
                         data = cell.summary, se = TRUE, verbose = FALSE)
      },
      error = function(cond) {
        cat("glmmTMB error:", conditionMessage(cond), "\n")
        return(NULL)
      },
      warning = function(cond) {
        cat("glmmTMB warning:", conditionMessage(cond), "\n")
        return(NULL)
      })
      
      if (is.null(fit)) cat("Model failed to converge for", x$SpeciesGroup[1], "\n")
      
      fit
    })
  }
)

summary(all.region.fits[["Okanagan Region"]][["MYLU"]])

all.region.slope.summary.singlet <- plyr::ldply(all.region.fits, function(region_fits) {
  plyr::ldply(region_fits, function(x) {
    if (is.null(x)) return(NULL)
    
    summary.table <- summary(x)
    slope <- as.data.frame(summary.table$coefficients$cond[2, , drop = FALSE])
    slope$intercept <- summary.table$coefficients$cond[1, , drop = FALSE]
    slope$n <- nobs(x)
    slope
  }, .id = "SpeciesGroup")
}, .id = "Region")

write.csv(all.region.slope.summary.singlet,"Data/Analyzed/ActivityTrends/region.slope.summary.singlet.csv") #this is what will be added to the report

