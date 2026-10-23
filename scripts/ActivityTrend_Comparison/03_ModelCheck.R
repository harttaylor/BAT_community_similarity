library(DHARMa)
library(performance)
library(insight)

# =============================================================================
# Runs the same diagnostics on two model types:
#   Stationary: total.detect ~ YearS + log(quad.nights) + (1|GRTS_Cell_ID/Quadrant) + covs
#   Transect:   total.detect ~ YearS + log(transect.length) + (1|GRTS_Cell_ID) + covs
# The grouping structure is detected per model: if Quadrant is in the model it
# is treated as nested within cell; otherwise cell is the only grouping level.
# =============================================================================

# ---- SETTINGS (shared) -------------------------------------------------------
n_sim    <- 1000          # DHARMa simulations (250 is faster for a first pass)
seed     <- 123

response <- "total.detect"
year_var <- "YearS"
cell_var <- "GRTS_Cell_ID"
quad_var <- "Quadrant"
lon_var  <- "Longitude"
lat_var  <- "Latitude"

# Coordinates for the spatial test (NULL = skip the spatial test).
# Stationary: one row per GRTS cell x quadrant, with columns
#   GRTS_Cell_ID, Quadrant, Longitude, Latitude. Example:
#   raw <- readxl::read_excel("AllYearActivitybyNight_2016_2025_BC_withALLmetadata.xlsx",
#                             sheet = "OccupancyKSP_NEAR")
#   coords_stationary <- aggregate(cbind(Longitude, Latitude) ~ GRTS_Cell_ID + Quadrant,
#                                  data = raw, FUN = mean)
# Transect: one row per GRTS cell, with columns GRTS_Cell_ID, Longitude, Latitude
#   (e.g. transect midpoint). If it also has a Quadrant column, it is averaged
#   to one point per cell automatically.
coords_stationary <- NULL
coords_transect   <- NULL
# -----------------------------------------------------------------------------

# ---- Helpers -----------------------------------------------------------------

is_model <- function(x) inherits(x, c("glmmTMB", "merMod"))

# Pull every model object out of a (possibly nested) list
flatten_models <- function(x, prefix = NULL) {
  if (is_model(x)) return(setNames(list(x), if (is.null(prefix)) "model" else prefix))
  if (!is.list(x)) return(list())
  nms <- names(x)
  if (is.null(nms)) nms <- as.character(seq_along(x))
  nms[nms == ""] <- as.character(which(nms == ""))
  out <- list()
  for (i in seq_along(x)) {
    nm  <- if (is.null(prefix)) nms[i] else paste(prefix, nms[i], sep = "_")
    out <- c(out, flatten_models(x[[i]], nm))
  }
  out
}

safe <- function(expr, label = "") {
  tryCatch(expr, error = function(e) {
    message("    ! ", label, " failed: ", conditionMessage(e))
    NULL
  })
}

ht <- function(x, prefix) {
  if (is.null(x)) return(setNames(list(NA_real_, NA_real_), paste0(prefix, c("_stat", "_p"))))
  setNames(list(unname(x$statistic[1]), unname(x$p.value)), paste0(prefix, c("_stat", "_p")))
}

bind <- function(lst) {
  if (!length(lst)) return(NULL)
  cols <- unique(unlist(lapply(lst, names)))
  do.call(rbind, lapply(lst, function(d) { d[setdiff(cols, names(d))] <- NA; d[cols] }))
}

# ---- Main function -----------------------------------------------------------

run_diagnostics <- function(model_list, out_dir, coords_df = NULL) {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  models <- flatten_models(model_list)
  message("Found ", length(models), " model(s)")
  
  summary_rows <- list(); cat_rows <- list(); vif_rows <- list()
  
  for (nm in names(models)) {
    m <- models[[nm]]
    message("\n== ", nm, " ==")
    set.seed(seed)
    
    mf <- model.frame(m)   # rows used in the fit, incl. log(quad.nights) / log(transect.length)
    
    # ---- Grouping structure: detected per model ---------------------------------
    has_cell <- cell_var %in% names(mf)
    has_quad <- has_cell && quad_var %in% names(mf)     # TRUE = stationary (nested)
    
    cell_f <- if (has_cell) factor(mf[[cell_var]]) else NULL
    cellquad_f <- if (has_quad) factor(paste(mf[[cell_var]], mf[[quad_var]], sep = "_")) else NULL
    quads_per_cell <- if (has_quad) {
      tapply(as.character(mf[[quad_var]]), cell_f, function(q) length(unique(q)))
    } else NULL
    
    # Spatial unit: cell x quadrant for stationary, cell for transect
    site_f <- if (has_quad) cellquad_f else cell_f
    
    message("    structure: ", if (has_quad) "cell / quadrant (nested)" else "cell only")
    
    pdf(file.path(out_dir, paste0("diagnostics_", make.names(nm), ".pdf")),
        width = 11, height = 7)
    
    # ---- 1. Simulated residuals, main plots, uniformity, outliers -------------
    res <- safe(simulateResiduals(m, n = n_sim, plot = FALSE), "simulateResiduals")
    if (is.null(res)) {
      dev.off()
      summary_rows[[nm]] <- data.frame(model = nm, error = "simulateResiduals failed")
      next
    }
    
    plot(res, title = paste("DHARMa residuals:", nm))
    uni <- safe(testUniformity(res, plot = FALSE), "uniformity")
    out <- safe(testOutliers(res, plot = FALSE,
                             type = if (nrow(mf) > 500) "binomial" else "bootstrap"),
                "outliers")
    
    # Residuals vs every fixed-effect predictor (YearS, log effort term,
    # species-specific covariates). Excludes the response and grouping factors.
    fixed_cols <- setdiff(names(mf), c(response, cell_var, quad_var,
                                       "(weights)", "(offset)"))
    fixed_cols <- fixed_cols[!grepl("^offset\\(", fixed_cols)]
    for (p in fixed_cols) {
      x <- mf[[p]]
      if (is.character(x)) x <- factor(x)
      if (is.matrix(x)) next   # e.g. poly() terms; their source variable is plotted instead
      safe(plotResiduals(res, form = x, xlab = p, main = paste("Residuals vs", p)),
           paste("plotResiduals", p))
    }
    
    # ---- 2. Dispersion ----------------------------------------------------------
    disp <- safe(testDispersion(res, plot = TRUE), "dispersion")
    
    # ---- 3. Zero inflation ------------------------------------------------------
    zi       <- safe(testZeroInflation(res, plot = TRUE), "zero-inflation")
    obs_zero <- mean(mf[[response]] == 0, na.rm = TRUE)
    
    # ---- 4a. Temporal autocorrelation (aggregated by year) --------------------
    temp <- NULL; n_years <- NA
    if (year_var %in% names(mf)) {
      yr      <- mf[[year_var]]
      n_years <- length(unique(yr))
      res_t   <- safe(recalculateResiduals(res, group = yr), "group by year")
      if (!is.null(res_t)) {
        temp <- safe(testTemporalAutocorrelation(res_t, time = sort(unique(yr)), plot = TRUE),
                     "temporal autocorrelation")
      }
    } else {
      message("    ! '", year_var, "' not in model frame: temporal test skipped")
    }
    
    # ---- 4b. Spatial autocorrelation (by cell x quadrant, or by cell) ---------
    spat <- NULL
    if (!is.null(coords_df) && !is.null(site_f)) {
      cd <- as.data.frame(coords_df)
      if (has_quad) {
        ck <- paste(cd[[cell_var]], cd[[quad_var]], sep = "_")
      } else {
        # one point per cell (averages quadrant coordinates if present)
        cd <- aggregate(cd[c(lon_var, lat_var)], by = list(cell = cd[[cell_var]]), FUN = mean)
        ck <- as.character(cd$cell)
      }
      res_s <- safe(recalculateResiduals(res, group = site_f), "group by site")
      if (!is.null(res_s)) {
        xy <- cd[match(levels(site_f), ck), c(lon_var, lat_var)]
        if (anyNA(xy)) {
          message("    ! ", sum(!complete.cases(xy)),
                  " site(s) without coordinates: spatial test skipped")
        } else {
          dup <- duplicated(xy)
          if (any(dup)) xy[[lon_var]][dup] <- xy[[lon_var]][dup] + runif(sum(dup), -1e-5, 1e-5)
          spat <- safe(testSpatialAutocorrelation(res_s, x = xy[[lon_var]], y = xy[[lat_var]],
                                                  plot = TRUE),
                       "spatial autocorrelation")
        }
      }
    } else if (is.null(coords_df)) {
      message("    - no coordinates supplied: spatial test skipped")
    }
    
    # ---- 5. Random effects --------------------------------------------------------
    grp_list <- list()
    if (has_cell) grp_list$cell <- cell_f
    if (has_quad) grp_list$cell_quadrant <- cellquad_f
    for (g in names(grp_list)) {
      grp <- grp_list[[g]]
      ct  <- safe(testCategorical(res, catPred = grp, plot = TRUE), paste("categorical", g))
      if (!is.null(ct)) {
        cat_rows[[paste(nm, g)]] <- data.frame(
          model = nm, group = g, n_levels = nlevels(grp),
          n_levels_nonuniform_p05 = sum(ct$uniformity$p.value < 0.05, na.rm = TRUE),
          levene_F = unname(ct$homogeneity$`F value`[1]),
          levene_p = unname(ct$homogeneity$`Pr(>F)`[1])
        )
      }
    }
    
    # QQ plots of the random intercepts ("GRTS_Cell_ID" and, for stationary
    # models, "Quadrant:GRTS_Cell_ID")
    re <- safe({ r <- lme4::ranef(m); if (inherits(m, "glmmTMB")) r$cond else r }, "ranef")
    re_sd <- c(cell = NA, cell_quadrant = NA)
    if (!is.null(re)) {
      for (g in names(re)) {
        v     <- re[[g]][["(Intercept)"]]
        label <- if (grepl(":", g)) "cell_quadrant" else if (g == cell_var) "cell" else g
        re_sd[label] <- sd(v)
        qqnorm(v, main = paste0("Random intercept QQ: ", g, " (", length(v), " levels)"))
        qqline(v, col = "red")
      }
    }
    
    # ---- 6. R2, singularity, collinearity ---------------------------------------
    # r2_nakagawa() returns a plain NA (with a warning) instead of a list when it
    # can't compute R2, e.g. when a random-effect variance is ~0
    r2  <- safe(r2_nakagawa(m), "r2_nakagawa")
    if (!is.list(r2)) {
      if (!is.null(r2)) message("    ! r2_nakagawa returned NA (often a near-zero random-effect variance)")
      r2 <- NULL
    }
    sng <- safe(check_singularity(m), "singularity")
    vif <- safe(check_collinearity(m), "collinearity")
    # check_collinearity() returns an empty table (0 rows) when there are fewer
    # than two fixed-effect terms to compare
    if (!is.data.frame(vif) || !"VIF" %in% names(vif) || nrow(vif) == 0) {
      if (!is.null(vif)) message("    - collinearity: fewer than 2 fixed-effect terms, VIF not computed")
      vif <- NULL
    }
    if (!is.null(vif)) vif_rows[[nm]] <- data.frame(model = nm, as.data.frame(vif))
    
    dev.off()
    
    # ---- Collect results -----------------------------------------------------------
    summary_rows[[nm]] <- data.frame(
      model               = nm,
      structure           = if (has_quad) "cell/quadrant" else "cell",
      n_obs               = nrow(mf),
      n_years             = n_years,
      n_cells             = if (has_cell) nlevels(cell_f) else NA,
      n_cell_quadrants    = if (has_quad) nlevels(cellquad_f) else NA,
      cells_with_1_quad   = if (has_quad) sum(quads_per_cell == 1) else NA,
      ht(uni,  "uniformity_KS"),
      outliers_n          = sum(outliers(res)),
      outliers_p          = if (!is.null(out)) out$p.value else NA,
      ht(disp, "dispersion"),
      ht(zi,   "zeroinfl_ratio"),
      observed_zero_prop  = obs_zero,
      ht(temp, "temporal_DW"),
      ht(spat, "spatial_MoranI"),
      RE_sd_cell          = unname(re_sd["cell"]),
      RE_sd_cell_quadrant = unname(re_sd["cell_quadrant"]),
      R2_marginal         = if (!is.null(r2)) unname(r2$R2_marginal)    else NA,
      R2_conditional      = if (!is.null(r2)) unname(r2$R2_conditional) else NA,
      singular            = if (!is.null(sng)) sng else NA,
      max_VIF             = if (!is.null(vif)) max(vif$VIF, na.rm = TRUE) else NA,
      error               = NA_character_,
      check.names = FALSE
    )
    message("    done")
  }
  
  # ---- Write tables ----------------------------------------------------------------
  diag_summary <- bind(summary_rows); rownames(diag_summary) <- NULL
  cat_summary  <- bind(cat_rows)
  vif_summary  <- bind(vif_rows)
  
  write.csv(diag_summary, file.path(out_dir, "diagnostics_summary.csv"), row.names = FALSE)
  if (!is.null(cat_summary)) write.csv(cat_summary, file.path(out_dir, "categorical_tests.csv"), row.names = FALSE)
  if (!is.null(vif_summary)) write.csv(vif_summary, file.path(out_dir, "collinearity.csv"),     row.names = FALSE)
  
  print(diag_summary, digits = 3)
  message("\nAll outputs written to: ", normalizePath(out_dir))
  invisible(diag_summary)
}

# ---- RUN ---------------------------------------------------------------------
# Replace the model list names with your own objects.

diag_stationary <- run_diagnostics(
  all.region.fits,
  out_dir   = "Data/Analyzed/ActivityTrends/model_diagnostics/ProvRegional",
  coords_df = coords_stationary
)

diag_transect <- run_diagnostics(
  transect_fits,                                   # <- your transect model list
  out_dir   = "Data/Analyzed/ActivityTrends/model_diagnostics/Transect",
  coords_df = coords_transect
)