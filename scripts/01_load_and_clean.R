# =================================================================================================================
# Load detector files, assign species, keep paired site-nights, and build count
# tables (with zeros) at three time scales.
#
# First step is to take individual excel files from Cami for each BC region and merge them all together into a csv
# saved as combined_raw_bat_detections.csv
#
# OUTPUT OBJECTS (used by scripts 02-05)
#   dat           one row per bat pass, paired site-nights only
#   night_counts  passes per site x night x detector x species
#   counts_60     passes per site x 1-hour bin x detector x species
#   counts_10     passes per site x 10-minute bin x detector x species
#   spp_keep      species with enough data to analyse
#
# =================================================================================================================
#
# It builds two datasets, one per classifier, and saves each:
#   data/processed/bat_data_kaleidoscope.rds   <- main analysis
#   data/processed/bat_data_sonobat.rds        <- classifier comparison
# Each .rds file is a list containing:
#   dat           one row per bat pass, paired site-nights only
#   night_counts  passes per site x night x detector x species
#   counts_60     passes per site x 1-hour bin x detector x species
#   counts_10     passes per site x 10-minute bin x detector x species
#   spp_keep      species with enough data to analyse
#
# =============================================================================

library(tidyverse)
library(readxl)
library(lubridate)

select <- dplyr::select

# Run from the PROJECT folder (the one that contains data/raw), not from data/
setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity")
dir.create("data/processed", showWarnings = FALSE)


# -----------------------------------------------------------------------------
# LOAD
# -----------------------------------------------------------------------------

xlsx_files <- list.files("data/raw", pattern = "\\.xlsx$", full.names = TRUE)
raw <- map_dfr(xlsx_files, read_excel)
print(length(xlsx_files))   # number of files
print(nrow(raw))            # number of rows (recordings)


# -----------------------------------------------------------------------------
# DETECTORS
# -----------------------------------------------------------------------------
# detector_type: Chorus vs Legacy (Legacy = the WA unit at that site).
# Legacy is the reference level, so model ratios are Chorus / Legacy.

raw$detector_model <- case_when(
  raw$Model == "Chorus" ~ "Anabat_Chorus",
  raw$Model == "SM2"    ~ "SM2"
)   # any other model (e.g. Swift) becomes NA and is removed in the filter
raw$detector_type <- ifelse(raw$detector_model == "Anabat_Chorus", "Chorus", "Legacy")


# -----------------------------------------------------------------------------
# CLASSIFIER CHECK
# -----------------------------------------------------------------------------
#   Species Auto ID          = SonoBat
#   WA|Kaleidoscope|Auto ID  = Kaleidoscope
# Each detector model needs IDs from a classifier for that classifier's
# dataset to be built.

print(table(raw$detector_model, kaleidoscope_id = !is.na(raw$`WA|Kaleidoscope|Auto ID`)))
print(table(raw$detector_model, sonobat_id      = !is.na(raw$`Species Auto ID`)))

# Which classifier does the "Passed AutoID" column belong to? If it is only
# "True" where a SonoBat ID exists, it is a SonoBat field and should not be
# used to filter the Kaleidoscope dataset.
print(table(raw$`Passed AutoID`, sonobat_id = !is.na(raw$`Species Auto ID`), useNA = "ifany"))


# -----------------------------------------------------------------------------
# TIME
# -----------------------------------------------------------------------------

raw$Timestamp <- as.POSIXct(raw$Timestamp, tz = "America/Vancouver")
raw$night     <- as.Date(raw$Timestamp - 3600 * 12)          # noon-to-noon night
raw$hour_bin  <- floor_date(raw$Timestamp, unit = "1 hour")
raw$bin_10min <- floor_date(raw$Timestamp, unit = "10 minutes")
raw$site_id   <- paste(raw$CellName, raw$Quadrant, sep = "_")


# -----------------------------------------------------------------------------
# SPECIES CODES (both classifiers -> 4-letter codes)
# -----------------------------------------------------------------------------
# SonoBat codes ("Mylu") and Kaleidoscope codes ("MYOLUC"). Anything else
# (NoID, Noise, etc.) becomes NA and is dropped.

spp_crosswalk <- c(
  "Anpa" = "ANPA", "ANTPAL" = "ANPA",   # Pallid bat
  "Coto" = "COTO", "CORTOW" = "COTO",   # Townsend's big-eared bat
  "Epfu" = "EPFU", "EPTFUS" = "EPFU",   # Big brown bat
  "Euma" = "EUMA", "EUDMAC" = "EUMA",   # Spotted bat
  "Labo" = "LABO", "LASBOR" = "LABO",   # Eastern red bat
  "Laci" = "LACI", "LASCIN" = "LACI",   # Hoary bat
  "Lano" = "LANO", "LASNOC" = "LANO",   # Silver-haired bat
  "Myca" = "MYCA", "MYOCAL" = "MYCA",   # California myotis
  "Myci" = "MYCI", "MYOCIL" = "MYCI",   # Western small-footed myotis
  "Myev" = "MYEV", "MYOEVO" = "MYEV",   # Long-eared myotis
  "Mylu" = "MYLU", "MYOLUC" = "MYLU",   # Little brown myotis
  "Myse" = "MYSE", "MYOSEP" = "MYSE",   # Northern long-eared myotis
  "Myth" = "MYTH", "MYOTHY" = "MYTH",   # Fringed myotis
  "Myvo" = "MYVO", "MYOVOL" = "MYVO",   # Long-legged myotis
  "Myyu" = "MYYU", "MYOYUM" = "MYYU",   # Yuma myotis
  "Pahe" = "PAHE", "PARHES" = "PAHE"    # Canyon bat
)


# =============================================================================
# BUILD ONE DATASET PER CLASSIFIER
# =============================================================================
# Everything below runs twice: once with Kaleidoscope IDs, once with SonoBat.

for (classifier in c("kaleidoscope", "sonobat")) {
  
  print(classifier)
  
  id_column <- if (classifier == "kaleidoscope") "WA|Kaleidoscope|Auto ID" else "Species Auto ID"
  raw$spp <- unname(spp_crosswalk[raw[[id_column]]])
  
  # Skip this classifier if any detector model has no IDs from it
  ids_per_model <- tapply(!is.na(raw$spp), raw$detector_model, sum)
  print(ids_per_model)
  if (any(ids_per_model == 0)) {
    print(paste("Skipping", classifier, "- a detector model has no IDs from it"))
    next
  }
  
  # ---- Filter ---------------------------------------------------------------
  dat <- raw %>%
    filter(format(Timestamp, "%Y") == "2025",
           !grepl("transect", filepath, ignore.case = TRUE),
           !Quadrant %in% c("Not Found", "Transects"),
           !is.na(detector_model),
           `Passed AutoID` == "True",    # see CLASSIFIER CHECK above
           !is.na(spp))
  print(nrow(dat))   # bat passes after filtering
  
  # ---- Paired site-nights ---------------------------------------------------
  paired_nights <- dat %>%
    group_by(site_id, night) %>%
    summarise(n_types = n_distinct(detector_type), .groups = "drop") %>%
    filter(n_types == 2) %>%
    select(site_id, night)
  dat <- inner_join(dat, paired_nights, by = c("site_id", "night"))
  
  # Which legacy unit is at each site (SM2)
  site_legacy <- dat %>%
    filter(detector_type == "Legacy") %>%
    distinct(site_id, legacy_model = detector_model)
  print(count(site_legacy, site_id) %>% filter(n > 1))   # should be empty
  dat <- left_join(dat, site_legacy, by = "site_id")
  
  # ---- Species with enough data ---------------------------------------------
  spp_keep <- dat %>%
    group_by(spp) %>%
    summarise(passes = n(), site_nights = n_distinct(site_id, night)) %>%
    filter(passes >= 50) %>%
    pull(spp)
  print(spp_keep)
  
  # ---- Count tables with zeros ----------------------------------------------
  # Same three steps at each scale:
  #   1. every time bin with any bat activity at the site (detector NOT used)
  #   2. give that bin both detectors and every species
  #   3. fill in the observed counts; everything else is 0
  
  obs <- dat %>%
    filter(spp %in% spp_keep) %>%
    count(site_id, night, detector_type, spp, name = "n_calls")
  night_counts <- dat %>%
    distinct(site_id, legacy_model, night) %>%
    crossing(detector_type = c("Legacy", "Chorus"), spp = spp_keep) %>%
    left_join(obs, by = c("site_id", "night", "detector_type", "spp")) %>%
    mutate(n_calls = replace_na(n_calls, 0),
           detector_type = factor(detector_type, levels = c("Legacy", "Chorus")))
  
  obs <- dat %>%
    filter(spp %in% spp_keep) %>%
    count(site_id, night, hour_bin, detector_type, spp, name = "n_calls")
  counts_60 <- dat %>%
    distinct(site_id, legacy_model, night, hour_bin) %>%
    crossing(detector_type = c("Legacy", "Chorus"), spp = spp_keep) %>%
    left_join(obs, by = c("site_id", "night", "hour_bin", "detector_type", "spp")) %>%
    mutate(n_calls = replace_na(n_calls, 0),
           detector_type = factor(detector_type, levels = c("Legacy", "Chorus")))
  
  obs <- dat %>%
    filter(spp %in% spp_keep) %>%
    count(site_id, night, bin_10min, detector_type, spp, name = "n_calls")
  counts_10 <- dat %>%
    distinct(site_id, legacy_model, night, bin_10min) %>%
    crossing(detector_type = c("Legacy", "Chorus"), spp = spp_keep) %>%
    left_join(obs, by = c("site_id", "night", "bin_10min", "detector_type", "spp")) %>%
    mutate(n_calls = replace_na(n_calls, 0),
           detector_type = factor(detector_type, levels = c("Legacy", "Chorus")))
  
  # ---- Checks and numbers for the methods -----------------------------------
  # Grid check: Legacy and Chorus must have the SAME number of rows
  print(table(night_counts$detector_type))
  print(table(counts_60$detector_type))
  print(table(counts_10$detector_type))
  
  # Passes per detector model
  print(count(dat, detector_model))
  
  # Sites and site-nights per legacy model
  print(night_counts %>%
          distinct(site_id, night, legacy_model) %>%
          group_by(legacy_model) %>%
          summarise(sites = n_distinct(site_id), site_nights = n()))
  
  # Passes per species per detector model
  print(dat %>%
          filter(spp %in% spp_keep) %>%
          count(spp, detector_model) %>%
          pivot_wider(names_from = detector_model, values_from = n, values_fill = 0))
  
  # Grids (CellName) and sites (grid + quadrant) used, with effort per site
  sites_used <- dat %>%
    group_by(CellName, Quadrant, site_id, legacy_model) %>%
    summarise(nights        = n_distinct(night),
              first_night   = min(night),
              last_night    = max(night),
              passes_chorus = sum(detector_type == "Chorus"),
              passes_legacy = sum(detector_type == "Legacy"),
              .groups = "drop") %>%
    arrange(legacy_model, CellName, Quadrant)
  print(sites_used, n = Inf)
  
  # Number of grids and sites per legacy model
  print(sites_used %>%
          group_by(legacy_model) %>%
          summarise(grids = n_distinct(CellName), sites = n(), nights = sum(nights)))
  
  write.csv(sites_used, paste0("data/processed/sites_used_", classifier, ".csv"), row.names = FALSE)
  
  # ---- Save -----------------------------------------------------------------
  write.csv(dat, paste0("data/processed/bat_detections_", classifier, ".csv"), row.names = FALSE)
  saveRDS(list(dat = dat, night_counts = night_counts, counts_60 = counts_60,
               counts_10 = counts_10, spp_keep = spp_keep),
          paste0("data/processed/bat_data_", classifier, ".rds"))
}

