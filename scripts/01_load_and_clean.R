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
library(tidyverse)
library(readxl)
library(lubridate)

select <- dplyr::select
setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity/data/raw")


# -----------------------------------------------------------------------------
# LOAD
# -----------------------------------------------------------------------------

xlsx_files <- list.files(pattern = "\\.xlsx$", full.names = TRUE)
raw <- map_dfr(xlsx_files, read_excel)
print(length(xlsx_files))   # number of files
print(nrow(raw))            # number of rows (recordings)


# -----------------------------------------------------------------------------
# DETECTORS
# -----------------------------------------------------------------------------
# detector_model: the actual unit (Anabat_Chorus, SM2, Swift)
# detector_type:  Chorus vs Legacy (Legacy = whichever WA unit was at the site)
# Legacy is the reference level, so model ratios are Chorus / Legacy.

raw$detector_model <- case_when(
  raw$Model == "Chorus" ~ "Anabat_Chorus",
  raw$Model == "SM2"    ~ "SM2",
  raw$Model == "Swift"  ~ "Swift"
)
raw$detector_type <- ifelse(raw$detector_model == "Anabat_Chorus", "Chorus", "Legacy")


# -----------------------------------------------------------------------------
# CLASSIFIER CHECK
# -----------------------------------------------------------------------------
# Which classifier columns are filled for each detector model?
#   Species Auto ID          = SonoBat
#   WA|Kaleidoscope|Auto ID  = Kaleidoscope
# Every detector model needs Kaleidoscope IDs. If a model shows only FALSE in
# the Kaleidoscope table, its files must be run through Kaleidoscope first.

print(table(raw$detector_model, kaleidoscope_id = !is.na(raw$`WA|Kaleidoscope|Auto ID`)))
print(table(raw$detector_model, sonobat_id      = !is.na(raw$`Species Auto ID`)))

# "Passed AutoID" filter (used below): check which classifier it belongs to.
# If it is only ever "True" where a SonoBat ID exists, it is a SonoBat field
# and should be removed from the filter so that only Kaleidoscope is used.
print(table(raw$`Passed AutoID`, sonobat_id = !is.na(raw$`Species Auto ID`), useNA = "ifany"))

no_kal <- raw %>%
  filter(!is.na(detector_model)) %>%
  group_by(detector_model) %>%
  summarise(n_kaleidoscope = sum(!is.na(`WA|Kaleidoscope|Auto ID`))) %>%
  filter(n_kaleidoscope == 0)
if (nrow(no_kal) > 0) {
  stop(paste("No Kaleidoscope IDs for:", paste(no_kal$detector_model, collapse = ", ")))
}


# -----------------------------------------------------------------------------
# SPECIES (Kaleidoscope codes -> 4-letter codes)
# -----------------------------------------------------------------------------
# Anything not in this list (NoID, Noise, etc.) becomes NA and is dropped.

spp_crosswalk <- c(
  "ANTPAL" = "ANPA",   # Pallid bat
  "CORTOW" = "COTO",   # Townsend's big-eared bat
  "EPTFUS" = "EPFU",   # Big brown bat
  "EUDMAC" = "EUMA",   # Spotted bat
  "LASBOR" = "LABO",   # Eastern red bat
  "LASCIN" = "LACI",   # Hoary bat
  "LASNOC" = "LANO",   # Silver-haired bat
  "MYOCAL" = "MYCA",   # California myotis
  "MYOCIL" = "MYCI",   # Western small-footed myotis
  "MYOEVO" = "MYEV",   # Long-eared myotis
  "MYOLUC" = "MYLU",   # Little brown myotis
  "MYOSEP" = "MYSE",   # Northern long-eared myotis
  "MYOTHY" = "MYTH",   # Fringed myotis
  "MYOVOL" = "MYVO",   # Long-legged myotis
  "MYOYUM" = "MYYU",   # Yuma myotis
  "PARHES" = "PAHE"    # Canyon bat
)
raw$spp <- unname(spp_crosswalk[raw$`WA|Kaleidoscope|Auto ID`])

# =============================================================================
# TIME
# =============================================================================
# --- Create timestamp and biological night -----------------------------------
# "Biological night" means we assign all recordings before noon to the
# previous night. So a recording at 02:00 on June 17 happened during the
# night that started on June 16 — we want to keep those together.

raw$Timestamp <- as.POSIXct(raw$Timestamp, tz = "America/Vancouver")
raw$night     <- as.Date(raw$Timestamp - 3600 * 12)          # noon-to-noon night
raw$hour_bin  <- floor_date(raw$Timestamp, unit = "1 hour")
raw$bin_10min <- floor_date(raw$Timestamp, unit = "10 minutes")
raw$site_id   <- paste(raw$CellName, raw$Quadrant, sep = "_")



# =============================================================================
# FILTER
# =============================================================================

dat <- raw %>%
  filter(format(Timestamp, "%Y") == "2025",
         !grepl("transect", filepath, ignore.case = TRUE),
         !Quadrant %in% c("Not Found", "Transects"),
         !is.na(detector_model),
         `Passed AutoID` == "True",    # see CLASSIFIER CHECK above
         !is.na(spp))

print(nrow(dat))   # bat passes after filtering


# =============================================================================
# KEEP ONLY PAIRED SITES AND PAIRED NIGHTS
# =============================================================================

paired_nights <- dat %>%
  group_by(site_id, night) %>%
  summarise(n_types = n_distinct(detector_type), .groups = "drop") %>%
  filter(n_types == 2) %>%
  select(site_id, night)

dat <- inner_join(dat, paired_nights, by = c("site_id", "night"))

# Which legacy unit is at each site (SM2 or Swift)
site_legacy <- dat %>%
  filter(detector_type == "Legacy") %>%
  distinct(site_id, legacy_model = detector_model)
print(count(site_legacy, site_id) %>% filter(n > 1))   # should be empty

dat <- left_join(dat, site_legacy, by = "site_id")

write.csv(dat, "bat_detections.csv", row.names = FALSE)


# =============================================================================
# BUILD NIGHTLY COUNTS WITH ZEROS
# =============================================================================

# Species with enough data 
spp_keep <- dat %>%
  group_by(spp) %>%
  summarise(passes = n(), site_nights = n_distinct(site_id, night)) %>%
  filter(passes >= 50, site_nights >= 10) %>%
  pull(spp)
print(spp_keep)

# Count tables with zeros
# Same three steps at each scale:
#   1. every time bin with any bat activity at the site (detector NOT used here)
#   2. give that bin both detectors and every species
#   3. fill in the observed counts; everything else is 0

# Nightly
obs <- dat %>%
  filter(spp %in% spp_keep) %>%
  count(site_id, night, detector_type, spp, name = "n_calls")
night_counts <- dat %>%
  distinct(site_id, legacy_model, night) %>%
  crossing(detector_type = c("Legacy", "Chorus"), spp = spp_keep) %>%
  left_join(obs, by = c("site_id", "night", "detector_type", "spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

# 1-hour bins
obs <- dat %>%
  filter(spp %in% spp_keep) %>%
  count(site_id, night, hour_bin, detector_type, spp, name = "n_calls")
counts_60 <- dat %>%
  distinct(site_id, legacy_model, night, hour_bin) %>%
  crossing(detector_type = c("Legacy", "Chorus"), spp = spp_keep) %>%
  left_join(obs, by = c("site_id", "night", "hour_bin", "detector_type", "spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

# 10-minute bins
obs <- dat %>%
  filter(spp %in% spp_keep) %>%
  count(site_id, night, bin_10min, detector_type, spp, name = "n_calls")
counts_10 <- dat %>%
  distinct(site_id, legacy_model, night, bin_10min) %>%
  crossing(detector_type = c("Legacy", "Chorus"), spp = spp_keep) %>%
  left_join(obs, by = c("site_id", "night", "bin_10min", "detector_type", "spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

# Factor levels: Legacy is the reference
night_counts$detector_type <- factor(night_counts$detector_type, levels = c("Legacy", "Chorus"))
counts_60$detector_type    <- factor(counts_60$detector_type,    levels = c("Legacy", "Chorus"))
counts_10$detector_type    <- factor(counts_10$detector_type,    levels = c("Legacy", "Chorus"))


# =============================================================================
# CHECKS AND NUMBERS FOR THE METHODS
# =============================================================================

# Grid check: Legacy and Chorus must have the SAME number of rows at each scale
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

# Passes per species per detector (paired data)
print(dat %>% filter(spp %in% spp_keep) %>% count(spp, detector_model) %>%
        pivot_wider(names_from = detector_model, values_from = n, values_fill = 0))



head(night_counts)
write.csv(night_counts, "night_counts.csv")
