
library(tidyverse)   # data wrangling and plotting
library(readxl)      # reading Excel files
library(MASS)        # glm.nb() — negative binomial models
library(lme4)        # glmer.nb() — negative binomial mixed models

# MASS overwrites dplyr's select() function — this line fixes that
select <- dplyr::select

setwd()

# =============================================================================
# LOAD AND CLEAN THE DATA
# =============================================================================

# --- Read all xlsx files from your working directory -------------------------
# Make sure all 9 xlsx files are in the same folder as this script,
# and that folder is set as your working directory (Session > Set Working Dir)

xlsx_files <- list.files(pattern = "\\.xlsx$", full.names = TRUE)
cat("Found", length(xlsx_files), "files\n")

raw <- map_dfr(xlsx_files, read_excel)
cat("Total rows:", nrow(raw), "\n")

#write.csv(raw, "combined_raw_bat_detections.csv")

# --- Label detector types ----------------------------------------------------
# The two detector types we want to compare are:
#   Anabat Chorus (Titley Scientific)  -> "Anabat_Chorus"
#   SM2 or Swift  (Wildlife Acoustics) -> "Other_FS"
# Everything else (Walkabout, Echo Meter Touch) gets NA and is dropped later.

raw$detector_type <- NA_character_
raw$detector_type[raw$Model == "Chorus"]               <- "Anabat_Chorus"
raw$detector_type[raw$Model %in% c("SM2", "Swift")]    <- "Other_FS"

# Make it a factor with Other_FS as the reference level.
# This means model results will tell us how Chorus compares TO Other_FS.
raw$detector_type <- factor(raw$detector_type, levels = c("Other_FS", "Anabat_Chorus"))


# --- Standardise species names -----------------------------------------------
# SonoBat (one classifier) uses informal codes like "Mylu" or "Epfu".
# Kaleidoscope (another classifier) uses codes like "MYOLUC" or "EPTFUS".
# These are the same species — we map everything to standard 4-letter codes.
# Anything not in this list (Noise, NoID) becomes NA and gets dropped later.

spp_crosswalk <- c(
  "Anpa" = "ANPA",  "ANTPAL" = "ANPA",   # Pallid bat
  "Coto" = "COTO",  "CORTOW" = "COTO",   # Townsend's big-eared bat
  "Epfu" = "EPFU",  "EPTFUS" = "EPFU",   # Big brown bat
  "Euma" = "EUMA",  "EUDMAC" = "EUMA",   # Spotted bat
  "Labo" = "LABO",  "LASBOR" = "LABO",   # Silver-haired bat
  "Laci" = "LACI",  "LASCIN" = "LACI",   # Hoary bat
  "Lano" = "LANO",  "LASNOC" = "LANO",   # Silver-haired bat
  "Myca" = "MYCA",  "MYOCAL" = "MYCA",   # California myotis
  "Myci" = "MYCI",  "MYOCIL" = "MYCI",   # Western small-footed myotis
  "Myev" = "MYEV",  "MYOEVO" = "MYEV",   # Long-eared myotis
  "Mylu" = "MYLU",  "MYOLUC" = "MYLU",   # Little brown myotis
  "Myse" = "MYSE",  "MYOSEP" = "MYSE",   # Northern long-eared myotis
  "Myth" = "MYTH",  "MYOTHY" = "MYTH",   # Fringed myotis
  "Myvo" = "MYVO",  "MYOVOL" = "MYVO",   # Long-legged myotis
  "Myyu" = "MYYU",  "MYOYUM" = "MYYU",   # Yuma myotis
  "Pahe" = "PAHE",  "PARHES" = "PAHE"    # Canyon bat
)

# Use SonoBat species ID where available, otherwise use Kaleidoscope ID
raw$spp_raw <- coalesce(raw$`Species Auto ID`, raw$`WA|Kaleidoscope|Auto ID`)

# Look up each value in the crosswalk — anything not in the list becomes NA
raw$spp <- spp_crosswalk[raw$spp_raw]


# --- Create timestamp and biological night -----------------------------------
# "Biological night" means we assign all recordings before noon to the
# previous night. So a recording at 02:00 on June 17 happened during the
# night that started on June 16 — we want to keep those together.

raw$Timestamp <- as.POSIXct(raw$Timestamp, tz = "America/Vancouver")
raw$night     <- as.Date(raw$Timestamp - 3600 * 12)   # subtract 12 hours, take the date
raw$site_id   <- paste(raw$CellName, raw$Quadrant, sep = "_")


# --- Filter to what we need --------------------------------------------------

dat <- raw %>%
  filter(
    format(Timestamp, "%Y") == "2025",                          # 2025 only
    !grepl("Transect|transect", filepath, ignore.case = TRUE),  # no driving transects
    !Quadrant %in% c("Not Found", "Transects"),                 # no unmatched sites
    !is.na(detector_type),                                      # known detector only
    `Passed AutoID` == "True",                                  # passed quality filter
    !is.na(spp)                                                 # identified species only
  )

cat("Bat passes after filtering:", nrow(dat), "\n")


# --- Keep only paired sites (both detectors present) -------------------------
# Some sites only had one detector deployed. We drop those because we can't
# compare detectors at a site where only one was running.

paired_sites <- dat %>%
  group_by(site_id) %>%
  summarise(
    has_chorus = any(detector_type == "Anabat_Chorus"),
    has_other  = any(detector_type == "Other_FS")
  ) %>%
  filter(has_chorus & has_other) %>%
  pull(site_id)

dat <- filter(dat, site_id %in% paired_sites)

# Within paired sites, keep only nights where BOTH detectors were active.
# (Some nights are edge nights — one detector was deployed or collected mid-night.)

paired_nights <- dat %>%
  group_by(site_id, night) %>%
  summarise(n_det_types = n_distinct(detector_type), .groups = "drop") %>%
  filter(n_det_types == 2) %>%
  select(site_id, night)

dat <- inner_join(dat, paired_nights, by = c("site_id", "night"))

cat("Paired sites:", n_distinct(dat$site_id), "\n")
cat("Paired site-nights:", nrow(paired_nights), "\n")

write.csv(dat, "bat_detections.csv")
# --- Build nightly count data with zeros -------------------------------------
# We need to count how many times each species was detected each night by each
# detector. Crucially, we also need the ZEROS — nights when a detector was
# running but didn't detect a species. Without zeros, our model would think
# every night had some detections, which would overestimate abundance.



# Step A: all active site-night-detector combinations
all_combos <- dat %>%
  distinct(site_id, CellName, Quadrant, detector_type, night)

# Step B: decide which species to model.
# We need enough data to fit a model reliably
# (rare species with only a handful of records won't give reliable results)

spp_keep <- dat %>% 
  group_by(spp) %>%
  summarise(
    total = n(),
    site_nights = n_distinct(paste(site_id, night))
  ) %>%
  filter(total >= 50, site_nights >= 10) %>%
  arrange(desc(total)) %>%
  pull(spp)

print(spp_keep)

# count observed detections per species per site-night-detector
obs_counts <- dat %>%
  group_by(site_id, CellName, Quadrant, detector_type, night, spp) %>%
  summarise(n_calls = n(), .groups = "drop")

# make the complete grid and fill unobserved cells with zero
night_counts <- crossing(all_combos, spp = spp_keep) %>%
  left_join(obs_counts,
            by = c("site_id","CellName","Quadrant","detector_type","night","spp")) %>%
  mutate(n_calls = replace_na(n_calls, 0))

cat("Total rows in model dataset (including zeros):", nrow(night_counts), "\n")

# Quick check — what does this look like?
# Each row is one night, one detector, one species, at one site.
head(night_counts)
