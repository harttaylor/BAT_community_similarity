
library(tidyverse)   # data wrangling and plotting
library(readxl)      # reading Excel files

# MASS overwrites dplyr's select() function — this line fixes that
select <- dplyr::select

setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity/data/raw")

# =============================================================================
# LOAD AND CLEAN THE DATA
# =============================================================================

# --- Read all xlsx files from your working directory -------------------------

xlsx_files <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

raw <- map_dfr(xlsx_files, read_excel)
cat("Total rows:", nrow(raw), "\n")

#write.csv(raw, "combined_raw_bat_detections.csv")

# --- Label detector models ----------------------------------------------------
# The two detector types we want to compare are:
#   Anabat Chorus (Titley Scientific)  -> "Anabat_Chorus"
#   SM2 or Swift  (Wildlife Acoustics) -> "Other_FS"
# Everything else (Walkabout, Echo Meter Touch) gets NA and is dropped later.
#raw <- read.csv("data/raw/combined_raw_bat_detections.csv")
raw$detector_model <- NA_character_
raw$detector_model[raw$Model == "Chorus"]      <- "Anabat_Chorus"
raw$detector_model[raw$Model == "SM2"]         <- "SM2"
raw$detector_model[raw$Model == "Swift"]       <- "Swift"

# Make it a factor with Other_FS as the reference level.
# This means model results will tell us how Chorus compares to the legacy units (SM2 and swift).
raw$detector_model <- factor(raw$detector_model, levels = c("SM2", "Swift", "Anabat_Chorus"))

# But now add a two-level grouping, using original name Other_FS 
# Reference level is still the legacy units, so a positive coefficient still 
# means "Chorus records more"

raw$detector_type <- case_when(
  raw$detector_model == "Anabat_Chorus" ~ "Anabat_Chorus",
  raw$detector_model %in% c("SM2", "Swift") ~ "Other_FS",
  TRUE ~ NA_character_
)
raw$detector_type <- factor(raw$detector_type, 
                            levels = c("Other_FS", "Anabat_Chorus"))


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
  "Labo" = "LABO",  "LASBOR" = "LABO",   # Eastern red bat  <- see note below
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


# =============================================================================
# FILTER
# =============================================================================

dat <- raw %>%
  filter(
    format(Timestamp, "%Y") == "2025",
    !grepl("Transect|transect", filepath, ignore.case = TRUE),
    !Quadrant %in% c("Not Found", "Transects"),
    !is.na(detector_model),
    `Passed AutoID` == "True",
    !is.na(spp)
  )

cat("Bat passes after filtering:", nrow(dat), "\n")


# =============================================================================
# KEEP ONLY PAIRED SITES AND PAIRED NIGHTS
# =============================================================================

paired_sites <- dat %>%
  group_by(site_id) %>%
  summarise(
    has_chorus = any(detector_type == "Anabat_Chorus"),
    has_legacy = any(detector_type == "Other_FS"),
    .groups = "drop"
  ) %>%
  filter(has_chorus & has_legacy) %>%
  pull(site_id)

dat <- filter(dat, site_id %in% paired_sites)

paired_nights <- dat %>%
  group_by(site_id, night) %>%
  summarise(n_det_types = n_distinct(detector_type), .groups = "drop") %>%
  filter(n_det_types == 2) %>%
  select(site_id, night)

dat <- inner_join(dat, paired_nights, by = c("site_id", "night"))

cat("Paired sites:", n_distinct(dat$site_id), "\n")
cat("Paired site-nights:", nrow(paired_nights), "\n")

write.csv(dat, "bat_detections.csv", row.names = FALSE)


# =============================================================================
# BUILD NIGHTLY COUNTS WITH ZEROS
# =============================================================================

all_combos <- dat %>%
  distinct(site_id, CellName, Quadrant, detector_model, detector_type, night)

spp_keep <- dat %>%
  group_by(spp) %>%
  summarise(total = n(),
            site_nights = n_distinct(paste(site_id, night)),
            .groups = "drop") %>%
  filter(total >= 50, site_nights >= 10) %>%
  arrange(desc(total)) %>%
  pull(spp)

print(spp_keep)

obs_counts <- dat %>%
  group_by(site_id, CellName, Quadrant, detector_model, detector_type, night, spp) %>%
  summarise(n_calls = n(), .groups = "drop")

night_counts <- crossing(all_combos, spp = spp_keep) %>%
  left_join(
    obs_counts,
    by = c("site_id", "CellName", "Quadrant",
           "detector_model", "detector_type", "night", "spp")
  ) %>%
  mutate(n_calls = replace_na(n_calls, 0))

cat("Total rows in model dataset (including zeros):", nrow(night_counts), "\n")


# =============================================================================
# ADD LEGACY MODEL LABEL  
# =============================================================================

site_legacy <- all_combos %>%
  filter(detector_model != "Anabat_Chorus") %>%
  distinct(site_id, detector_model) %>%
  arrange(site_id, detector_model) %>%
  group_by(site_id) %>%
  summarise(legacy_model    = paste(as.character(detector_model), collapse = " + "),
            n_legacy_models = n(),
            .groups = "drop")

mixed_sites <- filter(site_legacy, n_legacy_models > 1)
if (nrow(mixed_sites) > 0) {
  warning("These sites had more than one legacy model, which the design did not ",
          "anticipate. Check them before modelling:\n  ",
          paste(mixed_sites$site_id, collapse = ", "))
}

night_counts <- night_counts %>%
  left_join(select(site_legacy, site_id, legacy_model), by = "site_id") %>%
  mutate(legacy_model = factor(legacy_model),
         site_night   = paste(site_id, night, sep = "|"))


# =============================================================================
# DESIGN DIAGNOSTICS
# =============================================================================

cat("\n-- Bat passes by recorder model --\n")
print(table(dat$detector_model))

cat("\n-- Sites per legacy model --\n")
print(site_legacy %>% count(legacy_model, name = "n_sites"))

cat("\n-- Site-nights per legacy model --\n")
print(
  all_combos %>%
    distinct(site_id, night) %>%
    left_join(select(site_legacy, site_id, legacy_model), by = "site_id") %>%
    count(legacy_model, name = "n_site_nights")
)

cat("\n-- Sites x recorder model (deployment matrix) --\n")
print(table(all_combos$site_id, all_combos$detector_model))

n_legacy_types <- nlevels(droplevels(factor(site_legacy$legacy_model)))

cat("\n")
if (n_legacy_types >= 2) {
  cat(">> Both legacy models are represented. Script 02 will estimate a\n",
      "   Chorus/legacy correction factor and test whether that factor differs\n",
      "   between SM2 sites and Swift sites.\n")
} else {
  cat(">> Only one legacy model appears in the paired data. No SM2-vs-Swift\n",
      "   check is possible — report the correction factor as specific to that\n",
      "   unit rather than to legacy recorders in general.\n")
}

head(night_counts)
