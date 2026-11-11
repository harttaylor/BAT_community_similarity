# =============================================================================
# QUESTION 3: At the community level, do the Chorus and the SM2 record the
# same bats?
#
# Kaleidoscope IDs only.
#
# OUTPUTS (results/ folder)
#   fig5_rank_abundance.png/.pdf Fig. 5: rank-abundance curves
#   figS1_nmds_composition.png   Fig. S1: NMDS of proportions
#   figS_nmds_raw.png            NMDS of raw counts (for reference; not in the paper)
#   Test results are printed to the console.
#
# METHODS DETAILS
#   Diversity: richness and Pielou's evenness (J) per site-night, averaged per
#     site and unit; paired Wilcoxon tests with sites as pairs.
#   Species overlap: Sorensen per site-night (descriptive). Jaccard is still
#     calculated but not reported (it is a direct transformation of Sorensen).
#   Composition: PERMANOVA on Bray-Curtis of proportions (relative activity)
#     and of raw counts (composition + activity). Passes summed to one
#     community per unit per site. Model: dissimilarity ~ site_id +
#     detector_type, terms tested sequentially, 999 permutations restricted
#     within site. Report detector R2, F, p and site R2.
#   Dispersion: betadisper + permutest (p > 0.05 = assumption met).
#   NMDS: each site's two units joined by a line (short line = similar).
# =============================================================================

library(tidyverse)
library(vegan)
library(permute)
library(ggrepel)

setwd("C:/Users/hartt/Documents/Projects/BAT_detector_differences/BAT_community_similarity")
dir.create("results", showWarnings = FALSE)

bats <- readRDS("data/processed/bat_data_kaleidoscope.rds")
night_counts <- bats$night_counts
spp_keep     <- bats$spp_keep

# Detector colours ("Legacy" is the SM2)
det_colours <- c("Chorus" = "#5D3A9B", "Legacy" = "#E66100")
det_labels  <- c("Chorus" = "Chorus", "Legacy" = "SM2")


# -----------------------------------------------------------------------------
# WIDE TABLE: one row per site-night per detector, one column per species
# -----------------------------------------------------------------------------

wide <- night_counts %>%
  pivot_wider(id_cols = c(site_id, night, detector_type),
              names_from = spp, values_from = n_calls)
mat <- as.matrix(select(wide, all_of(spp_keep)))


# -----------------------------------------------------------------------------
# DIVERSITY: RICHNESS AND EVENNESS
# -----------------------------------------------------------------------------

wide$richness <- specnumber(mat)
wide$evenness <- diversity(mat, index = "shannon") / log(wide$richness)
wide$evenness[wide$richness < 2] <- NA            # evenness undefined for < 2 species

site_div <- wide %>%
  group_by(site_id, detector_type) %>%
  summarise(richness = mean(richness), evenness = mean(evenness, na.rm = TRUE),
            .groups = "drop") %>%
  pivot_wider(names_from = detector_type, values_from = c(richness, evenness))

print(summary(select(site_div, -site_id)))
print(wilcox.test(site_div$richness_Chorus, site_div$richness_Legacy, paired = TRUE))
print(wilcox.test(site_div$evenness_Chorus, site_div$evenness_Legacy, paired = TRUE))


# -----------------------------------------------------------------------------
# SPECIES OVERLAP, per site-night
# -----------------------------------------------------------------------------
# Sorensen = 2 x shared / (species on Chorus + species on SM2)
# Jaccard  = shared / species recorded by either unit (not reported)

overlap <- night_counts %>%
  mutate(present = n_calls > 0) %>%
  select(site_id, night, spp, detector_type, present) %>%
  pivot_wider(names_from = detector_type, values_from = present) %>%
  group_by(site_id, night) %>%
  summarise(shared   = sum(Chorus & Legacy),
            either   = sum(Chorus | Legacy),
            jaccard  = shared / either,
            sorensen = 2 * shared / (sum(Chorus) + sum(Legacy)),
            .groups = "drop")

print(summary(select(overlap, sorensen, jaccard)))   # min/max = range for the results


# -----------------------------------------------------------------------------
# ONE COMMUNITY PER DETECTOR PER SITE (for PERMANOVA and NMDS)
# -----------------------------------------------------------------------------

agg <- wide %>%
  group_by(site_id, detector_type) %>%
  summarise(across(all_of(spp_keep), sum), .groups = "drop")
agg_mat <- as.matrix(select(agg, all_of(spp_keep)))
perm <- how(nperm = 999, blocks = factor(agg$site_id))

d_rel <- vegdist(decostand(agg_mat, method = "total"), method = "bray")   # proportions
d_raw <- vegdist(agg_mat, method = "bray")                                # raw counts


# -----------------------------------------------------------------------------
# PERMANOVA, proportions (relative activity) + dispersion check
# -----------------------------------------------------------------------------

print(adonis2(d_rel ~ site_id + detector_type, data = agg, permutations = perm, by = "terms"))
print(permutest(betadisper(d_rel, agg$detector_type), permutations = 999))   # p > 0.05 = assumption met


# -----------------------------------------------------------------------------
# PERMANOVA, raw counts (composition + activity) + dispersion check
# -----------------------------------------------------------------------------

print(adonis2(d_raw ~ site_id + detector_type, data = agg, permutations = perm, by = "terms"))
print(permutest(betadisper(d_raw, agg$detector_type), permutations = 999))   # p > 0.05 = assumption met


# -----------------------------------------------------------------------------
# FIGURE 5: rank-abundance (share of all passes per species, by detector)
# -----------------------------------------------------------------------------

rank_ab <- night_counts %>%
  group_by(detector_type, spp) %>%
  summarise(total = sum(n_calls), .groups = "drop") %>%
  filter(total > 0) %>%
  group_by(detector_type) %>%
  mutate(rank = rank(-total, ties.method = "first"),
         prop = total / sum(total)) %>%
  ungroup()

fig5 <- ggplot(rank_ab, aes(rank, prop, colour = detector_type)) +
  geom_line(linewidth = 0.5) +
  geom_point(size = 2) +
  geom_text_repel(data = filter(rank_ab, detector_type == "Chorus"), aes(label = spp),
                  colour = "grey25", size = 2.8, seed = 1, nudge_y = 0.15,
                  segment.colour = "grey75", max.overlaps = Inf, show.legend = FALSE) +
  scale_y_log10() +
  scale_x_continuous(breaks = seq(2, 16, 2)) +
  scale_colour_manual(values = det_colours, labels = det_labels,
                      breaks = c("Chorus", "Legacy"), name = NULL) +
  labs(x = "Species rank", y = "Proportion of total passes") +
  theme_classic(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.85, 0.9))

ggsave("results/fig5_rank_abundance.png", fig5, width = 7, height = 5.5, dpi = 600)
ggsave("results/fig5_rank_abundance.pdf", fig5, width = 7, height = 5.5)


# -----------------------------------------------------------------------------
# FIGURE S1: NMDS, proportions, paired lines (Kamloops_NW labelled)
# -----------------------------------------------------------------------------

set.seed(1)
nmds_rel <- metaMDS(d_rel, k = 2, trymax = 100, trace = 0)
print(nmds_rel$stress)   # report in the figure caption

scores_rel <- as_tibble(scores(nmds_rel, display = "sites")) %>%
  bind_cols(select(agg, site_id, detector_type))

fig_nmds <- ggplot(scores_rel, aes(NMDS1, NMDS2)) +
  geom_line(aes(group = site_id), colour = "grey70", linewidth = 0.4) +
  geom_point(aes(colour = detector_type), size = 2.6) +
  geom_text(data = ~ filter(.x, site_id == "Kamloops_NW", detector_type == "Chorus"),
            label = "Kamloops_NW", hjust = -0.15, size = 3, colour = "grey30") +
  scale_colour_manual(values = det_colours, labels = det_labels,
                      breaks = c("Chorus", "Legacy"), name = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "inside", legend.position.inside = c(0.9, 0.92), aspect.ratio = 1)

ggsave("results/figS1_nmds_composition.png", fig_nmds, width = 6, height = 6, dpi = 600)
ggsave("results/figS1_nmds_composition.pdf", fig_nmds, width = 6, height = 6)


# -----------------------------------------------------------------------------
# NMDS, raw counts (composition + activity): for reference, not in the paper
# -----------------------------------------------------------------------------

set.seed(1)
nmds_raw <- metaMDS(d_raw, k = 2, trymax = 100, trace = 0)
print(nmds_raw$stress)

scores_raw <- as_tibble(scores(nmds_raw, display = "sites")) %>%
  bind_cols(select(agg, site_id, detector_type))

ggsave("results/figS_nmds_raw.png", fig_nmds %+% scores_raw, width = 6, height = 6, dpi = 600)


# -----------------------------------------------------------------------------
# CHECK: the long NMDS line - which site and why
# -----------------------------------------------------------------------------

scores_rel %>%
  pivot_wider(id_cols = site_id, names_from = detector_type, values_from = c(NMDS1, NMDS2)) %>%
  mutate(length = sqrt((NMDS1_Chorus - NMDS1_Legacy)^2 + (NMDS2_Chorus - NMDS2_Legacy)^2)) %>%
  arrange(desc(length))

night_counts %>% filter(site_id == "Kamloops_NW") %>%
  group_by(spp, detector_type) %>% summarise(passes = sum(n_calls), .groups = "drop") %>%
  pivot_wider(names_from = detector_type, values_from = passes)

