library(tidyverse)
library(vegan)

# ---- 1. Reshape to wide: one row per site x night x detector ----
spp_cols <- night_counts %>% distinct(spp) %>% pull(spp) %>% sort()

wide <- night_counts %>%
  mutate(sample_id = paste(site_id, night, sep = "_")) %>%
  pivot_wider(id_cols = c(sample_id, site_id, detector_type),
              names_from = spp, values_from = n_calls, values_fill = 0)

# Split by detector
chorus <- wide %>% filter(detector_type == "Anabat_Chorus")
sm2    <- wide %>% filter(detector_type == "Other_FS")

# ---- 2. Paired Jaccard & Sørensen similarity ----
paired <- inner_join(chorus, sm2, by = "sample_id", suffix = c("_ch", "_sm"))

similarity <- map_dfr(1:nrow(paired), function(i) {
  a <- as.numeric(paired[i, paste0(spp_cols, "_ch")] > 0)
  b <- as.numeric(paired[i, paste0(spp_cols, "_sm")] > 0)
  shared <- sum(a & b)
  either <- sum(a | b)
  tibble(
    sample_id = paired$sample_id[i],
    site_id   = paired$site_id_ch[i],
    jaccard   = ifelse(either > 0, shared / either, NA),
    sorensen  = ifelse((sum(a) + sum(b)) > 0,
                       2 * shared / (sum(a) + sum(b)), NA)
  )
})

# Quick summary
summary(similarity$jaccard)
summary(similarity$sorensen)

# Plot similarity distributions
similarity %>%
  pivot_longer(cols = c(jaccard, sorensen), names_to = "index") %>%
  ggplot(aes(x = index, y = value)) +
  geom_boxplot(width = 0.4) +
  geom_jitter(width = 0.1, alpha = 0.3) +
  labs(y = "Similarity (0 = no overlap, 1 = identical)",
       x = NULL, title = "Species overlap between Chorus & SM2-Swift") +
  theme_minimal()

# Optionally: similarity by site
ggplot(similarity, aes(x = reorder(site_id, jaccard, median), y = jaccard)) +
  geom_boxplot() +
  coord_flip() +
  labs(x = NULL, y = "Jaccard similarity",
       title = "Jaccard similarity by site") +
  theme_minimal()

# ---- Evenness (Pielou's J) by detector type ----
spp_matrix_ev <- wide %>% select(all_of(spp_cols)) %>% as.matrix()

evenness <- wide %>%
  select(sample_id, site_id, detector_type) %>%
  mutate(
    H = diversity(spp_matrix_ev, index = "shannon"),
    S = specnumber(spp_matrix_ev),
    J = ifelse(S > 1, H / log(S), NA)
  )

# Boxplot
ggplot(evenness, aes(x = detector_type, y = J, fill = detector_type)) +
  geom_boxplot(width = 0.4, alpha = 0.7) +
  labs(y = "Pielou's Evenness (J)", x = NULL,
       title = "Evenness of bat community by detector type") +
  theme_minimal() +
  theme(legend.position = "none")

# Paired Wilcoxon test
evenness_wide <- evenness %>%
  select(sample_id, detector_type, J) %>%
  pivot_wider(names_from = detector_type, values_from = J)

wilcox.test(evenness_wide$Anabat_Chorus, evenness_wide$Other_FS, paired = TRUE)


# ---- 4. NMDS ordination ----
spp_matrix <- wide %>% select(all_of(spp_cols)) %>% as.matrix()
bc_dist <- vegdist(spp_matrix, method = "bray")

nmds <- metaMDS(bc_dist, k = 2, trymax = 100)

nmds_df <- as_tibble(scores(nmds, display = "sites")) %>%
  bind_cols(wide %>% select(sample_id, site_id, detector_type))

ggplot(nmds_df, aes(x = NMDS1, y = NMDS2, color = detector_type)) +
  geom_point(size = 2, alpha = 0.7) +
  stat_ellipse(level = 0.95) +
  labs(title = paste("NMDS of bat communities (Stress:", 
                     round(nmds$stress, 3), ")"),
       color = "Detector") +
  theme_minimal()

# ---- 5. PERMANOVA ----
library(permute)

# Restrict permutations to within site_id (respects paired design)
perm <- how(nperm = 999, blocks = as.factor(wide$site_id))

adonis2(bc_dist ~ detector_type,
        data = wide, permutations = perm)


rank_abund <- night_counts %>%
  group_by(detector_type, spp) %>%
  summarise(total_calls = sum(n_calls), .groups = "drop") %>%
  group_by(detector_type) %>%
  mutate(
    rank = rank(-total_calls, ties.method = "first"),
    prop = total_calls / sum(total_calls)
  )

ggplot(rank_abund, aes(x = rank, y = prop, colour = detector_type)) +
  geom_line() + geom_point() +
  geom_text(aes(label = spp), size = 2.5, hjust = -0.1, show.legend = FALSE) +
  scale_y_log10() +
  labs(x = "Rank", y = "Proportion of calls (log scale)",
       title = "Rank-abundance by detector type") +
  theme_minimal()




