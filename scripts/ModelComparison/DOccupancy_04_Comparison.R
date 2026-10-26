# Compare species-specific occupancy trends: Chorus vs SM2 data, BC (all sites)
# Input:  ALL_SPECIES_trends_Chorus.csv, ALL_SPECIES_trends_SM2.csv
# Output: fig_trend_BC_all_sites.png
#         Left panel  = annual slope in occupancy (psi per year) +/- 95% CrI
#         Right panel = overall change in occupancy (psi_last - psi_first) +/- 95% CrI

library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)

chorus_path <- "Data/Analyzed/DOccupancy/ChorusOutputs/tables/ALL_SPECIES_trends_Chorus.csv"
sm2_path    <- "Data/Analyzed/DOccupancy/SM2Outputs/tables/ALL_SPECIES_trends_SM2.csv"
focal_scope <- "BC (all sites)"

det_cols <- c(Chorus = "#1B7F79", SM2 = "#D9822B") #change as needed to make it fit what we are using in other sections of the study

trends <- bind_rows(
  read_csv(chorus_path, show_col_types = FALSE) |> mutate(detector = "Chorus"),
  read_csv(sm2_path,    show_col_types = FALSE) |> mutate(detector = "SM2")
) |>
  filter(scope == focal_scope, n_detections > 0)

# long format: one row per species x detector x metric (slope, change)
trends_long <- trends |>
  select(species, detector,
         slope_est = slope, slope_lo95, slope_hi95,
         change_est = change, change_lo95, change_hi95) |>
  pivot_longer(-c(species, detector),
               names_to = c("metric", ".value"),
               names_sep = "_") |>
  mutate(
    metric = factor(metric, levels = c("slope", "change"),
                    labels = c("Annual trend (slope of occupancy per year)",
                               "Overall change in occupancy (last - first year)")),
    clear = lo95 > 0 | hi95 < 0      # 95% CrI excludes zero
  )

# order species by mean slope across both detectors
sp_order <- trends |>
  group_by(species) |>
  summarise(m = mean(slope), .groups = "drop") |>
  arrange(m) |>
  pull(species)
trends_long <- trends_long |> mutate(species = factor(species, levels = sp_order))

pd <- position_dodge(width = 0.6)

p <- ggplot(trends_long, aes(x = est, y = species, colour = detector, group = detector)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(aes(xmin = lo95, xmax = hi95),
                 height = 0, linewidth = 0.8, position = pd) +
  geom_point(aes(shape = clear), size = 2.2, position = pd,
             fill = "white", stroke = 0.9) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                     labels = c(`TRUE`  = "95% CI excludes 0",
                                `FALSE` = "95% CI overlaps 0"),
                     breaks = c("TRUE", "FALSE"), name = NULL) +
  scale_colour_manual(values = det_cols, name = "Data source") +
  facet_wrap(~ metric, ncol = 2, scales = "free_x") +
  labs(x = "Estimate (95% CI)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "grey92"),
        strip.background = element_rect(fill = "grey95"),
        legend.position = "bottom")

ggsave("Figures/fig_trend_BC_all_sites.png", p, width = 9, height = 6, dpi = 600)