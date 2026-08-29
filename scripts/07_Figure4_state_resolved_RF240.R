# =========================================================
# 07_Figure4_state_resolved_RF240.R
# Purpose:
#   Generate an alternative Figure 4 with Panel D reduced to top 15
#   RF240 module-state enrichment patterns.
#
# This script does NOT change the analysis.
# It reuses 60B biological-label tables and regenerates:
#   - Panel D top 15
#   - Combined Figure 4 top 15
#   - Support unclipped Panel D top 15
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(scales)
  library(patchwork)
})

# =========================================================
# =========================================================
# 1. Paths and audited derived inputs
# =========================================================

script_file <- tryCatch(
  normalizePath(
    sys.frame(1)$ofile,
    winslash = "/",
    mustWork = TRUE
  ),
  error = function(e) NA_character_
)

if (is.na(script_file)) {
  stop("Run this script with source().")
}

script_dir <- dirname(script_file)

project_root <- normalizePath(
  file.path(script_dir, ".."),
  winslash = "/",
  mustWork = TRUE
)

dir_inputs <- file.path(
  project_root,
  "data",
  "longitudinal",
  "projected_states"
)

dir_outputs <- file.path(
  project_root,
  "results",
  "figures",
  "Figure4"
)

dir_plots <- dir_outputs

dir_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Figure4_projected_states"
)

dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

path_panelA_values <- file.path(dir_inputs,
  "Fig4A_plot_values_biological_labels.csv")
path_panelB_values <- file.path(dir_inputs,
  "Fig4B_heatmap_values_biological_labels.csv")
path_panelB_annot <- file.path(dir_inputs,
  "Fig4B_annotations_biological_labels.csv")
path_panelC_values <- file.path(dir_inputs,
  "Fig4C_RF240_global_enrichment_values_biological_labels.csv")
path_panelC_annot <- file.path(dir_inputs,
  "Fig4C_RF240_global_enrichment_annotations_biological_labels.csv")
path_panelD_values <- file.path(dir_inputs,
  "Fig4D_top_module_state_values_biological_labels_clipped.csv")
path_panelD_annot <- file.path(dir_inputs,
  "Fig4D_annotations_biological_labels.csv")
path_panelD_rank <- file.path(dir_inputs,
  "Fig4D_top_module_state_ranking_biological_labels.csv")
path_identity_map <- file.path(dir_inputs,
  "Fig4_projected_state_biological_labels_used_60B.csv")

input_check <- tibble(
  object = c(
    "panelA_values", "panelB_values", "panelB_annot",
    "panelC_values", "panelC_annot", "panelD_values",
    "panelD_annot", "panelD_rank", "identity_map"
  ),
  path = c(
    path_panelA_values, path_panelB_values, path_panelB_annot,
    path_panelC_values, path_panelC_annot, path_panelD_values,
    path_panelD_annot, path_panelD_rank, path_identity_map
  ),
  exists = file.exists(path)
)

print(input_check)

if (any(!input_check$exists)) {
  stop("Missing one or more audited Figure 4 derived inputs.")
}

# 3. Visual configuration
# =========================================================

time_response_levels <- c(
  "Base\nNR", "Base\nR",
  "PD1\nNR", "PD1\nR",
  "RTPD1\nNR", "RTPD1\nR"
)

heatmap_low <- "#2166AC"
heatmap_mid <- "white"
heatmap_high <- "#B2182B"

composition_low <- "white"
composition_high <- "#08306B"

panelD_clip_limits <- c(-0.75, 0.75)
max_panelD_rows_to_plot <- 15

# =========================================================
# 4. Helper functions
# =========================================================

save_panel <- function(plot_obj, filename, width, height, dpi = 600) {
  
  pdf_file <- file.path(dir_plots, paste0(filename, ".pdf"))
  png_file <- file.path(dir_plots, paste0(filename, ".png"))
  tiff_file <- file.path(dir_plots, paste0(filename, ".tiff"))
  
  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    device = "pdf"
  )
  
  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  
  ggsave(
    filename = tiff_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw"
  )
  
  message("Saved: ", basename(pdf_file))
  message("Saved: ", basename(png_file))
  message("Saved: ", basename(tiff_file))
}

theme_clean <- function(base_size = 10) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", color = "grey70"),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      axis.text.x = element_text(color = "black"),
      axis.text.y = element_text(color = "black")
    )
}

# =========================================================
# 5. Load 60B tables
# =========================================================

cat("\n[1] Loading 60B output tables...\n")

panelA_values <- read_csv(path_panelA_values, show_col_types = FALSE)
panelB_values <- read_csv(path_panelB_values, show_col_types = FALSE)
panelB_annot  <- read_csv(path_panelB_annot, show_col_types = FALSE)
panelC_values <- read_csv(path_panelC_values, show_col_types = FALSE)
panelC_annot  <- read_csv(path_panelC_annot, show_col_types = FALSE)
panelD_values <- read_csv(path_panelD_values, show_col_types = FALSE)
panelD_annot  <- read_csv(path_panelD_annot, show_col_types = FALSE)
panelD_rank   <- read_csv(path_panelD_rank, show_col_types = FALSE)
identity_map  <- read_csv(path_identity_map, show_col_types = FALSE)

# =========================================================
# 6. Factor handling and top 15 Panel D selection
# =========================================================

state_label_levels <- identity_map$state_label

cluster_colors <- setNames(
  grDevices::hcl.colors(length(state_label_levels), palette = "Dark 3"),
  state_label_levels
)

panelA_values <- panelA_values %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels),
    state_label = factor(state_label, levels = state_label_levels)
  )

panelB_values <- panelB_values %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels)
  )

panelB_annot <- panelB_annot %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels)
  )

panelC_values <- panelC_values %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels)
  )

panelC_annot <- panelC_annot %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels)
  )

# Select top 15 using the original ranked table
top_module_states <- panelD_rank %>%
  slice_head(n = max_panelD_rows_to_plot) %>%
  pull(module_state_bio)

panelD_values_top15 <- panelD_values %>%
  filter(module_state_bio %in% top_module_states) %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels),
    module_state_bio = factor(module_state_bio, levels = rev(top_module_states)),
    module_state_bio_wrapped = stringr::str_wrap(as.character(module_state_bio), width = 52),
    module_state_bio_wrapped = factor(
      module_state_bio_wrapped,
      levels = stringr::str_wrap(rev(top_module_states), width = 52)
    ),
    median_module_enrichment_clipped = pmin(
      pmax(median_module_enrichment, panelD_clip_limits[1]),
      panelD_clip_limits[2]
    )
  )

panelD_annot_top15 <- panelD_annot %>%
  filter(module_state_bio %in% top_module_states) %>%
  mutate(
    x_label = factor(x_label, levels = time_response_levels),
    module_state_bio = factor(module_state_bio, levels = rev(top_module_states)),
    module_state_bio_wrapped = stringr::str_wrap(as.character(module_state_bio), width = 52),
    module_state_bio_wrapped = factor(
      module_state_bio_wrapped,
      levels = stringr::str_wrap(rev(top_module_states), width = 52)
    )
  )

# Export top15 tables
write_csv(
  panelD_values_top15,
  file.path(dir_tables, "Fig4D_top15_module_state_values_biological_labels_clipped.csv")
)

write_csv(
  panelD_annot_top15,
  file.path(dir_tables, "Fig4D_top15_annotations_biological_labels.csv")
)

write_csv(
  panelD_rank %>% slice_head(n = max_panelD_rows_to_plot),
  file.path(dir_tables, "Fig4D_top15_module_state_ranking_biological_labels.csv")
)

# =========================================================
# 7. Panel D clipping audit top 15
# =========================================================

panelD_clipping_audit_top15 <- panelD_values_top15 %>%
  summarise(
    n_values = n(),
    min_raw = min(median_module_enrichment, na.rm = TRUE),
    max_raw = max(median_module_enrichment, na.rm = TRUE),
    n_values_below_clip = sum(median_module_enrichment < panelD_clip_limits[1], na.rm = TRUE),
    n_values_above_clip = sum(median_module_enrichment > panelD_clip_limits[2], na.rm = TRUE),
    percent_values_clipped = 100 * (
      n_values_below_clip + n_values_above_clip
    ) / n_values,
    clip_lower = panelD_clip_limits[1],
    clip_upper = panelD_clip_limits[2]
  )

panelD_extreme_values_top15 <- panelD_values_top15 %>%
  filter(
    median_module_enrichment < panelD_clip_limits[1] |
      median_module_enrichment > panelD_clip_limits[2]
  ) %>%
  arrange(desc(abs(median_module_enrichment))) %>%
  select(
    module,
    projected_cluster_label,
    state_label,
    module_state_bio,
    timepoint,
    response,
    x_label,
    median_module_enrichment,
    median_module_enrichment_clipped,
    support_flag
  )

write_csv(
  panelD_clipping_audit_top15,
  file.path(dir_tables, "Fig4D_top15_clipping_audit.csv")
)

write_csv(
  panelD_extreme_values_top15,
  file.path(dir_tables, "Fig4D_top15_values_affected_by_clipping.csv")
)

# =========================================================
# 8. Recreate Panels A-C for combined figure
# =========================================================

cat("\n[2] Plotting Panels A-C...\n")

pA <- ggplot(
  panelA_values,
  aes(x = x_label, y = mean_prop, fill = state_label)
) +
  geom_col(width = 0.85, color = "grey20", linewidth = 0.15) +
  scale_fill_manual(values = cluster_colors, name = "Projected state") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(
    title = "A. Projected T-cell state composition across therapy",
    x = NULL,
    y = "Mean proportion"
  ) +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 8),
    legend.position = "right",
    legend.text = element_text(size = 7.5),
    legend.key.size = unit(0.42, "cm")
  )

pB <- ggplot(
  panelB_values,
  aes(x = x_label, y = state_label_wrapped, fill = median_prop)
) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  geom_text(
    data = panelB_annot,
    aes(x = x_label, y = state_label_wrapped, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3
  ) +
  scale_fill_gradient(
    low = composition_low,
    high = composition_high,
    labels = scales::percent_format(accuracy = 1),
    name = "Median\nproportion"
  ) +
  labs(
    title = "B. Median projected-state proportions",
    x = NULL,
    y = NULL
  ) +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 7.2)
  )

pC <- ggplot(
  panelC_values,
  aes(x = x_label, y = state_label_wrapped, fill = median_global_enrichment)
) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  geom_text(
    data = panelC_annot,
    aes(x = x_label, y = state_label_wrapped, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3
  ) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "Median RF240\nenrichment"
  ) +
  labs(
    title = "C. Global RF240 enrichment within projected states",
    x = NULL,
    y = NULL
  ) +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 7.2)
  )

# Save A-C again for this version
save_panel(
  pA,
  "Fig4A_projected_state_mean_composition_biological_labels",
  width = 10.4,
  height = 4.6
)

save_panel(
  pB,
  "Fig4B_projected_state_median_proportion_heatmap_biological_labels",
  width = 8.8,
  height = 4.9
)

save_panel(
  pC,
  "Fig4C_RF240_global_enrichment_projected_states_biological_labels",
  width = 8.8,
  height = 4.9
)

# =========================================================
# 9. Panel D top 15 clipped
# =========================================================

cat("\n[3] Plotting Panel D top 15...\n")

pD <- ggplot(
  panelD_values_top15,
  aes(x = x_label, y = module_state_bio_wrapped, fill = median_module_enrichment)
) +
  geom_tile(color = "grey85", linewidth = 0.30) +
  geom_text(
    data = panelD_annot_top15,
    aes(x = x_label, y = module_state_bio_wrapped, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 2.6
  ) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    limits = panelD_clip_limits,
    oob = scales::squish,
    name = paste0(
      "Median module\nenrichment\nclipped [",
      panelD_clip_limits[1],
      ", ",
      panelD_clip_limits[2],
      "]"
    )
  ) +
  labs(
    title = "D. Top RF240 module-state enrichment patterns",
    x = NULL,
    y = NULL
  ) +
  theme_clean(9) +
  theme(
    axis.text.x = element_text(face = "bold", size = 7.5),
    axis.text.y = element_text(face = "bold", size = 6.7),
    legend.position = "right"
  )

save_panel(
  pD,
  "Fig4D_top15_RF240_module_state_enrichment_patterns_biological_labels_clipped",
  width = 11.8,
  height = 6.4
)

# =========================================================
# 10. Combined Figure 4 top 15
# =========================================================

cat("\n[4] Building combined Figure 4 top 15...\n")

combined_fig4_top15 <- (pA / (pB | pC) / pD) +
  patchwork::plot_layout(heights = c(0.95, 1.05, 1.30))

save_panel(
  combined_fig4_top15,
  "Fig4_combined_reworked_biological_labels_clippedD_top15",
  width = 14.0,
  height = 15.0
)

# =========================================================
# 11. Support Panel D top 15, unclipped
# =========================================================

cat("\n[5] Plotting support unclipped Panel D top 15...\n")

pD_unclipped <- ggplot(
  panelD_values_top15,
  aes(x = x_label, y = module_state_bio_wrapped, fill = median_module_enrichment)
) +
  geom_tile(color = "grey85", linewidth = 0.30) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "Median module\nenrichment"
  ) +
  labs(
    title = "Support. Top RF240 module-state enrichment patterns, unclipped scale",
    x = NULL,
    y = NULL
  ) +
  theme_clean(9) +
  theme(
    axis.text.x = element_text(face = "bold", size = 7.5),
    axis.text.y = element_text(face = "bold", size = 6.7),
    legend.position = "right"
  )

save_panel(
  pD_unclipped,
  "Support_Fig4D_top15_RF240_module_state_enrichment_patterns_unclipped",
  width = 11.8,
  height = 6.4
)

# =========================================================
# 12. Manifest and summary
# =========================================================

cat("\n[6] Building manifest...\n")

manifest_tbl <- tibble(
  file = c(
    list.files(dir_tables, full.names = TRUE),
    list.files(dir_plots, full.names = TRUE)
  )
) %>%
  mutate(
    relative_file = fs::path_rel(file, start = dir_outputs),
    bytes = file.info(file)$size
  ) %>%
  arrange(relative_file)

write_csv(
  manifest_tbl,
  file.path(dir_tables, "07_Figure4_manifest_outputs.csv")
)

cat("\n\n=========================================================\n")
cat("07 DONE - Figure 4 projected-state RF240 analysis\n")
cat("Output directory:\n")
cat(dir_outputs, "\n")
cat("=========================================================\n")

cat("\n[Panel D top15 clipping audit]\n")
print(panelD_clipping_audit_top15)

cat("\n[Panel D top15 values affected by clipping]\n")
print(panelD_extreme_values_top15, n = 80)

cat("\n[Top 15 module-state rows shown in Panel D]\n")
print(
  panelD_rank %>%
    slice_head(n = max_panelD_rows_to_plot) %>%
    select(
      module,
      projected_cluster_label,
      state_label,
      module_state_bio,
      max_abs_median_enrichment,
      dynamic_range,
      mean_abs_median_enrichment
    ),
  n = 30
)

cat("\n[Generated plot files]\n")
print(list.files(dir_plots, full.names = FALSE))

cat("\n[Generated table files]\n")
print(list.files(dir_tables, full.names = FALSE))

cat("\n=========================================================\n")
cat("END OF 60C SUMMARY\n")
cat("=========================================================\n")
