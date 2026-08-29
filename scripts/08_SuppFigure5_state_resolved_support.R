# =========================================================
# 08_SuppFigure5_state_resolved_support.R
# Purpose:
#   Reproduce/update Supplementary Figure 5 associated with Figure 4.
#
# Concept:
#   Extended state-resolved longitudinal context for RF240 organization.
#
# This script uses audited Block 4 outputs from 60A/60B/60C.
#
# Panels:
#   A. Patient-level projected T-cell state composition.
#   B. Projected-state proportion deltas from baseline.
#   C. RF240 global score within projected T-cell states.
#   D. RF240 global enrichment within projected T-cell states.
#   E. RF240 global enrichment deltas from baseline.
#   F. Full RF240 module-state enrichment heatmap.
#
# Notes:
#   - Supplementary Figure 4 remains associated with Figure 3.
#   - Supplementary Figure 5 is associated with Figure 4.
#   - Projected states use biological labels confirmed from Figure 1 identity map.
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
  "supplementary_figures",
  "Supplementary_Figure5"
)

dir_plots <- dir_outputs

dir_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Supplementary_Figure5"
)

dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

path_composition <- file.path(dir_inputs,
  "Fig4_composition_zero_completed_standardized.csv")
path_rf240 <- file.path(dir_inputs,
  "Fig4_RF240_projected_enrichment_standardized.csv")
path_composition_delta <- file.path(dir_inputs,
  "Fig4_projected_cluster_proportion_delta_from_baseline.csv")
path_rf240_global_delta <- file.path(dir_inputs,
  "Fig4_RF240_global_enrichment_delta_from_baseline.csv")
path_module_delta <- file.path(dir_inputs,
  "Fig4_RF240_module_enrichment_delta_from_baseline.csv")
path_identity_map <- file.path(dir_inputs,
  "Fig4_projected_state_biological_labels_used_60B.csv")

input_check <- tibble(
  object = c(
    "composition", "rf240", "composition_delta",
    "rf240_global_delta", "module_delta", "identity_map"
  ),
  path = c(
    path_composition, path_rf240, path_composition_delta,
    path_rf240_global_delta, path_module_delta, path_identity_map
  ),
  exists = file.exists(path)
)

print(input_check)

if (any(!input_check$exists)) {
  stop("Missing one or more audited Supplementary Figure 5 derived inputs.")
}

# 2. Visual configuration
# =========================================================

timepoint_levels <- c("Base", "PD1", "RTPD1")
response_levels <- c("NR", "R")
delta_type_levels <- c("PD1 - Base", "RTPD1 - Base")

time_response_levels <- c(
  "Base\nNR", "Base\nR",
  "PD1\nNR", "PD1\nR",
  "RTPD1\nNR", "RTPD1\nR"
)

response_colors <- c(
  "NR" = "#E69F00",
  "R" = "#0072B2"
)

heatmap_low <- "#2166AC"
heatmap_mid <- "white"
heatmap_high <- "#B2182B"

module_heatmap_clip_limits <- c(-0.75, 0.75)

# =========================================================
# 3. Helper functions
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
# 4. Load audited tables
# =========================================================

cat("\n[1] Loading audited tables...\n")

composition <- read_csv(path_composition, show_col_types = FALSE)
rf240 <- read_csv(path_rf240, show_col_types = FALSE)
composition_delta <- read_csv(path_composition_delta, show_col_types = FALSE)
rf240_global_delta <- read_csv(path_rf240_global_delta, show_col_types = FALSE)
module_delta <- read_csv(path_module_delta, show_col_types = FALSE)
identity_map <- read_csv(path_identity_map, show_col_types = FALSE)

state_label_levels <- identity_map$state_label
state_label_wrapped_levels <- identity_map$state_label_wrapped

cluster_colors <- setNames(
  grDevices::hcl.colors(length(state_label_levels), palette = "Dark 3"),
  state_label_levels
)

# =========================================================
# 5. Add labels and standardize variables
# =========================================================

cat("\n[2] Adding biological labels...\n")

composition_labeled <- composition %>%
  mutate(
    projected_cluster_label = as.character(projected_cluster_label),
    timepoint = factor(as.character(timepoint), levels = timepoint_levels),
    response = factor(as.character(response), levels = response_levels),
    patient_id = as.character(patient_id)
  ) %>%
  left_join(
    identity_map %>%
      select(projected_cluster_label, state_label, state_label_wrapped),
    by = "projected_cluster_label"
  ) %>%
  mutate(
    state_label = factor(state_label, levels = state_label_levels),
    state_label_wrapped = factor(state_label_wrapped, levels = state_label_wrapped_levels),
    time_response = factor(
      paste(timepoint, response, sep = "\n"),
      levels = time_response_levels
    )
  )

composition_delta_labeled <- composition_delta %>%
  mutate(
    projected_cluster_label = as.character(projected_cluster_label),
    delta_type = factor(as.character(delta_type), levels = delta_type_levels),
    response = factor(as.character(response), levels = response_levels),
    patient_id = as.character(patient_id)
  ) %>%
  left_join(
    identity_map %>%
      select(projected_cluster_label, state_label, state_label_wrapped),
    by = "projected_cluster_label"
  ) %>%
  mutate(
    state_label = factor(state_label, levels = state_label_levels),
    state_label_wrapped = factor(state_label_wrapped, levels = state_label_wrapped_levels)
  )

rf240_labeled <- rf240 %>%
  mutate(
    projected_cluster_label = as.character(projected_cluster_label),
    timepoint = factor(as.character(timepoint), levels = timepoint_levels),
    response = factor(as.character(response), levels = response_levels),
    patient_id = as.character(patient_id),
    module = as.character(module)
  ) %>%
  left_join(
    identity_map %>%
      select(projected_cluster_label, state_label, state_label_wrapped),
    by = "projected_cluster_label"
  ) %>%
  mutate(
    state_label = factor(state_label, levels = state_label_levels),
    state_label_wrapped = factor(state_label_wrapped, levels = state_label_wrapped_levels),
    time_response = factor(
      paste(timepoint, response, sep = "\n"),
      levels = time_response_levels
    )
  )

rf240_global <- rf240_labeled %>%
  filter(module == "RF240 global")

rf240_global_delta_labeled <- rf240_global_delta %>%
  mutate(
    projected_cluster_label = as.character(projected_cluster_label),
    delta_type = factor(as.character(delta_type), levels = delta_type_levels),
    response = factor(as.character(response), levels = response_levels),
    patient_id = as.character(patient_id)
  ) %>%
  left_join(
    identity_map %>%
      select(projected_cluster_label, state_label, state_label_wrapped),
    by = "projected_cluster_label"
  ) %>%
  mutate(
    state_label = factor(state_label, levels = state_label_levels),
    state_label_wrapped = factor(state_label_wrapped, levels = state_label_wrapped_levels)
  )

module_delta_labeled <- module_delta %>%
  mutate(
    projected_cluster_label = as.character(projected_cluster_label),
    delta_type = factor(as.character(delta_type), levels = delta_type_levels),
    response = factor(as.character(response), levels = response_levels),
    patient_id = as.character(patient_id),
    module = as.character(module)
  ) %>%
  left_join(
    identity_map %>%
      select(projected_cluster_label, state_label, state_label_wrapped),
    by = "projected_cluster_label"
  ) %>%
  mutate(
    state_label = factor(state_label, levels = state_label_levels),
    state_label_wrapped = factor(state_label_wrapped, levels = state_label_wrapped_levels)
  )

# Export labeled source tables for Supp Figure 5
write_csv(
  composition_labeled,
  file.path(dir_tables, "SuppFig5_composition_labeled.csv")
)

write_csv(
  composition_delta_labeled,
  file.path(dir_tables, "SuppFig5_composition_delta_labeled.csv")
)

write_csv(
  rf240_labeled,
  file.path(dir_tables, "SuppFig5_RF240_state_resolved_labeled.csv")
)

write_csv(
  rf240_global_delta_labeled,
  file.path(dir_tables, "SuppFig5_RF240_global_enrichment_delta_labeled.csv")
)

write_csv(
  module_delta_labeled,
  file.path(dir_tables, "SuppFig5_RF240_module_delta_labeled.csv")
)

# =========================================================
# 6. Panel A: patient-level projected-state composition
# =========================================================

cat("\n[3] Plotting Panel A...\n")

panelA_data <- composition_labeled %>%
  mutate(
    patient_time_label = paste(patient_id, timepoint, response, sep = " | ")
  ) %>%
  group_by(timepoint, response) %>%
  mutate(
    patient_time_label = factor(
      patient_time_label,
      levels = unique(patient_time_label[order(patient_id)])
    )
  ) %>%
  ungroup()

pA <- ggplot(
  panelA_data,
  aes(x = patient_time_label, y = prop_cluster, fill = state_label)
) +
  geom_col(width = 0.9, color = NA) +
  facet_grid(response ~ timepoint, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = cluster_colors, name = "Projected state") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "A. Patient-level projected-state composition",
    x = NULL,
    y = "Projected-state proportion"
  ) +
  theme_clean(9) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.35, "cm")
  )

save_panel(
  pA,
  "SuppFig5A_patient_level_projected_state_composition",
  width = 13.0,
  height = 5.2
)

# =========================================================
# 7. Panel B: state proportion deltas from baseline
# =========================================================

cat("\n[4] Plotting Panel B...\n")

pB <- ggplot(
  composition_delta_labeled,
  aes(x = response, y = delta_prop_cluster, fill = response)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.75) +
  facet_grid(state_label_wrapped ~ delta_type, scales = "free_y") +
  scale_fill_manual(values = response_colors, guide = "none") +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "B. Projected-state proportion changes from baseline",
    x = NULL,
    y = "Delta projected-state proportion"
  ) +
  theme_clean(8.5) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 6.5)
  )

save_panel(
  pB,
  "SuppFig5B_projected_state_proportion_delta_from_baseline",
  width = 10.8,
  height = 10.0
)

# =========================================================
# 8. Panel C: global RF240 score within projected states
# =========================================================

cat("\n[5] Plotting Panel C...\n")

pC <- ggplot(
  rf240_global,
  aes(x = response, y = module_score_cluster, fill = response)
) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 0.75, alpha = 0.70) +
  facet_grid(state_label_wrapped ~ timepoint, scales = "free_y") +
  scale_fill_manual(values = response_colors, guide = "none") +
  labs(
    title = "C. Global RF240 score within projected states",
    x = NULL,
    y = "RF240 global score"
  ) +
  theme_clean(8.5) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 6.5)
  )

save_panel(
  pC,
  "SuppFig5C_RF240_global_score_within_projected_states",
  width = 10.8,
  height = 10.0
)

# =========================================================
# 9. Panel D: global RF240 enrichment within projected states
# =========================================================

cat("\n[6] Plotting Panel D...\n")

pD <- ggplot(
  rf240_global,
  aes(x = response, y = module_enrichment_delta, fill = response)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 0.75, alpha = 0.70) +
  facet_grid(state_label_wrapped ~ timepoint, scales = "free_y") +
  scale_fill_manual(values = response_colors, guide = "none") +
  labs(
    title = "D. Global RF240 enrichment within projected states",
    x = NULL,
    y = "RF240 enrichment\nstate score - patient-timepoint average"
  ) +
  theme_clean(8.5) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 6.5)
  )

save_panel(
  pD,
  "SuppFig5D_RF240_global_enrichment_within_projected_states",
  width = 10.8,
  height = 10.0
)

# =========================================================
# 10. Panel E: RF240 global enrichment deltas from baseline
# =========================================================

cat("\n[7] Plotting Panel E...\n")

pE <- ggplot(
  rf240_global_delta_labeled,
  aes(x = response, y = delta_global_RF240_enrichment, fill = response)
) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.35) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.75) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.75) +
  facet_grid(state_label_wrapped ~ delta_type, scales = "free_y") +
  scale_fill_manual(values = response_colors, guide = "none") +
  labs(
    title = "E. RF240 enrichment changes from baseline",
    x = NULL,
    y = "Delta RF240 enrichment"
  ) +
  theme_clean(8.5) +
  theme(
    axis.text.x = element_text(face = "bold"),
    axis.text.y = element_text(size = 6),
    strip.text.y = element_text(size = 6.5)
  )

save_panel(
  pE,
  "SuppFig5E_RF240_global_enrichment_delta_from_baseline",
  width = 10.8,
  height = 10.0
)

# =========================================================
# 11. Panel F: full module-state enrichment heatmap
# =========================================================

cat("\n[8] Plotting Panel F...\n")

module_heatmap_data <- rf240_labeled %>%
  filter(module != "RF240 global") %>%
  group_by(module, state_label, timepoint, response, time_response) %>%
  summarise(
    median_module_enrichment = median(module_enrichment_delta, na.rm = TRUE),
    median_n_cells_cluster = median(n_cells_cluster, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    module_state = paste0(module, " | ", state_label),
    module_state = stringr::str_wrap(module_state, width = 55)
  )

module_state_order <- module_heatmap_data %>%
  group_by(module_state) %>%
  summarise(
    max_abs_enrichment = max(abs(median_module_enrichment), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(max_abs_enrichment) %>%
  pull(module_state)

module_heatmap_data <- module_heatmap_data %>%
  mutate(
    module_state = factor(module_state, levels = module_state_order),
    median_module_enrichment_clipped = pmin(
      pmax(median_module_enrichment, module_heatmap_clip_limits[1]),
      module_heatmap_clip_limits[2]
    )
  )

pF <- ggplot(
  module_heatmap_data,
  aes(x = time_response, y = module_state, fill = median_module_enrichment)
) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    limits = module_heatmap_clip_limits,
    oob = scales::squish,
    name = paste0(
      "Median module\nenrichment\nclipped [",
      module_heatmap_clip_limits[1],
      ", ",
      module_heatmap_clip_limits[2],
      "]"
    )
  ) +
  labs(
    title = "F. Full RF240 module enrichment across projected states",
    x = NULL,
    y = NULL
  ) +
  theme_clean(8.2) +
  theme(
    axis.text.x = element_text(face = "bold", size = 7),
    axis.text.y = element_text(size = 5.2)
  )

save_panel(
  pF,
  "SuppFig5F_full_RF240_module_enrichment_projected_states",
  width = 11.8,
  height = 13.5
)

write_csv(
  module_heatmap_data,
  file.path(dir_tables, "SuppFig5F_full_module_enrichment_heatmap_values.csv")
)

# =========================================================
# 12. Combined supplementary figure
# =========================================================

cat("\n[9] Building combined Supplementary Figure 5...\n")

supp_fig5_combined <- (pA / pB / pC / pD / pE / pF) +
  patchwork::plot_layout(
    heights = c(0.75, 1.35, 1.35, 1.35, 1.35, 1.80)
  )

save_panel(
  supp_fig5_combined,
  "Supplementary_Figure_5_reworked_state_resolved_RF240_extended_context",
  width = 14.0,
  height = 40.0
)

# =========================================================
# 13. Manifest and summary
# =========================================================

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
  file.path(dir_tables, "08_SuppFigure5_manifest_outputs.csv")
)

cat("\n\n=========================================================\n")
cat("08 DONE - Supplementary Figure 5 state-resolved RF240 support\n")
cat("Output directory:\n")
cat(dir_outputs, "\n")
cat("=========================================================\n")

cat("\n[Input dimensions]\n")
cat("composition_labeled: "); print(dim(composition_labeled))
cat("composition_delta_labeled: "); print(dim(composition_delta_labeled))
cat("rf240_labeled: "); print(dim(rf240_labeled))
cat("rf240_global: "); print(dim(rf240_global))
cat("rf240_global_delta_labeled: "); print(dim(rf240_global_delta_labeled))
cat("module_delta_labeled: "); print(dim(module_delta_labeled))

cat("\n[Biological labels used]\n")
print(identity_map, n = 20)

cat("\n[Generated plot files]\n")
print(list.files(dir_plots, full.names = FALSE))

cat("\n[Generated table files]\n")
print(list.files(dir_tables, full.names = FALSE))

cat("\n=========================================================\n")
cat("END OF 60D SUMMARY\n")
cat("=========================================================\n")
