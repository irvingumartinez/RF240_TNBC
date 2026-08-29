# =========================================================
# 10_SuppFigure6_repertoire_structure.R
# Purpose:
#   Generate a cleaner, publication-oriented version of
#   candidate Supplementary Figure 6.
#
# Input:
#   Uses corrected outputs from 61B_2b plus original source tables.
#
# Main decision:
#   Statistical tests were performed in 61B_2b/61B_2c.
#   No comparison reached BH-adjusted q < 0.10.
#   Therefore, this clean version does not display significance symbols.
#
# Output:
#   Supplementary Figure 6 candidate, clean readable version.
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(patchwork)
})

# =========================================================
# =========================================================
# 1. Paths
# =========================================================

# Repository root is inferred from this script location.

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

dir_source_tables <- file.path(
  project_root,
  "data",
  "tcr",
  "supplementary_figure6"
)

# Corrected 61B_2b zero-completed tables are stored
# together with the other audited Supplementary Figure 6
# checkpoints in the public repository.
dir_tables_v2 <- dir_source_tables

dir_output_plots <- file.path(
  project_root,
  "results",
  "supplementary_figures",
  "Supplementary_Figure6"
)

dir_output_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Supplementary_Figure6"
)

dir_create(dir_output_plots)
dir_create(dir_output_tables)

# 2. Load data
# =========================================================

public_pt_zero <- read_csv(
  file.path(dir_tables_v2, "SuppFig6AB_public_fraction_zero_completed_patient_time.csv"),
  show_col_types = FALSE
)

public_state_zero <- read_csv(
  file.path(dir_tables_v2, "SuppFig6D_public_fraction_zero_completed_patient_time_state.csv"),
  show_col_types = FALSE
)

clonal_diversity <- read_csv(
  file.path(dir_source_tables, "clonal_diversity_patient_time.csv"),
  show_col_types = FALSE
)

rf240_public_private <- read_csv(
  file.path(dir_source_tables, "RF240_module_scores_public_private_patient_time_cluster.csv"),
  show_col_types = FALSE
)

rf240_expansion_modules <- read_csv(
  file.path(dir_source_tables, "RF240_module_scores_by_expansion_category_patient_time_cluster.csv"),
  show_col_types = FALSE
)

stats_summary <- read_csv(
  file.path(dir_tables_v2, "SuppFig6_v2_statistics_summary_by_panel.csv"),
  show_col_types = FALSE
)

# =========================================================
# 3. Labels
# =========================================================

time_order <- c("Base", "PD1", "RTPD1")
response_order <- c("NR", "R")

cluster_identity_full_map <- c(
  "1" = "P1 | IFN-high cytotoxic",
  "2" = "P2 | Activated CD4/Treg-like",
  "3" = "P3 | APC-like residual",
  "4" = "P4 | Ambiguous residual",
  "5" = "P5 | Cytotoxic CD8 effector-like",
  "6" = "P6 | Activated helper/memory IEG-like",
  "7" = "P7 | gdT/NK-like terminal cytotoxic",
  "8" = "P8 | Helper/memory basal-like"
)

cluster_identity_short_map <- c(
  "1" = "P1 IFN-high\ncytotoxic",
  "2" = "P2 Activated\nCD4/Treg-like",
  "3" = "P3 APC-like\nresidual",
  "4" = "P4 Ambiguous\nresidual",
  "5" = "P5 Cytotoxic CD8\neffector-like",
  "6" = "P6 Helper/memory\nIEG-like",
  "7" = "P7 gdT/NK-like\nterminal cytotoxic",
  "8" = "P8 Helper/memory\nbasal-like"
)

module_label_map <- c(
  "RF240_global" = "RF240 global",
  "IFN_cytotoxic_arm" = "IFN/cytotoxic",
  "Innate_like_arm" = "Innate-like",
  "NR_connected_hubs" = "NR hubs",
  "R_connected_hubs" = "R hubs",
  "Regulatory_adaptive_arm" = "Regulatory/adaptive",
  "Trafficking_RNA_proteostasis_arm" = "Trafficking/RNA/\nproteostasis"
)

expansion_category_order <- c(
  "no_clonotype",
  "singleton",
  "expanded_ge2",
  "expanded_ge5",
  "large_gt10"
)

expansion_category_label_map <- c(
  "no_clonotype" = "No TCR",
  "singleton" = "Sing.",
  "expanded_ge2" = "Exp. ≥2",
  "expanded_ge5" = "Exp. ≥5",
  "large_gt10" = "Large >10"
)

pal_response <- c(
  "NR" = "#D55E00",
  "R" = "#0072B2"
)

pal_public_private <- c(
  "Private" = "#969696",
  "Public" = "#009E73"
)

extract_cluster_num <- function(x) {
  str_extract(as.character(x), "[1-8]")
}

save_plot_multi <- function(plot, filename_base, width, height) {
  
  ggsave(
    file.path(dir_output_plots, paste0(filename_base, ".pdf")),
    plot,
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    file.path(dir_output_plots, paste0(filename_base, ".png")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  
  ggsave(
    file.path(dir_output_plots, paste0(filename_base, ".tiff")),
    plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
}

theme_supp <- function(base_size = 8.5) {
  theme_bw(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 1),
      strip.text = element_text(face = "bold", size = base_size - 1),
      legend.position = "top",
      panel.grid.minor = element_blank()
    )
}

# =========================================================
# 4. Harmonize
# =========================================================

public_pt_zero <- public_pt_zero %>%
  mutate(
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order)
  )

public_state_zero <- public_state_zero %>%
  mutate(
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order),
    cluster_num = extract_cluster_num(cluster_identity_full),
    cluster_identity_short = recode(cluster_num, !!!cluster_identity_short_map),
    cluster_identity_short = factor(
      cluster_identity_short,
      levels = unname(cluster_identity_short_map)
    )
  )

clonal_diversity <- clonal_diversity %>%
  mutate(
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order)
  )

rf240_public_private <- rf240_public_private %>%
  mutate(
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order),
    cluster_num = extract_cluster_num(cluster_projected),
    cluster_identity_short = recode(cluster_num, !!!cluster_identity_short_map),
    cluster_identity_short = factor(
      cluster_identity_short,
      levels = unname(cluster_identity_short_map)
    ),
    public_private_status = factor(public_private_status, levels = c("private", "public")),
    public_private_label = recode(
      as.character(public_private_status),
      "private" = "Private",
      "public" = "Public"
    ),
    public_private_label = factor(public_private_label, levels = c("Private", "Public")),
    module_label = recode(
      as.character(module),
      !!!module_label_map,
      .default = as.character(module)
    ),
    module_label = factor(
      module_label,
      levels = unname(module_label_map)
    )
  )

rf240_expansion_modules <- rf240_expansion_modules %>%
  mutate(
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order),
    cluster_num = extract_cluster_num(cluster_projected),
    cluster_identity_short = recode(cluster_num, !!!cluster_identity_short_map),
    cluster_identity_short = factor(
      cluster_identity_short,
      levels = unname(cluster_identity_short_map)
    ),
    clone_expansion_category = factor(
      clone_expansion_category,
      levels = expansion_category_order
    ),
    clone_expansion_label = recode(
      as.character(clone_expansion_category),
      !!!expansion_category_label_map
    ),
    clone_expansion_label = factor(
      clone_expansion_label,
      levels = unname(expansion_category_label_map)
    ),
    module_label = recode(
      as.character(module),
      !!!module_label_map,
      .default = as.character(module)
    ),
    module_label = factor(
      module_label,
      levels = unname(module_label_map)
    )
  )

# =========================================================
# 5. Panel A
# =========================================================

pA <- ggplot(
  public_pt_zero,
  aes(
    x = timepoint,
    y = prop_clonotypes,
    fill = response_pcr
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    width = 0.62,
    alpha = 0.78
  ) +
  geom_point(
    aes(color = response_pcr),
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.75
    ),
    size = 1.1,
    alpha = 0.75
  ) +
  scale_fill_manual(values = pal_response, drop = FALSE) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  labs(
    title = "A. Public clonotype fraction",
    x = NULL,
    y = "Public clonotype fraction",
    fill = "Response",
    color = "Response"
  ) +
  theme_supp(8.5)

# =========================================================
# 6. Panel B
# =========================================================

pB <- ggplot(
  public_pt_zero,
  aes(
    x = timepoint,
    y = prop_tcr_cells,
    fill = response_pcr
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    width = 0.62,
    alpha = 0.78
  ) +
  geom_point(
    aes(color = response_pcr),
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.75
    ),
    size = 1.1,
    alpha = 0.75
  ) +
  scale_fill_manual(values = pal_response, drop = FALSE) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  labs(
    title = "B. Public TCR-cell fraction",
    x = NULL,
    y = "Public TCR-cell fraction",
    fill = "Response",
    color = "Response"
  ) +
  theme_supp(8.5)

# =========================================================
# 7. Panel C
# =========================================================

# Final audited diversity metrics selected explicitly.
diversity_cols_main <- c(
  "normalized_shannon",
  "shannon",
  "simpson"
)

missing_diversity_cols <- setdiff(
  diversity_cols_main,
  names(clonal_diversity)
)

if (length(missing_diversity_cols) > 0) {
  stop(
    "Missing expected diversity metrics: ",
    paste(
      missing_diversity_cols,
      collapse = ", "
    )
  )
}

diversity_long <- clonal_diversity %>%
  select(patient_id, timepoint, response_pcr, all_of(diversity_cols_main)) %>%
  pivot_longer(
    cols = all_of(diversity_cols_main),
    names_to = "diversity_metric",
    values_to = "diversity_value"
  ) %>%
  mutate(
    diversity_metric = recode(
      diversity_metric,
      "shannon_entropy" = "Shannon",
      "simpson_index" = "Simpson",
      "clonotype_richness" = "Richness",
      .default = diversity_metric
    )
  )

pC <- ggplot(
  diversity_long,
  aes(
    x = timepoint,
    y = diversity_value,
    fill = response_pcr
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    width = 0.62,
    alpha = 0.78
  ) +
  geom_point(
    aes(color = response_pcr),
    position = position_jitterdodge(
      jitter.width = 0.10,
      dodge.width = 0.75
    ),
    size = 1.0,
    alpha = 0.75
  ) +
  facet_wrap(~ diversity_metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal_response, drop = FALSE) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  labs(
    title = "C. Clonal diversity",
    x = NULL,
    y = "Diversity metric",
    fill = "Response",
    color = "Response"
  ) +
  theme_supp(8.2) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

# =========================================================
# 8. Panel D
# =========================================================

panelD_values <- public_state_zero %>%
  group_by(cluster_identity_short, timepoint, response_pcr) %>%
  summarise(
    median_public_cell_fraction = median(prop_cells_cluster, na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  mutate(
    time_response = paste(timepoint, response_pcr, sep = " "),
    time_response = factor(
      time_response,
      levels = as.vector(outer(time_order, response_order, paste, sep = " "))
    )
  )

write_csv(
  panelD_values,
  file.path(dir_output_tables, "SuppFig6D_clean_public_cell_fraction_by_state_heatmap_values.csv")
)

pD <- ggplot(
  panelD_values,
  aes(
    x = time_response,
    y = cluster_identity_short,
    fill = median_public_cell_fraction
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "white",
    high = "#009E73",
    name = "Median\nfraction"
  ) +
  labs(
    title = "D. Public clonotype contribution by state",
    x = NULL,
    y = NULL
  ) +
  theme_supp(8.2) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.8),
    axis.text.y = element_text(size = 6.8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# =========================================================
# 9. Panel E
# =========================================================

rf240_score_col <- "median_module_score"

rf240_pp_global <- rf240_public_private %>%
  filter(module == "RF240_global") %>%
  filter(!is.na(.data[[rf240_score_col]]))

pE <- ggplot(
  rf240_pp_global,
  aes(
    x = public_private_label,
    y = .data[[rf240_score_col]],
    fill = public_private_label
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.60,
    alpha = 0.78
  ) +
  geom_point(
    aes(color = public_private_label),
    position = position_jitter(width = 0.10),
    size = 0.85,
    alpha = 0.60
  ) +
  facet_grid(response_pcr ~ timepoint) +
  scale_fill_manual(values = pal_public_private, drop = FALSE) +
  scale_color_manual(values = pal_public_private, drop = FALSE) +
  labs(
    title = "E. RF240 global score by public/private status",
    x = NULL,
    y = "RF240 global score",
    fill = "Status",
    color = "Status"
  ) +
  theme_supp(8.0) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  )

# =========================================================
# 10. Panel F
# =========================================================

panelF_values <- rf240_pp_global %>%
  group_by(cluster_identity_short, timepoint, response_pcr, public_private_label) %>%
  summarise(
    median_RF240_global = median(.data[[rf240_score_col]], na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  mutate(
    condition = paste(
      recode(as.character(timepoint), "Base" = "B", "PD1" = "PD1", "RTPD1" = "RTP"),
      response_pcr,
      recode(as.character(public_private_label), "Private" = "Priv", "Public" = "Pub"),
      sep = "_"
    ),
    condition = factor(
      condition,
      levels = as.vector(
        expand.grid(
          timepoint = c("B", "PD1", "RTP"),
          response_pcr = response_order,
          status = c("Priv", "Pub")
        ) %>%
          mutate(condition = paste(timepoint, response_pcr, status, sep = "_")) %>%
          pull(condition)
      )
    )
  )

write_csv(
  panelF_values,
  file.path(dir_output_tables, "SuppFig6F_clean_RF240_global_public_private_state_heatmap_values.csv")
)

pF <- ggplot(
  panelF_values,
  aes(
    x = condition,
    y = cluster_identity_short,
    fill = median_RF240_global
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "Median\nRF240"
  ) +
  labs(
    title = "F. State-resolved RF240 global score by public/private status",
    x = NULL,
    y = NULL
  ) +
  theme_supp(8.0) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 6.2),
    axis.text.y = element_text(size = 6.8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# =========================================================
# 11. Panel G
# =========================================================

rf240_pp_modules <- rf240_public_private %>%
  filter(module != "RF240_global") %>%
  filter(!is.na(.data[[rf240_score_col]]))

panelG_values <- rf240_pp_modules %>%
  group_by(module_label, timepoint, response_pcr, public_private_label) %>%
  summarise(
    median_module_score = median(.data[[rf240_score_col]], na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  mutate(
    condition = paste(
      recode(as.character(timepoint), "Base" = "B", "PD1" = "PD1", "RTPD1" = "RTP"),
      response_pcr,
      recode(as.character(public_private_label), "Private" = "Priv", "Public" = "Pub"),
      sep = "_"
    ),
    condition = factor(
      condition,
      levels = as.vector(
        expand.grid(
          timepoint = c("B", "PD1", "RTP"),
          response_pcr = response_order,
          status = c("Priv", "Pub")
        ) %>%
          mutate(condition = paste(timepoint, response_pcr, status, sep = "_")) %>%
          pull(condition)
      )
    )
  )

write_csv(
  panelG_values,
  file.path(dir_output_tables, "SuppFig6G_clean_RF240_modules_public_private_heatmap_values.csv")
)

pG <- ggplot(
  panelG_values,
  aes(
    x = condition,
    y = module_label,
    fill = median_module_score
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Median\nscore"
  ) +
  labs(
    title = "G. RF240 module scores by public/private status",
    x = NULL,
    y = NULL
  ) +
  theme_supp(8.0) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 6.2),
    axis.text.y = element_text(size = 6.8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# =========================================================
# 12. Panel H
# =========================================================

rf240_exp_modules <- rf240_expansion_modules %>%
  filter(module != "RF240_global") %>%
  filter(!is.na(.data[[rf240_score_col]]))

panelH_values <- rf240_exp_modules %>%
  group_by(module_label, timepoint, response_pcr, clone_expansion_label) %>%
  summarise(
    median_module_score = median(.data[[rf240_score_col]], na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  mutate(
    condition = paste(
      recode(as.character(timepoint), "Base" = "B", "PD1" = "PD1", "RTPD1" = "RTP"),
      response_pcr,
      clone_expansion_label,
      sep = "_"
    ),
    condition = factor(
      condition,
      levels = as.vector(
        expand.grid(
          timepoint = c("B", "PD1", "RTP"),
          response_pcr = response_order,
          expansion = unname(expansion_category_label_map)
        ) %>%
          mutate(condition = paste(timepoint, response_pcr, expansion, sep = "_")) %>%
          pull(condition)
      )
    )
  )

write_csv(
  panelH_values,
  file.path(dir_output_tables, "SuppFig6H_clean_RF240_modules_expansion_category_heatmap_values.csv")
)

pH <- ggplot(
  panelH_values,
  aes(
    x = condition,
    y = module_label,
    fill = median_module_score
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    name = "Median\nscore"
  ) +
  labs(
    title = "H. RF240 module scores by clonal expansion category",
    x = NULL,
    y = NULL
  ) +
  theme_supp(7.8) +
  theme(
    axis.text.x = element_text(angle = 60, hjust = 1, size = 5.8),
    axis.text.y = element_text(size = 6.8),
    panel.grid = element_blank(),
    legend.position = "right"
  )

# =========================================================
# 13. Save individual panels
# =========================================================

save_plot_multi(pA, "SuppFig6A_clean_public_clonotype_fraction", 4.6, 3.8)
save_plot_multi(pB, "SuppFig6B_clean_public_TCR_cell_fraction", 4.6, 3.8)
save_plot_multi(pC, "SuppFig6C_clean_clonal_diversity", 7.2, 3.8)
save_plot_multi(pD, "SuppFig6D_clean_public_contribution_by_state", 6.3, 4.9)
save_plot_multi(pE, "SuppFig6E_clean_RF240_global_public_private", 7.0, 4.9)
save_plot_multi(pF, "SuppFig6F_clean_state_RF240_global_public_private", 11.0, 5.1)
save_plot_multi(pG, "SuppFig6G_clean_RF240_modules_public_private", 10.5, 4.2)
save_plot_multi(pH, "SuppFig6H_clean_RF240_modules_expansion_category", 13.0, 4.2)

# =========================================================
# 14. Combined clean figure
# =========================================================

caption_note <- paste0(
  "Note: Public clonotype absence was retained as zero in panels A, B and D. ",
  "Statistical comparisons were performed where group sizes allowed; ",
  "no comparison reached BH-adjusted q < 0.10."
)

combined_clean <- (pA | pB | pC) /
  (pD | pE) /
  pF /
  pG /
  pH +
  plot_layout(
    heights = c(0.78, 0.96, 1.03, 0.80, 0.80)
  ) +
  plot_annotation(
    title = "Supplementary Figure 6. Public/private clonotype structure, clonal diversity, and RF240 module organization",
    subtitle = caption_note
  )

save_plot_multi(
  combined_clean,
  "SuppFig6_combined_public_private_clonotype_RF240_clean_readable",
  width = 16,
  height = 22
)

# =========================================================
# 15. Optional two-page version
# =========================================================
# This is often easier to inspect and can be useful if the final
# journal format allows supplementary figures as multi-page PDFs.

page1 <- (pA | pB | pC) /
  (pD | pE) +
  plot_layout(heights = c(0.85, 1.0)) +
  plot_annotation(
    title = "Supplementary Figure 6. Public/private clonotype structure and clonal diversity",
    subtitle = caption_note
  )

page2 <- pF / pG / pH +
  plot_layout(heights = c(1.05, 0.80, 0.80)) +
  plot_annotation(
    title = "Supplementary Figure 6 continued. RF240 organization by public/private and expansion status",
    subtitle = "Heatmaps show median RF240 global or module scores across longitudinal T-cell contexts."
  )

save_plot_multi(
  page1,
  "SuppFig6_page1_repertoire_clean_readable",
  width = 14,
  height = 9.5
)

save_plot_multi(
  page2,
  "SuppFig6_page2_RF240_clean_readable",
  width = 14,
  height = 13.5
)

# =========================================================
# 16. Summary table
# =========================================================

supp6_clean_summary <- tibble(
  panel = LETTERS[1:8],
  analysis = c(
    "Public clonotype fraction",
    "Public TCR-cell fraction",
    "Clonal diversity",
    "Public clonotype contribution by projected T-cell state",
    "RF240 global score by public/private clonotype status",
    "State-resolved RF240 global score by public/private clonotype status",
    "RF240 module scores by public/private clonotype status",
    "RF240 module scores by clonal expansion category"
  ),
  display_stat_symbols = FALSE,
  reason = "No comparison reached BH-adjusted q < 0.10 in the corrected statistical summary."
)

write_csv(
  supp6_clean_summary,
  file.path(dir_output_tables, "SuppFig6_clean_readable_panel_summary.csv")
)

write_csv(
  stats_summary,
  file.path(dir_output_tables, "SuppFig6_clean_readable_statistics_summary_reference.csv")
)

cat("\n=========================================================\n")
cat("10 DONE - Supplementary Figure 6 repertoire structure - Clean readable Supplementary Figure 6 generated\n")
cat("=========================================================\n")

cat("\n[Plots generated]\n")
print(list.files(dir_output_plots, full.names = FALSE))

cat("\n[Tables generated]\n")
print(list.files(dir_output_tables, full.names = FALSE))

cat("\n[Inspect]\n")
cat("- SuppFig6_combined_public_private_clonotype_RF240_clean_readable.pdf\n")
cat("- SuppFig6_page1_repertoire_clean_readable.pdf\n")
cat("- SuppFig6_page2_RF240_clean_readable.pdf\n")
cat("=========================================================\n")
