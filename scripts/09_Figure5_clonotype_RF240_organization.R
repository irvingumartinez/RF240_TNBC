# =========================================================
# 09_Figure5_clonotype_RF240_organization.R
# Purpose:
#   Final readability adjustment for Figure 5 v4.
#
# Builds from:
#   61A_4c_Figure5_restore_state_labels.R
#
# Main changes:
#   - Shorter expansion labels.
#   - Shorter heatmap condition labels.
#   - Preserve biological state annotations.
#   - Preserve statistics integrated into biological panels.
#
# No change to:
#   - source data
#   - statistical tests
#   - biological interpretation
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
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
  "figure5"
)

dir_output_plots <- file.path(
  project_root,
  "results",
  "figures",
  "Figure5"
)

dir_output_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Figure5_clonotype_RF240"
)

dir_create(dir_output_plots)
dir_create(dir_output_tables)

# 2. Load source tables
# =========================================================

clono <- read_csv(
  file.path(
    dir_source_tables,
    "clonotype_expansion_summary_patient_time_cluster.csv"
  ),
  show_col_types = FALSE
)

rf240_expanded <- read_csv(
  file.path(
    dir_source_tables,
    "RF240_module_scores_expanded_ge2_patient_time_cluster.csv"
  ),
  show_col_types = FALSE
)

# =========================================================
# 3. Orders, labels and palettes
# =========================================================

time_order <- c("Base", "PD1", "RTPD1")
response_order <- c("NR", "R")

cluster_num_order <- as.character(1:8)

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

cluster_identity_compact_map <- c(
  "1" = "P1 IFN-high\ncytotoxic",
  "2" = "P2 Activated\nCD4/Treg-like",
  "3" = "P3 APC-like\nresidual",
  "4" = "P4 Ambiguous\nresidual",
  "5" = "P5 Cytotoxic CD8\neffector-like",
  "6" = "P6 Helper/memory\nIEG-like",
  "7" = "P7 gdT/NK-like\nterminal cytotoxic",
  "8" = "P8 Helper/memory\nbasal-like"
)

expansion_order <- c("no_clonotype", "singleton", "expanded_ge2")

expansion_label_map <- c(
  "no_clonotype" = "No TCR",
  "singleton" = "Sing.",
  "expanded_ge2" = "Exp. >=2"
)

expansion_pretty_map <- c(
  "No TCR" = "No TCR",
  "Sing." = "Sing.",
  "Exp. >=2" = "Exp. ≥2"
)

# Same temporary palette as v3.
# Exact manuscript palette can be harmonized in final global figure pass.
pal_response <- c(
  "NR" = "#D55E00",
  "R" = "#0072B2"
)

pal_expansion <- c(
  "No TCR" = "#D9D9D9",
  "Sing." = "#969696",
  "Exp. >=2" = "#252525"
)

# =========================================================
# 4. Helper functions
# =========================================================

extract_cluster_num <- function(x) {
  x_chr <- as.character(x)
  str_extract(x_chr, "[1-8]")
}

format_q_compact <- function(q) {
  case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    q < 0.10 ~ ".",
    TRUE ~ ""
  )
}

safe_wilcox_two_groups <- function(df, group_col, value_col) {
  
  group_vec <- df[[group_col]]
  value_vec <- df[[value_col]]
  
  keep <- !is.na(group_vec) & !is.na(value_vec)
  group_vec <- droplevels(as.factor(group_vec[keep]))
  value_vec <- value_vec[keep]
  
  group_counts <- table(group_vec)
  
  if (length(group_counts) != 2 || any(group_counts < 2)) {
    return(
      tibble(
        n_groups = length(group_counts),
        group_1 = ifelse(length(group_counts) >= 1, names(group_counts)[1], NA_character_),
        group_2 = ifelse(length(group_counts) >= 2, names(group_counts)[2], NA_character_),
        n_group_1 = ifelse(length(group_counts) >= 1, as.integer(group_counts[1]), NA_integer_),
        n_group_2 = ifelse(length(group_counts) >= 2, as.integer(group_counts[2]), NA_integer_),
        median_group_1 = NA_real_,
        median_group_2 = NA_real_,
        delta_group_2_minus_group_1 = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  g1 <- names(group_counts)[1]
  g2 <- names(group_counts)[2]
  
  x1 <- value_vec[group_vec == g1]
  x2 <- value_vec[group_vec == g2]
  
  p <- suppressWarnings(
    wilcox.test(x2, x1, exact = FALSE)$p.value
  )
  
  tibble(
    n_groups = 2,
    group_1 = g1,
    group_2 = g2,
    n_group_1 = length(x1),
    n_group_2 = length(x2),
    median_group_1 = median(x1, na.rm = TRUE),
    median_group_2 = median(x2, na.rm = TRUE),
    delta_group_2_minus_group_1 =
      median(x2, na.rm = TRUE) - median(x1, na.rm = TRUE),
    p_value = p
  )
}

safe_kruskal <- function(df, group_col, value_col) {
  
  group_vec <- df[[group_col]]
  value_vec <- df[[value_col]]
  
  keep <- !is.na(group_vec) & !is.na(value_vec)
  group_vec <- droplevels(as.factor(group_vec[keep]))
  value_vec <- value_vec[keep]
  
  group_counts <- table(group_vec)
  
  if (length(group_counts) < 2 || any(group_counts < 2)) {
    return(
      tibble(
        n_groups = length(group_counts),
        n_total = length(value_vec),
        p_value = NA_real_
      )
    )
  }
  
  p <- suppressWarnings(
    kruskal.test(value_vec ~ group_vec)$p.value
  )
  
  tibble(
    n_groups = length(group_counts),
    n_total = length(value_vec),
    p_value = p
  )
}

save_plot_multi <- function(plot, filename_base, width, height) {
  
  ggsave(
    filename = file.path(dir_output_plots, paste0(filename_base, ".pdf")),
    plot = plot,
    width = width,
    height = height,
    units = "in"
  )
  
  ggsave(
    filename = file.path(dir_output_plots, paste0(filename_base, ".png")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300
  )
  
  ggsave(
    filename = file.path(dir_output_plots, paste0(filename_base, ".tiff")),
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    compression = "lzw"
  )
}

# =========================================================
# 5. Harmonize labels
# =========================================================

clono <- clono %>%
  mutate(
    cluster_num = extract_cluster_num(cluster_projected),
    cluster_identity_full = recode(cluster_num, !!!cluster_identity_full_map),
    cluster_identity_compact = recode(cluster_num, !!!cluster_identity_compact_map),
    cluster_identity_full = factor(
      cluster_identity_full,
      levels = unname(cluster_identity_full_map)
    ),
    cluster_identity_compact = factor(
      cluster_identity_compact,
      levels = unname(cluster_identity_compact_map)
    ),
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order)
  )

rf240_expanded <- rf240_expanded %>%
  mutate(
    cluster_num = extract_cluster_num(cluster_projected),
    cluster_identity_full = recode(cluster_num, !!!cluster_identity_full_map),
    cluster_identity_compact = recode(cluster_num, !!!cluster_identity_compact_map),
    cluster_identity_full = factor(
      cluster_identity_full,
      levels = unname(cluster_identity_full_map)
    ),
    cluster_identity_compact = factor(
      cluster_identity_compact,
      levels = unname(cluster_identity_compact_map)
    ),
    timepoint = factor(timepoint, levels = time_order),
    response_pcr = factor(response_pcr, levels = response_order),
    expansion_binary_ge2 = factor(expansion_binary_ge2, levels = expansion_order),
    expansion_label = recode(
      as.character(expansion_binary_ge2),
      !!!expansion_label_map
    ),
    expansion_label = factor(
      expansion_label,
      levels = unname(expansion_label_map)
    ),
    expansion_label_pretty = recode(
      as.character(expansion_label),
      !!!expansion_pretty_map
    ),
    expansion_label_pretty = factor(
      expansion_label_pretty,
      levels = unname(expansion_pretty_map)
    )
  )

cluster_qc <- bind_rows(
  clono %>%
    distinct(source = "clono", cluster_projected, cluster_num, cluster_identity_full),
  rf240_expanded %>%
    distinct(source = "rf240_expanded", cluster_projected, cluster_num, cluster_identity_full)
) %>%
  arrange(source, cluster_num)

write_csv(
  cluster_qc,
  file.path(dir_output_tables, "Fig5_v4_cluster_label_mapping_QC.csv")
)

cat("\n[Cluster label mapping QC]\n")
print(cluster_qc, n = 50)

# =========================================================
# 6. PANEL A
# =========================================================

panelA_value <- "prop_tcr_cells_expanded_ge2"

panelA_stats <- clono %>%
  group_by(cluster_identity_compact, timepoint) %>%
  group_modify(
    ~ safe_wilcox_two_groups(
      df = .x,
      group_col = "response_pcr",
      value_col = panelA_value
    )
  ) %>%
  ungroup() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    stat_label = format_q_compact(q_value)
  )

panelA_ypos <- clono %>%
  group_by(cluster_identity_compact, timepoint) %>%
  summarise(
    y_pos = max(.data[[panelA_value]], na.rm = TRUE) + 0.05,
    .groups = "drop"
  )

panelA_stats_annot <- panelA_stats %>%
  left_join(panelA_ypos, by = c("cluster_identity_compact", "timepoint"))

write_csv(
  panelA_stats,
  file.path(dir_output_tables, "Fig5A_v4_R_vs_NR_expansion_by_state_time_stats.csv")
)

pA <- ggplot(
  clono,
  aes(
    x = timepoint,
    y = .data[[panelA_value]],
    fill = response_pcr
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.78
  ) +
  geom_point(
    aes(color = response_pcr),
    position = position_jitterdodge(
      jitter.width = 0.11,
      dodge.width = 0.75
    ),
    size = 1.15,
    alpha = 0.75
  ) +
  geom_text(
    data = panelA_stats_annot,
    aes(
      x = timepoint,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 3.0,
    vjust = 0
  ) +
  facet_wrap(~ cluster_identity_compact, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = pal_response, drop = FALSE) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  labs(
    title = "A. Clonal expansion across projected T-cell states",
    subtitle = "Expansion = clonotypes represented by ≥2 cells; symbols show BH-adjusted R vs NR comparisons: * q<0.05, . q<0.10",
    x = NULL,
    y = "Proportion of TCR cells in expanded clonotypes",
    fill = "Response",
    color = "Response"
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 7),
    strip.text = element_text(face = "bold", size = 7),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

save_plot_multi(
  pA,
  "Fig5A_v4_clonal_expansion_by_state_time_response",
  width = 11.2,
  height = 7.1
)

# =========================================================
# 7. PANEL B
# =========================================================

panelB_values <- clono %>%
  group_by(cluster_identity_full, timepoint, response_pcr) %>%
  summarise(
    median_prop_expanded_ge2 = median(prop_tcr_cells_expanded_ge2, na.rm = TRUE),
    n_profiles = n(),
    .groups = "drop"
  ) %>%
  mutate(
    time_response = paste(timepoint, response_pcr, sep = "_"),
    time_response = factor(
      time_response,
      levels = as.vector(outer(time_order, response_order, paste, sep = "_"))
    )
  )

panelB_stats_annot <- panelA_stats %>%
  mutate(
    cluster_identity_full = recode(
      as.character(cluster_identity_compact),
      !!!setNames(unname(cluster_identity_full_map), unname(cluster_identity_compact_map))
    ),
    cluster_identity_full = factor(
      cluster_identity_full,
      levels = unname(cluster_identity_full_map)
    ),
    time_response = paste(timepoint, "R", sep = "_"),
    time_response = factor(time_response, levels = levels(panelB_values$time_response)),
    stat_label_heatmap = stat_label
  )

write_csv(
  panelB_values,
  file.path(dir_output_tables, "Fig5B_v4_median_clonal_expansion_heatmap_values.csv")
)

pB <- ggplot(
  panelB_values,
  aes(
    x = time_response,
    y = cluster_identity_full,
    fill = median_prop_expanded_ge2
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(
    data = panelB_stats_annot,
    aes(
      x = time_response,
      y = cluster_identity_full,
      label = stat_label_heatmap
    ),
    inherit.aes = FALSE,
    size = 3.0
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "Median\nexpanded ≥2"
  ) +
  labs(
    title = "B. Median clonal expansion by projected T-cell state",
    subtitle = "Symbols inside R columns show BH-adjusted R vs NR comparisons from panel A",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

save_plot_multi(
  pB,
  "Fig5B_v4_median_clonal_expansion_heatmap_annotated_states",
  width = 8.7,
  height = 5.4
)

# =========================================================
# 8. PANEL C
# =========================================================

rf240_global <- rf240_expanded %>%
  filter(module == "RF240_global") %>%
  filter(!is.na(median_module_score))

panelC_stats <- rf240_global %>%
  group_by(cluster_identity_compact, timepoint) %>%
  group_modify(
    ~ safe_kruskal(
      df = .x,
      group_col = "expansion_label",
      value_col = "median_module_score"
    )
  ) %>%
  ungroup() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    stat_label = format_q_compact(q_value)
  )

panelC_ypos <- rf240_global %>%
  group_by(cluster_identity_compact, timepoint) %>%
  summarise(
    y_pos = max(median_module_score, na.rm = TRUE) + 0.06,
    .groups = "drop"
  )

panelC_stats_annot <- panelC_stats %>%
  left_join(panelC_ypos, by = c("cluster_identity_compact", "timepoint"))

write_csv(
  panelC_stats,
  file.path(dir_output_tables, "Fig5C_v4_RF240_global_expansion_status_by_state_time_kruskal_stats.csv")
)

pC <- ggplot(
  rf240_global,
  aes(
    x = expansion_label,
    y = median_module_score,
    fill = expansion_label
  )
) +
  geom_boxplot(
    outlier.shape = NA,
    width = 0.62,
    alpha = 0.78
  ) +
  geom_point(
    aes(color = response_pcr),
    position = position_jitter(width = 0.11),
    size = 0.9,
    alpha = 0.65
  ) +
  geom_text(
    data = panelC_stats_annot,
    aes(
      x = 2,
      y = y_pos,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 2.8,
    vjust = 0
  ) +
  facet_grid(cluster_identity_compact ~ timepoint, scales = "free_y") +
  scale_fill_manual(
    values = pal_expansion,
    labels = unname(expansion_pretty_map),
    drop = FALSE
  ) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  scale_x_discrete(labels = unname(expansion_pretty_map)) +
  labs(
    title = "C. RF240 global score by clonotype expansion status within projected states",
    subtitle = "Points are colored by response; symbols show BH-adjusted Kruskal-Wallis tests across expansion categories within each state/timepoint",
    x = NULL,
    y = "RF240 global score",
    fill = "Expansion status",
    color = "Response"
  ) +
  theme_bw(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 6.5),
    axis.text.y = element_text(size = 6.5),
    strip.text.y = element_text(face = "bold", size = 6.5),
    strip.text.x = element_text(face = "bold", size = 7),
    plot.title = element_text(face = "bold"),
    legend.position = "top"
  )

save_plot_multi(
  pC,
  "Fig5C_v4_RF240_global_by_expansion_annotated_state_time",
  width = 11.2,
  height = 11.2
)

# =========================================================
# 9. PANEL D
# =========================================================

panelD_values <- rf240_global %>%
  group_by(cluster_identity_full, timepoint, response_pcr, expansion_label) %>%
  summarise(
    median_RF240_global = median(median_module_score, na.rm = TRUE),
    n_profiles = n(),
    median_n_cells = median(n_cells, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    time_abbrev = recode(
      as.character(timepoint),
      "Base" = "B",
      "PD1" = "PD1",
      "RTPD1" = "RTPD1"
    ),
    expansion_abbrev = recode(
      as.character(expansion_label),
      "No TCR" = "No",
      "Sing." = "S",
      "Exp. >=2" = "E"
    ),
    condition_short = paste(time_abbrev, response_pcr, expansion_abbrev, sep = "-"),
    condition_short = factor(
      condition_short,
      levels = as.vector(
        expand.grid(
          time_abbrev = c("B", "PD1", "RTPD1"),
          response_pcr = response_order,
          expansion_abbrev = c("No", "S", "E")
        ) %>%
          mutate(condition_short = paste(time_abbrev, response_pcr, expansion_abbrev, sep = "-")) %>%
          pull(condition_short)
      )
    )
  )

panelD_stats <- rf240_global %>%
  filter(expansion_binary_ge2 %in% c("singleton", "expanded_ge2")) %>%
  group_by(cluster_identity_full, timepoint, response_pcr) %>%
  group_modify(
    ~ safe_wilcox_two_groups(
      df = .x,
      group_col = "expansion_label",
      value_col = "median_module_score"
    )
  ) %>%
  ungroup() %>%
  mutate(
    q_value = p.adjust(p_value, method = "BH"),
    stat_label = format_q_compact(q_value),
    time_abbrev = recode(
      as.character(timepoint),
      "Base" = "B",
      "PD1" = "PD1",
      "RTPD1" = "RTPD1"
    ),
    expansion_abbrev = "E",
    condition_short = paste(time_abbrev, response_pcr, expansion_abbrev, sep = "-"),
    condition_short = factor(condition_short, levels = levels(panelD_values$condition_short))
  )

write_csv(
  panelD_values,
  file.path(dir_output_tables, "Fig5D_v4_RF240_global_state_time_response_expansion_heatmap_values.csv")
)

write_csv(
  panelD_stats,
  file.path(dir_output_tables, "Fig5D_v4_RF240_global_expanded_vs_singleton_by_state_time_response_stats.csv")
)

pD <- ggplot(
  panelD_values,
  aes(
    x = condition_short,
    y = cluster_identity_full,
    fill = median_RF240_global
  )
) +
  geom_tile(color = "white", linewidth = 0.25) +
  geom_text(
    data = panelD_stats,
    aes(
      x = condition_short,
      y = cluster_identity_full,
      label = stat_label
    ),
    inherit.aes = FALSE,
    size = 2.7
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "Median\nRF240 global"
  ) +
  labs(
    title = "D. RF240 global score by expansion status, state and timepoint",
    subtitle = "No = no TCR, S = singleton, E = expanded ≥2; symbols inside E columns show BH-adjusted expanded vs singleton comparisons",
    x = NULL,
    y = NULL
  ) +
  theme_bw(base_size = 8.5) +
  theme(
    axis.text.x = element_text(angle = 55, hjust = 1, size = 6.4),
    axis.text.y = element_text(size = 7),
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )

save_plot_multi(
  pD,
  "Fig5D_v4_RF240_global_heatmap_annotated_state_time_response_expansion",
  width = 12.5,
  height = 6.0
)

# =========================================================
# 10. Combined Figure 5 candidate v4
# =========================================================

if (requireNamespace("patchwork", quietly = TRUE)) {
  
  combined_fig5_v4 <- (pA | pB) / pC / pD +
    patchwork::plot_layout(
      heights = c(0.9, 1.45, 0.85)
    ) +
    patchwork::plot_annotation(
      title = "Figure 5. State-resolved clonal expansion and RF240 organization",
      subtitle = "Reworked candidate with biological state annotations, shorter expansion labels, and integrated statistics"
    )
  
  save_plot_multi(
    combined_fig5_v4,
    "Fig5_combined_state_resolved_candidate_v4_final_readability",
    width = 15,
    height = 20
  )
  
} else {
  
  cat("\n[Note]\n")
  cat("Package 'patchwork' is not installed. Individual panels were saved, but combined figure was not generated.\n")
}

# =========================================================
# 11. Summary
# =========================================================

fig5_v4_summary <- tibble(
  panel = c("Figure 5A", "Figure 5B", "Figure 5C", "Figure 5D"),
  key_analysis = c(
    "Clonal expansion by annotated projected state, timepoint and response",
    "Median clonal expansion heatmap by annotated projected state, timepoint and response",
    "RF240_global by expansion status within annotated projected states and timepoints",
    "Median RF240_global heatmap by annotated projected state, timepoint, response and expansion status"
  ),
  statistic_integrated = c(
    "R vs NR Wilcoxon within state/timepoint; BH across panel",
    "Same R vs NR statistic from panel A shown in R columns",
    "Kruskal-Wallis across expansion categories within state/timepoint; BH across panel",
    "Expanded >=2 vs singleton Wilcoxon within state/timepoint/response; BH across panel"
  ),
  readability_change = c(
    "Compact annotated facet labels retained",
    "Full annotated state labels retained",
    "Expansion labels shortened",
    "Condition labels shortened"
  )
)

write_csv(
  fig5_v4_summary,
  file.path(dir_output_tables, "Block5_Figure5_v4_final_readability_summary.csv")
)

cat("\n=========================================================\n")
cat("09 DONE - Figure 5 clonotype RF240 organization\n")
cat("=========================================================\n")

cat("\n[Plots generated]\n")
print(list.files(dir_output_plots, full.names = FALSE))

cat("\n[Tables generated]\n")
print(list.files(dir_output_tables, full.names = FALSE))

cat("\n[Inspect first]\n")
cat("- Fig5_combined_state_resolved_candidate_v4_final_readability.pdf\n")
cat("- Fig5A_v4_clonal_expansion_by_state_time_response.pdf\n")
cat("- Fig5B_v4_median_clonal_expansion_heatmap_annotated_states.pdf\n")
cat("- Fig5C_v4_RF240_global_by_expansion_annotated_state_time.pdf\n")
cat("- Fig5D_v4_RF240_global_heatmap_annotated_state_time_response_expansion.pdf\n")
cat("=========================================================\n")
