# =========================================================
# 05_Figure3_longitudinal_RF240_dynamics.R
# Purpose:
#   Replot Figure 3 from audited Block 3 input tables.
#
# Panels:
#   A. Global RF240 score across timepoints by response
#   B. Compact heatmap of RF240 module scores
#   C. Compact heatmap of RF240 module deltas from baseline
#   D. Global RF240 delta from baseline by response
#
# Inputs:
#   Methodological_audit_and_reproducibility/
#     59A_Figure3_block_inputs/source_tables/
#
# Outputs:
#   Methodological_audit_and_reproducibility/
#     59A_Figure3_reworked_outputs/
#       plots/
#       tables/
#
# Notes:
#   - Significance labels use FDR-adjusted q-values.
#   - q < 0.05 is shown as asterisks.
#   - 0.05 <= q < 0.10 is shown explicitly as q = ...
# =========================================================


# =========================================================
# 0. Packages
# =========================================================

required_pkgs <- c(
  "tidyverse",
  "fs",
  "scales",
  "patchwork"
)

missing_pkgs <- required_pkgs[
  !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_pkgs) > 0) {
  stop(
    "Install missing packages before running:\n",
    paste(missing_pkgs, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(scales)
  library(patchwork)
})


# =========================================================
# =========================================================
# 1. Paths
# =========================================================

# Repository root is inferred from this script location:
#   RF240_TNBC/scripts/05_Figure3_longitudinal_RF240_dynamics.R

script_file <- tryCatch(
  normalizePath(
    sys.frame(1)$ofile,
    winslash = "/",
    mustWork = TRUE
  ),
  error = function(e) NA_character_
)

if (is.na(script_file)) {
  stop(
    "Run this script with source() so the repository root can be detected."
  )
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
  "longitudinal"
)

dir_outputs <- file.path(
  project_root,
  "results"
)

dir_plots <- file.path(
  dir_outputs,
  "figures",
  "Figure3"
)

dir_tables <- file.path(
  dir_outputs,
  "tables",
  "Figure3_longitudinal_RF240"
)

dir.create(dir_outputs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)

path_global_scores <- file.path(
  dir_inputs,
  "RF240_longitudinal_global_score_patient_time.csv"
)

path_module_scores <- file.path(
  dir_inputs,
  "RF240_longitudinal_module_scores_patient_time.csv"
)

path_gene_matrix <- file.path(
  dir_inputs,
  "RF240_longitudinal_patient_time_gene_matrix.csv"
)

path_global_delta <- file.path(
  dir_inputs,
  "RF240_global_score_delta_table.csv"
)

path_module_delta <- file.path(
  dir_inputs,
  "RF240_module_score_delta_long_table.csv"
)

path_module_delta_stats <- file.path(
  dir_inputs,
  "RF240_module_delta_wilcox_R_vs_NR.csv"
)

input_check_tbl <- tibble(
  object = c(
    "global_scores",
    "module_scores",
    "gene_matrix",
    "global_delta",
    "module_delta",
    "module_delta_stats"
  ),
  path = c(
    path_global_scores,
    path_module_scores,
    path_gene_matrix,
    path_global_delta,
    path_module_delta,
    path_module_delta_stats
  ),
  exists = file.exists(
    c(
      path_global_scores,
      path_module_scores,
      path_gene_matrix,
      path_global_delta,
      path_module_delta,
      path_module_delta_stats
    )
  )
)

print(input_check_tbl)

if (any(!input_check_tbl$exists)) {
  stop("Missing required Figure 3 input files. Check input_check_tbl above.")
}

# 2. Visual configuration
# =========================================================

response_levels <- c("NR", "R")
timepoint_levels <- c("Base", "PD1", "RTPD1")
delta_type_levels <- c("PD1 - Base", "RTPD1 - Base")

response_colors <- c(
  "NR" = "#E69F00",
  "R"  = "#0072B2"
)

heatmap_low <- "#2166AC"
heatmap_mid <- "white"
heatmap_high <- "#B2182B"

module_order <- c(
  "IFN/cytotoxic",
  "Regulatory/adaptive",
  "Trafficking/RNA/proteostasis",
  "Innate-like",
  "Responder-connected hubs",
  "Non-responder-connected hubs"
)


# =========================================================
# 3. Helper functions
# =========================================================

standardize_response <- function(x) {
  x_chr <- as.character(x)
  x_low <- stringr::str_to_lower(stringr::str_trim(x_chr))
  
  out <- dplyr::case_when(
    x_low %in% c("r", "resp", "responder", "responders", "pcr", "1") ~ "R",
    x_low %in% c("nr", "nonresponder", "non-responder", "non responder", "non-pcr", "0") ~ "NR",
    TRUE ~ x_chr
  )
  
  factor(out, levels = response_levels)
}

standardize_timepoint <- function(x) {
  x_chr <- as.character(x)
  x_low <- stringr::str_to_lower(stringr::str_trim(x_chr))
  
  out <- dplyr::case_when(
    stringr::str_detect(x_low, "rtpd1|rt_pd1|rt-pd1|rt\\+pd1") ~ "RTPD1",
    stringr::str_detect(x_low, "base|baseline|pre") ~ "Base",
    stringr::str_detect(x_low, "pd1") ~ "PD1",
    TRUE ~ x_chr
  )
  
  factor(out, levels = timepoint_levels)
}

standardize_delta_type <- function(x) {
  x_chr <- as.character(x)
  
  # normalize different dash types
  x_norm <- x_chr %>%
    stringr::str_replace_all("−", "-") %>%
    stringr::str_replace_all("–", "-") %>%
    stringr::str_replace_all("—", "-") %>%
    stringr::str_squish()
  
  x_low <- stringr::str_to_lower(x_norm)
  
  out <- dplyr::case_when(
    stringr::str_detect(x_low, "rtpd1") & stringr::str_detect(x_low, "base") ~ "RTPD1 - Base",
    stringr::str_detect(x_low, "delta_rtpd1_minus_base") ~ "RTPD1 - Base",
    stringr::str_detect(x_low, "pd1") & stringr::str_detect(x_low, "base") & !stringr::str_detect(x_low, "rtpd1") ~ "PD1 - Base",
    stringr::str_detect(x_low, "delta_pd1_minus_base") ~ "PD1 - Base",
    TRUE ~ x_norm
  )
  
  factor(out, levels = delta_type_levels)
}

standardize_module_label <- function(x) {
  x_chr <- as.character(x)
  x_low <- stringr::str_to_lower(stringr::str_trim(x_chr))
  
  dplyr::case_when(
    stringr::str_detect(x_low, "ifn") | stringr::str_detect(x_low, "cytotoxic") ~ "IFN/cytotoxic",
    stringr::str_detect(x_low, "reg") | stringr::str_detect(x_low, "adaptive") ~ "Regulatory/adaptive",
    stringr::str_detect(x_low, "traffic") |
      stringr::str_detect(x_low, "rna") |
      stringr::str_detect(x_low, "proteo") ~ "Trafficking/RNA/proteostasis",
    stringr::str_detect(x_low, "innate") ~ "Innate-like",
    stringr::str_detect(x_low, "nr") |
      stringr::str_detect(x_low, "non.*responder") |
      stringr::str_detect(x_low, "non.*connected") ~ "Non-responder-connected hubs",
    stringr::str_detect(x_low, "r_connected") |
      stringr::str_detect(x_low, "responder") ~ "Responder-connected hubs",
    TRUE ~ x_chr
  )
}

format_q_label <- function(q) {
  dplyr::case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    q < 0.10 ~ paste0("q=", scales::number(q, accuracy = 0.01)),
    TRUE ~ ""
  )
}

safe_wilcox_2group <- function(dat, value_col, group_col = "response") {
  dat2 <- dat %>%
    filter(
      !is.na(.data[[value_col]]),
      !is.na(.data[[group_col]])
    )
  
  group_counts <- dat2 %>%
    count(.data[[group_col]], name = "n")
  
  if (nrow(group_counts) < 2) return(NA_real_)
  if (any(group_counts$n < 2)) return(NA_real_)
  
  tryCatch(
    wilcox.test(
      dat2[[value_col]] ~ dat2[[group_col]],
      exact = FALSE
    )$p.value,
    error = function(e) NA_real_
  )
}

save_panel <- function(plot_obj, filename, width, height, dpi = 600) {
  
  # PDF export
  # Use the standard PDF device for better compatibility across RStudio/macOS setups.
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
  
  # PNG export
  ggsave(
    filename = png_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  
  # TIFF export
  ggsave(
    filename = tiff_file,
    plot = plot_obj,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw"
  )
  
  # Explicit confirmation in console
  message("Saved: ", basename(pdf_file))
  message("Saved: ", basename(png_file))
  message("Saved: ", basename(tiff_file))
}

theme_clean <- function(base_size = 11) {
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
# 4. Load and standardize input tables
# =========================================================

cat("\n[1] Loading input tables...\n")

global_raw <- readr::read_csv(path_global_scores, show_col_types = FALSE)
module_raw <- readr::read_csv(path_module_scores, show_col_types = FALSE)
gene_matrix_raw <- readr::read_csv(path_gene_matrix, show_col_types = FALSE)
gdelta_raw <- readr::read_csv(path_global_delta, show_col_types = FALSE)
mdelta_raw <- readr::read_csv(path_module_delta, show_col_types = FALSE)
mdelta_stats_raw <- readr::read_csv(path_module_delta_stats, show_col_types = FALSE)

global_df <- global_raw %>%
  transmute(
    patient_id = as.character(patient_id),
    timepoint = standardize_timepoint(timepoint),
    response = standardize_response(response_pcr),
    n_cells = as.numeric(n_cells),
    global_score = as.numeric(RF240_global_score)
  ) %>%
  filter(!is.na(patient_id), !is.na(timepoint), !is.na(response), !is.na(global_score))

module_df <- module_raw %>%
  transmute(
    patient_id = as.character(patient_id),
    timepoint = standardize_timepoint(timepoint),
    response = standardize_response(response_pcr),
    n_cells = as.numeric(n_cells),
    module_raw = as.character(module),
    module = standardize_module_label(module),
    n_module_genes_present = n_module_genes_present,
    module_score = as.numeric(module_score)
  ) %>%
  mutate(
    module = factor(module, levels = module_order)
  ) %>%
  filter(!is.na(patient_id), !is.na(timepoint), !is.na(response), !is.na(module), !is.na(module_score))

# Global delta table is wide; convert to long.
gdelta_df <- gdelta_raw %>%
  transmute(
    patient_id = as.character(patient_id),
    response = standardize_response(response_pcr),
    `PD1 - Base` = as.numeric(delta_PD1_minus_Base),
    `RTPD1 - Base` = as.numeric(delta_RTPD1_minus_Base)
  ) %>%
  pivot_longer(
    cols = c(`PD1 - Base`, `RTPD1 - Base`),
    names_to = "delta_type",
    values_to = "delta_RF240"
  ) %>%
  mutate(
    delta_type = factor(delta_type, levels = delta_type_levels)
  ) %>%
  filter(!is.na(patient_id), !is.na(response), !is.na(delta_type), !is.na(delta_RF240))

# Module delta table is already long.
mdelta_df <- mdelta_raw %>%
  transmute(
    patient_id = as.character(patient_id),
    response = standardize_response(response_pcr),
    module_raw = as.character(module),
    module = standardize_module_label(module),
    delta_type = standardize_delta_type(delta_type),
    delta_score = as.numeric(delta_score)
  ) %>%
  mutate(
    module = factor(module, levels = module_order)
  ) %>%
  filter(!is.na(patient_id), !is.na(response), !is.na(module), !is.na(delta_type), !is.na(delta_score))

# Export standardized tables for reproducibility.
readr::write_csv(global_df, file.path(dir_tables, "Fig3_global_scores_standardized.csv"))
readr::write_csv(module_df, file.path(dir_tables, "Fig3_module_scores_standardized.csv"))
readr::write_csv(gdelta_df, file.path(dir_tables, "Fig3_global_delta_standardized.csv"))
readr::write_csv(mdelta_df, file.path(dir_tables, "Fig3_module_delta_standardized.csv"))
readr::write_csv(gene_matrix_raw, file.path(dir_tables, "Fig3_patient_time_gene_matrix_source_copy.csv"))
readr::write_csv(mdelta_stats_raw, file.path(dir_tables, "Fig3_module_delta_wilcox_source_copy.csv"))


# =========================================================
# 5. Audit summaries
# =========================================================

cat("\n[2] Building audit summaries...\n")

global_count_tbl <- global_df %>%
  count(timepoint, response, name = "n_patient_timepoints")

module_count_tbl <- module_df %>%
  count(module, timepoint, response, name = "n_rows")

gdelta_count_tbl <- gdelta_df %>%
  count(delta_type, response, name = "n_paired_deltas")

mdelta_count_tbl <- mdelta_df %>%
  count(module, delta_type, response, name = "n_paired_deltas")

readr::write_csv(global_count_tbl, file.path(dir_tables, "Fig3_global_score_patient_time_counts.csv"))
readr::write_csv(module_count_tbl, file.path(dir_tables, "Fig3_module_score_patient_time_counts.csv"))
readr::write_csv(gdelta_count_tbl, file.path(dir_tables, "Fig3_global_delta_counts.csv"))
readr::write_csv(mdelta_count_tbl, file.path(dir_tables, "Fig3_module_delta_counts.csv"))


# =========================================================
# 6. Statistics
# =========================================================

cat("\n[3] Computing statistics...\n")

# ---------------------------------------------------------
# Helper functions for robust grouped Wilcoxon statistics
# ---------------------------------------------------------

calc_r_vs_nr_stats <- function(df, group_vars, value_col) {
  
  df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~ {
      
      dat <- .x %>%
        filter(
          !is.na(.data[[value_col]]),
          !is.na(response)
        )
      
      n_R <- sum(dat$response == "R")
      n_NR <- sum(dat$response == "NR")
      
      median_R <- median(dat[[value_col]][dat$response == "R"], na.rm = TRUE)
      median_NR <- median(dat[[value_col]][dat$response == "NR"], na.rm = TRUE)
      
      p_value <- NA_real_
      
      if (n_R >= 2 && n_NR >= 2) {
        p_value <- tryCatch(
          wilcox.test(
            dat[[value_col]] ~ dat$response,
            exact = FALSE
          )$p.value,
          error = function(e) NA_real_
        )
      }
      
      tibble(
        n_R = n_R,
        n_NR = n_NR,
        median_R = median_R,
        median_NR = median_NR,
        p_value = p_value
      )
    }) %>%
    ungroup() %>%
    mutate(
      p_adj_fdr = p.adjust(p_value, method = "fdr"),
      sig_label = format_q_label(p_adj_fdr)
    )
}

# ---------------------------------------------------------
# Panel A: global RF240 score, R vs NR within each timepoint
# ---------------------------------------------------------

stats_A <- calc_r_vs_nr_stats(
  df = global_df,
  group_vars = c("timepoint"),
  value_col = "global_score"
)

# ---------------------------------------------------------
# Panel B: module score, R vs NR within each module/timepoint
# ---------------------------------------------------------

stats_B <- calc_r_vs_nr_stats(
  df = module_df,
  group_vars = c("module", "timepoint"),
  value_col = "module_score"
)

# ---------------------------------------------------------
# Panel C: module delta, R vs NR within each module/delta type
# ---------------------------------------------------------

stats_C <- calc_r_vs_nr_stats(
  df = mdelta_df,
  group_vars = c("module", "delta_type"),
  value_col = "delta_score"
)

# ---------------------------------------------------------
# Panel D: global delta, R vs NR within each delta type
# ---------------------------------------------------------

stats_D <- calc_r_vs_nr_stats(
  df = gdelta_df,
  group_vars = c("delta_type"),
  value_col = "delta_RF240"
)

# ---------------------------------------------------------
# Export stats
# ---------------------------------------------------------

readr::write_csv(
  stats_A,
  file.path(dir_tables, "Fig3A_global_score_R_vs_NR_stats.csv")
)

readr::write_csv(
  stats_B,
  file.path(dir_tables, "Fig3B_module_score_R_vs_NR_stats.csv")
)

readr::write_csv(
  stats_C,
  file.path(dir_tables, "Fig3C_module_delta_R_vs_NR_stats.csv")
)

readr::write_csv(
  stats_D,
  file.path(dir_tables, "Fig3D_global_delta_R_vs_NR_stats.csv")
)

cat("\n[Panel A stats]\n")
print(stats_A)

cat("\n[Panel B stats]\n")
print(stats_B %>% arrange(module, timepoint), n = 30)

cat("\n[Panel C stats]\n")
print(stats_C %>% arrange(module, delta_type), n = 30)

cat("\n[Panel D stats]\n")
print(stats_D)

# =========================================================
# 7. Plotting tables
# =========================================================

cat("\n[4] Building plotting tables...\n")

# Panel B heatmap:
# Use module-wise z-scored median scores across time-response cells.
# This reduces saturation caused by module scale differences.
panelB_values_raw <- module_df %>%
  group_by(module, timepoint, response) %>%
  summarise(
    median_score = median(module_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = factor(
      paste(timepoint, response, sep = "\n"),
      levels = c(
        "Base\nNR", "Base\nR",
        "PD1\nNR", "PD1\nR",
        "RTPD1\nNR", "RTPD1\nR"
      )
    )
  )

panelB_values <- panelB_values_raw %>%
  group_by(module) %>%
  mutate(
    median_score_z_by_module = as.numeric(scale(median_score))
  ) %>%
  ungroup()

panelB_annot <- stats_B %>%
  mutate(
    x_label = factor(
      paste(timepoint, "R", sep = "\n"),
      levels = levels(panelB_values$x_label)
    )
  ) %>%
  filter(sig_label != "") %>%
  select(module, timepoint, x_label, p_value, p_adj_fdr, sig_label)

# Panel C heatmap:
# Use raw median deltas centered at zero.
panelC_values <- mdelta_df %>%
  group_by(module, delta_type, response) %>%
  summarise(
    median_delta = median(delta_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    x_label = factor(
      paste(delta_type, response, sep = "\n"),
      levels = c(
        "PD1 - Base\nNR", "PD1 - Base\nR",
        "RTPD1 - Base\nNR", "RTPD1 - Base\nR"
      )
    )
  )

panelC_annot <- stats_C %>%
  mutate(
    x_label = factor(
      paste(delta_type, "R", sep = "\n"),
      levels = levels(panelC_values$x_label)
    )
  ) %>%
  filter(sig_label != "") %>%
  select(module, delta_type, x_label, p_value, p_adj_fdr, sig_label)

readr::write_csv(panelB_values, file.path(dir_tables, "Fig3B_module_score_heatmap_values.csv"))
readr::write_csv(panelB_annot, file.path(dir_tables, "Fig3B_module_score_heatmap_annotations.csv"))
readr::write_csv(panelC_values, file.path(dir_tables, "Fig3C_module_delta_heatmap_values.csv"))
readr::write_csv(panelC_annot, file.path(dir_tables, "Fig3C_module_delta_heatmap_annotations.csv"))


# =========================================================
# 8. Panel A
# =========================================================

cat("\n[5] Plotting Panel A...\n")

annot_A <- global_df %>%
  group_by(timepoint) %>%
  summarise(
    y_pos = max(global_score, na.rm = TRUE) +
      0.08 * diff(range(global_score, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(stats_A, by = "timepoint") %>%
  filter(sig_label != "") %>%
  mutate(x_pos = 1.5)

pA <- ggplot(global_df, aes(x = response, y = global_score, fill = response)) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.08, size = 1.7, alpha = 0.75, color = "black") +
  geom_text(
    data = annot_A,
    aes(x = x_pos, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.6
  ) +
  facet_wrap(~ timepoint, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = response_colors, drop = FALSE) +
  labs(
    title = "A. Global RF240 score across treatment timepoints",
    x = NULL,
    y = "Global RF240 score"
  ) +
  theme_clean(11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

save_panel(
  pA,
  "Fig3A_global_RF240_score_by_time_response",
  width = 8.6,
  height = 3.8
)


# =========================================================
# 9. Panel B
# =========================================================

cat("\n[6] Plotting Panel B...\n")

pB <- ggplot(
  panelB_values,
  aes(x = x_label, y = module, fill = median_score_z_by_module)
) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  geom_text(
    data = panelB_annot,
    aes(x = x_label, y = module, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3
  ) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "Median score\nmodule-scaled z"
  ) +
  labs(
    title = "B. RF240 module scores across time and response",
    x = NULL,
    y = NULL
  ) +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "right"
  )

save_panel(
  pB,
  "Fig3B_compact_heatmap_module_scores",
  width = 8.4,
  height = 4.4
)


# =========================================================
# 10. Panel C
# =========================================================

cat("\n[7] Plotting Panel C...\n")

pC <- ggplot(
  panelC_values,
  aes(x = x_label, y = module, fill = median_delta)
) +
  geom_tile(color = "grey85", linewidth = 0.35) +
  geom_text(
    data = panelC_annot,
    aes(x = x_label, y = module, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3
  ) +
  scale_fill_gradient2(
    low = heatmap_low,
    mid = heatmap_mid,
    high = heatmap_high,
    midpoint = 0,
    name = "Median\ndelta"
  ) +
  labs(
    title = "C. RF240 module deltas relative to baseline",
    x = NULL,
    y = NULL
  ) +
  theme_clean(10) +
  theme(
    axis.text.x = element_text(face = "bold", size = 8),
    axis.text.y = element_text(face = "bold", size = 9),
    legend.position = "right"
  )

save_panel(
  pC,
  "Fig3C_compact_heatmap_module_deltas",
  width = 8.0,
  height = 4.2
)


# =========================================================
# 11. Panel D
# =========================================================

cat("\n[8] Plotting Panel D...\n")

annot_D <- gdelta_df %>%
  group_by(delta_type) %>%
  summarise(
    y_pos = max(delta_RF240, na.rm = TRUE) +
      0.08 * diff(range(delta_RF240, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  left_join(stats_D, by = "delta_type") %>%
  filter(sig_label != "") %>%
  mutate(x_pos = 1.5)

pD <- ggplot(gdelta_df, aes(x = response, y = delta_RF240, fill = response)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_violin(trim = FALSE, alpha = 0.35, color = NA) +
  geom_boxplot(width = 0.18, alpha = 0.85, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.08, size = 1.7, alpha = 0.75, color = "black") +
  geom_text(
    data = annot_D,
    aes(x = x_pos, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    fontface = "bold",
    size = 3.6
  ) +
  facet_wrap(~ delta_type, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = response_colors, drop = FALSE) +
  labs(
    title = "D. Global RF240 longitudinal deltas",
    x = NULL,
    y = "Delta global RF240 score"
  ) +
  theme_clean(11) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold")
  )

save_panel(
  pD,
  "Fig3D_global_RF240_delta_by_response",
  width = 7.8,
  height = 3.8
)


# =========================================================
# 12. Combined figure
# =========================================================

cat("\n[9] Building combined figure...\n")

combined_fig3 <- (pA / (pB | pC) / pD) +
  patchwork::plot_layout(heights = c(1.0, 1.15, 1.0))

save_panel(
  combined_fig3,
  "Fig3_combined_reworked",
  width = 12.0,
  height = 13.5
)


# =========================================================
# 13. Manifest and summary
# =========================================================

cat("\n[10] Building manifest...\n")

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

readr::write_csv(
  manifest_tbl,
  file.path(dir_outputs, "05_Figure3_manifest_outputs.csv")
)

cat("\n\n=========================================================\n")
cat("05 DONE - Figure 3 longitudinal RF240 panels generated from audited derived inputs\n")
cat("Output directory:\n")
cat(dir_outputs, "\n")
cat("=========================================================\n")

cat("\n[Input table dimensions]\n")
cat("global_raw: "); print(dim(global_raw))
cat("module_raw: "); print(dim(module_raw))
cat("gene_matrix_raw: "); print(dim(gene_matrix_raw))
cat("gdelta_raw: "); print(dim(gdelta_raw))
cat("mdelta_raw: "); print(dim(mdelta_raw))

cat("\n[Standardized table dimensions]\n")
cat("global_df: "); print(dim(global_df))
cat("module_df: "); print(dim(module_df))
cat("gdelta_df: "); print(dim(gdelta_df))
cat("mdelta_df: "); print(dim(mdelta_df))

cat("\n[Patient/timepoint counts]\n")
print(global_count_tbl)

cat("\n[Global delta counts]\n")
print(gdelta_count_tbl)

cat("\n[Module delta counts]\n")
print(mdelta_count_tbl, n = 60)

cat("\n[Panel A stats]\n")
print(stats_A)

cat("\n[Panel B stats]\n")
print(stats_B %>% arrange(module, timepoint), n = 30)

cat("\n[Panel C stats]\n")
print(stats_C %>% arrange(module, delta_type), n = 30)

cat("\n[Panel D stats]\n")
print(stats_D)

cat("\n[Generated plot files]\n")
print(list.files(dir_plots, full.names = FALSE))

cat("\n[Generated table files]\n")
print(list.files(dir_tables, full.names = FALSE))

cat("\n=========================================================\n")
cat("END OF PUBLIC FIGURE 3 SUMMARY\n")
cat("=========================================================\n")
