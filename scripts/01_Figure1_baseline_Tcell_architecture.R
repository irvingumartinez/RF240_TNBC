# =========================================================
# 01_Figure1_baseline_Tcell_architecture.R
# Purpose:
#   Regenerate Figure 1 panels B-F and Supplementary Figure 1 panels A-D
#   from the curated baseline T-cell reanalysis object.
#
# This script DOES NOT reload h5ad or repeat QC/clustering.
# It uses:
#   - sce_tbase_curated_reanalysis.rds
#   - meta_reanalysis_baseline_Tcells.csv
#   - panelB_RF240_genes_present.csv
#
# Outputs:
#   - PDF, PNG, TIFF panels
#   - support CSV/XLSX tables
#   - manifest
# =========================================================


# =========================================================
# 0. Packages
# =========================================================

required_pkgs <- c(
  "tidyverse",
  "fs",
  "scales",
  "patchwork",
  "writexl",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "scater",
  "Matrix"
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
  library(writexl)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(scater)
  library(Matrix)
})

theme_set(theme_bw(base_size = 11))


# =========================================================
# =========================================================
# 1. Paths
# =========================================================

# Repository root is inferred from this script location:
#   RF240_TNBC/scripts/01_Figure1_baseline_Tcell_architecture.R

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
  "baseline"
)

dir_metadata <- file.path(
  project_root,
  "metadata"
)

dir_outputs <- file.path(
  project_root,
  "results"
)

dir_tables <- file.path(
  dir_outputs,
  "tables",
  "baseline_reanalysis"
)

dir_plots <- file.path(
  dir_outputs,
  "figures",
  "Figure1"
)

dir_rds <- file.path(
  dir_outputs,
  "derived_objects"
)

dir.create(dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_plots, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rds, recursive = TRUE, showWarnings = FALSE)

path_sce <- file.path(
  dir_inputs,
  "sce_tbase_curated_reanalysis.rds"
)

path_meta <- file.path(
  dir_inputs,
  "meta_reanalysis_baseline_Tcells.csv"
)

path_rf240_genes <- file.path(
  dir_metadata,
  "RF240_genes.csv"
)

path_cluster_identity <- file.path(
  dir_metadata,
  "baseline_cluster_identity_map.csv"
)

path_check_tbl <- tibble(
  object = c("path_sce", "path_meta", "path_rf240_genes", "path_cluster_identity"),
  path = c(
    path_sce,
    path_meta,
    path_rf240_genes,
    path_cluster_identity
  ),
  exists = file.exists(
    c(
      path_sce,
      path_meta,
      path_rf240_genes,
      path_cluster_identity
    )
  )
)

print(path_check_tbl)

if (any(!path_check_tbl$exists)) {
  stop("Missing required input files. Check path_check_tbl above.")
}

# 2. Visual standards
# =========================================================

pal_response <- c(
  "R"  = "#0072B2",
  "NR" = "#D55E00"
)

pal_cluster_identity_master <- c(
  "IFN-high cytotoxic" = "#1B9E77",
  "Activated CD4/Treg-like" = "#7570B3",
  "APC-like residual" = "#E7298A",
  "Ambiguous residual" = "#BDBDBD",
  "Cytotoxic CD8 effector-like" = "#D95F02",
  "Activated helper/memory IEG-like" = "#E6AB02",
  "gdT/NK-like terminal cytotoxic" = "#66A61E",
  "Helper/memory basal-like" = "#A6761D"
)

heatmap_div_low  <- "#2166AC"
heatmap_div_mid  <- "white"
heatmap_div_high <- "#B2182B"


# =========================================================
# 3. Helper functions
# =========================================================

save_panel <- function(plot_obj, filename, width, height, out_dir = dir_plots, dpi = 600) {
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".pdf")),
    plot = plot_obj,
    width = width,
    height = height,
    device = "pdf"
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".png")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi
  )
  
  ggsave(
    filename = file.path(out_dir, paste0(filename, ".tiff")),
    plot = plot_obj,
    width = width,
    height = height,
    dpi = dpi,
    compression = "lzw"
  )
}

p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 1e-4 ~ "****",
    p < 1e-3 ~ "***",
    p < 1e-2 ~ "**",
    p < 5e-2 ~ "*",
    TRUE ~ "ns"
  )
}

p_label_or_star <- function(p) {
  dplyr::case_when(
    is.na(p) ~ "",
    p < 5e-2 ~ p_to_stars(p),
    p < 1e-1 ~ paste0("FDR=", signif(p, 2)),
    TRUE ~ ""
  )
}

clean_response <- function(x) {
  x_chr <- as.character(x)
  
  dplyr::case_when(
    x_chr %in% c("R", "r", "Responder", "responder", "1", "pCR", "PCR") ~ "R",
    x_chr %in% c("NR", "nr", "Non-responder", "non-responder", "0", "non-pCR") ~ "NR",
    TRUE ~ NA_character_
  )
}

make_cluster_palette <- function(cluster_identity_levels) {
  pal <- pal_cluster_identity_master[cluster_identity_levels]
  
  missing_levels <- setdiff(
    cluster_identity_levels,
    names(pal_cluster_identity_master)
  )
  
  if (length(missing_levels) > 0) {
    new_cols <- scales::hue_pal()(length(missing_levels))
    names(new_cols) <- missing_levels
    pal <- c(pal, new_cols)
  }
  
  pal[cluster_identity_levels]
}

safe_wilcox <- function(df, value_col, group_col) {
  dsub <- df %>%
    filter(
      !is.na(.data[[value_col]]),
      !is.na(.data[[group_col]])
    )
  
  if (nrow(dsub) == 0) return(NA_real_)
  if (length(unique(dsub[[group_col]])) < 2) return(NA_real_)
  
  grp_sizes <- table(dsub[[group_col]])
  if (any(grp_sizes < 2)) return(NA_real_)
  
  tryCatch(
    wilcox.test(dsub[[value_col]] ~ dsub[[group_col]])$p.value,
    error = function(e) NA_real_
  )
}

mean_expr_by_group <- function(mat, groups) {
  groups <- factor(groups)
  
  out <- vapply(levels(groups), function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx])
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE])
    }
  }, numeric(nrow(mat)))
  
  out <- as.matrix(out)
  rownames(out) <- rownames(mat)
  colnames(out) <- levels(groups)
  out
}

pct_positive_by_group <- function(mat, groups, threshold = 0) {
  groups <- factor(groups)
  
  out <- vapply(levels(groups), function(g) {
    idx <- which(groups == g)
    if (length(idx) == 1) {
      as.numeric(mat[, idx] > threshold) * 100
    } else {
      Matrix::rowMeans(mat[, idx, drop = FALSE] > threshold) * 100
    }
  }, numeric(nrow(mat)))
  
  out <- as.matrix(out)
  rownames(out) <- rownames(mat)
  colnames(out) <- levels(groups)
  out
}


# =========================================================
# 4. Load curated object and metadata
# =========================================================

cat("\n[1] Loading curated SCE and metadata...\n")

sce <- readRDS(path_sce)
meta_csv <- readr::read_csv(path_meta, show_col_types = FALSE)

cat("\nSCE dimensions genes x cells:\n")
print(dim(sce))

cat("\nAssays:\n")
print(assayNames(sce))

cat("\nReducedDims:\n")
print(reducedDimNames(sce))

if (!"logcounts" %in% assayNames(sce)) {
  stop("The curated SCE does not contain a logcounts assay.")
}

if (!"UMAP" %in% reducedDimNames(sce)) {
  stop("The curated SCE does not contain reducedDim 'UMAP'.")
}

meta <- as_tibble(as.data.frame(colData(sce))) %>%
  mutate(
    cell_id = as.character(cell_id),
    patient_id = as.character(patient_id),
    response = factor(clean_response(response_pcr), levels = c("R", "NR")),
    timepoint = as.character(timepoint),
    original_cluster = as.character(subcluster_exact_original),
    cluster_reanalysis = as.character(cluster_reanalysis)
  )

# Harmonize cluster labels as C1-C8
meta <- meta %>%
  mutate(
    re_cluster = if_else(
      stringr::str_detect(cluster_reanalysis, "^C"),
      cluster_reanalysis,
      paste0("C", cluster_reanalysis)
    ),
    re_cluster = factor(re_cluster, levels = paste0("C", sort(unique(as.integer(stringr::str_remove(re_cluster, "^C"))))))
  )

cat("\nResponse distribution:\n")
print(table(meta$response, useNA = "ifany"))

cat("\nPatients by response:\n")
print(meta %>% distinct(patient_id, response) %>% count(response))

cat("\nReanalyzed cluster sizes:\n")
print(meta %>% count(re_cluster, name = "n_cells") %>% mutate(prop = n_cells / sum(n_cells)))

# Ensure order in SCE and meta is identical
if (!identical(meta$cell_id, colnames(sce))) {
  warning("meta$cell_id and colnames(sce) are not identical. Reordering metadata.")
  meta <- meta[match(colnames(sce), meta$cell_id), ]
}

stopifnot(identical(meta$cell_id, colnames(sce)))


# =========================================================
# 5. Cluster identity map
# =========================================================

cat("\n[2] Resolving cluster identity map...\n")

cluster_levels <- levels(meta$re_cluster)

if (file.exists(path_cluster_identity)) {
  cluster_identity_map_tbl <- readr::read_csv(path_cluster_identity, show_col_types = FALSE) %>%
    mutate(
      re_cluster = as.character(re_cluster),
      cluster_identity = as.character(cluster_identity)
    )
} else {
  cluster_identity_map_tbl <- tibble(
    re_cluster = cluster_levels,
    cluster_identity = cluster_levels
  )
  
  message(
    "No cluster_identity_map.csv found. Using C1-C8 labels temporarily.\n",
    "After reviewing marker and RF240 plots, create cluster_identity_map.csv."
  )
}

meta <- meta %>%
  left_join(cluster_identity_map_tbl, by = "re_cluster") %>%
  mutate(
    cluster_identity = if_else(
      is.na(cluster_identity),
      as.character(re_cluster),
      cluster_identity
    ),
    cluster_identity = factor(
      cluster_identity,
      levels = unique(cluster_identity_map_tbl$cluster_identity)
    )
  )

colData(sce)$re_cluster <- meta$re_cluster
colData(sce)$cluster_identity <- meta$cluster_identity

cluster_identity_levels <- levels(meta$cluster_identity)
pal_clusters <- make_cluster_palette(cluster_identity_levels)

readr::write_csv(
  cluster_identity_map_tbl,
  file.path(dir_tables, "Fig1_cluster_identity_map_used.csv")
)


# =========================================================
# 6. Figure 1B and Figure 1D UMAPs
# =========================================================

cat("\n[3] Generating Figure 1B and Figure 1D...\n")

umap_mat <- reducedDim(sce, "UMAP")

umap_df <- tibble(
  UMAP1 = umap_mat[, 1],
  UMAP2 = umap_mat[, 2]
) %>%
  bind_cols(meta)

p_fig1b <- ggplot(umap_df, aes(UMAP1, UMAP2, color = cluster_identity)) +
  geom_point(size = 0.25, alpha = 0.85) +
  scale_color_manual(values = pal_clusters, drop = FALSE) +
  coord_equal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Reanalyzed baseline T-cell clusters",
    color = "Cluster identity",
    x = "UMAP1",
    y = "UMAP2"
  )

p_fig1d <- ggplot(umap_df, aes(UMAP1, UMAP2, color = response)) +
  geom_point(size = 0.25, alpha = 0.85) +
  scale_color_manual(values = pal_response, drop = FALSE, na.value = "grey80") +
  coord_equal() +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(
    title = "Response distribution across the baseline T-cell landscape",
    color = "Response",
    x = "UMAP1",
    y = "UMAP2"
  )

save_panel(p_fig1b, "Fig1B_reanalyzed_clusters_umap", width = 6.4, height = 5.6)
save_panel(p_fig1d, "Fig1D_response_umap", width = 6.2, height = 5.4)


# =========================================================
# 7. Supplementary Figure 1A-B: original vs reanalyzed overlap
# =========================================================

cat("\n[4] Generating Supplementary Figure 1A-B...\n")

overlap_tbl <- meta %>%
  filter(!is.na(original_cluster)) %>%
  count(original_cluster, cluster_identity, name = "n") %>%
  group_by(original_cluster) %>%
  mutate(prop_within_original = n / sum(n)) %>%
  ungroup()

overlap_tbl_rev <- meta %>%
  filter(!is.na(original_cluster)) %>%
  count(cluster_identity, original_cluster, name = "n") %>%
  group_by(cluster_identity) %>%
  mutate(prop_within_reanalyzed = n / sum(n)) %>%
  ungroup()

readr::write_csv(overlap_tbl, file.path(dir_tables, "SuppFig1A_overlap_original_to_reanalyzed.csv"))
readr::write_csv(overlap_tbl_rev, file.path(dir_tables, "SuppFig1B_overlap_reanalyzed_to_original.csv"))

p_supp1a <- ggplot(
  overlap_tbl,
  aes(x = cluster_identity, y = original_cluster, fill = prop_within_original)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "grey20", labels = percent_format(accuracy = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Original cluster redistribution into reanalyzed clusters",
    x = "Reanalyzed cluster identity",
    y = "Original cluster",
    fill = "Within-original\nproportion"
  )

p_supp1b <- ggplot(
  overlap_tbl_rev,
  aes(x = original_cluster, y = cluster_identity, fill = prop_within_reanalyzed)
) +
  geom_tile(color = "white") +
  scale_fill_gradient(low = "white", high = "grey20", labels = percent_format(accuracy = 1)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Reanalyzed cluster composition by original annotation",
    x = "Original cluster",
    y = "Reanalyzed cluster identity",
    fill = "Within-reanalyzed\nproportion"
  )

save_panel(p_supp1a, "SuppFig1A_original_vs_reanalyzed_overlap", width = 6.8, height = 5.6)
save_panel(p_supp1b, "SuppFig1B_original_cluster_redistribution", width = 6.8, height = 5.6)


# =========================================================
# 8. Figure 1C: directed marker dot plot
# =========================================================

cat("\n[5] Generating Figure 1C...\n")

marker_panel_tbl <- tribble(
  ~marker_group,              ~gene,
  "T_lineage",                "CD3D",
  "T_lineage",                "CD3E",
  "T_lineage",                "TRAC",
  "T_lineage",                "TRBC1",
  "T_lineage",                "TRBC2",
  "Helper_memory",            "IL7R",
  "Helper_memory",            "LTB",
  "Helper_memory",            "TCF7",
  "Helper_memory",            "ICOS",
  "Helper_memory",            "MALAT1",
  "Cytotoxic_CD8",            "CD8A",
  "Cytotoxic_CD8",            "CCL5",
  "Cytotoxic_CD8",            "GZMK",
  "Cytotoxic_CD8",            "GZMB",
  "Cytotoxic_CD8",            "PRF1",
  "Cytotoxic_CD8",            "CXCR4",
  "Innate_gdT_NK_like",       "TRDC",
  "Innate_gdT_NK_like",       "NKG7",
  "Innate_gdT_NK_like",       "GNLY",
  "Innate_gdT_NK_like",       "XCL1",
  "Innate_gdT_NK_like",       "XCL2",
  "Innate_gdT_NK_like",       "KLRD1",
  "Innate_gdT_NK_like",       "TYROBP",
  "Innate_gdT_NK_like",       "HOPX",
  "Treg_checkpoint",          "FOXP3",
  "Treg_checkpoint",          "CTLA4",
  "Treg_checkpoint",          "TIGIT",
  "Treg_checkpoint",          "TNFRSF18",
  "Treg_checkpoint",          "TNFRSF4",
  "Treg_checkpoint",          "BATF",
  "APC_like",                 "CD74",
  "APC_like",                 "HLA-DRA",
  "APC_like",                 "HLA-DPA1",
  "APC_like",                 "IFI30",
  "APC_like",                 "APOE",
  "APC_like",                 "CST3",
  "IFN_RF240",                "ISG15",
  "IFN_RF240",                "GBP2",
  "IFN_RF240",                "NR3C1",
  "IFN_RF240",                "CNOT6L",
  "IFN_RF240",                "CD244"
)

marker_genes_present <- intersect(marker_panel_tbl$gene, rownames(sce))

marker_panel_use <- marker_panel_tbl %>%
  filter(gene %in% marker_genes_present)

marker_mat <- logcounts(sce)[marker_genes_present, , drop = FALSE]
cluster_group <- meta$cluster_identity

marker_mean <- mean_expr_by_group(marker_mat, cluster_group)
marker_pct <- pct_positive_by_group(marker_mat, cluster_group)

marker_mean_z <- t(scale(t(marker_mean)))
marker_mean_z[is.na(marker_mean_z)] <- 0

marker_dot_tbl <- as_tibble(marker_mean_z, rownames = "gene") %>%
  pivot_longer(-gene, names_to = "cluster_identity", values_to = "z_mean") %>%
  left_join(
    as_tibble(marker_pct, rownames = "gene") %>%
      pivot_longer(-gene, names_to = "cluster_identity", values_to = "pct_positive"),
    by = c("gene", "cluster_identity")
  ) %>%
  left_join(marker_panel_use, by = "gene") %>%
  mutate(
    cluster_identity = factor(cluster_identity, levels = cluster_identity_levels),
    gene = factor(gene, levels = rev(unique(marker_panel_use$gene))),
    marker_group = factor(marker_group, levels = unique(marker_panel_tbl$marker_group))
  )

readr::write_csv(marker_dot_tbl, file.path(dir_tables, "Fig1C_directed_marker_dotplot_data.csv"))
writexl::write_xlsx(list(marker_dotplot = marker_dot_tbl), file.path(dir_tables, "Fig1C_directed_marker_dotplot_data.xlsx"))

p_fig1c <- ggplot(marker_dot_tbl, aes(x = cluster_identity, y = gene)) +
  geom_point(aes(size = pct_positive, color = z_mean)) +
  scale_size(range = c(0.5, 7), limits = c(0, 100)) +
  scale_color_gradient2(
    low = heatmap_div_low,
    mid = heatmap_div_mid,
    high = heatmap_div_high,
    midpoint = 0
  ) +
  facet_grid(marker_group ~ ., scales = "free_y", space = "free_y") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    strip.background = element_rect(fill = "grey95"),
    strip.text.y = element_text(angle = 0)
  ) +
  labs(
    title = "Directed marker-based annotation of reanalyzed clusters",
    x = "Cluster identity",
    y = NULL,
    color = "Mean expr.\n(z-score)",
    size = "% positive"
  )

save_panel(p_fig1c, "Fig1C_directed_marker_dotplot", width = 8.2, height = 11.0)


# =========================================================
# =========================================================
# 9. Figure 1E: RF240 heatmap
#
# Final audited ordering:
#   1) dominant cluster defined by maximal z-scored mean
#      expression;
#   2) descending dominant-cluster z-score;
#   3) descending dominant-cluster mean expression.
# =========================================================

cat("\n[6] Generating Figure 1E...\n")

rf240_tbl_raw <- readr::read_csv(
  path_rf240_genes,
  show_col_types = FALSE
)

rf240_gene_col <- names(rf240_tbl_raw)[1]

rf240_genes <- unique(
  as.character(
    rf240_tbl_raw[[rf240_gene_col]]
  )
)

rf240_genes <- rf240_genes[
  !is.na(rf240_genes) &
    rf240_genes != ""
]

rf240_genes_present <- intersect(
  rf240_genes,
  rownames(sce)
)

if (length(rf240_genes_present) == 0) {
  stop("No RF240 genes present in the curated SCE.")
}

rf240_mat <- logcounts(sce)[
  rf240_genes_present,
  ,
  drop = FALSE
]

rf240_mean <- mean_expr_by_group(
  rf240_mat,
  cluster_group
)

rf240_pct <- pct_positive_by_group(
  rf240_mat,
  cluster_group
)

rf240_mean_z <- t(
  scale(
    t(rf240_mean)
  )
)

rf240_mean_z[is.na(rf240_mean_z)] <- 0

rf240_stats_tbl <- as_tibble(
  rf240_mean_z,
  rownames = "gene"
) %>%
  pivot_longer(
    -gene,
    names_to = "cluster_identity",
    values_to = "z_mean"
  ) %>%
  left_join(
    as_tibble(
      rf240_mean,
      rownames = "gene"
    ) %>%
      pivot_longer(
        -gene,
        names_to = "cluster_identity",
        values_to = "mean_expr"
      ),
    by = c("gene", "cluster_identity")
  ) %>%
  left_join(
    as_tibble(
      rf240_pct,
      rownames = "gene"
    ) %>%
      pivot_longer(
        -gene,
        names_to = "cluster_identity",
        values_to = "pct_positive"
      ),
    by = c("gene", "cluster_identity")
  ) %>%
  mutate(
    cluster_identity = factor(
      cluster_identity,
      levels = cluster_identity_levels
    )
  )

readr::write_csv(
  rf240_stats_tbl,
  file.path(
    dir_tables,
    "Fig1E_RF240_cluster_stats.csv"
  )
)

# Identify the dominant cluster using the maximal z-score.
rf240_gene_order_tbl <- rf240_stats_tbl %>%
  group_by(gene) %>%
  slice_max(
    order_by = z_mean,
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  rename(
    dominant_cluster = cluster_identity,
    dominant_z_mean = z_mean,
    dominant_mean_expr = mean_expr,
    dominant_pct_positive = pct_positive
  ) %>%
  mutate(
    dominant_cluster = factor(
      dominant_cluster,
      levels = cluster_identity_levels
    )
  ) %>%
  arrange(
    dominant_cluster,
    desc(dominant_z_mean),
    desc(dominant_mean_expr)
  ) %>%
  mutate(
    gene_order = row_number()
  )

rf240_rep_genes <- rf240_gene_order_tbl %>%
  group_by(dominant_cluster) %>%
  slice_max(
    order_by = dominant_z_mean,
    n = 4,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  arrange(
    dominant_cluster,
    desc(dominant_z_mean),
    desc(dominant_mean_expr)
  ) %>%
  pull(gene)

rf240_rep_gene_order <- rf240_gene_order_tbl %>%
  filter(
    gene %in% rf240_rep_genes
  ) %>%
  arrange(
    dominant_cluster,
    desc(dominant_z_mean),
    desc(dominant_mean_expr)
  ) %>%
  pull(gene)

rf240_rep_tbl <- rf240_stats_tbl %>%
  filter(
    gene %in% rf240_rep_genes
  ) %>%
  mutate(
    gene = factor(
      gene,
      levels = rev(rf240_rep_gene_order)
    ),
    cluster_identity = factor(
      cluster_identity,
      levels = cluster_identity_levels
    )
  )

readr::write_csv(
  rf240_gene_order_tbl,
  file.path(
    dir_tables,
    "Fig1E_RF240_representative_ordered_gene_support.csv"
  )
)

readr::write_csv(
  rf240_rep_tbl,
  file.path(
    dir_tables,
    "Fig1E_RF240_representative_ordered_heatmap_data.csv"
  )
)

p_fig1e <- ggplot(
  rf240_rep_tbl,
  aes(
    x = cluster_identity,
    y = gene,
    fill = z_mean
  )
) +
  geom_tile(color = NA) +
  scale_fill_gradient2(
    low = heatmap_div_low,
    mid = heatmap_div_mid,
    high = heatmap_div_high,
    midpoint = 0
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    axis.text.y = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Representative RF240 genes across reanalyzed clusters",
    x = "Cluster identity",
    y = "Representative RF240 genes ordered by dominant cluster",
    fill = "Mean expr.\n(z-score)"
  )

save_panel(
  p_fig1e,
  "Fig1E_RF240_representative_heatmap_ordered",
  width = 8.0,
  height = 7.6
)

# 10. Figure 1F and Supplementary Figure 1C
# =========================================================

cat("\n[7] Generating Figure 1F and Supplementary Figure 1C...\n")

cell_level_comp_tbl <- meta %>%
  count(response, cluster_identity, name = "n_cells") %>%
  group_by(response) %>%
  mutate(prop_cells = n_cells / sum(n_cells)) %>%
  ungroup()

readr::write_csv(cell_level_comp_tbl, file.path(dir_tables, "SuppFig1C_cell_level_cluster_composition_response.csv"))

overall_cluster_prop_tbl <- meta %>%
  count(cluster_identity, name = "n_total_cells") %>%
  mutate(
    overall_prop = n_total_cells / sum(n_total_cells)
  )

cell_bias_tbl <- cell_level_comp_tbl %>%
  select(response, cluster_identity, prop_cells) %>%
  pivot_wider(
    names_from = response,
    values_from = prop_cells,
    values_fill = 0
  ) %>%
  mutate(
    response_bias = R - NR,
    dominant_response = if_else(response_bias >= 0, "R", "NR")
  ) %>%
  left_join(overall_cluster_prop_tbl, by = "cluster_identity") %>%
  arrange(response_bias)

readr::write_csv(
  cell_bias_tbl,
  file.path(dir_tables, "Fig1F_cluster_identity_response_bias.csv")
)

p_fig1f <- ggplot(
  cell_bias_tbl,
  aes(
    x = response_bias,
    y = fct_reorder(cluster_identity, response_bias)
  )
) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_segment(
    aes(
      x = 0,
      xend = response_bias,
      y = cluster_identity,
      yend = cluster_identity
    ),
    linewidth = 0.8,
    color = "grey70"
  ) +
  geom_point(
    aes(size = overall_prop, color = dominant_response),
    alpha = 0.95
  ) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  scale_size_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Cluster identity and descriptive response bias",
    x = "Cell-level response bias (R proportion - NR proportion)",
    y = NULL,
    color = "Bias direction",
    size = "Overall\nproportion"
  )

p_supp1c <- ggplot(
  cell_level_comp_tbl,
  aes(x = response, y = prop_cells, fill = cluster_identity)
) +
  geom_col(width = 0.75, color = "white") +
  scale_fill_manual(values = pal_clusters, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "Cell-level cluster composition by response",
    x = "Response",
    y = "Cell-level proportion",
    fill = "Cluster identity"
  )

save_panel(
  p_fig1f,
  "Fig1F_cluster_identity_response_bias",
  width = 7.6,
  height = 5.2
)

save_panel(
  p_supp1c,
  "SuppFig1C_cell_level_cluster_composition_response",
  width = 6.2,
  height = 5.1
)


# =========================================================
# 11. Supplementary Figure 1D: patient-level proportions with FDR
# =========================================================

cat("\n[8] Generating Supplementary Figure 1D...\n")

patient_response_tbl <- meta %>%
  distinct(patient_id, response)

patient_level_prop_tbl <- meta %>%
  count(patient_id, response, cluster_identity, name = "n_cells") %>%
  group_by(patient_id, response) %>%
  mutate(prop_cells = n_cells / sum(n_cells)) %>%
  ungroup() %>%
  right_join(
    tidyr::expand_grid(
      patient_id = unique(meta$patient_id),
      cluster_identity = cluster_identity_levels
    ),
    by = c("patient_id", "cluster_identity")
  ) %>%
  left_join(
    patient_response_tbl,
    by = "patient_id",
    suffix = c("", ".patient")
  ) %>%
  mutate(
    response = coalesce(response, response.patient),
    n_cells = replace_na(n_cells, 0),
    prop_cells = replace_na(prop_cells, 0),
    cluster_identity = factor(cluster_identity, levels = cluster_identity_levels),
    response = factor(response, levels = c("R", "NR"))
  ) %>%
  select(-response.patient)

patient_stats_tbl <- patient_level_prop_tbl %>%
  group_by(cluster_identity) %>%
  group_modify(~ {
    tibble(
      n_R = sum(.x$response == "R"),
      n_NR = sum(.x$response == "NR"),
      median_R = median(.x$prop_cells[.x$response == "R"], na.rm = TRUE),
      median_NR = median(.x$prop_cells[.x$response == "NR"], na.rm = TRUE),
      p_value = safe_wilcox(.x, value_col = "prop_cells", group_col = "response")
    )
  }) %>%
  ungroup() %>%
  mutate(
    p_adj_fdr = p.adjust(p_value, method = "fdr"),
    sig_label = p_label_or_star(p_adj_fdr)
  )

readr::write_csv(
  patient_level_prop_tbl,
  file.path(dir_tables, "SuppFig1D_patient_level_cluster_proportions.csv")
)

readr::write_csv(
  patient_stats_tbl,
  file.path(dir_tables, "SuppFig1D_patient_level_cluster_proportions_wilcoxon_FDR.csv")
)

writexl::write_xlsx(
  list(
    patient_level_proportions = patient_level_prop_tbl,
    statistics = patient_stats_tbl
  ),
  path = file.path(dir_tables, "SuppFig1D_patient_level_cluster_proportions_support.xlsx")
)

patient_anno_tbl <- patient_level_prop_tbl %>%
  group_by(cluster_identity) %>%
  summarise(
    y_pos = max(prop_cells, na.rm = TRUE) * 1.10 + 1e-6,
    .groups = "drop"
  ) %>%
  left_join(patient_stats_tbl, by = "cluster_identity")

p_supp1d <- ggplot(
  patient_level_prop_tbl,
  aes(x = response, y = prop_cells, fill = response)
) +
  geom_boxplot(
    width = 0.60,
    outlier.shape = NA,
    alpha = 0.25,
    color = "grey25"
  ) +
  geom_jitter(
    aes(color = response),
    width = 0.12,
    size = 1.6,
    alpha = 0.9,
    show.legend = FALSE
  ) +
  facet_wrap(~ cluster_identity, scales = "free_y", ncol = 4) +
  geom_text(
    data = patient_anno_tbl %>% filter(sig_label != ""),
    aes(x = 1.5, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    size = 3.4
  ) +
  scale_fill_manual(values = pal_response, drop = FALSE) +
  scale_color_manual(values = pal_response, drop = FALSE) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme(
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey95")
  ) +
  labs(
    title = "Patient-level cluster proportions by response",
    subtitle = "* FDR < 0.05; exact FDR shown when 0.05 <= FDR < 0.10",
    x = "Response",
    y = "Proportion within patient",
    fill = "Response"
  )

save_panel(
  p_supp1d,
  "SuppFig1D_patient_level_cluster_proportions_with_FDR",
  width = 10.2,
  height = 7.6
)


# =========================================================
# 12. Cluster annotation support table
# =========================================================

cat("\n[9] Exporting cluster annotation support table...\n")

logmat_full <- logcounts(sce)
cluster_vec <- meta$cluster_identity

cluster_mean_expr_full <- mean_expr_by_group(logmat_full, cluster_vec)
cluster_pct_expr_full <- pct_positive_by_group(logmat_full, cluster_vec)

cluster_annotation_support_tbl <- as_tibble(
  cluster_mean_expr_full,
  rownames = "gene"
) %>%
  pivot_longer(
    -gene,
    names_to = "cluster_identity",
    values_to = "mean_expr"
  ) %>%
  left_join(
    as_tibble(cluster_pct_expr_full, rownames = "gene") %>%
      pivot_longer(
        -gene,
        names_to = "cluster_identity",
        values_to = "pct_positive"
      ),
    by = c("gene", "cluster_identity")
  ) %>%
  group_by(cluster_identity) %>%
  arrange(cluster_identity, desc(mean_expr), desc(pct_positive)) %>%
  slice_head(n = 30) %>%
  ungroup()

readr::write_csv(
  cluster_annotation_support_tbl,
  file.path(dir_tables, "Fig1_cluster_annotation_support_top30_by_cluster.csv")
)


# =========================================================
# 13. Save curated copy and manifest
# =========================================================

cat("\n[10] Saving curated object copy and manifest...\n")

saveRDS(
  sce,
  file.path(dir_rds, "sce_tbase_curated_reanalysis_used_for_Figure1.rds")
)

manifest_tbl <- tibble(
  file = c(
    list.files(dir_tables, full.names = TRUE),
    list.files(dir_plots, full.names = TRUE),
    list.files(dir_rds, full.names = TRUE)
  )
) %>%
  mutate(
    relative_file = fs::path_rel(file, start = dir_outputs),
    bytes = file.info(file)$size
  ) %>%
  arrange(relative_file)

readr::write_csv(
  manifest_tbl,
  file.path(dir_outputs, "57A_manifest_outputs_from_curated.csv")
)


# =========================================================
# 14. Copy/paste summary for review
# =========================================================

cat("\n\n=========================================================\n")
cat("57A DONE - Figure 1 and Supplementary Figure 1 regenerated from curated object\n")
cat("Output directory:\n")
cat(dir_outputs, "\n")
cat("=========================================================\n")

cat("\n[Input SCE dimensions]\n")
print(dim(sce))

cat("\n[Assays]\n")
print(assayNames(sce))

cat("\n[ReducedDims]\n")
print(reducedDimNames(sce))

cat("\n[Final cluster sizes]\n")
cluster_size_tbl <- meta %>%
  count(cluster_identity, name = "n_cells") %>%
  mutate(prop = n_cells / sum(n_cells))
print(cluster_size_tbl)

cat("\n[Patients by response]\n")
patient_summary_tbl <- meta %>%
  distinct(patient_id, response) %>%
  count(response, name = "n_patients")
print(patient_summary_tbl)

cat("\n[Cell-level response bias]\n")
print(cell_bias_tbl)

cat("\n[Patient-level cluster statistics]\n")
print(patient_stats_tbl)

cat("\n[Generated plot files]\n")
print(list.files(dir_plots, full.names = FALSE))

cat("\n[Generated table files]\n")
print(list.files(dir_tables, full.names = FALSE))

cat("\n=========================================================\n")
cat("END OF 57A SUMMARY\n")
cat("=========================================================\n")

