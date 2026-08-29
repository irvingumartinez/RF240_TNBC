# =========================================================
# 02_SuppFigure2_RF240_spatial_support.R
# Purpose:
#   Regenerate Supplementary Figure 2 from the curated baseline
#   T-cell reanalysis object.
#
# Panels:
#   SuppFig2A NR3C1 UMAP
#   SuppFig2B CNOT6L UMAP
#   SuppFig2C IFN score UMAP
#   SuppFig2D Regulatory/checkpoint score UMAP
#   SuppFig2E Full RF240 heatmap across reanalyzed clusters
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
  "writexl",
  "SingleCellExperiment",
  "SummarizedExperiment",
  "Matrix",
  "viridis"
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
  library(writexl)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(Matrix)
  library(viridis)
})

theme_set(theme_bw(base_size = 11))


# =========================================================
# =========================================================
# 1. Paths
# =========================================================

# Repository root is inferred from this script location:
#   RF240_TNBC/scripts/02_SuppFigure2_RF240_spatial_support.R

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
  "supplementary_figures",
  "Supplementary_Figure2"
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
  object = c(
    "path_sce",
    "path_meta",
    "path_rf240_genes",
    "path_cluster_identity"
  ),
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

calculate_simple_module_score <- function(sce, genes, assay_name = "logcounts") {
  genes_present <- intersect(genes, rownames(sce))
  
  if (length(genes_present) == 0) {
    warning("No module genes found in object.")
    return(rep(NA_real_, ncol(sce)))
  }
  
  score <- Matrix::colMeans(
    assay(sce, assay_name)[genes_present, , drop = FALSE]
  )
  
  score
}

make_umap_feature_plot <- function(plot_df, feature_col, title, legend_title) {
  ggplot(plot_df, aes(x = UMAP1, y = UMAP2, color = .data[[feature_col]])) +
    geom_point(size = 0.25, alpha = 0.85) +
    scale_color_viridis_c(
      option = "magma",
      direction = 1,
      na.value = "grey90"
    ) +
    coord_equal() +
    theme(
      panel.grid = element_blank(),
      legend.position = "right"
    ) +
    labs(
      title = title,
      x = "UMAP1",
      y = "UMAP2",
      color = legend_title
    )
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

meta <- meta %>%
  mutate(
    re_cluster = if_else(
      stringr::str_detect(cluster_reanalysis, "^C"),
      cluster_reanalysis,
      paste0("C", cluster_reanalysis)
    ),
    re_cluster = factor(
      re_cluster,
      levels = paste0(
        "C",
        sort(unique(as.integer(stringr::str_remove(re_cluster, "^C"))))
      )
    )
  )

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
  cluster_identity_map_tbl <- readr::read_csv(
    path_cluster_identity,
    show_col_types = FALSE
  ) %>%
    mutate(
      re_cluster = as.character(re_cluster),
      cluster_identity = as.character(cluster_identity)
    )
} else {
  cluster_identity_map_tbl <- tibble(
    re_cluster = cluster_levels,
    cluster_identity = cluster_levels
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
  file.path(dir_tables, "SuppFig2_cluster_identity_map_used.csv")
)


# =========================================================
# 6. Define gene sets for spatial support
# =========================================================

cat("\n[3] Defining genes and module scores...\n")

feature_genes <- c("NR3C1", "CNOT6L")

ifn_score_genes <- c(
  "ISG15",
  "GBP2",
  "IFI6",
  "IFITM1",
  "STAT1",
  "CXCL10"
)

regulatory_score_genes <- c(
  "FOXP3",
  "CTLA4",
  "TIGIT",
  "TNFRSF18",
  "TNFRSF4",
  "BATF"
)

feature_gene_status_tbl <- tibble(
  feature = c(
    feature_genes,
    ifn_score_genes,
    regulatory_score_genes
  ),
  present = feature %in% rownames(sce),
  feature_type = case_when(
    feature %in% feature_genes ~ "single_gene",
    feature %in% ifn_score_genes ~ "IFN_score",
    feature %in% regulatory_score_genes ~ "regulatory_score",
    TRUE ~ "other"
  )
)

readr::write_csv(
  feature_gene_status_tbl,
  file.path(dir_tables, "SuppFig2_feature_gene_status.csv")
)

print(feature_gene_status_tbl)


# =========================================================
# 7. Build UMAP feature dataframe
# =========================================================

cat("\n[4] Building UMAP feature dataframe...\n")

umap_mat <- reducedDim(sce, "UMAP")

plot_df <- tibble(
  UMAP1 = umap_mat[, 1],
  UMAP2 = umap_mat[, 2],
  cell_id = colnames(sce),
  cluster_identity = meta$cluster_identity,
  response = meta$response,
  patient_id = meta$patient_id
)

for (gene in feature_genes) {
  if (gene %in% rownames(sce)) {
    plot_df[[gene]] <- as.numeric(logcounts(sce)[gene, ])
  } else {
    plot_df[[gene]] <- NA_real_
  }
}

plot_df$IFN_score <- calculate_simple_module_score(
  sce,
  ifn_score_genes,
  assay_name = "logcounts"
)

plot_df$Regulatory_checkpoint_score <- calculate_simple_module_score(
  sce,
  regulatory_score_genes,
  assay_name = "logcounts"
)

readr::write_csv(
  plot_df %>%
    select(
      cell_id,
      cluster_identity,
      response,
      patient_id,
      NR3C1,
      CNOT6L,
      IFN_score,
      Regulatory_checkpoint_score
    ),
  file.path(dir_tables, "SuppFig2_cell_level_feature_scores.csv")
)


# =========================================================
# 8. Generate Supplementary Figure 2A-D UMAP feature plots
# =========================================================

cat("\n[5] Generating Supplementary Figure 2A-D...\n")

p_supp2a <- make_umap_feature_plot(
  plot_df,
  feature_col = "NR3C1",
  title = "NR3C1 expression",
  legend_title = "log-normalized\nexpression"
)

p_supp2b <- make_umap_feature_plot(
  plot_df,
  feature_col = "CNOT6L",
  title = "CNOT6L expression",
  legend_title = "log-normalized\nexpression"
)

p_supp2c <- make_umap_feature_plot(
  plot_df,
  feature_col = "IFN_score",
  title = "IFN module score",
  legend_title = "mean log-normalized\nexpression"
)

p_supp2d <- make_umap_feature_plot(
  plot_df,
  feature_col = "Regulatory_checkpoint_score",
  title = "Regulatory/checkpoint module score",
  legend_title = "mean log-normalized\nexpression"
)

save_panel(
  p_supp2a,
  "SuppFig2A_NR3C1_umap",
  width = 6.0,
  height = 5.2
)

save_panel(
  p_supp2b,
  "SuppFig2B_CNOT6L_umap",
  width = 6.0,
  height = 5.2
)

save_panel(
  p_supp2c,
  "SuppFig2C_IFN_score_umap",
  width = 6.0,
  height = 5.2
)

save_panel(
  p_supp2d,
  "SuppFig2D_regulatory_checkpoint_score_umap",
  width = 6.0,
  height = 5.2
)


# =========================================================
# =========================================================
# 9. Generate Supplementary Figure 2E full RF240 heatmap
#
# Final audited ordering:
#   1) dominant cluster defined by maximal z-scored mean
#      expression;
#   2) descending dominant-cluster z-score;
#   3) descending dominant-cluster mean expression.
# =========================================================

cat("\n[6] Generating Supplementary Figure 2E full RF240 heatmap...\n")

rf240_tbl_raw <- readr::read_csv(
  path_rf240_genes,
  show_col_types = FALSE
)

rf240_gene_col <- names(
  rf240_tbl_raw
)[1]

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

rf240_genes_missing <- setdiff(
  rf240_genes,
  rownames(sce)
)

rf240_presence_tbl <- tibble(
  gene = rf240_genes,
  present = gene %in% rownames(sce)
)

readr::write_csv(
  rf240_presence_tbl,
  file.path(
    dir_tables,
    "SuppFig2E_RF240_gene_presence.csv"
  )
)

if (length(rf240_genes_present) == 0) {
  stop("No RF240 genes present in the curated SCE.")
}

rf240_mat <- logcounts(sce)[
  rf240_genes_present,
  ,
  drop = FALSE
]

cluster_group <- meta$cluster_identity

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

rf240_mean_z[
  is.na(rf240_mean_z)
] <- 0

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

# Identify dominant cluster by maximal z-scored mean expression.
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

ordered_genes_top_to_bottom <- rf240_gene_order_tbl$gene

rf240_full_heatmap_tbl <- rf240_stats_tbl %>%
  left_join(
    rf240_gene_order_tbl %>%
      select(
        gene,
        dominant_cluster,
        dominant_z_mean,
        dominant_mean_expr,
        dominant_pct_positive,
        gene_order
      ),
    by = "gene"
  ) %>%
  mutate(
    cluster_identity = factor(
      cluster_identity,
      levels = cluster_identity_levels
    ),
    dominant_cluster = factor(
      dominant_cluster,
      levels = cluster_identity_levels
    ),
    gene = factor(
      gene,
      levels = rev(
        ordered_genes_top_to_bottom
      )
    )
  ) %>%
  arrange(
    gene_order,
    cluster_identity
  )

readr::write_csv(
  rf240_stats_tbl,
  file.path(
    dir_tables,
    "SuppFig2E_RF240_full_cluster_stats.csv"
  )
)

readr::write_csv(
  rf240_gene_order_tbl,
  file.path(
    dir_tables,
    "SuppFig2E_RF240_full_gene_order_ordered.csv"
  )
)

readr::write_csv(
  rf240_full_heatmap_tbl,
  file.path(
    dir_tables,
    "SuppFig2E_RF240_full_ordered_heatmap_data.csv"
  )
)

writexl::write_xlsx(
  list(
    RF240_gene_presence = rf240_presence_tbl,
    RF240_cluster_stats = rf240_stats_tbl,
    RF240_gene_order = rf240_gene_order_tbl
  ),
  path = file.path(
    dir_tables,
    "SuppFig2E_RF240_full_heatmap_support.xlsx"
  )
)

p_supp2e <- ggplot(
  rf240_full_heatmap_tbl,
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
    axis.text.y = element_text(size = 3),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Full RF240 gene distribution across reanalyzed clusters",
    x = "Cluster identity",
    y = "RF240 genes ordered by dominant cluster and expression",
    fill = "Mean expr.\n(z-score)"
  )

save_panel(
  p_supp2e,
  "SuppFig2E_full_RF240_heatmap_ordered",
  width = 8.8,
  height = 18.0
)

# 10. Manifest and summary
# =========================================================

cat("\n[7] Building manifest and summary...\n")

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
  file.path(dir_outputs, "57B_manifest_outputs.csv")
)

cat("\n\n=========================================================\n")
cat("57B DONE - Supplementary Figure 2 regenerated from curated object\n")
cat("Output directory:\n")
cat(dir_outputs, "\n")
cat("=========================================================\n")

cat("\n[Input SCE dimensions]\n")
print(dim(sce))

cat("\n[Feature gene status]\n")
print(feature_gene_status_tbl)

cat("\n[RF240 gene presence]\n")
cat("RF240 genes input: ", length(rf240_genes), "\n")
cat("RF240 genes present: ", length(rf240_genes_present), "\n")
cat("RF240 genes missing: ", length(rf240_genes_missing), "\n")

cat("\n[Generated plot files]\n")
print(list.files(dir_plots, full.names = FALSE))

cat("\n[Generated table files]\n")
print(list.files(dir_tables, full.names = FALSE))

cat("\n=========================================================\n")
cat("END OF 57B SUMMARY\n")
cat("=========================================================\n")
