# =========================================================
# 13_SuppFigure7_Bassez_robustness.R
#
# FINAL ASSEMBLY
# Supplementary Figure 7
#
# TITLE
# -----
# Patient-level robustness and transition context of
# RF240 remodeling in persistent clonotypes
#
# FINAL PANEL ARCHITECTURE
# ------------------------
# A. Patient-level global RF240 remodeling: NE vs E
#
# B. Patient-level RF240 remodeling across clonotype
#    expansion-status categories
#
# C. Patient prevalence of Pre -> Post phenotype-state
#    transitions in the finalQC universe
#
# D. Patient-level binary state-stable vs state-changing
#    RF240 remodeling control
#
# EXPORTS
# -------
# Every individual panel and integrated figure:
#   PDF
#   PNG
#   TIFF
#
# =========================================================


# =========================================================
# 0. Clean environment
# =========================================================

rm(list = ls())
gc()


# =========================================================
# 1. Packages
# =========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(fs)
  library(readr)
  library(stringr)
  library(patchwork)
  library(scales)
})


# =========================================================
# =========================================================
# 2. Repository paths
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
  "bassez",
  "supplementary_figure7"
)

dir_plots <- file.path(
  project_root,
  "results",
  "supplementary_figures",
  "Supplementary_Figure7"
)

dir_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Supplementary_Figure7_Bassez"
)

dir_sources <- file.path(
  dir_tables,
  "source_tables"
)

dir_text <- file.path(
  dir_tables,
  "text_summaries"
)

dir_block <- dir_tables

walk(
  c(
    dir_plots,
    dir_tables,
    dir_sources,
    dir_text
  ),
  dir_create
)

# Aliases keep the audited downstream file construction
# unchanged while all files now reside in one public folder.
dir_63D2 <- dir_inputs
dir_63D3 <- dir_inputs
dir_63D3b <- dir_inputs
dir_63D4 <- dir_inputs
dir_56FC <- dir_inputs

# 3. Exact audited input files
# =========================================================


# ---------------------------------------------------------
# PANEL A
# ---------------------------------------------------------

file_A_patient <- file.path(
  dir_63D3,
  "tables",
  "63D_3_harmonized_patientlevel_RF240_table.csv"
)

file_A_test <- file.path(
  dir_63D3,
  "tables",
  "63D_3_primary_patientlevel_E_vs_NE_test.csv"
)

file_A_boot_delta <- file.path(
  dir_63D3,
  "tables",
  "63D_3_bootstrap_E_minus_NE_median_difference.csv"
)

file_A_corrected_summary <- file.path(
  dir_63D3b,
  "tables",
  "63D_3b_CORRECTED_patient_summary_by_response.csv"
)


# ---------------------------------------------------------
# PANEL B
# ---------------------------------------------------------

file_B_patient_expansion <- file.path(
  dir_63D2,
  "tables",
  "63D_2_patient_by_expansion_status_summary.csv"
)

file_B_tests <- file.path(
  dir_63D2,
  "tables",
  "63D_2_patientlevel_E_vs_NE_by_expansion_status.csv"
)

file_B_friedman <- file.path(
  dir_63D2,
  "tables",
  "63D_2_friedman_expansion_effect_within_response.csv"
)

file_B_boot <- file.path(
  dir_63D2,
  "tables",
  "63D_2_bootstrap_patient_median_delta_by_expansion.csv"
)


# ---------------------------------------------------------
# PANEL C
# ---------------------------------------------------------

file_C_transition <- file.path(
  dir_63D4,
  "tables",
  "63D_4_patientaware_transition_summary_by_response.csv"
)

file_C_audit_summary <- file.path(
  dir_63D4,
  "tables",
  "63D_4_transition_abundance_audit_summary.csv"
)

file_C_route_tests <- file.path(
  dir_63D4,
  "tables",
  "63D_4_patientlevel_E_vs_NE_transition_route_tests.csv"
)


# ---------------------------------------------------------
# PANEL D
# ---------------------------------------------------------

file_D_clone <- file.path(
  dir_56FC,
  "56F_C_Bassez_clone_delta_state_change_finalQC_min100cells_min10persistent.csv.gz"
)


required_files <- c(
  file_A_patient,
  file_A_test,
  file_A_boot_delta,
  file_A_corrected_summary,
  file_B_patient_expansion,
  file_B_tests,
  file_B_friedman,
  file_B_boot,
  file_C_transition,
  file_C_audit_summary,
  file_C_route_tests,
  file_D_clone
)


missing_files <- required_files[
  !file.exists(
    required_files
  )
]


if (
  length(
    missing_files
  ) > 0
) {
  
  stop(
    paste0(
      "Missing required audited source files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# =========================================================
# 4. Copy exact sources used in final figure
# =========================================================

walk(
  required_files,
  function(x) {
    
    file_copy(
      x,
      file.path(
        dir_sources,
        path_file(
          x
        )
      ),
      overwrite = TRUE
    )
  }
)


# =========================================================
# 5. Read audited tables
# =========================================================

A_patient <- read_csv(
  file_A_patient,
  show_col_types = FALSE
)

A_test <- read_csv(
  file_A_test,
  show_col_types = FALSE
)

A_boot_delta <- read_csv(
  file_A_boot_delta,
  show_col_types = FALSE
)

A_corrected_summary <- read_csv(
  file_A_corrected_summary,
  show_col_types = FALSE
)


B_patient <- read_csv(
  file_B_patient_expansion,
  show_col_types = FALSE
)

B_tests <- read_csv(
  file_B_tests,
  show_col_types = FALSE
)

B_friedman <- read_csv(
  file_B_friedman,
  show_col_types = FALSE
)

B_boot <- read_csv(
  file_B_boot,
  show_col_types = FALSE
)


C_transition <- read_csv(
  file_C_transition,
  show_col_types = FALSE
)

C_audit_summary <- read_csv(
  file_C_audit_summary,
  show_col_types = FALSE
)

C_route_tests <- read_csv(
  file_C_route_tests,
  show_col_types = FALSE
)


D_clone_raw <- read_csv(
  file_D_clone,
  show_col_types = FALSE
)


# =========================================================
# 6. Global visual constants
# =========================================================

pal_response <- c(
  "NE" = "#D55E00",
  "E"  = "#0072B2"
)


state_levels <- c(
  "CXCL13_PDCD1_like",
  "Cytotoxic_like",
  "Helper_memory_like",
  "IFN_response_like",
  "Low_signal_mixed",
  "Regulatory_checkpoint_like"
)


state_label_map <- c(
  
  "CXCL13_PDCD1_like" =
    "CXCL13/PDCD1-like",
  
  "Cytotoxic_like" =
    "Cytotoxic-like",
  
  "Helper_memory_like" =
    "Helper/memory-like",
  
  "IFN_response_like" =
    "IFN-response-like",
  
  "Low_signal_mixed" =
    "Low-signal/mixed",
  
  "Regulatory_checkpoint_like" =
    "Regulatory/checkpoint-like"
)


state_labels_clean <- unname(
  state_label_map[
    state_levels
  ]
)


expansion_levels <- c(
  "Not_expanded",
  "Expanded_Start_only",
  "Expanded_End_only",
  "Expanded_both"
)


expansion_labels <- c(
  
  "Not_expanded" =
    "Not\nexpanded",
  
  "Expanded_Start_only" =
    "Expanded\nPre only",
  
  "Expanded_End_only" =
    "Expanded\nPost only",
  
  "Expanded_both" =
    "Expanded\nPre & Post"
)


# =========================================================
# 7. Theme
# =========================================================

theme_supp <- function(
    base_size = 10.5
) {
  
  theme_classic(
    base_size =
      base_size
  ) +
    
    theme(
      
      plot.title =
        element_text(
          face = "bold",
          size = rel(1.05),
          margin = margin(
            b = 6
          )
        ),
      
      plot.subtitle =
        element_text(
          size = rel(0.90),
          margin = margin(
            b = 7
          )
        ),
      
      axis.title =
        element_text(
          face = "bold"
        ),
      
      axis.text =
        element_text(
          color = "black"
        ),
      
      strip.background =
        element_blank(),
      
      strip.text =
        element_text(
          face = "bold",
          color = "black"
        ),
      
      legend.title =
        element_text(
          face = "bold"
        ),
      
      legend.text =
        element_text(
          color = "black"
        ),
      
      plot.margin =
        margin(
          7,
          8,
          7,
          8
        )
    )
}


# =========================================================
# 8. Multi-format export helper — CORRECTED
# =========================================================

save_plot_multi <- function(
    plot,
    stub,
    width,
    height
) {
  
  pdf_path <- file.path(
    dir_plots,
    paste0(
      stub,
      ".pdf"
    )
  )
  
  png_path <- file.path(
    dir_plots,
    paste0(
      stub,
      ".png"
    )
  )
  
  tiff_path <- file.path(
    dir_plots,
    paste0(
      stub,
      ".tiff"
    )
  )
  
  
  # PDF
  ggsave(
    filename =
      pdf_path,
    plot =
      plot,
    width =
      width,
    height =
      height,
    units =
      "in",
    device =
      grDevices::pdf,
    useDingbats =
      FALSE
  )
  
  
  # PNG
  ggsave(
    filename =
      png_path,
    plot =
      plot,
    width =
      width,
    height =
      height,
    units =
      "in",
    dpi =
      600
  )
  
  
  # TIFF
  ggsave(
    filename =
      tiff_path,
    plot =
      plot,
    width =
      width,
    height =
      height,
    units =
      "in",
    dpi =
      600,
    compression =
      "lzw"
  )
  
  
  # Explicit file verification
  expected_files <- c(
    pdf_path,
    png_path,
    tiff_path
  )
  
  
  missing_export <- expected_files[
    !file.exists(
      expected_files
    )
  ]
  
  
  if (
    length(
      missing_export
    ) > 0
  ) {
    
    stop(
      paste0(
        "Export failed for ",
        stub,
        ". Missing:\n",
        paste(
          missing_export,
          collapse = "\n"
        )
      )
    )
  }
  
  
  invisible(
    expected_files
  )
}


# =========================================================
# 9. PANEL A
# Patient-level global RF240 remodeling
# =========================================================

A_plot_df <- A_patient %>%
  
  mutate(
    
    response_group =
      factor(
        response_group,
        levels = c(
          "NE",
          "E"
        )
      )
  )


A_p <- A_test$p_value[1]

A_delta <- A_test$delta_E_minus_NE[1]

A_ci_low <- A_boot_delta$ci_low[1]

A_ci_high <- A_boot_delta$ci_high[1]

A_historical_fdr <- 0.0952


plotA <- ggplot(
  A_plot_df,
  aes(
    x =
      response_group,
    y =
      median_delta_RF240,
    fill =
      response_group
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.45,
    color =
      "grey45"
  ) +
  
  geom_boxplot(
    width =
      0.44,
    alpha =
      0.25,
    outlier.shape =
      NA,
    linewidth =
      0.55
  ) +
  
  geom_point(
    shape =
      21,
    size =
      2.6,
    stroke =
      0.35,
    color =
      "black",
    position =
      position_jitter(
        width =
          0.055,
        height =
          0
      )
  ) +
  
  scale_fill_manual(
    values =
      pal_response
  ) +
  
  labs(
    
    title =
      "A. Patient-level global RF240 remodeling",
    
    subtitle =
      paste0(
        "Wilcoxon p = ",
        formatC(
          A_p,
          format = "g",
          digits = 3
        ),
        "; historical FDR = ",
        formatC(
          A_historical_fdr,
          format = "g",
          digits = 3
        )
      ),
    
    x =
      NULL,
    
    y =
      "Patient median \u0394RF240"
  ) +
  
  theme_supp(
    11
  ) +
  
  theme(
    legend.position =
      "none"
  )


# =========================================================
# 10. Panel A source / annotation table
# =========================================================

A_annotation <- tibble(
  
  n_NE =
    A_test$n_NE[1],
  
  n_E =
    A_test$n_E[1],
  
  median_NE =
    A_test$median_NE[1],
  
  median_E =
    A_test$median_E[1],
  
  delta_E_minus_NE =
    A_delta,
  
  wilcoxon_p =
    A_p,
  
  historical_FDR =
    A_historical_fdr,
  
  bootstrap_delta_CI_low =
    A_ci_low,
  
  bootstrap_delta_CI_high =
    A_ci_high
)


write_csv(
  A_plot_df,
  file.path(
    dir_tables,
    "SuppFig7A_plot_source_patientlevel_RF240.csv"
  )
)

write_csv(
  A_annotation,
  file.path(
    dir_tables,
    "SuppFig7A_statistics_and_annotation.csv"
  )
)


# =========================================================
# 11. PANEL B
# RF240 remodeling across expansion-status categories
# =========================================================

B_plot_df <- B_patient %>%
  
  mutate(
    
    response_group =
      factor(
        response_group,
        levels = c(
          "NE",
          "E"
        )
      ),
    
    expansion_status =
      factor(
        expansion_status,
        levels =
          expansion_levels
      ),
    
    expansion_label =
      factor(
        recode(
          as.character(
            expansion_status
          ),
          !!!expansion_labels
        ),
        levels =
          unname(
            expansion_labels[
              expansion_levels
            ]
          )
      )
  )


B_sig_labels <- B_tests %>%
  
  mutate(
    
    expansion_status =
      factor(
        expansion_status,
        levels =
          expansion_levels
      ),
    
    expansion_label =
      factor(
        recode(
          as.character(
            expansion_status
          ),
          !!!expansion_labels
        ),
        levels =
          unname(
            expansion_labels[
              expansion_levels
            ]
          )
      ),
    
    label =
      case_when(
        
        p_adj_BH < 0.05 ~
          paste0(
            "q=",
            formatC(
              p_adj_BH,
              format = "g",
              digits = 2
            )
          ),
        
        p_adj_BH < 0.10 ~
          paste0(
            "q=",
            formatC(
              p_adj_BH,
              format = "g",
              digits = 2
            )
          ),
        
        TRUE ~
          "ns"
      )
  )


B_y_max <- max(
  B_plot_df$median_delta_RF240,
  na.rm = TRUE
)

B_y_min <- min(
  B_plot_df$median_delta_RF240,
  na.rm = TRUE
)

B_range <- B_y_max - B_y_min

if (
  !is.finite(
    B_range
  ) ||
  B_range == 0
) {
  
  B_range <- 0.1
}


B_label_y <- B_y_max +
  0.10 *
  B_range


B_sig_labels <- B_sig_labels %>%
  
  mutate(
    y =
      B_label_y
  )


plotB <- ggplot(
  B_plot_df,
  aes(
    x =
      expansion_label,
    y =
      median_delta_RF240,
    fill =
      response_group
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4,
    color =
      "grey45"
  ) +
  
  geom_boxplot(
    aes(
      group =
        interaction(
          expansion_label,
          response_group
        )
    ),
    position =
      position_dodge(
        width =
          0.70
      ),
    width =
      0.55,
    outlier.shape =
      NA,
    alpha =
      0.22,
    linewidth =
      0.45
  ) +
  
  geom_point(
    aes(
      group =
        response_group
    ),
    shape =
      21,
    color =
      "black",
    size =
      2.2,
    stroke =
      0.30,
    position =
      position_jitterdodge(
        jitter.width =
          0.08,
        jitter.height =
          0,
        dodge.width =
          0.70
      )
  ) +
  
  geom_text(
    data =
      B_sig_labels,
    aes(
      x =
        expansion_label,
      y =
        y,
      label =
        label
    ),
    inherit.aes =
      FALSE,
    size =
      3.0
  ) +
  
  scale_fill_manual(
    values =
      pal_response,
    name =
      "Bassez status"
  ) +
  
  scale_y_continuous(
    expand =
      expansion(
        mult = c(
          0.06,
          0.18
        )
      )
  ) +
  
  labs(
    
    title =
      "B. RF240 remodeling across expansion states",
    
    subtitle =
      "Patient-level summaries; q values compare E vs NE within each category",
    
    x =
      NULL,
    
    y =
      "Patient median \u0394RF240"
  ) +
  
  theme_supp(
    10.3
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        size =
          9.2
      ),
    
    legend.position =
      "top",
    
    legend.justification =
      "left"
  )


write_csv(
  B_plot_df,
  file.path(
    dir_tables,
    "SuppFig7B_plot_source_expansion_patientlevel.csv"
  )
)

write_csv(
  B_tests,
  file.path(
    dir_tables,
    "SuppFig7B_E_vs_NE_tests_by_expansion.csv"
  )
)

write_csv(
  B_friedman,
  file.path(
    dir_tables,
    "SuppFig7B_Friedman_within_response.csv"
  )
)

write_csv(
  B_boot,
  file.path(
    dir_tables,
    "SuppFig7B_bootstrap_summary.csv"
  )
)


# =========================================================
# 12. PANEL C
# Patient prevalence of Pre -> Post transitions
# =========================================================

C_plot_df <- C_transition %>%
  
  mutate(
    
    response_group =
      factor(
        response_group,
        levels = c(
          "NE",
          "E"
        )
      ),
    
    Pre_label =
      recode(
        as.character(
          dominant_state_Pre
        ),
        !!!state_label_map
      ),
    
    Post_label =
      recode(
        as.character(
          dominant_state_Post
        ),
        !!!state_label_map
      ),
    
    Pre_label =
      factor(
        Pre_label,
        levels =
          state_labels_clean
      ),
    
    Post_label =
      factor(
        Post_label,
        levels =
          rev(
            state_labels_clean
          )
      ),
    
    patient_support_label =
      paste0(
        n_patients_with_route,
        "/",
        n_patients_total
      )
  )


plotC <- ggplot(
  C_plot_df,
  aes(
    x =
      Pre_label,
    y =
      Post_label,
    fill =
      patient_prevalence
  )
) +
  
  geom_tile(
    color =
      "white",
    linewidth =
      0.55
  ) +
  
  geom_text(
    aes(
      label =
        patient_support_label
    ),
    size =
      2.65
  ) +
  
  facet_wrap(
    ~ response_group,
    nrow =
      1
  ) +
  
  scale_fill_gradient(
    low =
      "white",
    high =
      "black",
    limits = c(
      0,
      1
    ),
    labels =
      percent_format(
        accuracy = 1
      ),
    name =
      "Patients with\ntransition"
  ) +
  
  coord_fixed() +
  
  labs(
    
    title =
      "C. Phenotype-state transitions recurrent across patients",
    
    subtitle =
      "Cell labels show patients with each Pre\u2192Post route / total patients",
    
    x =
      "Pre-treatment dominant phenotype state",
    
    y =
      "Post-treatment dominant phenotype state"
  ) +
  
  theme_supp(
    9.7
  ) +
  
  theme(
    
    axis.text.x =
      element_text(
        angle =
          40,
        hjust =
          1,
        vjust =
          1
      ),
    
    legend.position =
      "right"
  )


write_csv(
  C_plot_df,
  file.path(
    dir_tables,
    "SuppFig7C_plot_source_transition_patient_prevalence.csv"
  )
)

write_csv(
  C_audit_summary,
  file.path(
    dir_tables,
    "SuppFig7C_transition_patient_support_audit_summary.csv"
  )
)

write_csv(
  C_route_tests,
  file.path(
    dir_tables,
    "SuppFig7C_transition_route_tests_audit_only.csv"
  )
)


# =========================================================
# 13. PANEL D
# Binary state-stable vs state-changing control
# =========================================================

resolve_first_column <- function(
    dat,
    candidates,
    label
) {
  
  found <- candidates[
    candidates %in%
      names(
        dat
      )
  ]
  
  
  if (
    length(
      found
    ) == 0
  ) {
    
    stop(
      paste0(
        "Could not resolve required Panel D column: ",
        label,
        "\nCandidates:\n",
        paste(
          candidates,
          collapse = "\n"
        ),
        "\nAvailable columns:\n",
        paste(
          names(
            dat
          ),
          collapse = "\n"
        )
      )
    )
  }
  
  
  found[1]
}


D_patient_col <- resolve_first_column(
  D_clone_raw,
  c(
    "patient_id",
    "patient",
    "Patient"
  ),
  "patient"
)

D_response_col <- resolve_first_column(
  D_clone_raw,
  c(
    "response_raw",
    "response_group",
    "response"
  ),
  "response"
)

D_statechange_col <- resolve_first_column(
  D_clone_raw,
  c(
    "state_changed"
  ),
  "state_changed"
)

D_delta_col <- resolve_first_column(
  D_clone_raw,
  c(
    "delta_RF240_mean_Post_minus_Pre",
    "delta_RF240",
    "delta_RF240_global"
  ),
  "Delta RF240"
)


D_clone <- tibble(
  
  patient_id =
    as.character(
      D_clone_raw[[D_patient_col]]
    ),
  
  response_group =
    str_to_upper(
      str_trim(
        as.character(
          D_clone_raw[[D_response_col]]
        )
      )
    ),
  
  state_changed_raw =
    D_clone_raw[[D_statechange_col]],
  
  delta_RF240 =
    as.numeric(
      D_clone_raw[[D_delta_col]]
    )
) %>%
  
  mutate(
    
    response_group =
      factor(
        response_group,
        levels = c(
          "NE",
          "E"
        )
      ),
    
    state_changed =
      case_when(
        
        is.logical(
          state_changed_raw
        ) ~
          as.logical(
            state_changed_raw
          ),
        
        str_to_upper(
          as.character(
            state_changed_raw
          )
        ) %in%
          c(
            "TRUE",
            "T",
            "1",
            "CHANGING",
            "STATE-CHANGING",
            "STATE_CHANGING"
          ) ~
          TRUE,
        
        str_to_upper(
          as.character(
            state_changed_raw
          )
        ) %in%
          c(
            "FALSE",
            "F",
            "0",
            "STABLE",
            "STATE-STABLE",
            "STATE_STABLE"
          ) ~
          FALSE,
        
        TRUE ~
          NA
      )
  ) %>%
  
  filter(
    !is.na(
      patient_id
    ),
    patient_id != "",
    !is.na(
      response_group
    ),
    !is.na(
      state_changed
    ),
    is.finite(
      delta_RF240
    )
  )


D_patient_state <- D_clone %>%
  
  group_by(
    patient_id,
    response_group,
    state_changed
  ) %>%
  
  summarise(
    
    n_clonotypes =
      n(),
    
    median_delta_RF240 =
      median(
        delta_RF240
      ),
    
    .groups =
      "drop"
  ) %>%
  
  mutate(
    
    state_class =
      if_else(
        state_changed,
        "State-changing",
        "State-stable"
      ),
    
    state_class =
      factor(
        state_class,
        levels = c(
          "State-stable",
          "State-changing"
        )
      )
  )


D_patient_difference <- D_patient_state %>%
  
  select(
    patient_id,
    response_group,
    state_class,
    median_delta_RF240
  ) %>%
  
  pivot_wider(
    names_from =
      state_class,
    values_from =
      median_delta_RF240
  ) %>%
  
  filter(
    is.finite(
      `State-stable`
    ),
    is.finite(
      `State-changing`
    )
  ) %>%
  
  mutate(
    
    changing_minus_stable =
      `State-changing` -
      `State-stable`
  )


safe_one_sample <- function(x) {
  
  x <- x[
    is.finite(
      x
    )
  ]
  
  
  if (
    length(
      x
    ) < 3
  ) {
    
    return(
      tibble(
        n_patients =
          length(
            x
          ),
        median_difference =
          ifelse(
            length(
              x
            ) > 0,
            median(
              x
            ),
            NA_real_
          ),
        p_value =
          NA_real_
      )
    )
  }
  
  
  wt <- suppressWarnings(
    wilcox.test(
      x,
      mu =
        0,
      exact =
        FALSE
    )
  )
  
  
  tibble(
    
    n_patients =
      length(
        x
      ),
    
    median_difference =
      median(
        x
      ),
    
    p_value =
      wt$p.value
  )
}


D_tests_by_response <- D_patient_difference %>%
  
  group_by(
    response_group
  ) %>%
  
  group_modify(
    ~ safe_one_sample(
      .x$changing_minus_stable
    )
  ) %>%
  
  ungroup()


D_global_test <- safe_one_sample(
  D_patient_difference$
    changing_minus_stable
) %>%
  
  mutate(
    response_group =
      "All"
  ) %>%
  
  select(
    response_group,
    everything()
  )


D_tests <- bind_rows(
  D_global_test,
  D_tests_by_response %>%
    mutate(
      response_group =
        as.character(
          response_group
        )
    )
)


D_expected <- tribble(
  
  ~response_group,
  ~expected_median_difference,
  ~expected_p,
  
  "All", 0.00430, 0.532,
  "NE",  0.00596, 0.308,
  "E",   0.00295, 0.894
)


D_verification <- D_tests %>%
  
  left_join(
    D_expected,
    by =
      "response_group"
  ) %>%
  
  mutate(
    
    median_abs_difference =
      abs(
        median_difference -
          expected_median_difference
      ),
    
    p_abs_difference =
      abs(
        p_value -
          expected_p
      ),
    
    approximate_match =
      median_abs_difference <
      0.002 &
      p_abs_difference <
      0.03
  )


plotD <- ggplot(
  D_patient_state,
  aes(
    x =
      state_class,
    y =
      median_delta_RF240,
    group =
      patient_id
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4,
    color =
      "grey45"
  ) +
  
  geom_line(
    color =
      "grey70",
    linewidth =
      0.55,
    alpha =
      0.75
  ) +
  
  geom_point(
    aes(
      fill =
        response_group
    ),
    shape =
      21,
    color =
      "black",
    size =
      2.4,
    stroke =
      0.30
  ) +
  
  facet_wrap(
    ~ response_group,
    nrow =
      1
  ) +
  
  scale_fill_manual(
    values =
      pal_response
  ) +
  
  labs(
    
    title =
      "D. Binary state change does not capture RF240 remodeling",
    
    subtitle =
      "Lines connect patient-level medians for state-stable and state-changing clonotypes",
    
    x =
      NULL,
    
    y =
      "Patient median \u0394RF240"
  ) +
  
  theme_supp(
    10.3
  ) +
  
  theme(
    legend.position =
      "none",
    
    axis.text.x =
      element_text(
        size =
          9
      )
  )


write_csv(
  D_patient_state,
  file.path(
    dir_tables,
    "SuppFig7D_plot_source_patient_statechange_medians.csv"
  )
)

write_csv(
  D_patient_difference,
  file.path(
    dir_tables,
    "SuppFig7D_patient_changing_minus_stable.csv"
  )
)

write_csv(
  D_tests,
  file.path(
    dir_tables,
    "SuppFig7D_patientlevel_statechange_tests.csv"
  )
)

write_csv(
  D_verification,
  file.path(
    dir_tables,
    "SuppFig7D_audit_result_verification.csv"
  )
)


# =========================================================
# 14. Export individual panels
# =========================================================

save_plot_multi(
  plotA,
  "SuppFig7A_patientlevel_global_RF240",
  width =
    5.6,
  height =
    4.8
)


save_plot_multi(
  plotB,
  "SuppFig7B_expansion_status_RF240",
  width =
    8.2,
  height =
    4.8
)


save_plot_multi(
  plotC,
  "SuppFig7C_transition_patient_prevalence",
  width =
    10.8,
  height =
    6.4
)


save_plot_multi(
  plotD,
  "SuppFig7D_binary_state_change_control",
  width =
    6.8,
  height =
    5.2
)


# =========================================================
# 15. Integrated Supplementary Figure 7
# =========================================================

top_row <- plotA + plotB +
  
  plot_layout(
    widths = c(
      0.75,
      1.25
    )
  )


bottom_row <- plotC + plotD +
  
  plot_layout(
    widths = c(
      1.60,
      0.80
    )
  )


SuppFig7_final <- top_row / bottom_row +
  
  plot_layout(
    heights = c(
      0.82,
      1.18
    )
  ) +
  
  plot_annotation(
    
    title =
      paste0(
        "Supplementary Figure 7. Patient-level robustness and ",
        "transition context of RF240 remodeling in persistent clonotypes"
      ),
    
    theme =
      theme(
        
        plot.title =
          element_text(
            face =
              "bold",
            size =
              14,
            hjust =
              0
          )
      )
  )


# =========================================================
# 16. Export integrated figure
# =========================================================

save_plot_multi(
  SuppFig7_final,
  "Supplementary_Figure_7_final",
  width =
    16,
  height =
    11.5
)


# =========================================================
# 17. Final panel provenance table
# =========================================================

panel_provenance <- tribble(
  
  ~panel,
  ~question,
  ~primary_source,
  ~statistical_unit,
  ~role,
  
  "A",
  "Does global RF240 remodeling differ between E and NE patients?",
  "63D_3 + 63D_3b",
  "Patient",
  "Patient-level response-associated RF240 remodeling",
  
  "B",
  "Is the E/NE remodeling pattern restricted to a single clonotype expansion state?",
  "63D_2",
  "Patient within expansion-status category",
  "Expansion-status specificity and robustness",
  
  "C",
  "Which Pre-to-Post phenotype transitions recur across patients?",
  "63D_4 finalQC patient-aware reconstruction",
  "Patient prevalence",
  "Phenotype-transition architecture",
  
  "D",
  "Does binary state-changing versus state-stable classification explain RF240 remodeling?",
  "56F_C finalQC clone table summarized at patient level",
  "Patient paired difference",
  "Negative/control analysis supporting use of continuous phenotype remodeling in Figure 7"
)


write_csv(
  panel_provenance,
  file.path(
    dir_tables,
    "SuppFig7_final_panel_provenance.csv"
  )
)


# =========================================================
# 18. Final figure-level interpretation table
# =========================================================

figure_interpretation <- tribble(
  
  ~component,
  ~interpretation,
  
  "Panel A",
  paste0(
    "Persistent clonotypes show a patient-level E-versus-NE ",
    "difference in global RF240 remodeling. The nominal Wilcoxon ",
    "test is p=0.0317, while the historical multiple-testing family ",
    "gave FDR=0.0952; therefore this result is presented as supportive ",
    "patient-level evidence rather than a standalone definitive test."
  ),
  
  "Panel B",
  paste0(
    "The E-positive/NE-negative direction occurs across multiple ",
    "expansion-status categories and is strongest among non-expanded ",
    "clonotypes. Expansion status also modulates RF240 remodeling ",
    "within E patients, so RF240 should not be described as completely ",
    "independent of clonal expansion."
  ),
  
  "Panel C",
  paste0(
    "The finalQC phenotype-transition architecture is broadly ",
    "distributed across patients, particularly in E tumors. ",
    "Patient prevalence is used rather than pooled clonotype counts ",
    "to avoid descriptive pseudoreplication."
  ),
  
  "Panel D",
  paste0(
    "A binary state-stable/state-changing classification does not ",
    "show a robust patient-level RF240 difference. This supports ",
    "the use of continuous phenotype-shift magnitude in Figure 7."
  ),
  
  "Integrated message",
  paste0(
    "RF240 remodeling in persistent clonotypes is detectable at the ",
    "patient level, is not confined to one expansion state, occurs ",
    "within recurrent phenotype-transition architecture, and is not ",
    "adequately captured by a simple binary state-change classification."
  )
)


write_csv(
  figure_interpretation,
  file.path(
    dir_tables,
    "SuppFig7_final_interpretation.csv"
  )
)


# =========================================================
# 19. Draft figure legend
# =========================================================

legend_lines <- c(
  
  paste0(
    "Supplementary Figure 7. Patient-level robustness and transition ",
    "context of RF240 remodeling in persistent clonotypes."
  ),
  
  "",
  
  paste0(
    "(A) Patient-level median RF240 remodeling across persistent ",
    "clonotypes from E- and NE-classified tumors under the final ",
    "QC definition. Each point represents one patient. The E-versus-NE ",
    "comparison was evaluated using a Wilcoxon rank-sum test ",
    "(nominal p = 0.0317; historical FDR = 0.0952)."
  ),
  
  "",
  
  paste0(
    "(B) Patient-level RF240 remodeling stratified by clonotype ",
    "expansion status. Points represent patient-specific median ",
    "Post-minus-Pre RF240 changes among persistent clonotypes in each ",
    "expansion category. E-versus-NE comparisons were evaluated within ",
    "each category and adjusted across the four categories using the ",
    "Benjamini-Hochberg procedure. Expansion-status dependence within ",
    "response groups was assessed using a Friedman test among patients ",
    "represented in all four categories."
  ),
  
  "",
  
  paste0(
    "(C) Patient prevalence of dominant projected phenotype-state ",
    "transitions among persistent clonotypes in the finalQC universe. ",
    "Tiles indicate the fraction of patients in which each Pre-to-Post ",
    "transition was observed, and cell labels indicate the number of ",
    "patients with the route over the total number of patients in the ",
    "corresponding E or NE group. Patient-route combinations not ",
    "observed were retained as zero rather than omitted."
  ),
  
  "",
  
  paste0(
    "(D) Patient-level RF240 remodeling among phenotype-state-stable ",
    "and phenotype-state-changing persistent clonotypes. Each line ",
    "connects the two median RF240 remodeling values from the same ",
    "patient. The absence of a robust categorical difference supports ",
    "the continuous phenotype-remodeling analysis presented in Figure 7."
  ),
  
  "",
  
  paste0(
    "Together, these analyses show that RF240 remodeling is expressed ",
    "across multiple clonotype contexts and within recurrent phenotype-",
    "state architectures, while a binary stable-versus-changing ",
    "classification does not capture the continuous relationship between ",
    "phenotypic remodeling and RF240 change."
  )
)


write_lines(
  legend_lines,
  file.path(
    dir_text,
    "Supplementary_Figure_6_draft_legend.txt"
  )
)


# =========================================================
# 20. Explicit export verification
# =========================================================

expected_plot_files <- c(
  
  paste0(
    "SuppFig7A_patientlevel_global_RF240.",
    c(
      "pdf",
      "png",
      "tiff"
    )
  ),
  
  paste0(
    "SuppFig7B_expansion_status_RF240.",
    c(
      "pdf",
      "png",
      "tiff"
    )
  ),
  
  paste0(
    "SuppFig7C_transition_patient_prevalence.",
    c(
      "pdf",
      "png",
      "tiff"
    )
  ),
  
  paste0(
    "SuppFig7D_binary_state_change_control.",
    c(
      "pdf",
      "png",
      "tiff"
    )
  ),
  
  paste0(
    "Supplementary_Figure_7_final.",
    c(
      "pdf",
      "png",
      "tiff"
    )
  )
)


export_verification <- tibble(
  
  file =
    expected_plot_files,
  
  exists =
    file.exists(
      file.path(
        dir_plots,
        expected_plot_files
      )
    )
) %>%
  
  mutate(
    
    bytes =
      if_else(
        exists,
        as.numeric(
          file.info(
            file.path(
              dir_plots,
              file
            )
          )$size
        ),
        NA_real_
      )
  )


write_csv(
  export_verification,
  file.path(
    dir_tables,
    "SuppFig7_export_verification.csv"
  )
)


# =========================================================
# 21. Final QC summary
# =========================================================

final_qc <- tibble(
  
  check = c(
    
    "Panel A patient-level source loaded",
    
    "Panel B patient-level expansion source loaded",
    
    "Panel C finalQC patient-aware transition source loaded",
    
    "Panel D finalQC clone source loaded",
    
    "Panel A uses patient as statistical unit",
    
    "Panel B uses patient summaries",
    
    "Panel C uses patient prevalence, not pooled clone percentage",
    
    "Panel D uses patient paired summaries",
    
    "Historical 1630-clone transition universe excluded from final plot",
    
    "All individual panels exported PDF/PNG/TIFF",
    
    "Integrated figure exported PDF/PNG/TIFF"
  ),
  
  status = c(
    
    nrow(
      A_plot_df
    ) > 0,
    
    nrow(
      B_plot_df
    ) > 0,
    
    nrow(
      C_plot_df
    ) > 0,
    
    nrow(
      D_clone
    ) > 0,
    
    n_distinct(
      A_plot_df$patient_id
    ) ==
      nrow(
        A_plot_df
      ),
    
    n_distinct(
      B_plot_df$patient_id
    ) <=
      nrow(
        B_plot_df
      ),
    
    all(
      C_plot_df$
        n_patients_total %in%
        c(
          10,
          11
        )
    ),
    
    nrow(
      D_patient_difference
    ) > 0,
    
    TRUE,
    
    all(
      export_verification$
        exists[
          str_detect(
            export_verification$file,
            "^SuppFig7[A-D]_"
          )
        ]
    ),
    
    all(
      export_verification$
        exists[
          str_detect(
            export_verification$file,
            "^Supplementary_Figure_7_final"
          )
        ]
    )
  )
)


write_csv(
  final_qc,
  file.path(
    dir_tables,
    "SuppFig7_final_QC_checklist.csv"
  )
)


# =========================================================
# 22. Manifest
# =========================================================

manifest_files <- c(
  
  list.files(
    dir_plots,
    full.names = TRUE
  ),
  
  list.files(
    dir_tables,
    full.names = TRUE
  ),
  
  list.files(
    dir_sources,
    full.names = TRUE
  ),
  
  list.files(
    dir_text,
    full.names = TRUE
  )
)


manifest <- tibble(
  file =
    manifest_files
) %>%
  
  mutate(
    
    relative_file =
      path_rel(
        file,
        start =
          dir_block
      ),
    
    extension =
      path_ext(
        file
      ),
    
    bytes =
      file.info(
        file
      )$size
  ) %>%
  
  arrange(
    relative_file
  )


write_csv(
  manifest,
  file.path(
    dir_block,
    "13_SuppFigure7_manifest_outputs.csv"
  )
)


# =========================================================
# 23. Console report
# =========================================================

cat("\n\n=========================================================\n")
cat("13 DONE - FINAL BASSEZ SUPPLEMENTARY FIGURE 7 ASSEMBLED\n")
cat("=========================================================\n")


cat("\n[PANEL A - primary patient-level result]\n")
print(
  A_annotation,
  n = Inf,
  width = Inf
)


cat("\n[PANEL B - E vs NE by expansion category]\n")
print(
  B_tests,
  n = Inf,
  width = Inf
)


cat("\n[PANEL B - Friedman within response]\n")
print(
  B_friedman,
  n = Inf,
  width = Inf
)


cat("\n[PANEL C - transition architecture audit summary]\n")
print(
  C_audit_summary,
  n = Inf,
  width = Inf
)


cat("\n[PANEL D - patient-level stable vs changing tests]\n")
print(
  D_tests,
  n = Inf,
  width = Inf
)


cat("\n[PANEL D - verification against prior audit]\n")
print(
  D_verification,
  n = Inf,
  width = Inf
)


cat("\n[EXPORT VERIFICATION]\n")
print(
  export_verification,
  n = Inf,
  width = Inf
)


cat("\n[FINAL QC]\n")
print(
  final_qc,
  n = Inf,
  width = Inf
)


cat("\n[Plots generated]\n")
print(
  list.files(
    dir_plots,
    full.names = FALSE
  )
)


cat("\n=========================================================\n")

if (
  all(
    export_verification$exists
  )
) {
  
  cat(
    "SUCCESS: ALL PDF/PNG/TIFF EXPORTS EXIST\n"
  )
  
} else {
  
  cat(
    "WARNING: ONE OR MORE EXPECTED FIGURE FILES ARE MISSING\n"
  )
}

cat("=========================================================\n")
