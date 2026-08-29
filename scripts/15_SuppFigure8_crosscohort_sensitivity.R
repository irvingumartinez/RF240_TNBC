# =========================================================
# 15_SuppFigure8_crosscohort_sensitivity.R
#
# FINAL GRAPHICAL REVISION ONLY
#
# IMPORTANT
# ---------
# This V2 does NOT change:
# - source tables
# - statistical tests
# - effect estimates
# - biological interpretation
# - panel architecture
#
# Graphical changes only:
# 1. Restore explicit Figure 8B panel heading.
# 2. Replace problematic Unicode symbols with safe ASCII text.
# 3. Improve Zhang patient labels.
# 4. Improve Supplementary Figure 7F readability.
# 5. Harmonize typography and spacing.
#
# =========================================================


# =========================================================
# 0. Clean session
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
  library(ggplot2)
  library(patchwork)
})


# =========================================================
# 2. Global plotting constants
# =========================================================

COL_NR <- "#D55E00"
COL_R  <- "#0072B2"

response_colors <- c(
  "NR" = COL_NR,
  "R"  = COL_R
)


theme_final <- theme_classic(
  base_size = 11
) +
  
  theme(
    
    axis.title =
      element_text(
        size = 11
      ),
    
    axis.text =
      element_text(
        size = 9
      ),
    
    plot.title =
      element_text(
        size = 12,
        face = "bold",
        margin = margin(
          b = 3
        )
      ),
    
    plot.subtitle =
      element_text(
        size = 9.5,
        margin = margin(
          b = 8
        )
      ),
    
    plot.caption =
      element_text(
        size = 8,
        hjust = 0
      ),
    
    legend.title =
      element_blank(),
    
    legend.position =
      "top",
    
    strip.background =
      element_blank(),
    
    strip.text =
      element_text(
        size = 10,
        face = "bold"
      ),
    
    plot.margin =
      margin(
        8,
        8,
        8,
        8
      )
  )


# =========================================================
# =========================================================
# 3. Repository paths
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

dir_in <- file.path(
  project_root,
  "data",
  "longitudinal",
  "figure8_crosscohort"
)

dir_supp7 <- file.path(
  project_root,
  "results",
  "supplementary_figures",
  "Supplementary_Figure8"
)

dir_supp7_panels <- file.path(
  dir_supp7,
  "individual_panels"
)

dir_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Supplementary_Figure8_crosscohort"
)

# Compatibility aliases from the historical combined script.
dir_out <- dir_supp7
dir_fig8 <- dir_supp7
dir_fig8_panels <- dir_supp7_panels

walk(
  c(
    dir_supp7,
    dir_supp7_panels,
    dir_tables
  ),
  dir_create
)

# 4. Input files
# =========================================================

file_A <- file.path(
  dir_in,
  "64E_Fig8A_Zhang_patient_points.csv"
)

file_B <- file.path(
  dir_in,
  "64E_Fig8B_Zhang_robustness.csv"
)

file_C <- file.path(
  dir_in,
  "64E_Fig8C_Shiao_patient_points.csv"
)

file_D <- file.path(
  dir_in,
  "64E_Fig8D_crosscohort_effects_bootstrap.csv"
)


file_Z_state <- file.path(
  dir_in,
  "64E_Supp7_Zhang_statechange_patientlevel.csv"
)

file_Z_rho <- file.path(
  dir_in,
  "64E_Supp7_Zhang_patient_rho.csv"
)

file_S_qc <- file.path(
  dir_in,
  "64E_Supp7_Shiao_QC_sensitivity.csv"
)

file_S_state <- file.path(
  dir_in,
  "64E_Supp7_Shiao_statechange_patientlevel.csv"
)

file_S_rho <- file.path(
  dir_in,
  "64E_Supp7_Shiao_patient_rho.csv"
)

file_S_expansion <- file.path(
  dir_in,
  "64E_Supp7_Shiao_expansion_patientlevel.csv"
)


required_files <- c(
  file_A,
  file_B,
  file_C,
  file_D,
  file_Z_state,
  file_Z_rho,
  file_S_qc,
  file_S_state,
  file_S_rho,
  file_S_expansion
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
      "Missing 64E input tables:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}


# =========================================================
# 5. Read tables
# =========================================================

A_dat <- read_csv(
  file_A,
  show_col_types = FALSE
)

B_dat <- read_csv(
  file_B,
  show_col_types = FALSE
)

C_dat <- read_csv(
  file_C,
  show_col_types = FALSE
)

D_dat <- read_csv(
  file_D,
  show_col_types = FALSE
)


Z_state <- read_csv(
  file_Z_state,
  show_col_types = FALSE
)

Z_rho <- read_csv(
  file_Z_rho,
  show_col_types = FALSE
)

S_qc <- read_csv(
  file_S_qc,
  show_col_types = FALSE
)

S_state <- read_csv(
  file_S_state,
  show_col_types = FALSE
)

S_rho <- read_csv(
  file_S_rho,
  show_col_types = FALSE
)

S_expansion <- read_csv(
  file_S_expansion,
  show_col_types = FALSE
)


# =========================================================
# 6. Export helper
# =========================================================

export_plot <- function(
    plot_object,
    output_directory,
    filename_stub,
    width,
    height,
    dpi = 600
) {
  
  pdf_path <- file.path(
    output_directory,
    paste0(
      filename_stub,
      ".pdf"
    )
  )
  
  png_path <- file.path(
    output_directory,
    paste0(
      filename_stub,
      ".png"
    )
  )
  
  tiff_path <- file.path(
    output_directory,
    paste0(
      filename_stub,
      ".tiff"
    )
  )
  
  
  grDevices::pdf(
    file =
      pdf_path,
    width =
      width,
    height =
      height,
    useDingbats =
      FALSE
  )
  
  print(
    plot_object
  )
  
  grDevices::dev.off()
  
  
  ggsave(
    filename =
      png_path,
    plot =
      plot_object,
    width =
      width,
    height =
      height,
    units =
      "in",
    dpi =
      dpi
  )
  
  
  ggsave(
    filename =
      tiff_path,
    plot =
      plot_object,
    width =
      width,
    height =
      height,
    units =
      "in",
    dpi =
      dpi,
    compression =
      "lzw"
  )
  
  
  invisible(
    c(
      pdf_path,
      png_path,
      tiff_path
    )
  )
}


# =========================================================
# 12. SUPPLEMENTARY FIGURE 7A
#
# Zhang stable vs changing
# =========================================================

Z_state_long <- Z_state %>%
  
  pivot_longer(
    
    cols = c(
      `State-stable`,
      `State-changing`
    ),
    
    names_to =
      "state_class",
    
    values_to =
      "median_delta_RF240"
  ) %>%
  
  mutate(
    
    state_class =
      factor(
        state_class,
        levels = c(
          "State-stable",
          "State-changing"
        )
      ),
    
    response =
      factor(
        response,
        levels = c(
          "NR",
          "R"
        )
      )
  )


pS7A <- ggplot(
  Z_state_long,
  aes(
    x =
      state_class,
    y =
      median_delta_RF240,
    group =
      patient_short,
    color =
      response
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +
  
  geom_line(
    alpha =
      0.55,
    linewidth =
      0.5
  ) +
  
  geom_point(
    size =
      2.4
  ) +
  
  scale_color_manual(
    values =
      response_colors
  ) +
  
  labs(
    
    title =
      "A  Zhang: binary phenotype-state change",
    
    subtitle =
      "Patient-paired stable vs changing clonotypes; global p = 0.281",
    
    x =
      NULL,
    
    y =
      "Delta RF240 (patient median)"
  ) +
  
  theme_final


# =========================================================
# 13. SUPPLEMENTARY FIGURE 7B
#
# Zhang patient-specific phenotype/RF240 rho
# =========================================================

Z_rho <- Z_rho %>%
  
  mutate(
    
    response =
      factor(
        response,
        levels = c(
          "NR",
          "R"
        )
      )
  )


pS7B <- ggplot(
  Z_rho,
  aes(
    x =
      response,
    y =
      spearman_rho,
    color =
      response
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +
  
  geom_jitter(
    width =
      0.08,
    size =
      2.8
  ) +
  
  geom_text(
    aes(
      label =
        patient_short
    ),
    position =
      position_nudge(
        x = 0.12
      ),
    size =
      2.7,
    show.legend =
      FALSE
  ) +
  
  scale_color_manual(
    values =
      response_colors
  ) +
  
  labs(
    
    title =
      "B  Zhang: patient-specific phenotype/RF240 coupling",
    
    subtitle =
      paste0(
        "Median rho: NR = 0.010, R = -0.053; ",
        "R vs NR p = 0.80"
      ),
    
    x =
      NULL,
    
    y =
      "Within-patient Spearman rho"
  ) +
  
  theme_final +
  
  theme(
    legend.position =
      "none"
  )


# =========================================================
# 14. SUPPLEMENTARY FIGURE 7C
#
# Shiao QC sensitivity
# =========================================================

qc_order <- c(
  "no_filter",
  "min50cells_min5persistent",
  "min100cells_min10persistent",
  "min200cells_min20persistent",
  "min300cells_min30persistent"
)


S_qc <- S_qc %>%
  
  mutate(
    
    qc_level =
      factor(
        qc_level,
        levels =
          qc_order
      ),
    
    qc_label =
      factor(
        qc_level,
        levels =
          qc_order,
        labels = c(
          "No filter",
          ">=50 / >=5",
          ">=100 / >=10",
          ">=200 / >=20",
          ">=300 / >=30"
        )
      ),
    
    comparison_label =
      case_when(
        
        comparison ==
          "Base_vs_PD1" ~
          "Base - PD1",
        
        comparison ==
          "Base_vs_RTPD1" ~
          "Base - RTPD1",
        
        comparison ==
          "PD1_vs_RTPD1" ~
          "PD1 - RTPD1",
        
        TRUE ~
          comparison
      )
  )


pS7C <- ggplot(
  S_qc,
  aes(
    x =
      qc_label,
    y =
      delta_R_minus_NR,
    group =
      comparison_label,
    linetype =
      comparison_label,
    shape =
      comparison_label
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +
  
  geom_line(
    linewidth =
      0.65
  ) +
  
  geom_point(
    size =
      2.3
  ) +
  
  labs(
    
    title =
      "C  Shiao: QC sensitivity of response effects",
    
    subtitle =
      "Nested thresholds are robustness checks, not independent hypothesis tests",
    
    x =
      "QC threshold",
    
    y =
      "Delta RF240 effect (R - NR)",
    
    linetype =
      "Interval",
    
    shape =
      "Interval"
  ) +
  
  theme_final +
  
  theme(
    axis.text.x =
      element_text(
        angle = 32,
        hjust = 1,
        size = 8
      )
  )


# =========================================================
# 15. SUPPLEMENTARY FIGURE 7D
#
# Shiao binary state change
# =========================================================

S_state <- S_state %>%
  
  mutate(
    
    response =
      factor(
        response,
        levels = c(
          "NR",
          "R"
        )
      ),
    
    comparison_label =
      case_when(
        
        comparison ==
          "Base_vs_PD1" ~
          "Base - PD1",
        
        comparison ==
          "Base_vs_RTPD1" ~
          "Base - RTPD1",
        
        comparison ==
          "PD1_vs_RTPD1" ~
          "PD1 - RTPD1",
        
        TRUE ~
          comparison
      ),
    
    comparison_label =
      factor(
        comparison_label,
        levels = c(
          "Base - PD1",
          "Base - RTPD1",
          "PD1 - RTPD1"
        )
      )
  )


pS7D <- ggplot(
  S_state,
  aes(
    x =
      response,
    y =
      changing_minus_stable,
    color =
      response
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +
  
  geom_jitter(
    width =
      0.12,
    size =
      1.8,
    alpha =
      0.72
  ) +
  
  stat_summary(
    fun =
      median,
    geom =
      "crossbar",
    width =
      0.4,
    linewidth =
      0.6
  ) +
  
  facet_wrap(
    ~ comparison_label,
    nrow =
      1
  ) +
  
  scale_color_manual(
    values =
      response_colors
  ) +
  
  labs(
    
    title =
      "D  Shiao: binary phenotype-state change",
    
    subtitle =
      paste0(
        "No interval shows a patient-level difference between ",
        "changing and stable clonotypes"
      ),
    
    x =
      NULL,
    
    y =
      "Delta RF240 changing - stable"
  ) +
  
  theme_final +
  
  theme(
    legend.position =
      "none"
  )


# =========================================================
# 16. SUPPLEMENTARY FIGURE 7E
#
# Shiao patient-specific phenotype/RF240 rho
# =========================================================

S_rho <- S_rho %>%
  
  mutate(
    
    response =
      factor(
        response,
        levels = c(
          "NR",
          "R"
        )
      ),
    
    comparison_label =
      case_when(
        
        comparison ==
          "Base_vs_PD1" ~
          "Base - PD1",
        
        comparison ==
          "Base_vs_RTPD1" ~
          "Base - RTPD1",
        
        comparison ==
          "PD1_vs_RTPD1" ~
          "PD1 - RTPD1",
        
        TRUE ~
          comparison
      ),
    
    comparison_label =
      factor(
        comparison_label,
        levels = c(
          "Base - PD1",
          "Base - RTPD1",
          "PD1 - RTPD1"
        )
      )
  )


pS7E <- ggplot(
  S_rho,
  aes(
    x =
      response,
    y =
      spearman_rho,
    color =
      response
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +
  
  geom_jitter(
    width =
      0.12,
    size =
      1.8,
    alpha =
      0.72
  ) +
  
  stat_summary(
    fun =
      median,
    geom =
      "crossbar",
    width =
      0.4,
    linewidth =
      0.6
  ) +
  
  facet_wrap(
    ~ comparison_label,
    nrow =
      1
  ) +
  
  scale_color_manual(
    values =
      response_colors
  ) +
  
  labs(
    
    title =
      "E  Shiao: patient-specific phenotype/RF240 coupling",
    
    subtitle =
      paste0(
        "No interval shows an R vs NR difference in ",
        "patient-specific correlations"
      ),
    
    x =
      NULL,
    
    y =
      "Within-patient Spearman rho"
  ) +
  
  theme_final +
  
  theme(
    legend.position =
      "none"
  )


# =========================================================
# 17. SUPPLEMENTARY FIGURE 7F
#
# Shiao patient-aware expansion context
# =========================================================

expansion_order <- c(
  "Not_expanded",
  "Expanded_Start_only",
  "Expanded_End_only",
  "Expanded_both"
)


S_expansion <- S_expansion %>%
  
  mutate(
    
    response =
      factor(
        response,
        levels = c(
          "NR",
          "R"
        )
      ),
    
    expansion_status =
      factor(
        expansion_status,
        levels =
          expansion_order,
        labels = c(
          "Not expanded",
          "Start only",
          "End only",
          "Both"
        )
      ),
    
    comparison_label =
      case_when(
        
        comparison ==
          "Base_vs_PD1" ~
          "Base - PD1",
        
        comparison ==
          "Base_vs_RTPD1" ~
          "Base - RTPD1",
        
        comparison ==
          "PD1_vs_RTPD1" ~
          "PD1 - RTPD1",
        
        TRUE ~
          comparison
      ),
    
    comparison_label =
      factor(
        comparison_label,
        levels = c(
          "Base - PD1",
          "Base - RTPD1",
          "PD1 - RTPD1"
        )
      )
  )


# Median table for cleaner emphasis
S_expansion_median <- S_expansion %>%
  
  group_by(
    comparison_label,
    expansion_status,
    response
  ) %>%
  
  summarise(
    
    median_delta_RF240 =
      median(
        median_delta_RF240
      ),
    
    .groups =
      "drop"
  )


pS7F <- ggplot(
  S_expansion,
  aes(
    x =
      expansion_status,
    y =
      median_delta_RF240,
    color =
      response
  )
) +
  
  geom_hline(
    yintercept =
      0,
    linetype =
      "dashed",
    linewidth =
      0.4
  ) +
  
  geom_jitter(
    position =
      position_jitterdodge(
        jitter.width =
          0.10,
        dodge.width =
          0.42
      ),
    size =
      1.2,
    alpha =
      0.35
  ) +
  
  geom_point(
    data =
      S_expansion_median,
    aes(
      x =
        expansion_status,
      y =
        median_delta_RF240,
      color =
        response,
      group =
        response
    ),
    position =
      position_dodge(
        width =
          0.42
      ),
    size =
      3.0,
    shape =
      18,
    inherit.aes =
      FALSE
  ) +
  
  facet_wrap(
    ~ comparison_label,
    nrow =
      1
  ) +
  
  scale_color_manual(
    values =
      response_colors
  ) +
  
  labs(
    
    title =
      "F  Shiao: expansion-status context",
    
    subtitle =
      paste0(
        "No expansion category rescues a consistent ",
        "response-associated RF240 difference"
      ),
    
    x =
      "Expansion status",
    
    y =
      "Delta RF240 (patient median)"
  ) +
  
  theme_final +
  
  theme(
    
    axis.text.x =
      element_text(
        angle = 28,
        hjust = 1,
        size = 8
      ),
    
    legend.position =
      "top"
  )


# =========================================================
# 18. Assemble Supplementary Figure 7
# =========================================================

SuppFig7 <- (
  pS7A |
    pS7B
) /
  (
    pS7C |
      pS7D
  ) /
  (
    pS7E |
      pS7F
  ) +
  
  plot_layout(
    heights =
      c(
        1,
        1,
        1
      )
  ) +
  
  plot_annotation(
    
    title =
      paste0(
        "Supplementary Figure 8. Patient-level robustness and ",
        "negative-control analyses of longitudinal RF240 remodeling"
      ),
    
    theme =
      theme(
        
        plot.title =
          element_text(
            size = 14,
            face = "bold",
            margin = margin(
              b = 8
            )
          )
      )
  )


# =========================================================


# =========================================================
# 19. Export Supplementary Figure 8 panels
# =========================================================

export_plot(
  pS7A,
  dir_supp7_panels,
  "SuppFig8A_Zhang_statechange",
  width = 5.8,
  height = 4.5
)

export_plot(
  pS7B,
  dir_supp7_panels,
  "SuppFig8B_Zhang_patient_rho",
  width = 5.8,
  height = 4.5
)

export_plot(
  pS7C,
  dir_supp7_panels,
  "SuppFig8C_Shiao_QC_sensitivity",
  width = 7.6,
  height = 4.8
)

export_plot(
  pS7D,
  dir_supp7_panels,
  "SuppFig8D_Shiao_statechange",
  width = 8.6,
  height = 4.5
)

export_plot(
  pS7E,
  dir_supp7_panels,
  "SuppFig8E_Shiao_patient_rho",
  width = 8.6,
  height = 4.5
)

export_plot(
  pS7F,
  dir_supp7_panels,
  "SuppFig8F_Shiao_expansion",
  width = 9.8,
  height = 4.9
)


# =========================================================
# 20. Export integrated Supplementary Figure 8
# =========================================================

export_plot(
  SuppFig7,
  dir_supp7,
  "Supplementary_Figure8_FINAL_integrated",
  width = 14.2,
  height = 15.3
)


# =========================================================
# 21. Source map
# =========================================================

figure_source_map <- tribble(
  ~figure, ~panel, ~source_table, ~interpretation,
  "Supplementary Figure 8", "A",
  "64E_Supp7_Zhang_statechange_patientlevel.csv",
  "Binary phenotype-state negative control in Zhang",
  "Supplementary Figure 8", "B",
  "64E_Supp7_Zhang_patient_rho.csv",
  "Patient-specific phenotype/RF240 coupling in Zhang",
  "Supplementary Figure 8", "C",
  "64E_Supp7_Shiao_QC_sensitivity.csv",
  "Shiao nested-QC robustness",
  "Supplementary Figure 8", "D",
  "64E_Supp7_Shiao_statechange_patientlevel.csv",
  "Binary phenotype-state negative control in Shiao",
  "Supplementary Figure 8", "E",
  "64E_Supp7_Shiao_patient_rho.csv",
  "Patient-specific phenotype/RF240 coupling in Shiao",
  "Supplementary Figure 8", "F",
  "64E_Supp7_Shiao_expansion_patientlevel.csv",
  "Patient-aware expansion context"
)

write_csv(
  figure_source_map,
  file.path(
    dir_tables,
    "15_SuppFigure8_source_map.csv"
  )
)


# =========================================================
# 22. Export verification
# =========================================================

expected_stubs <- c(
  file.path(dir_supp7_panels, "SuppFig8A_Zhang_statechange"),
  file.path(dir_supp7_panels, "SuppFig8B_Zhang_patient_rho"),
  file.path(dir_supp7_panels, "SuppFig8C_Shiao_QC_sensitivity"),
  file.path(dir_supp7_panels, "SuppFig8D_Shiao_statechange"),
  file.path(dir_supp7_panels, "SuppFig8E_Shiao_patient_rho"),
  file.path(dir_supp7_panels, "SuppFig8F_Shiao_expansion"),
  file.path(dir_supp7, "Supplementary_Figure8_FINAL_integrated")
)

expected_files <- expand_grid(
  stub = expected_stubs,
  extension = c("pdf", "png", "tiff")
) %>%
  mutate(
    file = paste0(stub, ".", extension),
    exists = file.exists(file),
    size_bytes = ifelse(
      exists,
      file.info(file)$size,
      NA_real_
    )
  )

write_csv(
  expected_files,
  file.path(
    dir_tables,
    "15_SuppFigure8_export_verification.csv"
  )
)

cat("\n=========================================================\n")
cat("15 DONE - SUPPLEMENTARY FIGURE 8 CROSS-COHORT SENSITIVITY\n")
cat("=========================================================\n")
cat("Expected exports present:", sum(expected_files$exists),
    "/", nrow(expected_files), "\n")
cat("=========================================================\n")
