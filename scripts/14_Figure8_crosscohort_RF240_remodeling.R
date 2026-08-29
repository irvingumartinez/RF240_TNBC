# =========================================================
# 14_Figure8_crosscohort_RF240_remodeling.R
#
# FINAL GRAPHICAL FIX — FIGURE 8 ONLY
#
# No statistics are recalculated.
# No source tables are modified.
#
# Fixes:
# 1. Explicit panel B heading that survives patchwork nesting.
# 2. Move P005 label away from Zhang median bar.
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
# 2. Constants
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
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 9),
    plot.title = element_text(
      size = 12,
      face = "bold",
      margin = margin(b = 3)
    ),
    plot.subtitle = element_text(
      size = 9.5,
      margin = margin(b = 8)
    ),
    plot.caption = element_text(
      size = 8,
      hjust = 0
    ),
    legend.title = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(
      size = 10,
      face = "bold"
    ),
    plot.margin = margin(8, 8, 8, 8)
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

dir_out <- file.path(
  project_root,
  "results",
  "figures",
  "Figure8"
)

dir_panels <- file.path(
  dir_out,
  "individual_panels"
)

dir_create(dir_out)
dir_create(dir_panels)

# 4. Inputs
# =========================================================

A_dat <- read_csv(
  file.path(
    dir_in,
    "64E_Fig8A_Zhang_patient_points.csv"
  ),
  show_col_types = FALSE
)

B_dat <- read_csv(
  file.path(
    dir_in,
    "64E_Fig8B_Zhang_robustness.csv"
  ),
  show_col_types = FALSE
)

C_dat <- read_csv(
  file.path(
    dir_in,
    "64E_Fig8C_Shiao_patient_points.csv"
  ),
  show_col_types = FALSE
)

D_dat <- read_csv(
  file.path(
    dir_in,
    "64E_Fig8D_crosscohort_effects_bootstrap.csv"
  ),
  show_col_types = FALSE
)


# =========================================================
# 5. Export helper
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
    paste0(filename_stub, ".pdf")
  )
  
  png_path <- file.path(
    output_directory,
    paste0(filename_stub, ".png")
  )
  
  tiff_path <- file.path(
    output_directory,
    paste0(filename_stub, ".tiff")
  )
  
  
  grDevices::pdf(
    file = pdf_path,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  print(plot_object)
  
  grDevices::dev.off()
  
  
  ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )
  
  
  ggsave(
    filename = tiff_path,
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw"
  )
}


# =========================================================
# 6. FIGURE 8A — final
# =========================================================

A_dat <- A_dat %>%
  
  mutate(
    response = factor(
      response,
      levels = c("NR", "R")
    ),
    
    # Explicit label positions.
    # P005 moved downward/right to avoid median line.
    label_x = case_when(
      patient %in% c("P017", "P002") ~ 1.09,
      patient == "P005" ~ 1.13,
      patient %in% c("P012", "P019") ~ 2.09,
      TRUE ~ as.numeric(response) + 0.09
    ),
    
    label_y = case_when(
      patient == "P005" ~ median_delta_RF240 - 0.008,
      TRUE ~ median_delta_RF240
    )
  )


A_median <- A_dat %>%
  
  group_by(response) %>%
  
  summarise(
    median_delta_RF240 =
      median(median_delta_RF240),
    .groups = "drop"
  )


p8A <- ggplot(
  A_dat,
  aes(
    x = response,
    y = median_delta_RF240,
    color = response
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  
  geom_jitter(
    width = 0.06,
    height = 0,
    size = 3.2,
    alpha = 0.95
  ) +
  
  geom_crossbar(
    data = A_median,
    aes(
      x = response,
      y = median_delta_RF240,
      ymin = median_delta_RF240,
      ymax = median_delta_RF240
    ),
    width = 0.42,
    linewidth = 0.7,
    color = "black",
    inherit.aes = FALSE
  ) +
  
  geom_text(
    aes(
      x = label_x,
      y = label_y,
      label = patient
    ),
    size = 3,
    hjust = 0,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = response_colors
  ) +
  
  scale_x_discrete(
    expand = expansion(
      mult = c(0.15, 0.28)
    )
  ) +
  
  labs(
    title =
      "A  Zhang: patient-level RF240 remodeling",
    
    subtitle =
      paste0(
        "3 NR vs 2 R patients; complete directional separation; ",
        "exact Wilcoxon p = 0.20"
      ),
    
    x = NULL,
    
    y = "Delta RF240 (patient median)"
  ) +
  
  theme_final +
  
  theme(
    legend.position = "none"
  )


# =========================================================
# 7. FIGURE 8B — final
#
# Use a dedicated title plot rather than nested
# plot_annotation(), so the B heading cannot disappear.
# =========================================================

B_threshold <- B_dat %>%
  
  filter(
    robustness_type ==
      "Clone-support threshold"
  ) %>%
  
  arrange(order) %>%
  
  mutate(
    label_safe =
      paste0(
        ">=",
        str_extract(label, "[0-9]+"),
        " persistent clones"
      ),
    
    label_safe =
      factor(
        label_safe,
        levels = label_safe
      )
  )


B_loo <- B_dat %>%
  
  filter(
    robustness_type ==
      "Leave-one-patient-out"
  ) %>%
  
  arrange(order) %>%
  
  mutate(
    label_safe =
      factor(
        label,
        levels = label
      )
  )


p8B_title <- ggplot() +
  
  annotate(
    "text",
    x = 0,
    y = 1,
    label = "B  Zhang: directional robustness",
    hjust = 0,
    vjust = 1,
    fontface = "bold",
    size = 4.2
  ) +
  
  annotate(
    "text",
    x = 0,
    y = 0.25,
    label =
      paste0(
        "R - NR remains positive across clone-support thresholds ",
        "and removal of each patient"
      ),
    hjust = 0,
    vjust = 1,
    size = 3.2
  ) +
  
  xlim(0, 1) +
  ylim(0, 1) +
  
  theme_void() +
  
  theme(
    plot.margin =
      margin(
        0,
        5,
        0,
        5
      )
  )


p8B1 <- ggplot(
  B_threshold,
  aes(
    x = label_safe,
    y = delta_R_minus_NR,
    group = 1
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  
  geom_line(
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 2.8
  ) +
  
  labs(
    title = "Clone-support sensitivity",
    x = NULL,
    y = "Delta R - NR"
  ) +
  
  theme_final +
  
  theme(
    plot.title =
      element_text(
        size = 10.5,
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 32,
        hjust = 1,
        size = 8
      ),
    
    legend.position = "none"
  )


p8B2 <- ggplot(
  B_loo,
  aes(
    x = label_safe,
    y = delta_R_minus_NR,
    group = 1
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  
  geom_line(
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 2.8
  ) +
  
  labs(
    title = "Leave-one-patient-out",
    x = NULL,
    y = "Delta R - NR"
  ) +
  
  theme_final +
  
  theme(
    plot.title =
      element_text(
        size = 10.5,
        face = "bold"
      ),
    
    axis.text.x =
      element_text(
        angle = 32,
        hjust = 1,
        size = 8
      ),
    
    legend.position = "none"
  )


p8B <- p8B_title /
  p8B1 /
  p8B2 +
  
  plot_layout(
    heights = c(
      0.20,
      1,
      1
    )
  )


# =========================================================
# 8. FIGURE 8C — unchanged analytical content
# =========================================================

C_dat <- C_dat %>%
  
  mutate(
    response =
      factor(
        response,
        levels = c(
          "NR",
          "R"
        )
      ),
    
    comparison_label_safe =
      case_when(
        as.character(comparison_label) ==
          "Base → PD1" ~
          "Base - PD1",
        
        as.character(comparison_label) ==
          "Base → RTPD1" ~
          "Base - RTPD1",
        
        as.character(comparison_label) ==
          "PD1 → RTPD1" ~
          "PD1 - RTPD1",
        
        TRUE ~
          as.character(comparison_label)
      ),
    
    comparison_label_safe =
      factor(
        comparison_label_safe,
        levels = c(
          "Base - PD1",
          "Base - RTPD1",
          "PD1 - RTPD1"
        )
      )
  )


C_median <- C_dat %>%
  
  group_by(
    comparison_label_safe,
    response
  ) %>%
  
  summarise(
    median_delta_RF240 =
      median(median_delta_RF240),
    .groups = "drop"
  )


p8C <- ggplot(
  C_dat,
  aes(
    x = response,
    y = median_delta_RF240,
    color = response
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4
  ) +
  
  geom_jitter(
    width = 0.11,
    height = 0,
    size = 2,
    alpha = 0.78
  ) +
  
  geom_crossbar(
    data = C_median,
    aes(
      x = response,
      y = median_delta_RF240,
      ymin = median_delta_RF240,
      ymax = median_delta_RF240,
      color = response
    ),
    width = 0.40,
    linewidth = 0.65,
    inherit.aes = FALSE
  ) +
  
  facet_wrap(
    ~ comparison_label_safe,
    nrow = 1
  ) +
  
  scale_color_manual(
    values = response_colors
  ) +
  
  labs(
    title =
      "C  Shiao: RF240 remodeling across treatment intervals",
    
    subtitle =
      paste0(
        "No interval shows a reproducible patient-level R vs NR ",
        "difference (all BH q = 0.962)"
      ),
    
    x = NULL,
    
    y = "Delta RF240 (patient median)"
  ) +
  
  theme_final +
  
  theme(
    legend.position = "none"
  )


# =========================================================
# 9. FIGURE 8D — unchanged analytical content
# =========================================================

D_dat <- D_dat %>%
  
  mutate(
    comparison_label_safe =
      case_when(
        comparison_label ==
          "Zhang: Pre → Post" ~
          "Zhang: Pre - Post",
        
        comparison_label ==
          "Shiao: Base → PD1" ~
          "Shiao: Base - PD1",
        
        comparison_label ==
          "Shiao: Base → RTPD1" ~
          "Shiao: Base - RTPD1",
        
        comparison_label ==
          "Shiao: PD1 → RTPD1" ~
          "Shiao: PD1 - RTPD1",
        
        TRUE ~
          comparison_label
      ),
    
    comparison_label_safe =
      factor(
        comparison_label_safe,
        levels = rev(
          c(
            "Zhang: Pre - Post",
            "Shiao: Base - PD1",
            "Shiao: Base - RTPD1",
            "Shiao: PD1 - RTPD1"
          )
        )
      ),
    
    dataset =
      factor(
        dataset,
        levels = c(
          "Zhang",
          "Shiao"
        )
      )
  )


p8D <- ggplot(
  D_dat,
  aes(
    x = comparison_label_safe,
    y = observed_difference,
    color = dataset
  )
) +
  
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45
  ) +
  
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.18,
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3
  ) +
  
  coord_flip() +
  
  scale_color_manual(
    values = c(
      "Zhang" = "black",
      "Shiao" = "grey45"
    )
  ) +
  
  labs(
    title =
      "D  Cross-cohort response-effect summary",
    
    subtitle =
      paste0(
        "Patient-level median R - NR difference with ",
        "percentile bootstrap intervals"
      ),
    
    x = NULL,
    
    y = "Delta RF240 effect (R - NR)",
    
    caption =
      paste0(
        "Zhang bootstrap interval is descriptive only ",
        "(2 R and 3 NR patients)."
      )
  ) +
  
  theme_final +
  
  theme(
    legend.position = "none"
  )


# =========================================================
# 10. FINAL FIGURE 8
# =========================================================

Figure8_v3 <- (
  p8A |
    p8B
) /
  (
    p8C |
      p8D
  ) +
  
  plot_layout(
    widths = c(
      1,
      1.25
    ),
    heights = c(
      1,
      1
    )
  ) +
  
  plot_annotation(
    title =
      paste0(
        "Figure 8. Longitudinal breast cancer cohorts reveal ",
        "context-dependent RF240 remodeling in persistent clonotypes"
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
# 11. Export individual panels
# =========================================================

export_plot(
  p8A,
  dir_panels,
  "Figure8A_Zhang_patientlevel_RF240",
  width = 5.5,
  height = 4.4
)

export_plot(
  p8B,
  dir_panels,
  "Figure8B_Zhang_robustness",
  width = 6.3,
  height = 6.4
)

export_plot(
  p8C,
  dir_panels,
  "Figure8C_Shiao_patientlevel_longitudinal_RF240",
  width = 9.3,
  height = 4.5
)

export_plot(
  p8D,
  dir_panels,
  "Figure8D_crosscohort_effect_summary",
  width = 7.0,
  height = 4.8
)


# =========================================================
# 12. Export integrated Figure 8
# =========================================================

export_plot(
  Figure8_v3,
  dir_out,
  "Figure8_FINAL_integrated",
  width = 13.8,
  height = 10.4
)


# =========================================================
# 13. Verify 15 exports
# =========================================================

expected_stubs <- c(
  file.path(
    dir_panels,
    "Figure8A_Zhang_patientlevel_RF240"
  ),
  file.path(
    dir_panels,
    "Figure8B_Zhang_robustness"
  ),
  file.path(
    dir_panels,
    "Figure8C_Shiao_patientlevel_longitudinal_RF240"
  ),
  file.path(
    dir_panels,
    "Figure8D_crosscohort_effect_summary"
  ),
  file.path(
    dir_out,
    "Figure8_FINAL_integrated"
  )
)


expected_files <- expand_grid(
  stub = expected_stubs,
  extension = c(
    "pdf",
    "png",
    "tiff"
  )
) %>%
  
  mutate(
    file =
      paste0(
        stub,
        ".",
        extension
      ),
    
    exists =
      file.exists(file),
    
    size_bytes =
      ifelse(
        exists,
        file.info(file)$size,
        NA_real_
      )
  )


# =========================================================
# 14. QC
# =========================================================

qc_checks <- tibble(
  
  check = c(
    "5 Zhang patients retained",
    "3 NR retained",
    "2 R retained",
    "All Zhang robustness effects positive",
    "3 Shiao intervals retained",
    "4 cross-cohort effects retained",
    "All 15 PDF/PNG/TIFF exports exist",
    "All exports are non-empty"
  ),
  
  pass = c(
    nrow(A_dat) == 5,
    sum(A_dat$response == "NR") == 3,
    sum(A_dat$response == "R") == 2,
    all(B_dat$delta_R_minus_NR > 0),
    n_distinct(C_dat$comparison_label_safe) == 3,
    nrow(D_dat) == 4,
    all(expected_files$exists),
    all(
      expected_files$size_bytes > 0,
      na.rm = TRUE
    )
  )
)


write_csv(
  expected_files,
  file.path(
    dir_out,
    "14_Figure8_export_verification.csv"
  )
)

write_csv(
  qc_checks,
  file.path(
    dir_out,
    "14_Figure8_QC_checks.csv"
  )
)


# =========================================================
# 15. Console output
# =========================================================

cat("\n\n=========================================================\n")
cat("14 DONE - FIGURE 8 CROSS-COHORT RF240 REMODELING\n")
cat("=========================================================\n")


cat("\n[EXPORT VERIFICATION]\n")

print(
  expected_files,
  n = Inf,
  width = Inf
)


cat("\n[QC CHECKS]\n")

print(
  qc_checks,
  n = Inf,
  width = Inf
)


cat("\n=========================================================\n")


if (
  all(qc_checks$pass)
) {
  
  cat("SUCCESS: ALL PUBLIC FIGURE 8 QC CHECKS PASSED\n")
  
  cat(
    paste0(
      "SUCCESS: ",
      sum(expected_files$exists),
      "/",
      nrow(expected_files),
      " EXPECTED FIGURE 8 EXPORTS EXIST\n"
    )
  )
  
} else {
  
  cat("WARNING: ONE OR MORE PUBLIC FIGURE 8 QC CHECKS FAILED\n")
}


cat("=========================================================\n")
