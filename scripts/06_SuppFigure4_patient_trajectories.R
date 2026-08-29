# =========================================================
# 06_SuppFigure4_patient_trajectories.R
#
# RF240 TNBC code compendium
#
# Supplementary Figure 4:
# Patient-level longitudinal trajectories of RF240 global
# and module scores.
#
# Audited historical source:
#   38B_RF240_longitudinal_scores_plots.R
#
# Only the trajectory panels used in final Supplementary
# Figure 4 are retained here. Group-level boxplots,
# Wilcoxon tests, and delta analyses are represented by
# the audited Figure 3 workflow and are not duplicated.
# =========================================================


# =========================================================
# 1. Packages
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})


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

dir_output <- file.path(
  project_root,
  "results",
  "supplementary_figures",
  "Supplementary_Figure4"
)

dir_tables <- file.path(
  project_root,
  "results",
  "tables",
  "Supplementary_Figure4"
)

dir.create(
  dir_output,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  dir_tables,
  recursive = TRUE,
  showWarnings = FALSE
)

file_global <- file.path(
  dir_inputs,
  "RF240_longitudinal_global_score_patient_time.csv"
)

file_modules <- file.path(
  dir_inputs,
  "RF240_longitudinal_module_scores_patient_time.csv"
)

required_inputs <- c(
  file_global,
  file_modules
)

if (any(!file.exists(required_inputs))) {
  stop(
    "Missing Supplementary Figure 4 input tables."
  )
}


# =========================================================
# 3. Read audited longitudinal checkpoints
# =========================================================

global_score <- readr::read_csv(
  file_global,
  show_col_types = FALSE
)

module_scores <- readr::read_csv(
  file_modules,
  show_col_types = FALSE
)

time_levels <- c(
  "Base",
  "PD1",
  "RTPD1"
)

pal_response <- c(
  "NR" = "#D55E00",
  "R" = "#0072B2"
)

global_score <- global_score %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = time_levels
    ),
    response_pcr = factor(
      response_pcr,
      levels = c("NR", "R")
    )
  )

module_scores <- module_scores %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = time_levels
    ),
    response_pcr = factor(
      response_pcr,
      levels = c("NR", "R")
    ),
    module = factor(module)
  )


# =========================================================
# 4. Panel A: global RF240 patient trajectories
# =========================================================

p_supp4a <- ggplot(
  global_score,
  aes(
    x = timepoint,
    y = RF240_global_score,
    group = patient_id,
    color = response_pcr
  )
) +
  geom_line(
    alpha = 0.45,
    linewidth = 0.6
  ) +
  geom_point(
    size = 2,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = pal_response
  ) +
  theme_classic(
    base_size = 13
  ) +
  labs(
    title = "Patient-level RF240 longitudinal trajectories",
    x = NULL,
    y = "Mean RF240 expression",
    color = "Response"
  )


# =========================================================
# 5. Panel B: module-level patient trajectories
# =========================================================

p_supp4b <- ggplot(
  module_scores,
  aes(
    x = timepoint,
    y = module_score,
    group = patient_id,
    color = response_pcr
  )
) +
  geom_line(
    alpha = 0.35,
    linewidth = 0.5
  ) +
  geom_point(
    size = 1.5,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ module,
    scales = "free_y",
    ncol = 2
  ) +
  scale_color_manual(
    values = pal_response
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Patient-level RF240 module trajectories",
    x = NULL,
    y = "Module score",
    color = "Response"
  )


# =========================================================
# 6. Export helper
# =========================================================

save_panel <- function(
    plot_object,
    filename_stub,
    width,
    height,
    dpi = 600
) {
  ggplot2::ggsave(
    filename = file.path(
      dir_output,
      paste0(filename_stub, ".pdf")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in"
  )

  ggplot2::ggsave(
    filename = file.path(
      dir_output,
      paste0(filename_stub, ".png")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi
  )

  ggplot2::ggsave(
    filename = file.path(
      dir_output,
      paste0(filename_stub, ".tiff")
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    compression = "lzw"
  )
}


# =========================================================
# 7. Export panels
# =========================================================

save_panel(
  p_supp4a,
  "SupplementaryFigure4A_global_patient_trajectories",
  width = 7,
  height = 5
)

save_panel(
  p_supp4b,
  "SupplementaryFigure4B_module_patient_trajectories",
  width = 10,
  height = 8
)


# =========================================================
# 8. Support tables and session information
# =========================================================

readr::write_csv(
  global_score,
  file.path(
    dir_tables,
    "SupplementaryFigure4A_global_patient_trajectory_data.csv"
  )
)

readr::write_csv(
  module_scores,
  file.path(
    dir_tables,
    "SupplementaryFigure4B_module_patient_trajectory_data.csv"
  )
)

capture.output(
  sessionInfo(),
  file = file.path(
    dir_tables,
    "SupplementaryFigure4_sessionInfo.txt"
  )
)


# =========================================================
# 9. Summary
# =========================================================

cat("\n=========================================================\n")
cat("06 DONE - Supplementary Figure 4 patient trajectories\n")
cat("=========================================================\n")
cat("Global patient-time rows:", nrow(global_score), "\n")
cat("Module rows:", nrow(module_scores), "\n")
cat("Output directory:\n", dir_output, "\n")
cat("=========================================================\n")
