# =========================================================
# 00_setup.R
#
# Common configuration for the RF240 TNBC code compendium
# =========================================================


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
# 2. Project root
#
# Detect repository root from the location of this script.
#
# Expected location:
#   RF240_TNBC/scripts/00_setup.R
#
# IMPORTANT:
# This script should be run with source(), not line by line.
# =========================================================

setup_file <- tryCatch(
  normalizePath(
    sys.frame(1)$ofile,
    winslash = "/",
    mustWork = TRUE
  ),
  error = function(e) NA_character_
)


if (is.na(setup_file)) {
  
  stop(
    paste0(
      "Could not determine the location of 00_setup.R.\n\n",
      "Please run it using source(), for example:\n",
      'source("~/Documents/R studio/RF240_TNBC/scripts/00_setup.R")\n\n',
      "Do not execute 00_setup.R line by line."
    )
  )
}


dir_scripts_detected <- dirname(
  setup_file
)


project_root <- normalizePath(
  file.path(
    dir_scripts_detected,
    ".."
  ),
  winslash = "/",
  mustWork = TRUE
)


# Safety check:
# the repository root must contain the expected directories.

expected_root_dirs <- c(
  "scripts",
  "metadata",
  "data",
  "results",
  "documentation"
)


root_check <- file.path(
  project_root,
  expected_root_dirs
)


if (!all(dir.exists(root_check))) {
  
  stop(
    paste0(
      "00_setup.R was found, but its parent directory does not ",
      "look like the RF240_TNBC repository.\n\n",
      "Detected project root:\n",
      project_root,
      "\n"
    )
  )
}


# =========================================================
# 3. Repository directories
# =========================================================

dir_scripts <- file.path(
  project_root,
  "scripts"
)

dir_metadata <- file.path(
  project_root,
  "metadata"
)

dir_data <- file.path(
  project_root,
  "data"
)

dir_data_baseline <- file.path(
  dir_data,
  "baseline"
)

dir_data_longitudinal <- file.path(
  dir_data,
  "longitudinal"
)

dir_data_tcr <- file.path(
  dir_data,
  "tcr"
)

dir_data_bassez <- file.path(
  dir_data,
  "bassez"
)

dir_data_zhang <- file.path(
  dir_data,
  "zhang"
)

dir_data_shiao <- file.path(
  dir_data,
  "shiao"
)

dir_results <- file.path(
  project_root,
  "results"
)

dir_figures <- file.path(
  dir_results,
  "figures"
)

dir_supp_figures <- file.path(
  dir_results,
  "supplementary_figures"
)

dir_tables <- file.path(
  dir_results,
  "tables"
)

dir_documentation <- file.path(
  project_root,
  "documentation"
)


# =========================================================
# 4. Confirm required repository structure
# =========================================================

required_dirs <- c(
  dir_scripts,
  dir_metadata,
  dir_data,
  dir_data_baseline,
  dir_data_longitudinal,
  dir_data_tcr,
  dir_data_bassez,
  dir_data_zhang,
  dir_data_shiao,
  dir_results,
  dir_figures,
  dir_supp_figures,
  dir_tables,
  dir_documentation
)


missing_dirs <- required_dirs[
  !dir.exists(
    required_dirs
  )
]


if (
  length(
    missing_dirs
  ) > 0
) {
  
  stop(
    paste0(
      "Repository structure is incomplete.\n\n",
      "Missing directories:\n",
      paste(
        missing_dirs,
        collapse = "\n"
      )
    )
  )
}


# =========================================================
# 5. Figure-specific output directories
# =========================================================

figure_dirs <- setNames(
  
  file.path(
    dir_figures,
    paste0(
      "Figure",
      1:8
    )
  ),
  
  paste0(
    "Figure",
    1:8
  )
)


supp_figure_dirs <- setNames(
  
  file.path(
    dir_supp_figures,
    paste0(
      "Supplementary_Figure",
      1:8
    )
  ),
  
  paste0(
    "Supplementary_Figure",
    1:8
  )
)


# =========================================================
# 6. Response conventions
# =========================================================

response_levels <- c(
  "NR",
  "R"
)

response_colors <- c(
  "NR" = "#D55E00",
  "R"  = "#0072B2"
)


# Bassez terminology
expansion_response_levels <- c(
  "NE",
  "E"
)

expansion_response_colors <- c(
  "NE" = "#D55E00",
  "E"  = "#0072B2"
)


# =========================================================
# 7. Longitudinal timepoint conventions
# =========================================================

timepoint_levels <- c(
  "Base",
  "PD1",
  "RTPD1"
)


comparison_levels <- c(
  "Base_vs_PD1",
  "Base_vs_RTPD1",
  "PD1_vs_RTPD1"
)


# =========================================================
# 8. Statistical conventions
# =========================================================

alpha_significant <- 0.05

# Used only in specific analyses where a trend threshold
# is explicitly part of the reported analysis.
alpha_trend <- 0.15


# =========================================================
# 9. General plotting theme
# =========================================================

theme_rf240 <- function(
    base_size = 11
) {
  
  theme_classic(
    base_size = base_size
  ) +
    
    theme(
      
      axis.title =
        element_text(
          size = base_size
        ),
      
      axis.text =
        element_text(
          size = base_size - 2,
          color = "black"
        ),
      
      plot.title =
        element_text(
          size = base_size + 1,
          face = "bold"
        ),
      
      plot.subtitle =
        element_text(
          size = base_size - 1
        ),
      
      strip.background =
        element_blank(),
      
      strip.text =
        element_text(
          face = "bold"
        ),
      
      legend.title =
        element_blank(),
      
      legend.position =
        "top"
    )
}


# =========================================================
# 10. Multi-format figure export
# =========================================================

export_plot <- function(
    plot_object,
    output_directory,
    filename_stub,
    width,
    height,
    dpi = 600
) {
  
  dir_create(
    output_directory
  )
  
  
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
    file = pdf_path,
    width = width,
    height = height,
    useDingbats = FALSE
  )
  
  print(
    plot_object
  )
  
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
  
  
  expected_files <- c(
    pdf_path,
    png_path,
    tiff_path
  )
  
  
  missing_exports <- expected_files[
    !file.exists(
      expected_files
    )
  ]
  
  
  if (
    length(
      missing_exports
    ) > 0
  ) {
    
    stop(
      paste0(
        "Figure export failed.\n\n",
        "Missing files:\n",
        paste(
          missing_exports,
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
# 11. Safe two-group Wilcoxon helper
# =========================================================

safe_wilcox_two_group <- function(
    x,
    y,
    exact = FALSE
) {
  
  x <- x[
    is.finite(
      x
    )
  ]
  
  y <- y[
    is.finite(
      y
    )
  ]
  
  
  if (
    length(
      x
    ) < 2 ||
    length(
      y
    ) < 2
  ) {
    
    return(
      tibble(
        n_x = length(x),
        n_y = length(y),
        statistic = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  wt <- suppressWarnings(
    wilcox.test(
      x,
      y,
      exact = exact
    )
  )
  
  
  tibble(
    
    n_x =
      length(
        x
      ),
    
    n_y =
      length(
        y
      ),
    
    statistic =
      unname(
        wt$statistic
      ),
    
    p_value =
      wt$p.value
  )
}


# =========================================================
# 12. Safe one-sample Wilcoxon helper
# =========================================================

safe_wilcox_onesample <- function(
    x,
    mu = 0,
    exact = FALSE,
    min_n = 3
) {
  
  x <- x[
    is.finite(
      x
    )
  ]
  
  
  if (
    length(
      x
    ) < min_n ||
    n_distinct(
      x
    ) < 2
  ) {
    
    return(
      tibble(
        
        n =
          length(
            x
          ),
        
        median =
          ifelse(
            length(
              x
            ) > 0,
            median(
              x
            ),
            NA_real_
          ),
        
        statistic =
          NA_real_,
        
        p_value =
          NA_real_
      )
    )
  }
  
  
  wt <- suppressWarnings(
    wilcox.test(
      x,
      mu = mu,
      exact = exact
    )
  )
  
  
  tibble(
    
    n =
      length(
        x
      ),
    
    median =
      median(
        x
      ),
    
    statistic =
      unname(
        wt$statistic
      ),
    
    p_value =
      wt$p.value
  )
}


# =========================================================
# 13. Safe Spearman helper
# =========================================================

safe_spearman <- function(
    x,
    y,
    min_n = 5
) {
  
  keep <- is.finite(
    x
  ) &
    is.finite(
      y
    )
  
  
  x <- x[
    keep
  ]
  
  y <- y[
    keep
  ]
  
  
  if (
    length(
      x
    ) < min_n ||
    n_distinct(
      x
    ) < 2 ||
    n_distinct(
      y
    ) < 2
  ) {
    
    return(
      tibble(
        n = length(x),
        rho = NA_real_,
        p_value = NA_real_
      )
    )
  }
  
  
  ct <- suppressWarnings(
    cor.test(
      x,
      y,
      method = "spearman",
      exact = FALSE
    )
  )
  
  
  tibble(
    
    n =
      length(
        x
      ),
    
    rho =
      unname(
        ct$estimate
      ),
    
    p_value =
      ct$p.value
  )
}


# =========================================================
# 14. Percentile bootstrap median CI
# =========================================================

bootstrap_median_ci <- function(
    x,
    n_boot = 5000,
    seed = 240,
    conf = 0.95
) {
  
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
        n = length(x),
        median = ifelse(
          length(x) > 0,
          median(x),
          NA_real_
        ),
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }
  
  
  set.seed(
    seed
  )
  
  
  boot <- replicate(
    
    n_boot,
    
    median(
      sample(
        x,
        size = length(x),
        replace = TRUE
      )
    )
  )
  
  
  alpha <- (
    1 - conf
  ) / 2
  
  
  tibble(
    
    n =
      length(
        x
      ),
    
    median =
      median(
        x
      ),
    
    ci_low =
      unname(
        quantile(
          boot,
          probs = alpha,
          na.rm = TRUE
        )
      ),
    
    ci_high =
      unname(
        quantile(
          boot,
          probs = 1 - alpha,
          na.rm = TRUE
        )
      )
  )
}


# =========================================================
# 15. Initialization message
# =========================================================

cat(
  "\n=========================================================\n"
)

cat(
  "RF240 TNBC code compendium initialized successfully\n"
)

cat(
  "=========================================================\n"
)

cat(
  "Project root:\n",
  project_root,
  "\n\n"
)

cat(
  "Scripts directory:\n",
  dir_scripts,
  "\n\n"
)

cat(
  "Results directory:\n",
  dir_results,
  "\n"
)

cat(
  "=========================================================\n"
)