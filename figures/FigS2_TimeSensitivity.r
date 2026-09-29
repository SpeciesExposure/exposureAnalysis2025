# FigS2_TimeSensitivity.r
# 3-panel vertical figure comparing Fig 1A histogram across baseline end-years:
#   Panel A – V8       (baseline ends 2022, the main run)
#   Panel B – V8_2015  (baseline ends 2015)
#   Panel C – V8_2009  (baseline ends 2009)
#
# The summary data for each panel is cached as an RDS file so the figure
# can be rebuilt quickly without re-loading the large qs sources.
#
# Each baseline is a separate run of 5_CalcExposure.r and 6_MakeDF.r. See
# 1_Setup.r for the run settings.
#
# Run from the repository root:
#   source("figures/FigS2_TimeSensitivity.r")
# ============================================================================

if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))

# ---- Libraries and paths of the main run ----
RUN_VERSION <- "V8"
RUN_NAME    <- "V8"
END_HISTORICAL_YEAR <- 2022L
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

# ============================================================================
# Output locations
# ============================================================================
out_dir   <- file.path(file.path(PROJECT_PATHS$output_v8, "figures", "main"))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cache_dir <- file.path(file.path(PROJECT_PATHS$output_v8, "figures", "main"), "FigS2_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

out_fig   <- file.path(out_dir, "FigS2_TimeSensitivity_v8.pdf")

# ============================================================================
# Colour / label helpers (shared across panels)
# ============================================================================
VAR_LABELS <- c(
  "temp__3__max_up"   = "High temperature (warmest 3 months)",
  "temp__12_up"       = "High temperature (annual mean)",
  "temp__3__min_lo"   = "Low temperature (coldest 3 months)",
  "temp__12_lo"       = "Low temperature (annual mean)",
  "precip__3__min_lo" = "Low precipitation (driest 3 months)",
  "precip__3__max_up" = "High precipitation (wettest 3 months)",
  "precip__12_up"     = "High precipitation (annual total)",
  "precip__12_lo"     = "Low precipitation (annual total)"
)

VAR_COLS <- c(
  "temp__3__max_up"   = "#d73027",
  "temp__12_up"       = "#fc8d59",
  "temp__3__min_lo"   = "#bdbdbd",
  "temp__12_lo"       = "#737373",
  "precip__3__max_up" = "#4575b4",
  "precip__12_up"     = "#74add1",
  "precip__3__min_lo" = "#fee090",
  "precip__12_lo"     = "#fdbb63"
)

# ============================================================================
# Smart reader: tries qs2 first, falls back to qs (older sensitivity runs)
# ============================================================================
read_qs_any <- function(path) {
  tryCatch(
    qs_read(path),
    error = function(e) {
      if (grepl("qs-legacy|QS format not detected", conditionMessage(e), ignore.case = TRUE)) {
        qs::qread(path)
      } else {
        stop(e)
      }
    }
  )
}

# ============================================================================
# Helper: load ex11, deduplicate ties, summarise → year × var counts
# and cache the summary as RDS.
# ============================================================================
get_hist_df <- function(run_name, int_subdir = "Int_V8") {
  cache_file <- file.path(cache_dir, paste0("hist_df_", run_name, "_v1.rds"))

  if (file.exists(cache_file)) {
    message("Reading cached histogram data for ", run_name)
    return(readRDS(cache_file))
  }

  ex11_path <- file.path(
    PROJECT_PATHS$output_base,
    tolower(run_name),
    int_subdir,
    "singleMaxExposure_v1.qs"
  )
  if (!file.exists(ex11_path)) {
    stop("singleMaxExposure_v1.qs not found for ", run_name, ": ", ex11_path)
  }

  message("Loading ex11 for ", run_name, " ...")
  ex11 <- read_qs_any(ex11_path)

  # Deduplicate species-year ties (keep highest pExp, then first var alphabetically)
  ex11 <- ex11 %>%
    dplyr::arrange(spName, year, dplyr::desc(pExp), var) %>%
    dplyr::group_by(spName, year) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  df <- ex11 %>%
    dplyr::group_by(year, var) %>%
    dplyr::summarise(n_species = dplyr::n_distinct(spName), .groups = "drop")

  message("Saving cached histogram data for ", run_name)
  saveRDS(df, cache_file)
  df
}

# ============================================================================
# Helper: build the stacked-bar histogram panel (replicates Fig 1A style)
# ============================================================================
make_hist_panel <- function(df, panel_label, subtitle = NULL, show_legend = TRUE) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = year, y = n_species, fill = var)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_x_continuous(breaks = sort(unique(df$year))) +
    ggplot2::scale_fill_manual(
      labels   = VAR_LABELS,
      values   = VAR_COLS,
      na.value = "grey70"
    ) +
    ggplot2::labs(
      x        = "Year",
      y        = "Number of species",
      fill     = NULL,
      title    = panel_label,
      subtitle = subtitle
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x          = ggplot2::element_text(size = 9, color = "black",
                                                    angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y          = ggplot2::element_text(size = 9, color = "black"),
      axis.title           = ggplot2::element_text(size = 11, color = "black"),
      legend.text          = ggplot2::element_text(size = 8),
      legend.key.size      = ggplot2::unit(0.45, "cm"),
      legend.position      = if (show_legend) c(0.02, 0.98) else "none",
      legend.justification = c(0, 1),
      panel.border         = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.5),
      axis.line            = ggplot2::element_blank(),
      plot.title           = ggplot2::element_text(size = 12, face = "bold"),
      plot.subtitle        = ggplot2::element_text(size = 9, color = "grey30")
    )
  p
}

# ============================================================================
# Load / cache summary data for each run
# ============================================================================
df_v8      <- get_hist_df("V8",      int_subdir = "Int_V8")
df_v8_2015 <- get_hist_df("V8_2015", int_subdir = "Int_V8")
df_v8_2009 <- get_hist_df("V8_2009", int_subdir = "Int_V8")

# Align y-axis range across all panels for comparability
y_max <- max(
  tapply(df_v8$n_species,      df_v8$year,      sum, na.rm = TRUE),
  tapply(df_v8_2015$n_species, df_v8_2015$year, sum, na.rm = TRUE),
  tapply(df_v8_2009$n_species, df_v8_2009$year, sum, na.rm = TRUE)
)
common_y <- ggplot2::scale_y_continuous(limits = c(0, y_max * 1.05))

# ============================================================================
# Build panels
# ============================================================================
pA <- make_hist_panel(df_v8,
  panel_label = "(a) V8 \u2013 baseline ends 2022 (canonical)",
  show_legend = TRUE
) + common_y

pB <- make_hist_panel(df_v8_2015,
  panel_label = "(b) V8_2015 \u2013 baseline ends 2015",
  show_legend = FALSE
) + common_y

pC <- make_hist_panel(df_v8_2009,
  panel_label = "(c) V8_2009 \u2013 baseline ends 2009",
  show_legend = FALSE
) + common_y

# ============================================================================
# Combine and save
# ============================================================================
multi_fig <- gridExtra::arrangeGrob(pA, pB, pC, ncol = 1)

ggplot2::ggsave(
  filename = out_fig,
  plot     = multi_fig,
  width    = 280,
  height   = 330,
  units    = "mm"
)

cat("Wrote FigS2 to:", out_fig, "\n")

