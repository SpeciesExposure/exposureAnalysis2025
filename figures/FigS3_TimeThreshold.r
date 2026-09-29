# FigS3_TimeThreshold.r
# Sensitivity of the main timeline (dominant climate threat per species-year)
# to quantile definition (year.q / cell.q thresholds).
#
# Panel A: baseline  (year.q=0.99, cell.q=0.99)
# Panel B: stricter annual   (year.q=1.00, cell.q=0.99)
# Panel C: stricter cell     (year.q=0.99, cell.q=1.00)
# Panel D: looser cell       (year.q=0.99, cell.q=0.95)
#
# Run from the repository root:  Rscript figures/FigS3_TimeThreshold.r

# ── 1. Bootstrap ─────────────────────────────────────────────────────────────
if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))

if (!exists("RUN_VERSION"))         RUN_VERSION         <- "V8"
if (!exists("RUN_NAME"))            RUN_NAME            <- "V8"
if (!exists("END_HISTORICAL_YEAR")) END_HISTORICAL_YEAR <- 2022L

source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))
# dataDir, intDir, plotDir, mc.cores now defined.

# ── 2. Scenarios ──────────────────────────────────────────────────────────────
scenarios <- data.frame(
  scenario_tag = c("q099c099", "q100c099", "q099c100", "q099c095"),
  scenario_label = c(
    "Panel A (Main text Figure 1 baseline): year.q=0.99, cell.q=0.99",
    "Panel B: year.q=1.00, cell.q=0.99",
    "Panel C: year.q=0.99, cell.q=1.00",
    "Panel D: year.q=0.99, cell.q=0.95"
  ),
  year_q_up = c(0.99, 1.00, 0.99, 0.99),
  cell_q_up = c(0.99, 0.99, 1.00, 0.95),
  stringsAsFactors = FALSE
)

rebuild_ex11 <- FALSE

# ── 3. Build ex11 for one threshold scenario ──────────────────────────────────
build_ex11 <- function(year_q_up, cell_q_up, scenario_tag, rebuild = FALSE) {
  year_q_lo <- 1 - year_q_up
  cell_q_lo <- 1 - cell_q_up

  cache_dir  <- file.path(intDir, "sensitivity", "timeline_quantiles")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  cache_file <- file.path(cache_dir,
                           paste0("singleMaxExposure_", scenario_tag, "_v1.qs"))

  if (file.exists(cache_file) && !rebuild) {
    message("  Loading cached ex11 for ", scenario_tag)
    return(qs2::qs_read(cache_file))
  }

  message("  Building ex11 for ", scenario_tag,
          "  (year.up=", year_q_up, ", cell.up=", cell_q_up, ")")

  vars <- list.files(file.path(dataDir, "variable_outputs"))

  ex1 <- lapply(vars, function(vv) {
    stats_dir <- file.path(dataDir, "variable_outputs", vv,
                            paste0(vv, "_stats"))
    ff <- list.files(stats_dir, full.names = TRUE)
    if (length(ff) == 0) return(NULL)

    this_year_q <- if (grepl("up", vv)) year_q_up else year_q_lo
    this_cell_q <- if (grepl("up", vv)) cell_q_up else cell_q_lo

    parallel::mclapply(ff, function(x) {
      tryCatch({
        # Use qs::qread explicitly—stats files are in legacy qs1 format
        d <- qs::qread(x)
        d <- d[dplyr::near(d$year.q,      this_year_q) &
               d$isRounded == FALSE       &
               dplyr::near(d$cell.q,      this_cell_q) &
               d$quantileType == 7, ]
        if (nrow(d) == 0) return(NULL)

        d <- dplyr::select(d, -dplyr::contains("mar"),
                              -dplyr::contains("XnM"),
                              -expPre1990,
                              -dplyr::any_of("year.cell.q"),
                              -quantileIndex)
        year_cols <- grep("^X[0-9]{4}$", names(d), value = TRUE)
        year_cols <- year_cols[year_cols != "X2026"]

        d[year_cols] <- lapply(d[year_cols],
                               function(z) round(z / d$rangeSize, 3))
        d$spName <- tools::file_path_sans_ext(basename(x))
        d$var    <- vv

        long <- tidyr::pivot_longer(d,
                                    cols       = dplyr::all_of(year_cols),
                                    names_to   = "year",
                                    values_to  = "pExp")
        long <- long[long$pExp >= 0.25, ]
        long$year <- as.numeric(sub("X", "", long$year))
        long
      }, error = function(e) {
        message("    SKIP ", basename(x), ": ", conditionMessage(e))
        NULL
      })
    }, mc.cores = mc.cores) |> dplyr::bind_rows()
  }) |> dplyr::bind_rows()

  if (nrow(ex1) == 0) {
    warning("No rows after filtering for scenario ", scenario_tag,
            ". Check year.q / cell.q values.")
    return(ex1)
  }

  ex11 <- ex1 |>
    dplyr::arrange(spName, year, dplyr::desc(pExp), var) |>
    dplyr::group_by(spName, year) |>
    dplyr::slice(1) |>
    dplyr::ungroup()

  qs2::qs_save(ex11, cache_file)
  ex11
}

# ── 4. Make one panel plot ────────────────────────────────────────────────────
make_panel <- function(ex11, panel_label, y_max, guide_breaks) {
  guide_breaks <- sort(unique(guide_breaks[guide_breaks >= 0 &
                                             guide_breaks <= y_max]))
  var_labels <- c(
    "temp__3__max_up"   = "High temperature (warmest 3 months)",
    "temp__12_up"       = "High temperature (annual mean)",
    "temp__3__min_lo"   = "Low temperature (coldest 3 months)",
    "temp__12_lo"       = "Low temperature (annual mean)",
    "precip__3__max_up" = "High precipitation (wettest 3 months)",
    "precip__12_up"     = "High precipitation (annual total)",
    "precip__3__min_lo" = "Low precipitation (driest 3 months)",
    "precip__12_lo"     = "Low precipitation (annual total)"
  )
  var_colors <- c(
    "temp__3__max_up"   = "#d73027",
    "temp__12_up"       = "#fc8d59",
    "temp__3__min_lo"   = "#bdbdbd",
    "temp__12_lo"       = "#737373",
    "precip__3__max_up" = "#4575b4",
    "precip__12_up"     = "#74add1",
    "precip__3__min_lo" = "#fee090",
    "precip__12_lo"     = "#fdbb63"
  )

  ex11 |>
    dplyr::group_by(year, var) |>
    dplyr::summarise(n_species = dplyr::n_distinct(spName), .groups = "drop") |>
    ggplot2::ggplot(ggplot2::aes(x = year, y = n_species, fill = var)) +
    ggplot2::geom_hline(
      data = data.frame(yint = guide_breaks[guide_breaks > 0]),
      ggplot2::aes(yintercept = yint),
      inherit.aes = FALSE,
      linetype = "dashed", linewidth = 0.35, color = "grey55"
    ) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::scale_x_continuous(breaks = sort(unique(ex11$year))) +
    ggplot2::scale_y_continuous(limits = c(0, y_max), breaks = guide_breaks) +
    ggplot2::scale_fill_manual(values = var_colors, labels = var_labels) +
    ggplot2::labs(x = "Year", y = "Number of species",
                  fill = NULL, title = panel_label) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x  = ggplot2::element_text(size = 11, color = "black",
                                            angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y  = ggplot2::element_text(size = 11, color = "black"),
      axis.title   = ggplot2::element_text(size = 13, color = "black"),
      legend.text  = ggplot2::element_text(size = 10),
      legend.key.size = grid::unit(0.5, "cm"),
      legend.position = c(0.01, 0.98),
      legend.justification = c(0, 1),
      legend.background = ggplot2::element_rect(
        fill  = scales::alpha("white", 0.85),
        color = "grey60", linewidth = 0.25),
      panel.border = ggplot2::element_rect(color = "black", fill = NA,
                                            linewidth = 0.4),
      axis.line    = ggplot2::element_blank(),
      plot.title   = ggplot2::element_text(size = 12, face = "bold",
                                            hjust = 0.5)
    )
}

# ── 5. Run all scenarios ──────────────────────────────────────────────────────
message("=== Building ex11 data for all threshold scenarios ===")
ex11_list <- vector("list", nrow(scenarios))
for (i in seq_len(nrow(scenarios))) {
  message(" - ", scenarios$scenario_tag[i])
  ex11_list[[i]] <- build_ex11(
    year_q_up    = scenarios$year_q_up[i],
    cell_q_up    = scenarios$cell_q_up[i],
    scenario_tag = scenarios$scenario_tag[i],
    rebuild      = rebuild_ex11
  )
}

global_y_max <- max(vapply(ex11_list, function(d) {
  if (nrow(d) == 0) return(0)
  d |> dplyr::distinct(spName, year) |> dplyr::count(year) |>
    dplyr::summarise(m = max(n, na.rm = TRUE)) |> dplyr::pull(m)
}, numeric(1)), na.rm = TRUE)

common_breaks <- pretty(c(0, global_y_max), n = 6)

# ── 6. Build plots ────────────────────────────────────────────────────────────
message("=== Generating plots ===")
plot_list <- lapply(seq_len(nrow(scenarios)), function(i) {
  make_panel(ex11_list[[i]], scenarios$scenario_label[i],
             global_y_max, common_breaks)
})

combined <- gridExtra::arrangeGrob(grobs = plot_list, ncol = 1, nrow = 4)

# ── 7. Save ───────────────────────────────────────────────────────────────────
dir.create(plotDir, recursive = TRUE, showWarnings = FALSE)
ver_tag  <- tolower(gsub("[^A-Za-z0-9_]", "_", RUN_NAME))
pdf_file <- file.path(plotDir, paste0("FigS3_TimeThreshold_", ver_tag, ".pdf"))
png_file <- file.path(plotDir, paste0("FigS3_TimeThreshold_", ver_tag, ".png"))

ggplot2::ggsave(pdf_file, combined, width = 300, height = 380, units = "mm")
ggplot2::ggsave(png_file, combined, width = 300, height = 380, units = "mm",
                dpi = 300)

message("Saved PDF: ", pdf_file)
message("Saved PNG: ", png_file)
system(paste0("open ", pdf_file))
