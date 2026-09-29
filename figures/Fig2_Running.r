# Fig2_Running.r
# Figure 2: repeated exposure through time.
# Panel A: cumulative pool of species ever exposed, by number of exposure years.
# Panel B: species exposed in each year, by number of exposure years so far.
# Panel C: first-time exposures by variable.
#
# Inputs:  singleMaxExposure_v1.qs
# Outputs: figure in plotDir

if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

ex11_file <- file.path(intDir, 'singleMaxExposure_v1.qs')
if (!file.exists(ex11_file)) {
  stop(
    "Missing single-max exposure file for RUN_VERSION ", RUN_VERSION, ": ", ex11_file
  )
}

ex11_raw <- qs_read(ex11_file)

ex11 <- ex11_raw %>%
  arrange(spName, year, desc(pExp), var) %>%
  group_by(spName, year) %>%
  slice(1) %>%
  ungroup()

running_count_cap <- if (exists("running_count_cap_override")) running_count_cap_override else 6
overflow_label <- if (is.finite(running_count_cap)) paste0(">", as.integer(running_count_cap) - 1) else NA_character_

year_levels <- sort(unique(ex11$year))
species_levels <- sort(unique(ex11$spName))

species_year_cum <- tidyr::expand_grid(
  spName = species_levels,
  year = year_levels
) %>%
  left_join(
    ex11 %>% distinct(spName, year) %>% mutate(exposed_this_year = 1L),
    by = c("spName", "year")
  ) %>%
  mutate(exposed_this_year = replace_na(exposed_this_year, 0L)) %>%
  group_by(spName) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(exposure_count = cumsum(exposed_this_year)) %>%
  ungroup() %>%
  filter(exposure_count > 0)

max_count <- max(species_year_cum$exposure_count, na.rm = TRUE)

if (is.finite(running_count_cap)) {
  if (running_count_cap <= 1) {
    bin_levels <- overflow_label
  } else if (max_count < running_count_cap) {
    bin_levels <- as.character(seq_len(max_count))
  } else {
    bin_levels <- c(as.character(seq_len(as.integer(running_count_cap) - 1)), overflow_label)
  }
} else {
  bin_levels <- as.character(seq_len(max_count))
}

message("running_count_cap = ", running_count_cap)
message("running bins = ", paste(bin_levels, collapse = ", "))

base_cols_41 <- c(
  "#d73027", "#fc8d59", "#fdbb63", "#fee090",
  "#74add1", "#4575b4", "#313695", "#542788"
)
if (length(bin_levels) <= length(base_cols_41)) {
  bin_cols <- base_cols_41[seq_along(bin_levels)]
} else {
  bin_cols <- grDevices::colorRampPalette(base_cols_41)(length(bin_levels))
}
bin_cols <- rev(bin_cols)
bin_col_map <- stats::setNames(bin_cols, bin_levels)

running_df <- species_year_cum %>%
  mutate(
    exposure_count_bin = ifelse(
      is.finite(running_count_cap) & exposure_count >= running_count_cap,
      overflow_label,
      as.character(exposure_count)
    )
  ) %>%
  count(year, exposure_count_bin, name = "n_species") %>%
  tidyr::complete(year = year_levels, exposure_count_bin = bin_levels, fill = list(n_species = 0)) %>%
  mutate(exposure_count_bin = factor(exposure_count_bin, levels = bin_levels))

annual_df <- species_year_cum %>%
  filter(exposed_this_year == 1L) %>%
  mutate(
    exposure_count_bin = ifelse(
      is.finite(running_count_cap) & exposure_count >= running_count_cap,
      overflow_label,
      as.character(exposure_count)
    )
  ) %>%
  count(year, exposure_count_bin, name = "n_species") %>%
  tidyr::complete(year = year_levels, exposure_count_bin = bin_levels, fill = list(n_species = 0)) %>%
  mutate(exposure_count_bin = factor(exposure_count_bin, levels = bin_levels))

running_totals <- running_df %>%
  group_by(year) %>%
  summarise(total_species = sum(n_species), .groups = "drop")

new_threat_first_year <- ex11 %>%
  distinct(spName, var, year) %>%
  group_by(spName, var) %>%
  summarise(first_year = min(year), .groups = "drop")

var_label_map <- c(
  "temp__3__max_up" = "High temperature (warmest 3 months)",
  "temp__12_up" = "High temperature (annual mean)",
  "temp__3__min_lo" = "Low temperature (coldest 3 months)",
  "temp__12_lo" = "Low temperature (annual mean)",
  "precip__3__min_lo" = "Low precipitation (driest 3 months)",
  "precip__3__max_up" = "High precipitation (wettest 3 months)",
  "precip__12_up" = "High precipitation (annual total)",
  "precip__12_lo" = "Low precipitation (annual total)"
)

var_color_map_full <- c(
  "temp__3__max_up" = "#d73027",
  "temp__12_up" = "#fc8d59",
  "temp__3__min_lo" = "#bdbdbd",
  "temp__12_lo" = "#737373",
  "precip__3__max_up" = "#4575b4",
  "precip__12_up" = "#74add1",
  "precip__3__min_lo" = "#fee090",
  "precip__12_lo" = "#fdbb63"
)

var_levels <- names(var_color_map_full)[names(var_color_map_full) %in% unique(new_threat_first_year$var)]
if (length(var_levels) == 0) {
  var_levels <- sort(unique(new_threat_first_year$var))
}

new_threat_var_labels <- stats::setNames(
  vapply(var_levels, function(vv) {
    if (vv %in% names(var_label_map)) var_label_map[[vv]] else vv
  }, character(1)),
  var_levels
)

base_var_colors <- vapply(var_levels, function(vv) {
  if (vv %in% names(var_color_map_full)) var_color_map_full[[vv]] else NA_character_
}, character(1))

if (all(!is.na(base_var_colors))) {
  new_threat_var_colors <- stats::setNames(base_var_colors, var_levels)
} else {
  fallback_colors <- scales::hue_pal()(length(var_levels))
  base_var_colors[is.na(base_var_colors)] <- fallback_colors[is.na(base_var_colors)]
  new_threat_var_colors <- stats::setNames(base_var_colors, var_levels)
}

new_threat_df <- new_threat_first_year %>%
  count(year = first_year, var, name = "n_new_exposures") %>%
  tidyr::complete(year = year_levels, var = var_levels, fill = list(n_new_exposures = 0)) %>%
  group_by(var) %>%
  arrange(year, .by_group = TRUE) %>%
  mutate(n_first_time_exposures_raw = n_new_exposures) %>%
  mutate(cum_exposures = cumsum(n_new_exposures)) %>%
  ungroup() %>%
  mutate(var = factor(var, levels = var_levels))

new_threat_totals <- new_threat_df %>%
  group_by(year) %>%
  summarise(total_exposure_contributions = sum(cum_exposures), .groups = "drop")

if (any(diff(running_totals$total_species) < 0)) {
  warning("Running-total plot check failed: yearly totals are not monotonic.")
} else {
  message("Running-total plot check passed: yearly totals are non-decreasing.")
}

p.running <- ggplot(running_df, aes(x = year, y = n_species, fill = exposure_count_bin)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
  scale_x_continuous(breaks = year_levels) +
  scale_fill_manual(values = bin_col_map, drop = FALSE) +
  labs(
    x = "Year",
    y = NULL,
    fill = "Running exposure\ncount",
    title = "A. All species: cumulative exposures up to year"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.45, "cm"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

p.annual <- ggplot(annual_df, aes(x = year, y = n_species, fill = exposure_count_bin)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
  scale_x_continuous(breaks = year_levels) +
  scale_fill_manual(values = bin_col_map, drop = FALSE) +
  labs(
    x = "Year",
    y = NULL,
    fill = "Running exposure count",
    title = "B. Repeated exposure count"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.45, "cm"),
    legend.position = c(0.02, 0.98),
    legend.justification = c(0, 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

p.new_threat <- ggplot(new_threat_df, aes(x = year, y = n_first_time_exposures_raw, fill = var)) +
  geom_col(width = 0.85, position = position_stack(reverse = TRUE)) +
  scale_x_continuous(breaks = year_levels) +
  scale_fill_manual(values = new_threat_var_colors, labels = new_threat_var_labels, drop = FALSE) +
  labs(
    x = "Year",
    y = NULL,
    fill = "Variable",
    title = "C. Species first-time exposures by variable"
  ) +
  guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "bottom",
    legend.justification = "center",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

# ---- Mark El Niño (red up-triangle) and La Niña (blue down-triangle) years above the bars ----
source(file.path(PROJECT_PATHS$figures_src, "enso_strip.R"))
tops_of <- function(d, ycol) { t <- tapply(d[[ycol]], d$year, sum, na.rm = TRUE); data.frame(year = as.integer(names(t)), top = as.numeric(t)) }
yrs <- min(year_levels):max(year_levels)
p.running    <- p.running    + enso_marks(yrs, tops_of(running_df, "n_species")) + list(
    ggplot2::annotate("point", x = 2000, y = 8300, shape = 24, size = 2.4 * 1.6, fill = "#d73027", colour = "black", stroke = 0.3),
    ggplot2::annotate("text",  x = 2000.6, y = 8300, label = "Strong El Ni\u00f1o", hjust = 0, size = 3.4),
    ggplot2::annotate("point", x = 2007.5, y = 8300, shape = 24, size = 2.4, fill = "#d73027", colour = "black", stroke = 0.3),
    ggplot2::annotate("text",  x = 2008.1, y = 8300, label = "El Ni\u00f1o", hjust = 0, size = 3.4),
    ggplot2::annotate("point", x = 2012.5, y = 8300, shape = 25, size = 2.4, fill = "#2166ac", colour = "black", stroke = 0.3),
    ggplot2::annotate("text",  x = 2013.1, y = 8300, label = "La Ni\u00f1a", hjust = 0, size = 3.4)) + scale_y_continuous(limits = c(0, 8600), expand = expansion(mult = c(0.01, 0.02)))
p.annual     <- p.annual     + enso_marks(yrs, tops_of(annual_df, "n_species"))
p.new_threat <- p.new_threat + enso_marks(yrs, tops_of(new_threat_df, "n_first_time_exposures_raw"))
combined_plot <- gridExtra::arrangeGrob(
  p.running,
  p.annual,
  p.new_threat,
  ncol = 1,
  nrow = 3,
  heights = c(1, 1, 1.35),   # panel C carries the legend below it
  left = grid::textGrob("Number of species", rot = 90)
)
plot_suffix <- if (is.finite(running_count_cap)) paste0('_cap', as.integer(running_count_cap)) else ''
plot.f.run <- paste0(plotDir, '/timeline_running_exposure_stacked_2025_v1', plot_suffix, '.pdf')

ggsave(filename = plot.f.run, plot = combined_plot, width = 180, height = 270, units = "mm")
ggsave(filename = sub("\\.pdf$", ".png", plot.f.run), plot = combined_plot, width = 180, height = 270, units = "mm", dpi = 300)

cat('Wrote running-total plot to: ', plot.f.run, '\n')
