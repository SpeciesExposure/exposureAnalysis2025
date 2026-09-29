# Fig1_FigS1_Time.r
# Figure 1: number of exposed species per year, 1990-2025, coloured by the
# climate variable with the largest exposed range fraction.
# Figure S1: the same counts split by taxonomic group, IUCN category, body-size
# group and range-size group.
# A species is exposed in a year when at least 25% of its range exceeds the
# historical limit of at least one variable.
#
# Inputs:  singleMaxExposure_v1.qs, spAttributes_v8.qs, data/body_mass_quartiles.csv
# Outputs: figures in plotDir

if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

ex11_file <- file.path(intDir, 'singleMaxExposure_v1.qs')
if (!file.exists(ex11_file)) {
	stop(
		"Missing ex11 file for RUN_VERSION ", RUN_VERSION, ": ", ex11_file,
		"\nRun 6_MakeDF.r to create it."
	)
}

ex11 <- qs_read(ex11_file)

bottom_panels_as_proportion <- if (exists("bottom_panels_as_proportion")) bottom_panels_as_proportion else FALSE
message("bottom_panels_as_proportion = ", bottom_panels_as_proportion)

if (!"redlistCategory" %in% names(ex11)) {
	spMeta_iucn <- qs_read(KEY_FILES$sp_attributes_v8) %>% dplyr::select(spName, redlistCategory) %>% unique
	ex11 <- ex11 %>% left_join(spMeta_iucn, by = "spName")
}

bs_quart_file <- KEY_FILES$body_mass_quartiles
if (!file.exists(path.expand(bs_quart_file))) {
	stop("Body-size quartile file not found: ", bs_quart_file)
}
bs_quart <- read.csv(path.expand(bs_quart_file), stringsAsFactors = FALSE) %>%
	dplyr::select(species, body_mass_quartile_llm) %>%
	distinct()
ex11 <- ex11 %>% left_join(bs_quart, by = c("spName" = "species"))

range_meta_candidates <- c(
	file.path(metaDir, 'spAttributes_v8.qs')
)
range_meta_file <- range_meta_candidates[file.exists(path.expand(range_meta_candidates))][1]

if (!is.na(range_meta_file)) {
	range_meta <- qs_read(path.expand(range_meta_file))
	if (!"spName" %in% names(range_meta) && "sp" %in% names(range_meta)) {
		range_meta <- range_meta %>% rename(spName = sp)
	}
	if ("rangeSizeKm" %in% names(range_meta)) {
		range_meta <- range_meta %>%
			dplyr::select(spName, rangeSizeKm) %>%
			filter(!is.na(spName), !is.na(rangeSizeKm), rangeSizeKm > 0) %>%
			distinct(spName, .keep_all = TRUE) %>%
			mutate(range_size_quartile = dplyr::ntile(rangeSizeKm, 4L))
		if ("rangeSizeKm" %in% names(ex11)) ex11 <- ex11 %>% dplyr::select(-rangeSizeKm)
		if ("range_size_quartile" %in% names(ex11)) ex11 <- ex11 %>% dplyr::select(-range_size_quartile)
		ex11 <- ex11 %>% left_join(range_meta, by = "spName")
	} else {
		message("rangeSizeKm not found in metadata file: ", range_meta_file)
	}
} else {
	message("Range-size metadata file not found in metaDir")
}

# Enforce one max-exposure variable per species-year (handles older cached files with ties)
n_before <- nrow(ex11)
ex11 <- ex11 %>%
	arrange(spName, year, desc(pExp), var) %>%
	group_by(spName, year) %>%
	slice(1) %>%
	ungroup()
n_after <- nrow(ex11)
if (n_after < n_before) {
	message("Deduplicated species-year ties in ex11: removed ", n_before - n_after, " rows")
}

# Cache the histogram summary so FigS2_TimeSensitivity.r can read it without
# re-loading the full ex11 qs file.
hist_cache_dir  <- file.path(plotDir, "FigS2_cache")
dir.create(hist_cache_dir, showWarnings = FALSE, recursive = TRUE)
hist_cache_file <- file.path(hist_cache_dir, paste0("hist_df_", runTag, "_v1.rds"))

hist_df <- ex11 %>%
  group_by(year, var) %>%
  summarise(n_species = n_distinct(spName), .groups = "drop")

if (!file.exists(hist_cache_file)) {
  saveRDS(hist_df, hist_cache_file)
  message("Saved histogram summary cache: ", hist_cache_file)
}

p.single.hist <- hist_df %>%
  ggplot(aes(x = year, y = n_species, fill = var)) +
  geom_col(width = 0.8) +
  scale_x_continuous(breaks = unique(ex11$year)) +
  scale_fill_manual(
    labels = c(
			"temp__3__max_up"   = "High temperature (warmest 3 months)",
			"temp__12_up"       = "High temperature (annual mean)",
			"temp__3__min_lo"   = "Low temperature (coldest 3 months)",
			"temp__12_lo"       = "Low temperature (annual mean)",
			"precip__3__min_lo" = "Low precipitation (driest 3 months)",
			"precip__3__max_up" = "High precipitation (wettest 3 months)",
			"precip__12_up"     = "High precipitation (annual total)",
			"precip__12_lo"     = "Low precipitation (annual total)"
    ),
    values = c(
      "temp__3__max_up"   = "#d73027",
      "temp__12_up"       = "#fc8d59",
      "temp__3__min_lo"   = "#bdbdbd",
      "temp__12_lo"       = "#737373",
      "precip__3__max_up" = "#4575b4",
      "precip__12_up"     = "#74add1",
      "precip__3__min_lo" = "#fee090",
      "precip__12_lo"     = "#fdbb63"
    )
  ) +
  labs(x = "Year", y = "Number of species", fill = NULL) +
  theme_classic() +
  theme(
		axis.text.x      = element_text(size = 11, color = "black", angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y      = element_text(size = 11, color = "black"),
    axis.title       = element_text(size = 13, color = "black"),
		legend.text      = element_text(size = 10),
		legend.key.size  = unit(0.5, "cm"),
		legend.position  = c(0.02, 0.98),
		legend.justification = c(0, 1),
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.line        = element_blank()
  )

year_total_exposed <- ex11 %>%
	group_by(year) %>%
	summarise(n_total_exposed = n_distinct(spName), .groups = "drop")

# 2 x 2 small panel block below main panel
group_levels <- c("reptiles", "amphibians", "birds", "mammals")

y_var <- if (isTRUE(bottom_panels_as_proportion)) "prop_species" else "n_species"
y_lab <- if (isTRUE(bottom_panels_as_proportion)) "Proportion of species exposed" else "Species"
y_scale_small <- if (isTRUE(bottom_panels_as_proportion)) {
	scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 1))
} else {
	scale_y_continuous()
}

small1_df <- ex11 %>%
	mutate(group = tolower(trimws(as.character(group)))) %>%
	filter(group %in% group_levels) %>%
	filter(year >= 2010) %>%
	group_by(year, group) %>%
	summarise(n_species = n_distinct(spName), .groups = "drop") %>%
	left_join(year_total_exposed, by = "year") %>%
	tidyr::complete(
		year = seq(2010, max(ex11$year), by = 1),
		group = group_levels,
		fill = list(n_species = 0, n_total_exposed = 0)
	) %>%
	mutate(prop_species = ifelse(n_total_exposed > 0, n_species / n_total_exposed, 0)) %>%
	mutate(group = factor(group, levels = group_levels))

if (sum(small1_df$n_species, na.rm = TRUE) == 0) {
	warning("Taxonomic Groups panel has no matching rows after group/year filtering. Check ex11$group labels and selected version.")
}

p.small1 <- ggplot(small1_df, aes(x = year, y = .data[[y_var]], fill = group)) +
	geom_col(width = 0.8) +
	scale_x_continuous(breaks = seq(2010, max(ex11$year), by = 5)) +
	y_scale_small +
	scale_fill_manual(values = c(
		"reptiles" = "#BDD7E7",
		"amphibians" = "#6BAED6",
		"birds" = "#3182BD",
		"mammals" = "#08519C"
	), labels = c(
		"reptiles" = "Reptiles",
		"amphibians" = "Amphibians",
		"birds" = "Birds",
		"mammals" = "Mammals"
	)) +
	labs(x = "Year", y = y_lab, title = "Taxonomic Groups", fill = NULL) +
	theme_classic() +
	theme(
		axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
		axis.text.y = element_text(size = 7),
		axis.title = element_text(size = 8),
		legend.position = c(0.02, 0.98),
		legend.justification = c(0, 1),
		legend.text = element_text(size = 6),
		legend.key.size = unit(0.35, "cm"),
		plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
	)

small2_df <- ex11 %>%
	filter(year >= 2010) %>%
	mutate(redlistCategory = ifelse(is.na(redlistCategory), "Unknown", gsub("_", " ", redlistCategory))) %>%
	group_by(year, redlistCategory) %>%
	summarise(n_species = n_distinct(spName), .groups = "drop") %>%
	left_join(year_total_exposed, by = "year") %>%
	tidyr::complete(
		year = seq(2010, max(ex11$year), by = 1),
		redlistCategory,
		fill = list(n_species = 0, n_total_exposed = 0)
	) %>%
	mutate(prop_species = ifelse(n_total_exposed > 0, n_species / n_total_exposed, 0)) %>%
	{
		std_levels <- c(
			"Critically Endangered", "Endangered", "Vulnerable", "Near Threatened", "Least Concern",
			"Data Deficient", "Not Evaluated", "Extinct", "Extinct in the Wild", "Unknown"
		)
		obs <- unique(.$redlistCategory)
		lvl <- c(std_levels[std_levels %in% obs], sort(setdiff(obs, std_levels)))
		mutate(., redlistCategory = factor(redlistCategory, levels = lvl))
	}

iucn_legend_breaks <- setdiff(as.character(levels(small2_df$redlistCategory)), "Unknown")

	p.small2 <- ggplot(small2_df, aes(x = year, y = .data[[y_var]], fill = redlistCategory)) +
	geom_col(width = 0.8) +
	scale_x_continuous(breaks = seq(2010, max(ex11$year), by = 5)) +
	y_scale_small +
	scale_fill_manual(values = c(
		"Critically Endangered" = "#67000D",
		"Endangered" = "#7F0000",
		"Vulnerable" = "#B30000",
		"Near Threatened" = "#D7301F",
		"Least Concern" = "#EF6548",
		"Data Deficient" = "#FC8D59",
		"Not Evaluated" = "#FDBB84",
		"Extinct" = "#3B0A0A",
		"Extinct in the Wild" = "#525252",
		"Unknown" = "#969696"
	), breaks = iucn_legend_breaks) +
	labs(x = "Year", y = y_lab, title = "IUCN Categories", fill = NULL) +
	theme_classic() +
	theme(
		axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
		axis.text.y = element_text(size = 7),
		axis.title = element_text(size = 8),
		legend.position = c(0.02, 0.98),
		legend.justification = c(0, 1),
		legend.text = element_text(size = 6),
		legend.key.size = unit(0.35, "cm"),
		plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
	)

make_placeholder <- function(title_txt) {
	ggplot() +
		annotate("text", x = 0.5, y = 0.5, label = "placeholder", size = 4, color = "grey40") +
		xlim(0, 1) + ylim(0, 1) +
		labs(title = title_txt, x = NULL, y = NULL) +
		theme_classic() +
		theme(
			axis.text = element_blank(),
			axis.ticks = element_blank(),
			axis.line = element_blank(),
			plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
			panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
		)
}

small3_df <- ex11 %>%
	filter(year >= 2010) %>%
	mutate(
		body_mass_quartile_llm = case_when(
			is.na(body_mass_quartile_llm) ~ "Unknown",
			body_mass_quartile_llm == 1 ~ "Smallest body-size",
			body_mass_quartile_llm == 2 ~ "Lower-middle body-size",
			body_mass_quartile_llm == 3 ~ "Upper-middle body-size",
			body_mass_quartile_llm == 4 ~ "Largest body-size",
			TRUE ~ "Unknown"
		)
	) %>%
	group_by(year, body_mass_quartile_llm) %>%
	summarise(n_species = n_distinct(spName), .groups = "drop") %>%
	left_join(year_total_exposed, by = "year") %>%
	tidyr::complete(
		year = seq(2010, max(ex11$year), by = 1),
		body_mass_quartile_llm,
		fill = list(n_species = 0, n_total_exposed = 0)
	) %>%
	mutate(prop_species = ifelse(n_total_exposed > 0, n_species / n_total_exposed, 0)) %>%
	mutate(body_mass_quartile_llm = factor(body_mass_quartile_llm, levels = c("Smallest body-size", "Lower-middle body-size", "Upper-middle body-size", "Largest body-size", "Unknown")))

p.small3 <- ggplot(small3_df, aes(x = year, y = .data[[y_var]], fill = body_mass_quartile_llm)) +
	geom_col(width = 0.8) +
	scale_x_continuous(breaks = seq(2010, max(ex11$year), by = 5)) +
	y_scale_small +
	scale_fill_manual(values = c(
		"Smallest body-size" = "#E6E6E6",
		"Lower-middle body-size" = "#BDBDBD",
		"Upper-middle body-size" = "#737373",
		"Largest body-size" = "#252525",
		"Unknown" = "#000000"
	)) +
	labs(x = "Year", y = y_lab, title = "Body-size Quartiles", fill = NULL) +
	theme_classic() +
	theme(
		axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
		axis.text.y = element_text(size = 7),
		axis.title = element_text(size = 8),
		legend.position = c(0.02, 0.98),
		legend.justification = c(0, 1),
		legend.text = element_text(size = 6),
		legend.key.size = unit(0.35, "cm"),
		plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
	)

small4_df <- ex11 %>%
	filter(year >= 2010) %>%
	mutate(
		range_size_quartile = case_when(
			is.na(range_size_quartile) ~ "Unknown",
			range_size_quartile == 1 ~ "Smallest ranges",
			range_size_quartile == 2 ~ "Lower-middle ranges",
			range_size_quartile == 3 ~ "Upper-middle ranges",
			range_size_quartile == 4 ~ "Largest ranges",
			TRUE ~ "Unknown"
		)
	) %>%
	group_by(year, range_size_quartile) %>%
	summarise(n_species = n_distinct(spName), .groups = "drop") %>%
	left_join(year_total_exposed, by = "year") %>%
	tidyr::complete(
		year = seq(2010, max(ex11$year), by = 1),
		range_size_quartile,
		fill = list(n_species = 0, n_total_exposed = 0)
	) %>%
	mutate(prop_species = ifelse(n_total_exposed > 0, n_species / n_total_exposed, 0)) %>%
	mutate(range_size_quartile = factor(range_size_quartile, levels = c("Smallest ranges", "Lower-middle ranges", "Upper-middle ranges", "Largest ranges", "Unknown")))

p.small4 <- ggplot(small4_df, aes(x = year, y = .data[[y_var]], fill = range_size_quartile)) +
	geom_col(width = 0.8) +
	scale_x_continuous(breaks = seq(2010, max(ex11$year), by = 5)) +
	y_scale_small +
	scale_fill_manual(values = c(
		"Smallest ranges" = "#FDD49E",
		"Lower-middle ranges" = "#FDBB84",
		"Upper-middle ranges" = "#EF6548",
		"Largest ranges" = "#B63600",
		"Unknown" = "#5A5A5A"
	)) +
	labs(x = "Year", y = y_lab, title = "Range-size Quartiles", fill = NULL) +
	theme_classic() +
	theme(
		axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1),
		axis.text.y = element_text(size = 7),
		axis.title = element_text(size = 8),
		legend.position = c(0.02, 0.98),
		legend.justification = c(0, 1),
		legend.text = element_text(size = 6),
		legend.key.size = unit(0.35, "cm"),
		plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
		panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
	)

small_block <- gridExtra::arrangeGrob(
	p.small1, p.small2,
	p.small3, p.small4,
	ncol = 2
)

p.main <- p.single.hist +
	labs(title = "Max Exposure Type per Species-Year") +
	theme(
		legend.position = c(0.02, 0.98),
		legend.justification = c(0, 1),
		plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
	)

# ---- Mark El Niño (red up-triangle) and La Niña (blue down-triangle) years above the bars ----
source(file.path(PROJECT_PATHS$figures_src, "enso_strip.R"))
tops_all <- data.frame(year = year_total_exposed$year, top = year_total_exposed$n_total_exposed)
p.main <- p.main + enso_marks(min(ex11$year):max(ex11$year), tops_all) + enso_key(x = 2003, y = 5700, dy = 450)
tops_small <- tops_all[tops_all$year >= 2010, ]
p.small1 <- p.small1 + enso_marks(2010:max(ex11$year), tops_small, size = 1.8) + enso_key(x = 2016.2, y = 5800, dy = 900, size = 1.9, text_size = 2.7)
p.small2 <- p.small2 + enso_marks(2010:max(ex11$year), tops_small, size = 1.8)
p.small3 <- p.small3 + enso_marks(2010:max(ex11$year), tops_small, size = 1.8)
p.small4 <- p.small4 + enso_marks(2010:max(ex11$year), tops_small, size = 1.8)
small_block <- gridExtra::arrangeGrob(p.small1, p.small2, p.small3, p.small4, ncol = 2)
# Fig. 1 = main panel, Fig. S1 = small block
ggsave(paste0(plotDir, '/Fig1_main_ENSO_', runTag, '.png'), p.main, width = 300, height = 110, units = "mm", dpi = 300)
ggsave(paste0(plotDir, '/Fig1_main_ENSO_', runTag, '.pdf'), p.main, width = 300, height = 110, units = "mm")
ggsave(paste0(plotDir, '/FigS1_small_ENSO_', runTag, '.png'), small_block, width = 300, height = 110, units = "mm", dpi = 300)
ggsave(paste0(plotDir, '/FigS1_small_ENSO_', runTag, '.pdf'), small_block, width = 300, height = 110, units = "mm")
multi_panel_plot <- gridExtra::arrangeGrob(
	p.main,
	small_block,
	ncol = 1,
	heights = c(1, 1)
)

plot_suffix <- if (isTRUE(bottom_panels_as_proportion)) '_propBottom' else '_countBottom'
plot.f2 <- paste0(plotDir, '/timeline_single_stacked_2025_v1', plot_suffix, '.pdf')
ggsave(filename = plot.f2, plot = multi_panel_plot, width = 300, height = 220, units = "mm")

cat('Wrote plot to: ', plot.f2, '\n')
