# Fig5_Stressor.r
# Figure 5A: the dominant stressor among exposed species in each ecoregion in 2025.
# Figure 5B: how the dominant stressor changed between 2015-2024 and 2025.
#
# Each exposed species-year has one variable, the one with the largest exposed
# range fraction. Variables are grouped into heat, cold, extreme precipitation
# and drought. A species counts towards every ecoregion in which it has an
# exposed cell in 2025.
#
# Inputs:  singleMaxExposure_v1.qs, AllCellExposureSpXVar.qs, ecoregion raster and shapefile
# Outputs: two tables in figProposedMapDir, two figures in plotDir

suppressPackageStartupMessages({
  library(dplyr)
  library(sf)
  library(ggplot2)
  library(ggnewscale)
  library(grid)
})
suppressPackageStartupMessages({
  library(qs2)
  library(terra)
})

if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

run_figures_dir <- if (exists("figuresDir") && nzchar(figuresDir)) figuresDir else dirname(plotDir)
run_main_figures_dir <- plotDir
proposed_map_dir <- if (exists("figProposedMapDir") && nzchar(figProposedMapDir)) figProposedMapDir else file.path(run_figures_dir, "proposed_maps")

sf::sf_use_s2(FALSE)

dir.create(proposed_map_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(run_main_figures_dir, recursive = TRUE, showWarnings = FALSE)

map5_table <- file.path(proposed_map_dir, "map5_dominant_combined_stressor_ecoregion_2025_table.csv")
map8_table <- file.path(proposed_map_dir, "map8_dominant_stressor_transition_ecoregion_2015_2024_to_2025_table.csv")

single_max_exposure_file <- KEY_FILES$single_max_exposure_v1
sp_attributes_file <- KEY_FILES$sp_attributes_v8

if (is.na(single_max_exposure_file) || !file.exists(single_max_exposure_file)) {
  stop("Missing singleMaxExposure file: ", single_max_exposure_file)
}
if (is.na(sp_attributes_file) || !file.exists(sp_attributes_file)) {
  stop("Missing spAttributes file: ", sp_attributes_file)
}

var_to_cvar <- c(
  temp__12_up = "temp_heat",
  temp__3__max_up = "temp_heat",
  temp__12_lo = "temp_cold",
  temp__3__min_lo = "temp_cold",
  precip__12_up = "precip_flood",
  precip__3__max_up = "precip_flood",
  precip__12_lo = "precip_dry",
  precip__3__min_lo = "precip_dry"
)

cvar_labels <- c(
  temp_heat = "Heat",
  temp_cold = "Cold",
  precip_flood = "Pluvial",
  precip_dry = "Drought"
)

message("Building the stressor tables...")

ex11 <- qs_read(single_max_exposure_file) %>%
  arrange(spName, year, desc(pExp), var) %>%
  group_by(spName, year) %>%
  slice(1) %>%
  ungroup()

# Species-to-ecoregion links from the exposed cells, so that species contribute to
# every ecoregion where they are exposed
all_exp_file <- KEY_FILES$all_cell_exposure
if (is.null(all_exp_file) || !file.exists(all_exp_file)) {
  stop("Missing AllCellExposureSpXVar file: ", all_exp_file)
}

ecoregion_raster_file <- KEY_FILES$ecoregion_raster
if (!file.exists(ecoregion_raster_file)) {
  stop("Missing ecoregion raster: ", ecoregion_raster_file)
}

message("Loading all-cell exposures and joining to ECO_ID raster...")
all_exp <- qs_read(all_exp_file)
e_r <- terra::rast(ecoregion_raster_file)
e_r <- terra::crop(e_r, template)
eco_df <- tibble(
  cell = 1:terra::ncell(e_r),
  ECO_ID = terra::values(e_r)[, 1]
) %>%
  filter(!is.na(ECO_ID))

eco_link <- all_exp %>%
  filter(year == 2025) %>%
  distinct(spName, cell) %>%
  inner_join(eco_df, by = "cell") %>%
  distinct(spName, ECO_ID = as.integer(ECO_ID))

sp_2025 <- ex11 %>%
  filter(year == 2025) %>%
  transmute(spName, var, cvar = unname(var_to_cvar[var])) %>%
  filter(!is.na(cvar))

eco_cvar_counts <- eco_link %>%
  inner_join(sp_2025 %>% dplyr::select(spName, cvar), by = "spName") %>%
  distinct(ECO_ID, spName, cvar) %>%
  count(ECO_ID, cvar, name = "n_species")

map5 <- eco_cvar_counts %>%
  group_by(ECO_ID) %>%
  mutate(total_exposed_species = sum(n_species)) %>%
  arrange(ECO_ID, desc(n_species), cvar) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(dom_share_pct = 100 * n_species / total_exposed_species) %>%
  arrange(desc(total_exposed_species), ECO_ID)

sp_prev_decade <- ex11 %>%
  filter(year >= 2015, year <= 2024) %>%
  transmute(spName, cvar = unname(var_to_cvar[var])) %>%
  filter(!is.na(cvar)) %>%
  count(spName, cvar, name = "n_years") %>%
  arrange(spName, desc(n_years), cvar) %>%
  group_by(spName) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(spName, prev_cvar = cvar)

sp_curr_2025 <- ex11 %>%
  filter(year == 2025) %>%
  transmute(spName, cvar_2025 = unname(var_to_cvar[var])) %>%
  filter(!is.na(cvar_2025)) %>%
  distinct(spName, cvar_2025)

sp_transition <- sp_curr_2025 %>%
  inner_join(sp_prev_decade, by = "spName")

map8 <- eco_link %>%
  inner_join(sp_transition, by = "spName") %>%
  distinct(ECO_ID, spName, prev_cvar, cvar_2025) %>%
  count(ECO_ID, prev_cvar, cvar_2025, name = "n_species") %>%
  group_by(ECO_ID) %>%
  mutate(total_transition_species = sum(n_species)) %>%
  arrange(ECO_ID, desc(n_species), prev_cvar, cvar_2025) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    prev_label = unname(cvar_labels[prev_cvar]),
    curr_label = unname(cvar_labels[cvar_2025]),
    transition_label = ifelse(
      !is.na(prev_label) & !is.na(curr_label),
      paste0(prev_label, " -> ", curr_label),
      NA_character_
    ),
    transition_changed = prev_cvar != cvar_2025,
    dom_share_pct = 100 * n_species / total_transition_species
  ) %>%
  arrange(desc(total_transition_species), ECO_ID)

write.csv(map5, map5_table, row.names = FALSE)
write.csv(map8, map8_table, row.names = FALSE)

eco_shp <- KEY_FILES$ecoregion_shp
if (!file.exists(eco_shp)) stop("Missing Ecoregions2017 shapefile")

map5 <- read.csv(map5_table, stringsAsFactors = FALSE)
map8 <- read.csv(map8_table, stringsAsFactors = FALSE)

e_sf <- st_read(eco_shp, quiet = TRUE) %>%
  mutate(ECO_ID = suppressWarnings(as.integer(ECO_ID)))

if (nrow(e_sf) > 0) {
  e_sf <- st_simplify(e_sf, dTolerance = 0.03, preserveTopology = TRUE)
}

panelB <- map5 %>%
  transmute(ECO_ID = as.integer(ECO_ID), cvar = as.character(cvar), dom_share_pct = as.numeric(dom_share_pct)) %>%
  filter(!is.na(ECO_ID), !is.na(cvar))

panelC <- map8 %>%
  transmute(
    ECO_ID = as.integer(ECO_ID),
    transition_label = as.character(transition_label),
    transition_changed = as.logical(transition_changed),
    dom_share_pct = as.numeric(dom_share_pct)
  ) %>%
  filter(!is.na(ECO_ID), !is.na(transition_label))

mapB_sf <- e_sf %>% left_join(panelB, by = "ECO_ID")
mapC_sf <- e_sf %>% left_join(panelC, by = "ECO_ID")

cvar_levels <- c("temp_heat", "precip_flood", "precip_dry")
cvar_labels <- c(temp_heat = "Heat", precip_flood = "Extreme Precip.", precip_dry = "Drought")
cvar_palette <- c(temp_heat = "#d73027", precip_flood = "#4575b4", precip_dry = "#fdae61")

mapB_sf$cvar <- ifelse(mapB_sf$cvar %in% cvar_levels, mapB_sf$cvar, "other")
mapB_sf$cvar <- factor(mapB_sf$cvar, levels = c(cvar_levels, "other"))
cvar_palette_full <- c(cvar_palette, other = "grey70")
cvar_labels_full <- c(cvar_labels, other = "No exposure")

normalize_stressor_label <- function(x) {
  x_clean <- tolower(trimws(gsub("_", " ", x)))
  dplyr::case_when(
    grepl("heat|temp", x_clean) ~ "Heat",
    grepl("dry|drought", x_clean) ~ "Drought",
    grepl("flood|pluvial|wet", x_clean) ~ "Extreme Precip.",
    TRUE ~ NA_character_
  )
}

mapC_sf <- mapC_sf %>%
  mutate(
    transition_label_raw = as.character(transition_label),
    is_arrow = grepl("->", transition_label_raw),
    prev_raw = ifelse(is_arrow, trimws(sub("\\s*->.*$", "", transition_label_raw)), NA_character_),
    curr_raw = ifelse(is_arrow, trimws(sub("^.*->\\s*", "", transition_label_raw)), NA_character_),
    prev_stressor = normalize_stressor_label(prev_raw),
    curr_stressor = normalize_stressor_label(curr_raw),
    is_same = is_arrow & !is.na(prev_stressor) & !is.na(curr_stressor) & prev_stressor == curr_stressor,
    is_transition = is_arrow & !is.na(prev_stressor) & !is.na(curr_stressor) & prev_stressor != curr_stressor,
    transition_clean = ifelse(is_transition, paste0(prev_stressor, " \u2192 ", curr_stressor), NA_character_),
    same_clean = ifelse(is_same, curr_stressor, NA_character_)
  )

mapC_transition_sf <- mapC_sf %>% filter(is_transition, !is.na(transition_clean))
mapC_same_sf <- mapC_sf %>% filter(is_same, !is.na(same_clean))

pick_label <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_character_)
  val <- as.character(x[[1]])
  if (!nzchar(trimws(val))) return(NA_character_)
  val
}

run_name_label <- if (exists('RUN_NAME')) {
  pick_label(RUN_NAME)
} else if (exists('RUN_VERSION')) {
  pick_label(RUN_VERSION)
} else {
  NA_character_
}
if (is.na(run_name_label)) run_name_label <- 'V8'
run_name_label <- tolower(gsub('[^A-Za-z0-9]+', '_', trimws(run_name_label)))

trans_levels <- sort(unique(stats::na.omit(mapC_transition_sf$transition_clean)))
if (length(trans_levels) == 0) trans_levels <- "No observed transitions"
mapC_transition_sf$transition_clean <- factor(mapC_transition_sf$transition_clean, levels = trans_levels)
trans_palette <- setNames(grDevices::hcl.colors(length(trans_levels), palette = "Dynamic"), trans_levels)

same_levels <- c("Heat", "Drought", "Extreme Precip.")
mapC_same_sf$same_clean <- factor(mapC_same_sf$same_clean, levels = same_levels)
same_palette <- c(Heat = cvar_palette[["temp_heat"]], Drought = cvar_palette[["precip_dry"]], `Extreme Precip.` = cvar_palette[["precip_flood"]])

common_map_theme <- theme_minimal(base_size = 13) +
  theme(
    legend.background = element_blank(),
    legend.box.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "right",
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, 0),
    panel.background = element_rect(fill = "#d6eaf8", color = NA),
    panel.grid.major = element_line(color = "white", linewidth = 0.25),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.title = element_text(size = 14, face = "bold"),
    plot.margin = margin(0, 2, 0, 2)
  )

pB <- ggplot() +
  geom_sf(data = e_sf, fill = "grey92", color = "grey70", linewidth = 0.1) +
  geom_sf(data = mapB_sf, aes(fill = cvar, geometry = geometry), color = NA, linewidth = 0) +
  scale_fill_manual(values = cvar_palette_full, labels = cvar_labels_full, drop = FALSE, name = "Dominant\n2025 stressor") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  labs(title = "A. Dominant stressor (2025)") +
  common_map_theme

# Muted (lightened/desaturated) versions of same_palette for background ecoregions
lighten_hex <- function(hex, amount = 0.55) {
  rgb_mat <- col2rgb(hex) / 255
  rgb_light <- rgb_mat + (1 - rgb_mat) * amount
  grDevices::rgb(rgb_light[1, ], rgb_light[2, ], rgb_light[3, ])
}
same_palette_muted <- setNames(
  sapply(same_palette, lighten_hex, amount = 0.35),
  names(same_palette)
)

# Semantic palette: FROM Heat = reds, FROM Precip = blues, FROM Drought = yellows/golds
semantic_trans_palette <- c(
  "Heat \u2192 Drought"         = "#67000D",   # very dark crimson
  "Heat \u2192 Extreme Precip." = "#EF3B2C",   # bright vivid red
  "Extreme Precip. \u2192 Heat" = "#1A3A8F",   # deep blue
  "Extreme Precip. \u2192 Drought" = "#4393C3", # steel blue
  "Drought \u2192 Heat"         = "#FF9900",   # vivid orange
  "Drought \u2192 Extreme Precip." = "#FFE135"  # bright canary yellow
)
trans_palette_vivid <- setNames(
  sapply(trans_levels, function(lv) {
    if (lv %in% names(semantic_trans_palette)) semantic_trans_palette[[lv]]
    else "#555555"
  }),
  trans_levels
)

pC <- ggplot() +
  geom_sf(data = e_sf, fill = "grey92", color = "grey70", linewidth = 0.1) +
  geom_sf(data = mapC_same_sf, aes(fill = same_clean, geometry = geometry), color = NA, linewidth = 0) +
  scale_fill_manual(values = same_palette_muted, breaks = same_levels, drop = FALSE, name = "Stayed same") +
  ggnewscale::new_scale_fill() +
  geom_sf(data = mapC_transition_sf, aes(fill = transition_clean, geometry = geometry), color = "black", linewidth = 0.12) +
  scale_fill_manual(values = trans_palette_vivid, drop = FALSE, name = "Transitioned") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  labs(title = "B. Dominant stressor transition (2015-2024 to 2025)") +
  common_map_theme +
  theme(legend.box = "vertical")

map5_png <- file.path(run_main_figures_dir, paste0("Fig5_StressorMap_", run_name_label, ".png"))
figS5_png <- file.path(run_main_figures_dir, paste0("Fig5b_StressorSwitching_", run_name_label, ".png"))

# Figure 5A
png(map5_png, width = 2400, height = 1000, res = 300)
print(pB)
dev.off()
message("Saved Figure 5A: ", map5_png)

# Figure 5B
png(figS5_png, width = 2400, height = 1000, res = 300)
print(pC)
dev.off()
message("Saved Figure 5B: ", figS5_png)

