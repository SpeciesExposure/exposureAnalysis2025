# Fig4_FigS5_Maps.r
# Figure 4: number of species exposed in 2025, per ecoregion (panel A) and per
# 0.25 degree cell in three regions (panels B-D).
# Figure S5: number of chronically exposed species per ecoregion. A species is
# counted as chronically exposed here when it has at least one exposed cell in
# 6 or more years between 1990 and 2025.
#
# A species is counted in a cell or ecoregion when it has at least one exposed
# cell there. The 25% of range rule is not applied in these maps.
#
# Inputs:  AllCellExposureSpXVar.qs, landTemplate.tif, ecoregion raster and shapefile
# Outputs: rasters and tables in intDir/MapProducts_4_Maps, figures in plotDir

suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(qs2)
})

if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

# ---- user-tunable settings --------------------------------------------------
target_year <- 2025
chronic_min_years <- 6
chronic_year_min <- 1990
chronic_year_max <- target_year

# ---- inputs ----------------------------------------------------------------
all_exp_file <- if (exists("KEY_FILES") && !is.null(KEY_FILES$all_cell_exposure)) {
  KEY_FILES$all_cell_exposure
} else {
  file.path(intDir, "AllCellExposureSpXVar.qs")
}

template_file <- file.path(intDir, "landTemplate.tif")
ecoregion_raster_file <- KEY_FILES$ecoregion_raster

if (!file.exists(all_exp_file)) stop("Missing all-cell exposure file: ", all_exp_file)
if (!file.exists(template_file)) stop("Missing template raster: ", template_file)
if (!file.exists(ecoregion_raster_file)) stop("Missing ecoregion raster: ", ecoregion_raster_file)

# ---- outputs ---------------------------------------------------------------
out_dir <- file.path(intDir, "MapProducts_4_Maps")
out_dir_rasters <- file.path(out_dir, "rasters")
out_dir_tables <- file.path(out_dir, "tables")
out_dir_figures <- plotDir
dir.create(out_dir_rasters, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir_figures, recursive = TRUE, showWarnings = FALSE)

suffix <- paste0("_", target_year, "_", RUN_VERSION)

cell_raster_file <- file.path(out_dir_rasters, paste0("Exposure_CellRichness", suffix, ".tif"))
eco_raster_out_file <- file.path(out_dir_rasters, paste0("Exposure_EcoregionRichness", suffix, ".tif"))
chronic_cell_raster_file <- file.path(
  out_dir_rasters,
  paste0("Exposure_ChronicCellRichness_GE", chronic_min_years, "yrs", suffix, ".tif")
)
chronic_eco_raster_file <- file.path(
  out_dir_rasters,
  paste0("Exposure_ChronicEcoregionRichness_GE", chronic_min_years, "yrs", suffix, ".tif")
)

cell_table_file <- file.path(out_dir_tables, paste0("Exposure_CellRichness", suffix, ".csv"))
eco_table_file <- file.path(out_dir_tables, paste0("Exposure_EcoregionRichness", suffix, ".csv"))
chronic_cell_table_file <- file.path(
  out_dir_tables,
  paste0("Exposure_ChronicCellRichness_GE", chronic_min_years, "yrs", suffix, ".csv")
)
chronic_eco_table_file <- file.path(
  out_dir_tables,
  paste0("Exposure_ChronicEcoregionRichness_GE", chronic_min_years, "yrs", suffix, ".csv")
)

write_cell_count_raster <- function(counts_df, template_r, out_file, value_col = "n") {
  if (!all(c("cell", value_col) %in% names(counts_df))) {
    stop("counts_df must include columns: cell and ", value_col)
  }
  r <- template_r
  values(r)[] <- NA_real_
  values(r)[counts_df$cell] <- as.numeric(counts_df[[value_col]])
  writeRaster(r, out_file, overwrite = TRUE)
  invisible(r)
}

write_ecoregion_count_raster <- function(eco_counts_df, eco_r, out_file) {
  if (!all(c("ECO_ID", "n") %in% names(eco_counts_df))) {
    stop("eco_counts_df must include columns: ECO_ID and n")
  }
  id_to_count <- setNames(as.numeric(eco_counts_df$n), as.character(eco_counts_df$ECO_ID))
  eco_vals <- values(eco_r)[, 1]
  out_vals <- id_to_count[as.character(eco_vals)]

  r_out <- eco_r
  values(r_out) <- as.numeric(out_vals)
  writeRaster(r_out, out_file, overwrite = TRUE)
  invisible(r_out)
}

cat("Reading all-cell exposure events:\n", all_exp_file, "\n")
all_exp <- qs_read(all_exp_file)

required_cols <- c("spName", "year", "cell")
missing_cols <- setdiff(required_cols, names(all_exp))
if (length(missing_cols) > 0) {
  stop("AllCellExposureSpXVar is missing required columns: ", paste(missing_cols, collapse = ", "))
}

template_r <- rast(template_file)
eco_r <- rast(ecoregion_raster_file)
eco_r <- crop(eco_r, template_r)
if (!terra::compareGeom(template_r, eco_r, stopOnError = FALSE)) {
  stop("Template and ecoregion rasters are not geometrically aligned after crop().")
}

all_exp <- all_exp %>% mutate(cell = suppressWarnings(as.integer(cell)))
bad_cells <- all_exp %>%
  filter(is.na(cell) | cell < 1L | cell > ncell(template_r))
if (nrow(bad_cells) > 0) {
  stop(
    "Invalid cell ids in AllCellExposureSpXVar.qs. Bad rows: ",
    nrow(bad_cells),
    " (expected integer cell in [1, ", ncell(template_r), "])."
  )
}

eco_lookup <- tibble(
  cell = seq_len(ncell(eco_r)),
  ECO_ID = values(eco_r)[, 1]
) %>%
  filter(!is.na(ECO_ID))

all_exp_year <- all_exp %>%
  filter(year == target_year)

cat("Target year rows:", nrow(all_exp_year), "\n")
if (nrow(all_exp_year) == 0) {
  stop("No exposure rows found for target_year=", target_year)
}

# 1) Cell exposure richness (unique species per cell)
cell_counts <- all_exp_year %>%
  distinct(spName, cell) %>%
  count(cell, name = "n") %>%
  arrange(desc(n), cell)

write.csv(cell_counts, cell_table_file, row.names = FALSE)
write_cell_count_raster(cell_counts, template_r, cell_raster_file)

# 2) Ecoregion exposure richness (unique species per ecoregion)
eco_link <- all_exp_year %>%
  distinct(spName, cell) %>%
  inner_join(eco_lookup, by = "cell") %>%
  distinct(spName, ECO_ID = as.integer(ECO_ID))

eco_counts <- eco_link %>%
  count(ECO_ID, name = "n") %>%
  arrange(desc(n), ECO_ID)

write.csv(eco_counts, eco_table_file, row.names = FALSE)
write_ecoregion_count_raster(eco_counts, eco_r, eco_raster_out_file)

# 3) Chronic exposure richness at cell level
chronic_species <- all_exp %>%
  filter(year >= chronic_year_min, year <= chronic_year_max) %>%
  distinct(spName, year) %>%
  count(spName, name = "n_exposure_years") %>%
  filter(n_exposure_years >= chronic_min_years) %>%
  dplyr::select(spName)

cat("Chronic species (>=", chronic_min_years, " years): ", nrow(chronic_species), "\n", sep = "")

chronic_cell_counts <- all_exp_year %>%
  semi_join(chronic_species, by = "spName") %>%
  distinct(spName, cell) %>%
  count(cell, name = "n") %>%
  arrange(desc(n), cell)

write.csv(chronic_cell_counts, chronic_cell_table_file, row.names = FALSE)
write_cell_count_raster(chronic_cell_counts, template_r, chronic_cell_raster_file)

# 4) Chronic exposure richness at ecoregion level (analogous to 2025 eco calc)
chronic_eco_link <- all_exp_year %>%
  semi_join(chronic_species, by = "spName") %>%
  distinct(spName, cell) %>%
  inner_join(eco_lookup, by = "cell") %>%
  distinct(spName, ECO_ID = as.integer(ECO_ID))

chronic_eco_counts <- chronic_eco_link %>%
  count(ECO_ID, name = "n") %>%
  arrange(desc(n), ECO_ID)

write.csv(chronic_eco_counts, chronic_eco_table_file, row.names = FALSE)
write_ecoregion_count_raster(chronic_eco_counts, eco_r, chronic_eco_raster_file)

cat("\nSaved rasters:\n")
cat(" - ", cell_raster_file, "\n", sep = "")
cat(" - ", eco_raster_out_file, "\n", sep = "")
cat(" - ", chronic_cell_raster_file, "\n", sep = "")
cat(" - ", chronic_eco_raster_file, "\n", sep = "")

cat("\nSaved tables:\n")
cat(" - ", cell_table_file, "\n", sep = "")
cat(" - ", eco_table_file, "\n", sep = "")
cat(" - ", chronic_cell_table_file, "\n", sep = "")
cat(" - ", chronic_eco_table_file, "\n", sep = "")


# ---- plotting + validation------------------------------------------
suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(scales)
  library(tidyterra)
  library(cowplot)
  library(grid)
  library(rnaturalearth)
  library(magick)
})

sf::sf_use_s2(FALSE)

eco_shp_file <- KEY_FILES$ecoregion_shp
if (!file.exists(eco_shp_file)) stop("Missing ecoregion shapefile: ", eco_shp_file)

e_sf <- st_read(eco_shp_file, quiet = TRUE) %>%
  mutate(ECO_ID = suppressWarnings(as.integer(ECO_ID)))


cm.cols2 <- function(x, bias = 1) {
  colorRampPalette(c("steelblue4", "steelblue1", "gold", "red1", "red4"), bias = bias)(x)
}
my_palette <- cm.cols2(100, bias = 3.0)

common_map_theme <- theme_minimal(base_size = 12) +
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

# Validation: ensure rasters encode the same counts generated from all_exp.
r_cell <- rast(cell_raster_file)[[1]]
r_eco <- rast(eco_raster_out_file)[[1]]
r_chronic <- rast(chronic_cell_raster_file)[[1]]
r_chronic_eco <- rast(chronic_eco_raster_file)[[1]]
r_diff_cell_minus_chronic <- r_cell - r_chronic
names(r_cell) <- "value"
names(r_eco) <- "value"
names(r_chronic) <- "value"
names(r_chronic_eco) <- "value"
names(r_diff_cell_minus_chronic) <- "value"

cell_val_tbl <- tibble(
  cell = seq_len(ncell(r_cell)),
  n_r = values(r_cell)[, 1]
) %>%
  filter(!is.na(n_r))

validate_exact_key_value <- function(df_tbl, ras_tbl, key_col, val_df_col = "n_df", val_r_col = "n_r", label = "validation") {
  if (nrow(df_tbl) == 0 && nrow(ras_tbl) == 0) return(invisible(TRUE))
  if (nrow(df_tbl) == 0 || nrow(ras_tbl) == 0) {
    stop(label, " failed: one side is empty and the other is not.")
  }

  miss_in_r <- dplyr::anti_join(df_tbl, ras_tbl, by = key_col)
  extra_in_r <- dplyr::anti_join(ras_tbl, df_tbl, by = key_col)
  if (nrow(miss_in_r) > 0 || nrow(extra_in_r) > 0) {
    stop(label, " failed: key coverage mismatch between dataframe and raster.")
  }

  chk <- dplyr::left_join(df_tbl, ras_tbl, by = key_col)
  if (!all(chk[[val_df_col]] == chk[[val_r_col]])) {
    stop(label, " failed: value mismatch.")
  }
  invisible(TRUE)
}

cell_df_tbl <- cell_counts %>% transmute(cell, n_df = as.numeric(n))
cell_r_tbl <- cell_val_tbl %>% transmute(cell, n_r = as.numeric(n_r))
validate_exact_key_value(cell_df_tbl, cell_r_tbl, key_col = "cell", label = "Cell richness check")

chronic_val_tbl <- tibble(
  cell = seq_len(ncell(r_chronic)),
  n_r = values(r_chronic)[, 1]
) %>%
  filter(!is.na(n_r))

chronic_df_tbl <- chronic_cell_counts %>% transmute(cell, n_df = as.numeric(n))
chronic_r_tbl <- chronic_val_tbl %>% transmute(cell, n_r = as.numeric(n_r))
validate_exact_key_value(chronic_df_tbl, chronic_r_tbl, key_col = "cell", label = "Chronic cell richness check")

eco_val_tbl <- tibble(
  cell = seq_len(ncell(r_eco)),
  n_r = values(r_eco)[, 1]
) %>%
  inner_join(eco_lookup, by = "cell") %>%
  filter(!is.na(ECO_ID), !is.na(n_r)) %>%
  group_by(ECO_ID) %>%
  summarise(
    n_unique = n_distinct(n_r),
    n_r = first(n_r),
    .groups = "drop"
  )

if (any(eco_val_tbl$n_unique > 1)) {
  stop("Ecoregion raster contains multiple values within at least one ECO_ID.")
}

eco_val_tbl <- eco_val_tbl %>% dplyr::select(ECO_ID, n_r)

eco_df_tbl <- eco_counts %>% transmute(ECO_ID = as.integer(ECO_ID), n_df = as.numeric(n))
eco_r_tbl <- eco_val_tbl %>% transmute(ECO_ID = as.integer(ECO_ID), n_r = as.numeric(n_r))
validate_exact_key_value(eco_df_tbl, eco_r_tbl, key_col = "ECO_ID", label = "Ecoregion richness check")

chronic_eco_val_tbl <- tibble(
  cell = seq_len(ncell(r_chronic_eco)),
  n_r = values(r_chronic_eco)[, 1]
) %>%
  inner_join(eco_lookup, by = "cell") %>%
  filter(!is.na(ECO_ID), !is.na(n_r)) %>%
  group_by(ECO_ID) %>%
  summarise(
    n_unique = n_distinct(n_r),
    n_r = first(n_r),
    .groups = "drop"
  )

if (any(chronic_eco_val_tbl$n_unique > 1)) {
  stop("Chronic ecoregion raster contains multiple values within at least one ECO_ID.")
}

chronic_eco_val_tbl <- chronic_eco_val_tbl %>% dplyr::select(ECO_ID, n_r)

chronic_eco_df_tbl <- chronic_eco_counts %>% transmute(ECO_ID = as.integer(ECO_ID), n_df = as.numeric(n))
chronic_eco_r_tbl <- chronic_eco_val_tbl %>% transmute(ECO_ID = as.integer(ECO_ID), n_r = as.numeric(n_r))
validate_exact_key_value(chronic_eco_df_tbl, chronic_eco_r_tbl, key_col = "ECO_ID", label = "Chronic ecoregion richness check")

if (any(values(r_diff_cell_minus_chronic)[, 1] < 0, na.rm = TRUE)) {
  stop("Difference raster has negative values; chronic cell counts should not exceed 2025 cell counts.")
}

cat("Validation checks passed: dataframe-derived values match all map rasters.\n")

# Figure 4, panels B-D: regional zoom windows of the cell-level map
world_zoom <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

r_zoom <- r_cell
names(r_zoom) <- "exposure"
r_zoom_cbr <- r_zoom
vmin_zoom <- 0
vmax_zoom <- 920

zoom_break_vals <- c(0, 100, 300, 500, 700, 900)
zoom_breaks <- zoom_break_vals

if (!exists("cm.cols", mode = "function")) cm.cols <- cm.cols1
zoom_pal <- cm.cols(100)

zoom_theme <- theme_void() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(0, 0, 0, 0)
  )

make_zoom_panel <- function(rast_cbr, xlim, ylim, show_legend = FALSE, panel_label = NULL) {
  r_crop <- terra::crop(rast_cbr, terra::ext(xlim[1], xlim[2], ylim[1], ylim[2]))
  p <- ggplot() +
    geom_sf(data = world_zoom, fill = "#d3d3d3", color = NA, inherit.aes = FALSE) +
    tidyterra::geom_spatraster(data = r_crop, aes(fill = exposure)) +
    scale_fill_gradientn(
      colors = zoom_pal,
      name = "Exposed\nspecies\nper\npixel",
      na.value = "white",
      limits = c(vmin_zoom, vmax_zoom),
      oob = scales::squish,
      breaks = zoom_breaks,
      labels = zoom_break_vals,
      guide = guide_colorbar(title.position = "top", label.position = "bottom")
    ) +
    geom_sf(data = world_zoom, fill = NA, color = "grey30", linewidth = 0.2, inherit.aes = FALSE) +
    coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
    zoom_theme +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      panel.border = element_rect(fill = NA, color = "black", linewidth = 0.75),
      plot.background = element_rect(fill = "white", color = NA)
    )

  if (show_legend) {
    p <- p + theme(
      legend.position = c(0.02, 0.067),
      legend.justification = c(0, 0),
      legend.direction = "vertical",
      legend.title = element_text(size = 17, face = "bold"),
      legend.text = element_text(size = 15),
      legend.key.width = grid::unit(17, "pt"),
      legend.key.height = grid::unit(24, "pt"),
      legend.background = element_rect(fill = "white", color = NA),
      legend.margin = margin(3, 3, 3, 3),
      legend.box.margin = margin(0, 0, 0, 0)
    ) + guides(fill = guide_colorbar(barwidth = grid::unit(17, "pt"), barheight = grid::unit(153, "pt")))
  } else {
    p <- p + theme(legend.position = "none")
  }

  if (!is.null(panel_label) && nzchar(panel_label)) {
    p <- p +
      labs(tag = panel_label) +
      theme(
        plot.tag.position = c(0.05, 0.94),
        plot.tag = element_text(face = "bold", size = 20)
      )
  }

  p
}

render_zoom_panel <- function(plot_obj, file, width_px, height_px, dpi = 200) {
  ggsave(
    filename = file,
    plot = plot_obj + theme(plot.margin = margin(0, 0, 0, 0)),
    width = width_px / dpi,
    height = height_px / dpi,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

render_top_panel_png <- function(plot_obj, file, width_px, height_px, dpi = 200) {
  ggsave(
    filename = file,
    plot = plot_obj + theme(plot.margin = margin(0, 0, 0, 0)),
    width = width_px / dpi,
    height = height_px / dpi,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

left_xlim <- c(-82, -60)
left_ylim <- c(-46, 12)
bot_cd_xlim <- c(4, 36)
bot_cd_ylim <- c(-8.7, 13.7)
bot3_xlim <- c(113.6, 156.4)
bot3_ylim <- c(-20.25, 4.25)

total_w <- 3200
total_h <- 1800
left_w <- round(total_w * 0.20)
right_w <- total_w - left_w
top_h <- round(total_h * (2.7 / 3.7))
bot_h <- total_h - top_h
gap_px <- 6
bot_panel_w <- right_w - 2 * gap_px
bot_w1 <- floor(bot_panel_w / 3)
bot_w2 <- floor(bot_panel_w / 3)
bot_w3 <- bot_panel_w - bot_w1 - bot_w2
bot_w_cd <- bot_w1 + gap_px + bot_w2
bot_w_e <- bot_w3

fit_extent_to_ratio <- function(xlim, ylim, target_ratio) {
  xr <- diff(xlim)
  yr <- diff(ylim)
  cur_ratio <- xr / yr
  xmid <- mean(xlim)
  ymid <- mean(ylim)

  if (cur_ratio < target_ratio) {
    new_xr <- yr * target_ratio
    xlim <- c(xmid - new_xr / 2, xmid + new_xr / 2)
  } else if (cur_ratio > target_ratio) {
    new_yr <- xr / target_ratio
    ylim <- c(ymid - new_yr / 2, ymid + new_yr / 2)
  }

  list(xlim = xlim, ylim = ylim)
}

left_fit <- fit_extent_to_ratio(left_xlim, left_ylim, left_w / total_h)
bot_cd_fit <- fit_extent_to_ratio(bot_cd_xlim, bot_cd_ylim, bot_w_cd / bot_h)
bot3_fit <- fit_extent_to_ratio(bot3_xlim, bot3_ylim, bot_w_e / bot_h)

p_left_zoom <- make_zoom_panel(r_zoom_cbr, left_fit$xlim, left_fit$ylim, show_legend = TRUE, panel_label = "B")
p_bot_cd_zoom <- make_zoom_panel(r_zoom_cbr, bot_cd_fit$xlim, bot_cd_fit$ylim, show_legend = FALSE, panel_label = "C")
p_bot3_zoom <- make_zoom_panel(r_zoom_cbr, bot3_fit$xlim, bot3_fit$ylim, show_legend = FALSE, panel_label = "D")

tmp_left_zoom <- tempfile(fileext = ".png")
tmp_bot_cd_zoom <- tempfile(fileext = ".png")
tmp_bot3_zoom <- tempfile(fileext = ".png")

render_zoom_panel(p_left_zoom, tmp_left_zoom, left_w, total_h)
render_zoom_panel(p_bot_cd_zoom, tmp_bot_cd_zoom, bot_w_cd, bot_h)
render_zoom_panel(p_bot3_zoom, tmp_bot3_zoom, bot_w_e, bot_h)

img_left_zoom <- magick::image_read(tmp_left_zoom) |>
  magick::image_background("white", flatten = TRUE)
img_top_main <- magick::image_blank(width = right_w, height = top_h, color = "white") |>
  magick::image_background("white", flatten = TRUE)
img_bot_cd_zoom <- magick::image_read(tmp_bot_cd_zoom) |>
  magick::image_background("white", flatten = TRUE)
img_bot3_zoom <- magick::image_read(tmp_bot3_zoom) |>
  magick::image_background("white", flatten = TRUE)

gap_img <- magick::image_blank(width = gap_px, height = bot_h, color = "white")

bot_row_img <- magick::image_append(c(img_bot_cd_zoom, gap_img, img_bot3_zoom), stack = FALSE)
right_col_img <- magick::image_append(c(img_top_main, bot_row_img), stack = TRUE)
fig4_style_composite <- magick::image_append(c(img_left_zoom, right_col_img), stack = FALSE)
fig4_style_composite <- magick::image_border(fig4_style_composite, color = "white", geometry = "24x24")

subpanel_png <- file.path(out_dir_figures, paste0("Fig4_BCD_CellExposureRichness_", target_year, "_", RUN_VERSION, ".png"))
magick::image_write(fig4_style_composite, path = subpanel_png, format = "png")

cat("Saved Figure 4 panels B-D:\n")
cat(" - ", subpanel_png, "\n", sep = "")

# Figure 4, panel A: ecoregion map
world_map_sf <- if (exists("world.shp") && inherits(world.shp, "sf")) {
  suppressWarnings(st_crop(world.shp, xmin = -180, xmax = 180, ymin = -61, ymax = 86))
} else {
  warning("world.shp is not available; using rnaturalearth fallback for standalone ecoregion map.")
  rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
    st_crop(xmin = -180, xmax = 180, ymin = -61, ymax = 86)
}
world_map_sf <- suppressWarnings(st_make_valid(world_map_sf))
world_land_sf <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf") %>%
  st_crop(xmin = -180, xmax = 180, ymin = -61, ymax = 86)

# colour anchors stretch the yellow band up (about 30% of the scale) and pull dark red down
# (about 70%), ending in near-black so the top ecoregions stay distinguishable
pal_4a <- scales::gradient_n_pal(
  c("steelblue4", "steelblue1", "gold", "red1", "red4", "#33000A"),
  values = c(0, 0.07, 0.30, 0.50, 0.70, 1)
)(seq(0, 1, length.out = 100))

p_eco_2025_fig4_style <- ggplot() +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -60, ymax = 85),
    fill = "#cfe8f8",
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_sf(data = world_land_sf, fill = "grey92", color = "grey70", linewidth = 0.1) +
  tidyterra::geom_spatraster(data = r_eco, aes(fill = value), maxcell = 5e6, alpha = 0.92) +
  geom_sf(data = world_map_sf, fill = NA, color = "grey45", linewidth = 0.15) +
  scale_fill_gradientn(
    colors = pal_4a,
    name = "Exposed species\nper ecoregion",
    limits = c(0, max(1, terra::global(r_eco, "max", na.rm = TRUE)[1, 1])),
    oob = squish,
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  annotate("text", x = -178, y = 84.5, label = "A", hjust = 0, vjust = 1, fontface = "bold", size = 7) +
  labs(title = paste0("Ecoregion exposure richness (", target_year, ")")) +
  common_map_theme +
  theme(
    panel.background = element_rect(fill = "#cfe8f8", color = NA),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.ontop = FALSE
  )

map_b_fig4style_png <- file.path(out_dir_figures, paste0("Fig4_A_EcoregionExposureRichness_", target_year, "_", RUN_VERSION, ".png"))
map_b_fig4style_pdf <- sub("\\.png$", ".pdf", map_b_fig4style_png)

ggsave(map_b_fig4style_png, p_eco_2025_fig4_style, width = 11, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(map_b_fig4style_pdf, p_eco_2025_fig4_style, width = 11, height = 6.5, units = "in", bg = "white")

cat("Saved Figure 4 panel A:\n")
cat(" - ", map_b_fig4style_png, "\n", sep = "")
cat(" - ", map_b_fig4style_pdf, "\n", sep = "")

# Figure S5: chronic exposure per ecoregion
p_chronic_eco_2025_fig4_style <- ggplot() +
  geom_rect(
    aes(xmin = -180, xmax = 180, ymin = -60, ymax = 85),
    fill = "#cfe8f8",
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_sf(data = world_land_sf, fill = "grey92", color = "grey70", linewidth = 0.1) +
  tidyterra::geom_spatraster(data = r_chronic_eco, aes(fill = value), maxcell = 5e6, alpha = 0.92) +
  geom_sf(data = world_map_sf, fill = NA, color = "grey45", linewidth = 0.15) +
  scale_fill_gradientn(
    colors = my_palette,
    name = paste0("Chronic species\n(>=", chronic_min_years, " years)\nper ecoregion"),
    limits = c(0, max(1, terra::global(r_chronic_eco, "max", na.rm = TRUE)[1, 1])),
    oob = squish,
    na.value = "transparent"
  ) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 85), expand = FALSE) +
  annotate("text", x = -178, y = 84.5, label = "D", hjust = 0, vjust = 1, fontface = "bold", size = 7) +
  labs(title = paste0("Chronic ecoregion exposure richness (", target_year, ")")) +
  common_map_theme +
  theme(
    panel.background = element_rect(fill = "#cfe8f8", color = NA),
    panel.grid.major = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.ontop = FALSE
  )

map_d_fig4style_png <- file.path(out_dir_figures, paste0("FigS5_ChronicEcoregionRichness_GE", chronic_min_years, "yrs_", target_year, "_", RUN_VERSION, ".png"))
map_d_fig4style_pdf <- sub("\\.png$", ".pdf", map_d_fig4style_png)

ggsave(map_d_fig4style_png, p_chronic_eco_2025_fig4_style, width = 11, height = 6.5, units = "in", dpi = 300, bg = "white")
ggsave(map_d_fig4style_pdf, p_chronic_eco_2025_fig4_style, width = 11, height = 6.5, units = "in", bg = "white")

cat("Saved Figure S5:\n")
cat(" - ", map_d_fig4style_png, "\n", sep = "")
cat(" - ", map_d_fig4style_pdf, "\n", sep = "")
