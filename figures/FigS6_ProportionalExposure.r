# FigS6_ProportionalExposure.r
# Supplementary Figure S6: proportion of resident vertebrate species exposed per
# 0.25-degree cell in 2025:
#   (a) number of species exposed per cell (numerator; raster from Fig4_FigS5_Maps.r)
#   (b) proportion exposed = exposed count / resident vertebrate richness
#   (c) resident vertebrate richness (denominator; stacked expert range maps)
#
# Inputs:  Int_V8/MapProducts_4_Maps/rasters/Exposure_CellRichness_2025_V8.tif  (from Fig4_FigS5_Maps.r)
#          spRangeTables_Prepped/Expert/<species>.qs (vector of land-template cell ids)
#          metadata_V8/spAttributes_v8.qs (defines the 32,345-species study set)
# Outputs: Int_V8/vertRichness_V8.tif, Int_V8/propExposed_2025_V8.tif
#          figures/main/FigS6_ProportionalExposure_V8.{pdf,png}

suppressPackageStartupMessages({
  library(dplyr)
  library(terra)
  library(sf)
  library(qs2)
  library(parallel)
})

if (!exists("RUN_VERSION")) RUN_VERSION <- "V8"
if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

target_year <- 2025
mc_cores <- 10

# ---- inputs ------------------------------------------------------------------
exp_cell_file <- file.path(intDir, "MapProducts_4_Maps", "rasters",
                           paste0("Exposure_CellRichness_", target_year, "_", RUN_VERSION, ".tif"))
range_dir <- file.path(dataDir, "spRangeTables_Prepped", "Expert")
sp_attr_file <- KEY_FILES$sp_attributes_v8
template_file <- KEY_FILES$land_template
for (f in c(exp_cell_file, sp_attr_file, template_file)) if (!file.exists(f)) stop("Missing input: ", f)
if (!dir.exists(range_dir)) stop("Missing range directory: ", range_dir)

# ---- outputs -----------------------------------------------------------------
rich_file <- file.path(intDir, paste0("vertRichness_", RUN_VERSION, ".tif"))
prop_file <- file.path(intDir, paste0("propExposed_", target_year, "_", RUN_VERSION, ".tif"))
fig_pdf <- file.path(plotDir, paste0("FigS6_ProportionalExposure_", RUN_VERSION, ".pdf"))
fig_png <- file.path(plotDir, paste0("FigS6_ProportionalExposure_", RUN_VERSION, ".png"))

template <- rast(template_file)

# ---- 1. resident richness: stack the expert range maps of the study species ----
if (!file.exists(rich_file)) {
  sp_attr <- qs_read(sp_attr_file)
  study_sp <- unique(sp_attr$spName)
  range_files <- file.path(range_dir, paste0(study_sp, ".qs"))
  has_file <- file.exists(range_files)
  cat("Study species:", length(study_sp), " with range file:", sum(has_file), "\n")
  if (sum(has_file) != 32345) stop("Expected 32,345 species with range files, got ", sum(has_file))
  range_files <- range_files[has_file]

  n_cells <- ncell(template)   # plain integer: terra objects cannot be used inside forked workers
  chunks <- split(range_files, ceiling(seq_along(range_files) / 500))
  counts <- mclapply(chunks, function(fs) {
    cells <- unlist(lapply(fs, function(f) as.integer(qs::qread(f))))   # range tables are legacy qs format
    tabulate(cells, nbins = n_cells)
  }, mc.cores = mc_cores)
  bad <- !vapply(counts, is.numeric, logical(1))
  if (any(bad)) stop("Range-file reading failed in ", sum(bad), " chunk(s): ", as.character(counts[[which(bad)[1]]]))
  rich_vals <- Reduce(`+`, counts)

  rich <- template
  values(rich) <- as.numeric(rich_vals)
  rich <- mask(rich, template)          # keep land only
  names(rich) <- "vertRichness"
  writeRaster(rich, rich_file, overwrite = TRUE)
  cat("Saved richness raster:", rich_file, "\n")
} else {
  cat("Using existing richness raster:", rich_file, "\n")
}
rich <- rast(rich_file)

# ---- 2. exposed count and proportion -----------------------------------------
expo <- rast(exp_cell_file)
if (!compareGeom(expo, rich, stopOnError = FALSE)) stop("Exposure and richness rasters are not aligned")
expo <- mask(expo, template)
expo[is.na(expo) & !is.na(template)] <- 0   # land cells with no exposed species

prop <- expo / rich
prop[is.infinite(prop)] <- NA
n_over <- global(prop > 1, "sum", na.rm = TRUE)[1, 1]
prop[prop > 1] <- NA                         # as in 2023: edge/island disagreements only
names(prop) <- "propExposed"
writeRaster(prop, prop_file, overwrite = TRUE)
cat("Saved proportion raster:", prop_file, "  (cells with count > richness set NA:", n_over, ")\n")

# ---- 3. summary statistics for the caption/text ------------------------------
land_n <- global(!is.na(template), "sum")[1, 1]
v_exp <- values(expo)[, 1]; v_rich <- values(rich)[, 1]; v_prop <- values(prop)[, 1]
ok <- !is.na(v_prop) & v_rich > 0
cat(sprintf("Land cells: %d; with >=1 exposed species: %d (%.1f%%)\n",
            land_n, sum(v_exp > 0, na.rm = TRUE), 100 * sum(v_exp > 0, na.rm = TRUE) / land_n))
qs_p <- quantile(v_prop[ok & v_exp > 0], c(.5, .9, .99), na.rm = TRUE)
cat(sprintf("Proportion exposed, cells with exposure: median %.3f, 90th %.3f, 99th %.3f, max %.3f\n",
            qs_p[1], qs_p[2], qs_p[3], max(v_prop[ok], na.rm = TRUE)))
for (th in c(.05, .10, .25, .50))
  cat(sprintf("  cells with proportion > %.2f: %d (%.2f%% of land)\n", th, sum(v_prop[ok] > th), 100 * sum(v_prop[ok] > th) / land_n))
cat(sprintf("Spearman rho (count vs proportion, exposed cells): %.3f\n",
            cor(v_exp[ok & v_exp > 0], v_prop[ok & v_exp > 0], method = "spearman")))
cat(sprintf("Spearman rho (richness vs count, exposed cells): %.3f\n",
            cor(v_rich[ok & v_exp > 0], v_exp[ok & v_exp > 0], method = "spearman")))

# ---- 4. figure: three stacked panels, as in the 2023 Fig. 1 ------------------
cnt_breaks <- c(0, 1, 5, 10, 25, 50, 100, 200, 400, 800)
cnt_cols <- c("grey93", cm.cols1(length(cnt_breaks) - 2, bias = .7))
prop_breaks <- c(0, .005, .01, .025, .05, .1, .2, .35, .5, 1)
prop_cols <- c("grey93", cm.cols1(length(prop_breaks) - 2, bias = 1.2))
rich_breaks <- c(0, 1, 5, 10, 25, 50, 100, 200, 400, 800, 1200)
rich_cols <- c("grey93", cm.cols1(length(rich_breaks) - 2, bias = .8))

draw_panel <- function(r, breaks, cols, title, label) {
  plot(r, type = "interval", breaks = breaks, col = cols,
       axes = FALSE, box = FALSE, mar = c(0.2, 0.2, 0.2, 6),
       plg = list(title = title, title.cex = .9, cex = .8))
  plot(st_geometry(world.shp), add = TRUE, lwd = .3, border = "grey40")
  mtext(label, side = 3, adj = .02, cex = 1.3, line = -1.3)
}

for (dev in c("pdf", "png")) {
  if (dev == "pdf") pdf(fig_pdf, width = 10, height = 9) else png(fig_png, width = 3000, height = 2700, res = 300)
  par(mfrow = c(3, 1), oma = c(0, 0, 0, 0))
  draw_panel(expo, cnt_breaks, cnt_cols, "# species\nexposed", "a")
  draw_panel(prop, prop_breaks, prop_cols, "proportion\nexposed", "b")
  draw_panel(rich, rich_breaks, rich_cols, "vertebrate\nrichness", "c")
  dev.off()
}
cat("Saved figure:", fig_pdf, "\n            ", fig_png, "\n")
