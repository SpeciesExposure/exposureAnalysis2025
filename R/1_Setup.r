# 1_Setup.r
# Loads packages, reads the run settings and defines the directories used by
# every other script.
#
# Run settings. Set these before sourcing this script to change them:
#   END_HISTORICAL_YEAR  last year of the historical baseline (default 2022)
#   RUN_NAME             label for the run; names the output folders
#                        (default "V8", the run reported in the main text)
#
# The baseline sensitivity runs in Fig. S2 used
#   END_HISTORICAL_YEAR = 2015L; RUN_NAME = "V8_2015"
#   END_HISTORICAL_YEAR = 2009L; RUN_NAME = "V8_2009"
# All runs share the same climate tables, species ranges and species metadata.

library(tidyverse)
library(sf)
library(terra)
library(qs)
library(qs2)
library(parallel)
library(tools)
library(ggplot2)
library(gridExtra)
library(scales)
library(tidyterra)

# reader that accepts both the legacy qs format and the qs2 format
# Both start 0x0B 0x0E; byte 4 is 0xC1 for qs2, 0x0C for legacy qs
qread <- function(file, ...) {
  magic <- readBin(file, raw(), n = 4)
  if (length(magic) >= 4 && magic[4] == as.raw(0xC1)) {
    qs2::qs_read(file, ...)
  } else {
    qs::qread(file, ...)
  }
}
# always write new files in qs2 format
qsave <- function(x, file, ...) qs2::qs_save(x, file, ...)

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))

source(file.path(PROJECT_PATHS$src_r, '2_Functions.r'))
if (!exists('mc.cores')) mc.cores=max(1, parallel::detectCores() - 1)
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")

# ========== RUN SETTINGS ==========
RUN_VERSION <- 'V8'
if (!exists('END_HISTORICAL_YEAR')) END_HISTORICAL_YEAR <- 2022L
END_HISTORICAL_YEAR <- as.integer(END_HISTORICAL_YEAR)

runTag <- if (exists("RUN_NAME") && is.character(RUN_NAME) && length(RUN_NAME) == 1 && !is.na(RUN_NAME) && nzchar(trimws(RUN_NAME))) {
	trimws(RUN_NAME)
} else if (END_HISTORICAL_YEAR == 2022L) {
	RUN_VERSION
} else {
	paste0(RUN_VERSION, '_', END_HISTORICAL_YEAR)
}
RUN_NAME <- runTag
run_tag_lower <- tolower(runTag)

year='2025'

# ========== DIRECTORIES ==========
dataDir <- file.path(EXTERNAL_PATHS$exposure_2025, runTag)
metaDir <- EXTERNAL_PATHS$metadata_v8
envDir <- EXTERNAL_PATHS$env_v8

baseDir <- file.path(PROJECT_PATHS$output_base, run_tag_lower)
figuresDir <- file.path(baseDir, 'figures')
plotDir <- file.path(figuresDir, 'main')
figProposedMapDir <- file.path(figuresDir, 'proposed_maps')
intDir <- file.path(baseDir, paste0('Int_', RUN_VERSION))
tablesDir <- file.path(baseDir, 'tables')

ensure_dirs(
	dataDir,
	file.path(dataDir, 'variable_outputs'),
	metaDir,
	envDir,
	EXTERNAL_PATHS$ranges_expert,
	plotDir,
	figProposedMapDir,
	intDir,
	PROJECT_PATHS$int_v8,
	tablesDir
)

# Sensitivity runs use the same species ranges as the main run
is_sensitivity_run <- toupper(runTag) != 'V8'
if (is_sensitivity_run) {
	v8_range_root <- file.path(EXTERNAL_PATHS$exposure_v8, 'spRangeTables_Prepped')
	run_range_root <- file.path(dataDir, 'spRangeTables_Prepped')
	if (!dir.exists(run_range_root) && dir.exists(v8_range_root)) {
		ok <- file.symlink(v8_range_root, run_range_root)
		if (!ok) {
			stop('Could not create symlink for spRangeTables_Prepped: ', run_range_root, ' -> ', v8_range_root)
		}
	}
}

KEY_FILES$land_template <- file.path(intDir, 'landTemplate.tif')
KEY_FILES$single_max_exposure_v1 <- file.path(intDir, 'singleMaxExposure_v1.qs')
KEY_FILES$single_max_exposure_2025 <- file.path(intDir, 'singleMaxExposure2025_v1.qs')
KEY_FILES$all_cell_exposure <- file.path(intDir, 'AllCellExposureSpXVar.qs')

# ========== WORLD BORDERS ==========
sf::sf_use_s2(FALSE)
world.shp=st_read(INPUTS$world_borders, quiet = TRUE); world.shp <- st_make_valid(world.shp)
world.shp=st_crop(world.shp,xmin=-180,xmax=180,ymin=-61,ymax=86)
sf::sf_use_s2(TRUE)

# ========== LAND TEMPLATE ==========
# The 0.25 degree land grid. It is written by 3b_ClimatePrep.r for the main
# run. Sensitivity runs copy it from the main run.
template_file <- file.path(intDir, 'landTemplate.tif')
if (!file.exists(template_file)) {
	main_template <- file.path(PROJECT_PATHS$int_v8, 'landTemplate.tif')
	if (file.exists(main_template)) file.copy(main_template, template_file)
}
if (file.exists(template_file)) {
	template=rast(template_file)
} else {
	message('landTemplate.tif not found in ', intDir, '. Run 3b_ClimatePrep.r to create it.')
}
