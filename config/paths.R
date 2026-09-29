# ============================================================================
# PATHS
# ============================================================================
# Where the analysis finds its inputs and writes its outputs.
# Run every script from the repository root, so that getwd() is this repository.
#
# To use your own locations, create config/paths_local.R and redefine
# INPUTS and/or WORK_DIR there. That file is read at the end of this section
# and is not tracked by git.
# ============================================================================

ROOT_PATH <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(ROOT_PATH, "config", "paths.R"))) {
  stop("Run scripts from the repository root (the folder that contains config/paths.R).")
}

# Working directory for large files: climate tables, species range tables and
# per-species exposure outputs. Allow tens of GB.
WORK_DIR <- file.path(ROOT_PATH, "work")

# Input data. See README.md for where to obtain each one.
INPUTS <- list(
  # ERA5 monthly means: temp2m.nc (2 m temperature), precip.nc (total
  # precipitation), ERA5_lsm.nc (land-sea mask)
  era5_dir        = file.path(ROOT_PATH, "inputs", "ERA5"),
  # Area of Habitat maps (Lumbierres et al. 2022): folders AoH_Mammals,
  # AoH_Birds_Nonmigratory, AoH_Birds_Migratory
  aoh_dir         = file.path(ROOT_PATH, "inputs", "AoH"),
  # IUCN Red List range polygons: AMPHIBIANS/AMPHIBIANS_PART1.shp, ..., REPTILES/REPTILES_PART2.shp
  iucn_maps_dir   = file.path(ROOT_PATH, "inputs", "IUCN_maps"),
  # IUCN Red List tables: one folder per group (Amphibians, Birds1, Birds2,
  # Mammals, Reptiles), each with assessments.csv and taxonomy.csv
  iucn_tables_dir = file.path(ROOT_PATH, "inputs", "IUCN_tables"),
  # RESOLVE Ecoregions 2017: Ecoregions2017.shp, Ecoregions2017_v2.tif (the
  # ECO_ID rasterized to the 0.25 degree grid), ecoregionMetadata.rds (the
  # shapefile attribute table)
  ecoregions_dir  = file.path(ROOT_PATH, "inputs", "Ecoregions2017"),
  # World country borders used as a map overlay (TM_WORLD_BORDERS_SIMPL-0.3)
  world_borders   = file.path(ROOT_PATH, "inputs", "TM_WORLD_BORDERS_SIMPL-0.3", "TM_WORLD_BORDERS_SIMPL-0.3.shp")
)

if (file.exists(file.path(ROOT_PATH, "config", "paths_local.R"))) {
  source(file.path(ROOT_PATH, "config", "paths_local.R"))
}
WORK_DIR <- path.expand(WORK_DIR)
INPUTS <- lapply(INPUTS, path.expand)

# ============================================================================
# PROJECT PATHS (inside this repository)
# ============================================================================
PROJECT_PATHS <- list(
  root        = ROOT_PATH,
  config      = file.path(ROOT_PATH, "config"),
  src_r       = file.path(ROOT_PATH, "R"),
  figures_src = file.path(ROOT_PATH, "figures"),
  tables_src  = file.path(ROOT_PATH, "tables"),
  data        = file.path(ROOT_PATH, "data"),
  output_base = file.path(ROOT_PATH, "output"),
  output_v8   = file.path(ROOT_PATH, "output", "v8"),
  int_v8      = file.path(ROOT_PATH, "output", "v8", "Int_V8")
)

# ============================================================================
# EXTERNAL PATHS (inside WORK_DIR)
# ============================================================================
EXTERNAL_PATHS <- list(
  exposure_2025 = WORK_DIR,
  exposure_v8   = file.path(WORK_DIR, "V8"),
  env_v8        = file.path(WORK_DIR, "V8", "env_V8"),
  metadata_v8   = file.path(WORK_DIR, "V8", "metadata_V8"),
  ranges_expert = file.path(WORK_DIR, "V8", "spRangeTables_Prepped", "Expert")
)

# ============================================================================
# KEY FILES
# ============================================================================
KEY_FILES <- list(
  land_template            = file.path(PROJECT_PATHS$int_v8, "landTemplate.tif"),
  single_max_exposure_v1   = file.path(PROJECT_PATHS$int_v8, "singleMaxExposure_v1.qs"),
  single_max_exposure_2025 = file.path(PROJECT_PATHS$int_v8, "singleMaxExposure2025_v1.qs"),
  all_cell_exposure        = file.path(PROJECT_PATHS$int_v8, "AllCellExposureSpXVar.qs"),
  sp_attributes_v8         = file.path(EXTERNAL_PATHS$metadata_v8, "spAttributes_v8.qs"),
  temp_12                  = file.path(EXTERNAL_PATHS$env_v8, "temp__12.qs"),

  ecoregion_raster         = file.path(INPUTS$ecoregions_dir, "Ecoregions2017_v2.tif"),
  ecoregion_shp            = file.path(INPUTS$ecoregions_dir, "Ecoregions2017.shp"),
  ecoregion_metadata       = file.path(INPUTS$ecoregions_dir, "ecoregionMetadata.rds"),

  # Small inputs shipped in data/
  body_mass_quartiles      = file.path(PROJECT_PATHS$data, "body_mass_quartiles.csv"),
  oni                      = file.path(PROJECT_PATHS$data, "oni_v5_noaa_cpc_2026-09-15.txt")
)

ensure_dirs <- function(...) {
  dirs <- list(...)
  for (d in dirs) {
    if (!dir.exists(d)) {
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
    }
  }
  invisible(TRUE)
}
