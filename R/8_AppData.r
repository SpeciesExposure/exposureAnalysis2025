# 8_AppData.r
# Builds the data files shipped with the exposureApp Shiny app
# (https://github.com/SpeciesExposure/exposureApp).
# 
# propExposed is the fraction of a species' range exposed to any variable in a
# year. Each range cell is counted once even if several variables exceed their
# limit there. Species-years below 1% are left out of the app data.
# 
# Inputs:  AllCellExposureSpXVar.qs, spAttributes_v8.qs
# Outputs: files in app_extdata_dir

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))
if (!exists('intDir')) source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))

# folder that receives the app data
app_extdata_dir <- file.path(PROJECT_PATHS$output_base, 'exposureApp_extdata')

allExp_file <- file.path(intDir, 'AllCellExposureSpXVar.qs')
if (!file.exists(allExp_file)) {
  stop("Missing allExp file: ", allExp_file, ". Run 6_MakeDF.r to create it.")
}
allExp=qs_read(allExp_file)
allExp <- allExp %>%
  group_by(spName, year) %>%
  # count each range cell at most once per year
  mutate(propExposed = dplyr::n_distinct(cell) / rangeSize) %>%
  ungroup()

sp_meta_file <- KEY_FILES$sp_attributes_v8

spMeta=qs_read(sp_meta_file)  %>%
  dplyr::select(spName,orderName,familyName,redlistCategory,group) %>%
  unique
aa=allExp %>% left_join(spMeta,by='spName')
aa=aa %>%
  dplyr::select(cell, year, var, spName, group = group.y, propExposed, orderName, familyName, redlistCategory) %>%
  mutate(propExposed=round(propExposed,2))

# species that were never exposed are added so that every study species appears in the app
exposed_sp <- unique(aa$spName)
zero_exp_sp <- spMeta %>%
  filter(!spName %in% exposed_sp) %>%
  mutate(cell = NA_integer_, year = NA_integer_, var = NA_character_,
         propExposed = 0) %>%
  dplyr::select(cell, year, var, spName, group, propExposed, orderName, familyName, redlistCategory)
bb <- bind_rows(aa, zero_exp_sp)
qs_save(bb, file = file.path(intDir, 'allExpForShiny.qs'))

if (!dir.exists(app_extdata_dir)) dir.create(app_extdata_dir, recursive = TRUE, showWarnings = FALSE)
app_year_dir <- file.path(app_extdata_dir, 'allExpForShiny_by_year')
if (!dir.exists(app_year_dir)) dir.create(app_year_dir, recursive = TRUE, showWarnings = FALSE)
unlink(list.files(app_year_dir, pattern = '^allExpForShiny_.*\\.qs$', full.names = TRUE))

# Keep species-years whose range-exposed fraction is at least 1%, using the
# unrounded fraction
app_keep_sp_year <- allExp %>%
  distinct(spName, year, propExposed) %>%
  filter(!is.na(year), propExposed >= 0.01) %>%
  dplyr::select(spName, year)

app_filtered <- bb %>%
  semi_join(app_keep_sp_year, by = c('spName', 'year'))

# rangeSize = number of grid cells in the species' range
sp_range <- qs_read(sp_meta_file) %>%
  dplyr::select(spName, rangeSize) %>%
  filter(!is.na(rangeSize)) %>%
  distinct(spName, .keep_all = TRUE)
app_species_meta <- app_filtered %>%
  dplyr::select(spName, group, orderName, familyName, redlistCategory) %>%
  mutate(
    spName = as.character(spName),
    group = as.character(group),
    orderName = str_to_title(orderName),
    familyName = str_to_title(familyName),
    redlistCategory = as.character(redlistCategory)
  ) %>%
  distinct(spName, group, orderName, familyName, redlistCategory) %>%
  left_join(sp_range, by = 'spName')

app_rows <- app_filtered %>%
  filter(!is.na(year), !is.na(var), !is.na(cell)) %>%
  mutate(
    spName = as.character(spName),
    var = as.character(var),
    group = as.character(group),
    orderName = str_to_title(orderName),
    familyName = str_to_title(familyName),
    redlistCategory = as.character(redlistCategory),
    propExposed = as.numeric(propExposed)
  ) %>%
  arrange(year, spName, var, cell)

app_years <- sort(unique(app_rows$year))
app_year_files <- setNames(character(length(app_years)), as.character(app_years))
for (yr in app_years) {
  year_file <- file.path(app_year_dir, paste0('allExpForShiny_', yr, '.qs'))
  qs_save(app_rows %>% filter(year == yr), file = year_file)
  app_year_files[as.character(yr)] <- file.path('allExpForShiny_by_year', basename(year_file))
}

cell_trend_df <- app_rows %>%
  distinct(spName, cell, year) %>%
  count(cell, year, name = 'n_sp')

sp_trend_df <- app_rows %>%
  group_by(spName, var, year) %>%
  summarize(mean_prop = mean(propExposed, na.rm = TRUE), .groups = 'drop')

species_year_cells <- app_rows %>%
  distinct(year, spName, var, cell) %>%
  arrange(year, spName, var, cell) %>%
  group_by(year, spName, var) %>%
  summarize(cells = list(as.integer(cell)), .groups = 'drop')

qs_save(app_species_meta, file = file.path(app_extdata_dir, 'allExpForShiny_species.qs'))
qs_save(app_rows, file = file.path(app_extdata_dir, 'allExpForShiny.qs'))
qs_save(
  list(
    version = 1L,
    cell_trend_df = cell_trend_df,
    sp_trend_df = sp_trend_df
  ),
  file = file.path(app_extdata_dir, 'allExpForShiny_trends_cache_v1.qs')
)
qs_save(species_year_cells, file = file.path(app_extdata_dir, 'allExpForShiny_species_year_cells_v1.qs'))
qs_save(
  list(
    version = 1L,
    years_avail = app_years,
    year_files = as.list(app_year_files),
    species_file = 'allExpForShiny_species.qs',
    trends_cache_file = 'allExpForShiny_trends_cache_v1.qs',
    avail_cells_all = sort(unique(app_rows$cell))
  ),
  file = file.path(app_extdata_dir, 'allExpForShiny_manifest_v1.qs')
)

#=================================================
# Range cells of mammals and birds, used by the app as a map background
mb <- allExp |>
  filter(group %in% c("Mammals", "Birds"), !is.na(cell)) |>
  dplyr::select(spName, cell) |>
  distinct()
range_cells_mb <- split(as.integer(mb$cell), mb$spName)
range_cells_mb <- lapply(range_cells_mb, unique)
qs_save(range_cells_mb, file.path(app_extdata_dir, 'range_cells_mb.qs'))
