# 6_MakeDF.r
# Collects the per-species outputs of 5_CalcExposure.r into the data frames used
# by the figures and tables. Thresholds used throughout: quantile 0.99 over
# years then 0.99 over cells for upper limits, 0.01 and 0.01 for lower limits,
# quantile type 7.
#
# Data frame: ex11
# One row per species-year (after tie-breaking to keep the variable with max pExp).
# Includes:
# - spName: species name
# - year: exposure year
# - var: climate variable producing the max single-variable exposure in that year
# - pExp: proportion of species range exposed
# - range/threshold stats inherited from variable stats files
# - species metadata joined from spAttributes_<version>.qs (group, redlistCategory)

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))

if (!exists("RUN_VERSION") || !exists("metaDir") || !exists("dataDir") || !exists("intDir") || !exists("plotDir") || !exists("mc.cores")) {
	source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))
}

ex11_file <- file.path(intDir, "singleMaxExposure_v1.qs")
sp_attributes_file <- file.path(metaDir, paste0("spAttributes_", tolower(RUN_VERSION), ".qs"))
if (!file.exists(sp_attributes_file)) {
	stop("Missing required metadata file for RUN_VERSION ", RUN_VERSION, ": ", sp_attributes_file)
}

vars <- list.files(file.path(dataDir, 'variable_outputs')) %>% file_path_sans_ext
spMeta <- qread(sp_attributes_file) %>% dplyr::select(spName, group, redlistCategory) %>% unique

normalize_sp_name_col <- function(df) {
	if ("spName" %in% names(df)) return(df)
	if ("sp" %in% names(df)) return(df %>% rename(spName = sp))
	df
}

ensure_ex11_schema <- function(df) {
	df <- normalize_sp_name_col(df)
	defaults <- tibble(
		spName = character(),
		year = numeric(),
		pExp = numeric(),
		var = character(),
		rangeSize = numeric(),
		group = character(),
		redlistCategory = character()
	)
	for (nm in setdiff(names(defaults), names(df))) df[[nm]] <- defaults[[nm]]
	df
}

ensure_allexp_schema <- function(df) {
	df <- normalize_sp_name_col(df)
	defaults <- tibble(
		year = numeric(),
		cell = numeric(),
		rangeSize = numeric(),
		var = character(),
		spName = character(),
		group = character()
	)
	for (nm in setdiff(names(defaults), names(df))) df[[nm]] <- defaults[[nm]]
	df
}

ex1 <- lapply(vars, function(vv) {
print(vv)
ff <- list.files(file.path(dataDir, 'variable_outputs', vv, paste0(vv, '_stats')), full.names = TRUE)
d11 <- mclapply(ff, function(x) {
	a <- qread(x) %>%
	filter(near(year.q, ifelse(grepl('up', vv), .99, .01)) &
				!isRounded &
				near(cell.q, ifelse(grepl('up', vv), .99, .01)) &
				quantileType == 7) %>%
	dplyr::select(-contains('mar')) %>%
	dplyr::select(-contains('XnM'), -expPre1990) %>%
	mutate(across(X1990:X2025, function(z) { round(z / rangeSize, 3) })) %>%
	mutate(spName = basename(file_path_sans_ext(x)), var = vv) %>%
	pivot_longer(starts_with('X'), names_to = 'year', values_to = 'pExp') %>%
	filter(!year == 'X2026') %>%
	dplyr::select(-quantileIndex) %>%
	filter(pExp >= .25) %>%
	mutate(year = as.numeric(sub('X', '', year)))
	a
}, mc.cores = mc.cores) %>% Filter(is.data.frame, .) %>% bind_rows

d11 %>% left_join(spMeta, by = c('spName')) %>% dplyr::select(-thresh.val) %>% filter(!is.na(group))
}) %>% bind_rows %>%
	ensure_ex11_schema() %>%
	dplyr::select(-any_of('year.cell.q'))

ex11 <- ex1 %>%
arrange(spName, year, desc(pExp), var) %>%
group_by(spName, year) %>%
slice(1) %>%
ungroup()

qsave(ex11, file = ex11_file)

ex11 <- qread(ex11_file)

#=================================================
#=================================================
# Data frame: allExp
# One row per exposed cell event for a species-variable-year combination.
# Includes:
# - spName: species name
# - var: climate variable
# - year, cell: exposed year and cell id from exposedYearCells
# - rangeSize: species range size used for proportions
# - group: taxonomic group from spAttributes_<version>.qs
allExp_file <- file.path(intDir, "AllCellExposureSpXVar.qs")

vars_allexp <- list.files(file.path(dataDir, 'variable_outputs'))
spMeta_allexp <- qread(sp_attributes_file) %>% dplyr::select(spName, group) %>% unique

allExp <- lapply(vars_allexp, function(vv) {
print(vv)
ff <- list.files(paste0(dataDir, '/variable_outputs/', vv, '/', vv, '_marg'), full.names = TRUE)
d2 <- mclapply(ff, function(x) {
	aa <- qread(x) %>%
	filter(near(year.q, ifelse(grepl('up', vv), .99, .01)) &
				!isRounded &
				near(cell.q, ifelse(grepl('up', vv), .99, .01)) &
				quantileType == 7) %>%
	dplyr::select(exposedYearCells, rangeSize) %>%
	unnest(exposedYearCells) %>%
	mutate(var = vv, spName = basename(file_path_sans_ext(x)))

	if (nrow(aa) == 0) return(NULL)
	aa
}, mc.cores = mc.cores) %>% Filter(function(x) is.data.frame(x) && nrow(x) > 0, .) %>% bind_rows

d2 %>% left_join(spMeta_allexp, by = c('spName'))
}) %>% bind_rows %>%
	ensure_allexp_schema()

qsave(allExp, file = allExp_file)

allExp <- qread(allExp_file)

#=================================================
# Data frame: singleMaxExposure2025
# One row per species-variable-ecoregion combination for 2025 where exposure occurs.
# Includes ex11 2025 rows joined to species metadata and exposed ecoregions.
single_max_exposure_2025_file <- file.path(intDir, "singleMaxExposure2025_v1.qs")

ex12 <- ex11 %>% filter(year == 2025)

spMeta_full <- qread(sp_attributes_file)
required_spmeta_cols <- c("spName", "orderName", "familyName", "redlistCategory", "group")
missing_spmeta_cols <- setdiff(required_spmeta_cols, names(spMeta_full))
if (length(missing_spmeta_cols) > 0) {
	stop(
		"Missing required columns in ",
		sp_attributes_file,
		": ",
		paste(missing_spmeta_cols, collapse = ", ")
	)
}

first_non_na <- function(x) {
	x <- x[!is.na(x) & nzchar(trimws(as.character(x)))]
	if (length(x) == 0) return(NA)
	x[[1]]
}

spMeta_2025 <- spMeta_full %>%
	dplyr::select(all_of(required_spmeta_cols)) %>%
	group_by(spName) %>%
	summarise(
		orderName = first_non_na(orderName),
		familyName = first_non_na(familyName),
		redlistCategory = first_non_na(redlistCategory),
		group = first_non_na(group),
		.groups = "drop"
	)

a <- ex12 %>% left_join(spMeta_2025, by = "spName")

ecoregion_raster_file <- KEY_FILES$ecoregion_raster
if (!file.exists(ecoregion_raster_file)) {
	stop("Missing ecoregion raster: ", ecoregion_raster_file)
}

e.r1 <- rast(ecoregion_raster_file)
e.r1 <- crop(e.r1, template)
eco.df <- tibble(ECO_ID = values(e.r1)[, 1], cell = 1:ncell(e.r1)) %>% na.omit

allExp_2025 <- allExp %>%
	filter(year == 2025) %>%
	left_join(eco.df, by = "cell") %>%
	dplyr::select(spName, var, rangeSize, year, ECO_ID) %>%
	filter(!is.na(ECO_ID)) %>%
	distinct()

singleMaxExposure2025 <- a %>%
	left_join(allExp_2025, by = c("spName", "var", "rangeSize", "year")) %>%
	distinct(spName, var, year, ECO_ID, .keep_all = TRUE)

qsave(singleMaxExposure2025, file = single_max_exposure_2025_file)
