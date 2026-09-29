# Table1_Table2.r
# Table 1: the 50 family x ecoregion x variable combinations with the most
# species exposed in 2025.
# Table 2: species exposed for the first time in 2025.
# 
# Inputs:  singleMaxExposure2025_v1.qs, singleMaxExposure_v1.qs, spAttributes_v8.qs
# Outputs: CSV files in tablesDir

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))
if (!exists('intDir')) source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))

#=================================================
# Table 1

aa=qs_read(paste0(intDir,'/singleMaxExposure2025_v1.qs'))  %>%
  filter(year == 2025, pExp >= 0.25) %>%
  mutate(
    orderName  = str_to_title(orderName),
    familyName = str_to_title(familyName)
  )%>% mutate(var = recode(var,
    "temp__12_up"       = "Temp 12mo High",
    "temp__3__max_up"   = "Temp 3mo High",
    "temp__3__min_lo"   = "Temp 3mo Low",
    "precip__12_up"     = "Precip 12mo High",
    "precip__12_lo"     = "Precip 12mo Low",
    "precip__3__max_up" = "Precip 3mo High",
    "precip__3__min_lo" = "Precip 3mo Low"
  ))
e=st_read(KEY_FILES$ecoregion_shp)
# number of species per ecoregion, variable and family
bb2=aa %>%
  filter(!is.na(ECO_ID)) %>%
  group_by(ECO_ID,var,familyName) %>%
  summarize(Nspecies=n_distinct(spName), .groups = 'drop') %>%
  left_join(e) %>%
  arrange(desc(Nspecies)) %>%
  ungroup

family_eco_top50 <- bb2 %>%
  st_drop_geometry() %>%
  dplyr::select(familyName, ECO_NAME, var, Nspecies) %>%
  head(50)

#=================================================
# Newly exposed species in 2025 table
# First year of exposure is derived from singleMaxExposure_v1.qs

single_max_file <- file.path(intDir, "singleMaxExposure_v1.qs")
if (!file.exists(single_max_file)) {
  stop("Missing single-max exposure file: ", single_max_file, ". Run 6_MakeDF.r first.")
}

ex11_raw <- qs_read(single_max_file)
if (!"spName" %in% names(ex11_raw) && "sp" %in% names(ex11_raw)) {
  ex11_raw <- ex11_raw %>% rename(spName = sp)
}

ex11 <- ex11_raw %>%
  arrange(spName, year, desc(pExp), var) %>%
  group_by(spName, year) %>%
  slice(1) %>%
  ungroup()

new_2025 <- ex11 %>%
  group_by(spName) %>%
  summarise(first_year = min(year), .groups = "drop") %>%
  filter(first_year == 2025)

sp_attributes_file <- file.path(metaDir, paste0("spAttributes_", tolower(RUN_VERSION), ".qs"))
if (!file.exists(sp_attributes_file)) {
  stop("Missing metadata file: ", sp_attributes_file)
}

sp_attr <- qs_read(sp_attributes_file)
if (!"spName" %in% names(sp_attr) && "sp" %in% names(sp_attr)) {
  sp_attr <- sp_attr %>% rename(spName = sp)
}

if (!"REALM" %in% names(sp_attr) && "realm" %in% names(sp_attr)) {
  sp_attr <- sp_attr %>% mutate(REALM = realm)
}

summarize_one <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(trimws(x))]
  if (!length(x)) return(NA_character_)
  paste(sort(unique(x)), collapse = "; ")
}

format_var <- function(v) {
  dplyr::recode(
    trimws(as.character(v)),
    "temp__12_up" = "Temp annual max",
    "temp__3__max_up" = "Temp summer max",
    "temp__3__min_lo" = "Temp winter min",
    "precip__12_up" = "Precip annual max",
    "precip__12_lo" = "Precip annual min",
    "precip__3__max_up" = "Precip 3 mo max",
    "precip__3__min_lo" = "Precip 3 mo min",
    .default = gsub("__", " ", trimws(as.character(v)))
  )
}

is_least_concern <- function(x) {
  x_norm <- gsub("[^a-z]", "", tolower(trimws(as.character(x))))
  x_norm %in% c("leastconcern", "lowerriskleastconcern", "lrlc")
}

meta_2025 <- sp_attr %>%
  filter(spName %in% new_2025$spName) %>%
  group_by(spName) %>%
  summarise(
    `Taxon group` = dplyr::first(na.omit(group)),
    `IUCN status` = dplyr::first(na.omit(redlistCategory)),
    Realms = summarize_one(REALM),
    .groups = "drop"
  )

exp_2025 <- ex11 %>%
  filter(year == 2025, spName %in% new_2025$spName) %>%
  group_by(spName) %>%
  slice_max(pExp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    spName,
    `Prop. Exposed` = paste0(round(100 * pExp, 0), "%"),
    Variable = format_var(var)
  )

tab_2025 <- new_2025 %>%
  dplyr::select(spName) %>%
  left_join(meta_2025, by = "spName") %>%
  left_join(exp_2025, by = "spName") %>%
  mutate(
    Scientific = gsub("_", " ", spName),
    Realms = ifelse(is.na(Realms) | !nzchar(Realms), "Unknown", Realms),
    `Taxon group` = ifelse(is.na(`Taxon group`) | !nzchar(`Taxon group`), "Unknown", `Taxon group`),
    `IUCN status` = ifelse(is.na(`IUCN status`) | !nzchar(`IUCN status`), "Unknown", `IUCN status`),
    `Prop. Exposed` = ifelse(is.na(`Prop. Exposed`) | !nzchar(`Prop. Exposed`), "NA", `Prop. Exposed`),
    Variable = ifelse(is.na(Variable) | !nzchar(Variable), "NA", Variable),
    prop_exposed_num = suppressWarnings(as.numeric(gsub("%", "", `Prop. Exposed`)))
  ) %>%
  arrange(desc(prop_exposed_num), Scientific) %>%
  dplyr::select(Scientific, Realms, `Taxon group`, `IUCN status`, `Prop. Exposed`, Variable, prop_exposed_num) %>%
  dplyr::select(Scientific, Realms, `Taxon group`, `IUCN status`, `Prop. Exposed`, Variable)

tab_2025_non_lc <- tab_2025 %>%
  filter(!is_least_concern(`IUCN status`))

#=================================================
# Write tables

dir.create(tablesDir, recursive = TRUE, showWarnings = FALSE)
run_tag_safe <- gsub("[^A-Za-z0-9_]", "_", tolower(runTag))
write.csv(family_eco_top50, file.path(tablesDir, paste0("family_by_ecoregion_top50_table_", run_tag_safe, ".csv")), row.names = FALSE)
write.csv(tab_2025, file.path(tablesDir, paste0("newly_exposed_species_2025_table_", run_tag_safe, ".csv")), row.names = FALSE)
write.csv(tab_2025_non_lc, file.path(tablesDir, paste0("newly_exposed_species_2025_table_non_lc_", run_tag_safe, ".csv")), row.names = FALSE)
