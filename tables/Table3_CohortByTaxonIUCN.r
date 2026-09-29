# Table3_CohortByTaxonIUCN.r
# Builds main-text Table 3: the exposure cohorts reported in Results 3-3.3, broken down by
# taxonomic group and by pooled IUCN category, so the counts can leave the prose.
#
# Rows (cohorts), all from singleMaxExposure_v1.qs (one row per species-year with >=25% of range exposed):
#   Study species              all species with an expert range map (denominator)
#   Exposed in 2025            species with a 2025 row
#   Chronically exposed        >=6 distinct exposure years over 1990-2025
#   Exposed 2023-2025          exposed in each of 2023, 2024 and 2025
#   First exposed in 2025      earliest exposure year is 2025
#   New variable in 2025       exposed in 2025 to a variable never before recorded for that species,
#                              among species already exposed in an earlier year (i.e. excludes the row above)
# Columns: Amphibians, Reptiles, Mammals, Birds | Threatened (CR+EN+VU), Near Threatened, Least Concern,
#          Data Deficient. Species with no IUCN category count in the row total but in no IUCN column.
# Cells: n (% of the row total).
#
# Outputs: output/v8/tables/Table3_cohorts_by_taxon_iucn_V8.csv        (formatted cells)
#          output/v8/tables/Table3_cohorts_by_taxon_iucn_V8_counts.csv (raw counts)

suppressPackageStartupMessages({ library(dplyr); library(tidyr); library(qs2) })

if (!exists("RUN_VERSION")) RUN_VERSION <- "V8"
if (!exists("PROJECT_PATHS")) source(file.path("config", "paths.R"))
source(file.path(PROJECT_PATHS$src_r, "1_Setup.r"))

single_max_file <- KEY_FILES$single_max_exposure_v1
sp_attr_file <- KEY_FILES$sp_attributes_v8
range_dir <- file.path(dataDir, "spRangeTables_Prepped", "Expert")
out_dir <- file.path(baseDir, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (f in c(single_max_file, sp_attr_file)) if (!file.exists(f)) stop("Missing input: ", f)

x <- qs_read(single_max_file)
ex11 <- x %>% arrange(spName, year, desc(pExp), var) %>% group_by(spName, year) %>% slice(1) %>% ungroup()

study <- qs_read(sp_attr_file) %>%
  distinct(spName, group, redlistCategory) %>%
  filter(file.exists(file.path(range_dir, paste0(spName, ".qs"))))
if (anyDuplicated(study$spName)) stop("Species with more than one group/category row in spAttributes")
if (nrow(study) != 32345) stop("Expected 32,345 study species, got ", nrow(study))

study <- study %>% mutate(iucn4 = case_when(
  redlistCategory %in% c("Critically_Endangered", "Endangered", "Vulnerable") ~ "Threatened",
  redlistCategory == "Near_Threatened" ~ "Near Threatened",
  redlistCategory %in% c("Least_Concern", "Lower_Risk/least concern",
                         "Lower_Risk/conservation dependent", "Lower_Risk/near threatened") ~ "Least Concern",
  redlistCategory == "Data_Deficient" ~ "Data Deficient",
  TRUE ~ "Not assessed"))

sp_by_year <- split(ex11$spName, ex11$year)
first_year <- ex11 %>% group_by(spName) %>% summarise(first = min(year), .groups = "drop")
first_var_year <- ex11 %>% group_by(spName, var) %>% summarise(first = min(year), .groups = "drop")
new_var_2025 <- ex11 %>% filter(year == 2025) %>% distinct(spName, var) %>%
  inner_join(first_var_year, by = c("spName", "var")) %>% filter(first == 2025) %>% distinct(spName) %>% pull(spName)
first_2025 <- first_year %>% filter(first == 2025) %>% pull(spName)

cohorts <- list(
  "Study species" = study$spName,
  "Exposed in 2025" = unique(sp_by_year[["2025"]]),
  "Chronically exposed (>=6 years, 1990-2025)" = ex11 %>% distinct(spName, year) %>% count(spName) %>% filter(n >= 6) %>% pull(spName),
  "Exposed in each of 2023, 2024 and 2025" = Reduce(intersect, list(sp_by_year[["2023"]], sp_by_year[["2024"]], sp_by_year[["2025"]])),
  "First exposed in 2025" = first_2025,
  "Exposed to a new variable in 2025 (previously exposed species)" = setdiff(new_var_2025, first_2025))

taxa <- c("Amphibians", "Reptiles", "Mammals", "Birds")
iucn <- c("Threatened", "Near Threatened", "Least Concern", "Data Deficient")

counts <- bind_rows(lapply(names(cohorts), function(nm) {
  s <- study %>% filter(spName %in% cohorts[[nm]])
  if (nrow(s) != length(cohorts[[nm]])) stop("Cohort '", nm, "' contains species outside the study set")
  tibble(cohort = nm, total = nrow(s),
         !!!setNames(lapply(taxa, function(g) sum(s$group == g)), taxa),
         !!!setNames(lapply(iucn, function(g) sum(s$iucn4 == g)), iucn),
         not_assessed = sum(s$iucn4 == "Not assessed"))
}))

# expected numbers from the current manuscript text, for the cross-check
expected <- tribble(
  ~cohort, ~total, ~Amphibians, ~Reptiles, ~Mammals, ~Birds, ~Threatened, ~`Data Deficient`, ~`Least Concern`,
  "Exposed in 2025", 2662, 1060, 1032, 315, 255, 1083, 713, 587,
  "Chronically exposed (>=6 years, 1990-2025)", 1584, 694, 604, 163, 123, NA, 451, 244,
  "Exposed in each of 2023, 2024 and 2025", 1133, 460, 481, 89, 103, 508, 329, 194,
  "First exposed in 2025", 129, 28, 49, 28, 24, NA, 13, 62,
  "Exposed to a new variable in 2025 (previously exposed species)", 450, 185, 163, 66, 36, NA, 157, 93)
chk <- expected %>% pivot_longer(-cohort, names_to = "col", values_to = "text") %>% filter(!is.na(text)) %>%
  left_join(counts %>% pivot_longer(-cohort, names_to = "col", values_to = "data"), by = c("cohort", "col"))
bad <- chk %>% filter(text != data)
cat("Cross-check against manuscript text:", nrow(chk), "numbers compared,", nrow(bad), "mismatches\n")
if (nrow(bad) > 0) print(as.data.frame(bad))

fmt <- function(n, tot) ifelse(is.na(n), "", sprintf("%s (%.1f)", format(n, big.mark = ","), 100 * n / tot))
formatted <- counts %>% mutate(across(all_of(c(taxa, iucn)), ~ fmt(.x, total))) %>%
  mutate(total = format(total, big.mark = ",")) %>% dplyr::select(-not_assessed)

write.csv(counts, file.path(out_dir, paste0("Table3_cohorts_by_taxon_iucn_", RUN_VERSION, "_counts.csv")), row.names = FALSE)
write.csv(formatted, file.path(out_dir, paste0("Table3_cohorts_by_taxon_iucn_", RUN_VERSION, ".csv")), row.names = FALSE)
cat("Not assessed (no IUCN category) per row:", paste(counts$cohort, counts$not_assessed, sep = " = ", collapse = "; "), "\n")
print(as.data.frame(formatted))
