# Species exposure to climate extremes, 1941 to 2025: analysis code

This repository holds the R code behind "The 2025 Biodiversity Exposure Report" (Merow et al., in revision at *BioScience*). It computes, for 32,345 terrestrial vertebrate species and every year from 1941 to 2025, where in each species' range the climate was more extreme than the species had experienced during a historical baseline.

This is a fixed copy of the code for the paper. It is not under active development.

Related resources:

- Web dashboard: <https://speciesexposure.github.io/>
- Shiny app and the derived data: <https://github.com/SpeciesExposure/exposureApp>. Its file `AI_HANDOFF.md` documents the derived data in detail.

## The method in brief

1. **Climate.** ERA5 monthly 2 m temperature and total precipitation are summarized by water year on a 0.25 degree grid. Water year Y runs from October of Y-1 to September of Y.
2. **Variables.** Six climate tables give eight exposure variables.

   | Climate table | Meaning | Upper limit | Lower limit |
   |---|---|---|---|
   | `temp__12` | Mean temperature of the water year | yes | yes |
   | `temp__3__max` | Mean temperature of the warmest 3 consecutive months | yes | |
   | `temp__3__min` | Mean temperature of the coldest 3 consecutive months | | yes |
   | `precip__12` | Total precipitation of the water year | yes | yes |
   | `precip__3__max` | Precipitation of the wettest 3 consecutive months | yes | |
   | `precip__3__min` | Precipitation of the driest 3 consecutive months | | yes |

3. **Ranges.** Each species' range is the set of land cells it occupies. Mammals and birds use Area of Habitat maps. Amphibians and reptiles use IUCN range polygons.
4. **Historical limit.** For each species and variable, the 99th percentile across the baseline years is taken in every range cell, and the limit is the 99th percentile of those values across cells. Lower limits use the 1st percentile in the same way. The baseline is 1941 to 2022.
5. **Exposure.** A cell is exposed in a year when its value is beyond the limit. A species is exposed in a year when at least 25% of its range is exposed for at least one variable.

## Requirements

- R 4.2 or later
- R packages: `tidyverse`, `sf`, `terra`, `qs`, `qs2`, `gridExtra`, `scales`, `tidyterra`, `ggnewscale`, `patchwork`, `cowplot`, `magick`, `rnaturalearth`, `jsonlite`
- A machine with several cores and tens of GB of free disk space. The exposure calculation writes one file per species and variable.

## Input data

None of the input data is stored here. Download it and either place it under `inputs/` as listed in `config/paths.R`, or create `config/paths_local.R` and point `INPUTS` and `WORK_DIR` to your own locations.

| Input | Source |
|---|---|
| ERA5 monthly means, 1940 to 2025: 2 m temperature, total precipitation, land-sea mask | Copernicus Climate Data Store, <https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels-monthly-means> |
| Area of Habitat maps for mammals and birds | Lumbierres et al. (2022) |
| Range polygons for amphibians and reptiles | IUCN Red List spatial data, <https://www.iucnredlist.org/resources/spatial-data-download> |
| Assessment and taxonomy tables for all four groups | IUCN Red List |
| Ecoregions | RESOLVE Ecoregions 2017 |
| Country borders for map overlays | TM World Borders 0.3 |

The ecoregion raster `Ecoregions2017_v2.tif` is the ecoregion id, `ECO_ID`, of the RESOLVE shapefile on the 0.25 degree grid. `ecoregionMetadata.rds` is the attribute table of the same shapefile.

Three small inputs are stored in `data/`:

| File | Content |
|---|---|
| `oni_v5_noaa_cpc_2026-09-15.txt` | Oceanic Niño Index from the NOAA Climate Prediction Center, used to mark El Niño and La Niña years |
| `body_mass_quartiles.csv` | Body-mass quartile, 1 to 4, of each species within its group. The quartiles were estimated with a large language model. They are used only for the body-size panel of Fig. S1 |
| `species_missing_taxon_group_filled.csv` | Taxonomic group of three species that have none in the IUCN tables. The groups were proposed by a large language model |

## Running the analysis

Run every script from the repository root. `R/0_RunAll.r` lists the steps in order.

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 1 | `R/1_Setup.r` | Packages, run settings, directories | |
| 2 | `R/2_Functions.r` | Functions for the exposure calculation | |
| 3 | `R/3b_ClimatePrep.r` | Climate tables and the land grid | `landTemplate.tif`, six climate tables |
| 4 | `R/3a_SpeciesRanges.r` | Species ranges on the grid | one file of cell ids per species |
| 5 | `R/4_Metadata.r` | IUCN category, taxonomy, ecoregions and range size | `spAttributes_v8.qs` |
| 6 | `R/5_CalcExposure.r` | Historical limits and exposed cells for every species and variable | per-species files in `variable_outputs/` |
| 7 | `R/6_MakeDF.r` | Collects the results | `singleMaxExposure_v1.qs`, `AllCellExposureSpXVar.qs`, `singleMaxExposure2025_v1.qs` |
| 8 | `R/8_AppData.r` | Data files for the Shiny app | files in `output/exposureApp_extdata/` |

`singleMaxExposure_v1.qs` has one row per species and year in which the species was exposed, with the variable that had the largest exposed range fraction. The counts in the paper come from this file. `AllCellExposureSpXVar.qs` has one row per species, variable, year and exposed cell.

| Paper item | Script |
|---|---|
| Table 1, Table 2 | `tables/Table1_Table2.r` |
| Table 3 | `tables/Table3_CohortByTaxonIUCN.r` |
| Figure 1, Figure S1 | `figures/Fig1_FigS1_Time.r` |
| Figure 2 | `figures/Fig2_Running.r` |
| Figure 3 | `figures/Fig3_IUCN.r` |
| Figure 4, Figure S5 | `figures/Fig4_FigS5_Maps.r` |
| Figure 5 | `figures/Fig5_Stressor.r` |
| Figure S2 | `figures/FigS2_TimeSensitivity.r` |
| Figure S3 | `figures/FigS3_TimeThreshold.r` |
| Figure S4 | `figures/FigS4_IUCN_Chronic.r` |
| Figure S6 | `figures/FigS6_ProportionalExposure.r` |

`figures/enso_strip.R` holds the functions that mark El Niño and La Niña years. Run `figures/Fig4_FigS5_Maps.r` before `figures/FigS6_ProportionalExposure.r`, because the second uses a raster written by the first.

Not included: the example species figure in Box 2, the screenshots in Figure 6, and the taxon group column that was added to Table 1 for the paper.

### Baseline sensitivity runs

Figure S2 compares the main run with two runs that end the historical baseline earlier. Each is a separate run of steps 6 and 7 with different settings:

```r
END_HISTORICAL_YEAR <- 2015L; RUN_NAME <- "V8_2015"
source("R/1_Setup.r"); source("R/5_CalcExposure.r"); source("R/6_MakeDF.r")
```

Figure S3 compares percentile choices. It needs no extra run, because step 6 stores the limits for several percentiles.

## Checks made on this copy

The scripts here were cleaned from the working code: notes, unused code and machine-specific paths were removed. Each script was then run and its output compared with the output used for the paper.

| Script | Check | Result |
|---|---|---|
| `3b_ClimatePrep.r` | Land grid and climate tables rebuilt from the ERA5 files | Land grid and five tables identical. Warmest 3 months identical to within 0.00003 K, which is the rounding of the stored file |
| `3a_SpeciesRanges.r` | Ranges rebuilt for a sample of 69 species from all four groups | All 69 identical |
| `4_Metadata.r` | Species attributes rebuilt for all species | Group, IUCN category, taxonomy, range size and ecoregions identical for all 32,345 species |
| `5_CalcExposure.r` | Limits and exposed cells recomputed for 68 species and 8 variables | Identical for 7 variables. For the warmest 3 months the limits differ by at most 0.00003 K and no exposed-cell count changes |
| `6_MakeDF.r` | The three data frames rebuilt from the per-species outputs | All identical, including 75.5 million exposed cell records |
| `8_AppData.r` | App data rebuilt | Identical to the data shipped with the app |
| Tables 1 to 3 | Rebuilt | Identical |
| Figures 1 to 4, S1 to S4, S6 | Redrawn | Identical, pixel for pixel, to the existing figure files |
| Figure S5 | Redrawn | Identical to the existing file before one correction. The chronically exposed species are now defined with the 25% of range rule, as in the text, which gives 1,584 species |
| Figure 5 | Redrawn | The underlying tables are identical to the existing tables |

## Notes on reproducibility

- Run R in a UTF-8 locale, which is the default on most systems. `figures/Fig5_Stressor.r` uses an arrow character in its legend labels.

- `R/4_Metadata.r` fills missing taxonomy from the GBIF backbone and the Catalogue of Life through their web services. A rerun can return different values if those services have changed.
- Outputs are written in the `qs2` format. Files in the older `qs` format can be read with the `qread()` function defined in `R/1_Setup.r`.

## Citation

Merow, C., et al. The 2025 Biodiversity Exposure Report. In revision, *BioScience*.

## License

MIT. See [LICENSE.md](LICENSE.md).
