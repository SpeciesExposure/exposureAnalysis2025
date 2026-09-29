# 0_RunAll.r
# Runs the whole analysis in order. Run from the repository root.
#
# The exposure calculation, 5_CalcExposure.r, is the long step. It writes one
# file per species and skips files that already exist, so it can be stopped and
# restarted.

source(file.path('config', 'paths.R'))

# ---- main run: historical baseline 1941-2022 ----
END_HISTORICAL_YEAR <- 2022L
RUN_NAME <- 'V8'

source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))
source(file.path(PROJECT_PATHS$src_r, '3b_ClimatePrep.r'))    # climate tables and land grid
source(file.path(PROJECT_PATHS$src_r, '3a_SpeciesRanges.r'))  # species ranges on the grid
source(file.path(PROJECT_PATHS$src_r, '4_Metadata.r'))        # species attributes
source(file.path(PROJECT_PATHS$src_r, '5_CalcExposure.r'))    # historical limits and exposed cells
source(file.path(PROJECT_PATHS$src_r, '6_MakeDF.r'))          # data frames for figures and tables

# ---- baseline sensitivity runs for Fig. S2 ----
for (end_year in c(2015L, 2009L)) {
	END_HISTORICAL_YEAR <- end_year
	RUN_NAME <- paste0('V8_', end_year)
	source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))
	source(file.path(PROJECT_PATHS$src_r, '5_CalcExposure.r'))
	source(file.path(PROJECT_PATHS$src_r, '6_MakeDF.r'))
}

# ---- tables and figures, from the main run ----
END_HISTORICAL_YEAR <- 2022L
RUN_NAME <- 'V8'
source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))

source(file.path(PROJECT_PATHS$tables_src, 'Table1_Table2.r'))
source(file.path(PROJECT_PATHS$tables_src, 'Table3_CohortByTaxonIUCN.r'))

source(file.path(PROJECT_PATHS$figures_src, 'Fig1_FigS1_Time.r'))
source(file.path(PROJECT_PATHS$figures_src, 'Fig2_Running.r'))
source(file.path(PROJECT_PATHS$figures_src, 'Fig3_IUCN.r'))
source(file.path(PROJECT_PATHS$figures_src, 'Fig4_FigS5_Maps.r'))
source(file.path(PROJECT_PATHS$figures_src, 'Fig5_Stressor.r'))
source(file.path(PROJECT_PATHS$figures_src, 'FigS2_TimeSensitivity.r'))
source(file.path(PROJECT_PATHS$figures_src, 'FigS3_TimeThreshold.r'))
source(file.path(PROJECT_PATHS$figures_src, 'FigS4_IUCN_Chronic.r'))
source(file.path(PROJECT_PATHS$figures_src, 'FigS6_ProportionalExposure.r'))

# ---- data for the Shiny app ----
source(file.path(PROJECT_PATHS$src_r, '8_AppData.r'))
