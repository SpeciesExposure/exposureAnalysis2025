# 3b_ClimatePrep.r
# Turns ERA5 monthly means into the climate tables used by the exposure
# calculation. One table per variable, one row per 0.25 degree cell, one column
# per water year 1941-2025. A water year runs October to September: water year
# Y is October of Y-1 to September of Y.
#
# Inputs (INPUTS$era5_dir): temp2m.nc, precip.nc, ERA5_lsm.nc
# Outputs: landTemplate.tif and, in envDir, temp__12, temp__3__max,
# temp__3__min, precip__12, precip__3__max, precip__3__min

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))
if (!exists('envDir')) source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))

# --------------------------------------------
# settings for all variables
processedEnvDir=file.path(EXTERNAL_PATHS$exposure_2025, paste0('processedEnvData_',RUN_VERSION))
.mkdir(processedEnvDir)
rawEnvDir=INPUTS$era5_dir
rglob <- rast(nrows=720, ncols=1440, nlyrs=1, xmin=-180, xmax=180,ymin=-90, ymax=90)
myExt=ext(-180, 180, -61, 86)
edge_grid <- rast(
  xmin = -180, xmax = 180,
  ymin = -61,  ymax = 86,
  resolution = 0.25,
  crs = "EPSG:4326")

#=========================================================
# Land mask: the template raster. Land is where the ERA5 land fraction exceeds 0.5

lm=rast(file.path(rawEnvDir, 'ERA5_lsm.nc'))
lm=lm>.5
lm_rot <- rotate(lm)
lm_crop <- crop(lm_rot, myExt)
mycrs="GEOGCRS[\"WGS 84\",\n    ENSEMBLE[\"World Geodetic System 1984 ensemble\",\n        MEMBER[\"World Geodetic System 1984 (Transit)\"],\n        MEMBER[\"World Geodetic System 1984 (G730)\"],\n        MEMBER[\"World Geodetic System 1984 (G873)\"],\n        MEMBER[\"World Geodetic System 1984 (G1150)\"],\n        MEMBER[\"World Geodetic System 1984 (G1674)\"],\n        MEMBER[\"World Geodetic System 1984 (G1762)\"],\n        MEMBER[\"World Geodetic System 1984 (G2139)\"],\n        MEMBER[\"World Geodetic System 1984 (G2296)\"],\n        ELLIPSOID[\"WGS 84\",6378137,298.257223563,\n            LENGTHUNIT[\"metre\",1]],\n        ENSEMBLEACCURACY[2.0]],\n    PRIMEM[\"Greenwich\",0,\n        ANGLEUNIT[\"degree\",0.0174532925199433]],\n    CS[ellipsoidal,2],\n        AXIS[\"geodetic latitude (Lat)\",north,\n            ORDER[1],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n        AXIS[\"geodetic longitude (Lon)\",east,\n            ORDER[2],\n            ANGLEUNIT[\"degree\",0.0174532925199433]],\n    USAGE[\n        SCOPE[\"Horizontal component of 3D system.\"],\n        AREA[\"World.\"],\n        BBOX[-90,-180,90,180]],\n    ID[\"EPSG\",4326]]"
lm_sh <- project(lm_crop, mycrs, method = "near")
lm_sh <- resample(lm_sh, edge_grid, method = "near")
lm_sh[lm_sh == 0] <- NA
lm_sh[lm_sh == 1] <- 0
writeRaster(lm_sh,file=paste0(intDir,'/landTemplate.tif'),overwrite=T)
template=rast(paste0(intDir,'/landTemplate.tif'))

# one row for each cell of the template

t23=tibble(cell=1:ncell(template), isLand=values(template)[,1])	

# ---------------------------------------------------------------
# ---------------------------------------------------------------
# ---------------------------------------------------------------
#--------------------------------------------
# annual temp -------------------------------
# Read the full monthly series and summarize by water year
a=rast(paste0(rawEnvDir,'/temp2m.nc'))
ts_raw <- sub("t2m_valid_time=", "", names(a))
# convert to POSIXct
ts_date <- as.POSIXct(as.numeric(ts_raw), origin = "1970-01-01", tz = "UTC")
# format as YYYY-MM
new_names <- format(ts_date, "%Y-%m")
# assign back
names(a) <- new_names

# Extract year and month
dates <- names(a)
yr  <- as.integer(str_sub(dates, 1, 4))
mo  <- as.integer(str_sub(dates, 6, 7))
# Compute water year: WY starts in October
water_year <- ifelse(mo >= 10, yr + 1, yr)
# Assign names like "WY1941"
wy_groups <- paste0("WY", water_year)
# Now compute mean for each group
a_wy_mean <- tapp(a, wy_groups, fun = mean, na.rm = TRUE)
# The output layers will have names like "WY1941", "WY1942", ...

# the source grid is 0-360 and has slightly different dimensions
v1.1 <- a_wy_mean |>
  rotate () |>  crop(myExt)  %>% resample(edge_grid, method = "near") 
writeRaster(v1.1,file=paste0(processedEnvDir,'/temp__12.tif'))

# get data frame
d24=values(v1.1)
v1.2=t23[,1] %>% bind_cols(d24)
v1.2=v1.2 %>% dplyr::select(-WY1940,-WY2026)
names(v1.2) = c('cell',1941:2025)
qs_save(v1.2,file=paste0(envDir,'/temp__12.qs'))

#----------------------------------------------------------
# hot season ----------------------------------------------
dates <- as.Date(paste0(names(a), "-01"))
# Water year: Oct → Sep
water_year <- ifelse(format(dates, "%m") >= "10",
                     as.integer(format(dates, "%Y")) + 1,
                     as.integer(format(dates, "%Y")))

# Add water year to layers
df <- tibble(layer = seq_along(dates),
             date = dates,
             wy = water_year)

cellwise_3mo_extreme <- function(r_monthly, wy_index, mode = c('max', 'min')) {
  mode <- match.arg(mode)
  wy_valid <- wy_index[vapply(wy_index, nrow, integer(1)) >= 3]

  out_list <- lapply(wy_valid, function(df_wy) {
    layers <- df_wy$layer
    r <- r_monthly[[layers]]

    triplet_means <- lapply(seq_len(nlyr(r) - 2), function(i) {
      (r[[i]] + r[[i + 1]] + r[[i + 2]]) / 3
    })

    triplet_stack <- rast(triplet_means)
    if (mode == 'max') {
      app(triplet_stack, max, na.rm = TRUE)
    } else {
      app(triplet_stack, min, na.rm = TRUE)
    }
  })

  out <- rast(out_list)
  names(out) <- names(wy_valid)
  names(out) <- sub('^X(?=[0-9]{4}$)', '', names(out), perl = TRUE)
  out
}

# For each water year: find warmest 3-month consecutive window
wy_list <- split(df, df$wy)
warm3 <- cellwise_3mo_extreme(a, wy_list, mode = 'max')

# the source grid is 0-360 and has slightly different dimensions
v1.1a <- warm3 |>
  rotate () |> 
  crop(myExt)  %>% resample(edge_grid, method = "near") 
writeRaster(v1.1a,file=paste0(processedEnvDir,'/temp__3__max.tif'), overwrite=T)

# get data frame
d24=values(v1.1a)
v1.2=t23[,1] %>% bind_cols(d24)
v1.2=v1.2 %>% dplyr::select(-any_of('1940'))
qs_save(v1.2,file=paste0(envDir,'/temp__3__max.qs'))

#======================================================
# cold season
# For each water year: find the coldest 3-month consecutive window
warm3 <- cellwise_3mo_extreme(a, wy_list, mode = 'min')

# the source grid is 0-360 and has slightly different dimensions
v1.1 <- warm3 |>
  rotate () |>
  crop(myExt)  %>% resample(edge_grid, method = "near") 
writeRaster(v1.1,file=paste0(processedEnvDir,'/temp__3__min.tif'), overwrite=T)

# get data frame
d24=values(v1.1)
v1.2=t23[,1] %>% bind_cols(d24)
v1.2=v1.2 %>% dplyr::select(-any_of('1940'))
qs_save(v1.2,file=paste0(envDir,'/temp__3__min.qs'))

#==========================================================
#==========================================================
#==========================================================
#==========================================================
#==========================================================
#--------------------------------------------
# annual precip  -------------------------------
# Read the full monthly series and summarize by water year

a=rast(paste0(rawEnvDir,'/precip.nc'))
ts_raw <- sub("tp_valid_time=", "", names(a))
# convert to POSIXct
ts_date <- as.POSIXct(as.numeric(ts_raw), origin = "1970-01-01", tz = "UTC")
# format as YYYY-MM
new_names <- format(ts_date, "%Y-%m")
# assign back
names(a) <- new_names

# Extract year and month
dates <- names(a)
yr  <- as.integer(str_sub(dates, 1, 4))
mo  <- as.integer(str_sub(dates, 6, 7))
# Compute water year: WY starts in October
water_year <- ifelse(mo >= 10, yr + 1, yr)
# Assign names like "WY1941"
wy_groups <- paste0("WY", water_year)
# sum for each water year
a_wy_mean <- tapp(a, wy_groups, fun = sum, na.rm = TRUE)
# The output layers will have names like "WY1941", "WY1942", ...

names(a_wy_mean)=sub('WY','',names(a_wy_mean))

# resample
v1.1 <- a_wy_mean |>
  rotate () |>
	crop(myExt)  %>% resample(edge_grid, method = "near") 
writeRaster(v1.1,file=paste0(processedEnvDir,'/precip__12.tif'),overwrite=T )

# get data frame
d24=values(v1.1)
v1.2=t23[,1] %>% bind_cols(d24)
v1.2=v1.2 %>% dplyr::select(-`1940`,-`2026`)
qs_save(v1.2,file=paste0(envDir,'/precip__12.qs'))

#--------------------------------------------
# wet season ----------------------------------------------
dates <- as.Date(paste0(names(a), "-01"))
# Water year: Oct → Sep
water_year <- ifelse(format(dates, "%m") >= "10",
                     as.integer(format(dates, "%Y")) + 1,
                     as.integer(format(dates, "%Y")))
# Add water year to layers
df <- tibble(layer = seq_along(dates),
             date = dates,
             wy = water_year)

# For each water year: find the wettest 3-month consecutive window
wy_list <- split(df, df$wy)
warm3 <- cellwise_3mo_extreme(a, wy_list, mode = 'max')

# resample
v1.1 <- warm3 |>
  rotate () |>
  crop(myExt)  %>% resample(edge_grid, method = "near") 
writeRaster(v1.1,file=paste0(processedEnvDir,'/precip__3__max.tif'), overwrite=T)

# get data frame
d24=values(v1.1)
v1.2=t23[,1] %>% bind_cols(d24)
v1.2=v1.2 %>% dplyr::select(-any_of('1940'))
qs_save(v1.2,file=paste0(envDir,'/precip__3__max.qs'))

#=======================================
# dry season
warm3 <- cellwise_3mo_extreme(a, wy_list, mode = 'min')

# resample
v1.1 <- warm3 |> rotate () |> crop(myExt)  %>% resample(edge_grid, method = "near") 
writeRaster(v1.1,file=paste0(processedEnvDir,'/precip__3__min.tif'), overwrite=T)

# get data frame
d24=values(v1.1)
v1.2=t23[,1] %>% bind_cols(d24)
v1.2=v1.2 %>% dplyr::select(-any_of('1940'))
qs_save(v1.2,file=paste0(envDir,'/precip__3__min.qs'))

