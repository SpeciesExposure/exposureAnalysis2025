# 3a_SpeciesRanges.r
# Puts every species range on the 0.25 degree land grid. The output is one file
# per species holding the ids of the grid cells in its range.
#
# Mammals and birds: Area of Habitat maps (Lumbierres et al. 2022).
#   For migratory birds the breeding (B) and resident (R) maps are combined.
#   Non-breeding (N) maps are not used.
# Amphibians and reptiles: IUCN Red List range polygons.
#
# Inputs:  INPUTS$aoh_dir, INPUTS$iucn_maps_dir, landTemplate.tif
# Outputs: EXTERNAL_PATHS$ranges_expert/<Genus_species>.qs
#
# Run 3b_ClimatePrep.r first, because it writes landTemplate.tif.

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))
if (!exists('template')) source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))

tableOutDir=EXTERNAL_PATHS$ranges_expert

#=============================================================
# Mammals and non-migratory birds
spIn=c(list.files(file.path(INPUTS$aoh_dir,'AoH_Mammals'),recursive=T,full.names=T), list.files(file.path(INPUTS$aoh_dir,'AoH_Birds_Nonmigratory'),recursive=T,full.names=T))

outinfo = parallel::mclapply(seq_along(spIn),function(x){
	out_qsfile=paste0(tableOutDir,'/', basename(file_path_sans_ext(spIn[x])),'.qs')
	print(x)
	try({
	f=rast(spIn[x]) %>% aggregate(fact=16,fun='max',na.rm=T)
	pts <- crds(f, na.rm = TRUE)
	# assign the habitat cells to the template grid
	f <- rasterize(pts, template, field = 1)
	f=mask(f,template)
	notNACells=cells(f)
	qs_save(notNACells,file = out_qsfile)
	})
},mc.cores=mc.cores)

#=============================================================
# Migratory birds: combine breeding (B) and resident (R) maps
spInB=list.files(file.path(INPUTS$aoh_dir,'AoH_Birds_Migratory'),pattern='_B.tif',recursive=T,full.names=T)
spInR=list.files(file.path(INPUTS$aoh_dir,'AoH_Birds_Migratory'),pattern='_R.tif',recursive=T,full.names=T)
sp=unique(c(sub('_B.tif','',spInB),sub('_R.tif','',spInR)))
outinfo = parallel::mclapply(seq_along(sp),function(x){
	print(x)
	out_qsfile=paste0(tableOutDir,'/', basename(file_path_sans_ext(sp[x])),'.qs')
	try({
	xx=list(spInB[grep(paste0(sp[x],'_'),spInB)], spInR[grep(paste0(sp[x],'_'),spInR)])
	xxx=sapply(Filter(length,xx),function(y) rast(y) %>% aggregate(fact=16,fun='max',na.rm=T))
	if(length(xxx)==2) { f <- try(do.call(merge, xxx))} else {f=xxx[[1]]}
	if(class(f)=='try-error') {print(sp[x]); return()}
	pts <-f %>%  crds(na.rm = TRUE)
	f <- rasterize(pts, template, field = 1)
	f=mask(f,template)
	notNACells=cells(f)
	qs_save(notNACells,file = out_qsfile)
	})
},mc.cores=mc.cores)

#=============================================================
# Amphibians and reptiles: IUCN range polygons
pp1=st_read(file.path(INPUTS$iucn_maps_dir,'AMPHIBIANS','AMPHIBIANS_PART1.shp')) %>% bind_rows(st_read(file.path(INPUTS$iucn_maps_dir,'AMPHIBIANS','AMPHIBIANS_PART2.shp')) ) %>% mutate(sci_name=gsub('[[:space:]]','_',sci_name))
pp2=st_read(file.path(INPUTS$iucn_maps_dir,'REPTILES','REPTILES_PART1.shp')) %>% bind_rows(st_read(file.path(INPUTS$iucn_maps_dir,'REPTILES','REPTILES_PART2.shp')) ) %>% mutate(sci_name=gsub('[[:space:]]','_',sci_name))
# Polygons labelled Extinct are dropped. Polygons labelled presence uncertain are kept.
spIn=pp1 %>% bind_rows(pp2) %>% filter(!grepl("Extinct", legend, ignore.case = TRUE) ) %>% filter(terrestria=='true')
sp=unique(spIn$sci_name)

sf_use_s2(FALSE)
outinfo = parallel::mclapply(1:length(sp),function(x){
	print(x)
	out_qsfile=paste0(tableOutDir,'/', sp[x],'.qs')
	xxx=spIn[grep(sp[x],spIn$sci_name),]
	if(nrow(xxx)>1) { this.poly <- try(st_union(xxx))} else {this.poly=xxx}
	this.poly <- st_transform(this.poly, crs = crs(template))
	# a cell is in the range if the polygon touches it
	this.rast=this.poly %>% vect %>% terra::rasterize(template,touches=T)
	this.rast=mask(this.rast,template)
	(notNACellsm=cells(this.rast))
	qs_save(notNACellsm,file = out_qsfile)
},mc.cores=mc.cores)
sf_use_s2(TRUE)
