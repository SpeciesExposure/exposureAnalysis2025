# 5_CalcExposure.r
# For every species and climate variable: (1) the historical limit, (2) the
# cells and years beyond it, (3) and (4) per-species summaries.
# Upper limits are computed for temp__12, temp__3__max, precip__12 and
# precip__3__max. Lower limits are computed for temp__12, temp__3__min,
# precip__12 and precip__3__min.
#
# Set FORCE_OVERWRITE_EXPOSURE <- TRUE to recompute outputs that already exist.
if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))

if (!exists('RUN_VERSION') || !exists('envDir') || !exists('dataDir') || !exists('metaDir')) {
	source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))
}

all.sp.f=list.files(EXTERNAL_PATHS$ranges_expert,full.names=T)

overwrite <- isTRUE(get0('FORCE_OVERWRITE_EXPOSURE', envir = .GlobalEnv, inherits = TRUE, ifnotfound = FALSE))
message('5_CalcExposure overwrite mode: ', overwrite)
vars1=list.files(envDir) %>% file_path_sans_ext 
vars1=vars1[ !grepl("min", vars1, ignore.case = TRUE) ]

quants=c(.95,.99,1)
quantsYearCell=c(.995,.999,1)
types=c('up')
vars=paste0(vars1,'_',types)
END_HISTORICAL_YEAR=as.integer(END_HISTORICAL_YEAR)

for(v in seq_along(vars)){
	print('-------------------')
	print(v)
	# 1. Calculate thresholds
	v.dir=paste0(dataDir,'/variable_outputs/',vars[v]); .mkdir(v.dir)
	outDir=paste0(v.dir,'/',vars[v],'_Thresh'); .mkdir(outDir)
	# species with no land cells are copied here
	brokeDir=paste0(v.dir,'/',vars[v],'_broke'); .mkdir(brokeDir)
	clim=qread(paste0(envDir,'/',vars1[v],'.qs'))
	clim_lookup_thresh=.make_clim_lookup(clim=clim,endHistoricalPeriod=END_HISTORICAL_YEAR)
	clim_lookup_full=.make_clim_lookup(clim=clim)
	done=mclapply(all.sp.f, function(x) {
		out.f=paste0(outDir,'/',basename(x))
		if(!file.exists(out.f) | overwrite ) 	almostEveryQuantileEver2(x,clim,outDir,brokeDir,quants,varName=vars[v],quantsYearCell,whichQuantileMethods=7,
			endHistoricalPeriod=END_HISTORICAL_YEAR,clim_lookup=clim_lookup_thresh)
	},mc.cores=mc.cores)
	# 2. Calculate exposure and marginal cells 
		outDirMarg=paste0(v.dir,'/',vars[v],'_marg'); .mkdir(outDirMarg)
	marg.fun(
		all.sp.f=all.sp.f,
		clim=clim,
		varName=vars[v],
		type=types[1],
		outDirMarg=outDirMarg,
		outDir=outDir,
		mc.cores=mc.cores,
		overwrite=overwrite,
		clim_lookup=clim_lookup_full
	)
	# 3. Count exposed cells per year
	outDirStats=paste0(v.dir,'/',vars[v],'_stats'); .mkdir(outDirStats)
	outDirSp=paste0(v.dir,'/',vars[v],'_spStats'); .mkdir(outDirSp)
	all.sp.f2=list.files(outDirMarg,full.names=T)
	statsFun2024(
		all.sp.f2=all.sp.f2,
		outDirStats=outDirStats,
		outDirMarg=outDirMarg,
		outDirSp=outDirSp,
		mc.cores=mc.cores,
		yearMax=2025
	) 
	# 4. Same counts for the subset of thresholds kept for mapping
	qSubsets=tibble(quantileIndex=as.integer(c(1,5,6,8,9)))
	outDirSpXCell=paste0(v.dir,'/',vars[v],'_spXCell'); .mkdir(outDirSpXCell)
	otherStatsFun(
		all.sp.f2=all.sp.f2,
		qSubset=qSubsets,
		outDirSpXCell=outDirSpXCell,
		outDirMarg=outDirMarg,
		outDirSp=outDirSp,
		mc.cores=mc.cores,
		yearMax=2025
	)
}

# lower limits
quants=1-c(.95,.99,1)
quantsYearCell=1-c(.995,.999,1)
types=c('lo')
vars1=list.files(envDir) %>% file_path_sans_ext 
vars2=vars1[-grep('max',vars1)]
vars=paste0(vars2,'_',types)

# lower limits are not computed for the warmest or wettest 3 months
for(v in seq_along(vars)){
	print(v)
	# 1. Calculate thresholds
	v.dir=paste0(dataDir,'/variable_outputs/',vars[v]); .mkdir(v.dir)
	outDir=paste0(v.dir,'/',vars[v],'_Thresh'); .mkdir(outDir)
	# species with no land cells are copied here
	brokeDir=paste0(v.dir,'/',vars[v],'_broke'); .mkdir(brokeDir)
	clim=qread(paste0(envDir,'/',vars2[v],'.qs'))
	clim_lookup_thresh=.make_clim_lookup(clim=clim,endHistoricalPeriod=END_HISTORICAL_YEAR)
	clim_lookup_full=.make_clim_lookup(clim=clim)
	print('thresholds')
	done=mclapply(all.sp.f, function(x) {
		out.f=paste0(outDir,'/',basename(x))
		if(!file.exists(out.f) | overwrite ) 	almostEveryQuantileEver2(x,clim,outDir,brokeDir,quants,varName=vars[v],quantsYearCell,whichQuantileMethods=7,endHistoricalPeriod=END_HISTORICAL_YEAR,clim_lookup=clim_lookup_thresh)
	},mc.cores=mc.cores)
	# 2. Calculate exposure and marginal cells 
	print('margins')
	outDirMarg=paste0(v.dir,'/',vars[v],'_marg'); .mkdir(outDirMarg)
	marg.fun(
		all.sp.f=all.sp.f,
		clim=clim,
		varName=vars[v],
		type=types[1],
		outDirMarg=outDirMarg,
		outDir=outDir,
		mc.cores=mc.cores,
		overwrite=overwrite,
		verbose=T,
		clim_lookup=clim_lookup_full
	)
	# 3. Count exposed cells per year
 	print('stats1')
	outDirStats=paste0(v.dir,'/',vars[v],'_stats'); .mkdir(outDirStats)
	outDirSp=paste0(v.dir,'/',vars[v],'_spStats'); .mkdir(outDirSp)
	all.sp.f2=list.files(outDirMarg,full.names=T)
	statsFun2024(
		all.sp.f2=all.sp.f2,
		outDirStats=outDirStats,
		outDirMarg=outDirMarg,
		outDirSp=outDirSp,
		mc.cores=mc.cores,
		yearMax=2025
	)
	# 4. Same counts for the subset of thresholds kept for mapping
	print('stats2')
	qSubsets=tibble(quantileIndex=as.integer(c(1,5,6,8,9)))
	outDirSpXCell=paste0(v.dir,'/',vars[v],'_spXCell'); .mkdir(outDirSpXCell)
	otherStatsFun(
		all.sp.f2=all.sp.f2,
		qSubset=qSubsets,
		outDirSpXCell=outDirSpXCell,
		outDirMarg=outDirMarg,
		outDirSp=outDirSp,
		mc.cores=mc.cores,
		yearMax=2025
	)
}

