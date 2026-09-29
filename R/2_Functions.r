# 2_Functions.r
# Functions used by the exposure pipeline.

#==========================================
# File helpers
#==========================================

.mkdir=function(x){if(!file.exists(x)) dir.create(x)}

.be =function(f) {basename(file_path_sans_ext(f))}

#==========================================
# Climate extraction helpers
#==========================================

.make_clim_lookup=function(clim, endHistoricalPeriod=NULL){
	clim_cells=as.integer(clim$cell)
	year_cols=setdiff(names(clim),'cell')
	if(!is.null(endHistoricalPeriod)){
		year_num=suppressWarnings(as.integer(year_cols))
		year_cols=year_cols[!is.na(year_num) & year_num<=endHistoricalPeriod]
	}
	clim_mat=as.matrix(clim[,year_cols,drop=FALSE])
	list(clim_cells=clim_cells,year_cols=year_cols,clim_mat=clim_mat)
}

.subset_clim_to_long=function(clim, occ_cells, varName, endHistoricalPeriod=NULL, clim_lookup=NULL){
	if(is.null(clim_lookup)){
		clim_lookup=.make_clim_lookup(clim=clim,endHistoricalPeriod=endHistoricalPeriod)
	}
	clim_cells=clim_lookup$clim_cells
	year_cols=clim_lookup$year_cols
	n_occ=length(occ_cells)
	n_year=length(year_cols)
	if(n_occ==0 || n_year==0){
		out=tibble(cell=numeric(0),year=character(0),value=numeric(0))
		names(out)[3]=varName
		return(out)
	}
	clim_mat=clim_lookup$clim_mat
	idx=match(as.integer(occ_cells),clim_cells)
	vals=matrix(NA_real_,nrow=n_occ,ncol=n_year)
	good=!is.na(idx)
	if(any(good)) vals[good,]=clim_mat[idx[good],,drop=FALSE]
	out=tibble(cell=rep(as.numeric(occ_cells),each=n_year),year=rep(year_cols,times=n_occ),value=as.vector(t(vals)))
	names(out)[3]=varName
	out
}

#==========================================
# Map palette
#==========================================

cm.cols1=function(x,bias=1) { colorRampPalette(c('grey90','steelblue4','steelblue1','gold','red1','red4'),bias=bias)(x)
}

#==========================================
# Threshold and exposure pipeline helpers
#==========================================

almostEveryQuantileEver2=function(sp.f,clim,outDir,brokeDir,quants,varName,quantsYearCell,endHistoricalPeriod=2022,focalYear=2023,whichQuantileMethods=1:9,verbose=F,clim_lookup=NULL){ 
	if(verbose) print(basename(sp.f))
	sp1=qread(sp.f)
	if(length(sp1)==0) {
		file.copy(sp.f,paste0(brokeDir,'/',basename(sp.f)));
		warning(basename(sp.f))
		return()}
	occ.cells=as.integer(sp1)
	sp.clim=.subset_clim_to_long(clim=clim,occ_cells=occ.cells,varName=varName,endHistoricalPeriod=endHistoricalPeriod,clim_lookup=clim_lookup)
	# 1. mean climate of each cell
	cell.means=sp.clim %>% group_by(cell) %>% dplyr::summarize(mean.clim=mean(!!sym(varName)))
	
	# 2. quantile over years within each cell, then quantile of those values over cells
		mythresh1=lapply(whichQuantileMethods,function(ty){sp.clim %>%  group_by(cell) %>%  reframe(year.q.val = quantile(!!sym(varName), quants,type=ty,na.rm=T), year.q = quants) %>% mutate(quantileType=ty,isRounded=F)}) %>% bind_rows
		# quantile over cells
		mythresh2=lapply(whichQuantileMethods,function(ty){ 
			tmp=mythresh1%>% ungroup %>% group_by(year.q,isRounded) 
				# keep quantile types separate
			tmp1=tmp %>% filter(quantileType==ty) %>% reframe(thresh.val = quantile(year.q.val, quants,type=ty,na.rm=T), cell.q = quants) %>% mutate(quantileType=ty)
			}) %>% bind_rows %>% mutate(year.cell.q=NA)
		
	# 3. quantile over all year x cell values pooled
		mythresh3=lapply(whichQuantileMethods,function(ty){sp.clim %>%  reframe(thresh.val = quantile(!!sym(varName), quantsYearCell,type=ty,na.rm=T), year.cell.q = quantsYearCell) %>% mutate(quantileType=ty,isRounded=F)}) %>% bind_rows
			mythresh.b = mythresh3 %>% mutate(year.q=NA,cell.q=NA)

	#  4. organize all the outputs and store
	allMyThresh=mythresh2 %>% bind_rows(mythresh.b) 
	allMyThresh =	allMyThresh %>% mutate(quantileIndex=1:nrow(	allMyThresh))
	
	qsave(allMyThresh,file=paste0(outDir,'/',basename(sp.f)))
}

#----------------------------
marg.fun=function(all.sp.f,clim,varName,type='up',outDirMarg,outDir,sensVal=0,mc.cores=11,overwrite=FALSE,verbose=F,clim_lookup=NULL){
	marg=mclapply(all.sp.f, function(x) {
		tryCatch({
		if(verbose) print(basename(x))
		out.f=paste0(outDirMarg,'/',basename(x))
		tr.file=paste0(outDir,'/',basename(x))
		if(!file.exists(tr.file)) return()
		if(!overwrite){
			if(file.exists(out.f)) 	return()
		}
		# read in everything
		occ.cells=as.integer(qread(x))

		if(length(occ.cells)==0) {print(paste0('no cells: ',x)); return()}
		# all years are used here, including those after the historical period
		sp.clim=.subset_clim_to_long(clim=clim,occ_cells=occ.cells,varName=varName,clim_lookup=clim_lookup)
		cell.means=sp.clim %>% group_by(cell) %>% dplyr::summarize(mean.clim=mean(!!sym(varName)))
		allMyThresh=qread(tr.file)
		# optional offset to the thresholds, for sensitivity analysis
		allMyThresh$thresh.val=	allMyThresh$thresh.val + sensVal
	
			# marginal cells: cells whose mean climate is already beyond the threshold
			# split by threshold type: year then cell, or year x cell pooled
		t.y.c=allMyThresh %>% filter(!isRounded,is.na(year.cell.q))	
		t.yc=allMyThresh %>% filter(!isRounded,!is.na(year.cell.q))	
		
		# thresholds defined by year then cell
		# direction depends on whether this is an upper or lower limit
		if(type=='up'){
			t.y.c1a=t.y.c %>% mutate(marginalCells=lapply(1:nrow(t.y.c),function(y){
				cell.means %>% filter(mean.clim>t.y.c$thresh.val[y]) %>% pull(cell) })) %>% mutate(nMarginalCells=sapply(marginalCells,length))  %>% mutate(rangeSize=nrow(cell.means)) %>% mutate(exposedYearCells=lapply(1:nrow(t.y.c),function(y){
						sp.clim %>% filter(!!sym(varName)>t.y.c$thresh.val[y]) %>% mutate(mag=!!sym(varName)-t.y.c$thresh.val[y])%>% dplyr::select(cell,year,mag) %>% mutate(cell=as.integer(cell),year=as.integer(year))})) 
		}
		if(type=='lo'){
			t.y.c1a=t.y.c %>% mutate(marginalCells=lapply(1:nrow(t.y.c),function(y){
					cell.means %>% filter(mean.clim<t.y.c$thresh.val[y]) %>% pull(cell) }))  %>% mutate(nMarginalCells=sapply(marginalCells,length))  %>% mutate(rangeSize=nrow(cell.means)) %>% 
			 		mutate(exposedYearCells=lapply(1:nrow(t.y.c),function(y){
						sp.clim %>% filter(!!sym(varName)<t.y.c$thresh.val[y]) %>%  
					mutate(mag=!!sym(varName)-t.y.c$thresh.val[y])%>% dplyr::select(cell,year,mag) %>% mutate(cell=as.integer(cell),year=as.integer(year))})) 
		}	
			
		t.y.c1 = t.y.c1a %>% mutate(expEvents=sapply(exposedYearCells,nrow)) %>% mutate(expByYear=lapply(exposedYearCells,function(z) { z %>% count(year) }))
		qsave(t.y.c1,file=out.f)
			list(ok=TRUE,file=basename(x),path=x)
			}, error=function(e){
				list(
					ok=FALSE,
					file=basename(x),
					path=x,
					message=conditionMessage(e)
				)
			})
	},mc.cores=mc.cores)

		err.idx=which(vapply(marg,function(z){is.list(z) && !isTRUE(z$ok)},logical(1)))
		if(length(err.idx)>0){
			err.lines=vapply(marg[err.idx],function(z){paste0(z$file,' :: ',z$message)},character(1))
			stop(paste0('marg.fun failed for ',length(err.idx),' file(s):\n',paste(err.lines,collapse='\n')),call.=FALSE)
		}
		invisible(marg)
}

#==========================================
# Exposure summary/statistics helpers
#==========================================

statsFun2024=function(all.sp.f2,outDirStats,outDirMarg,outDirSp,mc.cores=11,yearMin=1990,yearMax=2024,verbose=F){
	stats=mclapply(all.sp.f2, function(x1) {
		if(verbose) print(basename(x1))
		out.f=paste0(outDirStats,'/',basename(x1))
		mar.file=paste0(outDirMarg,'/',basename(x1))
		if(file.exists(out.f) | !file.exists(mar.file)) 	return()
		# read in everything
		mar=qread(mar.file) 
			# note: rows are in a different order than in the threshold file
		# to store stats
		toKeep=mar %>% dplyr::select(-exposedYearCells,-marginalCells,-expByYear)
		#------------------------------
		# unnest count of exposure by year. combine all the years before 1990
			# can do years and yearcells at the same time
		tmp1=mar  %>%  dplyr::select(quantileIndex,expByYear) %>% unnest(expByYear) 
		tmp1.pre90=tmp1 %>% filter(year<yearMin) %>% group_by(quantileIndex) %>% dplyr::summarize(expPre1990=sum(n),.groups='drop')
		tmp1=tmp1 %>% filter(year>=yearMin) %>%  mutate(year=paste0('X',year))
		tmp2.qi=mar %>% dplyr::select(quantileIndex)
		tmp2.template=expand_grid(quantileIndex=tmp2.qi$quantileIndex,year=paste0('X',yearMin:yearMax))
		tmp2=tmp2.template %>% left_join(tmp1 %>% dplyr::select(quantileIndex,year,n),by=c('quantileIndex','year')) %>% mutate(n=replace_na(n,0)) %>% pivot_wider(names_from='year',values_from='n')
		tmp2=tmp2.qi %>% left_join(tmp2,by='quantileIndex') %>% left_join(tmp1.pre90,by='quantileIndex') %>% mutate(expPre1990=replace_na(expPre1990,0))
		toKeep = toKeep %>% left_join(tmp2,by='quantileIndex')
		#------------------------------
		# find marginal cells that were counted as exposed for year.q and cell.q. ONLY APPLIES TO year + cell
			# exposed years and cells
		tmp3a=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex,exposedYearCells) %>% unnest(exposedYearCells) 
			# marginal cells
		tmp3a.1=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex,marginalCells) %>% unnest(marginalCells) %>% rename(cell=marginalCells)
			#find when marginal contribute to exposure in each year
		new.vars=paste0('mar',yearMin:yearMax)
		tmp5=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex) %>% mutate(!!!setNames(rep(0, length(new.vars)), new.vars))
		if(nrow(tmp3a.1)>0){
			inBoth2=tmp3a %>% inner_join(tmp3a.1,by=c('quantileIndex','cell'))
			tmp4a=inBoth2 %>% filter(year>=yearMin) %>%  mutate(year=paste0('mar',year)) %>% group_by(quantileIndex,year) %>% count
			tmp5.qi=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex)
			tmp5.template=expand_grid(quantileIndex=tmp5.qi$quantileIndex,year=paste0('mar',yearMin:yearMax))
			tmp5=tmp5.template %>% left_join(tmp4a,by=c('quantileIndex','year')) %>% mutate(n=as.numeric(replace_na(n,0))) %>% pivot_wider(names_from='year',values_from='n')
		}
		toKeep.yc=toKeep %>% filter(!is.na(year.cell.q))
		toKeep.y.c = toKeep %>% filter(is.na(year.cell.q)) %>% left_join(tmp5,by='quantileIndex')
		expNotMarg = (toKeep.y.c %>% dplyr::select(X1990:X2024)) - (toKeep.y.c %>% dplyr::select(mar1990:mar2024))
		names(expNotMarg)=sub('X','XnM',names(expNotMarg))
		toKeep.y.c = 	toKeep.y.c %>% bind_cols(expNotMarg)
			# organize and store
		tk=toKeep.y.c %>% bind_rows(toKeep.yc)
		num_cols=names(tk)[grepl('^(X|mar|XnM)|^expPre1990$',names(tk))]
		if(length(num_cols)>0) tk[num_cols]=lapply(tk[num_cols],as.numeric)
		qsave(tk,file=out.f)
	
		# =======================
		# species-level summary statistics
		spStats=tk %>% dplyr::select(year.q:quantileIndex,rangeSize,X2015,X2016,X2019,X2020,X2023,X2024,XnM2015,XnM2016,XnM2019,XnM2020,XnM2023,) %>% mutate(pExp2015=round(X2015/rangeSize,3),pExp2016=round(X2016/rangeSize,3),pExp2019=round(X2019/rangeSize,3),pExp2020=round(X2020/rangeSize,3),pExp2023=round(X2023/rangeSize,3),pExp2024=round(X2024/rangeSize,3))
		qsave(spStats,file=paste0(outDirSp,'/',basename(x1)))

	},mc.cores=mc.cores)

}

#-------------------------------------
otherStatsFun=function(all.sp.f2,qSubset,outDirSpXCell,outDirMarg,outDirSp,mc.cores,yearMin=1990,yearMax=2024,verbose=F){ 
	stats=mclapply(all.sp.f2, function(x1) {
		if(verbose) print(basename(x1))
		out.f=paste0(outDirSpXCell,'/',basename(x1))
		mar.file=paste0(outDirMarg,'/',basename(x1))
		if(file.exists(out.f) | !file.exists(mar.file)) 	return()
		# read in everything
		mar=qread(mar.file) 
		# subset
		mar2=mar %>% right_join(qSubset) %>% filter(quantileType==7)
		# to store stats
		toKeep=mar2 %>% dplyr::select(-exposedYearCells,-marginalCells,-expByYear)
		#------------------------------
		# unnest count of exposure by year. combine all the years before 1990
			# can do years and yearcells at the same time
		tmp1=mar  %>%  dplyr::select(quantileIndex,expByYear) %>% unnest(expByYear) 
		tmp1.pre90=tmp1 %>% filter(year<yearMin) %>% group_by(quantileIndex) %>% dplyr::summarize(expPre1990=sum(n)) %>% pivot_longer(expPre1990,names_to='year',values_to='n') %>% right_join(expand_grid(quantileIndex=mar$quantileIndex,year='expPre1990'))%>% mutate(n=replace_na(n,0))
	
		tmp1=tmp1 %>% filter(year>=yearMin) %>%  mutate(year=paste0('X',year))
		# need to fill in empty years. make template and and join so the long format has all the years
		template.l=expand_grid(quantileIndex=mar$quantileIndex, year=paste0('X',yearMin:yearMax)) 
		tmp2=template.l %>% left_join(tmp1) %>% mutate(n=replace_na(n,0)) %>% bind_rows(	tmp1.pre90) %>% pivot_wider(names_from='year',values_from='n')
	
		toKeep = toKeep %>% left_join(tmp2,by='quantileIndex')
	
		#------------------------------
		# find marginal cells that were counted as exposed for year.q and cell.q. ONLY APPLIES TO year + cell
			# exposed years and cells
		tmp3a=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex,exposedYearCells) %>% unnest(exposedYearCells) 
			# marginal cells
		tmp3a.1=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex,marginalCells) %>% unnest(marginalCells) %>% rename(cell=marginalCells)
			#find when marginal contribute to exposure in each year
		new.vars=paste0('mar',yearMin:yearMax)
		tmp5=mar %>% filter(is.na(year.cell.q)) %>% dplyr::select(quantileIndex) %>% mutate(!!!setNames(rep(0, length(new.vars)), new.vars))
		if(nrow(tmp3a.1)>0){
			inBoth2=tmp3a %>% inner_join(tmp3a.1,by=c('quantileIndex','cell'))
			tmp4a=inBoth2 %>% filter(year>yearMin) %>%  mutate(year=paste0('mar',year)) %>% group_by(quantileIndex,year) %>% count
			templ=expand_grid(quantileIndex=mar$quantileIndex, year=paste0('mar',yearMin:yearMax)) 
			tmp5=templ%>% left_join(tmp4a) %>% mutate(n=replace_na(n,0)) %>% pivot_wider(names_from='year',values_from='n')
		}
		toKeep.y.c = toKeep %>% filter(is.na(year.cell.q)) %>% left_join(tmp5,by='quantileIndex')
		expNotMarg = (toKeep.y.c %>% dplyr::select(paste0('X',yearMin):paste0('X',yearMax))) - (toKeep.y.c %>% dplyr::select(paste0('mar',yearMin):paste0('mar',yearMax)))
		names(expNotMarg)=sub('X','XnM',names(expNotMarg))
		toKeep.y.c = 	toKeep.y.c %>% bind_cols(expNotMarg)
		tk=toKeep.y.c
		qsave(tk,file=out.f)
	
		# =======================
		# species-level summary statistics
		spStats=tk %>% dplyr::select(year.q:quantileIndex,rangeSize,X2015,X2016,X2019,X2020,X2023,XnM2015,XnM2016,XnM2019,XnM2020,XnM2023) %>% mutate(pExp2015=round(X2015/rangeSize,3),pExp2016=round(X2016/rangeSize,3),pExp2019=round(X2019/rangeSize,3),pExp2020=round(X2020/rangeSize,3),pExp2023=round(X2023/rangeSize,3))
		qsave(spStats,file=paste0(outDirSp,'/',basename(x1)))

	},mc.cores=mc.cores)
}

