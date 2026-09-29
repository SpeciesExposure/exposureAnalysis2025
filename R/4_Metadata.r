# 4_Metadata.r
# Builds the species attributes table, spAttributes_v8.qs: IUCN Red List
# category, ecoregions, taxonomy and range size for every species with a range.
# The table has one row per species and ecoregion.
#
# Missing taxonomy is looked up in the GBIF backbone and the Catalogue of Life
# through their web services, so a rerun can return different values if those
# services have changed.
#
# Inputs:  INPUTS$iucn_tables_dir, INPUTS$ecoregions_dir, species range files
# Outputs: metaDir/spAttributes_v8.qs, metaDir/rangeSize_v3.qs

if (!exists('PROJECT_PATHS')) source(file.path('config', 'paths.R'))
if (!exists('template')) source(file.path(PROJECT_PATHS$src_r, '1_Setup.r'))

#+++++++++++++++++++++++++++++++++++++++++++++++
# IUCN Status
f=list.files(INPUTS$iucn_tables_dir,pattern='assess',full.names=T,recursive=T)
# read the assessment tables of all groups
d1=lapply(f,function(x){
	data.frame(read.csv(x),group=basename(dirname(x)))
}) %>% bind_rows  %>% mutate(group=replace(group,group=='Birds1','Birds')) %>% mutate(group=replace(group,group=='Birds2','Birds'))
# relevant subset
d.iucn = d1 %>% mutate(sp=sub('[[:space:]]','_',scientificName),redlistCategory=sub('[[:space:]]','_',redlistCategory)) %>% dplyr::select(assessmentId:criteriaVersion,populationTrend,realm,sp,group) %>% dplyr::select(-scientificName)
hasMap=.be(list.files(EXTERNAL_PATHS$ranges_expert))

saveRDS(d.iucn,file=paste0(metaDir,'/iucnStatus.rds'))
d.iucn=readRDS(paste0(metaDir,'/iucnStatus.rds'))

# species attributes, first version
qs_save(d.iucn,paste0(metaDir,'/spAttributes_tmp.qs'))

#+++++++++++++++++++++++++++++++++++++++++++++++
# species by ecoregion  
e.r1=terra::rast(KEY_FILES$ecoregion_raster)
e.template=terra::rast(KEY_FILES$land_template)
eco_layer_idx=which(toupper(names(e.r1))=='ECO_ID')
if(!length(eco_layer_idx)) eco_layer_idx=grep('ECO',toupper(names(e.r1)))
if(!length(eco_layer_idx)) stop('Could not find ECO_ID layer in Ecoregions2017_v2.tif')
eco_id_r=e.r1[[eco_layer_idx[1]]]

ff=list.files(paste0(dataDir,'/spRangeTables_Prepped/Expert'),full.names=T)
out=mclapply(ff,function(x){
	print(x)
	s=qread(x) %>% unique %>% as.integer
	s=s[!is.na(s) & s>=1 & s<=terra::ncell(e.template)]
	if(!length(s)) return(NULL)
	xy=terra::xyFromCell(e.template,s)
	eco_vals=terra::extract(eco_id_r,xy)
	eco_col=setdiff(names(eco_vals),'ID')
	if(!length(eco_col)) return(NULL)
	tibble(ECO_ID=as.integer(eco_vals[[eco_col[1]]]),sp=basename(file_path_sans_ext(x))) %>% filter(!is.na(ECO_ID)) %>% unique
},mc.cores=mc.cores) %>% bind_rows

spXEco=out
e.d=readRDS(KEY_FILES$ecoregion_metadata) %>% dplyr::select(- (COLOR:LICENSE))

spXEco1= spXEco %>% left_join(e.d,by='ECO_ID') 

qs_save(spXEco1,file=paste0(metaDir,'/SpeciesXEcoregion.qs'))
spXEco=qread(paste0(metaDir,'/SpeciesXEcoregion.qs'))

# add ecoregions to the species attributes
spMeta=qread(paste0(metaDir,'/spAttributes_tmp.qs')) %>% full_join(spXEco1,by=c('sp'))

# patch: species still missing ecoregion metadata after initial join
sp_missing_eco=spMeta %>%
	group_by(sp) %>%
	summarise(has_eco=any(!is.na(ECO_ID) | (!is.na(ECO_NAME) & nzchar(trimws(ECO_NAME)))),.groups='drop') %>%
	filter(!has_eco) %>%
	pull(sp)

if(length(sp_missing_eco)){
	message('Species still missing ecoregion metadata after initial join: ',length(sp_missing_eco))
	if(!exists('e.poly')){
		e.poly=sf::st_read(KEY_FILES$ecoregion_shp,quiet=TRUE) %>%
			sf::st_make_valid() %>%
			sf::st_transform(4326) %>%
			dplyr::select(ECO_ID,geometry)
	}

	if(!'match_method' %in% names(spMeta)) spMeta$match_method=NA_character_
	if(!'assignment_confidence' %in% names(spMeta)) spMeta$assignment_confidence=NA_character_
	if(!'assignment_distance_m' %in% names(spMeta)) spMeta$assignment_distance_m=NA_real_
	if(!'REALM' %in% names(spMeta) && 'realm' %in% names(spMeta)) spMeta$REALM=spMeta$realm
	if(!'BIOME_NAME' %in% names(spMeta)) spMeta$BIOME_NAME=NA_character_

	merge_sp_assignments=function(spMeta,sp_assign){
		if(!nrow(sp_assign)) return(spMeta)
		merge_cols=intersect(names(sp_assign),names(spMeta))
		merge_cols=setdiff(merge_cols,'sp')
		recovered_sp=unique(sp_assign$sp)
		spMeta_base=spMeta %>% filter(sp %in% recovered_sp) %>% dplyr::select(-any_of(merge_cols)) %>% distinct
		bind_rows(
			spMeta %>% filter(!(sp %in% recovered_sp)),
			spMeta_base %>% inner_join(sp_assign,by='sp')
		)
	}

	spXEco_missing=mclapply(sp_missing_eco,function(spn){
		f.sp=file.path(dataDir,'spRangeTables_Prepped','Expert',paste0(spn,'.qs'))
		if(!file.exists(f.sp)) return(NULL)
		s=qread(f.sp) %>% unique %>% as.integer
		s=s[!is.na(s) & s>=1 & s<=terra::ncell(e.template)]
		if(!length(s)) return(NULL)
		xy=terra::xyFromCell(e.template,s)
		pts=sf::st_as_sf(data.frame(x=xy[,1],y=xy[,2]),coords=c('x','y'),crs=4326)
		hits=sf::st_intersects(pts,e.poly)
		if(!any(lengths(hits)>0)) return(NULL)
		eco_ids=unique(e.poly$ECO_ID[unique(unlist(hits))])
		tibble(sp=spn,ECO_ID=as.integer(eco_ids),match_method='direct_shp_overlap',assignment_confidence='high',assignment_distance_m=0)
	},mc.cores=mc.cores) %>% bind_rows

	if(nrow(spXEco_missing)){
		spXEco_missing=spXEco_missing %>% left_join(e.d,by='ECO_ID') %>% distinct(sp,ECO_ID,.keep_all=TRUE)
		write.csv(spXEco_missing %>% arrange(sp,ECO_ID),file=paste0(metaDir,'/species_missing_ecoregion_recovered_from_shp.csv'),row.names=FALSE)
		spXEco1=bind_rows(spXEco1,spXEco_missing) %>% distinct(sp,ECO_ID,.keep_all=TRUE)
		spMeta=merge_sp_assignments(spMeta,spXEco_missing)
		message('Recovered missing-ecoregion species from shp: ',length(unique(spXEco_missing$sp)))
	}else{
		message('Recovered missing-ecoregion species from shp: 0')
	}

	still_missing=spMeta %>%
		group_by(sp) %>%
		summarise(has_eco=any(!is.na(ECO_ID) | (!is.na(ECO_NAME) & nzchar(trimws(ECO_NAME)))),.groups='drop') %>%
		filter(!has_eco)

	# second-pass patch: assign nearest ecoregion if within 10 km
	if(nrow(still_missing)){
		e.poly_m=sf::st_transform(e.poly,6933)
		spXEco_near10=mclapply(still_missing$sp,function(spn){
			f.sp=file.path(dataDir,'spRangeTables_Prepped','Expert',paste0(spn,'.qs'))
			if(!file.exists(f.sp)) return(NULL)
			s=qread(f.sp) %>% unique %>% as.integer
			s=s[!is.na(s) & s>=1 & s<=terra::ncell(e.template)]
			if(!length(s)) return(NULL)
			xy=terra::xyFromCell(e.template,s)
			pts=sf::st_as_sf(data.frame(x=xy[,1],y=xy[,2]),coords=c('x','y'),crs=4326)
			pts_m=sf::st_transform(pts,6933)
			nn_idx=sf::st_nearest_feature(pts_m,e.poly_m)
			nn_dist=as.numeric(sf::st_distance(pts_m,e.poly_m[nn_idx,],by_element=TRUE))
			if(!length(nn_dist) || all(is.na(nn_dist))) return(NULL)
			i=which.min(nn_dist)
			if(!length(i) || is.na(nn_dist[i]) || nn_dist[i]>10000) return(NULL)
			tibble(sp=spn,ECO_ID=as.integer(e.poly$ECO_ID[nn_idx[i]]),assignment_distance_m=nn_dist[i],match_method='nearest_ecoregion_10km',assignment_confidence='medium')
		},mc.cores=mc.cores) %>% bind_rows

		if(nrow(spXEco_near10)){
			spXEco_near10=spXEco_near10 %>% left_join(e.d,by='ECO_ID') %>% distinct(sp,ECO_ID,.keep_all=TRUE)
			write.csv(spXEco_near10 %>% arrange(sp,ECO_ID),file=paste0(metaDir,'/species_missing_ecoregion_assigned_within10km.csv'),row.names=FALSE)
			spXEco1=bind_rows(spXEco1,spXEco_near10) %>% distinct(sp,ECO_ID,.keep_all=TRUE)
			spMeta=merge_sp_assignments(spMeta,spXEco_near10)
			message('Recovered missing-ecoregion species within 10 km: ',length(unique(spXEco_near10$sp)))
		}else{
			message('Recovered missing-ecoregion species within 10 km: 0')
		}

		still_missing=spMeta %>%
			group_by(sp) %>%
			summarise(has_eco=any(!is.na(ECO_ID) | (!is.na(ECO_NAME) & nzchar(trimws(ECO_NAME)))),.groups='drop') %>%
			filter(!has_eco)

		# third-pass patch: assign nearest ecoregion if within 30 km
		if(nrow(still_missing)){
			spXEco_near30=mclapply(still_missing$sp,function(spn){
				f.sp=file.path(dataDir,'spRangeTables_Prepped','Expert',paste0(spn,'.qs'))
				if(!file.exists(f.sp)) return(NULL)
				s=qread(f.sp) %>% unique %>% as.integer
				s=s[!is.na(s) & s>=1 & s<=terra::ncell(e.template)]
				if(!length(s)) return(NULL)
				xy=terra::xyFromCell(e.template,s)
				pts=sf::st_as_sf(data.frame(x=xy[,1],y=xy[,2]),coords=c('x','y'),crs=4326)
				pts_m=sf::st_transform(pts,6933)
				nn_idx=sf::st_nearest_feature(pts_m,e.poly_m)
				nn_dist=as.numeric(sf::st_distance(pts_m,e.poly_m[nn_idx,],by_element=TRUE))
				if(!length(nn_dist) || all(is.na(nn_dist))) return(NULL)
				i=which.min(nn_dist)
				if(!length(i) || is.na(nn_dist[i]) || nn_dist[i]>30000) return(NULL)
				tibble(sp=spn,ECO_ID=as.integer(e.poly$ECO_ID[nn_idx[i]]),assignment_distance_m=nn_dist[i],match_method='nearest_ecoregion_30km',assignment_confidence='low')
			},mc.cores=mc.cores) %>% bind_rows

			if(nrow(spXEco_near30)){
				spXEco_near30=spXEco_near30 %>% left_join(e.d,by='ECO_ID') %>% distinct(sp,ECO_ID,.keep_all=TRUE)
				write.csv(spXEco_near30 %>% arrange(sp,ECO_ID),file=paste0(metaDir,'/species_missing_ecoregion_assigned_within30km.csv'),row.names=FALSE)
				spXEco1=bind_rows(spXEco1,spXEco_near30) %>% distinct(sp,ECO_ID,.keep_all=TRUE)
				spMeta=merge_sp_assignments(spMeta,spXEco_near30)
				message('Recovered missing-ecoregion species within 30 km: ',length(unique(spXEco_near30$sp)))
			}else{
				message('Recovered missing-ecoregion species within 30 km: 0')
			}
		}
	}

	still_missing=spMeta %>%
		group_by(sp) %>%
		summarise(has_eco=any(!is.na(ECO_ID) | (!is.na(ECO_NAME) & nzchar(trimws(ECO_NAME)))),.groups='drop') %>%
		filter(!has_eco)

	# fallback for remaining metadata gaps: assign nearest BIOME/REALM using species coordinates
	sp_missing_biome_realm=spMeta %>%
		group_by(sp) %>%
		summarise(
			has_biome=any(!is.na(BIOME_NAME) & nzchar(trimws(BIOME_NAME))),
			has_realm=any(!is.na(REALM) & nzchar(trimws(REALM))),
			.groups='drop'
		) %>%
		filter(!has_biome | !has_realm)

	if(nrow(sp_missing_biome_realm)){
		sp_biome_realm=mclapply(sp_missing_biome_realm$sp,function(spn){
			f.sp=file.path(dataDir,'spRangeTables_Prepped','Expert',paste0(spn,'.qs'))
			if(!file.exists(f.sp)) return(NULL)
			s=qread(f.sp) %>% unique %>% as.integer
			s=s[!is.na(s) & s>=1 & s<=terra::ncell(e.template)]
			if(!length(s)) return(NULL)
			xy=terra::xyFromCell(e.template,s)
			pts=sf::st_as_sf(data.frame(x=xy[,1],y=xy[,2]),coords=c('x','y'),crs=4326)
			pts_m=sf::st_transform(pts,6933)
			nn_idx=sf::st_nearest_feature(pts_m,e.poly_m)
			nn_dist=as.numeric(sf::st_distance(pts_m,e.poly_m[nn_idx,],by_element=TRUE))
			if(!length(nn_dist) || all(is.na(nn_dist))) return(NULL)
			i=which.min(nn_dist)
			tibble(sp=spn,ECO_ID=as.integer(e.poly$ECO_ID[nn_idx[i]]),assignment_distance_m=nn_dist[i])
		},mc.cores=mc.cores) %>% bind_rows

		if(nrow(sp_biome_realm)){
			sp_biome_realm=sp_biome_realm %>%
				left_join(e.d %>% dplyr::select(any_of(c('ECO_ID','BIOME_NAME','REALM'))),by='ECO_ID') %>%
				dplyr::select(sp,BIOME_NAME,REALM,assignment_distance_m) %>%
				distinct(sp,.keep_all=TRUE) %>%
				mutate(match_method='biome_realm_only',assignment_confidence='coarse')

			spMeta=spMeta %>%
				left_join(sp_biome_realm %>% rename(BIOME_NAME_fallback=BIOME_NAME,REALM_fallback=REALM),by='sp') %>%
				mutate(
					BIOME_NAME=ifelse(is.na(BIOME_NAME) | !nzchar(trimws(BIOME_NAME)),BIOME_NAME_fallback,BIOME_NAME),
					REALM=ifelse(is.na(REALM) | !nzchar(trimws(REALM)),REALM_fallback,REALM),
					match_method=ifelse(sp %in% sp_biome_realm$sp & (is.na(match_method.x) | !nzchar(trimws(match_method.x))),match_method.y,match_method.x),
					assignment_confidence=ifelse(sp %in% sp_biome_realm$sp & (is.na(assignment_confidence.x) | !nzchar(trimws(assignment_confidence.x))),assignment_confidence.y,assignment_confidence.x),
					assignment_distance_m=ifelse(sp %in% sp_biome_realm$sp & is.na(assignment_distance_m.x),assignment_distance_m.y,assignment_distance_m.x)
				) %>%
				dplyr::select(-any_of(c('BIOME_NAME_fallback','REALM_fallback','match_method.x','match_method.y','assignment_confidence.x','assignment_confidence.y','assignment_distance_m.x','assignment_distance_m.y')))

			write.csv(sp_biome_realm %>% arrange(sp),file=paste0(metaDir,'/species_missing_ecoregion_biome_realm_only_fallback.csv'),row.names=FALSE)
			message('Assigned nearest-coordinate BIOME/REALM fallback for species: ',nrow(sp_biome_realm))
		}else{
			message('Assigned nearest-coordinate BIOME/REALM fallback for species: 0')
		}
	}

	still_missing=spMeta %>%
		group_by(sp) %>%
		summarise(has_eco=any(!is.na(ECO_ID) | (!is.na(ECO_NAME) & nzchar(trimws(ECO_NAME)))),.groups='drop') %>%
		filter(!has_eco)

	if(nrow(still_missing)){
		spMeta=spMeta %>% mutate(
			match_method=ifelse(sp %in% still_missing$sp & (is.na(match_method) | !nzchar(trimws(match_method))),'unresolved_no_ecoregion',match_method),
			assignment_confidence=ifelse(sp %in% still_missing$sp & (is.na(assignment_confidence) | !nzchar(trimws(assignment_confidence))),'none',assignment_confidence)
		)
	}

	write.csv(still_missing,file=paste0(metaDir,'/species_still_missing_ecoregion_after_all_fallbacks.csv'),row.names=FALSE)
	message('Species still missing ecoregion metadata after all fallbacks: ',nrow(still_missing))
}

# keep all species, including those with no assessment
qs_save(spMeta,paste0(metaDir,'/spAttributes_tmp.qs'))

#+++++++++++++++++++++++++++++++++++++++++++++++
# taxonomic group
tax.f=list.files(INPUTS$iucn_tables_dir,pattern='taxonomy',full.names=T,recursive=T)

d.tax=lapply(tax.f,function(x){
	data.frame(read.csv(x),group=basename(dirname(x)))
}) %>% bind_rows  %>% mutate(group=replace(group,group=='Birds1','Birds')) %>% mutate(group=replace(group,group=='Birds2','Birds')) %>% dplyr::select(-(infraType:taxonomicNotes)) %>% mutate(sp=sub('[[:space:]]','_',scientificName)) %>% dplyr::select(-scientificName,-group)
		
saveRDS(d.tax,file=paste0(metaDir,'/SpeciesXTaxonomy.rds'))
d.tax=readRDS(paste0(metaDir,'/SpeciesXTaxonomy.rds'))

# add taxonomy to the species attributes
spMeta=qread(paste0(metaDir,'/spAttributes_tmp.qs'))  %>% full_join(d.tax,by=c('sp','internalTaxonId'))

# taxon group for species that have none in the IUCN tables, read from a file in data/
llm_tax_input_file=paste0(metaDir,'/species_missing_taxon_group_v1.csv')
llm_tax_filled_file=file.path(PROJECT_PATHS$data,'species_missing_taxon_group_filled.csv')

normalize_taxon_group <- function(x){
	x=trimws(as.character(x))
	xl=tolower(x)
	out=ifelse(xl %in% c('bird','birds','aves'),'Birds',
		ifelse(xl %in% c('mammal','mammals','mammalia'),'Mammals',
			ifelse(xl %in% c('bat','bats','chiroptera'),'Mammals',
				ifelse(xl %in% c('reptile','reptiles','reptilia'),'Reptiles',
					ifelse(xl %in% c('amphibian','amphibians','amphibia'),'Amphibians',x)))))
	out[out=='']=NA
	out
}

missing_tax_for_llm=spMeta %>%
	filter(is.na(group) | trimws(group)=='') %>%
	group_by(sp) %>%
	summarise(
		internalTaxonId=dplyr::first(internalTaxonId),
		kingdomName=dplyr::first(kingdomName),
		phylumName=dplyr::first(phylumName),
		className=dplyr::first(className),
		orderName=dplyr::first(orderName),
		familyName=dplyr::first(familyName),
		genusName=dplyr::first(genusName),
		group_current=dplyr::first(group),
		.groups='drop'
	)
write.csv(missing_tax_for_llm,file=llm_tax_input_file,row.names=FALSE)
message('Wrote species missing taxon group: ',nrow(missing_tax_for_llm))

if(file.exists(llm_tax_filled_file)){
	llm_fill=read.csv(llm_tax_filled_file,stringsAsFactors=FALSE)
	sp_col=intersect(c('sp','spName','species'),names(llm_fill))[1]
	grp_col=intersect(c('group','group_llm','taxon_group','taxon_group_llm','group2'),names(llm_fill))[1]
	if(!is.na(sp_col) & !is.na(grp_col)){
		llm_fill=llm_fill %>%
			transmute(sp=.data[[sp_col]],group_llm=normalize_taxon_group(.data[[grp_col]])) %>%
			filter(!is.na(sp),!is.na(group_llm)) %>%
			distinct(sp,.keep_all=TRUE)

		pre_missing=sum(is.na(spMeta$group) | trimws(spMeta$group)=='')
		spMeta=spMeta %>%
			left_join(llm_fill,by='sp') %>%
			mutate(group=ifelse(is.na(group) | trimws(group)=='',group_llm,group)) %>%
			dplyr::select(-group_llm)
		post_missing=sum(is.na(spMeta$group) | trimws(spMeta$group)=='')
		message('Taxon-group rows filled from file: ',pre_missing-post_missing)
		message('Rows still missing taxon group after fill from file: ',post_missing)
	}else{
		warning('Taxon-group file found but expected columns missing (species + group).')
	}
}

# web fallback for missing family and higher taxonomy (GBIF backbone)
tax_cols_needed=c('familyName','orderName','className','phylumName','kingdomName')
tax_cols_needed=intersect(tax_cols_needed,names(spMeta))

if(length(tax_cols_needed)){
	missing_tax_sp=spMeta %>%
		group_by(sp) %>%
		summarise(missing_any=any(if_any(all_of(tax_cols_needed),~is.na(.) | trimws(as.character(.))=='')),.groups='drop') %>%
		filter(missing_any) %>%
		pull(sp)

	write.csv(data.frame(sp=missing_tax_sp),file=paste0(metaDir,'/species_missing_taxonomy_for_web_lookup_v1.csv'),row.names=FALSE)

	if(length(missing_tax_sp) && requireNamespace('jsonlite',quietly=TRUE)){
		lookup_gbif_taxonomy <- function(spn){
			qname=gsub('_',' ',spn)
			url=paste0('https://api.gbif.org/v1/species/match?name=',utils::URLencode(qname,reserved=TRUE))
			res=try(jsonlite::fromJSON(url),silent=TRUE)
			if(inherits(res,'try-error') || is.null(res)){
				return(data.frame(
					sp=spn,
					familyName=NA_character_,
					orderName=NA_character_,
					className=NA_character_,
					phylumName=NA_character_,
					kingdomName=NA_character_,
					gbif_usageKey=NA_integer_,
					gbif_matchType=NA_character_,
					gbif_status=NA_character_,
					stringsAsFactors=FALSE
				))
			}
			data.frame(
				sp=spn,
				familyName=ifelse(is.null(res$family),NA_character_,as.character(res$family)),
				orderName=ifelse(is.null(res$order),NA_character_,as.character(res$order)),
				className=ifelse(is.null(res$class),NA_character_,as.character(res$class)),
				phylumName=ifelse(is.null(res$phylum),NA_character_,as.character(res$phylum)),
				kingdomName=ifelse(is.null(res$kingdom),NA_character_,as.character(res$kingdom)),
				gbif_usageKey=ifelse(is.null(res$usageKey),NA_integer_,as.integer(res$usageKey)),
				gbif_matchType=ifelse(is.null(res$matchType),NA_character_,as.character(res$matchType)),
				gbif_status=ifelse(is.null(res$status),NA_character_,as.character(res$status)),
				stringsAsFactors=FALSE
			)
		}

		gbif_tax=lapply(missing_tax_sp,lookup_gbif_taxonomy) %>% bind_rows
		write.csv(gbif_tax,file=paste0(metaDir,'/taxonomy_web_lookup_gbif_v1.csv'),row.names=FALSE)

		gbif_fill=gbif_tax %>%
			filter(if_any(all_of(tax_cols_needed),~!is.na(.) & trimws(as.character(.))!='')) %>%
			distinct(sp,.keep_all=TRUE)

		pre_missing_tax=spMeta %>%
			group_by(sp) %>%
			summarise(missing_any=any(if_any(all_of(tax_cols_needed),~is.na(.) | trimws(as.character(.))=='')),.groups='drop') %>%
			filter(missing_any) %>%
			nrow

		spMeta=spMeta %>%
			left_join(gbif_fill %>% dplyr::select(sp,all_of(tax_cols_needed)),by='sp',suffix=c('','.gbif'))

		for(col in tax_cols_needed){
			col_gbif=paste0(col,'.gbif')
			spMeta[[col]]=ifelse(is.na(spMeta[[col]]) | trimws(as.character(spMeta[[col]]))=='',spMeta[[col_gbif]],spMeta[[col]])
		}

		spMeta=spMeta %>% dplyr::select(-any_of(paste0(tax_cols_needed,'.gbif')))

		post_missing_tax=spMeta %>%
			group_by(sp) %>%
			summarise(missing_any=any(if_any(all_of(tax_cols_needed),~is.na(.) | trimws(as.character(.))=='')),.groups='drop') %>%
			filter(missing_any) %>%
			nrow

		message('GBIF web taxonomy backfill species with any taxonomy found: ',nrow(gbif_fill))
		message('Species still missing family/higher taxonomy after GBIF backfill: ',post_missing_tax)
		message('Species reduced by GBIF backfill: ',pre_missing_tax-post_missing_tax)

		# second web fallback for remaining missing taxonomy (Catalogue of Life)
		remaining_missing_sp=spMeta %>%
			group_by(sp) %>%
			summarise(missing_any=any(if_any(all_of(tax_cols_needed),~is.na(.) | trimws(as.character(.))=='')),.groups='drop') %>%
			filter(missing_any) %>%
			pull(sp)

		if(length(remaining_missing_sp)){
			lookup_col_taxonomy <- function(spn){
				qname=gsub('_',' ',spn)
				url=paste0(
					'https://api.checklistbank.org/dataset/3LR/nameusage/search?q=',
					utils::URLencode(qname,reserved=TRUE),
					'&type=EXACT&rank=SPECIES'
				)
				res=try(jsonlite::fromJSON(url),silent=TRUE)
				if(inherits(res,'try-error') || is.null(res) || is.null(res$result) || !length(res$result)){
					return(data.frame(
						sp=spn,
						familyName=NA_character_,
						orderName=NA_character_,
						className=NA_character_,
						phylumName=NA_character_,
						kingdomName=NA_character_,
						col_usageId=NA_character_,
						col_status=NA_character_,
						stringsAsFactors=FALSE
					))
				}

				x=res$result[1,]
				cl=x$classification[[1]]
				get_rank_name <- function(rank_nm){
					if(is.null(cl) || !is.data.frame(cl) || !nrow(cl)) return(NA_character_)
					i=which(tolower(as.character(cl$rank))==tolower(rank_nm))
					if(!length(i)) return(NA_character_)
					as.character(cl$name[i[1]])
				}

				data.frame(
					sp=spn,
					familyName=get_rank_name('family'),
					orderName=get_rank_name('order'),
					className=get_rank_name('class'),
					phylumName=get_rank_name('phylum'),
					kingdomName=get_rank_name('kingdom'),
					col_usageId=ifelse(is.null(x$usage$id[[1]]),NA_character_,as.character(x$usage$id[[1]])),
					col_status=ifelse(is.null(x$usage$status[[1]]),NA_character_,as.character(x$usage$status[[1]])),
					stringsAsFactors=FALSE
				)
			}

			col_tax=lapply(remaining_missing_sp,lookup_col_taxonomy) %>% bind_rows
			write.csv(col_tax,file=paste0(metaDir,'/taxonomy_web_lookup_col_v1.csv'),row.names=FALSE)

			col_fill=col_tax %>%
				filter(if_any(all_of(tax_cols_needed),~!is.na(.) & trimws(as.character(.))!='')) %>%
				distinct(sp,.keep_all=TRUE)

			pre_missing_tax_col=spMeta %>%
				group_by(sp) %>%
				summarise(missing_any=any(if_any(all_of(tax_cols_needed),~is.na(.) | trimws(as.character(.))=='')),.groups='drop') %>%
				filter(missing_any) %>%
				nrow

			spMeta=spMeta %>%
				left_join(col_fill %>% dplyr::select(sp,all_of(tax_cols_needed)),by='sp',suffix=c('','.col'))

			for(col in tax_cols_needed){
				col_col=paste0(col,'.col')
				spMeta[[col]]=ifelse(is.na(spMeta[[col]]) | trimws(as.character(spMeta[[col]]))=='',spMeta[[col_col]],spMeta[[col]])
			}

			spMeta=spMeta %>% dplyr::select(-any_of(paste0(tax_cols_needed,'.col')))

			post_missing_tax_col=spMeta %>%
				group_by(sp) %>%
				summarise(missing_any=any(if_any(all_of(tax_cols_needed),~is.na(.) | trimws(as.character(.))=='')),.groups='drop') %>%
				filter(missing_any) %>%
				nrow

			message('Catalogue of Life taxonomy backfill species with any taxonomy found: ',nrow(col_fill))
			message('Species still missing family/higher taxonomy after CoL backfill: ',post_missing_tax_col)
			message('Species reduced by CoL backfill: ',pre_missing_tax_col-post_missing_tax_col)
		}
	}else if(length(missing_tax_sp)){
		warning('jsonlite not available; wrote species_missing_taxonomy_for_web_lookup_v1.csv for manual web fill.')
	}
}

# taxonomy-derived fallback for missing taxon group (uses CoL/GBIF-enriched taxonomy)
infer_group_from_taxonomy <- function(className,orderName,phylumName,kingdomName){
	cls=tolower(trimws(as.character(className)))
	ord=tolower(trimws(as.character(orderName)))
	phyl=tolower(trimws(as.character(phylumName)))
	king=tolower(trimws(as.character(kingdomName)))

	out=dplyr::case_when(
		cls %in% c('aves') ~ 'Birds',
		cls %in% c('mammalia') ~ 'Mammals',
		cls %in% c('reptilia') ~ 'Reptiles',
		cls %in% c('amphibia') ~ 'Amphibians',
		cls %in% c('squamata','testudines') | ord %in% c('squamata','testudines') ~ 'Reptiles',
		ord %in% c('chiroptera') ~ 'Mammals',
		phyl %in% c('chordata') & king %in% c('animalia') ~ NA_character_,
		TRUE ~ NA_character_
	)
	out[out=='']=NA_character_
	out
}

group_missing_before=sum(is.na(spMeta$group) | trimws(as.character(spMeta$group))=='')

spMeta=spMeta %>%
	mutate(group_missing_before_flag=is.na(group) | trimws(as.character(group))=='') %>%
	mutate(group_taxonomy=infer_group_from_taxonomy(className,orderName,phylumName,kingdomName)) %>%
	mutate(group=ifelse(group_missing_before_flag,group_taxonomy,group)) %>%
	mutate(group_filled_from_taxonomy_flag=group_missing_before_flag & !is.na(group) & trimws(as.character(group))!='')

group_missing_after=sum(is.na(spMeta$group) | trimws(as.character(spMeta$group))=='')
message('Taxonomy-derived taxon-group replacement rows filled: ',group_missing_before-group_missing_after)
message('Rows still missing taxon group after taxonomy-derived replacement: ',group_missing_after)

group_filled_from_taxonomy=spMeta %>%
	filter(group_filled_from_taxonomy_flag) %>%
	group_by(sp) %>%
	summarise(
		group=dplyr::first(group),
		className=dplyr::first(className),
		orderName=dplyr::first(orderName),
		phylumName=dplyr::first(phylumName),
		kingdomName=dplyr::first(kingdomName),
		.groups='drop'
	)
write.csv(group_filled_from_taxonomy,file=paste0(metaDir,'/species_taxon_group_filled_from_taxonomy_v1.csv'),row.names=FALSE)

group_still_missing=spMeta %>%
	filter(is.na(group) | trimws(as.character(group))=='') %>%
	group_by(sp) %>%
	summarise(
		className=dplyr::first(className),
		orderName=dplyr::first(orderName),
		phylumName=dplyr::first(phylumName),
		kingdomName=dplyr::first(kingdomName),
		.groups='drop'
	)
write.csv(group_still_missing,file=paste0(metaDir,'/species_still_missing_taxon_group_after_taxonomy_fill_v1.csv'),row.names=FALSE)

spMeta=spMeta %>% dplyr::select(-group_missing_before_flag,-group_taxonomy,-group_filled_from_taxonomy_flag)

qs_save(spMeta,paste0(metaDir,'/spAttributes_tmp2.qs'))

#+++++++++++++++++++++++++++++++++++++++++++++++
# range size 
# use the same range files as the exposure calculation
ff=list.files(paste0(dataDir,'/spRangeTables_Prepped/Expert'),full.names=T)

ar=cellSize(template)
rangeSize=mclapply(ff,function(x){ # x=ff[28954]
	print(x)
	r=qread(x)
	tibble(rangeSize=length(r),rangeSizeKm= round(sum(ar[r],na.rm=T)), spName=basename(file_path_sans_ext(x)))
},mc.cores=mc.cores) %>% bind_rows

qs_save(rangeSize,paste0(metaDir,'/rangeSize_v3.qs'))
spMeta=qread(paste0(metaDir,'/spAttributes_tmp2.qs'))  %>%  rename(spName = sp) %>% full_join(rangeSize,by=c('spName'))

# record and remove species with zero land cells in Expert range tables
sp_missing_zero_land=spMeta %>% filter(!is.na(rangeSize) & rangeSize==0)
write.csv(sp_missing_zero_land,file=paste0(metaDir,'/species_missing_zero_land_cells_fullinfo.csv'),row.names=FALSE)
qs_save(sp_missing_zero_land,file=paste0(metaDir,'/species_missing_zero_land_cells_fullinfo.qs'))

spMeta=spMeta %>% filter(is.na(rangeSize) | rangeSize>0)

# consolidate realm naming: keep canonical REALM in final output
if('realm' %in% names(spMeta)){
	if(!'REALM' %in% names(spMeta)) spMeta$REALM=NA_character_
	spMeta$REALM=ifelse(is.na(spMeta$REALM) | !nzchar(trimws(as.character(spMeta$REALM))),spMeta$realm,spMeta$REALM)
	spMeta=spMeta %>% dplyr::select(-any_of('realm'))
}

# Normalize taxonomy name columns to title case (IUCN returns ALL-CAPS; GBIF/CoL return title case)
tax_name_cols <- intersect(c('familyName','orderName','className','phylumName','kingdomName','genusName'), names(spMeta))
for (col in tax_name_cols) {
	spMeta[[col]] <- tools::toTitleCase(tolower(trimws(as.character(spMeta[[col]]))))
	spMeta[[col]][spMeta[[col]] == ''] <- NA_character_
}

qs_save(spMeta,paste0(metaDir,'/spAttributes_v8.qs'))

