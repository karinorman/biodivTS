##======================================================================
##	Before rarefaction we classify studies as having samples
##	from one or more locations in space: Single Location (SL) and
##	many location (ML). We use this classification, and the mean size of
##  SL studies to grid all other ML studies to a similar extent, using dggrid.
##  New cell centres are estimated based on rarefyID (study_cell).
##  Gridded studies are then further filtered based on coverage, omitting years
##  below the threshold for abundance data, and entire "rarefyID"s below the
##  threshold for presence and biomass data.

# Primary Code Authors: Shane Blowes and Sarah Supp
# Email: sablowes@gmail.com, sarah@weecology.org
# DISCLAIMER: code is under development
# 07/03/2017

# Input is the BioTIME database and metadata

# Coverage:
# Abundance-based coverage is calculated as C-hat corrected (Eqn 4a, Chao & Jost, Ecology)
# Incidence-based coverage is calculated as C-hat corrected, using years as samples
#   (Eqn Chat_sample(T) from Table 2, Chao et al. 2014)
# Note that the resolution of coverage (and omission) differs for abundance data vs. presence and biomass data

# Output is a gridded version of the BioTIME database, with a new column for rarefyID (the "studyID" to use henceforth)
#  that has been filtered for adequate abundance- or incidence-based coverage.
#   BioTIME_grid_filtered.Rdata (data.frame to be used for rarefaction in rarefy_griddedData_allmetrics.R)
##======================================================================

library(tidyverse)
library(dggridR)
library(scales)
library(vegan)
library(iNEXT)

##	get the data and metadata
##	Get the raw data locally
bt <- #read.csv("data/biotime_query.csv") %>%
  read.csv(system.file("extdata", "biotime/BioTIMEQuery02_04_2018 2.csv", package = "biodivTS")) %>%
  rename_with(toupper)
##	Get the meta data locally
meta <- read.csv(system.file("extdata", "biotime/BioTIMEMetadata_02_04_2018.csv", package = "biodivTS")) %>%
  rename_with(toupper)

##	join abundance records with the metadata
bt <- inner_join(meta, bt, by='STUDY_ID')

# Add columns for Species species name
bt <- bt %>% unite(col=Species, GENUS, SPECIES, remove=FALSE)

##======================================================================
## 	IDENTIFY SL & ML studies in the data:
##	Single Location [SL] have only 1 geographic coordinate for the whole study
##	Many Locations [ML] have many coordinates
meta <- meta %>%
  mutate(StudyMethod = ifelse(NUMBER_LAT_LONG == 1, "SL", NA))

bt <- bt %>%
  mutate(StudyMethod = ifelse(NUMBER_LAT_LONG == 1, 'SL', 'ML'))

##	change ML studies with small extent to SL
##	(i.e., ML studies with extent < mean(extent)+sd(extent) of SL studies)

##  First, look at frequency distribution of SL studies, and check for any extreme outliers
##	  SB: The extent of STUDY_ID==337 is 15259km2. All other SL have extent ≤ 500km2
##	  Including this study increases the mean by an order of magnitude (110km2 vs 11km2)
##	  and the sd goes from 53km2 to 1233km2. Our grid increases from ~32km2 to ~864km2
ggplot(meta[meta$StudyMethod=="SL",], aes(AREA_SQ_KM)) + geom_histogram()
ggplot(meta[meta$StudyMethod=="SL",], aes(AREA_SQ_KM)) + geom_histogram(aes(fill = as.numeric(AREA_SQ_KM) > 500 )) + scale_x_log10()
meta %>% filter(StudyMethod=='SL' & AREA_SQ_KM > 500) %>% select(TITLE, STUDY_ID, NUMBER_LAT_LONG, AREA_SQ_KM)

##  Calculate the extent and mean for SL studies, without the outlier
SL_extent_mean <- meta %>% filter(StudyMethod=='SL' & AREA_SQ_KM<=500) %>%
	summarise(extent_mean = mean(AREA_SQ_KM, na.rm=TRUE)) %>% .$extent_mean
SL_extent_sd <- meta %>% filter(StudyMethod=='SL' & AREA_SQ_KM<=500) %>%
	summarise(extent_sd = sd(AREA_SQ_KM, na.rm=TRUE)) %>% .$extent_sd

##	change the MLs into SLs that satisfy the criterion (< mean + sd)
bt$StudyMethod <- with(bt, ifelse(AREA_SQ_KM < (SL_extent_mean+SL_extent_sd), 'SL', StudyMethod))

##	We'll want to 'grid' SL and ML studies differently, add new coords to dataframe
## 	If StudyMethod=='SL', want to use the central lat and long, else
##	ML uses observed coords for each observation
bt <- bt %>%
	mutate(lon_to_grid = ifelse(StudyMethod=='SL', CENT_LONG, LONGITUDE),
		lat_to_grid = ifelse(StudyMethod=='SL', CENT_LAT, LATITUDE))

##======================================================================
##	other checks to filter data before gridding?
## 	studies with only ONE YEAR of data (since we do rarefaction on a yearly scale)
oneyear <- bt %>%
  group_by(STUDY_ID) %>%
  filter(max(YEAR)-min(YEAR)==0) %>%
  summarise() %>%
  collect %>% .[["STUDY_ID"]]
bt <- bt %>% filter(!(STUDY_ID %in% oneyear))

##======================================================================
##	create a global grid with cells approximately equal to extent +/- sd of the 'true' SL studies?

dgg <- dgconstruct(res=12)

## 	determine the resolution closest to our cutoff point for SL vs ML studies
res <- dg_closest_res_to_area(dgg, SL_extent_mean+SL_extent_sd)

##	set the resolution
dgg <- dgsetres(dgg, res)

##	get the corresponding grid cells for all observations
bt <- bt %>% mutate(cell = dgtransform(dgg, lat=lat_to_grid, lon=lon_to_grid))

##	what just happened?
check <- bt %>% group_by(StudyMethod, STUDY_ID) %>% summarise(n_cell = n_distinct(cell))

##	do all SL studies have one grid cell?
if (sum(dplyr::filter(check, StudyMethod=='SL') %>% .$n_cell != 1)==0) {
  print("all SL studies have 1 grid cell")
  } else { print("ERROR: some SL studies have > 1 grid cell") }

##	ok, how many cells/year/study are there? (e.g. how spread out were samples in a given study and year?)
check2 <- bt %>%
  group_by(StudyMethod, STUDY_ID, YEAR) %>%
  summarise(n_cell = n_distinct(cell))

range(check2$n_cell)

##	Some of these studies are spread out across MANY cells (>5000!?)
ggplot(filter(check2, StudyMethod=='ML')) +
	geom_histogram(aes(n_cell, fill = n_cell > 1000), binwidth=500) + xlab("Number of cells in a study per year") + scale_y_log10()

##======================================================================
##	in order to calculate new centres we need this new bt object
##	add rarefyID:
bt <- bt %>% unite(col=rarefyID, STUDY_ID, cell, sep="_", remove=FALSE)
#save(bt, dgg, file='biotime_cells.Rdata')
#load('~/Desktop/current/BioTime/data/biotime_cells_res11.Rdata')

##	want to calculate a single coordinate (i.e., the centre) for each rarefyID to place rarefyID's
##	within ecoregions....
##	also calculate cell_extent = area in convex hull (polygon) of points within a rarefyID. We can use this
##	to calculate a new 'extent' for each study AreaSum = sum(cell_extent) Issue #7 github
rarefyID_coords_nest <- ungroup(bt) %>%
	##	we don't need to do anything with the SL studies
	filter(StudyMethod!='SL') %>%
	##	select columns
	select(STUDY_ID, rarefyID, LONGITUDE, LATITUDE) %>%
	##	retain only uniqe locations within rarefyIDs (may want to change this if we want to weight the calculation of the centre)
	distinct(rarefyID, LONGITUDE, LATITUDE, .keep_all=TRUE) %>%
	##	there are also some rarefyID's with only one unique geographic
	group_by(rarefyID) %>%
	mutate(n_locations = n_distinct(LONGITUDE,LATITUDE)) %>%
	ungroup() %>%
	##	drop rarefyIDs with only one location
	filter(n_locations > 1) %>%
	##	drop our location counter
	select(-n_locations) %>%
	##	group & nest
	group_by(STUDY_ID, rarefyID) %>%
	nest()

##	i can't get purrr::map and chull to play nice so a loop it is! This takes awhile....
cell_extent <- numeric()
centre_rarefyID_x <- numeric()
centre_rarefyID_y <- numeric()
vertices_check <- data.frame()
for(i in 1:nrow(rarefyID_coords_nest)){
	##	sanity check
    print(paste('rarefyID', i, 'out of', length(unique(rarefyID_coords_nest$rarefyID))))
    ##	put a convex hull around the coords
	hull = chull(x=unlist(rarefyID_coords_nest$data[[i]][,'LONGITUDE']), y=unlist(rarefyID_coords_nest $data[[i]][,'LATITUDE']))
	##	get the vertices of the convex hull
	vertices = rarefyID_coords_nest$data[[i]][hull,c('LONGITUDE', 'LATITUDE')]
	##	put some metadata together for checking later
	info = cbind.data.frame(Realm=rep(rarefyID_coords_nest$STUDY_ID[i], times=nrow(vertices)), rarefyID=rep(rarefyID_coords_nest$rarefyID[i], times=nrow(vertices)), vertices)
	vertices_check = rbind.data.frame(vertices_check, info)	# this could be used to check all these convex hulls
	##	calculate the extent and centres (NB cell_extent==0 if there are only two points)
	cell_extent[i] = geosphere::areaPolygon(data.frame(x=vertices$LONGITUDE, y=vertices$LATITUDE))	# km2
	centre_rarefyID_x[i] = geosphere::geomean(cbind(x=vertices$LONGITUDE, y=vertices$LATITUDE))[1]
	centre_rarefyID_y[i] = geosphere::geomean(cbind(x=vertices$LONGITUDE, y=vertices$LATITUDE))[2]
}

##	combine STUDY_ID, rarefyID and the new cell_extent, and geographic centres
rarefyID_cell_centre <- cbind.data.frame(rarefyID_coords_nest[,1:2], cell_extent, rarefyID_x=centre_rarefyID_x, rarefyID_y=centre_rarefyID_y)
rarefyID_cell_centre <- as_tibble(rarefyID_cell_centre)

##	need to combine with the SL studies and ML with only 1 location per rarefyID
SL_coords <- ungroup(bt) %>%
	##	get the SL studies
	filter(StudyMethod=='SL') %>%
	##	select columns (use the central coords)
	select(STUDY_ID, rarefyID, CENT_LONG, CENT_LAT) %>%
	##	create cell_extent column, and rename coords for joining with rarefyID_cell_centre
	mutate(cell_extent = 0,
		rarefyID_x = CENT_LONG,
		rarefyID_y = CENT_LAT) %>%
	select(-CENT_LONG, -CENT_LAT)

ML_coords <- ungroup(bt) %>%
	filter(StudyMethod!='SL') %>%
	##	select columns
	select(STUDY_ID, rarefyID, LONGITUDE, LATITUDE) %>%
	##	retain only uniqe locations within rarefyIDs (may want to change this if we want to weight the calculation of the centre)
	distinct(rarefyID, LONGITUDE, LATITUDE, .keep_all=TRUE) %>%
	##	there are also some rarefyID's with only one unique geographic
	group_by(rarefyID) %>%
	mutate(n_locations = n_distinct(LONGITUDE,LATITUDE)) %>%
	ungroup() %>%
	##	retain rarefyIDs with only one location
	filter(n_locations == 1) %>%
	##	create cell_extent column, and rename coords for joining with rarefyID_cell_centre
	mutate(cell_extent = 0,
		rarefyID_x = LONGITUDE,
		rarefyID_y = LATITUDE) %>%
	select(-LONGITUDE, -LATITUDE, -n_locations)

##	put them together
rarefyID_cell_centre <- bind_rows(rarefyID_cell_centre, SL_coords, ML_coords)
##	not sure why but I have multiple entries for each rarefyID!
rarefyID_cell_centre <- rarefyID_cell_centre %>% distinct(STUDY_ID, rarefyID, cell_extent, rarefyID_x, rarefyID_y)

##======================================================================
##	reduce to data required for rarefying, and rename a couple of columns
##	NB: we will rarefy to the smallest number of ObsEventID's within studies
bt_grid <- bt %>%
  dplyr::select(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, BIOMASS_TYPE, STUDY_ID, YEAR, PLOT,
         cell, Species, Abundance = SUM.ALLRAWDATA.ABUNDANCE, Biomass = SUM.ALLRAWDATA.BIOMASS)

# add column for rarefyID (STUDY_ID + cell = samples within a study that are spatially grouped)
bt_grid <- bt_grid %>% unite(col=rarefyID, STUDY_ID, cell, sep="_", remove=FALSE)

# add column for ObsEventID (STUDY_ID + PLOT + YEAR = discrete sampling events within a year)
bt_grid <- bt_grid %>% unite(col=ObsEventID, rarefyID, PLOT, YEAR, sep="_", remove=FALSE)

bt_grid_collate <- bt_grid %>%
  group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, BIOMASS_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
  summarise(
    Abundance = sum(as.numeric(Abundance), na.rm=TRUE),
    Biomass = sum(as.numeric(Biomass), na.rm=TRUE)) %>%
  ungroup() %>%
  group_by(rarefyID) %>%
  filter(n_distinct(YEAR) > 1)

#save locally
usethis::use_data(bt_grid_collate)

##==================================== EXPLORE AND FILTER ON COVERAGE ========================================
# ##	first collate taxa within cells, to calculate coverage for count data (see below for incidence)
# bt_grid_collate <- bt_grid %>%
#   group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
#   summarise(
#     Abundance = sum(as.numeric(Abundance), na.rm=TRUE),
#     Biomass = sum(as.numeric(Biomass), na.rm=TRUE)) %>%
#   ungroup()
#
# abund_coverage <- bt_grid_collate %>%
#   # get only rows representing count data (abundance)
#   filter(ABUNDANCE_TYPE!='Presence/Absence' & ABUNDANCE_TYPE!='<NA>') %>%
#   # remove zeroes and NAs(some studies are designated count, but record Biomass estimates)
#   filter(Abundance > 0 & !is.na(Abundance)) %>%
#   group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, rarefyID, STUDY_ID, YEAR, cell) %>%
#   summarise(
#     # how many singletons
#     singletons = sum(Abundance==1),
#     # how many doubletons
#     doubletons = sum(Abundance==2),
#     # how many individuals in total sample
#     N = sum(Abundance),
#     # eqn 4a in Chao & Jost 2012 Ecology (==eqn 12 in Chao et al 2014 Ecol Monogr)
#     Chat = 1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)),
#     # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
#     # correction for communities with no doubletons (this prevents NaN results for either singletons==0 or doubletons==0)
#     Chat_corrected = ifelse(doubletons>0,
#                             1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2*doubletons)),
#                             1 - (singletons/N) * (((N-1)*singletons)/((N-1)*singletons + 2))),
#     # the iNEXT coverage calculation has some extras in it that I don't understand
#     # e.g., it corrects using f0.hat, the number of unobserved species (Chao1 estimator)? Why? Exactly how? Ref?
#     f1_iNEXT = DataInfo(Abundance)$f1,
#     f2_iNEXT = DataInfo(Abundance)$f2,
#     #n_iNEXT = DataInfo(Abundance)$n,
#     coverage = DataInfo(Abundance)$SC) %>%
#   ungroup()
#
# # are my estimates of singletons and doubletons equal to the iNEXT estimates?
# sum(abund_coverage$singletons!=abund_coverage$f1_iNEXT)		# yes, good
# sum(abund_coverage$doubletons!=abund_coverage$f2_iNEXT)		# yes, good
#
# #png('count_coverage_comparison_cell_level', width=800, height=600)
# with(abund_coverage, plot(x=Chat, y=coverage, type='n', xlim=c(0,1), ylim=c(0, 1),
#                                 xlab='Chat or Chat_corrected', ylab='Coverage (iNEXT)', main='Count data coverage calculation check'))
# with(abund_coverage, points(x=Chat, y=coverage))
# with(abund_coverage, points(x=Chat_corrected, y=coverage, col=2))
# abline(c(0,1), lty=2)
# legend('bottomright', legend=c('eqn4', 'eqn4 corrected for f2=0'), col=c(1,2), pch=1, lty=1)
# #dev.off()
#
# mn=mean(abund_coverage$Chat_corrected, na.rm=TRUE)
# sd=sd(abund_coverage$Chat_corrected, na.rm=TRUE)
# ggplot(abund_coverage, aes(Chat_corrected)) + geom_histogram(aes(fill=Chat_corrected>=mn-sd), binwidth=0.05) + geom_vline(xintercept=mn)
#
# #-------------------
# ## PRESENCE/ABSENCE coverage (for INCIDENCE DATA)... by rarefyID (year as samples)
# bt_grid_collate_incidence <- bt_grid %>%
#   # incidence data and remove NAs
#   filter(ABUNDANCE_TYPE=='Presence/Absence' & ABUNDANCE_TYPE!='<NA>') %>%
#   filter(!is.na(Abundance)) %>%
#   #took out ObsEventID and PLOT as groups b/c no distinct plots for P/A data (Each year has only 1 sample)
#   group_by(CLIMATE, REALM, TAXA, StudyMethod, ABUNDANCE_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
#   summarise(
#     Abundance = sum(Abundance, na.rm=TRUE)) %>%
#   ungroup() %>%
#   dplyr::select(-CLIMATE, -REALM, -TAXA, -StudyMethod, -ABUNDANCE_TYPE) %>%
#   # create incidence column
#   mutate(incidence = ifelse(Abundance==0, 0, 1))
#
# pa_coverage <- data.frame()
# for (r in unique(bt_grid_collate_incidence$rarefyID)){
#   dat <- bt_grid_collate_incidence[bt_grid_collate_incidence$rarefyID == r,]
#   if(length(unique(dat$YEAR))<2){ next }
#   dat_pa <- dat %>%
#     dplyr::select(YEAR, Species, incidence) %>%
#     complete(YEAR, Species, fill = list(incidence = 0)) %>%
#     #each row is a unique cell and Species, each column is a year
#     spread(YEAR, incidence)
#   # how many singletons (present in only 1 year)    #FIXME: This is a *really* ugly way to count
#   singletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==1))
#   # how many doubletons (present in only 2 years)
#   doubletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==2))
#   # how many incidences in the whole matrix
#   N <- length(which(dat_pa[,2:ncol(dat_pa)] == 1))
#   # how many samples (years) in the matrix
#   T <- ncol(dat_pa) - 1
#   # Chat for incidence: eqn Chat_sample(T) from Table 2, Chao et al. 2014, similar to eqn 4a in Chao & Jost Ecology
#   Chat_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2*doubletons))
#   # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
#   # correction for communities with no doubletons (this prevents NaN results for
#   # either singletons==0 or doubletons==0)
#   if (doubletons > 0) {
#     Chat_corrected_i <- Chat_i  }
#   else {
#     Chat_corrected_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2))  }
#   calcs <- data.frame(rarefyID=r, singletons, doubletons, N, Chat_i, Chat_corrected_i)
#   pa_coverage <- rbind(pa_coverage, calcs)
#   rm(dat_pa)
# }
#
# #what is the mean and sd of the presence data?
# mean(pa_coverage$Chat_corrected_i)
# sd(pa_coverage$Chat_corrected_i)
# #plot the presence coverage data using the abundance-based mean and standard deviation as a guide
# ggplot(pa_coverage, aes(Chat_corrected_i)) + geom_histogram(aes(fill=Chat_corrected_i>=mn-sd), binwidth=0.05) + geom_vline(xintercept=mn)
#
# #-------------------
# # BIOMASS DATA coverage (as INCIDENCE)... by rarefyID (year as samples)
# bt_grid_collate_biomass <- bt_grid %>%
#   # incidence data and remove NAs
#   filter(is.na(ABUNDANCE_TYPE)) %>%
#   filter(!is.na(Biomass)) %>%
#   #took out ObsEventID and PLOT as groups b/c no distinct plots for biomass data (Few years have > 1 sample)
#   group_by(CLIMATE, REALM, TAXA, StudyMethod, BIOMASS_TYPE, rarefyID, STUDY_ID, YEAR, cell, Species) %>%
#   summarise(
#     Biomass = sum(Biomass, na.rm=TRUE)) %>%
#   ungroup() %>%
#   dplyr::select(-CLIMATE, -REALM, -TAXA, -StudyMethod, -BIOMASS_TYPE) %>%
#   # create incidence column
#   mutate(incidence = ifelse(Biomass==0, 0, 1))
#
# bm_coverage <- data.frame()
# for (r in unique(bt_grid_collate_biomass$rarefyID)){
#   dat <- bt_grid_collate_biomass[bt_grid_collate_biomass$rarefyID == r,]
#   if(length(unique(dat$YEAR))<2){ next }
#   dat_pa <- dat %>%
#     dplyr::select(YEAR, Species, incidence) %>%
#     complete(YEAR, Species, fill = list(incidence = 0)) %>%
#     #each row is a unique cell and Species, each column is a year
#     spread(YEAR, incidence)
#   # how many singletons (present in only 1 year)    #FIXME: This is a *really* ugly way to count
#   singletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==1))
#   # how many doubletons (present in only 2 years)
#   doubletons <- length(which(rowSums(dat_pa[,2:ncol(dat_pa)])==2))
#   # how many incidences in the whole matrix
#   N <- length(which(dat_pa[,2:ncol(dat_pa)] == 1))
#   # how many samples (years) in the matrix
#   T <- ncol(dat_pa) - 1
#   # Chat for incidence: eqn Chat_sample(T) from Table 2, Chao et al. 2014, similar to eqn 4a in Chao & Jost Ecology
#   Chat_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2*doubletons))
#   # with fix from Chao et al 2014 Subroutine for iNEXT code (appendix)
#   # correction for communities with no doubletons (this prevents NaN results for
#   # either singletons==0 or doubletons==0)
#   if (doubletons > 0) {
#     Chat_corrected_i <- Chat_i  }
#   else {
#     Chat_corrected_i <- 1 - (singletons/N) * (((T-1)*singletons)/((T-1)*singletons + 2))  }
#   calcs <- data.frame(rarefyID=r, singletons, doubletons, N, Chat_i, Chat_corrected_i)
#   bm_coverage <- rbind(bm_coverage, calcs)
#   rm(dat_pa)
# }
#
# mean(bm_coverage$Chat_corrected_i)
# sd(bm_coverage$Chat_corrected_i)
# #plot the presence coverage data using the abundance-based mean and standard deviation as a guide
# ggplot(bm_coverage, aes(Chat_corrected_i)) + geom_histogram(aes(fill=Chat_corrected_i>=mn-sd), binwidth=0.05)
#
#
# ## Make a list of the rarefyIDs to keep for analysis based on coverage threshold (>= mn-sd of count data)
# ##  NOTE: Count data will drop individual years that don't meet the criteria, and Presence and Biomass data will drop entire rarefyIDs
# countkeep <- unique(abund_coverage[abund_coverage$Chat_corrected>=mn-sd,c('rarefyID', 'YEAR')])
# countkeep <- unite(countkeep, col=keep, rarefyID, YEAR, sep="_")
# countkeep <- as.vector(countkeep$keep)
#
# pakeep <- unique(pa_coverage[pa_coverage$Chat_corrected_i>=mn-sd,'rarefyID']) # presence
# bmkeep <- unique(bm_coverage[bm_coverage$Chat_corrected_i>=mn-sd,'rarefyID']) # biomass
#
#
# # Filter gridded studies for the coverage cutoff, prior to rarefaction
# bt_count_filtered <- bt_grid %>%
#   unite(col=keep, rarefyID, YEAR, sep="_", remove=FALSE) %>%
#   filter(keep %in% countkeep) %>%
#   dplyr::select(-keep)
#
# bt_pabm_filtered <- bt_grid %>%
#   filter(rarefyID %in% pakeep | rarefyID %in% bmkeep)
#
# # DATASET filtered for coverage (>= mean-sd of count data) to be used for rarefaction
# #   Note: count data was filtered at rarefyID + YEAR, and presence and biomass data was filtered at rarefyID
# bt_grid_filtered <- rbind(bt_count_filtered, bt_pabm_filtered)
#
# #save locally
# usethis::use_data(dgg)
# usethis::use_data(bt_grid_filtered)
