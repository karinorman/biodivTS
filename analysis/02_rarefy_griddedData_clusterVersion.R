##======================================================================
##	23 Feb 2017: code modified to run on iDiv eve cluster
##	returns all rarefied metrics for 1 resamples; to be combined and means
##	calculated in separate script
##======================================================================
# Primary Code Authors: Shane Blowes and Sarah Supp
# Email: sablowes@gmail.com, sarah@weecology.org
##======================================================================
rm(list=ls())

##	load packages
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(lazyeval)
library(vegan)
library(betapart)
library(furrr)

##==========================================
##	Get the gridded data locally
#load('data/BioTIME_grid_filtered.Rdata')
bt_grid_filtered <- pins::pin_get("bt-traitfiltered", board = "github") %>%
  rename_with(toupper) %>%
  rename(Species = SPECIES, StudyMethod = STUDYMETHOD, ObsEventID = OBSEVENTID,
         rarefyID = RAREFYID, cell = CELL, Abundance = ABUNDANCE, Biomass = BIOMASS)

##	use the job and task id as a counter for resampling and to set the seed
##	of the random number generator
uniq_id <- runif(1, min = 1, max = 2000)
seed <- uniq_id
#seed <- substr(seed, 1, 4)
set.seed(seed)

#=================================FUNCTION TO RAREFY DATA============================

rarefy_diversity <- function(grid, type=c("count", "presence", "biomass"), resamples=100, trimsamples=FALSE, ...){

  #	CALCULATE RAREFIED METRICS for each study for all years
  #	restrict calculations to where there is abundance>0 AND
  #	following removal of NAs and 0's there are still more than 2 years

  # Check if the data is count abundance: if yes, calculate all rarefied metrics
  # Check if the data is presence or biomass: if yes, calculate only S and Jaccards
  # If is.na(ABUNDANCE_TYPE), then should calculate on the Biomass column
  # grid: input dataset
  # type: count, presence, or biomass
  # resamples: the number of bootstrap resampling events desired (default is 100)
  # trimsamples: TRUE means that years with < 1/2 the average number of samples should be removed to avoid excessive information loss. Default is FALSE.
  # calculate the number of sampling events per year, and find the minimum
  # resample the data to rarefy the diversity metrics
  # output a new dataframe

  if(type == "count" | type == "presence") { field = "Abundance"
  } else { field = "Biomass" }

  # Get the sample-size to rarefy to. How many sampling events per cell per year?
  # This is handled differently depending on the data type


  # define a filter for 'field' to be >0 and not NA
  zero_NA_filter <- interp(~y > x & !is.na(y),
                           .values = list(y = as.name(field), x = 0))

  #	define a function to calculate the sum(field) for use on the rarefied sample
  sum_field <- interp(~sum(as.numeric(var), na.rm=T),
                      var= as.name(field))


  nsamples <- ungroup(grid) %>%
    group_by(rarefyID, YEAR) %>%
    # remove 0's or NA's (this is to catch any places where abundance wasn't actually recorded
    # or the species is indicated as absent)
    filter_(.dots=zero_NA_filter) %>%
    # calculate how many observations per year per study
    dplyr::summarise(nsamples = n_distinct(ObsEventID))

  # Check if you wanted to remove years with especially low samples (< 1/2 the average number of samples)
  if(trimsamples) {
    # Calculate the mean number of samples per cell
    mean_samp <- ungroup(nsamples) %>% group_by(rarefyID) %>%
      mutate(mean_samp = mean(nsamples),
             lower_bound = mean(nsamples)/2) %>%
      filter(nsamples >= lower_bound)

    min_samp <- ungroup(mean_samp) %>% group_by(rarefyID) %>%
      mutate(min_samp = min(nsamples)) %>%
      # retain only the rows with the minimum sample size for a given cell
      filter(nsamples==min_samp) %>%
      distinct(rarefyID, min_samp, .keep_all=TRUE)

    # join the data and filter out years with  nsamples less than 1/2 the mean number of samples
    grid <- inner_join(grid, dplyr::select(min_samp, - YEAR, -nsamples)) %>%
      filter(n_samps >= lower_bound)
    rm(mean_samp,min_samp)
    print("trimsamples==TRUE: Removing years with < 1/2 the average number of samples for a given ID")

  } else {
    # Calculate the minimum number of samples per cell
    min_samp <- ungroup(nsamples) %>% group_by(rarefyID) %>%
      mutate(min_samp = min(nsamples)) %>%
      # retain only the rows with the minimum sample size for a given cell
      filter(nsamples==min_samp) %>%
      distinct(rarefyID, min_samp, .keep_all=TRUE)

    #	Add the min_samp to the data and tidy a little
    grid <- inner_join(grid, dplyr::select(min_samp, - YEAR, -nsamples))
    rm(min_samp)
  }

  # Re-calculate metadata
  new_meta <- ungroup(grid) %>%
    # remove 0's or NA's (this is to catch any places where abundance wasn't actually recorded
    # or the species is indicated as absent)
    filter_(.dots=zero_NA_filter) %>%
    group_by(rarefyID) %>%
    summarise(
      # total number of species in rarefyID time-series
      SamplePool = n_distinct(Species),
      # total number of individuals
      SampleN = ifelse(type=='count', sum(as.numeric(Abundance)),
                       NA),
      # number of years sampled
      num_years = n_distinct(YEAR),
      # duration of time series, start and end points
      duration = max(YEAR) - min(YEAR) + 1,
      startYear = min(YEAR),
      endYear = max(YEAR))

  #	Create dataframe where unique observations (i.e., the data of an
  #	ObsEventID's [individual species abundances])
  #	are nested within cells within years within studies
  bt_grid_nest <- ungroup(grid) %>%
    group_by(rarefyID, ObsEventID, cell, YEAR, min_samp) %>%
    # remove 0's or NA's (to catch any places where abundance wasn't actually recorded or the species is indicated as absent)
    filter_(.dots=zero_NA_filter) %>%
    # depending on type: nest(Species, Abundance) OR nest(Species, Biomass)
    nest(data=c("Species", field)) %>%
    # reduce to studies that have more than two time points for a given cell
    group_by(rarefyID) %>%
    #keeps all study_cells with 2 or more years of data
    filter(n_distinct(YEAR)>=2) %>%
    ungroup()

  ##	initialise df to store all biochange metrics
  rarefied_metrics <- data.frame()

  ## create FD function that returns an empty dataframe instead of erroring if FD can't be calculated
  get_FD_safe <- possibly(get_FD, otherwise = data_frame())

  ##	rarefy rarefy_resamps times
  for(i in 1:resamples){

    ## loop to do rarefaction for each study
    for(j in 1:length(unique(bt_grid_nest$rarefyID))){
      print(paste('rarefaction', i, 'out of', resamples, 'for study_cell', j, '(', unique(bt_grid_nest$rarefyID)[j], ')',  'in', length(unique(bt_grid_nest$rarefyID))))

      ##	get the jth study_cell
      study <- bt_grid_nest %>%
        filter(rarefyID==unique(bt_grid_nest$rarefyID)[j])

      # get minimum sample size for rarefaction
      min_samp <- study %>% distinct(min_samp) %>% .$min_samp

      # check that there is only one cell represented (This shouldn't be a problem)
      if(length(unique(study$cell))>1) {
        stop(paste0("ERROR: ", unique(study$rarefyID), " contains more than one grid cell")) }

      # check there there is more than one year in the cell
      if(length(unique(study$YEAR))<2) {
        print(paste0("ERROR: ", unique(study$rarefyID), " does not have more than one year"))
        next }

      rare_samp <- study %>%
        # rarefy to min_samp
        group_by(rarefyID, YEAR) %>%
        sample_n(size=min_samp) %>%
        # unpack and collate taxa from rarefied sample
        unnest(data) %>%
        # add unique counter for a resampling event
        mutate(rarefy_resamp = uniq_id) %>%
        # collate species within cells within years
        group_by(rarefyID, YEAR, cell, rarefy_resamp, Species) %>%
        dplyr::summarise_(
          Abundance=sum_field) %>%  #
        ungroup()

      # create community matrix of rarefied sample
      rare_comm  <- ungroup(rare_samp) %>%
        arrange(Species) %>% #order species so the trait matrix can in the same order
        spread(Species, Abundance, fill=0)

      #get vector of years in the same order as the community matrix to tack back on to FD dataframe later
      years <- rare_comm$YEAR

      rare_comm <- rare_comm %>%
        column_to_rownames("YEAR") %>%
        dplyr::select(-rarefyID, -cell, -rarefy_resamp)

      # betapart requires presence/absence matrix for Jaccard calculations of turnover/nestedness
      rare_comm_binary <- with(rare_comm, ifelse(rare_comm > 0, 1, 0))

      # get trait matrix for the species in the sample
      traits <- get_traitMat(rare_samp$Species, trait_ref)

      # initialise matrices for calculating turnover
      simbaseline <- data.frame(array(NA, dim=c(length(unique(rare_samp$YEAR)), 14)))
      names(simbaseline)<-c('YEAR', 'cell', 'Jaccard_base','Horn_base','Chao_base','Pearson_base','Gains_base','Losses_base', 'Jbeta_base', 'Jtu_base', 'Jne_base', 'Jbeta_base_func', 'Jtu_base_func', 'Jne_base_func')

      simnext <- data.frame(array(NA, dim=c(length(unique(rare_samp$YEAR)), 14)))
      names(simnext)<-c('YEAR', 'cell', 'Jaccard_next','Horn_next','Chao_next','Pearson_next','Gains_next','Losses_next', 'Jbeta_next', 'Jtu_next', 'Jne_next', 'Jbeta_next_func', 'Jtu_next_func', 'Jne_next_func')

      simhind <- data.frame(array(NA, dim=c(length(unique(rare_samp$YEAR)), 14)))
      names(simhind)<-c('YEAR', 'cell', 'Jaccard_hind','Horn_hind','Chao_hind','Pearson_hind','Gains_hind','Losses_hind', 'Jbeta_hind', 'Jtu_hind', 'Jne_hind', 'Jbeta_hind_func', 'Jtu_hind_func', 'Jne_hind_func')

      counter2 <- 1


      if(type=="count"){

        # calculating between year similarities (NOT DISTANCE!) with Jaccard, Morisita-Horn, Chao and Pearson correlations
        Pearsoncor <- cor(t(log(rare_comm+1)), method='pearson')
        Jacsim <- as.matrix(1-vegdist(rare_comm, method='jaccard', binary=TRUE))
        Hornsim <- as.matrix(1-vegdist(rare_comm, method='horn'))
        Chaosim <- as.matrix(1-vegdist(rare_comm, method='chao'))
        Gainssim <- as.matrix(designdist(rare_comm, method = "B-J", terms = "binary",  abcd = FALSE, alphagamma = FALSE, "gains"))
        Lossessim <- as.matrix(designdist(rare_comm, method = "A-J", terms = "binary",  abcd = FALSE, alphagamma = FALSE, "losses"))
        # two steps for Jaccard components (so as calculation is done only once)
        J_components <- beta.pair(rare_comm_binary, index.family='jaccard')	# distance
        Jbeta <- as.matrix(J_components$beta.jac)
        Jtu <- as.matrix(J_components$beta.jtu)
        Jne <- as.matrix(J_components$beta.jne)
        n <- length(years)

        # calculated functional diversity
        FD_mets <- get_FD_safe(species_mat = rare_comm, trait_mat = traits, year_list = years,
                               data_id = rare_samp$rarefyID, samp_id = uniq_id,
                               w.abun = TRUE, m = "min", corr = "cailliez")
        if (is.null(dim(FD_mets))){
          J_func_components <- functional.beta.pair(x = rare_comm_binary, traits = FD_mets$pca_traits, index.family='jaccard')	# distance
          Jbeta_func <- as.matrix(J_func_components$funct.beta.jac)
          Jtu_func <- as.matrix(J_func_components$funct.beta.jtu)
          Jne_func <- as.matrix(J_func_components$funct.beta.jne)
        } else {
          size <- dim(rare_comm_binary)
          Jbeta_func <- matrix(nrow = size, ncol = size)
          Jtu_func <- matrix(nrow = size, ncol = size)
          Jne_func <- matrix(nrow = size, ncol = size)
        }

        # calculate univariate metrics
        uni_metrics <- ungroup(rare_samp) %>%
          group_by(rarefyID, YEAR, cell, rarefy_resamp) %>%
          summarise(
            N = sum(as.numeric(Abundance)),
            S = n_distinct(Species),
            PIE = diversity(Abundance, index='simpson'),
            ENSPIE = diversity(Abundance, index='invsimpson'),
            pielou = diversity(Abundance)/log(S),
            Hill1 = renyi(Abundance,scales = 1,hill = T)[1]) %>%
          #Hill2 = renyi(Abundance,scales = 2,hill = T)[1]
          ungroup()

      }
      else {
        # ONLY the metrics we need for biomass and presence (S, and Jaccard's)
        # For presence data, first convert rare_comm to a binary species matrix (0,1)   #FIXME: Also do this for Biomass? Because of limiting calcs, the only thing this would affect is Pearson cor?
        if(type=="presence" | type=='biomass'){
          rare_comm[rare_comm >0 ] <- 1 }

        # calculating between year similarities (NOT DISTANCE!) with Jaccard, Morisita-Horn, Chao and Pearson correlations
        n <- length(unique(rare_samp$YEAR))
        Pearsoncor <- cor(t(log(rare_comm+1)), method='pearson')
        Jacsim <- as.matrix(1-vegdist(rare_comm, method='jaccard', binary=TRUE))
        Hornsim <- matrix(nrow=n, ncol=n, NA)
        Chaosim <- matrix(nrow=n, ncol=n, NA)
        Gainssim <- matrix(nrow=n, ncol=n, NA)
        Lossessim <- matrix(nrow=n, ncol=n, NA)
        # two steps for Jaccard components (so as calculation is done only once)
        J_components <- beta.pair(rare_comm, index.family='jaccard')	# distance
        Jbeta <- as.matrix(J_components$beta.jac)
        Jtu <- as.matrix(J_components$beta.jtu)
        Jne <- as.matrix(J_components$beta.jne)

        # calculated functional diversity
        FD_mets <- get_FD_safe(species_mat = rare_comm, trait_mat = traits, year_list = years,
                               data_id = rare_samp$rarefyID, samp_id = uniq_id,
                               w.abun = FALSE, m = "min", corr = "cailliez")
        if (is.null(dim(FD_mets))){
          J_func_components <- functional.beta.pair(x = rare_comm_binary, traits = FD_mets$pca_traits, index.family='jaccard')	# distance
          Jbeta_func <- as.matrix(J_func_components$funct.beta.jac)
          Jtu_func <- as.matrix(J_func_components$funct.beta.jtu)
          Jne_func <- as.matrix(J_func_components$funct.beta.jne)
        } else{
          size <- dim(rare_comm_binary)
          Jbeta_func <- matrix(nrow = size, ncol = size)
          Jtu_func <- matrix(nrow = size, ncol = size)
          Jne_func <- matrix(nrow = size, ncol = size)
        }

        # calculate univariate metrics
        uni_metrics <- ungroup(rare_samp) %>%
          group_by(rarefyID, YEAR, cell, rarefy_resamp) %>%
          summarise(
            N = NA,
            S = n_distinct(Species),
            PIE = NA,
            ENSPIE = NA,
            pielou = NA,
            Hill1 = NA) %>%
          #Hill2 = NA
          ungroup()
      }

      #Collate metrics
      # How baseline is calculated.
      # NB: we do not compare first year with itself here (to be added before model fitting later)
      simbaseline[counter2:(counter2+n-2),] <- cbind(
        unique(rare_samp$YEAR)[2:n],
        unique(rare_samp$cell),
        Jacsim[2:n],
        Hornsim[2:n],
        Chaosim[2:n],
        Pearsoncor[2:n],
        Gainssim[2:n],
        Lossessim[2:n],
        Jbeta[2:n],
        Jtu[2:n],
        Jne[2:n],
        Jbeta_func[2:n],
        Jtu_func[2:n],
        Jne_func[2:n])

      # How consecutive is calculated.
      simnext[counter2:(counter2+n-2),] <- cbind(
        unique(rare_samp$YEAR)[2:n],
        unique(rare_samp$cell),
        Jacsim[row(Jacsim)-col(Jacsim)==1],
        Hornsim[row(Hornsim)-col(Hornsim)==1],
        Chaosim[row(Chaosim)-col(Chaosim)==1],
        Pearsoncor[row(Pearsoncor)-col(Pearsoncor)==1],
        Gainssim[row(Gainssim)-col(Gainssim)==1],
        Lossessim[row(Lossessim)-col(Lossessim)==1],
        Jbeta[row(Jbeta)-col(Jbeta)==1],
        Jtu[row(Jtu)-col(Jtu)==1],
        Jne[row(Jne)-col(Jne)==1],
        Jbeta_func[row(Jbeta_func)-col(Jbeta_func)==1],
        Jtu_func[row(Jtu_func)-col(Jtu_func)==1],
        Jne_func[row(Jne_func)-col(Jne_func)==1])

      # How hindcasting is calculated.
      simhind[counter2:(counter2+n-2),] <- cbind(
        unique(rare_samp$YEAR)[1:(n-1)],
        unique(rare_samp$cell),
        Jacsim[row(Jacsim)%in%1:(max(row(Jacsim))-1) & col(Jacsim)==max(col(Jacsim))],
        Hornsim[row(Hornsim)%in%1:(max(row(Hornsim))-1) & col(Hornsim)==max(col(Hornsim))],
        Chaosim[row(Chaosim)%in%1:(max(row(Chaosim))-1) & col(Chaosim)==max(col(Chaosim))],
        Pearsoncor[row(Pearsoncor)%in%1:(max(row(Pearsoncor))-1) & col(Pearsoncor)==max(col(Pearsoncor))],
        Gainssim[row(Gainssim)%in%1:(max(row(Gainssim))-1) & col(Gainssim)==max(col(Gainssim))],
        Lossessim[row(Lossessim)%in%1:(max(row(Lossessim))-1) & col(Lossessim)==max(col(Lossessim))],
        Jbeta[row(Jbeta)%in%1:(max(row(Jbeta))-1) & col(Jbeta)==max(col(Jbeta))],
        Jtu[row(Jtu)%in%1:(max(row(Jtu))-1) & col(Jtu)==max(col(Jtu))],
        Jne[row(Jne)%in%1:(max(row(Jne))-1) & col(Jne)==max(col(Jne))],
        Jbeta_func[row(Jbeta_func)%in%1:(max(row(Jbeta_func))-1) & col(Jbeta_func)==max(col(Jbeta_func))],
        Jtu_func[row(Jtu_func)%in%1:(max(row(Jtu_func))-1) & col(Jtu_func)==max(col(Jtu_func))],
        Jne_func[row(Jne_func)%in%1:(max(row(Jne_func))-1) & col(Jne_func)==max(col(Jne_func))])


      # combine univariate and turnover metrics
      biochange_metrics <- full_join(uni_metrics, simbaseline[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
        full_join(simnext[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
        full_join(simhind[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell'))

      if(is.null(dim(FD_mets))){
        biochange_metrics <- full_join(biochange_metrics, FD_mets$fd_met, by = c("YEAR" = "year"))}

      # add to dataframe for all studies
      rarefied_metrics <- bind_rows(rarefied_metrics, biochange_metrics)

    }	# rarefyID loop (STUDY_CELL ID)
  }	# rarefaction loop
  rarefied_metrics <- as_tibble(rarefied_metrics)
  # combine with the new metadata
  rarefied_metrics <- inner_join(new_meta, rarefied_metrics, by='rarefyID')
  return(rarefied_metrics)
} # END function

##================================== Calculate mean rarefied diversity for each data type=======================================
get_samp_metrics <- function(data, samp_number){
  library(biodivTS)
  ## Separate true abundance (count) data from presence and biomass/cover data
  bt_grid_abund <- data %>%
    filter(ABUNDANCE_TYPE %in% c("Count", "Density", "MeanCount"))

  bt_grid_pres <- data %>%
    filter(ABUNDANCE_TYPE == "Presence/Absence")

  #We don't have any data with only biomass measurements, so don't need this
  # bt_grid_bmass <- data %>%
  #   filter(is.na(ABUNDANCE_TYPE)) #only want to calculate on biomass data when abundance data is also not available

  #-----------------------------Alternate for trimming years with low samples
  # rarefy diversity after trimming years with excessively low samples (sensu Dornelas et al. 2014, methods, Fig S10)

  ## Get rarefied sample for each measurement type (trimmed data)
  #rarefy_abund_trim <- rarefy_diversity(grid=bt_grid_abund, type="count", resamples=1, trimsamples=TRUE)
  #rarefy_pres_trim <- rarefy_diversity(grid=bt_grid_pres, type="presence", resamples=1, trimsamples=TRUE)
  #rarefy_biomass_trim <- rarefy_diversity(grid=bt_grid_bmass, type="biomass", resamples=1, trimsamples=TRUE)
  #-----------------------------
  ## Get rarefied resample for each type of measurement (all data)
  rarefy_abund <- rarefy_diversity(grid=bt_grid_abund, type="count", resamples=1)
  rarefy_pres <- rarefy_diversity(grid=bt_grid_pres, type="presence", resamples=1)
  #rarefy_biomass <- rarefy_diversity(grid=bt_grid_bmass, type="biomass", resamples=1)

  ##	save rarefied sample (for later collation and calculation of mean)
  #save(rarefy_abund_trim,rarefy_pres_trim, rarefy_biomass_trim,
  dir <- here::here("data", "rarefied_metrics")
  dir.create(dir)
  save(rarefy_abund, rarefy_pres, #rarefy_biomass,
       file=paste0(dir, "/resample", samp_number, ".rda"))
}

##### Set Up for rarefying via parallel mapping ############
# furrr can't detect functions loaded via devtools::load_all(), so we have to source them manually
source(here::here("R/get_FD.R"))
source(here::here("R/dbFD_mine.R"))
source(here::here("R/get_traitMat.R"))

# create a file for writing out the rarefyID's that error out for the FD calculation

plan("multiprocess")
future_map(1:200, get_samp_metrics, data = bt_grid_filtered)
