##======================================================================
##	23 Feb 2017: code modified to run on iDiv eve cluster
##	returns all rarefied metrics for 1 resamples; to be combined and means
##	calculated in separate script
##======================================================================
# Primary Code Authors: Shane Blowes and Sarah Supp
# Email: sablowes@gmail.com, sarah@weecology.org
##======================================================================

##	load packages
library(dplyr)
library(purrr)
library(tibble)
library(tidyr)
library(lazyeval)
library(vegan)
library(betapart)
library(doParallel)
library(foreach)
library(FD)
devtools::load_all()

##==========================================
pins::board_register_github(repo = "karinorman/biodivTS", branch = "master")

##	Get the gridded data locally
#load('data/BioTIME_grid_filtered.Rdata')
bt_grid_filtered <- pins::pin_get("bt-traitfiltered", board = "github") %>%
  rename_with(toupper) %>%
  rename(Species = SPECIES, StudyMethod = STUDYMETHOD, ObsEventID = OBSEVENTID,
         rarefyID = RAREFYID, cell = CELL, Abundance = ABUNDANCE, Biomass = BIOMASS)

trait_ref <- pins::pin_get("trait-ref", board = "github")

#=================================FUNCTION TO RAREFY DATA============================

rarefy_diversity <- function(grid, type=c("count", "presence", "biomass"), resamples=100,
                             parallel, core_num, ...){

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

  # Calculate the minimum number of samples per cell
  min_samp <- ungroup(nsamples) %>% group_by(rarefyID) %>%
    mutate(min_samp = min(nsamples), max_samp = max(nsamples)) %>%
    #check for studies that have consistent sampling and don't need to be rarefied
    mutate(rarefy = if_else(min_samp != max_samp, TRUE, FALSE)) %>%
    # retain only the rows with the minimum sample size for a given cell
    select(-max_samp, -nsamples, -YEAR) %>%
    distinct()

  #	Add the min_samp to the data and tidy a little
  grid <- inner_join(grid, min_samp)
  rm(min_samp)

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
      endYear = max(YEAR),
      rarefied = unique(rarefy))

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

  ## loop to do rarefaction for each study
  for(i in 1:length(unique(bt_grid_nest$rarefyID))){
    #for(j in 1:3){

    #get ID for the study and filter data
    filter_id <- unique(bt_grid_nest$rarefyID)[i]

    #if(filter_id %in% filter(new_meta, rarefied == TRUE)$rarefyID){browser()}

    study <- bt_grid_nest %>%
      filter(rarefyID==filter_id)

    # get minimum sample size for rarefaction
    min_samp <- study %>% distinct(min_samp) %>% .$min_samp

    # check that there is only one cell represented (This shouldn't be a problem)
    if(length(unique(study$cell))>1) {
      stop(paste0("ERROR: ", filter_id, " contains more than one grid cell")) }

    # check there there is more than one year in the cell
    if(length(unique(study$YEAR))<2) {
      print(paste0("ERROR: ", filter_id, " does not have more than one year"))
      next }

    #check if we need to rarify
    if(isFALSE(unique(study$rarefy))){
      study_resamples <- 1
    }else{study_resamples <- resamples}

    ##	initialise df to store all biochange metrics
    rarefied_metrics <- data.frame()

    #set up parallelization
    # if(isTRUE(parallel) & study_resamples > 1){
    #   # set up parallel loop
    #   cl <- parallel::makeForkCluster(core_num) #not to overload your computer
    #   registerDoParallel(cl)
    # }else{
    #   registerDoSEQ()
    # }

    ##	rarefy rarefy_resamps times
    #rarefied_metrics <- foreach(j = 1:study_resamples, .combine=rbind) %dopar% {
    for(j in 1:study_resamples){
      print(paste('rarefaction', j, 'out of', study_resamples, 'for study_cell', i, '(', filter_id, ')',  'in', length(unique(bt_grid_nest$rarefyID))))

      rare_samp <- study %>%
        # rarefy to min_samp
        group_by(rarefyID, YEAR) %>%
        sample_n(size=min_samp) %>%
        # unpack and collate taxa from rarefied sample
        unnest(data) %>%
        # add unique counter for a resampling event
        mutate(rarefy_resamp = j) %>%
        # collate species within cells within years
        group_by(rarefyID, YEAR, cell, rarefy_resamp, Species) %>%
        dplyr::summarise_(
          Abundance=sum_field) %>%  #
        ungroup()

      # create community matrix of rarefied sample
      rare_comm  <- ungroup(rare_samp) %>%
        arrange(Species) %>% #order species so the trait matrix can in the same order
        spread(Species, Abundance, fill=0)

      #save out sample species matrix
      rare_comm_save <- rare_comm %>% mutate(type = type)
      path <- here::here("data", "rarefied_samples")
      dir.create(path)
      save(rare_comm_save, file=paste0(path, "/", type, "_cell", unique(study$rarefyID), "_sample", j, ".rda"))


      #get vector of years in the same order as the community matrix to tack back on to FD dataframe later
      years <- rare_comm$YEAR

      rare_comm <- rare_comm %>%
        column_to_rownames("YEAR") %>%
        dplyr::select(-rarefyID, -cell, -rarefy_resamp)

      # betapart requires presence/absence matrix for Jaccard calculations of turnover/nestedness
      rare_comm_binary <- with(rare_comm, ifelse(rare_comm > 0, 1, 0))

      # initialise matrices for calculating turnover
      simbaseline <- data.frame(array(NA, dim=c(length(unique(rare_samp$YEAR)), 11)))
      names(simbaseline)<-c('YEAR', 'cell', 'Jaccard_base','Horn_base','Chao_base','Pearson_base','Gains_base','Losses_base', 'Jbeta_base', 'Jtu_base', 'Jne_base')#, 'Jbeta_base_func', 'Jtu_base_func', 'Jne_base_func')

      simnext <- data.frame(array(NA, dim=c(length(unique(rare_samp$YEAR)), 11)))
      names(simnext)<-c('YEAR', 'cell', 'Jaccard_next','Horn_next','Chao_next','Pearson_next','Gains_next','Losses_next', 'Jbeta_next', 'Jtu_next', 'Jne_next')#, 'Jbeta_next_func', 'Jtu_next_func', 'Jne_next_func')

      simhind <- data.frame(array(NA, dim=c(length(unique(rare_samp$YEAR)), 11)))
      names(simhind)<-c('YEAR', 'cell', 'Jaccard_hind','Horn_hind','Chao_hind','Pearson_hind','Gains_hind','Losses_hind', 'Jbeta_hind', 'Jtu_hind', 'Jne_hind')#, 'Jbeta_hind_func', 'Jtu_hind_func', 'Jne_hind_func')

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
        Jne[2:n]#,
        # Jbeta_func[2:n],
        # Jtu_func[2:n],
        # Jne_func[2:n]
      )

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
        Jne[row(Jne)-col(Jne)==1]#,
        # Jbeta_func[row(Jbeta_func)-col(Jbeta_func)==1],
        # Jtu_func[row(Jtu_func)-col(Jtu_func)==1],
        # Jne_func[row(Jne_func)-col(Jne_func)==1]
      )

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
        Jne[row(Jne)%in%1:(max(row(Jne))-1) & col(Jne)==max(col(Jne))]#,
        # Jbeta_func[row(Jbeta_func)%in%1:(max(row(Jbeta_func))-1) & col(Jbeta_func)==max(col(Jbeta_func))],
        # Jtu_func[row(Jtu_func)%in%1:(max(row(Jtu_func))-1) & col(Jtu_func)==max(col(Jtu_func))],
        # Jne_func[row(Jne_func)%in%1:(max(row(Jne_func))-1) & col(Jne_func)==max(col(Jne_func))]
      )


      # combine univariate and turnover metrics
      biochange_metrics <- full_join(uni_metrics, simbaseline[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
        full_join(simnext[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
        full_join(simhind[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell'))

      # add to dataframe for all studies
      rarefied_metrics <- bind_rows(rarefied_metrics, biochange_metrics)


    }	# rarefyID loop (STUDY_CELL ID)

    #if(exists("cl")){stopCluster(cl)}

    rarefied_metrics <- inner_join(new_meta, rarefied_metrics) %>%
      mutate(type = type)

    #save out the files
    path <- here::here("data", "rarefied_metrics")
    dir.create(path)
    save(rarefied_metrics, file=paste0(path, "/rarefyID_", filter_id, type, ".rda"))

  }	# rarefaction loop
} # END function

##================================== Calculate mean rarefied diversity for each data type=======================================
## Separate true abundance (count) data from presence and biomass/cover data
bt_grid_abund <- bt_grid_filtered %>%
  filter(ABUNDANCE_TYPE %in% c("Count", "Density", "MeanCount"))

bt_grid_pres <- bt_grid_filtered %>%
  filter(ABUNDANCE_TYPE == "Presence/Absence")

## Get rarefied resample for each type of measurement (all data)
rarefy_abund <- rarefy_diversity(grid=bt_grid_abund, trait_data = trait_ref, type="count", trait_axes = "min", parallel = TRUE, resamples=200, core_num = 44)
rarefy_pres <- rarefy_diversity(grid=bt_grid_pres, trait_data = trait_ref, type="presence", trait_axes = "min", parallel = TRUE, resamples=200, core_num = 20)

calc_FD_mets <- function(file_path, traits){
  rare_comm_save <- loadRData(file_path)

  years <- rare_comm_save$YEAR

  species_mat <- rare_comm_save %>%
    column_to_rownames("YEAR") %>%
    dplyr::select(-c(rarefyID, cell, type, rarefy_resamp))

  # get trait matrix for the species in the sample
  traits <- get_traitMat(colnames(species_mat), trait_data = traits)

  # create FD function that returns an empty dataframe instead of erroring if FD can't be calculated
  get_FD_safe <- possibly(get_FD, otherwise = data_frame())
  betapart_safe <- possibly(functional.beta.pair, otherwise = data.frame())

  # calculated functional diversity
  if(unique(rare_comm_save$type) == "count"){
    FD_mets <- get_FD_safe(species_mat = species_mat, trait_mat = traits, year_list = years,
                           data_id = unique(rare_comm_save$rarefyID), samp_id = unique(rare_comm_save$rarefy_resamp),
                           w.abun = TRUE, m = "min", corr = "cailliez")
  } else{
    FD_mets <- get_FD_safe(species_mat = species_mat, trait_mat = traits, year_list = years,
                           data_id = unique(rare_comm_save$rarefyID), samp_id = unique(rare_comm_save$rarefy_resamp),
                           w.abun = FALSE, m = "min", corr = "cailliez")
  }

  #check if FD metrics could be calculated, otherwise return empty dataframe
  if (is.null(dim(FD_mets))){
    species_mat_binary <- with(species_mat, ifelse(species_mat > 0, 1, 0))
    J_func_components <- betapart_safe(x = species_mat_binary, traits = FD_mets$pca_traits, index.family='jaccard')	# distance
    #check if the betapart metrics could be calculated, otherwise return dataframe of just FD metrics
    if(is.null(dim(J_func_components))){
      Jbeta_func <- as.matrix(J_func_components$funct.beta.jac)
      Jtu_func <- as.matrix(J_func_components$funct.beta.jtu)
      Jne_func <- as.matrix(J_func_components$funct.beta.jne)

      n <- length(unique(rare_comm_save$YEAR))

      # initialise matrices for calculating turnover
      simbaseline <- data.frame(array(NA, dim=c(n, 5)))
      names(simbaseline)<-c('YEAR', 'cell', 'Jbeta_base_func', 'Jtu_base_func', 'Jne_base_func')

      simnext <- data.frame(array(NA, dim=c(n, 5)))
      names(simnext)<-c('YEAR', 'cell', 'Jbeta_next_func', 'Jtu_next_func', 'Jne_next_func')

      simhind <- data.frame(array(NA, dim=c(n, 5)))
      names(simhind)<-c('YEAR', 'cell', 'Jbeta_hind_func', 'Jtu_hind_func', 'Jne_hind_func')

      counter2 <- 1

      #Collate metrics
      # How baseline is calculated.
      # NB: we do not compare first year with itself here (to be added before model fitting later)
      simbaseline[counter2:(counter2+n-2),] <- cbind(
        unique(rare_comm_save$YEAR)[2:n],
        unique(rare_comm_save$cell),
        Jbeta_func[2:n],
        Jtu_func[2:n],
        Jne_func[2:n]
      )

      # How consecutive is calculated.
      simnext[counter2:(counter2+n-2),] <- cbind(
        unique(rare_comm_save$YEAR)[2:n],
        unique(rare_comm_save$cell),
        Jbeta_func[row(Jbeta_func)-col(Jbeta_func)==1],
        Jtu_func[row(Jtu_func)-col(Jtu_func)==1],
        Jne_func[row(Jne_func)-col(Jne_func)==1]
      )

      # How hindcasting is calculated.
      simhind[counter2:(counter2+n-2),] <- cbind(
        unique(rare_comm_save$YEAR)[1:(n-1)],
        unique(rare_comm_save$cell),
        Jbeta_func[row(Jbeta_func)%in%1:(max(row(Jbeta_func))-1) & col(Jbeta_func)==max(col(Jbeta_func))],
        Jtu_func[row(Jtu_func)%in%1:(max(row(Jtu_func))-1) & col(Jtu_func)==max(col(Jtu_func))],
        Jne_func[row(Jne_func)%in%1:(max(row(Jne_func))-1) & col(Jne_func)==max(col(Jne_func))]
      )

      # combine univariate and turnover metrics
      biochange_metrics <- full_join(simbaseline[-n,],
                                     FD_mets$fd_met %>% dplyr::select(-c(nbsp, sing.sp)),
                                     by=c('YEAR' = 'year')) %>%
        mutate(cell = unique(na.omit(cell))) %>%
        full_join(simnext[-n,], by=c('YEAR', 'cell')) %>%
        full_join(simhind[-n,], by=c('YEAR', 'cell')) %>%
        mutate(rarefy_resamp = unique(rare_comm_save$rarefy_resamp))

    } else{
      biochange_metrics <- FD_mets$fd_met %>%
        select(-c(nbsp, sing.sp))
    }
  }else{ biochange_metrics <- tibble() }

  path <- here::here("data", "rarefied_metrics", "fd")
  dir.create(path)
  save(biochange_metrics, file=paste0(path, "/", unique(rare_comm_save$type),
                                      "_cell", unique(rare_comm_save$rarefyID),
                                      "_sample", unique(rare_comm_save$rarefy_resamp), ".rda"))
}

path <- here::here("data", "rarefied_samples")
files <- dir(path, "*.rda") %>%
  paste0(path, "/", .)

plan("multiprocess", workers = 4)
furrr::future_map(files[4968:10473], calc_FD_mets, traits = trait_ref)

