##============================================================
##	script to combine rarefied resamples and calculate mean
##	for analysis
library(dplyr)
library(tibble)
##============================================================
##	set the pattern to load the files to be compiled
path <- here::here("data", "rarefied_metrics")
filelist <-  dir(path, "*.rda") %>%
  paste0(path, "/", .)

species_metrics <- purrr::map_dfr(filelist, function(x) {load(x)
  return(rarefied_metrics)
	})

path <- here::here("data", "rarefied_metrics", "fd")
filelist <-  dir(path, "*.rda") %>%
  paste0(path, "/", .)

fd_metrics <- purrr::map_dfr(filelist, function(x) {load(x)
  return(biochange_metrics)
})

##	put them all together
rarefied_metrics <-

##	pull out new metadata
new_meta <- rarefied_metrics %>%
	distinct(rarefyID, SamplePool, SampleN, num_years, duration, startYear, endYear)


##	calculate the medians for all the metrics
rarefied_medians <- ungroup(rarefied_metrics) %>%
  group_by(rarefyID, YEAR, cell) %>%
  dplyr::summarise(
    N = median(N),
    N_int = round(median(N)),
    S = median(S),	# as above?
    S_int = round(median(S)),
    PIE = median(PIE),
    ENSPIE = median(ENSPIE),
    ENSPIE_int = round(median(ENSPIE)),
    pielou = median(pielou),
    Hill1 = median(Hill1),
    #MSA = median(MSA),
    Jaccard_base = median(Jaccard_base),
    Horn_base = median(Horn_base),
    Chao_base = median(Chao_base),
    Pearson_base = median(Pearson_base),
    Gains_base = median(Gains_base),
    Losses_base = median(Losses_base),
    Jbeta_base = median(Jbeta_base),
    Jtu_base = median(Jtu_base),
    Jne_base = median(Jne_base),
    Jbeta_base_func = median(Jbeta_base),
    Jtu_base_func = median(Jtu_base),
    Jne_base_func = median(Jne_base),
    Jaccard_next = median(Jaccard_next),
    Horn_next = median(Horn_next),
    Chao_next = median(Chao_next),
    Pearson_next = median(Pearson_next),
    Gains_next = median(Gains_next),
    Losses_next = median(Losses_next),
    Jbeta_next = median(Jbeta_next),
    Jtu_next = median(Jtu_next),
    Jne_next = median(Jne_next),
    Jbeta_next_func = median(Jbeta_next),
    Jtu_next_func = median(Jtu_next),
    Jne_next_func = median(Jne_next),
    Jaccard_hind = median(Jaccard_hind),
    Horn_hind = median(Horn_hind),
    Chao_hind = median(Chao_hind),
    Pearson_hind = median(Pearson_hind),
    Gains_hind = median(Gains_hind),
    Losses_hind = median(Losses_hind),
    Jbeta_hind = median(Jbeta_hind),
    Jtu_hind = median(Jtu_hind),
    Jne_hind = median(Jne_hind),
    Jbeta_hind_func = median(Jbeta_hind),
    Jtu_hind_func = median(Jtu_hind),
    Jne_hind_func = median(Jne_hind),
    FRic = median(FRic),
    FEve = median(FEve),
    FDiv = median(FDiv),
    FDis = median(FDis),
    RaoQ = median(RaoQ)) %>%
  ungroup()
##	recombine with new metadata
rarefied_medians <- inner_join(new_meta, rarefied_medians, by='rarefyID')

##	save
save(rarefied_medians, file=Sys.getenv('OFILE'))
