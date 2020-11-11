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
rarefied_metrics <- full_join(species_metrics, fd_metrics)

##	pull out new metadata
new_meta <- rarefied_metrics %>%
	distinct(rarefyID, SamplePool, SampleN, num_years, duration, startYear, endYear)


##	calculate the medians for all the metrics
rarefied_medians <- ungroup(rarefied_metrics) %>%
  select(-c(SamplePool, SampleN, num_years, duration, startYear, endYear, rarefy_resamp)) %>%
  group_by(rarefyID, YEAR, cell) %>%
  dplyr::summarise(across(.cols = everything(), median)) %>%
  ungroup()

rarefied_ints <- ungroup(rarefied_metrics) %>%
  select(-c(SamplePool, SampleN, num_years, duration, startYear, endYear, rarefy_resamp)) %>%
  group_by(rarefyID, YEAR, cell) %>%
  dplyr::summarise(N_int = round(median(N)),
                   S_int = round(median(S)),
                   ENSPIE_int = round(median(ENSPIE))) %>%
  ungroup()

##	recombine with new metadata
rarefied_medians <- inner_join(new_meta, rarefied_medians, by='rarefyID') %>%
  inner_join(rarefied_ints)

##	save
save(rarefied_medians, file=Sys.getenv('OFILE'))
