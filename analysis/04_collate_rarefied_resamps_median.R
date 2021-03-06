##============================================================
##	script to combine rarefied resamples and calculate mean
##	for analysis
library(dplyr)
library(tibble)
library(stringr)
##============================================================

##	load species based metrics
path <- here::here("data", "rarefied_metrics")
filelist <-  dir(path, "*.rda") %>%
  paste0(path, "/", .)

species_metrics <- purrr::map_dfr(filelist, function(x) {load(x)
  return(rarefied_metrics)
	}) %>% distinct()

## load fd based metrics
path <- here::here("data", "rarefied_metrics", "fd")
filelist <-  dir(path, "*.rda") %>%
  paste0(path, "/", .)

fd_metrics <- purrr::map_dfr(filelist, function(x) {
    load(x)
    rarefyID <- str_match(x,"_cell\\s*(.*?)\\s*_sample")[,2]

    biochange_metrics <- biochange_metrics %>%
      mutate(rarefyID = rarefyID)
    return(biochange_metrics)
  }) %>% distinct()

## get null models
null_table <- pins::pin_get("null-table", board = "github")

### Merge rarefied species and fd metrics and null models
##	put them all together
rarefied_metrics <- full_join(species_metrics, fd_metrics,
                              c("YEAR", "cell", "rarefy_resamp", "rarefyID")) %>%
  filter(S != 1) %>% # remove samples that had one species as they result in strange metric values
  select(-qual.FRic)


##	pull out new metadata
new_meta <- rarefied_metrics %>%
	distinct(rarefyID, SamplePool, SampleN, num_years, duration, startYear, endYear)

##	calculate the medians for all the metrics for studies that were rarefied
metric_cols <- ungroup(rarefied_metrics) %>%
  select(-c(SamplePool, SampleN, num_years, duration, startYear, endYear, CWM, type)) %>%
  pivot_longer(-c(rarefyID, rarefied, YEAR, cell, rarefy_resamp), names_to = "metric") %>%
  left_join(null_table %>% select(-c(richness, se, lowerCI, upperCI)),
            by = c("YEAR" = "year", "rarefyID" = "rarefyid", "rarefy_resamp" = "rarefy_resamp", "metric" = "metric"))

# get a dataframe of rarefyID's and rarefication samples that we don't have null_model for, and a count of the number of resamples that are missing
no_null <- metric_cols %>%
  filter(metric %in% c("FRic", "FEve", "FDiv", "FDis"), is.na(mean)) %>%
  select(rarefyID, rarefy_resamp) %>%
  distinct()

missing_nulls <- null_table %>%
  right_join(no_null, by = c("rarefyid" = "rarefyID","rarefy_resamp")) %>%
  group_by(rarefyid) %>%
  summarise(n_missing_nulls = n_distinct(rarefy_resamp)) %>%
  mutate(commplete_null_samps = "FALSE")

# calculate SES and incorporate into metric columns
metric_cols <- metric_cols %>%
  mutate(SES = (value - mean)/sd) %>%
  #for years where species richness = richness of the species pool, FRic sd = 0, so use the original value with no SES correction
  mutate(SES = case_when(
    sd == 0 ~ 0,
    TRUE ~ SES)) %>%
  select(-c(mean, sd)) %>%
  #the next three lines go back to wide format for both SES and unadjusted values, then to long format incorporating SES values
  pivot_wider(names_from = "metric", values_from = c("value", "SES")) %>%
  rename_at(vars(starts_with("value")), ~str_replace(., "value_", "")) %>%
  pivot_longer(-c(rarefyID, rarefied, YEAR, cell, rarefy_resamp), names_to = "metric") %>%
  filter(!is.na(value))

# rarefied_medians <- metric_cols %>%
#   filter(rarefied == TRUE) %>%
#   group_by(rarefyID, YEAR, cell, rarefied, type) %>%
#   dplyr::summarise(across(.cols = everything(), ~median(.x, na.rm = TRUE))) %>%
#   ungroup() %>%
#   bind_rows(metric_cols %>% filter(rarefied == FALSE))

rarefied_medians <- metric_cols %>%
  group_by(rarefyID, YEAR, cell, rarefied, metric) %>%
  summarise(value = median(value, na.rm = TRUE))


# rarefied_ints <- ungroup(rarefied_metrics) %>%
#   select(-c(SamplePool, SampleN, num_years, duration, startYear, endYear, rarefy_resamp)) %>%
#   group_by(rarefyID, YEAR, cell) %>%
#   dplyr::summarise(N_int = round(median(N)),
#                    S_int = round(median(S)),
#                    ENSPIE_int = round(median(ENSPIE))) %>%
#   ungroup()

##	recombine with new metadata
rarefied_metrics <- inner_join(new_meta, rarefied_medians, by='rarefyID') %>%
  left_join(missing_nulls, by = c("rarefyID" = "rarefyid"))
  #inner_join(rarefied_ints)

##	save
usethis::use_data(rarefied_metrics)

pins::board_register_github(name = "github", repo = "karinorman/biodivTS_data", branch = "master")
pins::pin(rarefied_metrics, board = "github")
