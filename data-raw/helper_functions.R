
get_match_counts <- function(data, species_col){
  species_col <- lazyeval::as_name(species_col)
  data %>%
    select(!!species_col) %>%
    drop_na() %>%
    distinct() %>%
    mutate(
      gbif = get_ids(!!species_col, "gbif"),
      col = get_ids(!!species_col, "col" ),
      itis = get_ids(!!species_col, "itis"),
      ncbi = get_ids(!!species_col, "ncbi"),
      wd = get_ids(!!species_col, "wd"  ),
      iucn = get_ids(!!species_col, "iucn"),
      ott = get_ids(!!species_col, "ott" )
    ) %>%
    #select(!!call("-", species_col)) %>%
    purrr::map_dbl(function(x) sum(!is.na(x)))
}


get_dupe_ids <- function(data, orig_id){
  orig_id <- lazyeval::as_name(orig_id)

  data %>%
    distinct(acceptedNameUsageID, !!orig_id) %>%
    group_by(!!orig_id) %>%
    filter(n()>1) %>%
    ungroup()
}

#Check if synonyms from other providers have matches in the original database
match_providers <- function(data, original_provider, common = FALSE) {
  unmatched <- data %>% filter(is.na(id)) %>% pull(sourceName) %>% unique()
  providers <- c("itis", "ncbi", "col", "gbif", "fb", "wd", "ott", "iucn")

  #match all the un-ID'd names to other providers, check if synonyms given by those providers match known ITIS species
  alt_names <-map_df(providers[providers != original_provider], function(name) taxadb::synonyms(unmatched, name)) %>%
    drop_na(acceptedNameUsageID) %>%
    #new matches can be in either the synonym or acceptedNameUsage column, so let's combine them
    gather(type, altProviderName, synonym, acceptedNameUsage) %>%
    select(altProviderName, input) %>%
    distinct() %>%
    drop_na(altProviderName)

  #match accepted names from
  syn_res <- by_name(na.omit(unique(alt_names$altProviderName)), original_provider) %>%
    drop_na(acceptedNameUsageID) %>%
    rename(altProviderName = input) %>%
    left_join(alt_names, na_matches = "never", by = "altProviderName")

  if (common == TRUE){
    common_providers <- c("col", "fb", "gbif", "itis", "iucn", "ncbi", "slb")
    unmatched <- unmatched[unmatched %in% syn_res$input]

    common_match <- map_df(common_providers[common_providers != original_provider],
           function(name) by_common(unmatched, name)) %>%
      drop_na(acceptedNameUsageID)

    sci_match <- by_name(unique(common_match$scientificName), "itis") %>%
      rename(altProviderName = input) %>%
      drop_na(acceptedNameUsageID) %>%
      left_join(common_match %>% select(scientificName, input),
                by = c("altProviderName" = "scientificName"))

    return(bind_rows(syn_res, sci_match))
  } else {return(syn_res)}
}

#Identify ids that have matched to multiple entries (for trait databases), and average within ID so there is only one entry
undupe_ids <- function(data) {
  #get ids that have more than one entry
  dupe_ids <- data %>%
    select(id, sourceName) %>%
    drop_na(id) %>%
    distinct() %>%
    group_by(id) %>%
    filter(n() > 1) %>%
    pull(id)

  #get data associated with ids
  dupe_data <- data %>%
    filter(id %in% dupe_ids)

  #resolve so that traits are averaged for species that need to be combined
  resolved <- dupe_data %>%
    select(-scientificName,-sourceName) %>%
    group_by(id) %>%
    summarise_all( ~ if (all(is.na(.))) {NA} else {mean(., na.rm = TRUE)}) %>%
    mutate(scientificName = get_names(id))

  #remove unresolved data and replace with new resolved
  resolved_data <- data %>%
    anti_join(dupe_data) %>%
    bind_rows(resolved)
}
