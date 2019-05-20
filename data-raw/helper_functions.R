
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
    distinct(acceptedNameUsageID, input, !!orig_id) %>%
    group_by(!!orig_id) %>%
    filter(n()>1) %>%
    ungroup()
}
