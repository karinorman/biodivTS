
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

  synonyms <-

  #match all the un-ID'd names to other providers, check if synonyms given by those providers match known ITIS species
  alt_names <-map_df(providers[providers != original_provider], function(name) synonyms(unmatched, name)) %>%
    drop_na(synonym) %>%
    #new matches can be in either the synonym or acceptedNameUsage column, so let's combine them
    gather(type, altProviderName, synonym, acceptedNameUsage) %>%
    select(altProviderName, input) %>%
    distinct() %>%
    drop_na(altProviderName)

  #match accepted names from
  syn_res <- filter_name(na.omit(unique(alt_names$altProviderName)), original_provider) %>%
    drop_na(acceptedNameUsageID) %>%
    rename(altProviderName = input) %>%
    left_join(alt_names, na_matches = "never", by = "altProviderName")

  if (common == TRUE){
    common_providers <- c("col", "fb", "gbif", "itis", "iucn", "ncbi", "slb")
    unmatched <- unmatched[unmatched %in% syn_res$input]

    common_match <- map_df(common_providers[common_providers != original_provider],
           function(name) filter_common(unmatched, name)) %>%
      drop_na(acceptedNameUsageID)

    sci_match <- filter_name(unique(common_match$scientificName), "itis") %>%
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

#' synonyms
#'
#' Resolve provided list of names against all known synonyms
#' @inheritParams filter_name
#' @importFrom dplyr left_join
#' @export
#' @examples
#' \donttest{
#'   \dontshow{
#'    ## All examples use a temporary directory
#'    Sys.setenv(TAXADB_HOME=tempdir())
#'   }
#'
#' sp <- c("Trochalopteron henrici gucenense",
#'         "Trochalopteron elliotii")
#' synonyms(sp)
#'
#' }
#'
synonyms <- function(name,
                     provider = getOption("taxadb_default_provider", "itis"),
                     #version = latest_version(),
                     collect = TRUE,
                     db = td_connect()){


  the_id_table <- filter_name(name, provider = provider, db = db)

  ## Get both accepted names & synonyms for anything with an acceptedNameUsageID

    syn <-
      taxa_tbl(provider = provider, db = db) %>%
      safe_right_join(the_id_table %>%
                        select("acceptedNameUsageID"),
                      by = "acceptedNameUsageID",
                      copy = TRUE) %>%
      dplyr::select("scientificName", "acceptedNameUsageID",
                    "taxonomicStatus", "taxonRank") %>%
      syn_table()

  ## Join that back onto the id table
  out <- the_id_table %>%
    dplyr::select("scientificName", "sort", "acceptedNameUsageID", "input") %>%
    dplyr::left_join(syn, by = "acceptedNameUsageID", copy = TRUE) %>%
    dplyr::select("acceptedNameUsage", "synonym", "taxonRank",
                  "acceptedNameUsageID", "input") %>%
    dplyr::distinct()

  if (collect && inherits(out, "tbl_lazy")) {
    return( dplyr::collect(out) )
  }

  out

}


globalVariables(c("taxonomicStatus", "scientificName", "taxonID",
                  "taxonRank", "acceptedNameUsageID", "synonym", "input"))
## A mapping in which synonym and accepted names are listed in the same row

#' @importFrom dplyr full_join filter select
syn_table <- function(taxon, accepted = "accepted"){

  dplyr::full_join(
    taxon %>%
      dplyr::filter(taxonomicStatus != "accepted") %>%
      dplyr::select(synonym = scientificName,
                    acceptedNameUsageID),
    taxon %>%
      dplyr::filter(taxonomicStatus == "accepted") %>%
      dplyr::select(acceptedNameUsage = scientificName,
                    acceptedNameUsageID,
                    taxonRank),
    by = "acceptedNameUsageID")

}


## Manually copy query into DB, since RSQLite lacks right_join,
## and dplyr `copy` can only copy table "y"
#' @importFrom dbplyr remote_con
#' @importFrom DBI dbWriteTable
#' @importFrom dplyr left_join tbl
safe_right_join <- function(x, y, by = NULL, copy = FALSE, ...){

  if(copy){
    con <- dbplyr::remote_con(x)
    if(!is.null(con)){ ## only attempt on remote tables!
      tmpname <-  paste0(sample(letters, 10, replace = TRUE), collapse = "")
      DBI::dbWriteTable(con, tmpname, y, temporary = TRUE)
      y <- dplyr::tbl(con, tmpname)
    }
  }
  dplyr::left_join(y, x, by = by, ...)
}
