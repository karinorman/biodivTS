get_traitMat <- function(species_list, trait_data, scale = TRUE){

  trait_mat <- trait_data %>%
    filter(Species %in% species_list) %>%
    select_if(~ length(unique(na.omit(.))) > 1)  #remove columns that have the same value for all columns

  #get binary variables so they can be excluded from rescaling
  bin_vars <- map(trait_mat, ~ all(na.omit(.) %in% 0:1))

  if (isTRUE(scale)){
    trait_mat <- trait_mat %>%
      mutate_at(vars(-c(Species, names(bin_vars[bin_vars == TRUE])), -ends_with("5cat")),
                list(~as.numeric(scale(.)))) %>% #rescale variables, not binary or categorical
      arrange(Species) %>%
      column_to_rownames("Species")
  }else{
    trait_mat <- trait_mat %>%
      arrange(Species) %>%
      column_to_rownames("Species")
  }
  return(trait_mat)
}
