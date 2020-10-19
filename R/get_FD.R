get_FD <- function(species_mat, trait_mat, year_list, data_id, samp_id, ...){
  #the convexhull c code creates an output file of vertices, in order to not have multiple threads writing
  #to the same file we have to create individual working directories
  tmp_dir <- paste0(here::here("tmp"), "/tmp_", data_id, "_", samp_id)
  fs::dir_create(tmp_dir)
  setwd(tmp_dir)

  fd_data <- dbFD_mine(trait_mat, species_mat, w.abun = TRUE, ...)
  fd_out <- as_tibble(fd_data$fd_met)[,1:8] %>%
    cbind(year = year_list, .)

  setwd(here::here())

  return(list(fd_met = fd_out, pca_traits = fd_data$pca_traits))
}
