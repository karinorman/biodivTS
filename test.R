
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
  FD_mets <- get_FD(species_mat = rare_comm, trait_mat = traits, year_list = years,
                    data_id = rare_samp$rarefyID, samp_id = uniq_id,
                    w.abun = TRUE, m = "min", corr = "cailliez")
  J_func_components <- functional.beta.pair(x = rare_comm_binary, traits = FD_mets$pca_traits, index.family='jaccard')	# distance
  Jbeta_func <- as.matrix(J_func_components$funct.beta.jac)
  Jtu_func <- as.matrix(J_func_components$funct.beta.jtu)
  Jne_func <- as.matrix(J_func_components$funct.beta.jne)

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

  # combine univariate and turnover metrics
  biochange_metrics <- full_join(uni_metrics, simbaseline[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
    full_join(simnext[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
    full_join(simhind[-length(unique(rare_samp$YEAR)),], by=c('YEAR', 'cell')) %>%
    full_join(FD_mets$fd_met, by = c("YEAR" = "year"))

  # add to dataframe for all studies
  rarefied_metrics <- bind_rows(rarefied_metrics, biochange_metrics)
