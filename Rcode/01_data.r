# clean body mass data ----

if(! dir.exists("data_output/cleaned_bodymass_data/")){
  generate_cleaned_CSV_files(DataDir = "data_raw/bodymass_data/",
                             species_path = "data_raw/species_table.csv", 
                             output_dir = "data_output/cleaned_bodymass_data/", 
                             FromRecord = 1, ToRecord = 198)
}

species = read_csv("data_raw/ReducedListSpecies.csv") # 175 sp

# Elton traits ----
if(file.exists("data_output/Elton_axes_bird.csv") &
   file.exists("data_output/Elton_axes_mammal.csv")){
  elton_axes_bird = read_csv("data_output/Elton_axes_bird.csv")
  elton_axes_mammal = read_csv("data_output/Elton_axes_mammal.csv")
} else {
  # modified from Ben Baiser's code; thus omit data explore 
  bird_traits = read_csv("data_raw/BirdFuncDat.csv") # bird traits from Elton Database
  mammal_traits = read_csv("data_raw/MamFuncDat.csv") # mammal traits from Elton Database
  
  # see if any species are not in Elton Data base
  not.in.all = setdiff(species$Species, c(bird_traits$Scientific, mammal_traits$Scientific))
  # 11 species not in either Elton database
  
  #synonyms use in Elton for the 11 taxa in question
  synnonyms = c("Spermophilus lateralis", "Carpodacus purpureus", 
                "Spermophilus tridecemlineatus", "Zoothera naevia",
                "Pipilo crissalis", "Pipilo fuscus", "Spermophilus beecheyi", 
                "Parus atricapillus", "Parus gambeli", "Carduelis psaltria", 
                "Carduelis tristis")
  not.in.all.df = data_frame(sp_new = as.character(not.in.all), 
                             Scientific = synnonyms)
  
  # replace synonyms in Elton for birds
  bird_traits$Scientific2 = left_join(bird_traits, not.in.all.df, by = "Scientific") %>% 
    mutate(sp_new2 = ifelse(is.na(sp_new), Scientific, sp_new)) %>% 
    pull(sp_new2)
  
  # replace synonyms in Elton for mammals
  mammal_traits$Scientific2 = left_join(mammal_traits, not.in.all.df, by = "Scientific") %>% 
    mutate(sp_new2 = ifelse(is.na(sp_new), Scientific, sp_new)) %>% 
    pull(sp_new2)
  
  
  # Trophic traits
  # extract trophic traits for the species in the bergmann's study
  bird_reduced = filter(bird_traits, Scientific2 %in% species$Species) %>% 
    select(Scientific2, `Diet-Inv`:`Diet-PlantO`) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "Scientific2")
  
  mammal_reduced = filter(mammal_traits, Scientific2%in% species$Species) %>% 
    select(Scientific2, `Diet-Inv`:`Diet-PlantO`) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "Scientific2")
  
  ## PCoA for Birds
  bird_dist = vegan::vegdist(bird_reduced, "bray") # bray-curtis distance
  cmd = cmdscale(bird_dist, k = 10, eig = TRUE) #run PCoA
  
  
  #Output Table
  eigenvalues <- cmd$eig[1:10]
  propVar <- eigenvalues / sum(eigenvalues)
  cumVar <- cumsum(propVar)
  PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  PCoA_Table
  
  # Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
  
  # plot the first two PCoA axes:
  bird_trophic_1 <- cmd$points[, 1]
  bird_trophic_2 <- cmd$points[, 2]
  
  plot(bird_trophic_1, bird_trophic_2, type = "n", xlab = "PCoA 1 (37%)", 
       ylab = "PCoA 2 (27%)", main = "Bird Tophic Traits",
       xlim = range(bird_trophic_1) * 1.2,
       ylim = range(bird_trophic_2) * 1.2)
  text(bird_trophic_1, bird_trophic_2, labels = rownames(cmd$points), cex = .9)
  
  # another way to plot:
  ordiplot(scores(cmd)[, c(1, 2)], type = "n", cex = 1, xlab = "PCoA 1 (37%)",
           ylab = "PCoA 2 (27%)", main = "Bird Tophic Traits")
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  col.vars <- wascores(cmd$points[, 1:2], bird_reduced)
  text(col.vars, rownames(col.vars), cex = .7, col = "red")
  
  ## PCoA for mammals
  mammal_dist = vegdist(mammal_reduced, "bray") # bray-curtis distance
  cmd <- cmdscale(mammal_dist, k=10, eig=TRUE) #run PCoA
  
  #Output Table
  eigenvalues<-cmd$eig[1:10]
  propVar<-eigenvalues/sum(eigenvalues)
  cumVar<-cumsum(propVar)
  PCoA_Table<-cbind(eigenvalues,propVar,cumVar)
  PCoA_Table
  
  # Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
  
  # plot the first two PCoA axes:
  
  mammal_trophic_1 <- cmd$points[, 1]
  mammal_trophic_2 <- cmd$points[, 2]
  plot(mammal_trophic_1, mammal_trophic_2, xlab = "Coordinate 1", ylab = "Coordinate 2",
       xlim = range(mammal_trophic_1) * 1.2, ylim = range(mammal_trophic_2) * 1.2,
       type = "n")
  text(mammal_trophic_1, mammal_trophic_2, labels = rownames(cmd$points), cex = .9)
  
  #another way to plot:
  ordiplot(scores(cmd)[, c(1, 2)], type = "n", cex = 1, xlab = "PCoA 1 (44%)",
           ylab = "PCoA 2  (25%)", main = " Mammal Trophic Traits")
  abline(h = 0, lty = 3)
  abline(v = 0, lty = 3)
  col.vars <- wascores(cmd$points[, 1:2], mammal_reduced)
  text(col.vars, rownames(col.vars), cex = .7, col = "red")
  
  # Habitat traits
  # extract trophic traits for bird species in bergmann's list
  bird_reduced <- filter(bird_traits, Scientific2 %in% species$Species) %>% 
    select(Scientific2, `ForStrat-watbelowsurf`:`PelagicSpecialist`) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "Scientific2")
  
  #for mammals, break out the foraging strata column into seperate columns for each level
  mammal_habitat = mutate(mammal_traits, d = 1, ForStrat.Value = paste0("ForStrat_", `ForStrat-Value`)) %>% 
    tidyr::spread(key = ForStrat.Value, value = d, fill = 0) %>% 
    as.tibble()
  
  #extract trophic traits for mammal species in bergmann's list
  mammal_reduced = filter(mammal_habitat, Scientific2%in% species$Species) %>% 
    select(Scientific2, ForStrat_A, ForStrat_Ar, ForStrat_G, ForStrat_S) %>% 
    as.data.frame() %>% 
    column_to_rownames(var = "Scientific2")
  # "Scientific2" "ForStrat_A"  "ForStrat_Ar" "ForStrat_G"  "ForStrat_S"
  
  ##PCoA for Birds
  bird_dist <- vegdist(bird_reduced, "bray")#bray-curtis distance
  cmd <- cmdscale(bird_dist, k = 10, eig = TRUE) #run PCoA
  
  #Output Table
  eigenvalues <- cmd$eig[1:10]
  propVar <- eigenvalues / sum(eigenvalues)
  cumVar <- cumsum(propVar)
  PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  PCoA_Table
  
  #Scree plot:
  plot(eigenvalues)
  lines(lowess(eigenvalues))
  
  # plot the first two PCoA axes:
  bird_habitat_1<-cmd$points[,1]
  bird_habitat_2<-cmd$points[,2]
  
  ## PCoA for mammals
  mammal_dist <- vegdist(mammal_reduced, "bray")#bray-curtis distance
  cmd <- cmdscale(mammal_dist, k = 10, eig = TRUE) #run PCoA
  
  #Output Table
  eigenvalues <- cmd$eig[1:10]
  propVar <- eigenvalues / sum(eigenvalues)
  cumVar <- cumsum(propVar)
  PCoA_Table <- cbind(eigenvalues, propVar, cumVar)
  PCoA_Table
  
  mammal_habitat_1<-cmd$points[,1]
  mammal_habitat_2<-cmd$points[,2]
  
  elton_axes_bird <- cbind(bird_habitat_1, bird_habitat_2, bird_trophic_1, bird_trophic_2) %>% 
    as.data.frame() %>% 
    setNames(c("habitat_axis1", "habitat_axis2","trophic_axis1", "trophic_axis2")) %>% 
    rownames_to_column("Species") %>% 
    as.tibble()
  write_csv(elton_axes_bird, "data_output/Elton_axes_bird.csv")
  
  elton_axes_mammal <- cbind(mammal_habitat_1, mammal_habitat_2, mammal_trophic_1, mammal_trophic_2) %>% 
    as.data.frame() %>% 
    setNames(c("habitat_axis1", "habitat_axis2","trophic_axis1", "trophic_axis2")) %>% 
    rownames_to_column("Species") %>% 
    as.tibble()
  write_csv(elton_axes_mammal, "data_output/Elton_axes_mammal.csv")
}

# combine all data ----
if(!file.exists('data_output/bergmann.rds')){
  sps = list.files('data_output/cleaned_bodymass_data/', pattern = 'csv$', full.names = F)
  sps = data_frame(files = as.character(1:length(sps)), fnames = sps)
  files = list.files('data_output/cleaned_bodymass_data', pattern = 'csv$', full.names = T)
  d = map_dfr(files, read_csv, .id = 'files') %>% 
    left_join(sps, by = 'files') %>% 
    rename(long = decimallongitude, lat = decimallatitude, sp = sciname_cor,
           length_mm = lengthinmm, mass_g = massing) %>% 
    select(-files, -Ids, -scientificname)
  # reduce to 175 sp
  sp2 = read_csv('data_raw/ReducedListSpecies.csv') %>% 
    mutate(Species = gsub(' ', '_', Species)) %>% 
    pull(Species)
  d = filter(d, sp %in% sp2) %>%
    rename(season = seasonality) %>% 
    mutate(mass_g_log10 = log10(mass_g),
           yr_80 = ifelse(year < 1980, 'pre_1980', 'post_1980'),
           yr_80 = factor(yr_80, levels = c('pre_1980', 'post_1980')),
           season = factor(season, levels = c('Spring', 'Summer', 'Fall', 'Winter')))
  saveRDS(d, file = 'data_output/bergmann.rds')
} else {
  d = readRDS('data/bergmann.rds')
}