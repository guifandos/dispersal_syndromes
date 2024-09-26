#' ####################################################################### #
#' PROJECT: Dispersal syndromes on European birds
#' CONTENTS: 
#'  - This code is to analyze interspecific patterns of covariation between dispersal and other traits (‘dispersal syndromes’) considering shared evolutionary history 
#'  DEPENDENCIES:
#'  - Code documents needed to execute this document
#'  - Data files: species dispersal distances and trait databases
#'  
#' AUTHOR: Guillermo Fandos
#' ####################################################################### #

# PREAMBLE ================================================================
rm(list=ls())

## Directories ------------------------------------------------------------
### Define dicrectories in relation to project directory
Dir.Base <- getwd()
Dir.Data <- file.path(Dir.Base, "data")
Dir.Exports <- file.path(Dir.Base, "Exports")
### Create directories which aren't present yet
Dirs <- c(Dir.Data, Dir.Exports)
CreateDir <- sapply(Dirs, function(x) if(!dir.exists(x)) dir.create(x))

## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr","tidyr", "tidyverse", "mgcv", "MuMIn", "purrr", "reshape2", "lattice", "car", "ape", "geiger", "phytools", "nlme", "raster", "ggplot2", "sjPlot", "ggtree", "MCMCglmm", "plotMCMC", "tidybayes" , "plotMCMC", "loo", "brms",
  "mice", "performance"# names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)


myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

# DATA ====================================================================

## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 


## Manipulation -----------------------------------------------------------
# Check to join the dispersal dataset with the phylo tree

# Select total dispersal, excluding breeding and natal

dispersal_traits <- dispersal_traits_total %>% 
  filter(type== "average") %>% 
  dplyr::rename(label= species) %>% 
  distinct(label, .keep_all = TRUE) 

distance_total_functions_join <- distance_total_functions %>%
  filter(type== "average") %>% 
  dplyr::select(species, median,upper_distance, function_id) %>% 
  myspread(function_id, c(median,upper_distance)) %>% 
  dplyr::rename(label= species)
dispersal_traits$label <- gsub(" ", "_", dispersal_traits$label)
distance_total_functions_join$label <- gsub(" ", "_", distance_total_functions_join$label)

distance_total_functions_join2 <- distance_total_functions_join %>% 
  mutate(Exponential_median_log= log(Exponential_median +1),
         Exponential_upper_distance_log= log(Exponential_upper_distance +1),
         Gamma_median_log= log(Gamma_median +1),
         Gamma_upper_distance_log= log(Gamma_upper_distance +1),
         Hcauchy_median_log= log(Hcauchy_median +1),
         Hcauchy_upper_distance_log= log(Hcauchy_upper_distance +1),
         Weibull_median_log= log(Weibull_median +1),
         Weibull_upper_distance_log= log(Weibull_upper_distance +1),)

variables_s <- distance_total_functions_join %>% 
  dplyr::select(Exponential_median, Exponential_upper_distance, Gamma_median, Gamma_upper_distance,Hcauchy_median, Hcauchy_upper_distance,
                Weibull_median, Weibull_upper_distance) %>% 
  # scale(.) %>% 
  as.data.frame(.) 

distance_total_functions_join <- distance_total_functions_join2
#distance_total_functions_join <- cbind(distance_total_functions_join$label, distance_total_functions_join2)
colnames(distance_total_functions_join)[colnames(distance_total_functions_join) == 'distance_total_functions_join$label'] <- 'label'

dispersal_traits <-  dispersal_traits %>%
  mutate(label = dplyr::recode(label, 
                               "Apus_melba" = "Tachymarptis_melba",
                               "Chlidonias_hybridus" =  "Chlidonias_hybrida"  ,          
                               "Delichon_urbica" =  "Delichon_urbicum",           
                               "Mergus_albellus" = "Mergellus_albellus",
                               "Saxicola_torquata" = "Saxicola_torquatus",
                               "Tetrao tetrix" = "Lyrurus tetrix",             
                               "Stercorarius_skua"= "Catharacta_skua" 
  ))
row.names(dispersal_traits) <- dispersal_traits$label

name_dispersal <- unique(dispersal_traits$label)
names_tree <- rf.tree$tip.label
species_name <- base::intersect(name_dispersal,names_tree)
dispersal_traits <- dispersal_traits %>% 
  filter(label %in% species_name)
dispersal_tree <- keep.tip(rf.tree, species_name)
dispersal_tree <- compute.brlen(dispersal_tree, method = "Grafen")
name.check(dispersal_tree, dispersal_traits) 
################################################
################################################

dispersal_analysis <- dispersal_traits %>% 
  dplyr::rename("habita_for"= "Habitat.niche.position.along.forest.open.area.gradient",
                "habitat_niche_breadth"="Habitat.niche.breadth",
                "humid_grad"= "Position.along.humidity.gradient",
                "human_set"="Position.along.humidity.gradient",
                "diet"= "Diet.niche.position",
                "LHS"= "Life.history.strategy",
                "climatic_pos"="Climatic.niche.position...C.",
                "climatic_breadth"="Climatic.niche.breadth...C.",
                "range_size"="Breeding.range.size..km2.",
                "distance_mig"="Migration.distance..km.",
                "body_mass"= "Body.mass..log.") %>% 
  mutate(dispersal_distance= log(median +1),
         lon_dispersal_distance= log(upper_distance +1) ) %>% 
  #dplyr::select(dispersal_distance, lon_dispersal_distance, body_mass, HWI, LHS, climatic_pos, climatic_breadth, range_size, distance_mig, Life.span, Age.of.first.breeding, Age.of.independence, TarsusU_MEAN, LengthU_MEAN, WeightU_MEAN, Clutch_MEAN) %>% 
  #dplyr::select(dispersal_distance, lon_dispersal_distance, body_mass, HWI, LHS, climatic_pos, climatic_breadth, range_size, distance_mig, PC1, PC2, Latitude, AnnualTemp, TempRange, AnnualPrecip, PrecipRange) %>% 
  dplyr::select(dispersal_distance, lon_dispersal_distance, body_mass, HWI, distance_mig, PC1, PC2, Latitude, AnnualTemp, TempRange, AnnualPrecip, PrecipRange) 

dispersal_analysis <- dispersal_analysis %>% 
  mutate(log_body_mass= log(body_mass +1),
         log_HWI= log(HWI +1)) %>% 
  modify_if(., is.character, as.numeric) %>% 
  scale(.) %>% 
  as.data.frame(.) 


label <- dispersal_traits$label
names(dispersal_analysis)

dispersal_analysis <- cbind(dispersal_analysis, dispersal_traits[, c("Territoriality.x", "Migration.1", "Diet.niche.position", "Habitat.niche.position.along.forest.open.area.gradient", "Habitat.niche.breadth")]) %>% 
  dplyr::rename("habita_for"= "Habitat.niche.position.along.forest.open.area.gradient",
                "habitat_niche_breadth"="Habitat.niche.breadth",
                "diet"= "Diet.niche.position", 
                "territoriality" = "Territoriality.x")

# Run bayesian analysis we do not need to delete NA
dispersal_analysis <- cbind(label, dispersal_analysis)

### Delete NA ####
dispersal_analysis_partial <- dispersal_analysis %>% 
  drop_na()
setdiff(dispersal_tree$tip.label, dispersal_analysis_partial$label) # Data that it is in the tree but not in the data
name.check(dispersal_tree, dispersal_analysis_partial) # All data that it is not matching
matches <- match(dispersal_analysis_partial$label, dispersal_tree$tip.label, nomatch = 0)
dispersal_analysis_partial <- subset(dispersal_analysis_partial, matches != 0)
row.names(dispersal_analysis_partial) <- dispersal_analysis_partial$label
name.check(dispersal_tree, dispersal_analysis_partial) # All data that it is not matching

ff <- name.check(dispersal_tree, dispersal_analysis_partial)
to_drop <- ff$tree_not_data
dispersal_tree_partial <-  drop.tip(dispersal_tree, to_drop)
name.check(dispersal_tree_partial, dispersal_analysis_partial) # All data that it is not matching
# ANALYSIS ================================================================


dispersal_analysis <- left_join(dispersal_analysis,distance_total_functions_join, by= c("label") )

dispersal_analysis <- dispersal_analysis %>% 
  drop_na()

dispersal_analysis$Weibull_median <- as.integer(dispersal_analysis$Weibull_median)
dispersal_analysis$Weibull_upper_distance <- as.integer(dispersal_analysis$Weibull_upper_distance)
dispersal_analysis$Weibull_median_log <- as.integer(dispersal_analysis$Weibull_median_log)
dispersal_analysis$Weibull_upper_distance_log <- as.integer(dispersal_analysis$Weibull_upper_distance_log)


A <- ape::vcv.phylo(dispersal_tree)
phylo <- dispersal_tree


# Primero correr brms_dispersal_24042021 hasta linea

order_sp <- dispersal_traits %>% 
  dplyr::select(label, Order.x, Family.x)

dispersal_analysis <- left_join(dispersal_analysis, order_sp, by= "label")
table(dispersal_analysis$Order.x)
table(dispersal_analysis$Family.x)

####### WITHIN-CLADE #######

# PAsserines

passerines <- dispersal_analysis %>% 
  filter(Order.x == "Passeriformes")
table(passerines$Family.x)


kfold <- function(x, k) {
  
  sample(nrow(x)) %>% split(1:nrow(x) %% k)
  
}

folds <-  kfold(passerines, 5)
names(folds) <- paste0("fold_", 1:5)
map_int(folds, length)


passerines_prediction <- imap_dfr(folds, function(x, y){
  
  test_ids <- x
  train_ids <- setdiff(1:nrow(passerines), test_ids)
  
  train <- passerines %>% slice(train_ids)
  
  test <- passerines %>% slice(test_ids)
  
  m <- brm(
    Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
    data = train, family = negbinomial(), 
    data2 = list(A = A),
    chains = 2, cores = 2, iter = 4000,
    control = list(adapt_delta = 0.95)
  )
  
  predictions <- predict(m, newdata = test, type = "response",
                         family = negbinomial(), allow_new_levels= TRUE,
                         #data2 = list(A = A),
                         chains = 2, cores = 2, iter = 4000,
                         control = list(adapt_delta = 0.95))
  
  predictions <- as.data.frame(predictions)
  v_res <- test$Weibull_median 
  model_test <- lm(v_res ~ predictions[, 1])
  ff <- model_performance(model_test)
  tibble(model = list(m), validation = ff, fold = y)
  
})

passerines_prediction$order <- "Passeriformes"

#### Accipitriformes ######
Accipitriformes <- dispersal_analysis %>% 
  filter(Order.x == "Accipitriformes")
table(Accipitriformes$Family.x)


kfold <- function(x, k) {
  
  sample(nrow(x)) %>% split(1:nrow(x) %% k)
  
}

folds <-  kfold(Accipitriformes, 2)
names(folds) <- paste0("fold_", 1:2)
map_int(folds, length)


Accipitriformes_prediction <- imap_dfr(folds, function(x, y){
  
  test_ids <- x
  train_ids <- setdiff(1:nrow(Accipitriformes), test_ids)
  
  train <- Accipitriformes %>% slice(train_ids)
  
  test <- Accipitriformes %>% slice(test_ids)
  
  m <- brm(
    Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
    data = train, family = negbinomial(), 
    data2 = list(A = A),
    chains = 2, cores = 2, iter = 4000,
    control = list(adapt_delta = 0.95)
  )
  
  predictions <- predict(m, newdata = test, type = "response",
                         family = negbinomial(), allow_new_levels= TRUE,
                         #data2 = list(A = A),
                         chains = 2, cores = 2, iter = 4000,
                         control = list(adapt_delta = 0.95))
  
  predictions <- as.data.frame(predictions)
  v_res <- test$Weibull_median 
  model_test <- lm(v_res ~ predictions[, 1])
  ff <- model_performance(model_test)
  tibble(model = list(m), validation = ff, fold = y)
  
})

Accipitriformes_prediction$order <- "Accipitriformes"
###

#### Anseriformes ######
Anseriformes <- dispersal_analysis %>% 
  filter(Order.x == "Anseriformes")
table(Anseriformes$Family.x)


kfold <- function(x, k) {
  
  sample(nrow(x)) %>% split(1:nrow(x) %% k)
  
}

folds <-  kfold(Anseriformes, 2)
names(folds) <- paste0("fold_", 1:2)
map_int(folds, length)


Anseriformes_prediction <- imap_dfr(folds, function(x, y){
  
  test_ids <- x
  train_ids <- setdiff(1:nrow(Anseriformes), test_ids)
  
  train <- Anseriformes %>% slice(train_ids)
  
  test <- Anseriformes %>% slice(test_ids)
  
  m <- brm(
    Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
    data = train, family = negbinomial(), 
    data2 = list(A = A),
    chains = 2, cores = 2, iter = 4000,
    control = list(adapt_delta = 0.95)
  )
  
  predictions <- predict(m, newdata = test, type = "response",
                         family = negbinomial(), allow_new_levels= TRUE,
                         #data2 = list(A = A),
                         chains = 2, cores = 2, iter = 4000,
                         control = list(adapt_delta = 0.95))
  
  predictions <- as.data.frame(predictions)
  v_res <- test$Weibull_median 
  model_test <- lm(v_res ~ predictions[, 1])
  ff <- model_performance(model_test)
  tibble(model = list(m), validation = ff, fold = y)
  
})
Anseriformes_prediction$order <- "Anseriformes"
#### Charadriiformes ######
Charadriiformes <- dispersal_analysis %>% 
  filter(Order.x == "Charadriiformes")
table(Charadriiformes$Family.x.x)


kfold <- function(x, k) {
  
  sample(nrow(x)) %>% split(1:nrow(x) %% k)
  
}

folds <-  kfold(Charadriiformes, 2)
names(folds) <- paste0("fold_", 1:2)
map_int(folds, length)


Charadriiformes_prediction <- imap_dfr(folds, function(x, y){
  
  test_ids <- x
  train_ids <- setdiff(1:nrow(Charadriiformes), test_ids)
  
  train <- Charadriiformes %>% slice(train_ids)
  
  test <- Charadriiformes %>% slice(test_ids)
  
  m <- brm(
    Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
    data = train, family = negbinomial(), 
    data2 = list(A = A),
    chains = 2, cores = 2, iter = 4000,
    control = list(adapt_delta = 0.95)
  )
  
  predictions <- predict(m, newdata = test, type = "response",
                         family = negbinomial(), allow_new_levels= TRUE,
                         #data2 = list(A = A),
                         chains = 2, cores = 2, iter = 4000,
                         control = list(adapt_delta = 0.95))
  
  predictions <- as.data.frame(predictions)
  v_res <- test$Weibull_median 
  model_test <- lm(v_res ~ predictions[, 1])
  ff <- model_performance(model_test)
  tibble(model = list(m), validation = ff, fold = y)
  
})
Charadriiformes_prediction$order <- "Charadriiformes"

within_clade_tibble <- bind_rows(passerines_prediction, Accipitriformes_prediction, Anseriformes_prediction, Charadriiformes_prediction)
within_clade_tibble$type <- "within"

#### BETWEEN- CLADE #####


#### Split dispersal analysis in dataset
order_list <- dispersal_analysis %>% 
  group_split(Order.x)
nrow_list <- keep(order_list, function(x) nrow(x)>5 )  
map_int(nrow_list, nrow)

names_list <- lapply(nrow_list, function(x) x[1, "Order.x"])
names(nrow_list) <- unlist(lapply(names_list, '[[', 1))  

dispersal_analysis_order <- dispersal_analysis %>% 
  filter(Order.x %in% names(nrow_list) )

names(nrow_list)
#1
train <- nrow_list %>% pluck(1)
test_y <- nrow_list %>% purrr::list_modify("Accipitriformes" = NULL)
test <- do.call("rbind", test_y)
m <- brm(
  Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
  data = test, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

predictions <- predict(m, newdata = train, type = "response",
                       family = negbinomial(), allow_new_levels= TRUE,
                       #data2 = list(A = A),
                       chains = 2, cores = 2, iter = 4000,
                       control = list(adapt_delta = 0.95))

predictions <- as.data.frame(predictions)
v_res <- train$Weibull_median 
model_test <- lm(v_res ~ predictions[, 1])
ff <- model_performance(model_test)
sp_table_total <- tibble(model = list(m), validation = ff, order = "Accipitriformes")

##2
train <- nrow_list %>% pluck(2)
test_y <- nrow_list %>% purrr::list_modify("Anseriformes" = NULL)
test <- do.call("rbind", test_y)
m <- brm(
  Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
  data = test, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

predictions <- predict(m, newdata = train, type = "response",
                       family = negbinomial(), allow_new_levels= TRUE,
                       #data2 = list(A = A),
                       chains = 2, cores = 2, iter = 4000,
                       control = list(adapt_delta = 0.95))

predictions <- as.data.frame(predictions)
v_res <- train$Weibull_median 
model_test <- lm(v_res ~ predictions[, 1])
ff <- model_performance(model_test)
sp_table <- tibble(model = list(m), validation = ff, order = "Anseriformes")
sp_table_total <- bind_rows(sp_table_total, sp_table)

## 3
train <- nrow_list %>% pluck(3)
test_y <- nrow_list %>% purrr::list_modify("Charadriiformes" = NULL)
test <- do.call("rbind", test_y)
m <- brm(
  Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
  data = test, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

predictions <- predict(m, newdata = train, type = "response",
                       family = negbinomial(), allow_new_levels= TRUE,
                       #data2 = list(A = A),
                       chains = 2, cores = 2, iter = 4000,
                       control = list(adapt_delta = 0.95))

predictions <- as.data.frame(predictions)
v_res <- train$Weibull_median 
model_test <- lm(v_res ~ predictions[, 1])
ff <- model_performance(model_test)
sp_table <- tibble(model = list(m), validation = ff, order = "Charadriiformes")
sp_table_total <- bind_rows(sp_table_total, sp_table)

## 4
train <- nrow_list %>% pluck(4)
test_y <- nrow_list %>% purrr::list_modify("Passeriformes" = NULL)
test <- do.call("rbind", test_y)
m <- brm(
  Weibull_median ~ body_mass + (1|gr(label, cov = A)), 
  data = test, family = negbinomial(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

predictions <- predict(m, newdata = train, type = "response",
                       family = negbinomial(), allow_new_levels= TRUE,
                       #data2 = list(A = A),
                       chains = 2, cores = 2, iter = 4000,
                       control = list(adapt_delta = 0.95))

predictions <- as.data.frame(predictions)
v_res <- train$Weibull_median 
model_test <- lm(v_res ~ predictions[, 1])
ff <- model_performance(model_test)
sp_table <- tibble(model = list(m), validation = ff, order = "Passeriformes")
sp_table_total <- bind_rows(sp_table_total, sp_table)

between_clade_tibble <- sp_table_total

between_clade_tibble$type <- "betweeen"

between_clade_tibble$fold <- "total"


prediction_total <- bind_rows(within_clade_tibble, between_clade_tibble)
names(prediction_total)
prediction_total$complexity <- "only_body_mass"

saveRDS(prediction_total, "prediction_between_within_clade_body_mass_weibull.rds")
