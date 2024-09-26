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
  "readr","dplyr", "caret","tidyr", "tidyverse","mgcv", "MuMIn", "purrr", "reshape2", "lattice", "car", "ape", "geiger", "phytools", "nlme", "raster", "ggplot2", "sjPlot", "MCMCglmm", "plotMCMC", "tidybayes" , "plotMCMC", "loo", "brms",
  "mice", "projpred","geiger", "caper", "phylolm", "knitr", "ggmice", "picante", "broom", "performance", 
  "DHARMa", "DHARMa.helpers", "kableExtra"# names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

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
  filter(type== "natal") %>% 
  dplyr::rename(label= species) %>% 
  distinct(label, .keep_all = TRUE) 

distance_total_functions_join <- distance_total_functions %>%
  filter(type== "natal") %>% 
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

##########################
# Select only passerines
##########################
names(dispersal_traits)
unique(dispersal_traits$Order)
dispersal_analysis_passerines <- dispersal_traits %>% 
  filter(Order== "Passeriformes")

################################################

dispersal_analysis <- dispersal_analysis_passerines %>% 
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


label <- dispersal_analysis_passerines$label
names(dispersal_analysis)

dispersal_analysis <- cbind(dispersal_analysis, dispersal_analysis_passerines[, c("Territoriality.x", "Migration.1", "Diet.niche.position", "Habitat.niche.position.along.forest.open.area.gradient", "Habitat.niche.breadth")]) %>% 
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


####################################

dispersal_analysis_partial <- left_join(dispersal_analysis_partial,distance_total_functions_join, by= c("label") )



##################################
# BRMS model ####################
# Use the partial dataset

names(dispersal_analysis_partial)

st_full <- c("HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude", "body_mass:habita_for", "body_mass:PC1", "body_mass:diet")
#st1 <- paste(st, collapse = " + ")
#st2 <- paste0(paste(st, collapse = " + "), " + (1|gr(label, cov = A))")
st_full_log <- c("log_HWI", "log_body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude", "log_body_mass:habita_for", "log_body_mass:PC1", "log_body_mass:diet",  "distance_mig:Latitude")
st <- c("HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")


data <- dispersal_analysis_partial
A <- ape::vcv.phylo(dispersal_tree_partial)
phylo <- dispersal_tree_partial

data$Weibull_median <- as.integer(data$Weibull_median)
data$Weibull_upper_distance <- as.integer(data$Weibull_upper_distance)
data$Weibull_median_log <- as.integer(data$Weibull_median_log)
data$Weibull_upper_distance_log <- as.integer(data$Weibull_upper_distance_log)


#### MEDIAN DISPERSAL ######

model_dispersal_gauss_log_complete <- brm(
  paste0("Weibull_median_log~", paste(st_full_log, collapse = " + ")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

model_dispersal_gauss_log_complete <- add_criterion(model_dispersal_gauss_log_complete, "loo", moment_match = TRUE)
simres <- dh_check_brms(model_dispersal_gauss_log_complete, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/natal/median_gauss_complete_dharma.png')
plot(simres)
dev.off()


#### Model selection median dispersal ######

# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
ncores <- parallel::detectCores(logical = FALSE)
doParallel::registerDoParallel(ncores)
# Final cv_varsel() run:
cvvs <- cv_varsel(
  model_dispersal_gauss_log_complete,
  cv_method = "loo",
  ### Only for the sake of speed (not recommended in general):
  method = "forward",
  nclusters_pred = 20,
  ###
  nterms_max = 11,
  parallel = TRUE,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = TRUE
  ### 
)
# Tear down the CV parallelization setup:
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()
plot(cvvs, stats = "mlpd", deltas = TRUE)
dev.copy(png,'revision_analysis/results/weibull/natal/median_variable_selection.png')
plot(cvvs)
dev.off()
size_decided <- suggest_size(cvvs, stat = "mlpd")

smmry <- summary(cvvs, stats = "mlpd", type = c("mean", "lower", "upper"),
                 deltas = TRUE)
print(smmry, digits = 1)
rk <- ranking(cvvs)
( pr_rk <- cv_proportions(rk) )
rk[["fulldata"]]
plot(pr_rk)
( predictors_final <- head(rk[["fulldata"]], size_decided) )
plot(cv_proportions(rk, cumulate = TRUE))

# Model selection without cv

vs1 <- varsel(model_dispersal_gauss_log_complete, method="forward", cv_method="loo", refit_prj = TRUE)
modsize_decided <- suggest_size(vs1)
soltrms <- solution_terms(vs1)
ranking_terms <- ranking(vs1, nterms_max = modsize_decided)
soltrms_final <- head(soltrms, modsize_decided) 
variable_sel_table <- kable(vs1$predictor_ranking)

## Run with predictors selected


model_dispersal_gauss_log_short2 <- brm(
  paste0("Weibull_median_log~", paste(soltrms_final, collapse = " + ")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

model_dispersal_gauss_log_short2 <- add_criterion(model_dispersal_gauss_log_short2, "loo", moment_match = TRUE)
simres <- dh_check_brms(model_dispersal_gauss_log_short2, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/natal/median_gauss_short_dharma.png')
plot(simres)
dev.off()

model_dispersal_HWI <- brm(
  paste0("Weibull_median_log~",paste0(paste("HWI", collapse = " + "))),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

model_dispersal_HWI <- add_criterion(model_dispersal_HWI, "loo", moment_match = TRUE)
simres <- dh_check_brms(model_dispersal_HWI, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/natal/median_HWI_dharma.png')
plot(simres)
dev.off()

null_median_model <- brm(
  paste0("Weibull_median_log~", paste0(paste("1", collapse = " + "))),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)
null_median_model <- add_criterion(null_median_model, "loo", moment_match = TRUE)


loo_compare(model_dispersal_gauss_log_short2,  model_dispersal_gauss_log_complete,  criterion = "loo")

pp <-loo_compare( model_dispersal_gauss_log_short2, model_dispersal_gauss_log_complete, model_dispersal_HWI,null_median_model ,  criterion = "loo")



model_dispersal_gauss_log_final <- brm(
  paste0("Weibull_median_log~", paste0(paste(soltrms_final, collapse = " + "), " + (1|gr(label, cov = A))")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

simres <- dh_check_brms(model_dispersal_gauss_log_final, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/natal/median_phylo_dharma.png')
plot(simres)
dev.off()

#### The best one is the complete model for the median
best_median_model <- model_dispersal_gauss_log_final

##### LONG DISTANCE DISPERSAL #######

model_long_dispersal_gauss_log_complete <- brm(
  paste0("Weibull_upper_distance_log~", paste(st_full_log, collapse = " + ")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

model_long_dispersal_gauss_log_complete <- add_criterion(model_long_dispersal_gauss_log_complete, "loo", moment_match = TRUE)

simres <- dh_check_brms(model_long_dispersal_gauss_log_complete, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/natal/long_gauss_complete_dharma.png')
plot(simres)
dev.off()

#### Model selection median dispersal ######


# For running projpred's CV in parallel (see cv_varsel()'s argument `parallel`):
ncores <- parallel::detectCores(logical = FALSE)
doParallel::registerDoParallel(ncores)
# Final cv_varsel() run:
cvvs_long <- cv_varsel(
  model_long_dispersal_gauss_log_complete,
  cv_method = "loo",
  ### Only for the sake of speed (not recommended in general):
  method = "forward",
  nclusters_pred = 20,
  ###
  nterms_max = 11,
  parallel = TRUE,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = TRUE
  ### 
)
# Tear down the CV parallelization setup:
doParallel::stopImplicitCluster()
foreach::registerDoSEQ()
plot(cvvs_long, stats = "mlpd", deltas = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/natal/variable_selection_long.png')
plot(cvvs_long)
dev.off()

size_decided_long <- suggest_size(cvvs_long, stat = "mlpd")

smmry_long <- summary(cvvs_long, stats = "mlpd", type = c("mean", "lower", "upper"),
                      deltas = TRUE)
print(smmry_long, digits = 1)
rk_long <- ranking(cvvs_long)
( pr_rk_long <- cv_proportions(rk_long) )
rk_long[["fulldata"]]
plot(pr_rk_long)
( predictors_final_long <- head(rk_long[["fulldata"]], size_decided_long) )
plot(cv_proportions(rk_long, cumulate = TRUE))

# Model selection without cv

vs1_long <- varsel(model_long_dispersal_gauss_log_complete, method="forward", cv_method="loo", refit_prj = TRUE)
plot(vs1_long, stats = "mlpd", deltas = TRUE)
modsize_decided_long <- suggest_size(vs1_long)
soltrms_long <- solution_terms(vs1_long)
ranking_terms_long <- ranking(vs1_long, nterms_max = modsize_decided_long)
soltrms_final_long <- head(soltrms_long, modsize_decided_long) 
variable_sel_table_long <- kable(vs1_long$predictor_ranking)

## Run with predictors selected

model_long_dispersal_gauss_log_short2 <- brm(
  paste0("Weibull_upper_distance_log~", paste(soltrms_final_long, collapse = " + ")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

model_long_dispersal_gauss_log_short2 <- add_criterion(model_long_dispersal_gauss_log_short2, "loo", moment_match = TRUE)
simres <- dh_check_brms(model_long_dispersal_gauss_log_short2, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/breeding/long_gauss_short_dharma.png')
plot(simres)
dev.off()

model_long_dispersal_HWI <- brm(
  paste0("Weibull_upper_distance_log~",paste0(paste("HWI", collapse = " + "))),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

model_long_dispersal_HWI <- add_criterion(model_long_dispersal_HWI, "loo", moment_match = TRUE)
simres <- dh_check_brms(model_long_dispersal_HWI, integer = TRUE)
dev.copy(png,'revision_analysis/results/weibull/natal/long_HWI_dharma.png')
plot(simres)
dev.off()

null_long_model <- brm(
  paste0("Weibull_upper_distance_log~", paste0(paste("1", collapse = " + "))),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)
null_long_model <- add_criterion(null_long_model, "loo", moment_match = TRUE)

compare_long <- loo_compare(model_long_dispersal_gauss_log_short2, model_long_dispersal_gauss_log_complete, null_long_model, model_long_dispersal_HWI,  criterion = "loo")


library(kableExtra)

#write.csv2(compare_long, "revision_analysis/results/weibull/natal/natal_comparison_long.csv", row.names= T)

model_long_dispersal_gauss_log_final <- brm(
  paste0("Weibull_upper_distance_log~", paste0(paste(soltrms_final_long, collapse = " + "), " + (1|gr(label, cov = A))")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

simres <- dh_check_brms(model_long_dispersal_gauss_log_final, integer = TRUE)
#dev.copy(png,'revision_analysis/results/weibull/breeding/long_phylo_dharma.png')
plot(simres)
dev.off()


#### The best one is the short model for the long
best_long_model <- model_long_dispersal_gauss_log_final


## Save results
## Save results

results_natal_weibull <- tibble(model = list(best_median_model, best_long_model), variable_selection= list(cvvs, cvvs_long), age = "natal", function_t= "best_model", type= c("median", "long") )
save(results_natal_weibull, file = "revision_analysis/results/weibull/natal/model_natal_weibull_passerines.RData")

results_complete_natal_models <- tibble(model = list(model_dispersal_gauss_log_complete, model_long_dispersal_gauss_log_complete), age = "natal", function_t= "complete_model", type= c("median", "long") )
save(results_complete_natal_models, file = "revision_analysis/results/weibull/natal/model_natal_complete_weibull_passerines.RData")
