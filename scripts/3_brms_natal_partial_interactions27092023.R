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
  "readr","dplyr", "caret","tidyr", "tidyverse", "mgcv", "MuMIn", "purrr", "reshape2", "lattice", "car", "ape", "geiger", "phytools", "nlme", "raster", "ggplot2", "sjPlot", "MCMCglmm", "plotMCMC", "tidybayes" , "plotMCMC", "loo", "brms",
  "mice", "projpred","geiger", "caper", "phylolm", "knitr", "ggmice", "picante", "broom"# names of the packages required placed here as character objects
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
rm(list = ls())
## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 


## Manipulation -----------------------------------------------------------
# Check to join the dispersal dataset with the phylo tree

# Select total dispersal, excluding natal and natal

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

distance_total_functions_join <- distance_total_functions_join %>% 
  mutate(Exponential_median= log(Exponential_median +1),
         Exponential_upper_distance= log(Exponential_upper_distance +1),
         Gamma_median= log(Gamma_median +1),
         Gamma_upper_distance= log(Gamma_upper_distance +1),
         Hcauchy_median= log(Hcauchy_median +1),
         Hcauchy_upper_distance= log(Hcauchy_upper_distance +1),
         Weibull_median= log(Weibull_median +1),
         Weibull_upper_distance= log(Weibull_upper_distance +1),)

variables_s <- distance_total_functions_join %>% 
  dplyr::select(Exponential_median, Exponential_upper_distance, Gamma_median, Gamma_upper_distance,Hcauchy_median, Hcauchy_upper_distance,
                Weibull_median, Weibull_upper_distance) %>% 
  scale(.) %>% 
  as.data.frame(.) 

distance_total_functions_join <- cbind(distance_total_functions_join$label, variables_s)
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
species_name <- intersect(name_dispersal,names_tree)
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
  #dplyr::select(dispersal_distance, lon_dispersal_distance, body_mass, HWI, LHS, climatic_pos, climatic_breadth, range_size, distance_mig, Life.span, Age.of.first.natal, Age.of.independence, TarsusU_MEAN, LengthU_MEAN, WeightU_MEAN, Clutch_MEAN) %>% 
  #dplyr::select(dispersal_distance, lon_dispersal_distance, body_mass, HWI, LHS, climatic_pos, climatic_breadth, range_size, distance_mig, PC1, PC2, Latitude, AnnualTemp, TempRange, AnnualPrecip, PrecipRange) %>% 
  dplyr::select(dispersal_distance, lon_dispersal_distance, body_mass, HWI, distance_mig, PC1, PC2, Latitude, AnnualTemp, TempRange, AnnualPrecip, PrecipRange) 

dispersal_analysis <- dispersal_analysis %>% 
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

### Inputate NA #####
md.pattern(dispersal_analysis)
summary(dispersal_analysis)
# visualize the incomplete data
ggmice(dispersal_analysis, ggplot2::aes(HWI, dispersal_distance)) + ggplot2::geom_point()
ggmice(dispersal_analysis, ggplot2::aes(habitat_niche_breadth, dispersal_distance)) + ggplot2::geom_point()

histogram(dispersal_analysis$HWI)
imp <- mice(dispersal_analysis, m = 8, print = FALSE)
ggmice(imp, ggplot2::aes(HWI, dispersal_distance)) + ggplot2::geom_point() 
ggmice(imp, ggplot2::aes(habitat_niche_breadth, dispersal_distance)) + ggplot2::geom_point() 
ggmice(imp, ggplot2::aes(x = HWI, group = .imp)) +
  ggplot2::geom_density() 
summary(dispersal_analysis)
summary(with(imp, mean(habita_for)))
plot(imp)
summary(complete(imp))
densityplot(imp)
dispersal_analysis_complete <- complete(imp)

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

#######################
#### Join with the rest of dispersal distance for the different distributions
#dispersal_analysis_complete <- left_join(dispersal_analysis_complete,distance_total_functions_join, by= c("label") )

dispersal_analysis_partial <- left_join(dispersal_analysis_partial,distance_total_functions_join, by= c("label") )


##################################
# BRMS model ####################
# Use the partial dataset



st <- c("HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude", "body_mass:Latitude", "body_mass:PC1", "body_mass:diet")
#st1 <- paste(st, collapse = " + ")
#st2 <- paste0(paste(st, collapse = " + "), " + (1|gr(label, cov = A))")


data <- dispersal_analysis_partial
A <- ape::vcv.phylo(dispersal_tree_partial)
phylo <- dispersal_tree_partial

# Median dispesal

model_dispersal <- brm(
  paste0("Weibull_median~", paste(st, collapse = " + ")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)
# Variable selection
vs1 <- varsel(model_dispersal, method="forward", cv_method="loo")
modsize_decided <- suggest_size(vs1)
soltrms <- solution_terms(vs1)
soltrms_final <- head(soltrms, modsize_decided) 
variable_sel_table <- kable(vs1$summary)
#st2 <- paste0(paste(soltrms_final, collapse = " + "), " + (1|gr(label, cov = A))")
# Final model median
model_dispersal_f <- brm(
  paste0("Weibull_median~", paste0(paste(soltrms_final, collapse = " + "), " + (1|gr(label, cov = A))")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

# Long dispersal

model_long_dispersal <- brm(
  paste0("Weibull_upper_distance~", paste(st, collapse = " + ")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)
# Variable selection
vs_long <- varsel(model_long_dispersal, method="forward", cv_method="loo")
modsize_decided_long <- suggest_size(vs_long)
soltrms_long <- solution_terms(vs_long)
soltrms_final_long <- head(soltrms_long, modsize_decided_long) 
variable_sel_table_long <- kable(vs_long$summary)
#st2 <- paste0(paste(soltrms_final, collapse = " + "), " + (1|gr(label, cov = A))")
# Final model median
model_long_dispersal_f <- brm(
  paste0("Weibull_upper_distance~", paste0(paste(soltrms_final, collapse = " + "), " + (1|gr(label, cov = A))")),
  data = data, family = gaussian(), 
  data2 = list(A = A),
  chains = 2, cores = 2, iter = 4000,
  control = list(adapt_delta = 0.95)
)

results_natal_weibull <- tibble(model = list(model_dispersal_f, model_long_dispersal_f), variable_selection= list(vs1, vs_long), type = "natal", function_t= "best_model" )
save(results_natal_weibull, file = "Exports/results/model_natal_weibull_interactions.RData")
