library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)   ## custom fonts

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
  "mice", "patchwork" ,"projpred","geiger", "caper", "phylolm", "knitr", "ggmice", "picante", "broom"# names of the packages required placed here as character objects
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

foo <- function(file_names) {
  lapply(file_names, load, environment())
  ls()
}

# DATA ====================================================================

## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data for all ages 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 

# Load partial models

load("./revision_analysis/results/weibull/best_model_weibull_review.RData")

partial_models <- results_total

##### Load dispersal  with trait data for all ages
dispersal_analysis_partial_age <- read_csv2("data/data_process/dispersal_analysis_partial_age.csv")

names(dispersal_analysis_partial_age)

summary(dispersal_analysis_partial_age)

## Manipulation -----------------------------------------------------------

age_type <- unique(dispersal_traits_total$type)
dispersal_total <- c()


for ( i in 1:length(age_type)) {
  dispersal_traits <- dispersal_traits_total %>% 
    filter(type== age_type[i]) %>% 
    dplyr::rename(label= species) %>% 
    distinct(label, .keep_all = TRUE) 
  
  distance_total_functions_join <- distance_total_functions %>%
    filter(type== age_type[i]) %>% 
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
    #scale(.) %>% 
    as.data.frame(.) 
  
  
  label <- dispersal_traits$label
  #names(dispersal_analysis)
  
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
  
  
  ####################################
  
  dispersal_analysis_partial <- left_join(dispersal_analysis_partial,distance_total_functions_join, by= c("label") )
  dispersal_analysis_partial$age <- age_type[i]
  
  dispersal_total <- rbind(dispersal_total, dispersal_analysis_partial)
  
  
}
  

age_type <- unique(dispersal_traits_total$type)
dispersal_total_complete <- c()


for ( i in 1:length(age_type)) {
  dispersal_traits <- dispersal_traits_total %>% 
    filter(type== age_type[i]) %>% 
    dplyr::rename(label= species) %>% 
    distinct(label, .keep_all = TRUE) 
  
  distance_total_functions_join <- distance_total_functions %>%
    filter(type== age_type[i]) %>% 
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
    #scale(.) %>% 
    as.data.frame(.) 
  
  
  label <- dispersal_traits$label
  #names(dispersal_analysis)
  
  dispersal_analysis <- cbind(dispersal_analysis, dispersal_traits[, c("Territoriality.x", "Migration.1", "Diet.niche.position", "Habitat.niche.position.along.forest.open.area.gradient", "Habitat.niche.breadth")]) %>% 
    dplyr::rename("habita_for"= "Habitat.niche.position.along.forest.open.area.gradient",
                  "habitat_niche_breadth"="Habitat.niche.breadth",
                  "diet"= "Diet.niche.position", 
                  "territoriality" = "Territoriality.x")
  
  # Run bayesian analysis we do not need to delete NA
  dispersal_analysis <- cbind(label, dispersal_analysis)
  
  dispersal_analysis_partial <- left_join(dispersal_analysis,distance_total_functions_join, by= c("label") )
  dispersal_analysis_partial$age <- age_type[i]
  
  dispersal_total_complete <- rbind(dispersal_total_complete, dispersal_analysis_partial)
  
  
}

# Dispersal total: Delete species with NA in some traits
# Dispersal total complete: All species even with NA in some traits


#### Plots -----------------------------------------------------------------
## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
theme_set(theme_tidybayes())


# Plot
library(ggthemes)

ggplot(data= dispersal_total_complete, aes(x = PC1, y= Weibull_median)) +
  geom_smooth(aes(fill = age, color = age),
              method = lm,
              alpha = 0.2, size = 0.75) +
  geom_point(aes(y = Weibull_median, color = age),
             size = 1.3, alpha= 0.4) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("Life history strategy (slow-fast)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") + 
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length=unit(.20, "cm"))
  

ggsave(filename = "./revision_analysis/results/weibull/figures/scatterplots/PC1.tiff", width = 8, height = 5, device='tiff', dpi=700)


ggplot(data= dispersal_total_complete, aes(x = log_HWI, y= Weibull_median)) +
  geom_smooth(aes(fill = age, color = age),
              method = lm,
              alpha = 0.2, size = 0.75) +
  geom_point(aes(y = Weibull_median, color = age),
             size = 1.3, alpha= 0.4) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("HWI (log)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length=unit(.20, "cm"))

ggsave(filename = "./revision_analysis/results/weibull/figures/scatterplots/HWI.tiff", width = 8, height = 5, device='tiff', dpi=700)



ggplot(data= dispersal_total_complete, aes(x = body_mass, y= Weibull_median)) +
  geom_smooth(aes(fill = age, color = age),
              method = lm,
              alpha = 0.2, size = 0.75) +
  geom_point(aes(y = Weibull_median, color = age),
             size = 1.3, alpha= 0.4) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("Body mass (log)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length=unit(.20, "cm"))

ggsave(filename = "./revision_analysis/results/weibull/figures/scatterplots/body_mass.tiff", width = 8, height = 5, device='tiff', dpi=700)





#############################################
# Annotate species name

Plot
library(ggthemes)
library(ggrepel)

x_limits <- c(-4, 2)
y_limits <- c(20, NA)

ggplot(data= dispersal_total_complete%>% filter(age== "average"), aes(x = PC1, y= Weibull_median)) +
  geom_smooth(aes(fill = age, color = age),
              method = lm,
              alpha = 0.2, size = 0.75) +
  geom_point(aes(y = Weibull_median, color = age),
             size = 1.3, alpha= 0.4) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("Life history strategy (slow-fast)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") + 
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length=unit(.20, "cm"))  +
  geom_text_repel(#data= dispersal_total_complete %>% filter(Weibull_median> 20),
    aes(label = label),
    family = "Poppins",
    size = 3,
    min.segment.length = 0.5, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = 1,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .08,
    nudge_y = .1,
    color = "grey50", 
    point.padding = 0.25,
    xlim= x_limits, 
    ylim= y_limits
  ) 


#ggsave(filename = "./revision_analysis/results/weibull/figures/scatterplots/PC1.tiff", width = 8, height = 5, device='tiff', dpi=700)

x_limits <- c(2, 4)
y_limits <- c(20, NA)
ggplot(data= dispersal_total_complete %>% filter(age== "average"), aes(x = log_HWI, y= Weibull_median)) +
  geom_smooth(aes(fill = age, color = age),
              method = lm,
              alpha = 0.2, size = 0.75) +
  geom_point(aes(y = Weibull_median, color = age),
             size = 1.3, alpha= 0.4) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("HWI (log)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length=unit(.20, "cm")) +
  geom_text_repel(#data= dispersal_total_complete %>% filter(Weibull_median> 20),
                  aes(label = label),
                  family = "Poppins",
                  size = 3,
                  min.segment.length = 0.5, 
                  seed = 42, 
                  box.padding = 0.5,
                  max.overlaps = 0,
                  arrow = arrow(length = unit(0.010, "npc")),
                  nudge_x = .08,
                  nudge_y = .1,
                  color = "grey50", 
                  point.padding = 0.25,
                  xlim= x_limits, 
                  ylim= y_limits
  ) 


#ggsave(filename = "./revision_analysis/results/weibull/figures/scatterplots/HWI.tiff", width = 8, height = 5, device='tiff', dpi=700)



ggplot(data= dispersal_total_complete%>% filter(age== "average"), aes(x = body_mass, y= Weibull_median)) +
  geom_smooth(aes(fill = age, color = age),
              method = lm,
              alpha = 0.2, size = 0.75) +
  geom_point(aes(y = Weibull_median, color = age),
             size = 1.3, alpha= 0.4) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("Body mass (log)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length=unit(.20, "cm")) +
  geom_text_repel(#data= dispersal_total_complete %>% filter(Weibull_median> 20),
    aes(label = label),
    family = "Poppins",
    size = 3,
    min.segment.length = 0.5, 
    seed = 42, 
    box.padding = 0.5,
    max.overlaps = 0,
    arrow = arrow(length = unit(0.010, "npc")),
    nudge_x = .08,
    nudge_y = .1,
    color = "grey50", 
    point.padding = 0.25,
    xlim= x_limits, 
    ylim= y_limits
  ) 

#ggsave(filename = "./revision_analysis/results/weibull/figures/scatterplots/body_mass.tiff", width = 8, height = 5, device='tiff', dpi=700)












#### MEDIAN DISPERSAL #####

model_average_weibull <- partial_models$model[[7]]
model_average_hcauchy <- partial_models$model[[9]]

model_breeding_weibull <- partial_models$model[[17]]
model_breeding_hcauchy <- partial_models$model[[19]]

model_natal_weibull <- partial_models$model[[27]]
model_natal_hcauchy <- partial_models$model[[29]]

# x-values and predictions based on the log(hp)-values
theme_set(theme_tidybayes())
plot_model(model_average_weibull, type = "pred", terms = "PC1") 

#### PRUEBA SOLO CON AVERAGE #####
## plot all conditional effects
plot(conditional_effects(model_average_weibull)) +
  geom_point(dispersal_analysis_partial, aes(x= PC1, y= Weibull_median))

me <- conditional_effects(model_average_weibull, "PC1")

plot(me, plot = FALSE)[[1]] +
  scale_color_grey() +
  scale_fill_grey()

# Only if we want to fit the model in new data
##### Load dispersal  with trait data for all ages
dispersal_analysis_partial_age <- read_csv2("data/data_process/dispersal_analysis_partial_age.csv")

dispersal_analysis_partial <- dispersal_analysis_partial_age %>% 
  filter(age== "average")
nd <- 
  tibble(body_mass = seq(from = min(dispersal_analysis_partial$body_mass), to = max(dispersal_analysis_partial$body_mass), length.out = 138),
         PC1 = seq(from = min(dispersal_analysis_partial$PC1), to = max(dispersal_analysis_partial$PC1), length.out = 138),
         Latitude= seq(from = min(dispersal_analysis_partial$Latitude), to = max(dispersal_analysis_partial$Latitude), length.out = 138),
         diet= seq(from = min(dispersal_analysis_partial$diet), to = max(dispersal_analysis_partial$diet), length.out = 138),
         habita_for= seq(from = min(dispersal_analysis_partial$habita_for), to = max(dispersal_analysis_partial$habita_for), length.out = 138),
         label= dispersal_analysis_partial$label) %>% 
  drop_na()

f_b7.1 <-
  fitted(model_average_weibull, newdata = nd,  allow_new_levels= TRUE) %>%
  as_tibble() %>%
  bind_cols(nd)

# If we want to fit the model in the data that it has been calibrated
f_b7.1 <-
  fitted(model_average_weibull, allow_new_levels= TRUE) %>%
  as_tibble() %>%
  bind_cols(dispersal_analysis_partial)

# Plot 

ggplot(data= dispersal_analysis_partial, aes(x = PC1)) +
  geom_smooth(data = f_b7.1, 
              aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
              stat = "smooth",
              alpha = 1/4, size = 0.75) +
  geom_point(aes(y = Weibull_median),
             size = 1.3) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("Life history strategy (slow-fast)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") +
  theme(text = element_text(family = "Times"),
        legend.position = "none") 

 ggplot(data= dispersal_analysis_partial, aes(x = PC1)) +
  geom_smooth(aes(y = Weibull_median), stat = "smooth",
              alpha = 1/4, size = 0.75) +
  geom_point(aes(y = Weibull_median),
             size = 1.3) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  scale_x_continuous("Life history strategy (fast-slow)", expand = c(0, 0)) +
  ylab("Median dispersal (Weibull)") +
  theme(text = element_text(family = "Times"),
        legend.position = "none") 

########################
 # The three ages #######
 
 f_weibull_average <-
   fitted(model_average_weibull, allow_new_levels= TRUE) %>%
   as_tibble() %>%
   bind_cols(dispersal_analysis_partial_age %>% filter(age== "average")) %>% 
   mutate(age= "average")
 
 f_weibull_breeding <-
   fitted(model_breeding_weibull, allow_new_levels= TRUE) %>%
   as_tibble() %>%
   bind_cols(dispersal_analysis_partial_age %>% filter(age== "breeding"))%>% 
   mutate(age= "breeding")
 
 f_weibull_natal <-
   fitted(model_natal_weibull, allow_new_levels= TRUE) %>%
   as_tibble() %>%
   bind_cols(dispersal_analysis_partial_age %>% filter(age== "natal"))%>% 
   mutate(age= "natal")

 f <-
   rbind(f_weibull_average, f_weibull_breeding, f_weibull_natal)  
 
 # Plot
 library(ggthemes)
 
 ggplot(data= dispersal_analysis_partial_age, aes(x = PC1)) +
   geom_smooth(data = f, 
               aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                   fill = age, color = age),
               stat = "smooth",
               alpha = 0.2, size = 0.75) +
   geom_point(aes(y = Weibull_median, color = age),
              size = 1.3, alpha= 0.4) +
   scale_fill_manual(values = my_pal) + 
   scale_color_manual(values = my_pal) +
   scale_x_continuous("Life history strategy (slow-fast)", expand = c(0, 0)) +
   ylab("Median dispersal (Weibull)") #+
   facet_wrap(~age)
   
   ggsave(filename = "Exports/PC1.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   
   
 ggplot(data= dispersal_analysis_partial_age, aes(x = body_mass)) +
     geom_smooth(data = f, 
                 aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                     fill = age, color = age),
                 stat = "smooth",
                 alpha = 0.2, size = 0.75) +
     geom_point(aes(y = Weibull_median, color = age),
                size = 1.3, alpha= 0.5) +
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_continuous("Body mass (log)", expand = c(0, 0)) +
     ylab("Median dispersal (Weibull)") #+
   facet_wrap(~age)
   ggsave(filename = "Exports/body_mass.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   ggplot(data= dispersal_analysis_partial_age, aes(x = HWI)) +
     geom_smooth(data = f, 
                 aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                     fill = age, color = age),
                 stat = "smooth",
                 alpha = 0.2, size = 0.75) +
     geom_point(aes(y = Weibull_median, color = age),
                size = 1.3, alpha= 0.5) +
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_continuous("HWI", expand = c(0, 0)) +
     ylab("Median dispersal (Weibull)") #+
   facet_wrap(~age)
   ggsave(filename = "Exports/HWI.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   ggplot(data= dispersal_analysis_partial_age, aes(x = diet)) +
     geom_smooth(data = f, 
                 aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                     fill = age, color = age),
                 stat = "smooth",
                 alpha = 0.2, size = 0.75) +
     geom_point(aes(y = Weibull_median, color = age),
                size = 1.3, alpha= 0.5) +
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_continuous("Diet", expand = c(0.1, 0.1)) +
     ylab("Median dispersal (Weibull)") #+
   facet_wrap(~age)
   ggsave(filename = "Exports/Diet.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   ggplot(data= dispersal_analysis_partial_age) +
     geom_boxplot(aes(x = as.factor(diet),y = Weibull_median,
                     fill = age),
                  outlier.shape = NA, alpha = .5, width = .8, colour = "black") +
     geom_point(aes(x = diet, y = Weibull_median, color = age, group=age),
                position = position_jitterdodge(jitter.width = .2,
                                                jitter.height = 0,
                                                dodge.width = 1),
       shape = 95, size = 3, alpha = .8
     ) + 
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_discrete("Diet (Herb-Carn)", expand = c(0.2, 0.2)) +
     ylab("Median dispersal (Weibull)") #+
   facet_wrap(~age)
   
   
   ggplot(data= dispersal_analysis_partial_age, aes(x = Latitude)) +
     geom_smooth(data = f, 
                 aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                     fill = age, color = age),
                 stat = "smooth",
                 alpha = 0.2, size = 0.75) +
     geom_point(aes(y = Weibull_upper_distance, color = age),
                size = 1.3, alpha= 0.5) +
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_continuous("Latitude", expand = c(0.1, 0.1)) +
     ylab("Median dispersal (Weibull)")
   ggsave(filename = "Exports/Latitude.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   ggplot(data= dispersal_analysis_partial_age, aes(x = habita_for)) +
     geom_smooth(data = f, 
                 aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                     fill = age, color = age),
                 stat = "smooth",
                 alpha = 0.2, size = 0.75) +
     geom_point(aes(y = Weibull_upper_distance, color = age),
                size = 1.3, alpha= 0.5) +
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_continuous("Habitat (forest-open)", expand = c(0.1, 0.1)) +
     ylab("Median dispersal (Weibull)")
   ggsave(filename = "Exports/Habitat.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   
   ggplot(data= dispersal_analysis_partial_age, aes(x = distance_mig)) +
     geom_smooth(data = f, 
                 aes(y = Estimate, ymin = Q2.5, ymax = Q97.5,
                     fill = age, color = age),
                 stat = "smooth",
                 alpha = 0.2, size = 0.75) +
     geom_point(aes(y = Weibull_median, color = age),
                size = 1.3, alpha= 0.5) +
     scale_fill_manual(values = my_pal) + 
     scale_color_manual(values = my_pal) +
     scale_x_continuous("Migration distance", expand = c(0, 0)) +
     ylab("Median dispersal (Weibull)") #+
   facet_wrap(~age)
   
   ggsave(filename = "Exports/migration_distance.tiff", width = 8, height = 5, device='tiff', dpi=700)
   
   
   
   