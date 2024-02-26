library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)   ## custom fonts
library(standardize)

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


##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 

# Load partial models
load("./Exports/results/brms_models/partial_models/interactions/model_weibull_partial_interactions.RData")

partial_models <- results_total

load("./Exports/results/brms_models/model_total_partial_full_null.RData")

partial_models_full <- results_total

partial_totals <- left_join(partial_models, partial_models_full, by= c("type", "function_t", "dispersal_mode"))
#save(partial_totals, file = "Exports/results/brms_models/partial_total.RData" )
# Evaluate models
best_model <- partial_models$model[[1]]
full_model <- partial_models_full$model_full[[1]]
null_model <- partial_models_full$model_null[[1]]



var_weib_average <- as.data.frame(partial_models$variable_selection[[1]]$summary) %>% 
  mutate(age= "average", function_id = "weibull", type= "median")
var_weib_breeding <- as.data.frame(partial_models$variable_selection[[3]]$summary) %>% 
  mutate(age= "breeding", function_id = "weibull", type= "median")
var_weib_natal <- as.data.frame(partial_models$variable_selection[[5]]$summary) %>% 
  mutate(age= "natal", function_id = "weibull", type= "median")

long_weib_average <- as.data.frame(partial_models$variable_selection[[2]]$summary)%>% 
  mutate(age= "average", function_id = "weibull", type= "long")

long_weib_breeding <- as.data.frame(partial_models$variable_selection[[4]]$summary)%>% 
  mutate(age= "breeding", function_id = "weibull", type= "long")

long_weib_natal <- as.data.frame(partial_models$variable_selection[[6]]$summary) %>% 
  mutate(age= "natal", function_id = "weibull", type= "long")

median_total_sel <- rbind(var_weib_average,
                          var_weib_breeding,
                          var_weib_natal)

median_total_sel$solution_terms <- as.factor(median_total_sel$solution_terms)

median_total_sel$variable <- 
  recode_factor(median_total_sel$solution_terms, "PC1"= "Life history", "Latitude"= "Latitude", "body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "HWI"= "HWI", "distance_mig"= "Distance migration",
                "body_mass:PC1"= "Body mass : PC1", "b_body_mass:diet"= "Body mass : Diet", "b_body_mass:habita_for"= "Body mass : Habitat")
write.csv(median_total_sel, "Exports/results/tables_results/median_variable_selection_interactions.csv")
write.csv(median_total_sel, "median_variable_selection_interactions.csv")

median_total_scale <- median_total_sel %>% 
  drop_na(solution_terms) %>% 
  mutate_at(c("diff", "diff.se"), ~(scale_by(. ~ age + function_id) %>% as.vector))

write.csv(median_total_scale, "Exports/results/tables_results/median_variable_selection_scale_interactions.csv")
write.csv(median_total_scale, "median_variable_selection_scale_interactions.csv")


long_total_sel <- rbind(long_weib_average,
                        long_weib_breeding,
                        long_weib_natal)

long_total_sel$variable <- 
  recode_factor(long_total_sel$solution_terms, "PC1"= "Life history", "Latitude"= "Latitude", "body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "HWI"= "HWI", "distance_mig"= "Distance migration",
                "body_mass:PC1"= "Body mass : PC1", "b_body_mass:diet"= "Body mass : Diet", "b_body_mass:habita_for"= "Body mass : Habitat")
write.csv(long_total_sel, "Exports/results/tables_results/long_variable_selection_interactions.csv")
write.csv(long_total_sel, "long_variable_selection_interactions.csv")

long_total_scale <- long_total_sel %>% 
  drop_na(solution_terms) %>% 
  mutate_at(c("diff", "diff.se"), ~(scale_by(. ~ age + function_id) %>% as.vector))
write.csv(long_total_scale, "Exports/results/tables_results/long_variable_selection_scale_interactions.csv")
write.csv(long_total_scale, "long_variable_selection_scale_interactions.csv")

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]



ggplot(median_total_sel %>% 
         #filter(function_id== "weibull") %>% 
         drop_na(solution_terms), aes(y = variable, x = diff *-1, xmin = (diff*-1) +diff.se, xmax = (diff*-1) - diff.se)) +
  geom_pointinterval(aes(colour= age, fill=age),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  #facet_wrap(~function_id) +
  xlab("Difference elpd") 
ggsave("Exports/results/figure/median_dispersal/median_variable_selection.png", plot = median_weibull, width = 8, height = 6)


ggplot(long_total_sel %>% 
         #filter(function_id== "weibull") %>% 
         drop_na(solution_terms), aes(y = variable, x = diff *-1, xmin = (diff*-1) +diff.se, xmax = (diff*-1) - diff.se)) +
  geom_pointinterval(aes(colour= age, fill=age),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~function_id) +
  xlab("Difference elpd") 


ggsave("Exports/results/figure/long_dispersal/long_variable_selection.png", plot = median_weibull, width = 8, height = 6)



ggplot(median_total_sel, aes(y = solution_terms, x = elpd, xmin = elpd-se, xmax = elpd + se)) +
  geom_pointinterval(
                     ) +
  coord_flip()
  #geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  #facet_wrap(~function_id) +
  #xlab("Estimates") 

plot(var_sel_1)

var_sel_2 <- partial_models$variable_selection[[17]]
as.data.frame(var_sel_2$summary)
plot(var_sel_2)
