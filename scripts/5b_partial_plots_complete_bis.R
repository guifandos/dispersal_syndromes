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
theme_set(theme_classic())
# DATA ====================================================================

## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 

# Load complete models

load("./revision_analysis/results/weibull/complete_model_weibull_review.RData")

partial_models <- results_total_complete

#### Plots

#### MEDIAN DISPERSAL #####

model_average_weibull <- partial_models$model[[1]]

model_breeding_weibull <- partial_models$model[[3]]

model_natal_weibull <- partial_models$model[[5]]


data_model_average_weibull <- model_average_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_log_body_mass, b_diet, b_habita_for, b_distance_mig, `b_log_body_mass:diet`, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`, b_log_HWI, `b_distance_mig:Latitude`) %>%
  summarise_draws() %>% 
  mutate(age= "average", function_id = "weibull")

plot_model(model_average_weibull)
plot_model(model_average_weibull, type = "pred", terms = c("PC1", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("habita_for", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("distance_mig", "Latitude"))

plot_model(model_natal_weibull, type = "pred", terms = c("PC1", "log_body_mass"))

data_model_breeding_weibull <- model_breeding_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_log_body_mass, b_diet, b_habita_for, b_distance_mig, `b_log_body_mass:diet`, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`, b_log_HWI, `b_distance_mig:Latitude`) %>%
  summarise_draws() %>% 
  mutate(age= "breeding", function_id = "weibull")

data_model_natal_weibull <- model_natal_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_log_body_mass, b_diet, b_habita_for, b_distance_mig, `b_log_body_mass:diet`, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`, b_log_HWI, `b_distance_mig:Latitude`) %>%
  summarise_draws() %>% 
  mutate(age= "natal", function_id = "weibull")


data_models_median <- rbind(data_model_average_weibull, 
                     data_model_breeding_weibull,
                     data_model_natal_weibull) 
data_models_median$variable <- as.factor(data_models_median$variable)
data_models_median$function_id <- as.factor(data_models_median$function_id)
data_models_median$variable <- 
  recode_factor(data_models_median$variable, "b_PC1"= "Life history", "b_Latitude"= "Latitude", "b_log_body_mass"= "Body mass",
                "b_diet"= "Diet", "b_habita_for"= "Habitat", "b_log_HWI"= "HWI", "b_distance_mig"= "Distance migration", "b_log_body_mass:PC1"= "Body mass : PC1",
                "b_PC1:body_mass"= "Body mass : PC1", "b_log_body_mass:diet"= "Body mass : Diet", "b_log_body_mass:habita_for"= "Body mass : Habitat", "b_distance_mig:Latitude"= "Distance migration : Latitude")

median_models <- data_models_median %>% 
  mutate(descriptor= "median")


write.csv2(data_models_median, "./revision_analysis/results/weibull/weibull_models_median_complete.csv")

ggplot(data_models_median, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, group= age, shape= function_id),) 

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]


median_weibull <- ggplot(median_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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
 # facet_wrap(~function_id) +
  xlab("Estimates") 

ggsave("./revision_analysis/results/weibull/figures/median_weibull_complete.png", plot = median_weibull, width = 8, height = 6)


#### LONG DISPERSAL #####

model_average_weibull <- partial_models$model[[2]]

model_breeding_weibull <- partial_models$model[[4]]

model_natal_weibull <- partial_models$model[[6]]

plot_model(model_average_weibull, type = "pred", terms = c("PC1", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("habita_for", "log_body_mass"))

plot_model(model_natal_weibull, type = "pred", terms = c("habita_for", "log_body_mass"))
plot_model(model_natal_weibull, type = "pred", terms = c("PC1", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "log_body_mass"))


data_model_average_weibull <- model_average_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_log_body_mass, b_diet, b_habita_for, b_distance_mig, `b_log_body_mass:diet`, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`, b_log_HWI,  `b_distance_mig:Latitude`) %>%
  summarise_draws() %>% 
  mutate(age= "average", function_id = "weibull")


data_model_breeding_weibull <- model_breeding_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_log_body_mass, b_diet, b_habita_for, b_distance_mig, `b_log_body_mass:diet`, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`, b_log_HWI,  `b_distance_mig:Latitude`) %>%
  summarise_draws() %>% 
  mutate(age= "breeding", function_id = "weibull")

data_model_natal_weibull <- model_natal_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_log_body_mass, b_diet, b_habita_for, b_distance_mig, `b_log_body_mass:diet`, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`, b_log_HWI,  `b_distance_mig:Latitude`) %>%
  summarise_draws() %>% 
  mutate(age= "natal", function_id = "weibull")


data_models_long <- rbind(data_model_average_weibull, 
                     data_model_breeding_weibull,
                     data_model_natal_weibull) 
data_models_long$variable <- as.factor(data_models_long$variable)
data_models_long$function_id <- as.factor(data_models_long$function_id)
data_models_long$variable <- 
  recode_factor(data_models_long$variable, "b_PC1"= "Life history", "b_Latitude"= "Latitude", "b_log_body_mass"= "Body mass",
                "b_diet"= "Diet", "b_habita_for"= "Habitat", "b_log_HWI"= "HWI", "b_distance_mig"= "Distance migration", "b_log_body_mass:PC1"= "Body mass : PC1",
                "b_PC1:body_mass"= "Body mass : PC1", "b_log_body_mass:diet"= "Body mass : Diet", "b_log_body_mass:habita_for"= "Body mass : Habitat", "b_distance_mig:Latitude"= "Distance migration : Latitude")

long_models <- data_models_long %>% 
  mutate(descriptor= "long")

write.csv2(long_models, "./revision_analysis/results/weibull/weibull_models_complete_long.csv")

ggplot(data_models_long, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, group= age, shape= function_id),) 

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]


long_weibull <- ggplot(data_models_long, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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
  # facet_wrap(~function_id) +
  xlab("Estimates") 

ggsave("./revision_analysis/results/weibull/figures/long_weibull_complete.png", plot = long_weibull, width = 8, height = 6)

#### BOTH descriptors together
data_models_median$descriptor <- "median dispersal"
data_models_long$descriptor <- "long dispersal"
data_models_both <- rbind(data_models_median, data_models_long)



both_weibull <-  ggplot(data_models_both %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age,  group= age),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~descriptor) +
  xlab("Estimates") +
  scale_x_continuous(limits = c(-2, 5))
  xlim(-5,10)
both_weibull
ggsave("./revision_analysis/results/weibull/figures/both_weibull_complete.png", plot = both_weibull, width = 8, height = 6)

