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

# Load partial models

load("./Exports/results/brms_models/partial_models/interactions/model_weibull_partial_interactions.RData")

partial_models <- results_total

#### Plots

#### MEDIAN DISPERSAL #####

model_average_weibull <- partial_models$model[[1]]

model_breeding_weibull <- partial_models$model[[3]]

model_natal_weibull <- partial_models$model[[5]]


data_model_average_weibull <- model_average_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_body_mass, b_diet, b_habita_for, b_distance_mig, `b_body_mass:diet`, `b_PC1:body_mass`, `b_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "average", function_id = "weibull")

plot_model(model_average_weibull)
plot_model(model_average_weibull, type = "pred", terms = c("PC1", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("habita_for", "body_mass"))

plot_model(model_natal_weibull, type = "pred", terms = c("habita_for", "body_mass"))
plot_model(model_natal_weibull, type = "pred", terms = c("PC1", "body_mass"))
plot_model(model_natal_weibull, type = "pred", terms = c("diet", "body_mass"))

data_model_breeding_weibull <- model_breeding_weibull %>%
  spread_draws(b_PC1, b_body_mass, b_diet, b_habita_for) %>%
  summarise_draws() %>% 
  mutate(age= "breeding", function_id = "weibull")

data_model_natal_weibull <- model_natal_weibull %>%
  spread_draws(b_body_mass, b_diet, b_PC1, b_HWI, b_habita_for, b_Latitude, `b_body_mass:PC1`, `b_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "natal", function_id = "weibull")


data_models <- rbind(data_model_average_weibull, 
                     data_model_breeding_weibull,
                     data_model_natal_weibull) 
data_models$variable <- as.factor(data_models$variable)
data_models$function_id <- as.factor(data_models$function_id)
data_models$variable <- 
  recode_factor(data_models$variable, "b_PC1"= "Life history", "b_Latitude"= "Latitude", "b_body_mass"= "Body mass",
                "b_diet"= "Diet", "b_habita_for"= "Habitat", "b_HWI"= "HWI", "b_distance_mig"= "Distance migration", "b_body_mass:PC1"= "Body mass : PC1",
                "b_PC1:body_mass"= "Body mass : PC1", "b_body_mass:diet"= "Body mass : Diet", "b_body_mass:habita_for"= "Body mass : Habitat")

median_models <- data_models %>% 
  mutate(descriptor= "median")

write.csv(data_models, "data_models.csv")
write.csv2(data_models, "data_models_interactions_median.csv")

ggplot(data_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, group= age, shape= function_id),) 

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]


ggplot(median_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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

# Only weibull

ggplot(data_models %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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
  xlab("Estimates") 



## Modify point size depending on the variable importance

#median_var_sel <- read.csv("Exports/results/tables_results/median_variable_selection_scale.csv")
median_var_sel <- read.csv("Exports/results/tables_results/median_variable_selection_interactions.csv")
median_var_sel <- median_var_sel %>%
  filter(!size==0) %>% 
  group_by(function_id, age) %>%
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 


average_var_sel <- median_var_sel %>% 
  filter(age== "average") %>% 
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 
  
breeding_var_sel <- median_var_sel %>% 
  filter(age== "breeding") %>% 
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 

natal_var_sel <- median_var_sel %>% 
  filter(age== "natal") %>% 
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 

median_var_sel_scaled <- rbind(average_var_sel, breeding_var_sel, natal_var_sel)

#mutate(scaled_elpd = scales::rescale((elpd*-1), to = c(0, 10)))
# Join with data models and variable importance

data_models_var_median <- left_join(median_models, median_var_sel_scaled, by= c("variable", "function_id", "age"))

median_weibull <- ggplot(data_models_var_median %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, size= scaled_elpd, group= age),
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
  xlab("Estimates") 
median_weibull
ggsave("Exports/results/figure/median_dispersal/median_weibull_variable_selection_interactions.png", plot = median_weibull, width = 8, height = 6)

ggplot(data_models_var_median %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, size= (elpd), group= age),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_size_continuous(limits = c(0, 10)) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~function_id) +
  xlab("Estimates") 



#### LONG DISPERSAL #####

model_average_weibull <- partial_models$model[[2]]

model_breeding_weibull <- partial_models$model[[4]]

model_natal_weibull <- partial_models$model[[6]]

plot_model(model_average_weibull, type = "pred", terms = c("PC1", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("habita_for", "body_mass"))

plot_model(model_natal_weibull, type = "pred", terms = c("habita_for", "body_mass"))
plot_model(model_natal_weibull, type = "pred", terms = c("PC1", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "body_mass"))


data_model_average_weibull <- model_average_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_body_mass, b_diet, b_habita_for, b_distance_mig, `b_body_mass:diet`, `b_PC1:body_mass`, `b_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "average", function_id = "weibull")


data_model_breeding_weibull <- model_breeding_weibull %>%
  spread_draws(b_PC1, b_body_mass, b_diet, b_habita_for) %>%
  summarise_draws() %>% 
  mutate(age= "breeding", function_id = "weibull")

data_model_natal_weibull <- model_natal_weibull %>%
  spread_draws(b_body_mass, b_diet, b_PC1, b_HWI, b_habita_for, b_Latitude, `b_body_mass:PC1`, `b_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "natal", function_id = "weibull")


data_models <- rbind(data_model_average_weibull, 
                     data_model_breeding_weibull,
                     data_model_natal_weibull) 
data_models$variable <- as.factor(data_models$variable)
data_models$function_id <- as.factor(data_models$function_id)
data_models$variable <- 
  recode_factor(data_models$variable, "b_PC1"= "Life history", "b_Latitude"= "Latitude", "b_body_mass"= "Body mass",
                "b_diet"= "Diet", "b_habita_for"= "Habitat", "b_HWI"= "HWI", "b_distance_mig"= "Distance migration", "b_body_mass:PC1"= "Body mass : PC1",
                "b_PC1:body_mass"= "Body mass : PC1", "b_body_mass:diet"= "Body mass : Diet", "b_body_mass:habita_for"= "Body mass : Habitat")

long_models <- data_models %>% 
  mutate(descriptor= "long")

write.csv2(data_models, "data_models_interactions_long.csv")

ggplot(data_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, group= age, shape= function_id),) 

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]


ggplot(data_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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

# Only weibull

ggplot(data_models %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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
  xlab("Estimates") 


## Modify point size depending on the variable importance

#long_var_sel <- read.csv("Exports/results/tables_results/long_variable_selection_scale.csv")
long_var_sel <- read.csv("Exports/results/tables_results/long_variable_selection_interactions.csv")

long_var_sel <- long_var_sel %>%
  filter(!size==0) %>% 
  group_by(function_id, age) %>%
  mutate(scaled_elpd = scales::rescale((elpd*-1)))



average_var_sel <- long_var_sel %>% 
  filter(age== "average") %>% 
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 

breeding_var_sel <- long_var_sel %>% 
  filter(age== "breeding") %>% 
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 

natal_var_sel <- long_var_sel %>% 
  filter(age== "natal") %>% 
  mutate(scaled_elpd = scales::rescale((elpd*-1))) 

long_var_sel_scaled <- rbind(average_var_sel, breeding_var_sel, natal_var_sel)

#mutate(scaled_elpd = scales::rescale((elpd*-1), to = c(0, 10)))
# Join with data models and variable importance

data_models_var_long <- left_join(long_models, long_var_sel_scaled, by= c("variable", "function_id", "age"))




ggplot(data_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, group= age, shape= function_id),) 

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]


ggplot(data_models_var_long, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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
  xlab("Estimates") 

# Only weibull
long_weibull <- ggplot(data_models_var_long %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, size= scaled_elpd, group= age),
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
  xlab("Estimates") 
ggsave("Exports/results/figure/median_dispersal/long_weibull_variable_selection_interactions.png", plot = long_weibull, width = 8, height = 6)
#### BOTH descriptors together
data_models_var_median$descriptor <- "median dispersal"
data_models_var_long$descriptor <- "long dispersal"
data_models_both <- rbind(data_models_var_median, data_models_var_long)

ggplot(data_models_both %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
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
  facet_wrap(~descriptor) +
  xlab("Estimates") 


both_weibull <-  ggplot(data_models_both %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, size= scaled_elpd, group= age),
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
  xlab("Estimates") 

ggsave("Exports/results/figure/median_dispersal/both_weibull_variable_selection_interactions.png", plot = both_weibull, width = 8, height = 6)

