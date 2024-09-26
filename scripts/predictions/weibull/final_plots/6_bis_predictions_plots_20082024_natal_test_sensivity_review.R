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
library(tidyverse)     ## data wrangling + ggplot2
library(colorspace)    ## adjust colors
library(rcartocolor)   ## Carto palettes
library(ggforce)       ## sina plots
library(ggdist)        ## halfeye plots
library(ggridges)      ## ridgeline plots
library(ggbeeswarm)    ## beeswarm plots
library(gghalves)      ## off-set jitter
library(systemfonts)   ## custom fonts
library(ggthemes)
library(ggrepel)

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

new_scale <- function(new_aes) {
  structure(ggplot2::standardise_aes_names(new_aes), class = "new_aes")
}

ggplot_add.new_aes <- function(object, plot, object_name) {
  plot$layers <- lapply(plot$layers, bump_aes, new_aes = object)
  plot$scales$scales <- lapply(plot$scales$scales, bump_aes, new_aes = object)
  plot$labels <- bump_aes(plot$labels, new_aes = object)
  plot
}

# DATA ====================================================================
## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 

##### Load dispersal  with trait data for all ages
dispersal_analysis_partial_age <- read_csv2("data/data_process/dispersal_analysis_partial_age.csv")



## Manipulation -----------------------------------------------------------

####################
# Load predictions files


### natal




prediction_complex_natal <- readRDS("revision_analysis/results/predictions/weibull/median/natal/prediction_between_within_clade_complex_Weibull_natal_interactions_test_sensitivity.rds")

prediction_complex_natal$complexity <- "total"


prediction_simple_b_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_body_mass_weibull_natal.rds")
prediction_simple_h_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_HWI_Weibull_natal.rds")
prediction_simple_l_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_latitude_weibull_natal.rds")
prediction_simple_d_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_distance_mig_Weibull_natal.rds")
prediction_simple_p_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_PC1_weibull_natal.rds")
prediction_simple_die_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_diet_weibull_natal.rds")
prediction_simple_hab_natal <- readRDS("Exports/results/prediction/weibull/natal/prediction_between_within_clade_habitat_Weibull_natal.rds")

prediction_natal <- bind_rows(prediction_complex_natal, prediction_simple_b_natal, prediction_simple_h_natal, prediction_simple_l_natal, prediction_simple_d_natal, prediction_simple_p_natal,
                              prediction_simple_die_natal,prediction_simple_hab_natal )



validation <- prediction_natal$validation 
validation$id <- c(1:nrow(validation)) 
prediction_natal$id <- c(1:nrow(prediction_natal)) 
validation_long <-
  prediction_natal %>% 
  dplyr::select(-validation) %>% 
  left_join(., validation, by= "id") %>% 
  mutate(R2_s = (R2 - min(R2)) / (max(R2) - min(R2)))


summary_type <- validation_long %>% 
  group_by(type) %>% 
  summarise(Mean=mean(R2), Max=max(R2), Min=min(R2), Median=median(R2), Std=sd(R2))

names(validation_long)
export_validation <- validation_long %>% 
  dplyr::select(-model) 
write.csv(export_validation, "prediction_median_natal_sensivity_test_train.csv")

###############################
# Plots

theme_set(theme_classic())
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]



ggplot(validation_long, aes(x = complexity, y = R2 , fill = type)) +
  ggdist::stat_halfeye(aes(fill = type),position = position_nudge(x = .1, y = 0), adjust = 1.5, trim = FALSE, alpha = .5, colour = NA,
                       .width = c(.5, .95))+
  stat_halfeye(
    aes(color = type), .width = c(0), slab_fill = NA
  ) +
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 0.75), clip = "on")  +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 95,
              size = 10,
              alpha = .2
  ) 

ggplot(validation_long, aes(x = complexity, y = R2 , fill = type)) +
  geom_boxplot(aes(x = complexity, y = R2, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  geom_point(aes(x = complexity, y = R2, colour = type, group= type),position = position_jitterdodge(jitter.width = .2), size = 1, shape = 20)+
  #scale_colour_brewer(palette = "Dark2")+
  #scale_fill_brewer(palette = "Dark2") +
  coord_cartesian(ylim = c(0, 0.75), clip = "on")  +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 95,
              size = 10,
              alpha = .2
  ) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) 

ggplot(validation_long, aes(x = complexity, y = R2 , fill = type)) +
  geom_boxplot(aes(x = complexity, y = R2, fill = type),outlier.shape = NA, alpha = .5, width = .8, colour = "black")+
  #geom_point(aes(x = complexity, y = R2, colour = type, group= type),position = position_jitterdodge(jitter.width = .2), size = 1, shape = 20)+
  coord_cartesian(ylim = c(0, 0.75), clip = "on")  +
  geom_point( aes(color = type),
              position = position_jitterdodge(jitter.width = 0),
              ## draw horizontal lines instead of points
              shape = 95,
              size = 8,
              alpha = .4
  ) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) 
  #scale_colour_brewer(palette = "Dark2")+
  #scale_fill_brewer(palette = "Dark2") 
  #geom_label(
  #  label=validation_long$order, 
  #  nudge_x = 0.25, nudge_y = 0.25, 
  #  check_overlap = T
  #) +
  #geom_text_repel(aes(label = order))
  
library(Rmisc)
sumrepdat <- summarySE(validation_long, measurevar = "R2", groupvars=c("type", "complexity", "order"))

summary_order <- validation_long %>% 
  group_by(type, complexity, order) %>% 
  dplyr::summarize(Mean_R2 = mean(R2, na.rm=TRUE))


ggplot(validation_long, aes(x = complexity, y = R2 , fill = type)) +
  geom_boxplot(aes(x = complexity, y = R2, fill = type),outlier.shape = NA, alpha = .5, width = .8, colour = "black")+
  #geom_point(aes(x = complexity, y = R2, colour = type, group= type),position = position_jitterdodge(jitter.width = .2), size = 1, shape = 20)+
  coord_cartesian(ylim = c(0, 0.75), clip = "on")  +
  geom_point( aes(color = type),
              position = position_jitterdodge(jitter.width = 0),
              ## draw horizontal lines instead of points
              shape = 95,
              size = 8,
              alpha = .4
  ) +
  ## median points
 # stat_halfeye(data = summary_order,
 #   aes(aes(x = complexity, y = Mean_R2, colour = type),
#        .width=c(0), slab_fill=NA)) +
  geom_point(data = summary_order, aes(x = complexity, y = Mean_R2, group = type, colour = type),
             position = position_dodge(width= .8), shape = 1, size=4) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  geom_text_repel(data=summary_order, aes(x = complexity, y = Mean_R2, label = order, group = type),
                  min.segment.length = 0, seed = 42, box.padding = 0.5, position = position_dodge(width= .8)) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) 
  
#scale_colour_brewer(palette = "Dark2")+
#scale_fill_brewer(palette = "Dark2") 
#geom_label(
#  label=validation_long$order, 
#  nudge_x = 0.25, nudge_y = 0.25, 
#  check_overlap = T
#) +
#geom_text_repel(aes(label = order))


table_prediction <- summary_order %>% 
  mutate(complexity= as.factor(complexity)) 

table_prediction$complexity <- 
  recode_factor(table_prediction$complexity, "only_PC1"= "Life history", "only_latitude"= "Latitude", "only_body_mass"= "Body mass",
                 "only_HWI"= "HWI", "only_distance_mig"= "Distance migration", "only_diet"= "Diet", "only_habitat_for"= "Habitat", "total"= "Dispersal syndrome")


ggplot(table_prediction, aes(x = complexity, y = Mean_R2 , fill = type)) +
  geom_boxplot(aes(x = complexity, y = Mean_R2, fill = type),outlier.shape = NA, alpha = .5, width = .8, colour = "black")+
  #geom_point(aes(x = complexity, y = R2, colour = type, group= type),position = position_jitterdodge(jitter.width = .2), size = 1, shape = 20)+
  coord_cartesian(ylim = c(0, 0.75), clip = "on")  +
  geom_point( aes(color = type),
              position = position_jitterdodge(jitter.width = 0),
              ## draw horizontal lines instead of points
              shape = 95,
              size = 8,
              alpha = .4
  ) +
  ## median points
  # stat_halfeye(data = summary_order,
  #   aes(aes(x = complexity, y = Mean_R2, colour = type),
  #        .width=c(0), slab_fill=NA)) +
  #geom_point(data = summary_order, aes(x = complexity, y = Mean_R2, group = type, colour = type),
  #           position = position_dodge(width= .8), shape = 1, size=4) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_text_repel(data=summary_order, aes(x = complexity, y = Mean_R2, label = order, group = type),
  #                min.segment.length = 0, seed = 42, box.padding = 0.5, position = position_dodge(width= .8)) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) 

#scale_colour_brewer(palette = "Dark2")+
#scale_fill_brewer(palette = "Dark2") 
#geom_label(
#  label=validation_long$order, 
#  nudge_x = 0.25, nudge_y = 0.25, 
#  check_overlap = T
#) +
#geom_text_repel(aes(label = order))

sumrepdat2 <- summarySE(table_prediction, measurevar = "Mean_R2", groupvars=c("type", "complexity"))




ggplot(sumrepdat2, aes(x = complexity, y = Mean_R2 , fill = type, group= type)) +
    geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 2,
              alpha = 1, 
              position=position_dodge(1)
  )  +
  geom_errorbar(aes(ymin=Mean_R2-sd, ymax=Mean_R2+sd, color = type), width=.2, 
                position=position_dodge(1))+
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(guide = guide_axis(n.dodge=2)) +
  coord_cartesian(ylim = c(0, 0.75), clip = "on")  + theme(text = element_text(size = 20))       
  


