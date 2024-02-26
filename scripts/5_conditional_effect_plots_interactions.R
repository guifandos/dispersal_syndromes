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
load("./Exports/results/brms_models/model_total_partial.RData")
partial_models <- results_total
load("./Exports/results/brms_models/model_total_partial_full_null.RData")

partial_models_full <- results_total

partial_totals <- left_join(partial_models, partial_models_full, by= c("type", "function_t", "dispersal_mode"))
#save(partial_totals, file = "Exports/results/brms_models/partial_total.RData" )
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data for all ages
dispersal_analysis_partial_age <- read_csv2("data/data_process/dispersal_analysis_partial_age.csv")



## Manipulation -----------------------------------------------------------

#### Plots
## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]
theme_set(theme_tidybayes())
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
   
   
   
   