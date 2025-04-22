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

#' Define a custom color scale
#'
#' @param pal Name of the color palette, as part of the
#'   original palette function name.
#' @param palette Palette type, as defined in the
#'   original palette function (optional).
#' @param n Number of (first) colors to fetch from the original palette.
#' @param order A vector of color index (optional).
#' @param alpha Transparency level.
#'
#' @return A custom color scale function.
scale_color_custom <- function(pal, palette, n, order, alpha = 1) {
  pal <- getFromNamespace(paste0("pal_", pal), "ggsci")
  
  colors <- if (missing(palette)) {
    pal(alpha = alpha)(n)
  } else {
    pal(palette = palette, alpha = alpha)(n)
  }
  
  if (length(order) > length(colors)) {
    stop("The length of order exceeds the number of colors.", call. = FALSE)
  }
  colors <- if (!missing(order)) colors[order]
  
  ggplot2::scale_color_manual(values = colors)
}

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
theme_set(theme_tidybayes())
# Set global options
options(
  ggplot2.discrete.colour = ggsci::scale_color_npg,
  ggplot2.discrete.fill = ggsci::scale_fill_npg
)

# DATA ====================================================================

## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 

# Load partial models

load("./results/models/best_model_weibull_partial_interactions.RData")

partial_models <- results_total

#### Plots

#### MEDIAN DISPERSAL #####

model_average_weibull <- partial_models$model[[1]]

model_breeding_weibull <- partial_models$model[[3]]

model_natal_weibull <- partial_models$model[[5]]


data_model_average_weibull <- model_average_weibull %>%
  spread_draws(b_PC1, b_Latitude, b_body_mass, b_diet, b_habita_for, b_distance_mig, b_log_HWI, `b_body_mass:diet`, `b_PC1:body_mass`, `b_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "average", function_id = "weibull")

plot_model(model_average_weibull)
plot_model(model_average_weibull, type = "pred", terms = c("PC1", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("habita_for", "log_body_mass"))

plot_model(model_natal_weibull, type = "pred", terms = c("PC1", "log_body_mass"))
summary(model_average_weibull)
summary(model_breeding_weibull)
summary(model_natal_weibull)

data_model_breeding_weibull <- model_breeding_weibull %>%
  spread_draws(b_PC1, b_log_body_mass, b_diet, b_habita_for, `b_log_body_mass:PC1`, `b_log_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "breeding", function_id = "weibull")
plot_model(model_breeding_weibull)

data_model_natal_weibull <- model_natal_weibull %>%
  spread_draws(b_log_body_mass, b_diet, b_PC1, b_log_HWI,  `b_log_body_mass:PC1`, b_habita_for) %>%
  summarise_draws() %>% 
  mutate(age= "natal", function_id = "weibull")

plot_model(model_natal_weibull)
data_models <- rbind(data_model_average_weibull, 
                     data_model_breeding_weibull,
                     data_model_natal_weibull) 
data_models$variable <- as.factor(data_models$variable)
data_models$function_id <- as.factor(data_models$function_id)
data_models$variable <- 
  recode_factor(data_models$variable, "b_PC1"= "Life history", "b_Latitude"= "Latitude", "b_log_body_mass"= "Body mass",
                "b_body_mass"= "Body mass", "b_PC1:body_mass"= "Body mass : Life history", 
                "b_body_mass:diet"= "Body mass : Diet", "b_body_mass:habita_for"= "Body mass : Habitat",
                "b_diet"= "Diet", "b_habita_for"= "Habitat", "b_log_HWI"= "HWI", "b_distance_mig"= "Distance migration", "b_log_body_mass:PC1"= "Body mass : Life history",
                "b_PC1:log_body_mass"= "Body mass : Life history", "b_log_body_mass:diet"= "Body mass : Diet", "b_log_body_mass:habita_for"= "Body mass : Habitat", "b_distance_mig:Latitude"= "Distance migration : Latitude",
                "b_habita_for:log_body_mass"= "Body mass : Habitat")

median_models <- data_models %>% 
  mutate(descriptor= "median")


write.csv2(data_models, "./results/weibull/weibull_models_median_best.csv")

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


## Modify point size depending on the variable importance

#median_var_sel <- read.csv("Exports/results/tables_results/median_variable_selection_scale.csv")
median_var_sel <- read.csv2("./results/weibull/variable_selection/median_variable_selection.csv")
median_var_sel <- median_var_sel %>%
  group_by( age) %>%
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 


average_var_sel <- median_var_sel %>% 
  filter(age== "average") %>% 
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 
  
breeding_var_sel <- median_var_sel %>% 
  filter(age== "breeding") %>% 
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 

natal_var_sel <- median_var_sel %>% 
  filter(age== "natal") %>% 
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 

median_var_sel_scaled <- rbind(average_var_sel, breeding_var_sel, natal_var_sel)

median_var_sel_scaled$variable <- as.factor(median_var_sel_scaled$variable)

median_var_sel_scaled$variable <- 
  recode_factor(median_var_sel_scaled$ranking_fulldata, "PC1"= "Life history", "Latitude"= "Latitude", "log_body_mass"= "Body mass",
                "body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "log_HWI"= "HWI", "distance_mig"= "Distance migration",
                "log_body_mass:PC1"= "Body mass : Life history", "log_body_mass:diet"= "Body mass : Diet", "log_body_mass:habita_for"= "Body mass : Habitat",
                "body_mass:PC1"= "Body mass : Life history", "body_mass:diet"= "Body mass : Diet", "body_mass:habita_for"= "Body mass : Habitat",
                "distance_mig:Latitude"= "Distance migration : Latitude",
                "log_HWI:Latitude"= "HWI : Latitude")

#mutate(scaled_elpd = scales::rescale((elpd*-1), to = c(0, 10)))
# Join with data models and variable importance



data_models_var_median <- left_join(median_models, median_var_sel_scaled, by= c("variable", "function_id", "age"))

median_weibull <- ggplot(data_models_var_median, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, size= cv_proportions_diag, fill=age,  group= age),
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
  xlab("Estimates") +
  theme_tidybayes()
median_weibull
ggsave("./results/weibull/figures/median_weibull_variable_selection_review2.png", plot = median_weibull, width = 8, height = 6)


#### LONG DISPERSAL #####

model_average_weibull <- partial_models$model[[2]]

model_breeding_weibull <- partial_models$model[[4]]

model_natal_weibull <- partial_models$model[[6]]

plot_model(model_average_weibull, type = "pred", terms = c("PC1", "log_body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("habita_for", "body_mass"))

plot_model(model_natal_weibull, type = "pred", terms = c("habita_for", "body_mass"))
plot_model(model_natal_weibull, type = "pred", terms = c("PC1", "body_mass"))
plot_model(model_average_weibull, type = "pred", terms = c("diet", "body_mass"))


data_model_average_weibull <- model_average_weibull %>%
  spread_draws(b_body_mass, b_diet) %>%
  summarise_draws() %>% 
  mutate(age= "average", function_id = "weibull")


data_model_breeding_weibull <- model_breeding_weibull %>%
  spread_draws(b_PC1, b_log_body_mass, b_diet, b_habita_for, b_Latitude, `b_log_body_mass:diet`, `b_log_body_mass:habita_for`) %>%
  summarise_draws() %>% 
  mutate(age= "breeding", function_id = "weibull")

data_model_natal_weibull <- model_natal_weibull %>%
  spread_draws(b_log_body_mass, b_diet) %>%
  summarise_draws() %>% 
  mutate(age= "natal", function_id = "weibull")


data_models <- rbind(data_model_average_weibull, 
                     data_model_breeding_weibull,
                     data_model_natal_weibull) 
data_models$variable <- as.factor(data_models$variable)
data_models$function_id <- as.factor(data_models$function_id)


data_models$variable <- 
  recode_factor(data_models$variable, "b_PC1"= "Life history", "b_Latitude"= "Latitude", "b_log_body_mass"= "Body mass",
                "b_body_mass"= "Body mass", "b_PC1:body_mass"= "Body mass : Life history", 
                "b_body_mass:diet"= "Body mass : Diet", "b_body_mass:habita_for"= "Body mass : Habitat",
                "b_diet"= "Diet", "b_habita_for"= "Habitat", "b_log_HWI"= "HWI", "b_distance_mig"= "Distance migration", "b_log_body_mass:PC1"= "Body mass : Life history",
                "b_PC1:log_body_mass"= "Body mass : Life history", "b_log_body_mass:diet"= "Body mass : Diet", "b_log_body_mass:habita_for"= "Body mass : Habitat", "b_distance_mig:Latitude"= "Distance migration : Latitude",
                "b_habita_for:log_body_mass"= "Body mass : Habitat")


long_models <- data_models %>% 
  mutate(descriptor= "long")

write.csv2(long_models, "./results/weibull/weibull_models_long_best.csv")

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


## Modify point size depending on the variable importance

#long_var_sel <- read.csv("Exports/results/tables_results/long_variable_selection_scale.csv")
long_var_sel <- read.csv2("./results/weibull/variable_selection/long_variable_selection.csv")

long_var_sel <- long_var_sel %>%
  filter(!size==0) %>% 
  group_by(age) %>%
  mutate(scaled_mlpd = scales::rescale((mlpd*-1)))



average_var_sel <- long_var_sel %>% 
  filter(age== "average") %>% 
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 

breeding_var_sel <- long_var_sel %>% 
  filter(age== "breeding") %>% 
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 

natal_var_sel <- long_var_sel %>% 
  filter(age== "natal") %>% 
  mutate(scaled_mlpd = scales::rescale((mlpd*-1))) 

long_var_sel_scaled <- rbind(average_var_sel, breeding_var_sel, natal_var_sel)

long_var_sel_scaled$ranking_fulldata <- as.factor(long_var_sel_scaled$ranking_fulldata)

long_var_sel_scaled$variable <- 
  recode_factor(long_var_sel_scaled$ranking_fulldata, "PC1"= "Life history", "Latitude"= "Latitude", "log_body_mass"= "Body mass",
                "body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "log_HWI"= "HWI", "distance_mig"= "Distance migration",
                "log_body_mass:PC1"= "Body mass : Life history", "log_body_mass:diet"= "Body mass : Diet", "log_body_mass:habita_for"= "Body mass : Habitat",
                "body_mass:PC1"= "Body mass : Life history", "body_mass:diet"= "Body mass : Diet", "body_mass:habita_for"= "Body mass : Habitat",
                "distance_mig:Latitude"= "Distance migration : Latitude",
                "log_HWI:Latitude"= "HWI : Latitude")




#mutate(scaled_elpd = scales::rescale((elpd*-1), to = c(0, 10)))
# Join with data models and variable importance


data_models_var_long <- left_join(long_models, long_var_sel_scaled, by= c("variable", "function_id","age"))


unique(long_var_sel_scaled$variable)
unique(long_models$variable)



ggplot(data_models, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, group= age, shape= function_id),) 

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]



long_weibull <- ggplot(data_models_var_long, aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, size= cv_proportions_diag, group= age),
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
long_weibull
ggsave("./results/weibull/figures/long_weibull_variable_selection.png", plot = long_weibull, width = 8, height = 6)


#### BOTH descriptors together
data_models_var_median$descriptor <- "Median dispersal"
data_models_var_long$descriptor <- "Long-distance dispersal"
data_models_both <- rbind(data_models_var_median, data_models_var_long)

ggplot(data_models_both %>% filter(age == "average"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  # Plot for both median and long-distance dispersal using descriptor for color
  geom_pointinterval(aes(colour = descriptor, fill = descriptor),
                     position = position_dodge(width = 0.5)) + 
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +
  
  # Custom color and fill using descriptor
  scale_fill_manual(values = c("Median dispersal" = "#1f77b4",  # Blue for median dispersal
                               "Long-distance dispersal" = "#ff7f0e")) +  # Orange for long-distance dispersal
  scale_color_manual(values = c("Median dispersal" = "#1f77b4",  # Blue for median dispersal
                                "Long-distance dispersal" = "#ff7f0e")) + 
  
  xlab("Estimates") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), 
        axis.ticks = element_line(linewidth = 1),
        axis.ticks.length = unit(.20, "cm"))


log(data_models_both$cv_proportions_diag)

both_weibull <-  ggplot(data_models_both %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= age, fill=age, size= cv_proportions_diag*10, group= age),
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
  theme_classic()
both_weibull
ggsave("./results/weibull/figures/both_weibull_variable_selection.png", plot = both_weibull, width = 8, height = 6)


both_weibull <-  ggplot(data_models_both %>% filter(function_id== "weibull"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= descriptor, fill=descriptor, size= cv_proportions_diag*10, group= descriptor),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~age) +
  xlab("Estimates") +
  theme_classic()
both_weibull

library(ggsci)

both_weibull_average <-  ggplot(data_models_both %>% filter(function_id== "weibull" & age== "average"), aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour= descriptor, fill=descriptor, size= cv_proportions_diag*10, group= descriptor),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  #scale_fill_manual(values = my_pal) + 
  #scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  #facet_wrap(~age) +
  xlab("Estimates") +
  theme_classic()
both_weibull_average
ggsave("./results/weibull/figures/both_weibull_variable_selection.png", plot = both_weibull_average, width = 8, height = 6)

library(ggplot2)

# Define colors for green and purple
purple_color <- "#7F3C8D"
blue <- "#3969AC"
orange <- "#E68310"

# Create the plot
both_weibull_average <- ggplot(data_models_both %>% filter(function_id == "weibull" & age == "average"), 
                               aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour = descriptor, fill = descriptor, size = cv_proportions_diag * 10, group = descriptor),
                     position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +  # Dashed line at x = 0
  scale_colour_manual(values = c("Long-distance dispersal" = blue, "Median dispersal" = orange)) +  # Assign green and purple to the descriptors
  scale_fill_manual(values = c("Long-distance dispersal" = blue, "Median dispersal" = orange)) +  # Fill color as well
  xlab("Estimates") +
  theme_classic() +
  theme(legend.title = element_blank())  # Optional: Remove legend title for simplicity

both_weibull_average

ggsave("./results/weibull/figures/long_median_weibull_variable_selection_color.png", plot = both_weibull_average, width = 8, height = 6)


both_natal_breeding_weibull_average <- ggplot(data_models_both %>% filter(! age == "average"), 
                               aes(y = variable, x = mean, xmin = q5, xmax = q95)) +
  geom_pointinterval(aes(colour = descriptor, fill = descriptor, size = cv_proportions_diag * 10, group = descriptor),
                     position = position_dodge(width = 0.5)) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +  # Dashed line at x = 0
  scale_colour_manual(values = c("Long-distance dispersal" = blue, "Median dispersal" = orange)) +  # Assign green and purple to the descriptors
  scale_fill_manual(values = c("Long-distance dispersal" = blue, "Median dispersal" = orange)) +  # Fill color as well
  scale_linetype_manual(values = c("Long-distance dispersal" = "dashed", "Median dispersal" = "solid")) +  # Líneas discontinua para long y continua para median
  xlab("Estimates") +
  theme_classic() +
  facet_wrap(~age)+
  theme(legend.title = element_blank())  # Optional: Remove legend title for simplicity +
both_natal_breeding_weibull_average

ggsave("./results/weibull/figures/natal_breeding_weibull_variable_selection.png", plot = both_natal_breeding_weibull_average, width = 8, height = 6)

