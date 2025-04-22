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
load("./results/models/best/best_model_weibull_partial_interactions.Rdata")

partial_models <- results_total

model_average_weibull <- partial_models$model[[1]]

#load("./Exports/results/brms_models/model_total_partial_full_null.RData")

#partial_models_full <- results_total

#partial_totals <- left_join(partial_models, partial_models_full, by= c("type", "function_t", "dispersal_mode"))
#save(partial_totals, file = "Exports/results/brms_models/partial_total.RData" )
# Evaluate models
#best_model <- partial_models$model[[1]]
#full_model <- partial_models_full$model_full[[1]]
#null_model <- partial_models_full$model_null[[1]]
library(posterior)

prueba <- partial_models$variable_selection[[1]]
smmry <- summary(prueba, stats = "mlpd", type = c("mean", "lower", "upper"),
                 deltas = TRUE)
smmry$perf_sub
as.data.frame(smmry$perf_sub)


smmry_average <- summary(partial_models$variable_selection[[1]], stats = "mlpd", type = c("mean", "lower", "upper"),
                         deltas = TRUE)
var_weib_average <- as.data.frame(smmry_average$perf_sub) %>% 
  mutate(age= "average", function_id = "weibull", type= "median")
smmry_breeding<- summary(partial_models$variable_selection[[3]], stats = "mlpd", type = c("mean", "lower", "upper"),
                         deltas = TRUE)
var_weib_breeding <- as.data.frame(smmry_breeding$perf_sub) %>% 
  mutate(age= "breeding", function_id = "weibull", type= "median")
smmry_natal <- summary(partial_models$variable_selection[[5]], stats = "mlpd", type = c("mean", "lower", "upper"),
                         deltas = TRUE)
var_weib_natal <- as.data.frame(smmry_natal$perf_sub) %>% 
  mutate(age= "natal", function_id = "weibull", type= "median")

smmry_long_average <- summary(partial_models$variable_selection[[2]], stats = "mlpd", type = c("mean", "lower", "upper"),
                         deltas = TRUE)
long_weib_average <- as.data.frame(smmry_long_average$perf_sub)%>% 
  mutate(age= "average", function_id = "weibull", type= "long")
smmry_long_breeding <- summary(partial_models$variable_selection[[4]], stats = "mlpd", type = c("mean", "lower", "upper"),
                              deltas = TRUE)
long_weib_breeding <- as.data.frame(smmry_long_breeding$perf_sub)%>% 
  mutate(age= "breeding", function_id = "weibull", type= "long")
smmry_long_natal <- summary(partial_models$variable_selection[[6]], stats = "mlpd", type = c("mean", "lower", "upper"),
                               deltas = TRUE)
long_weib_natal <- as.data.frame(smmry_long_natal$perf_sub) %>% 
  mutate(age= "natal", function_id = "weibull", type= "long")

median_total_sel <- rbind(var_weib_average,
                          var_weib_breeding,
                          var_weib_natal)
median_total_sel <- median_total_sel %>% 
  filter(!ranking_fulldata== "(Intercept)")

median_total_sel$ranking_fulldata <- as.factor(median_total_sel$ranking_fulldata)

median_total_sel$variable <- 
  recode_factor(median_total_sel$ranking_fulldata, "PC1"= "Life history", "Latitude"= "Latitude", "log_body_mass"= "Body mass",
                "body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "log_HWI"= "HWI", "distance_mig"= "Distance migration",
                "log_body_mass:PC1"= "Body mass : Life history", "log_body_mass:diet"= "Body mass : Diet", "log_body_mass:habita_for"= "Body mass : Habitat",
                "body_mass:PC1"= "Body mass : Life history", "body_mass:diet"= "Body mass : Diet", "body_mass:habita_for"= "Body mass : Habitat",
                "distance_mig:Latitude"= "Distance migration : Latitude",
                "log_HWI:Latitude"= "HWI : Latitude")
write.csv2(median_total_sel, "./results/weibull/median_variable_selection.csv")

library(standardize)
median_total_scale <- median_total_sel %>%
  mutate(mlpd_scale= mlpd) %>% 
  mutate_at(c("mlpd_scale"), ~(scale_by(. ~ age) %>% as.vector))

write.csv2(median_total_scale, "./results/weibull/median_variable_selection_scale.csv")


long_total_sel <- rbind(long_weib_average,
                        long_weib_breeding,
                        long_weib_natal)
long_total_sel <- long_total_sel %>% 
  filter(!ranking_fulldata== "(Intercept)")

long_total_sel$ranking_fulldata <- as.factor(long_total_sel$ranking_fulldata)


long_total_sel$variable <- 
  recode_factor(long_total_sel$ranking_fulldata, "PC1"= "Life history", "Latitude"= "Latitude", "log_body_mass"= "Body mass",
                "body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "log_HWI"= "HWI", "distance_mig"= "Distance migration",
                "log_body_mass:PC1"= "Body mass : Life history", "log_body_mass:diet"= "Body mass : Diet", "log_body_mass:habita_for"= "Body mass : Habitat",
                "body_mass:PC1"= "Body mass : Life history", "body_mass:diet"= "Body mass : Diet", "body_mass:habita_for"= "Body mass : Habitat",
                "distance_mig:Latitude"= "Distance migration : Latitude",
                "log_HWI:Latitude"= "HWI : Latitude")

write.csv2(long_total_sel, "./results/weibull/long_variable_selection.csv")

long_total_scale <- long_total_sel %>%
  mutate(mlpd_scale= mlpd) %>% 
  mutate_at(c("mlpd_scale"), ~(scale_by(. ~ age) %>% as.vector))

write.csv2(long_total_scale, "./results/weibull/long_variable_selection_scale.csv")

## custom colors
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]



median_weibull <- ggplot(median_total_sel, aes(y = variable, x = mlpd *-1, xmin = (mlpd.lower*-1), xmax = (mlpd.upper*-1))) +
  geom_pointinterval(aes(colour= age, fill=age),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  #geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~age) +
  xlab("MLPD *-1") 
ggsave("./results/weibull/figures/median_variable_selection.png", plot = median_weibull, width = 8, height = 6)


long_weibull <- ggplot(long_total_sel, aes(y = variable, x = mlpd *-1, xmin = (mlpd.lower*-1), xmax = (mlpd.upper*-1) )) +
  geom_pointinterval(aes(colour= age, fill=age),
                     position = position_dodge(
                       ## control randomness and range of jitter
                       width = 0.5
                     )
  ) +
  #geom_vline(xintercept=0, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")+
  facet_wrap(~age) +
  xlab("MLPD *-1") 


ggsave("./results/weibull/figures/long_variable_selection.png", plot = long_weibull, width = 8, height = 6)


cv_proportion_median <- ggplot(median_total_sel, aes(y = variable, x = cv_proportions_diag)) +
  geom_point(aes(colour= age, fill=age),
             position = position_dodge(
               ## control randomness and range of jitter
               width = 0.5
             )
  ) +
  geom_vline(xintercept=0.5, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~age) +
  xlab("CV Proportion") 

ggsave("./results/weibull/figures/median_variable_selection_cv_proportion.png", plot = cv_proportion_median, width = 8, height = 6)


cv_proportion_long <- ggplot(long_total_sel, aes(y = variable, x = cv_proportions_diag)) +
  geom_point(aes(colour= age, fill=age),
             position = position_dodge(
               ## control randomness and range of jitter
               width = 0.5
             )
  ) +
  geom_vline(xintercept=0.5, lty=2, size =0.5, col="grey") +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  #scale_fill_manual(values = my_pal, guide = "none")
  facet_wrap(~age) +
  xlab("CV Proportion") 
ggsave("./results/weibull/figures/long_variable_selection_cv_proportion.png", plot = cv_proportion_long, width = 8, height = 6)

##############
variable_selection_both <- rbind(median_total_sel, long_total_sel)
library(LaCroixColoR)

variable_selection_plot <- ggplot(variable_selection_both %>% mutate(lower= 0), aes(y = variable, x = cv_proportions_diag *100)) + 
  geom_vline(aes(xintercept = 50), color = "grey85") +
  geom_linerange(aes(xmin = lower, xmax = cv_proportions_diag *100, colour= type), position = position_dodge(
    ## control randomness and range of jitter
    width = 0.5 )) +
  geom_point(aes(colour= type, fill=type),
             position = position_dodge(
               ## control randomness and range of jitter
               width = 0.5
             )) +
  facet_wrap(~age) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  xlab("CV Proportion %") +
  theme_tidybayes()

  ggsave("./results/weibull/figures/variable_selection_plot_both.png", plot = variable_selection_plot, width = 8, height = 6)
  
x11()
variable_selection_plot

# Natal dispersal

natal_dispersal_variable_selection <- rbind(var_weib_natal,
                                            long_weib_natal)
natal_dispersal_variable_selection <- natal_dispersal_variable_selection %>% 
  filter(!ranking_fulldata== "(Intercept)")

natal_dispersal_variable_selection$ranking_fulldata <- as.factor(natal_dispersal_variable_selection$ranking_fulldata)

natal_dispersal_variable_selection$variable <- 
  recode_factor(natal_dispersal_variable_selection$ranking_fulldata, "PC1"= "Life history", "Latitude"= "Latitude", "log_body_mass"= "Body mass",
                "diet"= "Diet", "habita_for"= "Habitat", "log_HWI"= "HWI", "distance_mig"= "Distance migration",
                "log_body_mass:PC1"= "Body mass : PC1", "log_body_mass:diet"= "Body mass : Diet", "log_body_mass:habita_for"= "Body mass : Habitat",
                "distance_mig:Latitude"= "Distance migration : Latitude")

cv_proportion_natal <- ggplot(natal_dispersal_variable_selection %>% mutate(lower= 0), aes(y = variable, x = cv_proportions_diag *100)) +
  #geom_vline(aes(xintercept = 50), color = "grey85") +
  geom_linerange(aes(xmin = lower, xmax = cv_proportions_diag *100, colour= type), position = position_dodge(
    ## control randomness and range of jitter
    width = 0.5 )) +
  geom_point(aes(colour= type, fill=type),
             position = position_dodge(
               ## control randomness and range of jitter
               width = 0.5
             )) +
  facet_wrap(~age) +
  scale_x_continuous(limits = c(0, 100)) +
  scale_fill_manual(values = my_pal) + 
  scale_color_manual(values = my_pal) +
  xlab("CV Proportion %") +
  theme_tidybayes()
ggsave("./results/weibull/figures/median_natal_variable_selection_cv_proportion.png", plot = cv_proportion_natal, width = 8, height = 6)

########################
library(here)
library(rstanarm)
options(mc.cores = parallel::detectCores())
library(loo)
library(projpred)
library(ggplot2)
library(bayesplot)
theme_set(bayesplot::theme_default(base_family = "sans"))
library(corrplot)
library(knitr)

median_natal_validation <- partial_models$variable_selection[[5]]
median_natal_validation <- cvvs

natal_median_size <- suggest_size(median_natal_validation, stat = "mlpd")
rk_median_natal <- ranking(median_natal_validation)
( pr_rk_natal <- cv_proportions(rk_median_natal) )
( predictors_final_natal_median <- head(rk_median_natal[["fulldata"]], natal_median_size) )
plot(cv_proportions(rk_median_natal, cumulate = TRUE))
plot(median_natal_validation, stats = c('elpd', 'rmse'), deltas=FALSE)



smmry_natal <- summary(median_natal_validation, stats = "mlpd", type = c("mean", "lower", "upper"),
                       deltas = TRUE)
var_weib_natal <- as.data.frame(smmry_natal$perf_sub) %>% 
  mutate(age= "natal", function_id = "weibull", type= "median")



median_natal_validation_old <- partial_models$variable_selection_old[[5]]



smmry_natal <- summary(partial_models$variable_selection[[5]], stats = "elpd", type = c("mean", "lower", "upper"),
                       deltas = TRUE)
var_weib_natal <- as.data.frame(smmry_natal$perf_sub) %>% 
  mutate(age= "natal", function_id = "weibull", type= "median")



(nsel <- suggest_size(median_natal_validation, alpha=0.1))
(vsel <- solution_terms(median_natal_validation)[1:nsel])



projrhs <- project(median_natal_validation, nv = nsel, ns = 4000)
mcmc_areas(as.matrix(projrhs), pars = vsel)

size_decided <- suggest_size(median_natal_validation, stat = "mlpd")
rk <- ranking(median_natal_validation)
( pr_rk <- cv_proportions(rk) )
rk[["fulldata"]]
plot(pr_rk)
( predictors_final <- head(rk[["fulldata"]], size_decided) )


prj <- project(
  median_natal_validation,
  predictor_terms = predictors_final,
  ### In interactive use, we recommend not to deactivate the verbose mode:
  verbose = FALSE
  ###
)
prj_mat <- as.matrix(prj)
library(posterior)
prj_drws <- as_draws_matrix(prj_mat)
prj_smmry <- summarize_draws(
  prj_drws,
  "median", "mad", function(x) quantile(x, probs = c(0.025, 0.975))
)
# Coerce to a `data.frame` because pkgdown versions > 1.6.1 don't print the
# tibble correctly:
prj_smmry <- as.data.frame(prj_smmry)
print(prj_smmry, digits = 1)
library(bayesplot)
bayesplot_theme_set(ggplot2::theme_bw())
mcmc_intervals(prj_mat) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))


refm_mat <- as.matrix(partial_models$model[[5]])
mcmc_intervals(refm_mat, pars = colnames(prj_mat)) +
  ggplot2::coord_cartesian(xlim = c(-1.5, 1.6))



median_natal_validation_old <- partial_models$variable_selection_old[[5]]
plot(median_natal_validation_old, stats = "mlpd", deltas = TRUE)
