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
library(readr)
library(Rmisc)

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


################## adopted version with weighted Means --- weighted SE
weighted.summarySE <- function(data=NULL, measurevar,  groupvars=NULL, w, na.rm=FALSE,
                               conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  #weighted - SD function!
  w.sd <- function(x, w,na.rm=TRUE )  ( (sum(w*x*x, na.rm=na.rm)/sum(w, na.rm=na.rm)) - weighted.mean(x,w, na.rm=na.rm)^2 )^.5
  
  
  # This does the summary. For each group's data frame, return a vector with
  datac <- ddply(data, groupvars,
                 .fun = function(xx, col, w) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = weighted.mean(xx[[col]], xx[[w]], na.rm=na.rm),
                     sd   = w.sd(xx[[col]], xx[[w]], na.rm=na.rm)
                   )
                 },
                 measurevar, w
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
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


###########################
# Load number of species per order
average_data <- read_csv2("revision_analysis/data_review/dispersal_analysis_average.csv")
average_orders <- as.data.frame(table(average_data$Order.x)) %>% 
  filter(Freq > 5) %>% 
  mutate(age="average")
names(average_orders) <- c("order", "number_sp", "age")

breeding_data <- read_csv2("revision_analysis/data_review/dispersal_analysis_breeding.csv")
breeding_orders <- as.data.frame(table(breeding_data$Order.x)) %>% 
  filter(Freq > 5) %>% 
  mutate(age="breeding")
names(breeding_orders) <- c("order", "number_sp", "age")
natal_data <- read_csv2("revision_analysis/data_review/dispersal_analysis_natal.csv")
natal_orders <- as.data.frame(table(natal_data$Order.x)) %>% 
  filter(Freq > 5) %>% 
  mutate(age="natal")
names(natal_orders) <- c("order", "number_sp", "age")

prediction_orders <- rbind(average_orders, breeding_orders, natal_orders) 

## Manipulation -----------------------------------------------------------

####################
# Load table predictions files

prediction_average <- read_csv("revision_analysis/results/predictions/weibull/median/average/prediction_median_average.csv") %>% 
  mutate(age="average")

prediction_natal <- read_csv("revision_analysis/results/predictions/weibull/median/natal/prediction_median_natal_complete.csv") %>% 
  mutate(age="natal") %>% 
  filter(!order== "Accipitriformes")
prediction_breeding <- read_csv("revision_analysis/results/predictions/weibull/median/breeding/prediction_median_breeding_complete.csv") %>% 
  mutate(age="breeding") %>% 
  filter(!order== "Accipitriformes" )

prediction_total <- rbind(prediction_average, prediction_natal, prediction_breeding) #%>% 
drop_na(Sigma, R2) %>% 
  filter(!AICc== Inf)
###########################

prediction_total <- left_join(prediction_total, prediction_orders, by= c("order", "age"))




theme_set(theme_classic())
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]


summary_order <- prediction_total %>% 
  group_by(age, type, complexity, order, number_sp) %>% 
  dplyr::summarize(Mean_R2 = weighted.mean(R2, number_sp,na.rm=TRUE))



table_prediction <- summary_order %>% 
  mutate(complexity= as.factor(complexity)) 

table_prediction$complexity <- 
  recode_factor(table_prediction$complexity, "only_PC1"= "Life history", "only_latitude"= "Latitude", "only_body_mass"= "Body mass", "only_Latitude"= "Latitude",
                "only_HWI"= "HWI", "only_distance_mig"= "Distance migration", "only_diet"= "Diet", "only_habitat"= "Habitat",
                "only_habitat_for"= "Habitat", "only_habita_for"= "Habitat", "total"= "Dispersal syndrome")

sumrepdat2 <- weighted.summarySE(table_prediction, measurevar = "Mean_R2", w="number_sp", groupvars=c("age", "type", "complexity"))


# Plots



ggplot(sumrepdat2, aes(x = complexity, y = Mean_R2 , fill = type, group= type)) +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 2,
              alpha = 1, 
              position=position_dodge(0.5)
  )  +
  geom_errorbar(aes(ymin=Mean_R2-se, ymax=Mean_R2+se, color = type), width=.2, 
                position=position_dodge(0.5)
  )+
  geom_point(data = table_prediction, aes(x = complexity, y = Mean_R2, group = type, colour = type, shape= order), size= 3, alpha=0.5,
             position=position_dodge(0.5)) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_text_repel(data=table_prediction, aes(x = complexity, y = Mean_R2, label = order, group = type),
  #               min.segment.length = 0, seed = 42, box.padding = 0.5) +
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values=c(15, 17, 4, 8))+
  scale_x_discrete(guide = guide_axis(n.dodge=1)) +
  coord_cartesian(ylim = c(0, 1), clip = "on")  + theme(text = element_text(size = 20))  +
  facet_wrap(vars(age), scales = "free")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Fig_3_median_predictions", width = 8, height = 5, device='tiff', dpi=800)

###############################
ggplot(sumrepdat2, aes(x = complexity, y = Mean_R2 , fill = type, group= type)) +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 2,
              alpha = 1, 
              position=position_dodge(0.5)
  )  +
  geom_errorbar(aes(ymin=Mean_R2-se, ymax=Mean_R2+se, color = type), width=.2, 
                position=position_dodge(0.5)
  )+
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_x_discrete(guide = guide_axis(n.dodge=1)) +
  coord_cartesian(ylim = c(0, 1), clip = "on")  + theme(text = element_text(size = 20))  +
  facet_wrap(vars(age), scales = "free")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



# Average

ggplot(sumrepdat2 %>% filter(age== "average"), aes(x = complexity, y = Mean_R2 , fill = type, group= type)) +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 3,
              alpha = 1, 
              position=position_dodge(0.5)
  )  +
  geom_errorbar(aes(ymin=Mean_R2-se, ymax=Mean_R2+se, color = type), width=.2, 
                position=position_dodge(0.5)
  )+
  geom_point(data = table_prediction %>% filter(age== "average"), aes(x = complexity, y = Mean_R2, group = type, colour = type, shape= order), size= 2, alpha=0.5,
             position=position_dodge(0.5)) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_text_repel(data=table_prediction, aes(x = complexity, y = Mean_R2, label = order, group = type),
  #               min.segment.length = 0, seed = 42, box.padding = 0.5) +
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values=c(15, 17, 4, 8))+
  scale_x_discrete(guide = guide_axis(n.dodge=1)) +
  coord_cartesian(ylim = c(0, 1), clip = "on")  + theme(text = element_text(size = 20))  +
  #facet_wrap(vars(age), scales = "free")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Fig_3_median_average_prediction.tiff", width = 8, height = 5, device='tiff', dpi=700)


ggplot(sumrepdat2 %>% filter(age== "natal"), aes(x = complexity, y = Mean_R2 , fill = type, group= type)) +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 3,
              alpha = 1, 
              position=position_dodge(0.5)
  )  +
  geom_errorbar(aes(ymin=Mean_R2-se, ymax=Mean_R2+se, color = type), width=.2, 
                position=position_dodge(0.5)
  )+
  geom_point(data = table_prediction %>% filter(age== "natal"), aes(x = complexity, y = Mean_R2, group = type, colour = type, shape= order), size= 2, alpha=0.5,
             position=position_dodge(0.5)) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_text_repel(data=table_prediction, aes(x = complexity, y = Mean_R2, label = order, group = type),
  #               min.segment.length = 0, seed = 42, box.padding = 0.5) +
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values=c(17, 4, 8))+
  scale_x_discrete(guide = guide_axis(n.dodge=1)) +
  coord_cartesian(ylim = c(0, 1), clip = "on")  + theme(text = element_text(size = 20))  +
  #facet_wrap(vars(age), scales = "free")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Fig_3_median_natal_prediction.tiff", width = 8, height = 5, device='tiff', dpi=700)


ggplot(sumrepdat2 %>% filter(age== "breeding"), aes(x = complexity, y = Mean_R2 , fill = type, group= type)) +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 3,
              alpha = 1, 
              position=position_dodge(0.5)
  )  +
  geom_errorbar(aes(ymin=Mean_R2-se, ymax=Mean_R2+se, color = type), width=.2, 
                position=position_dodge(0.5)
  )+
  geom_point(data = table_prediction %>% filter(age== "breeding" ), aes(x = complexity, y = Mean_R2, group = type, colour = type, shape= order), size= 2, alpha=0.5,
             position=position_dodge(0.5)) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_text_repel(data=table_prediction, aes(x = complexity, y = Mean_R2, label = order, group = type),
  #               min.segment.length = 0, seed = 42, box.padding = 0.5) +
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values=c(17, 4, 8))+
  scale_x_discrete(guide = guide_axis(n.dodge=1)) +
  coord_cartesian(ylim = c(0, 1), clip = "on")  + theme(text = element_text(size = 20))  +
  #facet_wrap(vars(age), scales = "free")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Fig_3_median_breeding_prediction.tiff", width = 8, height = 5, device='tiff', dpi=700)


###########################
### Only passerines
###########################

prediction_total_passerines <- rbind(prediction_average, prediction_natal, prediction_breeding) %>% 
  filter(order== "Passeriformes")
###########################




theme_set(theme_classic())
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]

sumrepdat_passerines <- summarySE(prediction_total_passerines, measurevar = "R2", groupvars=c("age","type", "complexity"))

summary_passerines <- prediction_total_passerines %>% 
  group_by(age, type, complexity) %>% 
  dplyr::summarize(Mean_R2 = mean(R2, na.rm=TRUE))


table_prediction_passerines <- summary_passerines %>% 
  mutate(complexity= as.factor(complexity)) 

table_prediction_passerines$complexity <- 
  recode_factor(table_prediction_passerines$complexity, "only_PC1"= "Life history", "only_latitude"= "Latitude", "only_body_mass"= "Body mass", "only_Latitude"= "Latitude",
                "only_HWI"= "HWI", "only_distance_mig"= "Distance migration", "only_diet"= "Diet", "only_habitat"= "Habitat",
                "only_habitat_for"= "Habitat", "only_habita_for"= "Habitat", "total"= "Dispersal syndrome")

prediction_total_passerines$complexity <- 
  recode_factor(prediction_total_passerines$complexity, "only_PC1"= "Life history", "only_latitude"= "Latitude", "only_body_mass"= "Body mass", "only_Latitude"= "Latitude",
                "only_HWI"= "HWI", "only_distance_mig"= "Distance migration", "only_diet"= "Diet", "only_habitat"= "Habitat",
                "only_habitat_for"= "Habitat", "only_habita_for"= "Habitat", "total"= "Dispersal syndrome")

sumrepdat2_passerines <- summarySE(prediction_total_passerines, measurevar = "R2", groupvars=c("age", "type", "complexity"))

ggplot(sumrepdat2_passerines , aes(x = complexity, y = R2 , fill = type, group= type)) +
  geom_point( aes(color = type),
              ## draw horizontal lines instead of points
              shape = 16,
              size = 2,
              alpha = 1, 
              position=position_dodge(0.5)
  )  +
  geom_errorbar(aes(ymin=R2-se, ymax=R2+se, color = type), width=.2, 
                position=position_dodge(0.5)
  )+
  geom_point(data = table_prediction_passerines, aes(x = complexity, y = Mean_R2, group = type, colour = type), size= 3, alpha=0.5,
             position=position_dodge(0.5)) +
  #geom_point(summary_order, aes(x = complexity, y = R2, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_text_repel(data=table_prediction, aes(x = complexity, y = Mean_R2, label = order, group = type),
  #               min.segment.length = 0, seed = 42, box.padding = 0.5) +
  #geom_point(aes(x = as.numeric(R2)-.15, y = complexity, colour = type),position = position_jitter(width = .005), size = 1, shape = 20)+
  #geom_point(aes(colour = type),size = 1, shape = 20)+
  #geom_boxplot(aes(x = R2_s, y = complexity, fill = type),outlier.shape = NA, alpha = .5, width = .5, colour = "black")+
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2") +
  scale_shape_manual(values=c(15, 17, 4, 8))+
  scale_x_discrete(guide = guide_axis(n.dodge=1)) +
  coord_cartesian(ylim = c(0, 1), clip = "on")  + theme(text = element_text(size = 20))  +
  facet_wrap(vars(age), scales = "free")  +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(filename = "Fig_3_median_prediction_passerines.tiff", width = 8, height = 5, device='tiff', dpi=700)
