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
  "readr","dplyr", "caret","tidyr", "tidyverse","mgcv", "MuMIn", "purrr", "reshape2", "lattice", "car", "ape", "geiger", "phytools", "nlme", "raster", "ggplot2", "sjPlot", "MCMCglmm", "plotMCMC", "tidybayes" , "plotMCMC", "loo", "brms",
  "mice", "projpred","geiger", "caper", "phylolm", "knitr", "ggmice", "picante", "broom", "performance", 
  "DHARMa", "DHARMa.helpers"# names of the packages required placed here as character objects
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

## Loading ----------------------------------------------------------------
# Load tree data
rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")

##### Load dispersal  with trait data 
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
names(dispersal_traits_total)

distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 


## Manipulation -----------------------------------------------------------
# Check to join the dispersal dataset with the phylo tree

# Select total dispersal, excluding breeding and natal

dispersal_traits <- dispersal_traits_total %>% 
  filter(type== "average") %>% 
  dplyr::rename(label= species) %>% 
  distinct(label, .keep_all = TRUE) 

distance_total_functions_join <- distance_total_functions %>%
  filter(type== "average") %>% 
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

data <- dispersal_analysis_partial
A <- ape::vcv.phylo(dispersal_tree_partial)
phylo <- dispersal_tree_partial


# Run univariate models for body mass, HWI and PC1 -----

run_separate_distance_models <- function(data, variables, A, save_path) {
  
  # Store the models in a list
  models <- list()
  
  # Loop over each variable and fit both models (long-distance and median distance)
  for (var in variables) {
    # Formula for median distance model
    formula_median <- paste0("Weibull_median_log ~ ", var, " + (1|gr(label, cov = A))")
    
    # Formula for long-distance model
    formula_long <- paste0("Weibull_upper_distance_log ~ ", var, " + (1|gr(label, cov = A))")
    
    # Fit the median distance model
    model_median <- brm(
      formula = formula_median,
      data = data, family = gaussian(),
      data2 = list(A = A),
      chains = 2, cores = 2, iter = 4000,
      control = list(adapt_delta = 0.95),
      save_pars = save_pars(all = TRUE)
    )
    
    # Fit the long-distance model
    model_long <- brm(
      formula = formula_long,
      data = data, family = gaussian(),
      data2 = list(A = A),
      chains = 2, cores = 2, iter = 4000,
      control = list(adapt_delta = 0.95),
      save_pars = save_pars(all = TRUE)
    )
    
    # Save the models in the list
    models[[paste0(var, "_median")]] <- model_median
    models[[paste0(var, "_long")]] <- model_long
  }
  
  # Save the list of models as an RData file
  save(models, file = save_path)
  
  return(models)
}

# Example list of variables
st <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")
#st <- c("log_HWI", "body_mass")

# Define the path to save the models
save_path <- "results/weibull/average/separate_univariate_distance_models.RData"

# Run the models and save the results
models_list <- run_separate_distance_models(data, st, A, save_path)
models_list <- models

# The models are now saved to the file "separate_univariate_distance_models.RData"


# Make scatterplots of the predictions of univariate models and the data----

library(brms)
library(ggplot2)
library(dplyr)



# Función para obtener predicciones y graficar
plot_dispersal_predictions <- function(model_median, model_long, data, variable_name, x_axis_label) {
  
  # Obtener las predicciones del modelo ajustado para la variable pasada como argumento
  pred_median <- predict(model_median, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% 
    as.data.frame()
  pred_long <- predict(model_long, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% 
    as.data.frame()
  
  # Las predicciones incluyen los valores estimados y los intervalos de credibilidad (95%)
  data$predicted_median <- pred_median$Estimate
  data$lower_median <- pred_median$Q2.5 
  data$upper_median <- pred_median$Q97.5
  
  data$predicted_long <- pred_long$Estimate
  data$lower_long <- pred_long$Q2.5 
  data$upper_long <- pred_long$Q97.5
  
  # Graficar las predicciones en un scatterplot con intervalos de credibilidad al 95%
  ggplot(data = data, aes_string(x = variable_name)) +
    # Línea para dispersión mediana
    geom_line(aes(y = predicted_median, color = "Median Dispersal"), size = 1) +
    geom_point(aes(y = Weibull_median_log, color = "Median Dispersal"),
               size = 1.3, alpha = 0.6) +
    # Área sombreada para intervalo de confianza mediana
    #geom_ribbon(aes(ymin = lower_median, ymax = upper_median, fill = "95% CI Median Dispersal"),
     #           alpha = 0.2) +
    
    # Línea para dispersión larga
    geom_line(aes(y = predicted_long, color = "Long-Distance Dispersal"), size = 1, linetype = "dashed") +
    # Área sombreada para intervalo de confianza larga
    #geom_ribbon(aes(ymin = lower_long, ymax = upper_long, fill = "95% CI Long-Distance Dispersal"),
     #           alpha = 0.2) +
    geom_point(aes(y = Weibull_upper_distance_log, color = "Long-Distance Dispersal"),
               size = 1.3, alpha = 0.6, shape = 17) +
    
    # Personalización de colores y etiquetas
    #scale_fill_manual(values = c("95% CI Median Dispersal" = "#ff7f0e", "95% CI Long-Distance Dispersal" = "#1f77b4")) +
    scale_color_manual(values = c("Median Dispersal" = "#ff7f0e", "Long-Distance Dispersal" = "#1f77b4")) +
    
    scale_x_continuous(x_axis_label, expand = c(0, 0)) +
    ylab("Dispersal [log]") +
    theme_classic() +
    theme(axis.text = element_text(size = 12), 
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length = unit(.20, "cm")) +
    labs(color = "Dispersal Type", fill = "Dispersal Type")
}

# Uso de la función con body_mass y personalización del nombre del eje X
load("./results/weibull/average/separate_univariate_distance_models.RData")
models_list <- models


body_mass_plot <- plot_dispersal_predictions(models_list$body_mass_median, models_list$body_mass_long, data, "body_mass", "Body Mass [log]")


HWI_plot <- plot_dispersal_predictions(models_list$log_HWI_median, models_list$log_HWI_long, data, "log_HWI", "HWI [log]") +
   theme(legend.position = "none")



PC1_plot <- plot_dispersal_predictions(models_list$PC1_median, models_list$PC1_long, data, "PC1", "Life history strategy (slow-fast)")+
   theme(legend.position = "none")


library(patchwork)

all_plots <- PC1_plot + HWI_plot + body_mass_plot   
ggsave(all_plots, filename = "./results/weibull/figures/scatterplots/PC1_body_mass_HWI_median_long_prediction.tiff", width = 16, height = 5, device='tiff', dpi=700)

####################
# Para poner los nombres-----

library(ggthemes)
library(ggrepel)

x_limits <- c(-4, 2)
y_limits <- c(20, NA)



# Función para obtener predicciones y graficar con ggrepel para nombres de especies
plot_dispersal_predictions_label <- function(model_median, model_long, data, variable_name, x_axis_label) {
  
  # Obtener las predicciones del modelo ajustado para la variable pasada como argumento
  pred_median <- predict(model_median, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% 
    as.data.frame()
  pred_long <- predict(model_long, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% 
    as.data.frame()
  
  # Las predicciones incluyen los valores estimados y los intervalos de credibilidad (95%)
  data$predicted_median <- pred_median$Estimate
  data$lower_median <- pred_median$Q2.5 
  data$upper_median <- pred_median$Q97.5
  
  data$predicted_long <- pred_long$Estimate
  data$lower_long <- pred_long$Q2.5 
  data$upper_long <- pred_long$Q97.5
  
  # Graficar las predicciones en un scatterplot con intervalos de credibilidad al 95%
  ggplot(data = data, aes_string(x = variable_name)) +
    # Línea para dispersión mediana
    geom_line(aes(y = predicted_median, color = "Median Dispersal"), size = 1) +
    geom_point(aes(y = Weibull_median_log, color = "Median Dispersal"),
               size = 1.3, alpha = 0.6) +
    
    # Línea para dispersión larga
    geom_line(aes(y = predicted_long, color = "Long-Distance Dispersal"), size = 1, linetype = "dashed") +
    geom_point(aes(y = Weibull_upper_distance_log, color = "Long-Distance Dispersal"),
               size = 1.3, alpha = 0.6, shape = 17) +
    
    # Añadir las etiquetas de nombres de especies con ggrepel
    geom_text_repel(aes(y = Weibull_median_log, label = label), 
                    data = data, 
                    box.padding = 0.5, 
                    max.overlaps = 10, 
                    size = 3, 
                    color = "black", 
                    segment.size = 0.5, 
                    segment.color = "grey50", 
                    force = 2) +  # Adjust force for better placement
    
    # Personalización de colores y etiquetas
    scale_color_manual(values = c("Median Dispersal" = "#ff7f0e", "Long-Distance Dispersal" = "#1f77b4")) +
    
    scale_x_continuous(x_axis_label, expand = c(0, 0)) +
    ylab("Dispersal [log]") +
    theme_classic() +
    theme(axis.text = element_text(size = 12), 
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length = unit(.20, "cm")) +
    labs(color = "Dispersal Type", fill = "Dispersal Type")
}

# Example: Assume 'data' includes a 'label' column with species names


body_mass_plot_label <- plot_dispersal_predictions_label(models_list$body_mass_median, models_list$body_mass_long, data, "body_mass", "Body Mass [log]")
HWI_plot_label <- plot_dispersal_predictions_label(models_list$log_HWI_median, models_list$log_HWI_long, data, "log_HWI", "HWI [log]") +
  theme(legend.position = "none")
PC1_plot_label <- plot_dispersal_predictions_label(models_list$PC1_median, models_list$PC1_long, data, "PC1", "Life history strategy (slow-fast)") +
  theme(legend.position = "none")




library(patchwork)

all_plots <- PC1_plot_label + HWI_plot_label + body_mass_plot_label   
ggsave(all_plots, filename = "./results/weibull/figures/scatterplots/PC1_body_mass_HWI_median_long_prediction_label.tiff", width = 16, height = 5, device='tiff', dpi=700)


ggsave(body_mass_plot_label, filename = "./results/weibull/figures/scatterplots/body_mass_label.tiff", width = 8, height = 5, device='tiff', dpi=700)

# all variables------

# The type if it is signficant or not
plot_dispersal_predictions <- function(model_median, model_long, data, variable_name, x_axis_label) {
  
  # Obtener predicciones del modelo
  pred_median <- predict(model_median, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% as.data.frame()
  pred_long <- predict(model_long, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% as.data.frame()
  
  # Agregar predicciones a los datos
  data$predicted_median <- pred_median$Estimate
  data$lower_median <- pred_median$Q2.5 
  data$upper_median <- pred_median$Q97.5
  
  data$predicted_long <- pred_long$Estimate
  data$lower_long <- pred_long$Q2.5 
  data$upper_long <- pred_long$Q97.5
  
  # Extraer coeficientes y CIs correctamente usando fixef()
  fixef_median <- fixef(model_median)
  fixef_long <- fixef(model_long)
  
  # Verificar si la variable existe en los coeficientes
  if (!(variable_name %in% rownames(fixef_median))) {
    stop(paste("Variable", variable_name, "not found in median dispersal model coefficients"))
  }
  if (!(variable_name %in% rownames(fixef_long))) {
    stop(paste("Variable", variable_name, "not found in long-distance dispersal model coefficients"))
  }
  
  # Determinar si el coeficiente es significativo (IC al 95% no cruza 0)
  median_significant <- !(fixef_median[variable_name, "Q2.5"] < 0 & 
                            fixef_median[variable_name, "Q97.5"] > 0)
  
  long_significant <- !(fixef_long[variable_name, "Q2.5"] < 0 & 
                          fixef_long[variable_name, "Q97.5"] > 0)
  
  # Crear columnas de tipo de línea dentro de los datos
  data$line_type_median <- ifelse(median_significant, "solid", "dashed")
  data$line_type_long <- ifelse(long_significant, "solid", "dashed")
  
  # Generar el gráfico
  ggplot(data, aes_string(x = variable_name)) +
    geom_line(aes(y = predicted_median, color = "Median Dispersal", linetype = line_type_median), size = 1) +
    geom_point(aes(y = Weibull_median_log, color = "Median Dispersal"), size = 1.3, alpha = 0.6) +
    
    geom_line(aes(y = predicted_long, color = "Long-Distance Dispersal", linetype = line_type_long), size = 1) +
    geom_point(aes(y = Weibull_upper_distance_log, color = "Long-Distance Dispersal"), size = 1.3, alpha = 0.6, shape = 17) +
    
    scale_color_manual(values = c("Median Dispersal" = "#ff7f0e", "Long-Distance Dispersal" = "#1f77b4")) +
    scale_linetype_manual(values = c("solid" = "solid", "dashed" = "dashed")) + 
    
    scale_x_continuous(x_axis_label, expand = c(0, 0)) +
    ylab("Dispersal [log]") +
    theme_classic() +
    theme(axis.text = element_text(size = 12), 
          axis.ticks = element_line(linewidth = 1),
          axis.ticks.length = unit(.20, "cm"),
          legend.position = "right") + 
    labs(color = "Dispersal Type", linetype = "Significance")
}


library(ggplot2)
library(patchwork)

# Example list of variables
st <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")

st <- c("habita_for", "diet", "distance_mig", "Latitude")

st <- c("log_HWI", "body_mass", "PC1")

# Define a custom mapping for axis labels
axis_labels <- c(
  "log_HWI" = "HWI (log)",
  "body_mass" = "Body mass (log)",
  "habita_for" = "Habitat (forest-open)",
  "distance_mig" = "Migration distance",
  "Latitude" = "Latitude",
  "PC1" = "Life history trait (slow-fast)",
  "diet" = "Diet"
)


# Number of columns in the layout (adjust as needed)
ncol_layout <- 2  # For example, 3 columns in total
nrow_layout <- 2  # Adjust this based on your row count

# Loop through the plots and add theme to remove legends where necessary
all_plots <- lapply(1:length(st), function(i) {
  # Get the corresponding axis label
  x_label <- axis_labels[st[i]]
  
  plot <- plot_dispersal_predictions(
    model_median = models_list[[paste0(st[i], "_median")]],
    model_long = models_list[[paste0(st[i], "_long")]],
    data = data,
    variable_name = st[i],
    x_axis_label = x_label  # Use the updated axis label
  )
  
  # Calculate the row and column position in the layout
  row_position <- ceiling(i / ncol_layout)
  col_position <- i %% ncol_layout
  if (col_position == 1 && row_position == 3) {
    # Keep legend only for the plot in column 1, row 3
    plot <- plot + theme(legend.position = "right") + 
      guides(linetype = "none")  # Remove the linetype legend
  } else {
    # Remove legend for all other plots
    plot <- plot + theme(legend.position = "none") + 
      guides(linetype = "none")  # Remove the linetype legend
  }
  
  return(plot)
})

# Combine the plots using patchwork
combined_plot <- wrap_plots(all_plots, ncol = ncol_layout) + 
  plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Save the plot
ggsave(combined_plot, filename = "./results/weibull/figures/scatterplots/rest_median_long_prediction.tiff", 
       width = 15, height = 10, device='tiff', dpi=700)


ggsave(combined_plot, filename = "./results/weibull/figures/scatterplots/HWI_body_mass_PC1_median_long_prediction.tiff", 
       width = 16, height = 5, device='tiff', dpi=700)



ggsave(combined_plot, filename = "./results/weibull/figures/scatterplots/all_median_long_prediction.tiff", 
       width = 16, height = 5, device='tiff', dpi=700)







# Generate table-----
library(brms)
library(dplyr)
library(tidyr)

# Initialize an empty list to store model summaries
model_summary_list <- list()

# Loop through each variable in the models list
for (var in st) {
  # Get the median and long models for the variable
  model_median <- models_list[[paste0(var, "_median")]]
  model_long <- models_list[[paste0(var, "_long")]]
  
  # Dynamically create the parameter name for each variable
  param_name <- as.vector(paste0("b_", var))
  
  # Extract posterior samples using spread_draws (excluding intercept)
  # Dynamically access the parameter for each variable
  draws_median <- spread_draws(model_median, !!sym(param_name)) %>%
    summarise_draws()
  
  draws_long <- spread_draws(model_long,!!sym(param_name)) %>%
    summarise_draws()
 
  
  # Combine results for both models
  
  model_summary <- bind_rows(
    data.frame(Variable = var, draws_median, type = "Median"),
    data.frame(Variable = var, draws_long, type = "Long")
  )
  
  # Append the result to the list
  model_summary_list[[var]] <- model_summary
}

# Combine all model summaries into a single data frame
model_summary_df <- bind_rows(model_summary_list) %>% 
  dplyr::select(!variable)

model_summary_df$Variable <- dplyr::recode(model_summary_df$Variable,
                                    "log_HWI" = "HWI (log)",
                                    "body_mass" = "Body mass (log)",
                                    "habitat_for" = "Habitat (forest-open)",
                                    "distance_mig" = "Migration distance",
                                    "Latitude" = "Latitude",
                                    "PC1" = "Life history trait (slow-fast)",
                                    "diet" = "Diet"
)

# View the table
print(model_summary_df)

# Save the summary table to a CSV file
write.csv(model_summary_df, "./results/weibull/univariate_model_summary_results_no_intercept.csv", row.names = FALSE)
