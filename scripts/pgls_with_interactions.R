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
library(raster)
library(phylolm)
library(picante)
library(geiger)
library(caper)
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

##################################
# BRMS model ####################
# Use the partial dataset

names(dispersal_analysis_partial)

st_full <- c("HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude", "body_mass:habita_for", "body_mass:PC1", "body_mass:diet")
#st1 <- paste(st, collapse = " + ")
#st2 <- paste0(paste(st, collapse = " + "), " + (1|gr(label, cov = A))")
st_full_log <- c("log_HWI", "log_body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude", "log_body_mass:habita_for", "log_body_mass:PC1", "log_body_mass:diet", "distance_mig:Latitude")
st <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")



# prepare a comparative.data structure 
comp.data<-comparative.data(dispersal_tree_partial, dispersal_analysis_partial, names.col="label", vcv.dim=2, warn.dropped=TRUE)
modelo_pgls <- pgls(Weibull_median_log~ log_HWI + body_mass  + PC1 +  distance_mig + habita_for + Latitude + diet + body_mass:habita_for + body_mass:PC1 +  body_mass:diet + distance_mig:Latitude, 
                    data = comp.data, lambda="ML")

summary(modelo_pgls)
plot(modelo_pgls)
mod.l <- pgls.profile(modelo_pgls, 'lambda')
plot(mod.l)

# Extract residuals and fitted values
residuals_pgls <- residuals(modelo_pgls)
fitted_values_pgls <- fitted(modelo_pgls)

set.seed(123)  # For reproducibility
n_simulations <- 100  # Number of simulations

# Create a matrix to store simulated responses
sim_res_matrix <- matrix(NA, nrow = length(fitted_values_pgls), ncol = n_simulations)

for (i in 1:n_simulations) {
  # Sample residuals with replacement to create a new simulated response
  simulated_residuals <- sample(residuals_pgls, replace = TRUE)
  sim_res_matrix[, i] <- fitted_values_pgls + simulated_residuals
}

# Extract residuals from the PGLS model
residuals_pgls <- residuals(modelo_pgls)

# Plot histogram of residuals
hist(residuals_pgls, main = "Histogram of Residuals", xlab = "Residuals", col = "lightblue", border = "black")

# Check normality using the Shapiro-Wilk test
shapiro.test(residuals_pgls)

dharma_res <- createDHARMa(
  observedResponse = dispersal_analysis_partial$Weibull_median,
  simulatedResponse = sim_res_matrix
)

# Test for overdispersion
testOverdispersion(dharma_res)

qqnorm(residuals_pgls)
qqline(residuals_pgls, col = "red")

plot(fitted(modelo_pgls), residuals_pgls, 
     xlab = "Fitted Values", ylab = "Residuals",
     main = "Residuals vs Fitted Values")
abline(h = 0, col = "red")


# Calculate the residuals
residuals_pgls <- residuals(modelo_pgls)

# Calculate the residual variance
residual_variance <- var(residuals_pgls)

# Calculate the fitted values
fitted_values_pgls <- fitted(modelo_pgls)

# Calculate the variance of fitted values
fitted_variance <- var(fitted_values_pgls)

# If the ratio is significantly greater than 1, it indicates that the residual variance is higher than expected, suggesting overdispersion.
# Calculate the ratio of residual variance to fitted variance
residual_variance / fitted_variance

##############################################################
# Extract the coefficients (estimates) from the PGLS model
coef_estimates <- coef(modelo_pgls)
# View the summary of the model
model_summary <- summary(modelo_pgls)

# Extract the model coefficients and related statistics
coefficients_table <- model_summary$coefficients

# View the coefficients table
coefficients_table

# Create a dataframe with coefficients and CI
coef_df <- data.frame(
  term = rownames(coefficients_table),
  estimate = coefficients_table[, 1],
  se = coefficients_table[, 2],
  lower_ci = coefficients_table[, 1] - 1.96 * coefficients_table[, 2],
  upper_ci = coefficients_table[, 1] + 1.96 * coefficients_table[, 2]
)

# Load ggplot2 package
library(ggplot2)

# Plot the coefficients with 95% CI using ggplot2
ggplot(coef_df, aes(x = term, y = estimate)) +
  geom_point(color = "grey9", size = 3) +  # Points for the estimates
  geom_errorbar(aes(ymin = lower_ci, ymax = upper_ci), width = 0.2, color = "grey50") +  # CI bars
  coord_flip() +  # Flip coordinates for a horizontal plot
  labs(title = "Model Coefficients with 95% Confidence Intervals",
       x = "Predictor Variables",
       y = "Estimate") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))  # Optional: change text size


#####################################
# Model selection ####
model_selection <-dredge(modelo_pgls)

# Calculate variable importance
var_importance <- sw(model_selection)
# Convert the variable importance table to a data frame
var_importance_df1 <- as.data.frame(var_importance)
# Convert the variable importance result into a data frame
var_importance_df <- data.frame(
  variable = rownames(var_importance_df1),  # Extract row names as variable names
  importance = var_importance_df1[, 1]     # Extract importance values
)

# Plot variable importance
ggplot(var_importance_df, aes(x = reorder(variable, importance), y = importance)) +
  geom_bar(stat = "identity", fill = "#66c2a5") +  # Bar plot for importance
  coord_flip() +  # Flip coordinates to make it horizontal
  labs(title = "Variable Importance Based on Model Selection (AICc Weights)",
       x = "Variables", y = "Importance") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))  # Customize axis text size

# Model average


pp <- subset(model_selection, cumsum(results$weight) <= .95)



# Model average using all candidate models, always use revised.var = TRUE
MA.ests<-model.avg(pp, revised.var = TRUE, fit=TRUE)
summary(MA.ests)

library(data.table)

mA<-summary(MA.ests) #pulling out model averages
df1<-as.data.frame(mA$coefmat.full) #selecting full model coefficient averages

CI <- as.data.frame(confint(MA.ests, full=T)) # get confidence intervals for full model
df1$CI.min <-CI$`2.5 %` #pulling out CIs and putting into same df as coefficient estimates
df1$CI.max <-CI$`97.5 %`# order of coeffients same in both, so no mixups; but should check anyway
setDT(df1, keep.rownames = "coefficient") #put rownames into column
names(df1) <- gsub(" ", "", names(df1)) # remove spaces from column headers

ggplot(df1[-1,], aes(x = coefficient, y = Estimate)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + 
  geom_point(color = "grey9", size = 3) +  # Points for the estimates
  geom_errorbar(aes(ymin = Estimate-Std.Error, ymax = Estimate+Std.Error), width = 0.2, color = "grey50") +  # CI bars
  coord_flip() +  # Flip coordinates for a horizontal plot
  # Add red dashed vertical line at 0
  labs(title = "Model Coefficients with 95% Confidence Intervals",
       x = "Predictor Variables",
       y = "Estimate") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 12))  # Optional: change text size


######
# Plot scatterplot HWI with dispersal distance
# Extract the coefficients and confidence intervals
coef_summary <- summary(modelo_pgls)$coefficients
intercept <- coef_summary["(Intercept)", "Estimate"]
slope <- coef_summary["log_HWI", "Estimate"]
se_slope <- coef_summary["log_HWI", "Std. Error"]


# Calculate 95% confidence intervals
lower_bound <- slope - 1.96 * se_slope
upper_bound <- slope + 1.96 * se_slope

conf_int_pgls <- data_frame(lower_bound= lower_bound, upper_bound= upper_bound)

# Check if slope is significantly different from 0
significant_slope <- !(lower_bound <= 0 & upper_bound >= 0)

# Print results
cat("Slope:", slope, "\n95% CI:", lower_bound, "to", upper_bound, 
    "\nSignificant Slope:", significant_slope, "\n")

# Load ggplot2 for plotting
library(ggplot2)

# Create a sequence of predictor values for the ribbon
new_data <- data.frame(log_HWI = seq(min(dispersal_analysis_partial$log_HWI), 
                                     max(dispersal_analysis_partial$log_HWI), 
                                     length.out = 100), 
                       Weibull_median_log= seq(min(dispersal_analysis_partial$Weibull_median_log), 
                                               max(dispersal_analysis_partial$Weibull_median_log), 
                                               length.out = 100)
)

# Calculate fitted values
new_data$fit <- intercept + slope * new_data$log_HWI

# Calculate standard errors of the predictions
new_data$se_fit <- se_slope * new_data$log_HWI

# Compute 95% confidence intervals
new_data$upper <- new_data$fit + 1.96 * new_data$se_fit
new_data$lower <- new_data$fit - 1.96 * new_data$se_fit

# Scatterplot with PGLS regression line
ggplot(dispersal_analysis_partial, aes(x = log_HWI, y = Weibull_median_log)) +
  geom_point(color = "lightblue") +  # Scatterplot of data points
  geom_abline(intercept = intercept, 
              slope = slope, 
              color = "grey9", 
              linetype = ifelse(significant_slope, "solid", "dashed")) +  # Line type based on significance
   labs(title = "Relationship Between HWI and Dispersal Distance",
       subtitle = ifelse(significant_slope, "Significant Relationship", "No Significant Relationship"),
       x = "Log(HWI)", y = "Dispersal Distance (km)") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))


ggplot(dispersal_analysis_partial, aes(x = log_HWI, y = Weibull_median_log)) +
  geom_point(color = "lightblue") +  # Scatterplot of empirical data
  geom_line(data = new_data, aes(x = log_HWI, y = fit), 
            color = "grey9", 
            linetype = ifelse(significant_slope, "solid", "dashed")) +  # Regression line
  geom_ribbon(data = new_data, aes(x = log_HWI, ymin = lower, ymax = upper), 
              fill = "grey50", alpha = 0.2) +  # Confidence band
  labs(title = "Relationship Between HWI and Dispersal Distance",
       x = "Log(HWI)", y = "Log (Dispersal Distance[km])") +
  theme_classic() +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))



######### Loop over all predictors

library(ggplot2)

pgls_model <- modelo_pgls
data <- dispersal_analysis_partial

# Define the function
plot_pgls <- function(pgls_model, data) {
  
  # Extract coefficients and standard errors
  coef_summary <- summary(pgls_model)$coefficients
  intercept <- coef_summary["(Intercept)", "Estimate"]
  
  # Identify all predictors
  predictors <- rownames(coef_summary)[rownames(coef_summary) != "(Intercept)"]
  
  # Create an empty list to store plots
  plot_list <- list()
  
  # Loop through each predictor
  for (predictor in predictors) {
    # Extract slope and standard error for the current predictor
    slope <- coef_summary[predictor, "Estimate"]
    se_slope <- coef_summary[predictor, "Std. Error"]
    
    # Create prediction data
    pred_data <- data.frame(x = seq(min(data[[predictor]], na.rm = TRUE), 
                                    max(data[[predictor]], na.rm = TRUE), 
                                    length.out = 100))
    pred_data$fit <- intercept + slope * pred_data$x
    pred_data$se_fit <- se_slope * pred_data$x
    pred_data$upper <- pred_data$fit + 1.96 * pred_data$se_fit
    pred_data$lower <- pred_data$fit - 1.96 * pred_data$se_fit
    
    # Generate the plot
    p <- ggplot(data, aes_string(x = "body_mass", y = names(data)["Weibull_median_log"])) +  # Y-variable is the response
      geom_point(color = "lightblue") +  # Scatterplot of empirical data
      geom_line(data = pred_data, aes(x = x, y = fit), 
                color = "grey9", 
                linetype = ifelse(significant_slope, "solid", "dashed")) +  # Regression line
      geom_ribbon(data = pred_data, aes(x = x, ymin = lower, ymax = upper), 
                  fill = "grey50", alpha = 0.2) +  # Confidence band
      labs(title = "Relationship Between HWI and Dispersal Distance",
           x = "Log(HWI)", y = "Log (Dispersal Distance[km])") +
      theme_classic() +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
      
      
    p <- ggplot(data, aes(x = body_mass, y = Weibull_median_log)) +  # Y-variable is the response
      geom_point(color = "lightblue") +  # Scatterplot of empirical data
      geom_line(data = pred_data, aes(x = x, y = fit), 
                color = "grey9", 
                linetype = ifelse(significant_slope, "solid", "dashed")) +  # Regression line
      geom_ribbon(data = pred_data, aes(x = x, ymin = lower, ymax = upper), 
                  fill = "grey50", alpha = 0.2) +  # Confidence band
      labs(title = "Relationship Between HWI and Dispersal Distance",
           x = "Log(HWI)", y = "Log (Dispersal Distance[km])") +
      theme_classic() +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
    
    
       geom_point(color = "blue") +  # Scatterplot of data
      geom_line(data = pred_data, aes(x = x, y = fit), color = "red") +  # Regression line
      geom_ribbon(data = pred_data, aes(x = x, ymin = lower, ymax = upper), 
                  fill = "red", alpha = 0.2) +  # Confidence band
      labs(title = paste("Relationship Between", predictor, "and Response"),
           x = predictor, y = "Response") +
      theme_minimal() +
      theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
    
    # Add plot to the list
    plot_list[[predictor]] <- p
  }
  
  # Return the list of plots
  return(plot_list)
}

plot_list <- plot_pgls(modelo_pgls, dispersal_analysis_partial)

