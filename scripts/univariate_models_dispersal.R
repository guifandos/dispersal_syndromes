# Model evaluation and variable selection in a bayesian framework applying a horseshow prior. 

data <- dispersal_analysis_partial
A <- ape::vcv.phylo(dispersal_tree_partial)
phylo <- dispersal_tree_partial

#########################################
st <- c("log_HWI", "body_mass", "habita_for",  "diet", "distance_mig", "Latitude", "PC1")

st_full_log <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude", "body_mass:PC1", "body_mass:diet", "distance_mig:Latitude")

st_reduced<- c("log_HWI", "body_mass",  "diet", "distance_mig", "PC1")


# Parameters to fit univariate models for
parameters <- st  # Vector of parameter names to include individually

# Initialize a list to store models
univariate_models <- list()

# Loop through each parameter and fit the univariate model
for (param in parameters) {
  # Define the formula dynamically for the univariate model
  formula <- as.formula(
    paste0("Weibull_median ~ ", param, " + (1|gr(label, cov = A))")
  )
  
  # Fit the model
  univariate_models[[param]] <- brm(
    formula = formula,
    data = dispersal_analysis_partial,
    family = weibull(),
    data2 = list(A = A),
    prior = c(
      prior(normal(0, 2), class = "Intercept"),
      prior(horseshoe(), class = "b")
    ),
    chains = 2,
    cores = 2,
    iter = 4000,
    control = list(adapt_delta = 0.95)
  )
  
  # Optional: Print progress
  message("Fitted model for parameter: ", param)
}

# Check the results
univariate_models

# Include model assumptions

# Load necessary libraries
library(DHARMa)
# Load necessary libraries
library(DHARMa.helpers)

# Initialize a data frame to store slopes, significance, and DHARMa results
model_summary_data <- data.frame(
  Parameter = character(),
  Slope = numeric(),
  Lower_95_CI = numeric(),
  Upper_95_CI = numeric(),
  Significance = character(),
  Dispersion = numeric(),
  Uniformity = numeric(),
  ZeroInflation = numeric(),
  stringsAsFactors = FALSE
)

# Loop through each univariate model and extract coefficients, confidence intervals, and DHARMa diagnostics
for (param in parameters) {
  # Extract the summary of the fitted model
  summary_model <- summary(univariate_models[[param]])
  
  # Extract slope and confidence intervals
  slope <- summary_model$fixed[param, "Estimate"]
  lower_ci <- summary_model$fixed[param, "l-95% CI"]
  upper_ci <- summary_model$fixed[param, "u-95% CI"]
  
  # Determine if the slope is significant (if CI does not include 0)
  significance <- ifelse(lower_ci > 0 | upper_ci < 0, "Significant", "Not Significant")
  
  # Add the information to the summary data
  model_summary_data <- rbind(model_summary_data, data.frame(
    Parameter = param,
    Slope = slope,
    Lower_95_CI = lower_ci,
    Upper_95_CI = upper_ci,
    Significance = significance,
    Dispersion = NA,   # Placeholder for dispersion test
    Uniformity = NA,   # Placeholder for uniformity test
    ZeroInflation = NA # Placeholder for zero inflation test
  ))
  
  # DHARMa.helpers diagnostics
  simulated_residuals <- dh_check_brms(univariate_models[[param]])  # Using DHARMa.helpers for brms
  
  # Test for dispersion (check for overdispersion)
  dispersion <- testDispersion(simulated_residuals)$p.value
  
  # Test for uniformity (check for deviations from expected residual distribution)
  uniformity <- testUniformity(simulated_residuals)$p.value
  
  # Test for zero inflation (check if there are too many zeros in the residuals)
  zero_inflation <- testZeroInflation(simulated_residuals)$p.value
  
  # Add the DHARMa results to the model summary data
  model_summary_data[model_summary_data$Parameter == param, "Dispersion"] <- dispersion
  model_summary_data[model_summary_data$Parameter == param, "Uniformity"] <- uniformity
  model_summary_data[model_summary_data$Parameter == param, "ZeroInflation"] <- zero_inflation
}

# Print the results to see the summary of models with DHARMa diagnostics
print(model_summary_data)




# Export as a table

# Required library for exporting to Excel
library(writexl)

# Phylo signal 

hyp <- "sd_label__Intercept^2 / (sd_label__Intercept^2 + shape^2) = 0"
# Initialize a list to store model summaries
model_summaries <- list()

# Loop through each univariate model and extract the summary statistics
for (param in parameters) {
  # Extract the summary of the fitted model
  summary_model <- summary(univariate_models[[param]])
  
  # Extract relevant information (e.g., estimate, standard error, and confidence intervals)
  model_summary <- summary_model$fixed[, c("Estimate", "Est.Error", "l-95% CI", "u-95% CI")]
  
  # Remove intercept row (if present)
  model_summary <- model_summary[rownames(model_summary) != "Intercept", ]
  
  # Add the parameter name as a column for identification
  model_summary$Parameter <- param
  
  # Add phylogenetic signal hypothesis
 
  phylo_signal <- hypothesis(univariate_models[[param]], hyp, class = NULL)
  model_summary$Phylogenetic_Signal <- phylo_signal$hypothesis$Estimate[1]
  model_summary$Phylogenetic_Signal_CI.Lower <- phylo_signal$hypothesis$CI.Lower[1]
  model_summary$Phylogenetic_Signal_CI.Upper <- phylo_signal$hypothesis$CI.Upper[1]
  
  # Store the summary in the list
  model_summaries[[param]] <- model_summary
  
  # Optional: Print progress
  message("Extracted summary for parameter: ", param)
}

# Combine the summaries into a single data frame
univariate_results <- do.call(rbind, model_summaries)

# Reformat the table for better readability before exporting
univariate_results$`l-95% CI` <- round(univariate_results$`l-95% CI`, 2)
univariate_results$`u-95% CI` <- round(univariate_results$`u-95% CI`, 2)
univariate_results$Estimate <- round(univariate_results$Estimate, 2)
univariate_results$Est.Error <- round(univariate_results$Est.Error, 2)
univariate_results$Phylogenetic_Signal <- round(univariate_results$Phylogenetic_Signal, 2)
univariate_results$Phylogenetic_Signal_CI.Lower <- round(univariate_results$Phylogenetic_Signal_CI.Lower, 2)
univariate_results$Phylogenetic_Signal_CI.Upper <- round(univariate_results$Phylogenetic_Signal_CI.Upper, 2)
# Export to Excel
# Recode the parameters in the final table using the correct syntax
univariate_results$Parameter_Label <- dplyr::recode(univariate_results$Parameter,
                                             "log_HWI" = "Hand-wing index (log)", 
                                             "body_mass" = "Body mass (log)", 
                                             "habita_for" = "Habitat", 
                                             "diet" = "Diet", 
                                             "distance_mig" = "Migration distance", 
                                             "Latitude" = "Latitude"
)

# Print the updated model summary data with new parameter labels
print(univariate_results)


write_xlsx(univariate_results, path = "univariate_results.xlsx")



##### Plots of these univariate results

library(ggplot2)
library(tidybayes)
library(brms)
library(dplyr)

# Initialize a data frame to store slopes and significance information
model_summary_data <- data.frame(
  Parameter = character(),
  Slope = numeric(),
  Lower_95_CI = numeric(),
  Upper_95_CI = numeric(),
  Significance = character(),
  stringsAsFactors = FALSE
)

# Loop through each univariate model and extract coefficients and confidence intervals
for (param in parameters) {
  # Extract the summary of the fitted model
  summary_model <- summary(univariate_models[[param]])
  
  # Extract slope and confidence intervals
  slope <- summary_model$fixed[param, "Estimate"]
  lower_ci <- summary_model$fixed[param, "l-95% CI"]
  upper_ci <- summary_model$fixed[param, "u-95% CI"]
  
  # Determine if the slope is significant (if CI does not include 0)
  significance <- ifelse(lower_ci > 0 | upper_ci < 0, "Significant", "Not Significant")
  
  # Add the information to the summary data
  model_summary_data <- rbind(model_summary_data, data.frame(
    Parameter = param,
    Slope = slope,
    Lower_95_CI = lower_ci,
    Upper_95_CI = upper_ci,
    Significance = significance
  ))
}


# Recode the parameters in the final table using the correct syntax
model_summary_data$Parameter_Label <- dplyr::recode(model_summary_data$Parameter,
                                             "log_HWI" = "Hand-wing index (log)", 
                                             "body_mass" = "Body mass (log)", 
                                             "habita_for" = "Habitat", 
                                             "diet" = "Diet", 
                                             "distance_mig" = "Migration distance", 
                                             "Latitude" = "Latitude",
                                             "PC1" = "Life History"
)

# Print the updated model summary data with new parameter labels
print(model_summary_data)

# Add a column with parameter names for ggplot to use in facet_wrap
dispersal_analysis_partial_long <- dispersal_analysis_partial %>%
  gather(key = "Parameter", value = "Trait_value", all_of(parameters)) 

dispersal_analysis_partial_long$Parameter_Label <- dplyr::recode(dispersal_analysis_partial_long$Parameter,
                                                    "log_HWI" = "Hand-wing index (log)", 
                                                    "body_mass" = "Body mass (log)", 
                                                    "habita_for" = "Habitat", 
                                                    "diet" = "Diet", 
                                                    "distance_mig" = "Migration distance", 
                                                    "Latitude" = "Latitude",
                                                    "PC1" = "Life History"
)

# Plot the relationship between traits and dispersal distances for each parameter
ggplot(dispersal_analysis_partial_long, aes(x = Trait_value, y = Weibull_median )) +
  geom_point(color = "lightgrey") +  # Set point color to light grey
  # Add regression lines based on PGLS results
  geom_abline(data = model_summary_data, 
              aes(slope = Slope, intercept = 0, color = Significance),
              linetype = ifelse(model_summary_data$Significance == "Significant", "solid", "dashed")) +
    # Customize the plot with theme_tidybayes
  facet_wrap(~ Parameter_Label, scales = "free", ncol = 3) +  # Use free scales for both axes
  labs(x = "Trait Value",
       y = "Dispersal Distance (Weibull median)",
       title = "Relationship between Traits and Dispersal Distances",
       color = "Significance") +
  theme_tidybayes() +  # Apply the tidybayes theme
  theme(legend.position = "top")



# Merge significance data with the long-format dataset
dispersal_analysis_partial_long <- dispersal_analysis_partial %>%
  gather(key = "Parameter", value = "Trait_value", all_of(parameters)) %>%
  left_join(model_summary_data, by = "Parameter")  # Merge significance info

# Plot the relationship between traits and dispersal distances for each parameter
ggplot(dispersal_analysis_partial_long, aes(x = Trait_value, y = Weibull_median)) +
  geom_point(color = "lightgrey") +  # Set point color to light grey
  # Add regression lines with color and line type depending on significance
  geom_smooth(method = "lm", aes(color = Significance, linetype = ifelse(Significance == "Significant", "dashed", "solid")), 
              se = FALSE) +  # Remove confidence intervals
  # Customize the plot with facet_wrap
  facet_wrap(~ Parameter_Label, scales = "free", ncol = 3) +  # Use free scales for both axes
  labs(x = "Trait Value",
       y = "Dispersal Distance (Weibull median)",
       title = "Relationship between Traits and Dispersal Distances",
       color = "Significance",
       linetype = "Significance") +
  theme_tidybayes() +  # Apply the tidybayes theme
  theme(legend.position = "top")  + # Position legend at the top
  theme(legend.position = "none",  # Remove legend
      plot.title = element_blank())  # Remove title


ggsave(filename = "./results/weibull/figures/univariate_plots.tiff", width = 8, height = 5, device='tiff', dpi=700)
