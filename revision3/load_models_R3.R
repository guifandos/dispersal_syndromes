# MODEL COMBINATION AND ANALYSIS ============================================

library(tidyverse)
library(standardize)

# CONFIGURATION ==============================================================

# Directories
RESULTS_DIR <- "results/weibull"
OUTPUT_DIR <- "results/combined"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# LOAD INDIVIDUAL MODEL RESULTS =============================================

# Load reduced models
load(file.path(RESULTS_DIR, "average", "model_average_weibull_R3.RData"))
average_short <- results_short

load(file.path(RESULTS_DIR, "breeding", "model_breeding_weibull_R3.RData"))
breeding_short <- results_short

load(file.path(RESULTS_DIR, "natal", "model_natal_weibull_R3.RData"))
natal_short <- results_short

# Load complete models
load(file.path(RESULTS_DIR, "average", "model_average_complete_R3.RData"))
average_complete <- results_complete

load(file.path(RESULTS_DIR, "breeding", "model_breeding_complete_R3.RData"))
breeding_complete <- results_complete

load(file.path(RESULTS_DIR, "natal", "model_natal_complete_R3.RData"))
natal_complete <- results_complete

# COMBINE MODELS =============================================================

# Combine reduced models
all_short_models <- bind_rows(
  average_short %>% mutate(age = "average"),
  breeding_short %>% mutate(age = "breeding"),
  natal_short %>% mutate(age = "natal")
) %>%
  mutate(
    function_id = "weibull",
    dispersal_mode = type
  ) %>%
  select(-type)  # Remove redundant columns

# Combine complete models
all_complete_models <- bind_rows(
  average_complete %>% mutate(age = "average"),
  breeding_complete %>% mutate(age = "breeding"), 
  natal_complete %>% mutate(age = "natal")
) %>%
  mutate(
    function_id = "weibull",
    dispersal_mode = type
  ) %>%
  select(-type)  # Remove redundant columns

# SAVE COMBINED MODELS =======================================================

save(all_short_models, file = file.path(OUTPUT_DIR, "short_models_combined.RData"))
save(all_complete_models, file = file.path(OUTPUT_DIR, "complete_models_combined.RData"))

cat("Combined models saved:\n")
cat("- Short models:", file.path(OUTPUT_DIR, "short_models_combined.RData"), "\n")
cat("- Complete models:", file.path(OUTPUT_DIR, "complete_models_combined.RData"), "\n")


model_average_weibull <- all_short_models$model[[4]]
model_average_weibull <- all_short_models$model[[2]]
model_average_weibull <- all_complete_models$model[[2]]
model_average_weibull <- all_complete_models$model[[4]]
model_average_weibull <- all_short_models$model[[4]]


# VARIABLE SELECTION ANALYSIS ===============================================

extract_variable_selection <- function(models_data, response_type) {
  
  # Filter for specific response type
  models_subset <- models_data %>% filter(dispersal_mode == response_type)
  
  # Extract variable selection summaries
  var_selection_list <- list()
  
  for (i in 1:nrow(models_subset)) {
    
    # Get variable selection object
    var_sel <- models_subset$variable_selection[[i]]
    
    # Extract summary
    summary_result <- summary(var_sel, 
                              stats = "mlpd", 
                              type = c("mean", "lower", "upper"),
                              deltas = TRUE)
    
    # Convert to data frame and add metadata
    var_df <- as.data.frame(summary_result$perf_sub) %>%
      mutate(
        age = models_subset$age[i],
        function_id = models_subset$function_id[i],
        dispersal_mode = models_subset$dispersal_mode[i]
      ) %>%
      filter(ranking_fulldata != "(Intercept)")  # Remove intercept
    
    var_selection_list[[i]] <- var_df
  }
  
  # Combine all results
  combined_results <- bind_rows(var_selection_list)
  
  return(combined_results)
}

# Extract variable selection for median dispersal
median_var_selection <- extract_variable_selection(all_short_models, "median")

# Extract variable selection for long distance dispersal  
long_var_selection <- extract_variable_selection(all_short_models, "long")

# CLEAN VARIABLE NAMES =======================================================

create_clean_variable_names <- function(var_data) {
  var_data %>%
    mutate(
      ranking_fulldata = as.factor(ranking_fulldata),
      variable = recode_factor(
        ranking_fulldata,
        "PC1" = "Life history",
        "Latitude" = "Latitude", 
        "log_body_mass" = "Body mass",
        "body_mass" = "Body mass",
        "diet" = "Diet",
        "habita_for" = "Habitat",
        "log_HWI" = "HWI",
        "distance_mig" = "Distance migration",
        "log_body_mass:PC1" = "Body mass : Life history",
        "log_body_mass:diet" = "Body mass : Diet", 
        "log_body_mass:habita_for" = "Body mass : Habitat",
        "body_mass:PC1" = "Body mass : Life history",
        "body_mass:diet" = "Body mass : Diet",
        "body_mass:habita_for" = "Body mass : Habitat",
        "distance_mig:Latitude" = "Distance migration : Latitude"
      )
    )
}

# Apply clean names
median_var_clean <- create_clean_variable_names(median_var_selection)
long_var_clean <- create_clean_variable_names(long_var_selection)

# SCALE MLPD VALUES ==========================================================

scale_mlpd_by_age <- function(var_data) {
  var_data %>%
    mutate(
      mlpd_scaled = scale_by(mlpd ~ age) %>% as.vector()
    )
}

# Apply scaling
median_var_scaled <- scale_mlpd_by_age(median_var_clean)
long_var_scaled <- scale_mlpd_by_age(long_var_clean)

# SAVE VARIABLE SELECTION RESULTS ===========================================

# Save raw results
write_csv(median_var_clean, file.path(OUTPUT_DIR, "median_variable_selection.csv"))
write_csv(long_var_clean, file.path(OUTPUT_DIR, "long_variable_selection.csv"))

# Save scaled results
write_csv(median_var_scaled, file.path(OUTPUT_DIR, "median_variable_selection_scaled.csv"))
write_csv(long_var_scaled, file.path(OUTPUT_DIR, "long_variable_selection_scaled.csv"))

# SUMMARY STATISTICS =========================================================

# Create summary tables
create_selection_summary <- function(var_data, title) {
  
  cat("\n", title, "\n")
  cat(rep("=", nchar(title)), "\n")
  
  # Count by age and variable
  summary_table <- var_data %>%
    group_by(age, variable) %>%
    summarise(
      mean_mlpd = mean(mlpd, na.rm = TRUE),
      min_rank = min(as.numeric(as.character(ranking_fulldata)), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(age, desc(mean_mlpd))
  
  print(summary_table)
  
  # Top variables by age
  top_vars <- var_data %>%
    group_by(age) %>%
    slice_max(mlpd, n = 3) %>%
    select(age, variable, mlpd, ranking_fulldata)
  
  cat("\nTop 3 variables by age:\n")
  print(top_vars)
  
  return(summary_table)
}

# Generate summaries
median_summary <- create_selection_summary(median_var_clean, "MEDIAN DISPERSAL VARIABLE SELECTION")
long_summary <- create_selection_summary(long_var_clean, "LONG DISTANCE DISPERSAL VARIABLE SELECTION")

# DIAGNOSTIC PLOTS ===========================================================

create_selection_plot <- function(var_data, title, filename) {
  
  p <- var_data %>%
    ggplot(aes(x = reorder(variable, mlpd), y = mlpd, fill = age)) +
    geom_col(position = "dodge", alpha = 0.8) +
    coord_flip() +
    labs(
      title = title,
      x = "Variable",
      y = "MLPD (Mean Log Predictive Density)",
      fill = "Dispersal Age"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom"
    )
  
  ggsave(file.path(OUTPUT_DIR, filename), p, width = 10, height = 8)
  
  return(p)
}

# Create plots
median_plot <- create_selection_plot(
  median_var_clean, 
  "Variable Selection: Median Dispersal", 
  "median_variable_selection_plot.png"
)

long_plot <- create_selection_plot(
  long_var_clean, 
  "Variable Selection: Long Distance Dispersal", 
  "long_variable_selection_plot.png"
)

# Display plots
print(median_plot)
print(long_plot)

# FINAL SUMMARY ===============================================================

cat("\n", rep("=", 60), "\n")
cat("ANALYSIS COMPLETED\n")
cat(rep("=", 60), "\n")

cat("\nFiles created:\n")
cat("Models:\n")
cat("- Short models combined:", file.path(OUTPUT_DIR, "short_models_combined.RData"), "\n")
cat("- Complete models combined:", file.path(OUTPUT_DIR, "complete_models_combined.RData"), "\n")

cat("\nVariable Selection:\n")
cat("- Median selection (raw):", file.path(OUTPUT_DIR, "median_variable_selection.csv"), "\n")
cat("- Median selection (scaled):", file.path(OUTPUT_DIR, "median_variable_selection_scaled.csv"), "\n")
cat("- Long selection (raw):", file.path(OUTPUT_DIR, "long_variable_selection.csv"), "\n")
cat("- Long selection (scaled):", file.path(OUTPUT_DIR, "long_variable_selection_scaled.csv"), "\n")

cat("\nPlots:\n")
cat("- Median selection plot:", file.path(OUTPUT_DIR, "median_variable_selection_plot.png"), "\n")
cat("- Long selection plot:", file.path(OUTPUT_DIR, "long_variable_selection_plot.png"), "\n")

cat("\nDataset dimensions:\n")
cat("- Combined short models:", nrow(all_short_models), "rows,", ncol(all_short_models), "columns\n")
cat("- Combined complete models:", nrow(all_complete_models), "rows,", ncol(all_complete_models), "columns\n")
cat("- Median variable selection:", nrow(median_var_clean), "observations\n")
cat("- Long variable selection:", nrow(long_var_clean), "observations\n")

