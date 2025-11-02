# ==============================================================================
# BAYESIAN STACKING FOR COMPETITIVE MODELS - REDUCED MULTIVARIATE VERSION
# ==============================================================================
#
# DESCRIPTION:
# Implements Bayesian stacking using reduced multivariate models and univariate
# models. Uses the comprehensive model comparison table to identify competitive
# models (Δ LOOIC ≤ 2).
#
# KEY DIFFERENCES BETWEEN MODEL TYPES:
# 
# 1. UNIVARIATE MODELS:
#    - Single predictor + phylogeny
#    - No interactions
#    - Independent model for each variable
#
# 2. REDUCED MULTIVARIATE MODELS:
#    - Selected predictors via projection predictive selection
#    - May include interactions between predictors
#    - Account for multicollinearity and predictor relationships
#    - Can contain 1-7 predictors depending on variable selection
#
# 3. SPECIAL CASE - "Effectively Univariate" Reduced Models:
#    - When variable selection chooses only ONE predictor
#    - Still called "reduced_model" but functionally similar to univariate
#    - Difference: fitted within multivariate framework (different priors/setup)
#    - May have slightly different coefficients due to model fitting approach
#
# INPUTS REQUIRED:
# - Fitted univariate models
# - Fitted reduced multivariate models
# - results/comparison_R3/comprehensive_model_comparison.csv
#
# OUTPUTS:
# - Stacking weights for each competitive model group
# - Weighted predictions and coefficients
# - Summary tables with model weights and interpretations
#
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(loo)
  library(brms)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(stringr)
})

# LOAD REQUIRED DATA ========================================================

cat("=== LOADING COMPETITIVE MODELS DATA ===\n")

# Load comprehensive model comparison table
if (!file.exists("results/comparison_R3/comprehensive_model_comparison.csv")) {
  stop("Please ensure comprehensive_model_comparison.csv exists in results/comparison_R3/")
}

competitive_models <- read.csv("results/comparison_R3/comprehensive_model_comparison.csv")
cat("Total models loaded:", nrow(competitive_models), "models\n")

# Identify competitive models (delta_looic <= 2)
competitive_only <- competitive_models %>%
  filter(delta_looic <= 2) %>%
  arrange(dispersal_type, distance_type, looic_rank)

cat("Competitive models (Δ LOOIC ≤ 2):", nrow(competitive_only), "models\n")

# Identify which reduced models are effectively univariate
reduced_model_summary <- competitive_models %>%
  filter(model_type == "multivariate_reduced") %>%
  select(dispersal_type, distance_type, variable, model_id, n_params) %>%
  mutate(
    effectively_univariate = n_params <= 3,  # Intercept + 1 predictor + phylogeny
    model_complexity = case_when(
      n_params <= 3 ~ "Effectively univariate (1 predictor)",
      n_params <= 5 ~ "Simple multivariate (2-3 predictors)", 
      n_params <= 7 ~ "Moderate multivariate (4-5 predictors)",
      TRUE ~ "Complex multivariate (6+ predictors)"
    )
  )

cat("\nREDUCED MODEL COMPLEXITY SUMMARY:\n")
print(table(reduced_model_summary$model_complexity))

# Show which reduced models are effectively univariate
effectively_univariate <- reduced_model_summary %>%
  filter(effectively_univariate) %>%
  arrange(dispersal_type, distance_type)

if (nrow(effectively_univariate) > 0) {
  cat("\nEFFECTIVELY UNIVARIATE REDUCED MODELS:\n")
  for (i in 1:nrow(effectively_univariate)) {
    cat(sprintf("• %s %s: %s (n_params = %d)\n",
                str_to_title(effectively_univariate$dispersal_type[i]),
                effectively_univariate$distance_type[i], 
                effectively_univariate$variable[i],
                effectively_univariate$n_params[i]))
  }
  cat("\nNote: These differ from true univariate models due to:\n")
  cat("- Different fitting framework (multivariate vs univariate setup)\n")
  cat("- Different prior specifications\n") 
  cat("- Variable selection process vs direct fitting\n")
}

# LOAD MODEL OBJECTS ========================================================

cat("\n=== LOADING MODEL OBJECTS ===\n")

# Function to load models by dispersal type
load_models_by_type <- function(disp_type) {
  
  distance_types <- c("long", "median")
  
  # Load univariate models
  univariate_path <- paste0("results/weibull/", disp_type, "/univariate_models.RData")
  if (file.exists(univariate_path)) {
    load(univariate_path)
    univariate_models <- models
    cat("✓ Loaded", length(univariate_models), "univariate models for", disp_type, "\n")
  } else {
    stop("Univariate models not found for ", disp_type)
  }
  
  # Load reduced models
  reduced_path <- "results/combined/short_models_combined.RData" 
  if (file.exists(reduced_path) && !exists("all_short_models")) {
    load(reduced_path)
    reduced_models <- all_short_models
    cat("✓ Loaded reduced multivariate models\n")
  } else if (exists("all_short_models")) {
    reduced_models <- all_short_models
  } else {
    stop("Reduced models not found")
  }
  
  return(list(
    univariate = univariate_models,
    reduced = reduced_models
  ))
}

# Load all models
models_average <- load_models_by_type("average")
models_natal <- load_models_by_type("natal") 
models_breeding <- load_models_by_type("breeding")

# FUNCTION TO GET MODEL OBJECT BY ID ========================================

get_model_by_id <- function(dispersal_type, distance_type, variable, model_type) {
  
  if (model_type %in% c("multivariate_reduced", "reduced_model")) {
    # Reduced multivariate model
    if (dispersal_type == "average") {
      model_row <- models_average$reduced %>%
        filter(age == dispersal_type, dispersal_mode == distance_type)
    } else if (dispersal_type == "natal") {
      model_row <- models_natal$reduced %>%
        filter(age == dispersal_type, dispersal_mode == distance_type)  
    } else if (dispersal_type == "breeding") {
      model_row <- models_breeding$reduced %>%
        filter(age == dispersal_type, dispersal_mode == distance_type)
    }
    
    if (nrow(model_row) > 0) {
      return(model_row$model[[1]])
    }
    
  } else if (model_type == "univariate") {
    # Univariate model
    model_name <- paste0(variable, "_", distance_type)
    
    if (dispersal_type == "average") {
      return(models_average$univariate[[model_name]])
    } else if (dispersal_type == "natal") {
      return(models_natal$univariate[[model_name]])
    } else if (dispersal_type == "breeding") {
      return(models_breeding$univariate[[model_name]])
    }
  }
  
  return(NULL)
}

# FUNCTION TO CALCULATE STACKING WEIGHTS ====================================

calculate_stacking_weights <- function(group_data) {
  
  cat("Processing group:", unique(group_data$dispersal_type), 
      unique(group_data$distance_type), "\n")
  
  # Get model objects
  models_list <- list()
  loo_list <- list()
  
  for (i in 1:nrow(group_data)) {
    dispersal_type <- group_data$dispersal_type[i]
    distance_type <- group_data$distance_type[i]
    variable <- group_data$variable[i]
    model_type <- group_data$model_type[i]
    
    # Create unique model identifier
    if (model_type %in% c("multivariate_reduced", "reduced_model")) {
      model_id <- paste0(dispersal_type, "_", distance_type, "_reduced")
      
      # Check if this is effectively univariate
      n_params <- group_data$n_params[i]
      if (n_params <= 3) {
        model_id <- paste0(model_id, "_effectively_univariate")
      }
      
    } else {
      model_id <- paste0(dispersal_type, "_", distance_type, "_", variable, "_univariate")
    }
    
    model_obj <- get_model_by_id(dispersal_type, distance_type, variable, model_type)
    
    if (!is.null(model_obj)) {
      models_list[[model_id]] <- model_obj
      
      # Calculate LOO for this model using existing LOOIC values
      tryCatch({
        # Use pre-calculated LOOIC from comparison table
        looic_val <- group_data$looic[i]
        elpd_loo_val <- group_data$elpd_loo[i]
        se_looic_val <- group_data$se_looic[i]
        p_loo_val <- group_data$p_loo[i]
        
        # Create LOO object structure (simplified for stacking)
        # Note: For proper stacking, we'd need full LOO calculation
        # Here we use a workaround with available metrics
        
        cat("  ✓ Using pre-calculated LOO for", model_id, "(LOOIC =", round(looic_val, 2), ")\n")
        
        # Calculate LOO properly for stacking
        loo_obj <- loo(model_obj, refresh = 0)
        loo_list[[model_id]] <- loo_obj
        
      }, error = function(e) {
        cat("  ✗ Error calculating LOO for", model_id, ":", e$message, "\n")
      })
    } else {
      cat("  ⚠ Could not find model for", dispersal_type, distance_type, variable, model_type, "\n")
    }
  }
  
  # Calculate stacking weights if we have multiple valid models
  if (length(loo_list) > 1) {
    tryCatch({
      weights <- loo_model_weights(loo_list, method = "stacking")
      
      # Create results data frame
      results <- data.frame(
        model_id = names(weights),
        stacking_weight = as.numeric(weights),
        stringsAsFactors = FALSE
      )
      
      # Extract model information from model_id  
      results <- results %>%
        mutate(
          dispersal_type = sapply(strsplit(model_id, "_"), `[`, 1),
          distance_type = sapply(strsplit(model_id, "_"), `[`, 2),
          model_type = case_when(
            grepl("reduced", model_id) & grepl("effectively_univariate", model_id) ~ "reduced_effectively_univariate",
            grepl("reduced", model_id) ~ "multivariate_reduced", 
            TRUE ~ "univariate"
          ),
          variable = case_when(
            grepl("reduced", model_id) ~ "reduced_model",
            TRUE ~ sapply(strsplit(model_id, "_"), function(x) paste(x[3:(length(x)-1)], collapse = "_"))
          ),
          effectively_univariate = grepl("effectively_univariate", model_id)
        ) %>%
        arrange(desc(stacking_weight))
      
      cat("  ✓ Stacking weights calculated for", length(weights), "models\n")
      return(results)
      
    }, error = function(e) {
      cat("  ✗ Error in stacking calculation:", e$message, "\n")
      return(NULL)
    })
  } else {
    cat("  ⚠ Not enough valid models for stacking (need >1, have", length(loo_list), ")\n")
    return(NULL)
  }
}

# MAIN STACKING ANALYSIS ====================================================

cat("\n=== PERFORMING STACKING ANALYSIS ===\n")

# Group competitive models by dispersal_type and distance_type
competitive_groups <- competitive_only %>%
  group_by(dispersal_type, distance_type) %>%
  filter(n() > 1) %>%  # Only groups with multiple competitive models
  group_split()

cat("Found", length(competitive_groups), "groups with multiple competitive models\n")

# Show summary of what we're analyzing
cat("\nCOMPETITIVE GROUPS SUMMARY:\n")
for (i in seq_along(competitive_groups)) {
  group <- competitive_groups[[i]]
  group_name <- paste(unique(group$dispersal_type), unique(group$distance_type), sep = "_")
  
  n_univariate <- sum(group$model_type == "univariate")
  n_reduced <- sum(group$model_type == "multivariate_reduced")
  
  # Check if reduced model is effectively univariate
  reduced_n_params <- group$n_params[group$model_type == "multivariate_reduced"]
  effectively_univariate <- any(reduced_n_params <= 3)
  
  cat(sprintf("• %s: %d univariate + %d reduced%s (n_competitive = %d)\n",
              str_to_title(group_name),
              n_univariate,
              n_reduced,
              ifelse(effectively_univariate, " [effectively univariate]", ""),
              nrow(group)))
}

# Calculate stacking weights for each group
stacking_results <- list()

for (i in seq_along(competitive_groups)) {
  group <- competitive_groups[[i]]
  group_name <- paste(unique(group$dispersal_type), unique(group$distance_type), sep = "_")
  
  weights <- calculate_stacking_weights(group)
  
  if (!is.null(weights)) {
    stacking_results[[group_name]] <- weights
  }
}

cat("Successfully calculated stacking weights for", length(stacking_results), "groups\n")

# CREATE SUMMARY TABLES =====================================================

cat("\n=== CREATING STACKING SUMMARY TABLES ===\n")

# Combine all stacking results
all_stacking_weights <- do.call(rbind, lapply(names(stacking_results), function(name) {
  parts <- strsplit(name, "_")[[1]]
  dispersal_type <- parts[1]
  distance_type <- parts[2]
  
  result <- stacking_results[[name]]
  result$dispersal_type <- dispersal_type
  result$distance_type <- distance_type
  result$group_name <- name
  
  return(result)
}))

# Save detailed stacking weights
if (!dir.exists("results/tables")) dir.create("results/tables", recursive = TRUE)
write.csv(all_stacking_weights, "results/tables/stacking_weights_detailed_reduced.csv", row.names = FALSE)

# Create summary table
stacking_summary <- all_stacking_weights %>%
  group_by(dispersal_type, distance_type) %>%
  summarise(
    n_competitive_models = n(),
    best_model = variable[which.max(stacking_weight)],
    best_model_type = model_type[which.max(stacking_weight)],
    best_weight = max(stacking_weight),
    weight_concentration = max(stacking_weight),
    effective_models = sum(stacking_weight > 0.1),
    model_uncertainty = 1 - max(stacking_weight),
    has_effectively_univariate = any(effectively_univariate),
    n_reduced_models = sum(grepl("reduced", model_type)),
    n_univariate_models = sum(model_type == "univariate"),
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = case_when(
      weight_concentration > 0.8 ~ "Strong preference for best model",
      weight_concentration > 0.6 ~ "Moderate preference for best model", 
      weight_concentration > 0.4 ~ "Substantial model uncertainty",
      TRUE ~ "High model uncertainty"
    ),
    recommendation = case_when(
      model_uncertainty <= 0.2 ~ "Use best model with confidence",
      model_uncertainty <= 0.4 ~ "Consider best model with some caution",
      TRUE ~ "Use model averaging"
    )
  )

write.csv(stacking_summary, "results/tables/stacking_summary_reduced.csv", row.names = FALSE)

# FUNCTION TO CREATE WEIGHTED PREDICTIONS ===================================

create_weighted_predictions <- function(group_data, weights_data, newdata = NULL) {
  
  if (is.null(weights_data) || nrow(weights_data) == 0) {
    return(NULL)
  }
  
  # Use original data if newdata not provided
  if (is.null(newdata)) {
    dispersal_type <- unique(weights_data$dispersal_type)
    data_path <- paste0("data/processed/dispersal_", dispersal_type, "_complete.RData")
    
    if (file.exists(data_path)) {
      load(data_path)
      if (dispersal_type == "average") newdata <- data_average
      else if (dispersal_type == "natal") newdata <- data_natal
      else if (dispersal_type == "breeding") newdata <- data_breeding
    } else {
      cat("  Warning: Data file not found for", dispersal_type, "\n")
      return(NULL)
    }
  }
  
  weighted_pred <- NULL
  total_weight <- 0
  
  for (i in 1:nrow(weights_data)) {
    dispersal_type <- weights_data$dispersal_type[i]
    distance_type <- weights_data$distance_type[i]
    variable <- weights_data$variable[i]
    model_type_original <- ifelse(grepl("reduced", weights_data$model_type[i]), 
                                  "multivariate_reduced", "univariate")
    weight <- weights_data$stacking_weight[i]
    
    model_obj <- get_model_by_id(dispersal_type, distance_type, variable, model_type_original)
    
    if (!is.null(model_obj) && weight > 0.01) {
      tryCatch({
        pred <- predict(model_obj, newdata = newdata, re_formula = NA)
        
        if (is.null(weighted_pred)) {
          weighted_pred <- weight * pred[, "Estimate"]
        } else {
          weighted_pred <- weighted_pred + weight * pred[, "Estimate"]
        }
        
        total_weight <- total_weight + weight
        
      }, error = function(e) {
        cat("  Error in prediction for", dispersal_type, distance_type, variable, ":", e$message, "\n")
      })
    }
  }
  
  if (!is.null(weighted_pred) && total_weight > 0) {
    weighted_pred <- weighted_pred / total_weight
    return(weighted_pred)
  }
  
  return(NULL)
}

# GENERATE WEIGHTED PREDICTIONS =============================================

cat("\n=== GENERATING WEIGHTED PREDICTIONS ===\n")

prediction_results <- list()

for (i in seq_along(competitive_groups)) {
  group <- competitive_groups[[i]]
  group_name <- paste(unique(group$dispersal_type), unique(group$distance_type), sep = "_")
  
  if (group_name %in% names(stacking_results)) {
    weights_data <- stacking_results[[group_name]]
    
    cat("Creating weighted predictions for", group_name, "\n")
    weighted_pred <- create_weighted_predictions(group, weights_data)
    
    if (!is.null(weighted_pred)) {
      prediction_results[[group_name]] <- weighted_pred
      cat("  ✓ Weighted predictions created\n")
    }
  }
}

# CALCULATE WEIGHTED COEFFICIENTS ============================================

cat("\n=== CALCULATING WEIGHTED COEFFICIENTS ===\n")

calculate_weighted_coefficients <- function(group_data, weights_data) {
  
  if (is.null(weights_data) || nrow(weights_data) == 0) {
    return(NULL)
  }
  
  weighted_coefs <- NULL
  total_weight <- 0
  coefficient_names <- NULL
  
  for (i in 1:nrow(weights_data)) {
    dispersal_type <- weights_data$dispersal_type[i]
    distance_type <- weights_data$distance_type[i]
    variable <- weights_data$variable[i]
    model_type_original <- ifelse(grepl("reduced", weights_data$model_type[i]), 
                                  "multivariate_reduced", "univariate")
    weight <- weights_data$stacking_weight[i]
    
    model_obj <- get_model_by_id(dispersal_type, distance_type, variable, model_type_original)
    
    if (!is.null(model_obj) && weight > 0.01) {
      tryCatch({
        coefs <- fixef(model_obj)[, "Estimate"]
        coef_names <- names(coefs)
        
        if (is.null(weighted_coefs)) {
          weighted_coefs <- weight * coefs
          coefficient_names <- coef_names
        } else {
          # Align coefficients 
          aligned_coefs <- rep(0, length(coefficient_names))
          names(aligned_coefs) <- coefficient_names
          
          for (name in coef_names) {
            if (name %in% coefficient_names) {
              aligned_coefs[name] <- coefs[name]
            } else {
              # Add new coefficient
              coefficient_names <- c(coefficient_names, name)
              weighted_coefs <- c(weighted_coefs, 0)
              names(weighted_coefs) <- coefficient_names
              aligned_coefs <- c(aligned_coefs, coefs[name])
              names(aligned_coefs) <- coefficient_names
            }
          }
          
          weighted_coefs <- weighted_coefs + weight * aligned_coefs
        }
        
        total_weight <- total_weight + weight
        
      }, error = function(e) {
        cat("  Error extracting coefficients for", dispersal_type, distance_type, variable, ":", e$message, "\n")
      })
    }
  }
  
  if (!is.null(weighted_coefs) && total_weight > 0) {
    weighted_coefs <- weighted_coefs / total_weight
    
    results <- data.frame(
      coefficient = names(weighted_coefs),
      weighted_estimate = as.numeric(weighted_coefs),
      dispersal_type = unique(weights_data$dispersal_type),
      distance_type = unique(weights_data$distance_type),
      stringsAsFactors = FALSE
    )
    
    return(results)
  }
  
  return(NULL)
}

all_weighted_coefficients <- list()

for (i in seq_along(competitive_groups)) {
  group <- competitive_groups[[i]]
  group_name <- paste(unique(group$dispersal_type), unique(group$distance_type), sep = "_")
  
  if (group_name %in% names(stacking_results)) {
    weights_data <- stacking_results[[group_name]]
    
    cat("Calculating weighted coefficients for", group_name, "\n")
    weighted_coefs <- calculate_weighted_coefficients(group, weights_data)
    
    if (!is.null(weighted_coefs)) {
      all_weighted_coefficients[[group_name]] <- weighted_coefs
      cat("  ✓ Weighted coefficients calculated\n")
    }
  }
}

if (length(all_weighted_coefficients) > 0) {
  combined_weighted_coefficients <- do.call(rbind, all_weighted_coefficients)
  rownames(combined_weighted_coefficients) <- NULL
  
  write.csv(combined_weighted_coefficients, 
            "results/tables/weighted_coefficients_reduced.csv", 
            row.names = FALSE)
  
  cat("Weighted coefficients calculated for", length(all_weighted_coefficients), "groups\n")
}

# CREATE VISUALIZATIONS ======================================================

cat("\n=== CREATING STACKING VISUALIZATIONS ===\n")

if (!dir.exists("results/figures")) dir.create("results/figures", recursive = TRUE)

if (nrow(all_stacking_weights) > 0) {
  
  # 1. STACKING WEIGHTS PLOT =================================================
  
  plot_data <- all_stacking_weights %>%
    mutate(
      group_label = paste(str_to_title(dispersal_type), distance_type),
      model_label = case_when(
        model_type == "multivariate_reduced" ~ "Reduced Multivariate",
        model_type == "reduced_effectively_univariate" ~ "Reduced (Effectively Univariate)",
        model_type == "univariate" ~ paste("Univariate:", variable),
        TRUE ~ variable
      ),
      model_category = case_when(
        grepl("reduced", model_type) ~ "Reduced Multivariate",
        model_type == "univariate" ~ "Univariate",
        TRUE ~ "Other"
      )
    )
  
  p1 <- ggplot(plot_data, aes(x = reorder(model_label, stacking_weight), 
                              y = stacking_weight, 
                              fill = model_category)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_text(aes(label = paste0(round(stacking_weight * 100, 1), "%")), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    facet_wrap(~group_label, scales = "free_y") +
    scale_fill_manual(values = c("Reduced Multivariate" = "#2E86AB", 
                                 "Univariate" = "#A23B72",
                                 "Other" = "#F18F01")) +
    labs(
      title = "Bayesian Stacking Weights: Reduced vs Univariate Models",
      subtitle = "Higher weights indicate better predictive performance",
      x = "Model",
      y = "Stacking Weight",
      fill = "Model Type",
      caption = "Reduced models selected via projection predictive inference"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      strip.text = element_text(face = "bold", color = "black"),
      strip.background = element_rect(fill = "white", color = "black"),
      legend.position = "bottom",
      plot.caption = element_text(hjust = 0, color = "grey60"),
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.25),
      plot.title = element_text(color = "black"),
      plot.subtitle = element_text(color = "black")
    ) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1))
  
  ggsave("results/figures/stacking_weights_reduced.png", p1, 
         width = 12, height = 8, dpi = 300)
  
  cat("✓ Stacking weights plot saved\n")
  
  # 2. MODEL UNCERTAINTY PLOT ================================================
  
  uncertainty_data <- stacking_summary %>%
    mutate(
      group_label = paste(str_to_title(dispersal_type), distance_type),
      uncertainty_level = case_when(
        model_uncertainty <= 0.2 ~ "Low (≤20%)",
        model_uncertainty <= 0.4 ~ "Moderate (20-40%)",
        model_uncertainty <= 0.6 ~ "High (40-60%)",
        TRUE ~ "Very High (>60%)"
      ),
      uncertainty_level = factor(uncertainty_level, 
                                 levels = c("Low (≤20%)", "Moderate (20-40%)", 
                                            "High (40-60%)", "Very High (>60%)")),
      model_composition = paste0("R:", n_reduced_models, " U:", n_univariate_models)
    )
  
  p2 <- ggplot(uncertainty_data, aes(x = reorder(group_label, model_uncertainty), 
                                     y = model_uncertainty, 
                                     fill = uncertainty_level)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", alpha = 0.7, size = 0.8) +
    geom_text(aes(label = paste0(round(model_uncertainty * 100, 1), "%\n(", model_composition, ")")), 
              hjust = -0.1, size = 2.5, color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Low (≤20%)" = "#2E8B57", 
                                 "Moderate (20-40%)" = "#FFD700", 
                                 "High (40-60%)" = "#FF8C00", 
                                 "Very High (>60%)" = "#DC143C")) +
    labs(
      title = "Model Uncertainty: Reduced vs Univariate Competition",
      subtitle = "Red line: averaging threshold (40%). R = Reduced, U = Univariate",
      x = "Dispersal Group",
      y = "Model Uncertainty (1 - Best Weight)",
      fill = "Uncertainty Level"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.position = "bottom",
      legend.text = element_text(color = "black"),
      legend.title = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.25),
      plot.title = element_text(color = "black"),
      plot.subtitle = element_text(color = "black")
    ) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1))
  
  ggsave("results/figures/model_uncertainty_reduced.png", p2, 
         width = 10, height = 6, dpi = 300)
  
  cat("✓ Model uncertainty plot saved\n")
  
  # 3. WEIGHTED COEFFICIENTS PLOT ============================================
  
  if (exists("combined_weighted_coefficients") && nrow(combined_weighted_coefficients) > 0) {
    
    coef_plot_data <- combined_weighted_coefficients %>%
      filter(!grepl("Intercept", coefficient)) %>%
      mutate(
        group_label = paste(str_to_title(dispersal_type), distance_type),
        coefficient_clean = gsub(":", " × ", coefficient),
        coefficient_clean = gsub("_", " ", coefficient_clean),
        significant = abs(weighted_estimate) > 0.1
      )
    
    p3 <- ggplot(coef_plot_data, aes(x = weighted_estimate, 
                                     y = reorder(coefficient_clean, weighted_estimate),
                                     color = significant)) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "black", size = 0.8) +
      geom_point(size = 3, alpha = 0.8) +
      facet_wrap(~group_label, scales = "free") +
      scale_color_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#CCCCCC")) +
      labs(
        title = "Weighted Coefficients from Reduced vs Univariate Stacking",
        subtitle = "Coefficients averaged by stacking weights across competitive models",
        x = "Weighted Coefficient Estimate",
        y = "Predictor Variable",
        color = "Notable Effect\n(|coef| > 0.1)"
      ) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        strip.text = element_text(face = "bold", color = "black"),
        strip.background = element_rect(fill = "white", color = "black"),
        legend.position = "bottom",
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        panel.grid.major = element_line(color = "grey90", size = 0.5),
        panel.grid.minor = element_line(color = "grey95", size = 0.25),
        plot.title = element_text(color = "black"),
        plot.subtitle = element_text(color = "black")
      )
    
    ggsave("results/figures/weighted_coefficients_reduced.png", p3, 
           width = 12, height = 8, dpi = 300)
    
    cat("✓ Weighted coefficients plot saved\n")
    
  } else {
    cat("⚠ No weighted coefficients available for plotting\n")
  }
  
  # 4. ADDITIONAL PLOT: COMPARISON HEATMAP ===================================
  
  # Create a heatmap comparing reduced vs univariate performance
  if (exists("performance_comparison") && nrow(performance_comparison) > 0) {
    
    heatmap_data <- performance_comparison %>%
      select(dispersal_type, distance_type, weight_difference) %>%
      mutate(
        group_label = paste(str_to_title(dispersal_type), distance_type),
        advantage_category = case_when(
          weight_difference > 0.2 ~ "Strong Reduced Advantage",
          weight_difference > 0.05 ~ "Moderate Reduced Advantage", 
          weight_difference > -0.05 ~ "Similar Performance",
          weight_difference > -0.2 ~ "Moderate Univariate Advantage",
          TRUE ~ "Strong Univariate Advantage"
        ),
        advantage_category = factor(advantage_category, levels = c(
          "Strong Univariate Advantage", "Moderate Univariate Advantage",
          "Similar Performance", "Moderate Reduced Advantage", "Strong Reduced Advantage"
        ))
      )
    
    p4 <- ggplot(heatmap_data, aes(x = distance_type, y = dispersal_type, 
                                   fill = weight_difference)) +
      geom_tile(color = "white", size = 1) +
      geom_text(aes(label = paste0(ifelse(weight_difference > 0, "+", ""), 
                                   round(weight_difference * 100, 1), "%")),
                color = "white", fontface = "bold", size = 4) +
      scale_fill_gradient2(low = "#A23B72", mid = "white", high = "#2E86AB",
                           midpoint = 0, 
                           name = "Weight\nDifference",
                           labels = scales::percent_format()) +
      labs(
        title = "Reduced vs Univariate Model Performance Comparison",
        subtitle = "Positive values favor reduced models, negative favor univariate",
        x = "Distance Type",
        y = "Dispersal Type"
      ) +
      theme_bw() +
      theme(
        axis.text = element_text(color = "black"),
        axis.title = element_text(color = "black"),
        legend.position = "right",
        legend.text = element_text(color = "black"),
        legend.title = element_text(color = "black"),
        plot.title = element_text(color = "black"),
        plot.subtitle = element_text(color = "black"),
        panel.grid = element_blank()
      ) +
      coord_equal()
    
    ggsave("results/figures/performance_comparison_heatmap.png", p4, 
           width = 8, height = 6, dpi = 300)
    
    cat("✓ Performance comparison heatmap saved\n")
  }
  
  # 5. MODEL COMPLEXITY DISTRIBUTION =========================================
  
  # Show distribution of model complexities
  complexity_data <- all_stacking_weights %>%
    left_join(
      competitive_models %>% select(dispersal_type, distance_type, model_type, variable, n_params),
      by = c("dispersal_type", "distance_type", "model_type", "variable")
    ) %>%
    mutate(
      complexity_category = case_when(
        model_type == "univariate" ~ "Univariate",
        n_params <= 3 ~ "Effectively Univariate",
        n_params <= 5 ~ "Simple Multivariate",
        n_params <= 7 ~ "Moderate Multivariate", 
        TRUE ~ "Complex Multivariate"
      ),
      complexity_category = factor(complexity_category, levels = c(
        "Univariate", "Effectively Univariate", "Simple Multivariate",
        "Moderate Multivariate", "Complex Multivariate"
      ))
    )
  
  p5 <- ggplot(complexity_data, aes(x = complexity_category, y = stacking_weight,
                                    fill = complexity_category)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.6) +
    geom_jitter(width = 0.2, alpha = 0.4, size = 2) +
    scale_fill_viridis_d(name = "Model\nComplexity") +
    labs(
      title = "Stacking Weights by Model Complexity",
      subtitle = "Distribution of weights across complexity categories",
      x = "Model Complexity Category",
      y = "Stacking Weight"
    ) +
    theme_bw() +
    theme(
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(color = "black"),
      plot.subtitle = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90", size = 0.5),
      panel.grid.minor = element_line(color = "grey95", size = 0.25)
    ) +
    scale_y_continuous(labels = scales::percent_format())
  
  ggsave("results/figures/complexity_distribution.png", p5, 
         width = 10, height = 6, dpi = 300)
  
  cat("✓ Model complexity distribution plot saved\n")
  
} else {
  cat("⚠ No stacking weights available for visualization\n")
}

# SUMMARY OF GENERATED PLOTS =================================================

cat("\nVISUALIZATIONS GENERATED:\n")
cat("1. stacking_weights_reduced.png - Comparison of model weights by type\n")
cat("2. model_uncertainty_reduced.png - Model uncertainty with averaging threshold\n") 
cat("3. weighted_coefficients_reduced.png - Averaged coefficients across models\n")
cat("4. performance_comparison_heatmap.png - Heatmap of reduced vs univariate performance\n")
cat("5. complexity_distribution.png - Distribution of weights by model complexity\n")

cat("\nAll plots optimized for manuscript publication (300 DPI, clean themes)\n")

