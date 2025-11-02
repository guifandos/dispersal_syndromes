# COMPREHENSIVE MODEL COMPARISON ============================================
# Compare: Univariate vs Multivariate Complete vs Multivariate Selected vs Interactions
# ===========================================================================

library(brms)
library(loo)
library(dplyr)
library(tidyr)

# LOAD ALL MODEL TYPES =======================================================

load_all_models <- function() {
  
  cat("=== LOADING ALL MODEL TYPES ===\n")
  
  all_models <- list()
  
  # 1. UNIVARIATE MODELS -----------------------------------------------------
  cat("Loading univariate models...\n")
  
  dispersal_types <- c("average", "natal", "breeding")
  univariate_models <- list()
  
  for (type in dispersal_types) {
    file_path <- paste0("results/weibull/", type, "/univariate_models.RData")
    if (file.exists(file_path)) {
      load(file_path)  # Loads 'models' object
      univariate_models[[type]] <- models
      cat("  ✓", type, "univariate models loaded\n")
    } else {
      cat("  ✗", file_path, "not found\n")
    }
  }
  
  all_models$univariate <- univariate_models
  
  # 2. MULTIVARIATE COMPLETE MODELS -----------------------------------------
  cat("Loading multivariate complete models...\n")
  
  if (file.exists("results/combined/complete_models_combined.RData")) {
    load("results/combined/complete_models_combined.RData")
    all_models$multivariate_complete <- all_complete_models
    cat("  ✓ Complete multivariate models loaded\n")
  } else {
    cat("  ✗ results/combined/complete_models_combined.RData not found\n")
  }
  
  # 3. MULTIVARIATE SELECTED MODELS -----------------------------------------
  cat("Loading multivariate selected models...\n")
  
  if (file.exists("results/combined/short_models_combined.RData")) {
    load("results/combined/short_models_combined.RData")
    all_models$multivariate_selected <- all_short_models
    cat("  ✓ Selected multivariate models loaded\n")
  } else {
    cat("  ✗ results/combined/short_models_combined.RData not found\n")
  }
  
  # 4. INTERACTION MODELS ----------------------------------------------------
  cat("Loading interaction models...\n")
  
  if (file.exists("./results/models/best/best_model_weibull_partial_interactions.Rdata")) {
    load("./results/models/best/best_model_weibull_partial_interactions.Rdata")
    partial_models <- results_total
    all_models$interactions <- partial_models
    cat("  ✓ Interaction models loaded\n")
  } else {
    cat("  ✗ ./results/models/best/best_model_weibull_partial_interactions.Rdata not found\n")
  }
  
  cat("✓ All available models loaded\n\n")
  return(all_models)
}

# STANDARDIZE MODEL ACCESS ===================================================

get_model_standardized <- function(all_models, model_type, dispersal_type, distance_type, variable = NULL) {
  
  if (model_type == "univariate") {
    if (is.null(variable)) return(NULL)
    
    models <- all_models$univariate[[dispersal_type]]
    if (is.null(models)) return(NULL)
    
    model_name <- paste0(variable, "_", distance_type)
    return(models[[model_name]])
    
  } else if (model_type == "multivariate_complete") {
    models_df <- all_models$multivariate_complete
    if (is.null(models_df)) return(NULL)
    
    model_row <- models_df %>%
      filter(age == dispersal_type, dispersal_mode == distance_type)
    
    if (nrow(model_row) > 0) {
      return(model_row$model[[1]])
    }
    
  } else if (model_type == "multivariate_selected") {
    models_df <- all_models$multivariate_selected
    if (is.null(models_df)) return(NULL)
    
    model_row <- models_df %>%
      filter(age == dispersal_type, dispersal_mode == distance_type)
    
    if (nrow(model_row) > 0) {
      return(model_row$model[[1]])
    }
    
  } else if (model_type == "interactions") {
    # Handle both tibble format and list format
    if (is.data.frame(all_models$interactions) || "tbl" %in% class(all_models$interactions)) {
      # Tibble format (results_total)
      models_df <- all_models$interactions
      model_row <- models_df %>%
        filter(age == dispersal_type, type == distance_type)
      
      if (nrow(model_row) > 0) {
        return(model_row$model[[1]])
      }
    } else {
      # List format (from files)
      interaction_data <- all_models$interactions[[dispersal_type]]
      if (is.null(interaction_data)) return(NULL)
      
      if (distance_type == "median") {
        return(interaction_data$models$median_multi)
      } else {
        return(interaction_data$models$long_multi)
      }
    }
  }
  
  return(NULL)
}

# COMPREHENSIVE COMPARISON FUNCTION ==========================================

compare_all_model_types <- function(all_models, verbose = TRUE) {
  
  if (verbose) cat("=== COMPREHENSIVE MODEL COMPARISON ===\n")
  
  # Define comparisons to make
  dispersal_types <- c("average", "natal", "breeding")
  distance_types <- c("median", "long")
  univariate_variables <- c("body_mass", "PC1", "log_HWI", "diet", "habita_for", "distance_mig", "Latitude")
  
  # Results storage
  comparison_results <- data.frame(
    dispersal_type = character(),
    distance_type = character(),
    model_type = character(),
    variable = character(),
    delta_looic = numeric(),
    se_diff = numeric(),
    interpretation = character(),
    method = character(),
    converged = logical(),
    stringsAsFactors = FALSE
  )
  
  for (dispersal_type in dispersal_types) {
    for (distance_type in distance_types) {
      
      if (verbose) cat("\n--- Comparing", toupper(dispersal_type), toupper(distance_type), "dispersal ---\n")
      
      # Get reference model (try multivariate complete first)
      reference_model <- get_model_standardized(all_models, "multivariate_complete", dispersal_type, distance_type)
      reference_type <- "multivariate_complete"
      
      if (is.null(reference_model)) {
        reference_model <- get_model_standardized(all_models, "multivariate_selected", dispersal_type, distance_type)
        reference_type <- "multivariate_selected"
      }
      
      if (is.null(reference_model)) {
        if (verbose) cat("No reference model found for", dispersal_type, distance_type, "\n")
        next
      }
      
      if (verbose) cat("Using", reference_type, "as reference\n")
      
      # Add LOO to reference if needed
      if (!"loo" %in% names(reference_model$criteria)) {
        tryCatch({
          reference_model <- add_criterion(reference_model, "loo", moment_match = TRUE, reloo = TRUE)
        }, error = function(e) {
          reference_model <- add_criterion(reference_model, "waic")
        })
      }
      
      # Compare each model type
      model_types_to_compare <- c("multivariate_complete", "multivariate_selected", "interactions")
      
      # Add univariate models
      for (var in univariate_variables) {
        
        model <- get_model_standardized(all_models, "univariate", dispersal_type, distance_type, var)
        
        if (!is.null(model)) {
          
          if (verbose) cat("  Comparing univariate", var, "...\n")
          
          tryCatch({
            # Add criterion if needed
            if (!"loo" %in% names(model$criteria)) {
              model <- add_criterion(model, "loo", moment_match = TRUE, reloo = TRUE)
              method <- "loo"
            } else {
              method <- "loo"
            }
            
            # Compare
            comparison <- loo_compare(loo(reference_model), loo(model))
            delta <- comparison[2, "elpd_diff"]
            se_diff <- comparison[2, "se_diff"]
            
            interpretation <- case_when(
              delta > 2 * se_diff ~ "Strong improvement",
              delta > se_diff ~ "Moderate improvement",
              abs(delta) <= se_diff ~ "No difference",
              delta < -se_diff ~ "Worse performance",
              TRUE ~ "Inconclusive"
            )
            
            comparison_results <- rbind(comparison_results, data.frame(
              dispersal_type = dispersal_type,
              distance_type = distance_type,
              model_type = "univariate",
              variable = var,
              delta_looic = round(delta, 2),
              se_diff = round(se_diff, 2),
              interpretation = interpretation,
              method = method,
              converged = TRUE,
              stringsAsFactors = FALSE
            ))
            
            if (verbose) cat("    Δ LOOIC =", round(delta, 2), "- ", interpretation, "\n")
            
          }, error = function(e) {
            if (verbose) cat("    Error:", e$message, "\n")
            comparison_results <<- rbind(comparison_results, data.frame(
              dispersal_type = dispersal_type,
              distance_type = distance_type,
              model_type = "univariate",
              variable = var,
              delta_looic = NA,
              se_diff = NA,
              interpretation = "Failed",
              method = "failed",
              converged = FALSE,
              stringsAsFactors = FALSE
            ))
          })
        }
      }
      
      # Compare other multivariate models
      for (model_type in model_types_to_compare) {
        
        if (model_type == reference_type) next  # Skip reference
        
        model <- get_model_standardized(all_models, model_type, dispersal_type, distance_type)
        
        if (!is.null(model)) {
          
          if (verbose) cat("  Comparing", model_type, "...\n")
          
          tryCatch({
            # Add criterion if needed
            if (!"loo" %in% names(model$criteria)) {
              model <- add_criterion(model, "loo", moment_match = TRUE, reloo = TRUE)
              method <- "loo"
            } else {
              method <- "loo"
            }
            
            # Compare
            comparison <- loo_compare(loo(reference_model), loo(model))
            delta <- comparison[2, "elpd_diff"]
            se_diff <- comparison[2, "se_diff"]
            
            interpretation <- case_when(
              delta > 2 * se_diff ~ "Strong improvement",
              delta > se_diff ~ "Moderate improvement",
              abs(delta) <= se_diff ~ "No difference",
              delta < -se_diff ~ "Worse performance",
              TRUE ~ "Inconclusive"
            )
            
            comparison_results <- rbind(comparison_results, data.frame(
              dispersal_type = dispersal_type,
              distance_type = distance_type,
              model_type = model_type,
              variable = "full_model",
              delta_looic = round(delta, 2),
              se_diff = round(se_diff, 2),
              interpretation = interpretation,
              method = method,
              converged = TRUE,
              stringsAsFactors = FALSE
            ))
            
            if (verbose) cat("    Δ LOOIC =", round(delta, 2), "- ", interpretation, "\n")
            
          }, error = function(e) {
            if (verbose) cat("    Error:", e$message, "\n")
            comparison_results <<- rbind(comparison_results, data.frame(
              dispersal_type = dispersal_type,
              distance_type = distance_type,
              model_type = model_type,
              variable = "full_model",
              delta_looic = NA,
              se_diff = NA,
              interpretation = "Failed",
              method = "failed",
              converged = FALSE,
              stringsAsFactors = FALSE
            ))
          })
        }
      }
    }
  }
  
  return(comparison_results)
}

# CREATE COMPREHENSIVE SUMMARY TABLES =======================================

create_comprehensive_summaries <- function(comparison_results, all_models) {
  
  cat("\n=== CREATING COMPREHENSIVE SUMMARY TABLES ===\n")
  
  # 1. BEST MODELS BY AGE AND DISPERSAL TYPE --------------------------------
  cat("Creating best models by age and dispersal type table...\n")
  
  best_by_age_type <- comparison_results %>%
    filter(converged == TRUE) %>%
    group_by(dispersal_type, distance_type) %>%
    slice_max(delta_looic, n = 1, with_ties = FALSE) %>%
    arrange(dispersal_type, distance_type) %>%
    mutate(
      model_description = case_when(
        model_type == "univariate" ~ paste("Univariate:", variable),
        model_type == "multivariate_complete" ~ "Multivariate: All variables",
        model_type == "multivariate_selected" ~ "Multivariate: Selected variables", 
        model_type == "interactions" ~ "With interactions",
        TRUE ~ model_type
      ),
      performance_category = case_when(
        delta_looic > 4 ~ "Excellent (Δ > 4)",
        delta_looic > 2 ~ "Strong (Δ > 2)",
        delta_looic > 0 ~ "Moderate (Δ > 0)",
        delta_looic > -2 ~ "Similar (Δ ≈ 0)",
        TRUE ~ "Poor (Δ < -2)"
      )
    ) %>%
    select(dispersal_type, distance_type, model_type, variable, model_description, 
           delta_looic, performance_category, interpretation)
  
  # 2. OVERALL RANKING BY MODEL TYPE ----------------------------------------
  cat("Creating overall ranking table...\n")
  
  overall_ranking <- comparison_results %>%
    filter(converged == TRUE) %>%
    group_by(model_type) %>%
    summarise(
      n_contexts = n(),
      mean_delta = round(mean(delta_looic, na.rm = TRUE), 2),
      median_delta = round(median(delta_looic, na.rm = TRUE), 2),
      sd_delta = round(sd(delta_looic, na.rm = TRUE), 2),
      n_improvements = sum(grepl("improvement", interpretation)),
      n_worse = sum(interpretation == "Worse performance"),
      n_no_diff = sum(interpretation == "No difference"),
      best_performance = round(max(delta_looic, na.rm = TRUE), 2),
      worst_performance = round(min(delta_looic, na.rm = TRUE), 2),
      consistency = round(1 - (sd_delta / (abs(mean_delta) + 1)), 2),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_delta)) %>%
    mutate(
      overall_rank = row_number(),
      success_rate = round(n_improvements / n_contexts * 100, 1),
      recommendation = case_when(
        mean_delta > 2 & success_rate > 60 ~ "Highly Recommended",
        mean_delta > 0 & success_rate > 40 ~ "Recommended", 
        mean_delta > -1 & success_rate > 20 ~ "Consider with caution",
        TRUE ~ "Not recommended"
      )
    )
  
  # 3. DETAILED PERFORMANCE MATRIX ------------------------------------------
  cat("Creating detailed performance matrix...\n")
  
  performance_matrix <- comparison_results %>%
    filter(converged == TRUE) %>%
    select(dispersal_type, distance_type, model_type, variable, delta_looic) %>%
    mutate(
      model_id = case_when(
        model_type == "univariate" ~ paste(model_type, variable, sep = "_"),
        TRUE ~ model_type
      )
    ) %>%
    select(-model_type, -variable) %>%
    pivot_wider(
      names_from = c(dispersal_type, distance_type),
      values_from = delta_looic,
      names_sep = "_"
    ) %>%
    arrange(desc(rowMeans(select(., -model_id), na.rm = TRUE)))
  
  # 4. FILE LOCATIONS REFERENCE ---------------------------------------------
  cat("Creating file locations reference...\n")
  
  file_locations <- data.frame(
    model_type = c(
      "Univariate Models",
      "Univariate Models", 
      "Univariate Models",
      "Multivariate Complete",
      "Multivariate Selected",
      "Interaction Models"
    ),
    dispersal_age = c(
      "average", "natal", "breeding",
      "all", "all", "all"
    ),
    file_path = c(
      "results/weibull/average/univariate_models.RData",
      "results/weibull/natal/univariate_models.RData",
      "results/weibull/breeding/univariate_models.RData", 
      "results/combined/complete_models_combined.RData",
      "results/combined/short_models_combined.RData",
      "./results/models/best/best_model_weibull_partial_interactions.Rdata"
    ),
    object_name = c(
      "models", "models", "models",
      "all_complete_models", "all_short_models", "results_total"
    ),
    description = c(
      "Individual variable models for average dispersal",
      "Individual variable models for natal dispersal", 
      "Individual variable models for breeding dispersal",
      "Full multivariate models with all variables",
      "Variable-selected multivariate models",
      "Models with ecological interactions"
    ),
    stringsAsFactors = FALSE
  )
  
  # 5. MODEL CHARACTERISTICS SUMMARY ----------------------------------------
  cat("Creating model characteristics summary...\n")
  
  model_characteristics <- data.frame(
    model_type = c("univariate", "multivariate_complete", "multivariate_selected", "interactions"),
    n_variables = c("1", "7", "3-5", "7+"),
    complexity = c("Low", "High", "Medium", "Very High"),
    interpretability = c("High", "Medium", "High", "Low"),
    overfitting_risk = c("Low", "High", "Medium", "Very High"),
    computational_cost = c("Low", "Medium", "Medium", "High"),
    cross_validation_stability = c("High", "Medium", "Medium", "Low"),
    recommended_use = c(
      "Prediction, simple interpretation",
      "Understanding full syndromes", 
      "Balanced prediction and interpretation",
      "Testing specific ecological hypotheses"
    ),
    stringsAsFactors = FALSE
  )
  
  return(list(
    best_by_context = best_by_age_type,
    overall_ranking = overall_ranking,
    performance_matrix = performance_matrix,
    file_locations = file_locations,
    model_characteristics = model_characteristics,
    detailed_results = comparison_results
  ))
}

# MAIN EXECUTION FUNCTION ====================================================

run_comprehensive_comparison <- function(save_results = TRUE) {
  
  cat("=== COMPREHENSIVE MODEL COMPARISON ANALYSIS ===\n")
  cat("Comparing: Univariate vs Multivariate Complete vs Multivariate Selected vs Interactions\n\n")
  
  # Load all models
  all_models <- load_all_models()
  
  # Run comparisons
  comparison_results <- compare_all_model_types(all_models)
  
  # Create comprehensive summaries
  summaries <- create_comprehensive_summaries(comparison_results, all_models)
  
  # Print key results
  cat("\n=== BEST MODELS BY CONTEXT ===\n")
  print(summaries$best_by_context)
  
  cat("\n=== OVERALL MODEL RANKING ===\n")
  print(summaries$overall_ranking)
  
  cat("\n=== PERFORMANCE MATRIX ===\n")
  print(summaries$performance_matrix)
  
  # Save all results
  if (save_results) {
    if (!dir.exists("results/comprehensive_comparison")) {
      dir.create("results/comprehensive_comparison", recursive = TRUE)
    }
    
    # Save individual tables
    write.csv(summaries$best_by_context, "results/comprehensive_comparison/best_models_by_context.csv", row.names = FALSE)
    write.csv(summaries$overall_ranking, "results/comprehensive_comparison/overall_model_ranking.csv", row.names = FALSE)
    write.csv(summaries$performance_matrix, "results/comprehensive_comparison/performance_matrix.csv", row.names = FALSE)
    write.csv(summaries$file_locations, "results/comprehensive_comparison/file_locations_reference.csv", row.names = FALSE)
    write.csv(summaries$model_characteristics, "results/comprehensive_comparison/model_characteristics.csv", row.names = FALSE)
    write.csv(summaries$detailed_results, "results/comprehensive_comparison/detailed_comparisons.csv", row.names = FALSE)
    
    # Save combined RData
    save(all_models, comparison_results, summaries, 
         file = "results/comprehensive_comparison/comprehensive_comparison_complete.RData")
    
    cat("\n=== ALL TABLES SAVED ===\n")
    cat("- best_models_by_context.csv (Best model for each age × dispersal type)\n")
    cat("- overall_model_ranking.csv (Overall performance ranking)\n") 
    cat("- performance_matrix.csv (Δ LOOIC matrix across contexts)\n")
    cat("- file_locations_reference.csv (Where to find each model type)\n")
    cat("- model_characteristics.csv (Model properties comparison)\n")
    cat("- detailed_comparisons.csv (All individual comparisons)\n")
    cat("- comprehensive_comparison_complete.RData (Complete analysis)\n")
    cat("\nAll files in: results/comprehensive_comparison/\n")
  }
  
  return(list(
    models = all_models,
    comparisons = comparison_results,
    summaries = summaries
  ))
}

# EXECUTE ANALYSIS ============================================================

cat("Do you want to run the comprehensive comparison? (y/n): ")
response <- readline()

if (tolower(response) %in% c("y", "yes", "1")) {
  results <- run_comprehensive_comparison(save_results = TRUE)
  cat("\n✅ Comprehensive comparison completed!\n")
  
  # Show key findings
  cat("\n=== KEY FINDINGS ===\n")
  best_overall <- results$summaries$overall[1, ]
  cat("Best overall model type:", best_overall$model_type, "\n")
  cat("Mean Δ LOOIC:", best_overall$mean_delta, "\n")
  cat("Improvements vs reference:", best_overall$n_improvements, "/", best_overall$n_comparisons, "\n")
  
} else {
  cat("\nComparison not run. Use run_comprehensive_comparison() when ready.\n")
}