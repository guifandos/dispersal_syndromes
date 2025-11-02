# CLEAN MODEL COMPARISON TABLE ===============================================
# Compare: Univariate vs Full vs Reduced vs Interactions in single table
# =========================================================================

library(brms)
library(loo)
library(dplyr)
library(tidyr)

# HELPER FUNCTION FOR ROBUST LOO =============================================

ensure_robust_loo <- function(model, verbose = FALSE) {
  if (!"loo" %in% names(model$criteria)) {
    if (verbose) cat("Adding robust LOO criterion...\n")
    tryCatch({
      model <- add_criterion(model, "loo", moment_match = TRUE, reloo = TRUE)
    }, error = function(e) {
      if (verbose) cat("LOO failed, trying WAIC...\n")
      model <- add_criterion(model, "waic")
    })
  }
  return(model)
}

# LOAD ALL MODELS =============================================================

load_models_for_comparison <- function() {
  
  cat("=== LOADING MODELS FOR CLEAN COMPARISON ===\n")
  
  models <- list()
  
  # 1. Univariate models
  dispersal_types <- c("average", "natal", "breeding")
  models$univariate <- list()
  
  for (type in dispersal_types) {
    file_path <- paste0("results/weibull/", type, "/univariate_models.RData")
    if (file.exists(file_path)) {
      load(file_path)
      models$univariate[[type]] <- models
      cat("✓ Loaded", type, "univariate models\n")
    }
  }
  
  # 2. Full models (complete)
  if (file.exists("results/combined/complete_models_combined.RData")) {
    load("results/combined/complete_models_combined.RData")
    if (exists("all_complete_models") && !is.null(all_complete_models)) {
      models$full <- all_complete_models
      cat("✓ Loaded full models:", nrow(all_complete_models), "models\n")
    } else {
      cat("✗ all_complete_models is NULL or doesn't exist\n")
      models$full <- NULL
    }
  } else {
    cat("✗ complete_models_combined.RData not found\n")
    models$full <- NULL
  }
  
  # 3. Reduced models (selected variables)
  if (file.exists("results/combined/short_models_combined.RData")) {
    load("results/combined/short_models_combined.RData")
    if (exists("all_short_models") && !is.null(all_short_models)) {
      models$reduced <- all_short_models
      cat("✓ Loaded reduced models:", nrow(all_short_models), "models\n")
    } else {
      cat("✗ all_short_models is NULL or doesn't exist\n")
      models$reduced <- NULL
    }
  } else {
    cat("✗ short_models_combined.RData not found\n")
    models$reduced <- NULL
  }
  
  # 4. Interaction models
  if (file.exists("./results/models/best/best_model_weibull_partial_interactions.Rdata")) {
    load("./results/models/best/best_model_weibull_partial_interactions.Rdata")
    if (exists("results_total") && !is.null(results_total)) {
      models$interactions <- results_total
      cat("✓ Loaded interaction models:", nrow(results_total), "models\n")
    } else {
      cat("✗ results_total is NULL or doesn't exist\n")
      models$interactions <- NULL
    }
  } else {
    cat("✗ best_model_weibull_partial_interactions.Rdata not found\n")
    models$interactions <- NULL
  }
  
  cat("Model loading complete\n\n")
  return(models)
}

# GET MODEL FUNCTION ==========================================================

get_model <- function(models, model_type, age, distance_type, variable = NULL) {
  
  if (model_type == "univariate") {
    if (is.null(variable) || is.null(models$univariate) || is.null(models$univariate[[age]])) {
      return(NULL)
    }
    model_name <- paste0(variable, "_", distance_type)
    return(models$univariate[[age]][[model_name]])
    
  } else if (model_type == "full") {
    if (is.null(models$full)) return(NULL)
    
    tryCatch({
      model_row <- models$full %>%
        filter(age == !!age, dispersal_mode == distance_type)
      if (nrow(model_row) > 0) return(model_row$model[[1]])
    }, error = function(e) {
      cat("Error filtering full models:", e$message, "\n")
      return(NULL)
    })
    
  } else if (model_type == "reduced") {
    if (is.null(models$reduced)) return(NULL)
    
    tryCatch({
      model_row <- models$reduced %>%
        filter(age == !!age, dispersal_mode == distance_type)
      if (nrow(model_row) > 0) return(model_row$model[[1]])
    }, error = function(e) {
      cat("Error filtering reduced models:", e$message, "\n")
      return(NULL)
    })
    
  } else if (model_type == "interactions") {
    if (is.null(models$interactions)) return(NULL)
    
    tryCatch({
      model_row <- models$interactions %>%
        filter(age == !!age, type == distance_type)
      if (nrow(model_row) > 0) return(model_row$model[[1]])
    }, error = function(e) {
      cat("Error filtering interaction models:", e$message, "\n")
      return(NULL)
    })
  }
  
  return(NULL)
}

# FIND BEST UNIVARIATE MODEL =================================================

find_best_univariate <- function(models, age, distance_type, reference_model) {
  
  variables <- c("body_mass", "PC1", "log_HWI", "diet", "habita_for", "distance_mig", "Latitude")
  
  best_variable <- NA
  best_delta <- -Inf
  best_interpretation <- "Failed"
  
  for (var in variables) {
    
    model <- get_model(models, "univariate", age, distance_type, var)
    
    if (!is.null(model)) {
      tryCatch({
        model <- ensure_robust_loo(model)
        
        comparison <- loo_compare(loo(reference_model), loo(model))
        delta <- comparison[2, "elpd_diff"]
        
        if (delta > best_delta) {
          best_delta <- delta
          best_variable <- var
          
          se_diff <- comparison[2, "se_diff"]
          best_interpretation <- case_when(
            delta > 2 * se_diff ~ "Strong improvement",
            delta > se_diff ~ "Moderate improvement",
            abs(delta) <= se_diff ~ "No difference",
            delta < -se_diff ~ "Worse performance",
            TRUE ~ "Inconclusive"
          )
        }
        
      }, error = function(e) {
        # Skip this variable if comparison fails
      })
    }
  }
  
  return(list(
    variable = best_variable,
    delta = ifelse(is.finite(best_delta), round(best_delta, 2), NA),
    interpretation = best_interpretation
  ))
}

# COMPARE SINGLE CONTEXT =====================================================

compare_single_context <- function(models, age, distance_type, verbose = FALSE) {
  
  if (verbose) cat("Comparing", age, distance_type, "dispersal...\n")
  
  # Get reference model (full model)
  reference_model <- get_model(models, "full", age, distance_type)
  
  if (is.null(reference_model)) {
    if (verbose) cat("  No full model found, trying reduced model as reference...\n")
    reference_model <- get_model(models, "reduced", age, distance_type)
    
    if (is.null(reference_model)) {
      if (verbose) cat("  No reference model available\n")
      return(data.frame(
        age = age,
        distance_type = distance_type,
        best_univariate_var = "No reference",
        univariate_delta = NA,
        full_delta = NA,
        reduced_delta = NA,
        interactions_delta = NA,
        winner = "No comparison possible",
        univariate_interp = "No reference",
        reduced_interp = "No reference",
        interactions_interp = "No reference",
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Ensure reference has robust LOO
  tryCatch({
    reference_model <- ensure_robust_loo(reference_model, verbose)
  }, error = function(e) {
    if (verbose) cat("  Reference model LOO failed:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(reference_model)) return(NULL)
  
  # Find best univariate
  best_univariate <- find_best_univariate(models, age, distance_type, reference_model)
  
  # Compare reduced model
  reduced_delta <- NA
  reduced_interpretation <- "Not available"
  reduced_model <- get_model(models, "reduced", age, distance_type)
  
  if (!is.null(reduced_model)) {
    tryCatch({
      reduced_model <- ensure_robust_loo(reduced_model)
      comparison <- loo_compare(loo(reference_model), loo(reduced_model))
      reduced_delta <- round(comparison[2, "elpd_diff"], 2)
      
      se_diff <- comparison[2, "se_diff"]
      reduced_interpretation <- case_when(
        reduced_delta > 2 * se_diff ~ "Strong improvement",
        reduced_delta > se_diff ~ "Moderate improvement",
        abs(reduced_delta) <= se_diff ~ "No difference",
        reduced_delta < -se_diff ~ "Worse performance",
        TRUE ~ "Inconclusive"
      )
    }, error = function(e) {
      reduced_interpretation <- "Failed"
    })
  }
  
  # Compare interactions model
  interactions_delta <- NA
  interactions_interpretation <- "Not available"
  interactions_model <- get_model(models, "interactions", age, distance_type)
  
  if (!is.null(interactions_model)) {
    tryCatch({
      interactions_model <- ensure_robust_loo(interactions_model)
      comparison <- loo_compare(loo(reference_model), loo(interactions_model))
      interactions_delta <- round(comparison[2, "elpd_diff"], 2)
      
      se_diff <- comparison[2, "se_diff"]
      interactions_interpretation <- case_when(
        interactions_delta > 2 * se_diff ~ "Strong improvement",
        interactions_delta > se_diff ~ "Moderate improvement",
        abs(interactions_delta) <= se_diff ~ "No difference",
        interactions_delta < -se_diff ~ "Worse performance",
        TRUE ~ "Inconclusive"
      )
    }, error = function(e) {
      interactions_interpretation <- "Failed"
    })
  }
  
  # Determine winner
  deltas <- c(
    univariate = ifelse(is.na(best_univariate$delta), -Inf, best_univariate$delta),
    full = 0,  # Reference
    reduced = ifelse(is.na(reduced_delta), -Inf, reduced_delta),
    interactions = ifelse(is.na(interactions_delta), -Inf, interactions_delta)
  )
  
  winner <- names(deltas)[which.max(deltas)]
  
  # Format winner name
  winner_formatted <- case_when(
    winner == "univariate" ~ paste0("Univariate (", best_univariate$variable, ")"),
    winner == "full" ~ "Full model",
    winner == "reduced" ~ "Reduced model", 
    winner == "interactions" ~ "With interactions",
    TRUE ~ winner
  )
  
  return(data.frame(
    age = age,
    distance_type = distance_type,
    best_univariate_var = ifelse(is.na(best_univariate$variable), "None", best_univariate$variable),
    univariate_delta = best_univariate$delta,
    full_delta = 0,  # Reference
    reduced_delta = reduced_delta,
    interactions_delta = interactions_delta,
    winner = winner_formatted,
    univariate_interp = best_univariate$interpretation,
    reduced_interp = reduced_interpretation,
    interactions_interp = interactions_interpretation,
    stringsAsFactors = FALSE
  ))
}

# MAIN COMPARISON FUNCTION ===================================================

create_clean_comparison_table <- function(save_results = TRUE) {
  
  cat("=== CREATING CLEAN MODEL COMPARISON TABLE ===\n\n")
  
  # Load models
  models <- load_models_for_comparison()
  
  # Define contexts
  ages <- c("average", "natal", "breeding")
  distance_types <- c("median", "long")
  
  # Create comparison table
  comparison_table <- data.frame()
  
  for (age in ages) {
    for (distance_type in distance_types) {
      
      result <- compare_single_context(models, age, distance_type, verbose = TRUE)
      
      if (!is.null(result)) {
        comparison_table <- rbind(comparison_table, result)
      }
    }
  }
  
  # Create display table (clean version)
  display_table <- comparison_table %>%
    select(age, distance_type, best_univariate_var, univariate_delta, 
           full_delta, reduced_delta, interactions_delta, winner) %>%
    rename(
      "Age" = age,
      "Distance" = distance_type,
      "Best Univariate" = best_univariate_var,
      "Univariate Δ" = univariate_delta,
      "Full Δ" = full_delta,
      "Reduced Δ" = reduced_delta,
      "Interactions Δ" = interactions_delta,
      "Winner" = winner
    )
  
  # Create summary statistics
  summary_stats <- data.frame(
    model_type = c("Univariate", "Full (reference)", "Reduced", "Interactions"),
    n_wins = c(
      sum(grepl("Univariate", comparison_table$winner)),
      sum(grepl("Full", comparison_table$winner)),
      sum(grepl("Reduced", comparison_table$winner)),
      sum(grepl("interactions", comparison_table$winner))
    ),
    mean_delta = c(
      round(mean(comparison_table$univariate_delta, na.rm = TRUE), 2),
      0,
      round(mean(comparison_table$reduced_delta, na.rm = TRUE), 2),
      round(mean(comparison_table$interactions_delta, na.rm = TRUE), 2)
    ),
    best_performance = c(
      round(max(comparison_table$univariate_delta, na.rm = TRUE), 2),
      0,
      round(max(comparison_table$reduced_delta, na.rm = TRUE), 2),
      round(max(comparison_table$interactions_delta, na.rm = TRUE), 2)
    ),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(n_wins), desc(mean_delta))
  
  # Print results
  cat("\n=== CLEAN MODEL COMPARISON TABLE ===\n")
  print(display_table, row.names = FALSE)
  
  cat("\n=== SUMMARY STATISTICS ===\n")
  print(summary_stats, row.names = FALSE)
  
  # Interpretation
  cat("\n=== INTERPRETATION ===\n")
  cat("• Δ > 0: Better than full model (reference)\n")
  cat("• Δ < 0: Worse than full model\n")
  cat("• Winner: Best performing model for each context\n\n")
  
  # Key findings
  best_overall <- summary_stats$model_type[1]
  cat("KEY FINDINGS:\n")
  cat("• Overall best model type:", best_overall, "\n")
  cat("• Contexts where univariate wins:", sum(grepl("Univariate", comparison_table$winner)), "/", nrow(comparison_table), "\n")
  
  # Most successful univariate variable
  if (sum(!is.na(comparison_table$best_univariate_var)) > 0) {
    best_vars <- table(comparison_table$best_univariate_var[!is.na(comparison_table$best_univariate_var)])
    most_successful_var <- names(best_vars)[which.max(best_vars)]
    cat("• Most successful univariate variable:", most_successful_var, "\n")
  }
  
  # Save results
  if (save_results) {
    if (!dir.exists("results/clean_comparison")) {
      dir.create("results/clean_comparison", recursive = TRUE)
    }
    
    write.csv(display_table, "results/clean_comparison/model_comparison_table.csv", row.names = FALSE)
    write.csv(summary_stats, "results/clean_comparison/summary_statistics.csv", row.names = FALSE)
    write.csv(comparison_table, "results/clean_comparison/detailed_comparison.csv", row.names = FALSE)
    
    save(models, comparison_table, display_table, summary_stats,
         file = "results/clean_comparison/clean_comparison_results.RData")
    
    cat("\nResults saved in: results/clean_comparison/\n")
    cat("• model_comparison_table.csv (main table)\n")
    cat("• summary_statistics.csv (overall stats)\n") 
    cat("• detailed_comparison.csv (with interpretations)\n")
  }
  
  return(list(
    display_table = display_table,
    summary_stats = summary_stats,
    detailed = comparison_table
  ))
}

# EXECUTE ANALYSIS ============================================================

cat("Running clean model comparison...\n\n")
results <- create_clean_comparison_table(save_results = TRUE)

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Check results/clean_comparison/ for output files\n")