# DIAGNOSTIC LOO TESTING ====================================================

library(brms)
library(loo)
library(dplyr)

# FUNCTION TO TEST INDIVIDUAL MODEL LOO =====================================

test_individual_model_loo <- function(model, model_name) {
  
  cat("=== TESTING LOO FOR", toupper(model_name), "===\n")
  
  # Check model object
  cat("1. Model object check:\n")
  cat("   Class:", class(model), "\n")
  cat("   Family:", family(model)$family, "\n")
  cat("   Formula:", deparse1(formula(model)), "\n")
  cat("   N obs:", nobs(model), "\n")
  cat("   Chains:", model$fit@sim$chains, "\n")
  cat("   Iterations:", model$fit@sim$iter, "\n")
  
  # Check existing criteria
  cat("\n2. Existing criteria:\n")
  if (length(model$criteria) > 0) {
    cat("   Found:", names(model$criteria), "\n")
  } else {
    cat("   No criteria found\n")
  }
  
  # Test LOO calculation
  cat("\n3. Testing LOO calculation:\n")
  
  tryCatch({
    cat("   Attempting basic LOO...\n")
    loo_result <- loo(model, refresh = 0)
    cat("   ✓ Basic LOO successful\n")
    cat("   LOOIC:", round(loo_result$estimates["looic", "Estimate"], 2), "\n")
    cat("   SE:", round(loo_result$estimates["looic", "SE"], 2), "\n")
    
    # Check for problematic observations
    if (!is.null(loo_result$diagnostics)) {
      pareto_k <- loo_result$diagnostics$pareto_k
      high_k <- sum(pareto_k > 0.7, na.rm = TRUE)
      cat("   Problematic obs (k > 0.7):", high_k, "/", length(pareto_k), "\n")
      
      if (high_k > 0) {
        cat("   WARNING: High Pareto-k values detected\n")
        
        # Try moment matching
        cat("   Attempting LOO with moment matching...\n")
        tryCatch({
          loo_mm <- loo(model, moment_match = TRUE, refresh = 0)
          cat("   ✓ Moment matching successful\n")
          return(list(success = TRUE, method = "loo_mm", loo = loo_mm))
        }, error = function(e2) {
          cat("   ✗ Moment matching failed:", e2$message, "\n")
        })
        
        # Try reloo
        cat("   Attempting reloo...\n")
        tryCatch({
          loo_reloo <- loo(model, reloo = TRUE, refresh = 0)
          cat("   ✓ Reloo successful\n")
          return(list(success = TRUE, method = "reloo", loo = loo_reloo))
        }, error = function(e3) {
          cat("   ✗ Reloo failed:", e3$message, "\n")
        })
      }
    }
    
    return(list(success = TRUE, method = "basic", loo = loo_result))
    
  }, error = function(e1) {
    cat("   ✗ Basic LOO failed:", e1$message, "\n")
    
    # Try WAIC as fallback
    cat("   Attempting WAIC fallback...\n")
    tryCatch({
      waic_result <- waic(model)
      cat("   ✓ WAIC successful\n")
      cat("   WAIC:", round(waic_result$estimates["waic", "Estimate"], 2), "\n")
      return(list(success = TRUE, method = "waic", loo = waic_result))
    }, error = function(e2) {
      cat("   ✗ WAIC also failed:", e2$message, "\n")
      return(list(success = FALSE, method = "failed", error = e1$message))
    })
  })
}

# FUNCTION TO TEST MODEL COMPARISON ==========================================

test_model_comparison <- function(model1, model2, name1, name2) {
  
  cat("\n=== TESTING MODEL COMPARISON ===\n")
  cat("Comparing:", name1, "vs", name2, "\n")
  
  # Test both models first
  result1 <- test_individual_model_loo(model1, name1)
  result2 <- test_individual_model_loo(model2, name2)
  
  if (!result1$success || !result2$success) {
    cat("Cannot compare - one or both models failed LOO/WAIC\n")
    return(NULL)
  }
  
  # Attempt comparison
  cat("\nAttempting comparison...\n")
  tryCatch({
    
    if (result1$method == result2$method) {
      comparison <- loo_compare(result1$loo, result2$loo)
      cat("✓ Comparison successful using", result1$method, "\n")
      print(comparison)
      
      delta <- comparison[2, "elpd_diff"]
      se_diff <- comparison[2, "se_diff"]
      cat("Δ:", round(delta, 2), "± SE:", round(se_diff, 2), "\n")
      
      return(list(
        success = TRUE,
        delta = delta,
        se_diff = se_diff,
        method = result1$method,
        comparison = comparison
      ))
    } else {
      cat("✗ Cannot compare - different methods used\n")
      cat("Model 1 method:", result1$method, "\n")
      cat("Model 2 method:", result2$method, "\n")
      return(NULL)
    }
    
  }, error = function(e) {
    cat("✗ Comparison failed:", e$message, "\n")
    return(NULL)
  })
}

# MAIN DIAGNOSTIC FUNCTION ==================================================

diagnose_model_comparison_issues <- function() {
  
  cat("=== COMPREHENSIVE MODEL COMPARISON DIAGNOSIS ===\n\n")
  
  # Load models
  cat("Loading models for diagnosis...\n")
  
  # Univariate models
  load("results/weibull/average/univariate_models.RData")
  univariate_models <- models
  
  # Reference models
  load("results/combined/complete_models_combined.RData")
  complete_models <- all_complete_models
  
  cat("✓ Models loaded\n\n")
  
  # Find reference model
  ref_model_row <- complete_models %>%
    filter(age == "average", dispersal_mode == "median")
  
  if (nrow(ref_model_row) == 0) {
    cat("✗ No reference model found\n")
    return()
  }
  
  ref_model <- ref_model_row$model[[1]]
  cat("Reference model: Complete multivariate (average median)\n")
  
  # Test reference model
  ref_result <- test_individual_model_loo(ref_model, "reference_complete")
  
  if (!ref_result$success) {
    cat("\n❌ REFERENCE MODEL FAILS LOO - This is the root problem!\n")
    return()
  }
  
  cat("\n✓ Reference model LOO works\n")
  
  # Test univariate models
  cat("\n" + paste(rep("=", 60), collapse = "") + "\n")
  cat("TESTING UNIVARIATE MODELS\n")
  
  univariate_names <- names(univariate_models)
  median_models <- univariate_names[grepl("_median$", univariate_names)]
  
  for (model_name in median_models) {
    
    cat("\n" + paste(rep("-", 40), collapse = "") + "\n")
    
    model <- univariate_models[[model_name]]
    
    if (is.null(model)) {
      cat("Model", model_name, "is NULL\n")
      next
    }
    
    # Test individual model
    result <- test_individual_model_loo(model, model_name)
    
    if (result$success) {
      cat("\n✓ Model", model_name, "LOO works\n")
      
      # Test comparison with reference
      cat("Testing comparison with reference...\n")
      comp_result <- test_model_comparison(ref_model, model, "reference", model_name)
      
      if (!is.null(comp_result)) {
        cat("✓ Comparison successful!\n")
      }
    } else {
      cat("\n✗ Model", model_name, "LOO failed\n")
    }
  }
  
  cat("\n=== DIAGNOSIS COMPLETE ===\n")
}

# QUICK TEST FUNCTION ========================================================

quick_test_specific_model <- function(model_name = "body_mass_median") {
  
  cat("=== QUICK TEST:", toupper(model_name), "===\n")
  
  # Load univariate models
  load("results/weibull/average/univariate_models.RData")
  
  if (!model_name %in% names(models)) {
    cat("Available models:", paste(names(models), collapse = ", "), "\n")
    return()
  }
  
  model <- models[[model_name]]
  result <- test_individual_model_loo(model, model_name)
  
  if (result$success) {
    cat("\n✅ SUCCESS: Model", model_name, "LOO works fine\n")
    cat("Method used:", result$method, "\n")
    cat("This suggests the problem is in the comparison function, not the individual models\n")
  } else {
    cat("\n❌ FAILED: Model", model_name, "has LOO issues\n")
    cat("Error:", result$error, "\n")
  }
}

# EXECUTE TESTS ==============================================================

cat("Choose diagnostic test:\n")
cat("1. quick_test_specific_model('body_mass_median')\n")
cat("2. diagnose_model_comparison_issues()  # Full diagnosis\n\n")

# Run quick test by default
cat("Running quick test for body_mass_median...\n\n")
quick_test_specific_model("body_mass_median")