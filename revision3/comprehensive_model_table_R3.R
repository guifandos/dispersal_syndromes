# COMPREHENSIVE MODEL COMPARISON TABLE ======================================
# Create table like the provided one but including reduced models

library(brms)
library(loo)
library(dplyr)
library(tidyr)

# Increase memory limit for LOO calculations
options(future.globals.maxSize = 1000 * 1024^2)  # 1GB limit

cat("Memory limit set to 1GB for LOO calculations\n")

# CREATE RESULTS FOLDER ======================================================
if (!dir.exists("results/comparison_R3")) {
  dir.create("results/comparison_R3", recursive = TRUE)
}

# LOAD ALL MODELS =============================================================
cat("Loading all models for comprehensive comparison...\n")

# Univariate models
dispersal_types <- c("average", "natal", "breeding")
distance_types <- c("median", "long")

all_results <- data.frame()

for (disp_type in dispersal_types) {
  
  # Load univariate models (only if not already loaded)
  if (!exists("univariate_models") || !all(names(univariate_models) %in% paste0(
    c("body_mass", "PC1", "log_HWI", "diet", "habita_for", "distance_mig", "Latitude"),
    "_", rep(distance_types, each = 7)
  ))) {
    load(paste0("results/weibull/", disp_type, "/univariate_models.RData"))
    univariate_models <- models
  }
  
  # Load complete models only once
  if (!exists("complete_models")) {
    load("results/combined/complete_models_combined.RData")
    complete_models <- all_complete_models
  }
  
  # Load reduced models only once
  if (!exists("reduced_models")) {
    load("results/combined/short_models_combined.RData")
    reduced_models <- all_short_models
  }
  
  for (dist_type in distance_types) {
    
    cat("Processing", disp_type, dist_type, "dispersal...\n")
    
    # Get multivariate models
    complete_model <- complete_models %>%
      filter(age == disp_type, dispersal_mode == dist_type) %>%
      pull(model) %>% .[[1]]
    
    reduced_model <- reduced_models %>%
      filter(age == disp_type, dispersal_mode == dist_type) %>%
      pull(model) %>% .[[1]]
    
    # Process univariate models
    variables <- c("body_mass", "PC1", "log_HWI", "diet", "habita_for", "distance_mig", "Latitude")
    
    for (var in variables) {
      
      model_name <- paste0(var, "_", dist_type)
      model <- univariate_models[[model_name]]
      
      if (!is.null(model)) {
        
        # Add LOO if needed - with memory-efficient options
        if (!"loo" %in% names(model$criteria)) {
          tryCatch({
            model <- add_criterion(model, "loo")
          }, error = function(e1) {
            tryCatch({
              model <- add_criterion(model, "loo", moment_match = TRUE)
            }, error = function(e2) {
              cat("  LOO failed for", var, ", using WAIC\n")
              model <- add_criterion(model, "waic")
            })
          })
        }
        
        # Get LOO stats
        loo_result <- loo(model)
        elpd_loo <- loo_result$estimates["elpd_loo", "Estimate"]
        looic <- loo_result$estimates["looic", "Estimate"]
        se_looic <- loo_result$estimates["looic", "SE"]
        p_loo <- loo_result$estimates["p_loo", "Estimate"]
        
        all_results <- rbind(all_results, data.frame(
          model_type = "univariate",
          dispersal_type = disp_type,
          distance_type = dist_type,
          variable = var,
          model_id = paste0(disp_type, "_", dist_type, "_", var, "_uni"),
          elpd_loo = round(elpd_loo, 2),
          looic = round(looic, 2),
          se_looic = round(se_looic, 2),
          p_loo = round(p_loo, 2),
          n_params = 3,
          stringsAsFactors = FALSE
        ))
      }
    }
    
    # COMPLETE MULTIVARIATE ==================================================
    if (!is.null(complete_model)) {
      if (!"loo" %in% names(complete_model$criteria)) {
        tryCatch({
          complete_model <- add_criterion(complete_model, "loo")
        }, error = function(e1) {
          tryCatch({
            complete_model <- add_criterion(complete_model, "loo", moment_match = TRUE)
          }, error = function(e2) {
            cat("  LOO failed for complete model, using WAIC\n")
            complete_model <- add_criterion(complete_model, "waic")
          })
        })
      }
      
      loo_result <- loo(complete_model)
      elpd_loo <- loo_result$estimates["elpd_loo", "Estimate"]
      looic <- loo_result$estimates["looic", "Estimate"]
      se_looic <- loo_result$estimates["looic", "SE"]
      p_loo <- loo_result$estimates["p_loo", "Estimate"]
      
      all_results <- rbind(all_results, data.frame(
        model_type = "multivariate_complete",
        dispersal_type = disp_type,
        distance_type = dist_type,
        variable = "full_model",
        model_id = paste0(disp_type, "_", dist_type, "_complete_multi"),
        elpd_loo = round(elpd_loo, 2),
        looic = round(looic, 2),
        se_looic = round(se_looic, 2),
        p_loo = round(p_loo, 2),
        n_params = length(fixef(complete_model)[,1]),
        stringsAsFactors = FALSE
      ))
    }
    
    # REDUCED MULTIVARIATE ===================================================
    if (!is.null(reduced_model)) {
      if (!"loo" %in% names(reduced_model$criteria)) {
        tryCatch({
          reduced_model <- add_criterion(reduced_model, "loo")
        }, error = function(e1) {
          tryCatch({
            reduced_model <- add_criterion(reduced_model, "loo", moment_match = TRUE)
          }, error = function(e2) {
            cat("  LOO failed for reduced model, using WAIC\n")
            reduced_model <- add_criterion(reduced_model, "waic")
          })
        })
      }
      
      loo_result <- loo(reduced_model)
      elpd_loo <- loo_result$estimates["elpd_loo", "Estimate"]
      looic <- loo_result$estimates["looic", "Estimate"]
      se_looic <- loo_result$estimates["looic", "SE"]
      p_loo <- loo_result$estimates["p_loo", "Estimate"]
      
      all_results <- rbind(all_results, data.frame(
        model_type = "multivariate_reduced",
        dispersal_type = disp_type,
        distance_type = dist_type,
        variable = "reduced_model",
        model_id = paste0(disp_type, "_", dist_type, "_reduced_multi"),
        elpd_loo = round(elpd_loo, 2),
        looic = round(looic, 2),
        se_looic = round(se_looic, 2),
        p_loo = round(p_loo, 2),
        n_params = length(fixef(reduced_model)[,1]),
        stringsAsFactors = FALSE
      ))
    }
  }
}

# CALCULATE RANKINGS AND PERFORMANCE =========================================
cat("Calculating rankings and performance metrics...\n")

final_results <- all_results %>%
  group_by(dispersal_type, distance_type) %>%
  mutate(
    best_looic = min(looic, na.rm = TRUE),
    best_elpd = max(elpd_loo, na.rm = TRUE),
    delta_looic = round(looic - best_looic, 2),
    delta_elpd = round(elpd_loo - best_elpd, 2),
    looic_rank = rank(looic),
    elpd_rank = rank(-elpd_loo),
    se_units = round(abs(delta_elpd) / se_looic, 2)
  ) %>%
  mutate(
    performance = case_when(
      looic_rank == 1 ~ "Best model",
      delta_looic < 2 ~ "Competitive",
      delta_looic < 5 ~ "Worse",
      delta_looic < 10 ~ "Much worse",
      TRUE ~ "Very poor"
    )
  ) %>%
  ungroup() %>%
  arrange(dispersal_type, distance_type, looic_rank)

final_results <- final_results %>%
  select(model_type, dispersal_type, distance_type, variable, model_id,
         elpd_loo, looic, se_looic, p_loo, n_params,
         best_looic, best_elpd, delta_looic, delta_elpd,
         looic_rank, elpd_rank, performance, se_units)

# SAVE RESULTS ================================================================
out_file <- "results/comparison_R3/comprehensive_model_comparison.csv"
write.csv(final_results, out_file, row.names = FALSE)

cat("Comprehensive comparison table created!\n")
cat("Results saved as:", out_file, "\n")
