# ==============================================================================
# BAYESIAN STACKING FOR COMPETITIVE MODELS
# ==============================================================================
#
# DESCRIPTION:
# Implements Bayesian stacking for cases where multiple models show equivalent
# predictive performance (Δ LOOIC ≤ 2). Stacking optimizes leave-one-out 
# predictive performance and is more robust than traditional BMA.
#
# INPUTS REQUIRED:
# - Fitted models (univariate and multivariate)
# - LOO objects for each competitive model
# - loo_competitive_models.csv (from previous analysis)
#
# OUTPUTS:
# - Stacking weights for each competitive model group
# - Weighted predictions and coefficients
# - Summary tables with model weights
#
# DEPENDENCIES:
# library(loo)      # For stacking weights
# library(brms)     # For model objects
# library(dplyr)    # Data manipulation
#
# ==============================================================================

# Load libraries
suppressPackageStartupMessages({
  library(loo)
  library(brms)
  library(dplyr)
  library(tidyr)
})

# LOAD REQUIRED DATA ========================================================

cat("=== LOADING COMPETITIVE MODELS DATA ===\n")

# Load competitive models table
if (!file.exists("results/tables/complete/loo_competitive_models_complete.csv")) {
  stop("Please run LOO analysis first to generate competitive models table")
}

competitive_models <- read.csv("results/tables/complete/loo_competitive_models_complete.csv")
cat("Competitive models loaded:", nrow(competitive_models), "models\n")

# Load fitted models
load("results/weibull/average/univariate_models.RData")
models_average <- models
load("results/weibull/natal/univariate_models.RData")
models_natal <- models
load("results/weibull/breeding/univariate_models.RData")
models_breeding <- models

# Load complete multivariate models
load("results/combined/complete_models_combined.RData")

rm(models)  # Clean up

# FUNCTION TO GET MODEL OBJECT BY ID ========================================

get_model_by_id <- function(dispersal_type, distance_type, variable, model_type) {
  
  if (model_type == "multivariate" || model_type == "complete_multivariate") {
    # Multivariate model
    model_row <- all_complete_models %>%
      filter(age == dispersal_type, dispersal_mode == distance_type)
    
    if (nrow(model_row) > 0) {
      return(model_row$model[[1]])
    }
  } else {
    # Univariate model
    model_name <- paste0(variable, "_", distance_type)
    
    if (dispersal_type == "average") {
      return(models_average[[model_name]])
    } else if (dispersal_type == "natal") {
      return(models_natal[[model_name]])
    } else if (dispersal_type == "breeding") {
      return(models_breeding[[model_name]])
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
    if (model_type == "multivariate" || model_type == "complete_multivariate") {
      model_id <- paste0(dispersal_type, "_", distance_type, "_multivariate")
    } else {
      model_id <- paste0(dispersal_type, "_", distance_type, "_", variable, "_univariate")
    }
    
    model_obj <- get_model_by_id(dispersal_type, distance_type, variable, model_type)
    
    if (!is.null(model_obj)) {
      models_list[[model_id]] <- model_obj
      
      # Calculate LOO for this model
      tryCatch({
        loo_obj <- loo(model_obj, refresh = 0)
        loo_list[[model_id]] <- loo_obj
        cat("  ✓ LOO calculated for", model_id, "\n")
      }, error = function(e) {
        cat("  ✗ Error calculating LOO for", model_id, ":", e$message, "\n")
      })
    } else {
      cat("  ⚠ Could not find model for", dispersal_type, distance_type, variable, "\n")
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
          model_type = ifelse(grepl("multivariate", model_id), "multivariate", "univariate"),
          variable = ifelse(grepl("multivariate", model_id), "full_model", 
                            sapply(strsplit(model_id, "_"), function(x) paste(x[3:(length(x)-1)], collapse = "_")))
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

# FUNCTION TO CREATE WEIGHTED PREDICTIONS ===================================

create_weighted_predictions <- function(group_data, weights_data, newdata = NULL) {
  
  if (is.null(weights_data) || nrow(weights_data) == 0) {
    return(NULL)
  }
  
  # Use original data if newdata not provided
  if (is.null(newdata)) {
    dispersal_type <- unique(weights_data$dispersal_type)
    if (dispersal_type == "average") {
      load("data/processed/dispersal_average_complete.RData")
      newdata <- data_average
    } else if (dispersal_type == "natal") {
      load("data/processed/dispersal_natal_complete.RData")
      newdata <- data_natal
    } else if (dispersal_type == "breeding") {
      load("data/processed/dispersal_breeding_complete.RData")
      newdata <- data_breeding
    }
  }
  
  weighted_pred <- NULL
  total_weight <- 0
  
  for (i in 1:nrow(weights_data)) {
    dispersal_type <- weights_data$dispersal_type[i]
    distance_type <- weights_data$distance_type[i]
    variable <- weights_data$variable[i]
    model_type <- weights_data$model_type[i]
    weight <- weights_data$stacking_weight[i]
    
    model_obj <- get_model_by_id(dispersal_type, distance_type, variable, model_type)
    
    if (!is.null(model_obj) && weight > 0.01) {  # Only include models with meaningful weight
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
    # Normalize weights if they don't sum to 1
    weighted_pred <- weighted_pred / total_weight
    return(weighted_pred)
  }
  
  return(NULL)
}

# MAIN STACKING ANALYSIS ====================================================

cat("\n=== PERFORMING STACKING ANALYSIS ===\n")

# Group competitive models by dispersal_type and distance_type
competitive_groups <- competitive_models %>%
  group_by(dispersal_type, distance_type) %>%
  filter(n() > 1) %>%  # Only groups with multiple competitive models
  group_split()

cat("Found", length(competitive_groups), "groups with multiple competitive models\n")

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
write.csv(all_stacking_weights, "results/tables/stacking_weights_detailed.csv", row.names = FALSE)

# Create summary table
stacking_summary <- all_stacking_weights %>%
  group_by(dispersal_type, distance_type) %>%
  summarise(
    n_competitive_models = n(),
    best_model = variable[which.max(stacking_weight)],
    best_model_type = model_type[which.max(stacking_weight)],
    best_weight = max(stacking_weight),
    weight_concentration = max(stacking_weight),  # How concentrated weights are
    effective_models = sum(stacking_weight > 0.1),  # Models with >10% weight
    model_uncertainty = 1 - max(stacking_weight),   # 1 - max weight = uncertainty
    .groups = "drop"
  ) %>%
  mutate(
    interpretation = case_when(
      weight_concentration > 0.8 ~ "Strong preference for best model",
      weight_concentration > 0.6 ~ "Moderate preference for best model", 
      weight_concentration > 0.4 ~ "Substantial model uncertainty",
      TRUE ~ "High model uncertainty"
    )
  )

write.csv(stacking_summary, "results/tables/stacking_summary.csv", row.names = FALSE)

# GENERATE WEIGHTED PREDICTIONS =============================================

cat("\n=== GENERATING WEIGHTED PREDICTIONS ===\n")

# Create weighted predictions for each competitive group
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

# RESULTS SUMMARY ===========================================================

cat("\n=== STACKING ANALYSIS RESULTS ===\n")

print("STACKING SUMMARY:")
print(stacking_summary)

cat("\nDETAILED INTERPRETATION:\n")
for (i in 1:nrow(stacking_summary)) {
  row <- stacking_summary[i, ]
  cat(sprintf("• %s %s: %s (best: %s with %.1f%% weight)\n",
              stringr::str_to_title(row$dispersal_type),
              row$distance_type,
              row$interpretation,
              row$best_model,
              row$best_weight * 100))
}

# Model averaging recommendations
cat("\nRECOMMENDATIONS:\n")
high_uncertainty <- stacking_summary %>% filter(model_uncertainty > 0.4)
if (nrow(high_uncertainty) > 0) {
  cat("Groups with high model uncertainty (>40%) - consider model averaging:\n")
  for (i in 1:nrow(high_uncertainty)) {
    cat("• ", high_uncertainty$dispersal_type[i], high_uncertainty$distance_type[i], "\n")
  }
} else {
  cat("All groups show reasonable model certainty (uncertainty ≤ 40%)\n")
}

cat("\n=== FILES GENERATED ===\n")
cat("- results/tables/stacking_weights_detailed.csv (detailed weights for all models)\n")
cat("- results/tables/stacking_summary.csv (summary with interpretation)\n")

cat("\n=== STACKING ANALYSIS COMPLETED ===\n")

# EXAMPLE USAGE IN PAPER ====================================================

cat("\n=== FOR YOUR PAPER ===\n")
cat("Example text for manuscript:\n\n")

example_text <- 'When multiple models showed equivalent predictive performance 
(Δ LOOIC ≤ 2), we used Bayesian stacking to weight predictions based on 
leave-one-out cross-validation performance. Stacking weights reflect each 
model\'s contribution to optimal predictive accuracy while accounting for 
model uncertainty. For [specific case], the multivariate model received 
[X]% weight, with the remaining weight distributed among [Y] competitive 
univariate models, indicating [interpretation].'

cat(example_text)

cat("\n\nInterpretation guide:\n")
cat("• Weight > 80%: Strong evidence for single best model\n")
cat("• Weight 60-80%: Moderate preference with some uncertainty\n") 
cat("• Weight 40-60%: Substantial uncertainty, use averaging\n")
cat("• Weight < 40%: High uncertainty, multiple viable models\n")

# ==============================================================================
# EXTENSIÓN: COEFICIENTES PROMEDIADOS Y VISUALIZACIONES
# ==============================================================================

# FUNCIÓN PARA CALCULAR COEFICIENTES PROMEDIADOS =============================

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
    model_type <- weights_data$model_type[i]
    weight <- weights_data$stacking_weight[i]
    
    model_obj <- get_model_by_id(dispersal_type, distance_type, variable, model_type)
    
    if (!is.null(model_obj) && weight > 0.01) {
      tryCatch({
        # Extraer coeficientes fijos
        coefs <- fixef(model_obj)[, "Estimate"]
        coef_names <- names(coefs)
        
        # Inicializar en primera iteración
        if (is.null(weighted_coefs)) {
          weighted_coefs <- weight * coefs
          coefficient_names <- coef_names
        } else {
          # Alinear coeficientes (algunos modelos pueden tener diferentes variables)
          aligned_coefs <- rep(0, length(coefficient_names))
          names(aligned_coefs) <- coefficient_names
          
          # Añadir coeficientes existentes
          for (name in coef_names) {
            if (name %in% coefficient_names) {
              aligned_coefs[name] <- coefs[name]
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
    # Normalizar por peso total
    weighted_coefs <- weighted_coefs / total_weight
    
    # Crear data frame con resultados
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

# CALCULAR COEFICIENTES PROMEDIADOS PARA TODOS LOS GRUPOS ===================

cat("\n=== CALCULATING WEIGHTED COEFFICIENTS ===\n")

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

# Combinar todos los coeficientes promediados
if (length(all_weighted_coefficients) > 0) {
  combined_weighted_coefficients <- do.call(rbind, all_weighted_coefficients)
  rownames(combined_weighted_coefficients) <- NULL
  
  # Guardar tabla de coeficientes promediados
  write.csv(combined_weighted_coefficients, 
            "results/tables/weighted_coefficients.csv", 
            row.names = FALSE)
  
  cat("Weighted coefficients calculated for", length(all_weighted_coefficients), "groups\n")
} else {
  cat("No weighted coefficients calculated\n")
}

# CREAR VISUALIZACIONES ======================================================

cat("\n=== CREATING STACKING VISUALIZATIONS ===\n")

library(ggplot2)

# 1. GRÁFICO DE PESOS DE STACKING ============================================

if (nrow(all_stacking_weights) > 0) {
  
  # Preparar datos para visualización
  plot_data <- all_stacking_weights %>%
    mutate(
      group_label = paste(stringr::str_to_title(dispersal_type), distance_type),
      model_label = ifelse(model_type == "multivariate", 
                           "Multivariate", 
                           paste("Univariate:", variable)),
      model_label = factor(model_label, levels = unique(model_label))
    )
  
  # Gráfico de barras de pesos de stacking
  p1 <- ggplot(plot_data, aes(x = reorder(model_label, stacking_weight), 
                              y = stacking_weight, 
                              fill = model_type)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_text(aes(label = paste0(round(stacking_weight * 100, 1), "%")), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    facet_wrap(~group_label, scales = "free_y") +
    scale_fill_manual(values = c("multivariate" = "#2E86AB", "univariate" = "#A23B72")) +
    labs(
      title = "Bayesian Stacking Weights for Competitive Models",
      subtitle = "Higher weights indicate better predictive performance",
      x = "Model",
      y = "Stacking Weight",
      fill = "Model Type"
    ) +
    theme_bw() +
    theme(
      axis.text.y = element_text(size = 8, color = "black"),
      axis.text.x = element_text(color = "black"),
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
    ) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1))
  
  ggsave("results/figures/stacking_weights.png", p1, 
         width = 12, height = 8, dpi = 300)
  
  cat("✓ Stacking weights plot saved\n")
  
  # 2. GRÁFICO DE INCERTIDUMBRE DEL MODELO ===================================
  
  uncertainty_data <- stacking_summary %>%
    mutate(
      group_label = paste(stringr::str_to_title(dispersal_type), distance_type),
      uncertainty_level = case_when(
        model_uncertainty <= 0.2 ~ "Low (≤20%)",
        model_uncertainty <= 0.4 ~ "Moderate (20-40%)",
        model_uncertainty <= 0.6 ~ "High (40-60%)",
        TRUE ~ "Very High (>60%)"
      ),
      uncertainty_level = factor(uncertainty_level, 
                                 levels = c("Low (≤20%)", "Moderate (20-40%)", 
                                            "High (40-60%)", "Very High (>60%)"))
    )
  
  p2 <- ggplot(uncertainty_data, aes(x = reorder(group_label, model_uncertainty), 
                                     y = model_uncertainty, 
                                     fill = uncertainty_level)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "red", alpha = 0.7, size = 0.8) +
    geom_text(aes(label = paste0(round(model_uncertainty * 100, 1), "%")), 
              hjust = -0.1, size = 3, color = "black") +
    coord_flip() +
    scale_fill_manual(values = c("Low (≤20%)" = "#2E8B57", 
                                 "Moderate (20-40%)" = "#FFD700", 
                                 "High (40-60%)" = "#FF8C00", 
                                 "Very High (>60%)" = "#DC143C")) +
    labs(
      title = "Model Uncertainty Across Competitive Groups",
      subtitle = "Red line indicates threshold for model averaging (40%)",
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
  
  ggsave("results/figures/model_uncertainty.png", p2, 
         width = 10, height = 6, dpi = 300)
  
  cat("✓ Model uncertainty plot saved\n")
  
  # 3. TABLA VISUAL DE COEFICIENTES PROMEDIADOS ============================
  
  if (exists("combined_weighted_coefficients")) {
    
    coef_plot_data <- combined_weighted_coefficients %>%
      filter(!grepl("Intercept", coefficient)) %>%  # Remover intercepto
      mutate(
        group_label = paste(stringr::str_to_title(dispersal_type), distance_type),
        coefficient_clean = gsub(":", " × ", coefficient),  # Mejor formato para interacciones
        significant = abs(weighted_estimate) > 0.1  # Umbral para "significativo"
      )
    
    if (nrow(coef_plot_data) > 0) {
      p3 <- ggplot(coef_plot_data, aes(x = weighted_estimate, 
                                       y = reorder(coefficient_clean, weighted_estimate),
                                       color = significant)) +
        geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5, color = "black", size = 0.8) +
        geom_point(size = 3, alpha = 0.8) +
        facet_wrap(~group_label, scales = "free") +
        scale_color_manual(values = c("TRUE" = "#2E86AB", "FALSE" = "#CCCCCC")) +
        labs(
          title = "Weighted Model Coefficients from Stacking",
          subtitle = "Coefficients averaged across competitive models by stacking weights",
          x = "Weighted Coefficient Estimate",
          y = "Predictor Variable",
          color = "Notable Effect\n(|coef| > 0.1)"
        ) +
        theme_bw() +
        theme(
          axis.text = element_text(color = "black"),
          axis.title = element_text(color = "black"),
          legend.position = "bottom",
          legend.text = element_text(color = "black"),
          legend.title = element_text(color = "black"),
          strip.text = element_text(face = "bold", color = "black"),
          strip.background = element_rect(fill = "white", color = "black"),
          panel.grid.major = element_line(color = "grey90", size = 0.5),
          panel.grid.minor = element_line(color = "grey95", size = 0.25),
          plot.title = element_text(color = "black"),
          plot.subtitle = element_text(color = "black")
        )
      
      ggsave("results/figures/weighted_coefficients.png", p3, 
             width = 12, height = 8, dpi = 300)
      
      cat("✓ Weighted coefficients plot saved\n")
    }
  }
  
} else {
  cat("⚠ No stacking weights available for visualization\n")
}

# RESUMEN FINAL PARA EL MANUSCRITO ===========================================

cat("\n=== SUMMARY FOR MANUSCRIPT ===\n")

if (nrow(stacking_summary) > 0) {
  cat("\nSTACKING ANALYSIS SUMMARY:\n")
  
  # Estadísticas generales
  high_uncertainty <- sum(stacking_summary$model_uncertainty > 0.4)
  total_groups <- nrow(stacking_summary)
  
  cat("• Total competitive groups analyzed:", total_groups, "\n")
  cat("• Groups with high uncertainty (>40%):", high_uncertainty, "\n")
  cat("• Groups benefiting from model averaging:", 
      round(100 * high_uncertainty / total_groups, 1), "%\n")
  
  # Mejores modelos
  best_models_summary <- stacking_summary %>%
    count(best_model_type, name = "frequency") %>%
    arrange(desc(frequency))
  
  cat("\nMOST FREQUENT BEST MODELS:\n")
  for (i in 1:nrow(best_models_summary)) {
    cat("•", best_models_summary$best_model_type[i], ":", 
        best_models_summary$frequency[i], "times\n")
  }
  
  # Ejemplo de texto para manuscrito
  cat("\nEXAMPLE TEXT FOR MANUSCRIPT:\n")
  example_case <- stacking_summary[1, ]
  cat(sprintf('For %s %s dispersal, stacking assigned %.1f%% weight to the %s model, 
with remaining weight distributed among %d competitive models (model uncertainty = %.1f%%). 
This %s suggests %s.\n\n',
              example_case$dispersal_type,
              example_case$distance_type,
              example_case$best_weight * 100,
              example_case$best_model,
              example_case$n_competitive_models - 1,
              example_case$model_uncertainty * 100,
              tolower(example_case$interpretation),
              ifelse(example_case$model_uncertainty > 0.4, 
                     "model averaging is recommended", 
                     "the best model can be used with confidence")))
}

cat("\n=== FILES GENERATED ===\n")
cat("Tables:\n")
cat("- results/tables/weighted_coefficients.csv\n")
cat("Figures:\n")
cat("- results/figures/stacking_weights.png\n")
cat("- results/figures/model_uncertainty.png\n")
cat("- results/figures/weighted_coefficients.png\n")

cat("\n=== STACKING EXTENSION COMPLETED ===\n")
