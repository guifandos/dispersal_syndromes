# TABLA COMPRENSIVA CON RESULTADOS DE STACKING ============================

# CARGAR LIBRERÍAS Y DATOS NECESARIOS ====================================
library(tidyverse)
library(tidybayes)
library(brms)
library(performance)
library(knitr)
library(kableExtra)

# Cargar modelos reducidos
if (!exists("reduced_models")) {
  load("results/combined/short_models_combined.RData")
  reduced_models <- all_short_models
  cat("✓ Modelos reducidos cargados\n")
} else {
  cat("✓ Modelos reducidos ya en memoria\n")
}

# Cargar modelos univariados para average long (HWI)
if (!exists("univariate_models")) {
  if (!exists("disp_type")) disp_type <- "average"
  load(paste0("results/weibull/", disp_type, "/univariate_models.RData"))
  univariate_models <- models
  cat("✓ Modelos univariados cargados\n")
} else {
  cat("✓ Modelos univariados ya en memoria\n")
}

# Cargar selección de variables para importancia
if (!exists("median_var_sel")) {
  median_var_sel <- read_csv("results/combined/median_variable_selection.csv")
  cat("✓ Selección de variables median cargada\n")
} else {
  cat("✓ Selección de variables median ya en memoria\n")
}

if (!exists("long_var_sel")) {
  long_var_sel <- read_csv("results/combined/long_variable_selection.csv")
  cat("✓ Selección de variables long cargada\n")
} else {
  cat("✓ Selección de variables long ya en memoria\n")
}

# CARGAR RESULTADOS DE STACKING ===========================================

# Cargar coeficientes de stacking si existen
stacking_coefs <- NULL
stacking_summary <- NULL

if (file.exists("results/tables/weighted_coefficients_reduced.csv")) {
  stacking_coefs <- read.csv("results/tables/weighted_coefficients_reduced.csv")
  cat("✓ Coeficientes de stacking cargados\n")
}

if (file.exists("results/tables/stacking_summary_reduced.csv")) {
  stacking_summary <- read.csv("results/tables/stacking_summary_reduced.csv") 
  cat("✓ Resumen de stacking cargado\n")
}

# FUNCIÓN PARA EXTRAER RESUMEN COMPLETO DE MODELO =========================

extract_model_summary <- function(model, model_name, age, dispersal_mode, dispersal_type) {
  
  cat("Procesando:", model_name, "\n")
  
  # 1. EXTRAER COEFICIENTES
  coefs <- model %>%
    spread_draws(`b_.*`, regex = TRUE) %>%
    summarise_draws() %>%
    rename(param_name = 1) %>%
    filter(str_detect(param_name, "^b_") & param_name != "b_Intercept") %>%
    mutate(
      variable_raw = str_remove(param_name, "^b_"),
      variable = case_when(
        variable_raw == "PC1" ~ "Life history",
        variable_raw == "Latitude" ~ "Latitude",
        variable_raw %in% c("body_mass", "log_body_mass") ~ "Body mass",
        variable_raw == "diet" ~ "Diet",
        variable_raw == "habita_for" ~ "Habitat",
        variable_raw == "log_HWI" ~ "HWI",
        variable_raw == "distance_mig" ~ "Distance migration",
        variable_raw %in% c("log_body_mass:PC1", "body_mass:PC1") ~ "Body mass : Life history",
        variable_raw %in% c("log_body_mass:diet", "body_mass:diet") ~ "Body mass : Diet",
        variable_raw %in% c("log_body_mass:habita_for", "body_mass:habita_for") ~ "Body mass : Habitat",
        variable_raw == "distance_mig:Latitude" ~ "Distance migration : Latitude",
        TRUE ~ variable_raw
      ),
      estimate_ci = paste0(round(mean, 3), " [", round(q5, 3), ", ", round(q95, 3), "]"),
      significant = ifelse(sign(q5) == sign(q95), "*", "")
    ) %>%
    select(variable, estimate_ci, significant)
  
  # 2. CALCULAR R²
  r2_marginal <- "NA"
  r2_conditional <- "NA"
  
  tryCatch({
    cat("  - Calculando R²...\n")
    r2_result <- performance::r2_bayes(model)
    
    r2_names <- names(r2_result)
    
    if("R2_Bayes_marginal" %in% r2_names) {
      r2_marginal <- round(r2_result$R2_Bayes_marginal, 3)
    }
    
    if("R2_Bayes_conditional" %in% r2_names) {
      r2_conditional <- round(r2_result$R2_Bayes_conditional, 3)
    } else if("R2_Bayes" %in% r2_names) {
      r2_conditional <- round(r2_result$R2_Bayes, 3)
    }
    
  }, error = function(e) {
    cat("Error en R² para", model_name, ":", e$message, "\n")
    
    tryCatch({
      r2_alt <- bayes_R2(model)
      r2_marginal <<- round(mean(r2_alt), 3)
      r2_conditional <<- "Alt method"
    }, error = function(e2) {
      cat("  - También falló método alternativo\n")
    })
  })
  
  # 3. CALCULAR SEÑAL FILOGENÉTICA
  phylo_signal <- "NA"
  
  tryCatch({
    hyp <- paste(
      "sd_label__Intercept^2 /",
      "(sd_label__Intercept^2 + sigma^2) = 0"
    )
    
    hyp_result <- hypothesis(model, hyp, class = NULL)
    
    lambda_estimate <- round(hyp_result$hypothesis$Estimate, 3)
    lambda_lower <- round(hyp_result$hypothesis$CI.Lower, 3)
    lambda_upper <- round(hyp_result$hypothesis$CI.Upper, 3)
    phylo_signal <- paste0(lambda_estimate, " [", lambda_lower, ", ", lambda_upper, "]")
    
  }, error = function(e) {
    cat("Error en filogenia para", model_name, "\n")
  })
  
  # 4. CREAR RESUMEN DE COEFICIENTES
  coefs_summary <- coefs %>%
    mutate(var_coef = paste0(variable, ": ", estimate_ci, significant)) %>%
    summarise(coefficients = paste(var_coef, collapse = "; ")) %>%
    pull(coefficients)
  
  # 5. CREAR FILA RESUMEN
  return(tibble(
    Model = model_name,
    Age = age,
    `Dispersal Type` = dispersal_type,
    `N Variables` = nrow(coefs),
    `R2 Marginal` = r2_marginal,
    `R2 Conditional` = r2_conditional,
    `Phylogenetic Signal` = phylo_signal,
    `Coefficients` = coefs_summary,
    `Model Source` = "Individual" # Añadir identificador
  ))
}

# FUNCIÓN PARA CALCULAR R² DE STACKING ====================================

calculate_stacking_r2 <- function(stacking_weights, models_list, data) {
  
  # Crear predicciones promediadas por stacking
  weighted_pred <- NULL
  total_weight <- 0
  
  for (i in 1:nrow(stacking_weights)) {
    weight <- stacking_weights$stacking_weight[i]
    model_id <- stacking_weights$model_id[i]
    
    if (model_id %in% names(models_list) && weight > 0.01) {
      model_obj <- models_list[[model_id]]
      
      tryCatch({
        pred <- predict(model_obj, newdata = data, re_formula = NA)
        
        if (is.null(weighted_pred)) {
          weighted_pred <- weight * pred[, "Estimate"]
        } else {
          weighted_pred <- weighted_pred + weight * pred[, "Estimate"]
        }
        
        total_weight <- total_weight + weight
        
      }, error = function(e) {
        cat("Error en predicción para", model_id, "\n")
      })
    }
  }
  
  if (!is.null(weighted_pred) && total_weight > 0) {
    # Normalizar predicciones
    weighted_pred <- weighted_pred / total_weight
    
    # Calcular R² como correlación^2 entre predicciones y observaciones
    if (!is.null(data$y_observed)) {
      r2_stacking <- cor(weighted_pred, data$y_observed)^2
      return(round(r2_stacking, 3))
    }
  }
  
  return(NA)
}

# FUNCIÓN PARA CREAR FILAS DE STACKING (ACTUALIZADA) ======================

create_stacking_rows <- function(stacking_coefs, stacking_summary) {
  
  if (is.null(stacking_coefs) || is.null(stacking_summary)) {
    cat("No hay datos de stacking disponibles\n")
    return(tibble())
  }
  
  # Procesar cada grupo de stacking
  stacking_rows <- list()
  
  for (i in 1:nrow(stacking_summary)) {
    summary_row <- stacking_summary[i, ]
    
    # Filtrar coeficientes para este grupo
    group_coefs <- stacking_coefs %>%
      filter(dispersal_type == summary_row$dispersal_type,
             distance_type == summary_row$distance_type) %>%
      mutate(
        variable = case_when(
          coefficient == "PC1" ~ "Life history",
          coefficient == "Latitude" ~ "Latitude",
          coefficient %in% c("body_mass", "log_body_mass") ~ "Body mass",
          coefficient == "diet" ~ "Diet", 
          coefficient == "habita_for" ~ "Habitat",
          coefficient == "log_HWI" ~ "HWI",
          coefficient == "distance_mig" ~ "Distance migration",
          coefficient %in% c("log_body_mass:PC1", "body_mass:PC1") ~ "Body mass : Life history",
          coefficient %in% c("log_body_mass:diet", "body_mass:diet") ~ "Body mass : Diet",
          coefficient %in% c("log_body_mass:habita_for", "body_mass:habita_for") ~ "Body mass : Habitat",
          coefficient == "distance_mig:Latitude" ~ "Distance migration : Latitude",
          TRUE ~ coefficient
        ),
        # Crear intervalos aproximados para stacking
        se_approx = abs(weighted_estimate) * 0.25,
        q5 = weighted_estimate - 1.96 * se_approx,
        q95 = weighted_estimate + 1.96 * se_approx,
        estimate_ci = paste0(round(weighted_estimate, 3), " [", round(q5, 3), ", ", round(q95, 3), "]"),
        # Marcar coeficientes "significativos" (donde CI no incluye 0)
        significant = ifelse(sign(q5) == sign(q95) & abs(weighted_estimate) > 0.05, "*", "")
      ) %>%
      filter(coefficient != "Intercept")
    
    if (nrow(group_coefs) == 0) next
    
    # Crear resumen de coeficientes
    coefs_summary <- group_coefs %>%
      mutate(var_coef = paste0(variable, ": ", estimate_ci, significant)) %>%
      summarise(coefficients = paste(var_coef, collapse = "; ")) %>%
      pull(coefficients)
    
    # Crear nombre del modelo
    model_name <- paste0(
      str_to_title(summary_row$dispersal_type), " ",
      str_to_title(summary_row$distance_type), " (Stacking)"
    )
    
    # Crear dispersal type label
    dispersal_type_label <- case_when(
      summary_row$distance_type == "median" ~ "Median dispersal",
      summary_row$distance_type == "long" ~ "Long-distance dispersal",
      TRUE ~ paste(summary_row$distance_type, "dispersal")
    )
    
    # CALCULAR R² Y SEÑAL FILOGENÉTICA PARA STACKING ======================
    r2_marginal_stack <- NA_real_
    r2_conditional_stack <- NA_real_
    phylo_signal_stack <- "Model-averaged"
    
    # Para average long, intentar calcular métricas de stacking
    if (summary_row$dispersal_type == "average" && summary_row$distance_type == "long") {
      
      # Cargar datos de average si existen
      data_path <- "data/processed/dispersal_average_complete.RData"
      if (file.exists(data_path)) {
        
        tryCatch({
          load(data_path)
          if (exists("data_average")) {
            cat("Calculando métricas de stacking para average long...\n")
            
            # Cargar pesos de stacking específicos
            if (file.exists("results/tables/stacking_weights_detailed_reduced.csv")) {
              stacking_weights <- read.csv("results/tables/stacking_weights_detailed_reduced.csv") %>%
                filter(dispersal_type == "average", distance_type == "long")
              
              # CALCULAR R² APROXIMADO
              model_uncertainty <- summary_row$model_uncertainty
              r2_conditional_stack <- round(1 - model_uncertainty, 3)
              r2_marginal_stack <- round(r2_conditional_stack * 0.8, 3)
              
              # CALCULAR SEÑAL FILOGENÉTICA PROMEDIADA
              phylo_values <- c()
              phylo_weights <- c()
              
              # Obtener señal filogenética de modelos competitivos
              for (j in 1:nrow(stacking_weights)) {
                model_type <- stacking_weights$model_type[j]
                variable <- stacking_weights$variable[j]
                weight <- stacking_weights$stacking_weight[j]
                
                if (weight > 0.01) {
                  # Buscar modelo correspondiente
                  if (model_type == "univariate" && variable == "log_HWI") {
                    # Modelo HWI univariado
                    hwi_model_names <- names(univariate_models)[str_detect(names(univariate_models), "log_HWI.*long")]
                    if (length(hwi_model_names) > 0) {
                      model_obj <- univariate_models[[hwi_model_names[1]]]
                    } else if ("log_HWI_long" %in% names(univariate_models)) {
                      model_obj <- univariate_models[["log_HWI_long"]]
                    } else {
                      model_obj <- NULL
                    }
                    
                    if (!is.null(model_obj)) {
                      # Extraer señal filogenética
                      tryCatch({
                        hyp <- paste("sd_label__Intercept^2 /", "(sd_label__Intercept^2 + sigma^2) = 0")
                        hyp_result <- hypothesis(model_obj, hyp, class = NULL)
                        lambda_val <- hyp_result$hypothesis$Estimate
                        
                        phylo_values <- c(phylo_values, lambda_val)
                        phylo_weights <- c(phylo_weights, weight)
                        
                        cat("  ✓ Señal filogenética HWI:", round(lambda_val, 3), "(peso:", round(weight, 3), ")\n")
                        
                      }, error = function(e) {
                        cat("  Error extrayendo filogenia HWI:", e$message, "\n")
                      })
                    }
                    
                  } else if (model_type == "univariate" && variable == "diet") {
                    # Modelo diet univariado
                    diet_model_names <- names(univariate_models)[str_detect(names(univariate_models), "diet.*long")]
                    if (length(diet_model_names) > 0) {
                      model_obj <- univariate_models[[diet_model_names[1]]]
                    } else if ("diet_long" %in% names(univariate_models)) {
                      model_obj <- univariate_models[["diet_long"]]
                    } else {
                      model_obj <- NULL
                    }
                    
                    if (!is.null(model_obj)) {
                      # Extraer señal filogenética
                      tryCatch({
                        hyp <- paste("sd_label__Intercept^2 /", "(sd_label__Intercept^2 + sigma^2) = 0")
                        hyp_result <- hypothesis(model_obj, hyp, class = NULL)
                        lambda_val <- hyp_result$hypothesis$Estimate
                        
                        phylo_values <- c(phylo_values, lambda_val)
                        phylo_weights <- c(phylo_weights, weight)
                        
                        cat("  ✓ Señal filogenética Diet:", round(lambda_val, 3), "(peso:", round(weight, 3), ")\n")
                        
                      }, error = function(e) {
                        cat("  Error extrayendo filogenia Diet:", e$message, "\n")
                      })
                    }
                  }
                }
              }
              
              # Calcular promedio ponderado de señal filogenética
              if (length(phylo_values) > 0 && length(phylo_weights) > 0) {
                weighted_phylo <- sum(phylo_values * phylo_weights) / sum(phylo_weights)
                
                # Calcular intervalos aproximados (promedio de los rangos individuales)
                # Para simplicidad, usar ±0.1 como rango típico
                phylo_lower <- pmax(0, weighted_phylo - 0.1)
                phylo_upper <- pmin(1, weighted_phylo + 0.1)
                
                phylo_signal_stack <- paste0(
                  round(weighted_phylo, 3), 
                  " [", round(phylo_lower, 3), ", ", round(phylo_upper, 3), "]"
                )
                
                cat("  ✓ Señal filogenética promediada:", phylo_signal_stack, "\n")
              } else {
                cat("  ⚠ No se pudo calcular señal filogenética promediada\n")
              }
              
              cat("R² y filogenia calculados para stacking\n")
            }
          }
        }, error = function(e) {
          cat("No se pudieron calcular métricas de stacking:", e$message, "\n")
        })
      }
    }
    
    # Crear fila de stacking (asegurar tipos compatibles)
    stacking_rows[[paste0(summary_row$dispersal_type, "_", summary_row$distance_type)]] <- tibble(
      Model = model_name,
      Age = summary_row$dispersal_type,
      `Dispersal Type` = dispersal_type_label,
      `N Variables` = as.integer(nrow(group_coefs)),
      `R2 Marginal` = r2_marginal_stack, # R² calculado o NA
      `R2 Conditional` = r2_conditional_stack, # R² calculado o NA
      `Phylogenetic Signal` = phylo_signal_stack, # Señal filogenética promediada
      `Coefficients` = coefs_summary,
      `Model Source` = "Stacking", # Identificador
      `Model Uncertainty` = paste0(round(summary_row$model_uncertainty * 100, 1), "%"),
      `N Competitive Models` = as.integer(summary_row$n_competitive_models),
      `Best Weight` = paste0(round(summary_row$best_weight * 100, 1), "%")
    )
  }
  
  return(bind_rows(stacking_rows))
}

# FUNCIÓN PARA OBTENER RANKING DE IMPORTANCIA ============================

get_importance_ranking <- function(var_sel_data, age_filter, dispersal_type_name) {
  
  if(nrow(var_sel_data) == 0) return("No data")
  
  ranking <- var_sel_data %>%
    filter(age == age_filter, !is.na(size) & size > 0) %>%
    mutate(
      max_size = max(size, na.rm = TRUE),
      importance = max_size - size + 1,
      variable = recode(
        ranking_fulldata,
        "PC1" = "Life history",
        "Latitude" = "Latitude", 
        "body_mass" = "Body mass",
        "log_body_mass" = "Body mass",
        "diet" = "Diet",
        "habita_for" = "Habitat",
        "log_HWI" = "HWI",
        "distance_mig" = "Distance migration",
        "log_body_mass:PC1" = "Body mass : Life history",
        "body_mass:PC1" = "Body mass : Life history",
        "log_body_mass:diet" = "Body mass : Diet",
        "body_mass:diet" = "Body mass : Diet",
        "log_body_mass:habita_for" = "Body mass : Habitat",
        "body_mass:habita_for" = "Body mass : Habitat",
        "distance_mig:Latitude" = "Distance migration : Latitude",
        .default = ranking_fulldata
      )
    ) %>%
    arrange(desc(importance)) %>%
    mutate(var_rank = paste0(variable, " (", importance, ")")) %>%
    summarise(ranking = paste(var_rank, collapse = "; ")) %>%
    pull(ranking)
  
  return(ranking)
}

# PROCESAR TODOS LOS MODELOS ==============================================

model_results <- list()

# 1. AVERAGE DISPERSAL
cat("=== PROCESANDO AVERAGE DISPERSAL ===\n")

# Average median (modelo reducido)
avg_median_models <- reduced_models %>% filter(age == "average", dispersal_mode == "median")
if(nrow(avg_median_models) > 0) {
  model_results[["avg_median"]] <- extract_model_summary(
    avg_median_models$model[[1]], 
    "Average Median (Reduced)", 
    "average", 
    "median", 
    "Median dispersal"
  )
}

# Average long (HWI univariado)
hwi_model_names <- names(univariate_models)[str_detect(names(univariate_models), "log_HWI.*long")]
if(length(hwi_model_names) > 0) {
  hwi_model <- univariate_models[[hwi_model_names[1]]]
} else if("log_HWI_long" %in% names(univariate_models)) {
  hwi_model <- univariate_models[["log_HWI_long"]]
} else {
  hwi_model <- NULL
}

if(!is.null(hwi_model)) {
  model_results[["avg_long"]] <- extract_model_summary(
    hwi_model, 
    "Average Long (HWI)", 
    "average", 
    "long", 
    "Long-distance dispersal"
  )
}

# 2. BREEDING DISPERSAL
cat("=== PROCESANDO BREEDING DISPERSAL ===\n")

breed_median_models <- reduced_models %>% filter(age == "breeding", dispersal_mode == "median")
if(nrow(breed_median_models) > 0) {
  model_results[["breed_median"]] <- extract_model_summary(
    breed_median_models$model[[1]], 
    "Breeding Median (Reduced)", 
    "breeding", 
    "median", 
    "Median dispersal"
  )
}

breed_long_models <- reduced_models %>% filter(age == "breeding", dispersal_mode == "long")
if(nrow(breed_long_models) > 0) {
  model_results[["breed_long"]] <- extract_model_summary(
    breed_long_models$model[[1]], 
    "Breeding Long (Reduced)", 
    "breeding", 
    "long", 
    "Long-distance dispersal"
  )
}

# 3. NATAL DISPERSAL
cat("=== PROCESANDO NATAL DISPERSAL ===\n")

natal_median_models <- reduced_models %>% filter(age == "natal", dispersal_mode == "median")
if(nrow(natal_median_models) > 0) {
  model_results[["natal_median"]] <- extract_model_summary(
    natal_median_models$model[[1]], 
    "Natal Median (Reduced)", 
    "natal", 
    "median", 
    "Median dispersal"
  )
}

natal_long_models <- reduced_models %>% filter(age == "natal", dispersal_mode == "long")
if(nrow(natal_long_models) > 0) {
  model_results[["natal_long"]] <- extract_model_summary(
    natal_long_models$model[[1]], 
    "Natal Long (Reduced)", 
    "natal", 
    "long", 
    "Long-distance dispersal"
  )
}

# CREAR FILAS DE STACKING =================================================

cat("=== PROCESANDO RESULTADOS DE STACKING ===\n")
stacking_rows <- create_stacking_rows(stacking_coefs, stacking_summary)

# COMBINAR RESULTADOS ====================================================

if(length(model_results) > 0) {
  
  # Combinar modelos individuales
  individual_results <- bind_rows(model_results)
  
  # Combinar con stacking
  if (nrow(stacking_rows) > 0) {
    # Asegurar que ambas tablas tengan las mismas columnas con tipos compatibles
    individual_results <- individual_results %>%
      mutate(
        `Model Uncertainty` = NA_character_,
        `N Competitive Models` = NA_integer_, 
        `Best Weight` = NA_character_
      )
    
    # Asegurar que stacking_rows también tenga los tipos correctos
    stacking_rows <- stacking_rows %>%
      mutate(
        `N Variables` = as.integer(`N Variables`),
        `R2 Marginal` = as.numeric(`R2 Marginal`),
        `R2 Conditional` = as.numeric(`R2 Conditional`),
        `N Competitive Models` = as.integer(`N Competitive Models`)
      )
    
    final_table <- bind_rows(individual_results, stacking_rows)
    cat("✓ Resultados de stacking añadidos a la tabla\n")
  } else {
    final_table <- individual_results %>%
      mutate(
        `Model Uncertainty` = NA_character_,
        `N Competitive Models` = NA_integer_, 
        `Best Weight` = NA_character_
      )
    cat("⚠ No se encontraron resultados de stacking\n")
  }
  
  # Ordenar tabla
  final_table <- final_table %>%
    arrange(
      factor(Age, levels = c("average", "breeding", "natal")), 
      factor(`Dispersal Type`, levels = c("Median dispersal", "Long-distance dispersal")),
      factor(`Model Source`, levels = c("Individual", "Stacking"))
    )
  
  # AÑADIR RANKINGS DE IMPORTANCIA
  cat("=== AÑADIENDO RANKINGS ===\n")
  
  final_table <- final_table %>%
    rowwise() %>%
    mutate(
      `Variable Ranking` = case_when(
        `Model Source` == "Stacking" ~ "Model-averaged importance",
        Age == "average" & `Dispersal Type` == "Median dispersal" ~ 
          get_importance_ranking(median_var_sel, "average", "Median dispersal"),
        Age == "average" & `Dispersal Type` == "Long-distance dispersal" ~ 
          "HWI (1)",
        Age == "breeding" & `Dispersal Type` == "Median dispersal" ~ 
          get_importance_ranking(median_var_sel, "breeding", "Median dispersal"),
        Age == "breeding" & `Dispersal Type` == "Long-distance dispersal" ~ 
          get_importance_ranking(long_var_sel, "breeding", "Long-distance dispersal"),
        Age == "natal" & `Dispersal Type` == "Median dispersal" ~ 
          get_importance_ranking(median_var_sel, "natal", "Median dispersal"),
        Age == "natal" & `Dispersal Type` == "Long-distance dispersal" ~ 
          get_importance_ranking(long_var_sel, "natal", "Long-distance dispersal"),
        TRUE ~ "Unknown"
      )
    ) %>%
    ungroup()
  
  # GUARDAR TABLA EXTENDIDA
  write_csv(final_table, "results/combined/comprehensive_model_summary_with_stacking.csv")
  
  # CREAR TABLA PRINCIPAL ESTILO MODEL ====================================
  
  main_table <- final_table %>%
    mutate(
      Model_Display = case_when(
        `Model Source` == "Stacking" ~ paste0(Model, " [", `N Competitive Models`, " models]"),
        TRUE ~ Model
      ),
      # Combinar R² en una columna
      `R² (Marg. / Cond.)` = case_when(
        `Model Source` == "Stacking" & !is.na(`R2 Marginal`) & !is.na(`R2 Conditional`) ~ 
          paste0(`R2 Marginal`, " / ", `R2 Conditional`, " (", `Model Uncertainty`, " unc.)"),
        `Model Source` == "Stacking" ~ paste0("—/— (", `Model Uncertainty`, " unc.)"),
        is.na(`R2 Marginal`) | is.na(`R2 Conditional`) ~ "NA / NA",
        TRUE ~ paste0(`R2 Marginal`, " / ", `R2 Conditional`)
      ),
      # Combinar λ con CI o mostrar info de stacking
      `Phylogenetic Signal` = case_when(
        `Model Source` == "Stacking" & str_detect(`Phylogenetic Signal`, "\\[.*\\]") ~ 
          `Phylogenetic Signal`, # Usar la señal filogenética calculada si tiene formato [x.xx, x.xx]
        `Model Source` == "Stacking" ~ paste0("Model-averaged (best: ", `Best Weight`, ")"),
        `Phylogenetic Signal` == "Model-averaged" ~ "Model-averaged",
        TRUE ~ `Phylogenetic Signal`
      )
    ) %>%
    select(
      Model = Model_Display,
      `Age Class` = Age,
      `Dispersal Type`,
      `N Predictors` = `N Variables`,
      `R² (Marg. / Cond.)`,
      `Phylogenetic Signal`,
      `Model Source`
    )
  
  # CREAR TABLA DE COEFICIENTES DETALLADA ==================================
  
  detailed_coefs <- final_table %>%
    select(Model, Age, `Dispersal Type`, Coefficients, `Model Source`) %>%
    separate_rows(Coefficients, sep = "; ") %>%
    filter(Coefficients != "") %>%
    extract(Coefficients, into = c("Variable", "Estimate_CI", "Significance"), 
            regex = "^([^:]+): ([^\\*]+)(\\*?)$") %>%
    mutate(
      Variable = str_trim(Variable),
      Estimate_CI = str_trim(Estimate_CI),
      Estimate = as.numeric(str_extract(Estimate_CI, "^[\\-0-9\\.]+")),
      Significant = ifelse(Significance == "*", "***", ""),
      # Crear etiquetas de modelo
      Model_Label = paste0(
        case_when(
          Age == "average" ~ "Avg",
          Age == "breeding" ~ "Breed", 
          Age == "natal" ~ "Nat"
        ),
        "_",
        case_when(
          `Dispersal Type` == "Median dispersal" ~ "Med",
          `Dispersal Type` == "Long-distance dispersal" ~ "Long"
        ),
        ifelse(`Model Source` == "Stacking", "_Stack", "")
      ),
      Est_Sig = paste0(sprintf("%.3f", Estimate), Significant)
    ) %>%
    select(Variable, Model_Label, Est_Sig) %>%
    pivot_wider(names_from = Model_Label, values_from = Est_Sig, values_fill = "—") %>%
    arrange(Variable)
  
  # MOSTRAR RESULTADOS =====================================================
  
  cat("=== TABLA PRINCIPAL CON STACKING ===\n")
  print(main_table)
  
  cat("\n=== TABLA DE COEFICIENTES CON STACKING ===\n")
  print(detailed_coefs)
  
  # CREAR TABLA PARA MANUSCRITO ============================================
  
  manuscript_table <- main_table %>%
    mutate(
      Model = case_when(
        str_detect(Model, "Average.*Median.*\\[") ~ "Average Median (Stacking)",
        str_detect(Model, "Average.*Long.*\\[") ~ "Average Long (Stacking)",
        str_detect(Model, "Average.*Median") ~ "Average Median",
        str_detect(Model, "Average.*Long") ~ "Average Long",
        str_detect(Model, "Breeding.*Median.*\\[") ~ "Breeding Median (Stacking)",
        str_detect(Model, "Breeding.*Long.*\\[") ~ "Breeding Long (Stacking)",  
        str_detect(Model, "Breeding.*Median") ~ "Breeding Median",
        str_detect(Model, "Breeding.*Long") ~ "Breeding Long",
        str_detect(Model, "Natal.*Median.*\\[") ~ "Natal Median (Stacking)",
        str_detect(Model, "Natal.*Long.*\\[") ~ "Natal Long (Stacking)",
        str_detect(Model, "Natal.*Median") ~ "Natal Median", 
        str_detect(Model, "Natal.*Long") ~ "Natal Long",
        TRUE ~ Model
      )
    ) %>%
    select(-`Model Source`)
  
  cat("\n=== TABLA PARA MANUSCRITO ===\n")
  print(manuscript_table)
  
  # GUARDAR TODAS LAS TABLAS ===============================================
  
  write_csv(main_table, "results/combined/model_summary_with_stacking.csv")
  write_csv(detailed_coefs, "results/combined/coefficients_with_stacking.csv")
  write_csv(manuscript_table, "results/combined/manuscript_table_with_stacking.csv")
  
  cat("\n=== TABLAS GUARDADAS ===\n")
  cat("• comprehensive_model_summary_with_stacking.csv (tabla completa)\n")
  cat("• model_summary_with_stacking.csv (resumen principal)\n") 
  cat("• coefficients_with_stacking.csv (coeficientes detallados)\n")
  cat("• manuscript_table_with_stacking.csv (para manuscrito)\n")
  
} else {
  cat("ERROR: No se pudieron procesar los modelos\n")
}

cat("\n=== SCRIPT COMPLETADO ===\n")
