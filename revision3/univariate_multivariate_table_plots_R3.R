# ==============================================================================
# ANÁLISIS COMPLETO DE MODELOS UNIVARIANTES PARA DISPERSIÓN EN AVES
# ==============================================================================
#
# DESCRIPCIÓN:
# Script completo para correr modelos univariantes de dispersión para tres tipos:
# average (total), natal, y breeding dispersal. Genera tablas comparativas y 
# gráficos de predicciones.
#
# FUNCIONES PRINCIPALES:
# 1. run_univariate_models(): Corre modelos Bayesianos para cada variable
# 2. extract_univariate_coefficients(): Extrae coeficientes y estadísticos
# 3. plot_dispersal_predictions(): Genera scatterplots con predicciones
# 4. create_comparison_tables(): Crea tablas comparativas
#
# INPUTS REQUERIDOS:
# - data/processed/dispersal_average_complete.RData
# - data/processed/dispersal_natal_complete.RData  
# - data/processed/dispersal_breeding_complete.RData
#
# OUTPUTS GENERADOS:
# MODELOS:
# - results/weibull/average/univariate_models.RData
# - results/weibull/natal/univariate_models.RData
# - results/weibull/breeding/univariate_models.RData
#
# TABLAS:
# - results/tables/univariate_vs_multivariate_comparison.csv (completa)
# - results/tables/simple_univariate_vs_multivariate.csv (simplificada)
# - results/tables/direct_univariate_vs_multivariate.csv (comparación directa)
# - results/tables/significant_all_coefficients.csv (solo significativos)
# - results/tables/comprehensive_summary_stats.csv (estadísticos resumen)
#
# GRÁFICOS:
# - results/figures/univariate_plots_average.png
# - results/figures/univariate_plots_natal.png  
# - results/figures/univariate_plots_breeding.png
# - results/figures/univariate_comparison_grid.png
#
# DEPENDENCIAS:
# library(brms)         # Modelos Bayesianos
# library(dplyr)        # Manipulación de datos
# library(tidyr)        # Pivoteo de datos
# library(ggplot2)      # Gráficos
# library(patchwork)    # Combinación de plots
# library(kableExtra)   # Tablas formateadas
# library(tidybayes)    # Para extraer coeficientes multivariantes
# library(posterior)    # Para summarise_draws
# library(stringr)      # Para manipulación de strings
#
# VARIABLES ANALIZADAS:
# body_mass, log_HWI, habita_for, PC1, diet, distance_mig, Latitude
#
# TIEMPO ESTIMADO: ~2-4 horas (dependiendo del hardware)
#
# ==============================================================================

# Cargar librerías
suppressPackageStartupMessages({
  library(brms)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(kableExtra)
  library(tidybayes)
  library(posterior) 
  library(stringr)
})

# Crear directorios
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/weibull/natal", recursive = TRUE, showWarnings = FALSE)
dir.create("results/weibull/breeding", recursive = TRUE, showWarnings = FALSE)

# FUNCIÓN PRINCIPAL PARA MODELOS UNIVARIANTES =============================

run_univariate_models <- function(data, variables, A, dispersal_type = "average", save_path) {
  
  cat("\n=== CORRIENDO MODELOS UNIVARIANTES PARA", toupper(dispersal_type), "===\n")
  cat("Variables a analizar:", paste(variables, collapse = ", "), "\n")
  cat("N especies:", nrow(data), "\n\n")
  
  # Store the models in a list
  models <- list()
  
  # Loop over each variable and fit both models
  for (i in seq_along(variables)) {
    var <- variables[i]
    cat(sprintf("[%d/%d] Ajustando modelos para: %s\n", i, length(variables), var))
    
    # Formulas
    formula_median <- paste0("Weibull_median_log ~ ", var, " + (1|gr(label, cov = A))")
    formula_long <- paste0("Weibull_upper_distance_log ~ ", var, " + (1|gr(label, cov = A))")
    
    # Fit median distance model
    cat("  - Modelo median dispersal...")
    model_median <- brm(
      formula = formula_median,
      data = data, family = gaussian(),
      data2 = list(A = A),
      chains = 2, cores = 2, iter = 4000, warmup = 2000,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      save_pars = save_pars(all = TRUE),
      refresh = 0, silent = 2
    )
    cat(" ✓\n")
    
    # Fit long-distance model  
    cat("  - Modelo long-distance dispersal...")
    model_long <- brm(
      formula = formula_long,
      data = data, family = gaussian(),
      data2 = list(A = A),
      chains = 2, cores = 2, iter = 4000, warmup = 2000,
      control = list(adapt_delta = 0.95, max_treedepth = 12),
      save_pars = save_pars(all = TRUE),
      refresh = 0, silent = 2
    )
    cat(" ✓\n")
    
    # Save models
    models[[paste0(var, "_median")]] <- model_median
    models[[paste0(var, "_long")]] <- model_long
  }
  
  # Save all models
  save(models, file = save_path)
  cat("Modelos guardados en:", save_path, "\n")
  
  return(models)
}

# FUNCIÓN PARA EXTRAER COEFICIENTES UNIVARIANTES =========================

extract_univariate_coefficients <- function(models_list, dispersal_type) {
  
  results <- data.frame()
  
  for (model_name in names(models_list)) {
    model <- models_list[[model_name]]
    
    # Extract variable name and distance type
    parts <- strsplit(model_name, "_")[[1]]
    distance_type <- parts[length(parts)]  
    variable <- paste(parts[-length(parts)], collapse = "_")
    
    # Get model summary
    model_summary <- summary(model)
    coef_data <- model_summary$fixed
    
    # Extract coefficient for the variable (skip intercept)
    if (nrow(coef_data) >= 2) {
      var_row <- 2  # Second row is the variable
      
      result_row <- data.frame(
        dispersal_type = dispersal_type,
        distance_type = distance_type,
        variable = variable,
        estimate = coef_data[var_row, "Estimate"],
        se = coef_data[var_row, "Est.Error"],
        lower_ci = coef_data[var_row, "l-95% CI"],
        upper_ci = coef_data[var_row, "u-95% CI"],
        rhat = coef_data[var_row, "Rhat"],
        significant = !(coef_data[var_row, "l-95% CI"] <= 0 & coef_data[var_row, "u-95% CI"] >= 0),
        model_r2 = round(bayes_R2(model)[1, "Estimate"], 3)
      )
      
      results <- rbind(results, result_row)
    }
  }
  
  return(results)
}

# FUNCIÓN PARA EXTRAER COEFICIENTES MULTIVARIANTES ======================

extract_multivariate_coefficients_from_combined <- function() {
  
  cat("=== CARGANDO MODELOS MULTIVARIANTES DESDE ARCHIVO COMBINADO ===\n")
  
  # Cargar archivo combinado
  if (file.exists("results/combined/short_models_combined.RData")) {
    load("results/combined/short_models_combined.RData")
    cat("Archivo cargado exitosamente\n")
  } else {
    cat("ERROR: No se encuentra results/combined/short_models_combined.RData\n")
    return(data.frame())
  }
  
  # Función para extraer coeficientes de un modelo brms
  extract_model_coefs <- function(model) {
    if (!"tidybayes" %in% loadedNamespaces()) {
      suppressPackageStartupMessages(library(tidybayes))
    }
    if (!"posterior" %in% loadedNamespaces()) {
      suppressPackageStartupMessages(library(posterior))
    }
    if (!"stringr" %in% loadedNamespaces()) {
      suppressPackageStartupMessages(library(stringr))
    }
    
    model %>%
      spread_draws(`b_.*`, regex = TRUE) %>%
      summarise_draws() %>%
      rename(param_name = 1) %>%
      filter(str_detect(param_name, "^b_") & param_name != "b_Intercept") %>%
      mutate(
        variable_raw = str_remove(param_name, "^b_"),
        variable = case_when(
          variable_raw == "PC1" ~ "PC1",  # Mantener nombres consistentes
          variable_raw == "Latitude" ~ "Latitude",
          variable_raw %in% c("body_mass", "log_body_mass") ~ "body_mass",
          variable_raw == "diet" ~ "diet",
          variable_raw == "habita_for" ~ "habita_for",
          variable_raw == "log_HWI" ~ "log_HWI",
          variable_raw == "distance_mig" ~ "distance_mig",
          variable_raw %in% c("log_body_mass:PC1", "body_mass:PC1") ~ "body_mass_x_PC1",
          variable_raw %in% c("log_body_mass:diet", "body_mass:diet") ~ "body_mass_x_diet",
          variable_raw %in% c("log_body_mass:habita_for", "body_mass:habita_for") ~ "body_mass_x_habita_for",
          variable_raw == "distance_mig:Latitude" ~ "distance_mig_x_Latitude",
          TRUE ~ variable_raw
        )
      ) %>%
      select(variable, mean, q5, q95)
  }
  
  # Extraer coeficientes para todos los modelos
  model_coefs <- all_short_models %>%
    rowwise() %>%
    mutate(
      coefs = list(extract_model_coefs(model))
    ) %>%
    unnest(coefs, names_sep = "_") %>%
    select(-model, -variable_selection, -variable_selection_old, -function_t, -function_id) %>%
    select(
      dispersal_type = age, 
      distance_type = dispersal_mode,
      variable = coefs_variable,
      estimate = coefs_mean,
      lower_ci = coefs_q5,
      upper_ci = coefs_q95
    ) %>%
    mutate(
      # Calcular SE aproximado desde intervalos de credibilidad
      se = (upper_ci - lower_ci) / 3.92,  # 3.92 ≈ 2 * 1.96 para 90% CI
      rhat = NA,  # No disponible en este formato
      significant = !(lower_ci <= 0 & upper_ci >= 0),
      model_r2 = NA,  # Calcular después si es necesario
      model_type = "multivariate"
    )
  
  cat("Coeficientes multivariantes extraídos:", nrow(model_coefs), "\n")
  
  return(model_coefs)
}

# FUNCIÓN PARA CALCULAR R2 DE MODELOS MULTIVARIANTES ====================

calculate_multivariate_r2 <- function(multivariate_results) {
  
  if (!"all_short_models" %in% ls(envir = .GlobalEnv)) {
    return(multivariate_results)
  }
  
  r2_values <- data.frame()
  
  for (i in 1:nrow(all_short_models)) {
    model_info <- all_short_models[i, ]
    model <- model_info$model[[1]]
    
    r2_est <- try(bayes_R2(model)[1, "Estimate"], silent = TRUE)
    if (class(r2_est) != "try-error") {
      r2_row <- data.frame(
        dispersal_type = model_info$age,
        distance_type = model_info$dispersal_mode,
        model_r2 = round(r2_est, 3)
      )
      r2_values <- rbind(r2_values, r2_row)
    }
  }
  
  if (nrow(r2_values) > 0) {
    multivariate_results <- multivariate_results %>%
      select(-model_r2) %>%
      left_join(r2_values, by = c("dispersal_type", "distance_type"))
    
    # Rellenar NAs con el R2 del modelo correspondiente
    multivariate_results$model_r2[is.na(multivariate_results$model_r2)] <- 
      rep(r2_values$model_r2, times = table(paste(multivariate_results$dispersal_type, multivariate_results$distance_type)))[is.na(multivariate_results$model_r2)]
  }
  
  return(multivariate_results)
}

# FUNCIÓN PARA PLOTS DE PREDICCIONES ====================================

plot_dispersal_predictions <- function(model_median, model_long, data, variable_name, 
                                       x_axis_label, dispersal_type = "average") {
  
  # Predicciones
  pred_median <- predict(model_median, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% 
    as.data.frame()
  pred_long <- predict(model_long, newdata = data, re_formula = NA, probs = c(0.025, 0.975)) %>% 
    as.data.frame()
  
  # Añadir predicciones a data
  data$predicted_median <- pred_median$Estimate
  data$predicted_long <- pred_long$Estimate
  
  # Plot
  p <- ggplot(data = data, aes_string(x = variable_name)) +
    # Puntos observados
    geom_point(aes(y = Weibull_median_log, color = "Median Dispersal"),
               size = 1.5, alpha = 0.7) +
    geom_point(aes(y = Weibull_upper_distance_log, color = "Long-Distance Dispersal"),
               size = 1.5, alpha = 0.7, shape = 17) +
    # Líneas de predicción
    geom_line(aes(y = predicted_median, color = "Median Dispersal"), size = 1.2) +
    geom_line(aes(y = predicted_long, color = "Long-Distance Dispersal"), 
              size = 1.2, linetype = "dashed") +
    
    # Personalización
    scale_color_manual(values = c("Median Dispersal" = "#ff7f0e", 
                                  "Long-Distance Dispersal" = "#1f77b4")) +
    scale_x_continuous(x_axis_label) +
    ylab("Dispersal distance [log]") +
    #ggtitle(paste("Univariate model:", x_axis_label, "-", stringr::str_to_title(dispersal_type))) +
    theme_classic() +
    theme(
      axis.text = element_text(size = 11), 
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 13, hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank()
    )
  
  return(p)
}

# CARGAR DATOS PREPARADOS ===============================================

cat("=== CARGANDO DATOS PREPARADOS ===\n")

# Average
load("data/processed/dispersal_average_complete.RData")
cat("Average dispersal loaded:", nrow(data_average), "especies\n")

# Natal  
load("data/processed/dispersal_natal_complete.RData")
cat("Natal dispersal loaded:", nrow(data_natal), "especies\n")

# Breeding
load("data/processed/dispersal_breeding_complete.RData") 
cat("Breeding dispersal loaded:", nrow(data_breeding), "especies\n")

# Variables a analizar
variables <- c("body_mass", "log_HWI", "habita_for", "PC1", "diet", "distance_mig", "Latitude")

# CORRER MODELOS UNIVARIANTES ===========================================

# Average
models_average <- run_univariate_models(
  data = data_average, 
  variables = variables, 
  A = A_average, 
  dispersal_type = "average",
  save_path = "results/weibull/average/univariate_models.RData"
)
load("results/weibull/average/univariate_models.RData")
models_average <- models
# Natal
models_natal <- run_univariate_models(
  data = data_natal, 
  variables = variables, 
  A = A_natal, 
  dispersal_type = "natal",
  save_path = "results/weibull/natal/univariate_models.RData"
)
load("results/weibull/natal/univariate_models.RData")
models_natal <- models
# Breeding
models_breeding <- run_univariate_models(
  data = data_breeding, 
  variables = variables, 
  A = A_breeding, 
  dispersal_type = "breeding",
  save_path = "results/weibull/breeding/univariate_models.RData"
)
load("results/weibull/breeding/univariate_models.RData")
models_breeding <- models
# EXTRAER COEFICIENTES ==================================================

cat("\n=== EXTRAYENDO COEFICIENTES ===\n")

# Coeficientes univariantes
coef_average <- extract_univariate_coefficients(models_average, "average")
coef_natal <- extract_univariate_coefficients(models_natal, "natal") 
coef_breeding <- extract_univariate_coefficients(models_breeding, "breeding")

# Añadir columna model_type
coef_average$model_type <- "univariate"
coef_natal$model_type <- "univariate"
coef_breeding$model_type <- "univariate"

# Coeficientes multivariantes desde archivo combinado
multivariate_results <- extract_multivariate_coefficients_from_combined()

# Calcular R2 para modelos multivariantes si es posible
if (nrow(multivariate_results) > 0) {
  multivariate_results <- calculate_multivariate_r2(multivariate_results)
}

# Combinar todos los resultados
all_univariate_results <- rbind(coef_average, coef_natal, coef_breeding)

if (nrow(multivariate_results) > 0) {
  all_results <- rbind(all_univariate_results, multivariate_results)
  cat("Coeficientes multivariantes añadidos:", nrow(multivariate_results), "\n")
} else {
  all_results <- all_univariate_results
  cat("ADVERTENCIA: No se encontraron modelos multivariantes. Solo se usarán univariantes.\n")
}

# CREAR TABLAS COMPARATIVAS ============================================

cat("=== CREANDO TABLAS COMPARATIVAS ===\n")

# FUNCIÓN PARA CREAR TABLA UNIVARIANTE (R² por variable) ================

create_univariate_table <- function(all_results) {
  
  cat("Creando tabla univariante...\n")
  
  univariate_data <- all_results %>%
    filter(model_type == "univariate") %>%
    mutate(
      # Etiquetas claras
      disp_label = case_when(
        dispersal_type == "average" ~ "Total",
        dispersal_type == "natal" ~ "Natal", 
        dispersal_type == "breeding" ~ "Breeding"
      ),
      dist_label = case_when(
        distance_type == "median" ~ "Median",
        distance_type == "long" ~ "Long-distance"
      ),
      variable_clean = case_when(
        variable == "body_mass" ~ "Body mass",
        variable == "log_HWI" ~ "HWI",
        variable == "PC1" ~ "Life history",
        variable == "habita_for" ~ "Habitat openness",
        variable == "diet" ~ "Diet",
        variable == "distance_mig" ~ "Migration distance",
        variable == "Latitude" ~ "Latitude",
        TRUE ~ variable
      ),
      # Formato: coeficiente (IC) [R²] con *
      coef_formatted = paste0(
        sprintf("%.3f", estimate),
        " (", sprintf("%.3f", lower_ci), ", ", sprintf("%.3f", upper_ci), ")",
        ifelse(!is.na(model_r2), paste0(" [R² = ", sprintf("%.3f", model_r2), "]"), ""),
        ifelse(significant, "*", "")
      ),
      # Crear columna header
      col_header = paste(disp_label, dist_label, sep = " - ")
    )
  
  # Crear tabla pivotada
  univariate_table <- univariate_data %>%
    select(variable_clean, col_header, coef_formatted) %>%
    pivot_wider(
      names_from = col_header,
      values_from = coef_formatted
    ) %>%
    arrange(variable_clean)
  
  return(univariate_table)
}

# FUNCIÓN PARA CREAR TABLA MULTIVARIANTE (R² por modelo) ================

create_multivariate_table <- function(all_results) {
  
  cat("Creando tabla multivariante...\n")
  
  multivariate_data <- all_results %>%
    filter(model_type == "multivariate") %>%
    mutate(
      # Etiquetas claras
      disp_label = case_when(
        dispersal_type == "average" ~ "Total",
        dispersal_type == "natal" ~ "Natal", 
        dispersal_type == "breeding" ~ "Breeding"
      ),
      dist_label = case_when(
        distance_type == "median" ~ "Median",
        distance_type == "long" ~ "Long-distance"
      ),
      variable_clean = case_when(
        variable == "body_mass" ~ "Body mass",
        variable == "log_HWI" ~ "HWI",
        variable == "PC1" ~ "Life history",
        variable == "habita_for" ~ "Habitat openness",
        variable == "diet" ~ "Diet",
        variable == "distance_mig" ~ "Migration distance",
        variable == "Latitude" ~ "Latitude",
        variable == "body_mass_x_PC1" ~ "Body mass × Life history",
        variable == "body_mass_x_diet" ~ "Body mass × Diet", 
        variable == "body_mass_x_habita_for" ~ "Body mass × Habitat",
        variable == "distance_mig_x_Latitude" ~ "Migration × Latitude",
        TRUE ~ variable
      ),
      # Formato: solo coeficiente (IC) con *
      coef_formatted = paste0(
        sprintf("%.3f", estimate),
        " (", sprintf("%.3f", lower_ci), ", ", sprintf("%.3f", upper_ci), ")",
        ifelse(significant, "*", "")
      ),
      # Crear columna header
      col_header = paste(disp_label, dist_label, sep = " - ")
    )
  
  # Crear tabla pivotada
  multivariate_table <- multivariate_data %>%
    select(variable_clean, col_header, coef_formatted) %>%
    pivot_wider(
      names_from = col_header,
      values_from = coef_formatted
    ) %>%
    arrange(variable_clean)
  
  # Añadir R² del modelo como nota al pie
  r2_info <- all_results %>%
    filter(model_type == "multivariate", !is.na(model_r2)) %>%
    distinct(dispersal_type, distance_type, model_r2) %>%
    mutate(
      disp_label = case_when(
        dispersal_type == "average" ~ "Total",
        dispersal_type == "natal" ~ "Natal", 
        dispersal_type == "breeding" ~ "Breeding"
      ),
      dist_label = case_when(
        distance_type == "median" ~ "Median",
        distance_type == "long" ~ "Long-distance"
      ),
      model_info = paste0(disp_label, " - ", dist_label, ": R² = ", sprintf("%.3f", model_r2))
    )
  
  # Crear nota de R² 
  r2_note <- paste("Model R²:", paste(r2_info$model_info, collapse = "; "))
  
  # Añadir nota como atributo de la tabla
  attr(multivariate_table, "r2_note") <- r2_note
  
  return(multivariate_table)
}

# GENERAR TABLAS SEPARADAS ==============================================

# Verificar que tenemos ambos tipos de modelos
has_univariate <- any(all_results$model_type == "univariate")
has_multivariate <- any(all_results$model_type == "multivariate")

if (has_univariate) {
  # 1. Tabla univariante (R² por variable)
  univariate_table <- create_univariate_table(all_results)
  write.csv(univariate_table, "results/tables/univariate_models_table.csv", row.names = FALSE)
  
  cat("\n=== TABLA UNIVARIANTE CREADA ===\n")
  cat("Formato: Coeficiente (IC 95%) [R² individual] *\n")
  cat("* = Significativo\n")
} else {
  cat("No se encontraron modelos univariantes\n")
}

if (has_multivariate) {
  # 2. Tabla multivariante (R² del modelo como nota)
  multivariate_table <- create_multivariate_table(all_results)
  write.csv(multivariate_table, "results/tables/multivariate_models_table.csv", row.names = FALSE)
  
  # Guardar nota de R² en archivo separado
  r2_note <- attr(multivariate_table, "r2_note")
  if (!is.null(r2_note)) {
    writeLines(r2_note, "results/tables/multivariate_r2_note.txt")
    
    cat("\n=== TABLA MULTIVARIANTE CREADA ===\n")
    cat("Formato: Coeficiente (IC 95%) *\n")
    cat("* = Significativo\n")
    cat("R² del modelo:", r2_note, "\n")
  }
} else {
  cat("No se encontraron modelos multivariantes\n")
}

# FUNCIÓN PARA CREAR TABLA DE SOLO SIGNIFICATIVOS =======================

create_significant_only_table <- function(all_results) {
  
  cat("Creando tabla de significativos...\n")
  
  significant_data <- all_results %>%
    filter(significant == TRUE) %>%
    mutate(
      disp_label = case_when(
        dispersal_type == "average" ~ "Total",
        dispersal_type == "natal" ~ "Natal", 
        dispersal_type == "breeding" ~ "Breeding"
      ),
      dist_label = case_when(
        distance_type == "median" ~ "Median",
        distance_type == "long" ~ "Long-distance"
      ),
      model_label = case_when(
        model_type == "univariate" ~ "Uni",
        model_type == "multivariate" ~ "Multi"
      ),
      variable_clean = case_when(
        variable == "body_mass" ~ "Body mass",
        variable == "log_HWI" ~ "HWI",
        variable == "PC1" ~ "Life history",
        variable == "habita_for" ~ "Habitat openness",
        variable == "diet" ~ "Diet",
        variable == "distance_mig" ~ "Migration distance",
        variable == "Latitude" ~ "Latitude",
        TRUE ~ stringr::str_replace(variable, "_x_", " × ")
      ),
      model_combo = paste(disp_label, dist_label, model_label, sep = " - "),
      coef_clean = sprintf("%.3f [%.3f, %.3f]", estimate, lower_ci, upper_ci)
    ) %>%
    select(variable_clean, model_combo, coef_clean) %>%
    pivot_wider(names_from = model_combo, values_from = coef_clean) %>%
    arrange(variable_clean)
  
  return(significant_data)
}

# FUNCIÓN PARA RESUMEN DE R² POR MODELO =================================

create_r2_summary_table <- function(all_results) {
  
  cat("Creando resumen de R²...\n")
  
  r2_summary <- all_results %>%
    filter(!is.na(model_r2)) %>%
    distinct(dispersal_type, distance_type, model_type, model_r2) %>%
    mutate(
      disp_label = case_when(
        dispersal_type == "average" ~ "Total",
        dispersal_type == "natal" ~ "Natal", 
        dispersal_type == "breeding" ~ "Breeding"
      ),
      dist_label = case_when(
        distance_type == "median" ~ "Median",
        distance_type == "long" ~ "Long-distance"
      ),
      model_label = case_when(
        model_type == "univariate" ~ "Univariate",
        model_type == "multivariate" ~ "Multivariate"
      ),
      # Interpretación del R²
      r2_interpretation = case_when(
        model_r2 < 0.1 ~ "Small effect",
        model_r2 < 0.3 ~ "Moderate effect",
        model_r2 < 0.5 ~ "Large effect",
        model_r2 >= 0.5 ~ "Very large effect",
        TRUE ~ "Unknown"
      )
    ) %>%
    arrange(disp_label, dist_label, model_label) %>%
    select(`Dispersal Type` = disp_label,
           `Distance Type` = dist_label, 
           `Model Type` = model_label,
           `R²` = model_r2,
           `Interpretation` = r2_interpretation)
  
  return(r2_summary)
}

# GENERAR TABLAS LIMPIAS ================================================

# Verificar que tenemos ambos tipos de modelos
has_univariate <- any(all_results$model_type == "univariate")
has_multivariate <- any(all_results$model_type == "multivariate")

if (has_univariate) {
  # 1. Tabla univariante (R² por variable)
  univariate_table <- create_univariate_table(all_results)
  write.csv(univariate_table, "results/tables/univariate_models_table.csv", row.names = FALSE)
  
  cat("\n=== TABLA UNIVARIANTE CREADA ===\n")
  cat("Formato: Coeficiente (IC 95%) [R² individual] *\n")
  cat("* = Significativo\n")
} else {
  cat("No se encontraron modelos univariantes\n")
}

if (has_multivariate) {
  # 2. Tabla multivariante (R² del modelo como nota)
  multivariate_table <- create_multivariate_table(all_results)
  write.csv(multivariate_table, "results/tables/multivariate_models_table.csv", row.names = FALSE)
  
  # Guardar nota de R² en archivo separado
  r2_note <- attr(multivariate_table, "r2_note")
  if (!is.null(r2_note)) {
    writeLines(r2_note, "results/tables/multivariate_r2_note.txt")
    
    cat("\n=== TABLA MULTIVARIANTE CREADA ===\n")
    cat("Formato: Coeficiente (IC 95%) *\n")
    cat("* = Significativo\n")
    cat("R² del modelo:", r2_note, "\n")
  }
} else {
  cat("No se encontraron modelos multivariantes\n")
}

# 1. Tabla comparativa principal (tipo sjPlot) - MANTENEMOS LA ORIGINAL
main_comparison_table <- all_results %>%
  mutate(
    disp_label = case_when(
      dispersal_type == "average" ~ "Total",
      dispersal_type == "natal" ~ "Natal", 
      dispersal_type == "breeding" ~ "Breeding"
    ),
    dist_label = case_when(
      distance_type == "median" ~ "Median", 
      distance_type == "long" ~ "Long-distance"
    ),
    model_label = case_when(
      model_type == "univariate" ~ "Univariate",
      model_type == "multivariate" ~ "Multivariate"
    ),
    variable_clean = case_when(
      variable == "body_mass" ~ "Body mass",
      variable == "log_HWI" ~ "HWI",
      variable == "PC1" ~ "Life history",
      variable == "habita_for" ~ "Habitat openness",
      variable == "diet" ~ "Diet",
      variable == "distance_mig" ~ "Migration distance",
      variable == "Latitude" ~ "Latitude",
      TRUE ~ stringr::str_replace(variable, "_x_", " × ")
    ),
    col_header = paste(disp_label, dist_label, model_label, sep = " - "),
    coef_formatted = paste0(
      sprintf("%.3f", estimate),
      " (", sprintf("%.3f", lower_ci), ", ", sprintf("%.3f", upper_ci), ")",
      ifelse(significant, "*", "")
    )
  ) %>%
  select(variable_clean, col_header, coef_formatted) %>%
  pivot_wider(names_from = col_header, values_from = coef_formatted) %>%
  arrange(variable_clean)

write.csv(main_comparison_table, "results/tables/clean_univariate_vs_multivariate.csv", row.names = FALSE)

# 2. Tabla solo significativos
significant_only_table <- create_significant_only_table(all_results)
write.csv(significant_only_table, "results/tables/clean_significant_only.csv", row.names = FALSE)

# 3. Resumen de R² con interpretación
r2_summary_table <- create_r2_summary_table(all_results)
write.csv(r2_summary_table, "results/tables/r2_summary_interpretation.csv", row.names = FALSE)

# 4. Tabla comparativa original (mantenemos para compatibilidad)
comparison_table <- all_results %>%
  mutate(
    coef_text = paste0(
      round(estimate, 3), 
      " (", round(lower_ci, 3), ", ", round(upper_ci, 3), ")",
      ifelse(significant, "*", "")
    )
  ) %>%
  select(dispersal_type, distance_type, variable, coef_text, model_r2, model_type) %>%
  pivot_wider(
    names_from = c(dispersal_type, distance_type, model_type), 
    values_from = c(coef_text, model_r2),
    names_sep = "_"
  ) %>%
  arrange(variable)

write.csv(comparison_table, "results/tables/univariate_vs_multivariate_comparison.csv", row.names = FALSE)

# INTERPRETACIÓN DE R² ==================================================

cat("\n=== INTERPRETACIÓN DE R² ===\n")
cat("R² mide qué proporción de varianza explica el modelo:\n")
cat("• 0.00 - 0.10 = Efecto pequeño (0-10% de varianza explicada)\n")
cat("• 0.10 - 0.30 = Efecto moderado (10-30% de varianza explicada)\n") 
cat("• 0.30 - 0.50 = Efecto grande (30-50% de varianza explicada)\n")
cat("• 0.50+ = Efecto muy grande (>50% de varianza explicada)\n\n")

cat("Si R² multivariante > R² univariante:\n")
cat("→ El conjunto de variables explica más que las variables por separado\n")
cat("→ Hay sinergias/interacciones entre variables\n\n")

print("RESUMEN DE R² POR MODELO:")
print(r2_summary_table)

# Estadísticos resumen (incluyendo multivariantes)
summary_stats <- all_results %>%
  group_by(dispersal_type, distance_type, model_type) %>%
  summarise(
    n_variables = n(),
    n_significant = sum(significant),
    prop_significant = round(n_significant/n_variables, 2),
    mean_r2 = round(mean(model_r2, na.rm = TRUE), 3),
    mean_effect_size = round(mean(abs(estimate)), 3),
    .groups = "drop"
  )

write.csv(summary_stats, "results/tables/comprehensive_summary_stats.csv", row.names = FALSE)

# GENERAR GRÁFICOS ====================================================

cat("=== GENERANDO GRÁFICOS ===\n")

# Variables clave para plotear
key_vars <- c("body_mass", "log_HWI", "PC1")
var_labels <- c("Body Mass [log]", "HWI [log]", "Life History (slow-fast)")

# Función para crear plots para un tipo de dispersión
create_plots_for_type <- function(models, data, dispersal_type) {
  
  plots <- list()
  
  for (i in seq_along(key_vars)) {
    var <- key_vars[i]
    label <- var_labels[i]
    
    model_median <- models[[paste0(var, "_median")]]
    model_long <- models[[paste0(var, "_long")]]
    
    p <- plot_dispersal_predictions(
      model_median, model_long, data, var, label, dispersal_type
    )
    
    if (i > 1) p <- p + theme(legend.position = "none")
    plots[[var]] <- p
  }
  
  # Combinar plots
  combined <- plots[[1]] + plots[[2]] + plots[[3]]
  combined <- combined + plot_layout(ncol = 3)
  
  return(combined)
}

# Crear plots para cada tipo
plots_average <- create_plots_for_type(models_average, data_average, "average")
plots_natal <- create_plots_for_type(models_natal, data_natal, "natal") 
plots_breeding <- create_plots_for_type(models_breeding, data_breeding, "breeding")

# Guardar plots individuales
ggsave(plots_average, filename = "results/figures/univariate_plots_average.png", 
       width = 16, height = 5, dpi = 300)
ggsave(plots_natal, filename = "results/figures/univariate_plots_natal.png", 
       width = 16, height = 5, dpi = 300)
ggsave(plots_breeding, filename = "results/figures/univariate_plots_breeding.png", 
       width = 16, height = 5, dpi = 300)

# Plot combinado de todos los tipos
all_plots <- plots_average / plots_natal / plots_breeding
ggsave(all_plots, filename = "results/figures/univariate_comparison_grid.png", 
       width = 16, height = 12, dpi = 300)

# RESUMEN FINAL ========================================================

cat("\n=== RESUMEN DEL ANÁLISIS ===\n")
cat("Total coeficientes univariantes:", nrow(all_univariate_results), "\n")
cat("Total coeficientes multivariantes:", nrow(multivariate_results), "\n")
cat("Total coeficientes analizados:", nrow(all_results), "\n")
cat("Asociaciones significativas univariantes:", sum(all_univariate_results$significant), "\n")
if (nrow(multivariate_results) > 0) {
  cat("Asociaciones significativas multivariantes:", sum(multivariate_results$significant), "\n")
}

cat("\nPROPORCIÓN SIGNIFICATIVA POR TIPO DE MODELO:\n")
model_summary <- all_results %>%
  group_by(model_type) %>%
  summarise(
    total = n(),
    significant = sum(significant),
    prop_significant = round(significant/total, 2),
    .groups = "drop"
  )
print(model_summary)

cat("\nRESUMEN DETALLADO:\n")
print(summary_stats)

# Variables que cambian de significancia
if (nrow(multivariate_results) > 0) {
  cat("\n=== ANÁLISIS DE CAMBIOS EN SIGNIFICANCIA ===\n")
  
  # Variables que son significativas en univariante pero no en multivariante
  uni_sig <- all_univariate_results %>% 
    filter(significant) %>% 
    select(dispersal_type, distance_type, variable) %>%
    mutate(key = paste(dispersal_type, distance_type, variable, sep = "_"))
  
  multi_sig <- multivariate_results %>%
    filter(significant) %>%
    select(dispersal_type, distance_type, variable) %>%
    mutate(key = paste(dispersal_type, distance_type, variable, sep = "_"))
  
  # Perdieron significancia
  lost_significance <- setdiff(uni_sig$key, multi_sig$key)
  if (length(lost_significance) > 0) {
    cat("Variables que pierden significancia en multivariante:\n")
    for (var in lost_significance) {
      cat("- ", var, "\n")
    }
  }
  
  # Ganaron significancia  
  gained_significance <- setdiff(multi_sig$key, uni_sig$key)
  if (length(gained_significance) > 0) {
    cat("Variables que ganan significancia en multivariante:\n")
    for (var in gained_significance) {
      cat("- ", var, "\n")
    }
  }
}

cat("\nVariables más consistentes (>50% significant en univariantes):\n")
variable_consistency <- all_univariate_results %>%
  group_by(variable) %>%
  summarise(
    total_tests = n(),
    significant_tests = sum(significant),
    consistency = round(significant_tests/total_tests, 2),
    .groups = "drop"
  ) %>%
  filter(consistency > 0.5) %>%
  arrange(desc(consistency))

if (nrow(variable_consistency) > 0) {
  print(variable_consistency)
} else {
  cat("No hay variables consistentemente significativas.\n")
}

cat("\n=== ARCHIVOS GENERADOS ===\n")
cat("MODELOS:\n")
cat("- results/weibull/average/univariate_models.RData\n")
cat("- results/weibull/natal/univariate_models.RData\n") 
cat("- results/weibull/breeding/univariate_models.RData\n")
cat("\nTABLAS PRINCIPALES (RECOMENDADAS):\n")
cat("- results/tables/univariate_models_table.csv (R² individual por variable)\n")
cat("- results/tables/multivariate_models_table.csv (coeficientes + R² en nota)\n")
cat("- results/tables/multivariate_r2_note.txt (R² de cada modelo multivariante)\n")
cat("\nTABLAS COMPLEMENTARIAS:\n")
cat("- results/tables/clean_univariate_vs_multivariate.csv (comparación lado a lado)\n")
cat("- results/tables/clean_significant_only.csv (solo significativos)\n")
cat("- results/tables/r2_summary_interpretation.csv (R² con interpretación)\n")
cat("- results/tables/comprehensive_summary_stats.csv (estadísticos resumen)\n")
cat("\nGRÁFICOS:\n")
cat("- results/figures/univariate_plots_average.png\n")
cat("- results/figures/univariate_plots_natal.png\n")
cat("- results/figures/univariate_plots_breeding.png\n")
cat("- results/figures/univariate_comparison_grid.png\n")

cat("\n=== ANÁLISIS COMPLETADO ===\n")

# NOTAS PARA EL USUARIO ===============================================

cat("\n=== NOTAS IMPORTANTES ===\n")
cat("📊 TABLAS RECOMENDADAS:\n")
cat("   → univariate_models_table.csv (R² por variable individual)\n")
cat("   → multivariate_models_table.csv (R² del modelo completo)\n\n")

cat("🔍 FORMATO DE TABLAS:\n")
cat("   UNIVARIANTE: 0.450 (-0.120, 0.820) [R² = 0.234]*\n")
cat("                ├── Coeficiente (IC 95%) [R² individual] *\n")
cat("   MULTIVARIANTE: 0.123 (-0.045, 0.291)*\n")
cat("                  ├── Coeficientes (IC 95%) *\n")
cat("                  ├── R² del modelo completo en nota separada\n\n")

cat("⚠️  INTERPRETACIÓN:\n")
cat("   • R² univariante = Varianza explicada por esa variable sola\n")
cat("   • R² multivariante = Varianza explicada por todas las variables juntas\n")
cat("   • * = Significativo (IC no incluye 0)\n")
cat("   • Comparar ambas tablas para ver confounding/synergies\n")

