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

extract_multivariate_coefficients <- function(model_median, model_long, dispersal_type) {
  
  results <- data.frame()
  
  # Lista de modelos
  models_list <- list(
    "median" = model_median,
    "long" = model_long
  )
  
  for (distance_type in names(models_list)) {
    model <- models_list[[distance_type]]
    
    if (!is.null(model)) {
      # Get model summary
      model_summary <- summary(model)
      coef_data <- model_summary$fixed
      
      # Extraer coeficientes (saltando intercept)
      for (i in 2:nrow(coef_data)) {
        var_name <- rownames(coef_data)[i]
        
        # Limpiar nombre de variable (remover interacciones complejas)
        clean_var_name <- gsub(":", "_x_", var_name)
        
        result_row <- data.frame(
          dispersal_type = dispersal_type,
          distance_type = distance_type,
          variable = clean_var_name,
          estimate = coef_data[i, "Estimate"],
          se = coef_data[i, "Est.Error"],
          lower_ci = coef_data[i, "l-95% CI"],
          upper_ci = coef_data[i, "u-95% CI"],
          rhat = coef_data[i, "Rhat"],
          significant = !(coef_data[i, "l-95% CI"] <= 0 & coef_data[i, "u-95% CI"] >= 0),
          model_r2 = round(bayes_R2(model)[1, "Estimate"], 3),
          model_type = "multivariate"
        )
        
        results <- rbind(results, result_row)
      }
    }
  }
  
  return(results)
}

# FUNCIÓN PARA CARGAR MODELOS MULTIVARIANTES ===========================

load_multivariate_models <- function() {
  
  cat("=== CARGANDO MODELOS MULTIVARIANTES ===\n")
  
  multivariate_models <- list()
  
  # Rutas esperadas de modelos multivariantes
  paths <- list(
    average_median = "results/weibull/average/multivariate_median_model.RData",
    average_long = "results/weibull/average/multivariate_long_model.RData",
    natal_median = "results/weibull/natal/multivariate_median_model.RData", 
    natal_long = "results/weibull/natal/multivariate_long_model.RData",
    breeding_median = "results/weibull/breeding/multivariate_median_model.RData",
    breeding_long = "results/weibull/breeding/multivariate_long_model.RData"
  )
  
  for (name in names(paths)) {
    path <- paths[[name]]
    if (file.exists(path)) {
      cat("Cargando:", name, "\n")
      load(path)  # Asume que el modelo se guarda como 'model' o nombre específico
      
      # Intentar diferentes nombres de objeto comunes
      if (exists("model")) {
        multivariate_models[[name]] <- model
        rm(model)
      } else if (exists("final_model")) {
        multivariate_models[[name]] <- final_model  
        rm(final_model)
      } else if (exists("selected_model")) {
        multivariate_models[[name]] <- selected_model
        rm(selected_model)
      } else {
        # Si no encuentra el objeto, listar objetos disponibles
        objs <- ls()
        if (length(objs) > 0) {
          multivariate_models[[name]] <- get(objs[1])  # Tomar el primero
          cat("  Usando objeto:", objs[1], "\n")
        } else {
          cat("  ADVERTENCIA: No se encontró modelo en", path, "\n")
        }
      }
    } else {
      cat("ARCHIVO NO ENCONTRADO:", path, "\n")
    }
  }
  
  return(multivariate_models)
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
    ggtitle(paste("Univariate model:", x_axis_label, "-", stringr::str_to_title(dispersal_type))) +
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

# Natal
models_natal <- run_univariate_models(
  data = data_natal, 
  variables = variables, 
  A = A_natal, 
  dispersal_type = "natal",
  save_path = "results/weibull/natal/univariate_models.RData"
)

# Breeding
models_breeding <- run_univariate_models(
  data = data_breeding, 
  variables = variables, 
  A = A_breeding, 
  dispersal_type = "breeding",
  save_path = "results/weibull/breeding/univariate_models.RData"
)

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

# Coeficientes multivariantes
multivariate_models <- load_multivariate_models()

multivariate_results <- data.frame()

if (length(multivariate_models) > 0) {
  
  # Extraer para average
  if ("average_median" %in% names(multivariate_models) && "average_long" %in% names(multivariate_models)) {
    coef_avg_multi <- extract_multivariate_coefficients(
      multivariate_models$average_median, 
      multivariate_models$average_long, 
      "average"
    )
    multivariate_results <- rbind(multivariate_results, coef_avg_multi)
  }
  
  # Extraer para natal
  if ("natal_median" %in% names(multivariate_models) && "natal_long" %in% names(multivariate_models)) {
    coef_nat_multi <- extract_multivariate_coefficients(
      multivariate_models$natal_median, 
      multivariate_models$natal_long, 
      "natal"
    )
    multivariate_results <- rbind(multivariate_results, coef_nat_multi)
  }
  
  # Extraer para breeding
  if ("breeding_median" %in% names(multivariate_models) && "breeding_long" %in% names(multivariate_models)) {
    coef_breed_multi <- extract_multivariate_coefficients(
      multivariate_models$breeding_median, 
      multivariate_models$breeding_long, 
      "breeding"
    )
    multivariate_results <- rbind(multivariate_results, coef_breed_multi)
  }
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

# Tabla principal con todos los coeficientes (univariante vs multivariante)
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

# Tabla comparativa simplificada (solo coeficientes)
simple_comparison <- all_results %>%
  mutate(
    coef_text = paste0(
      round(estimate, 3), 
      ifelse(significant, "*", "")
    )
  ) %>%
  select(dispersal_type, distance_type, variable, coef_text, model_type) %>%
  unite("model_desc", dispersal_type, distance_type, model_type, sep = "_") %>%
  pivot_wider(
    names_from = model_desc,
    values_from = coef_text
  ) %>%
  arrange(variable)

write.csv(simple_comparison, "results/tables/simple_univariate_vs_multivariate.csv", row.names = FALSE)

# Tabla solo significativos (ambos tipos de modelo)
significant_table <- all_results %>%
  filter(significant == TRUE) %>%
  mutate(
    coef_text = paste0(round(estimate, 3), " [", round(lower_ci, 3), ", ", round(upper_ci, 3), "]")
  ) %>%
  select(dispersal_type, distance_type, variable, coef_text, model_r2, model_type) %>%
  arrange(model_type, dispersal_type, distance_type, variable)

write.csv(significant_table, "results/tables/significant_all_coefficients.csv", row.names = FALSE)

# Tabla de comparación directa univariante vs multivariante
if (nrow(multivariate_results) > 0) {
  
  # Variables en común entre ambos tipos de modelo
  common_vars <- intersect(
    unique(all_univariate_results$variable),
    unique(multivariate_results$variable)
  )
  
  if (length(common_vars) > 0) {
    
    direct_comparison <- all_results %>%
      filter(variable %in% common_vars) %>%
      mutate(
        model_key = paste(dispersal_type, distance_type, sep = "_"),
        coef_text = paste0(round(estimate, 3), ifelse(significant, "*", ""))
      ) %>%
      select(variable, model_key, coef_text, model_type) %>%
      pivot_wider(
        names_from = c(model_key, model_type),
        values_from = coef_text,
        names_sep = "_"
      ) %>%
      arrange(variable)
    
    write.csv(direct_comparison, "results/tables/direct_univariate_vs_multivariate.csv", row.names = FALSE)
    
    cat("Tabla de comparación directa creada con", length(common_vars), "variables\n")
  }
}

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
cat("\nTABLAS:\n")
cat("- results/tables/univariate_vs_multivariate_comparison.csv (tabla completa)\n")
cat("- results/tables/simple_univariate_vs_multivariate.csv (tabla simple)\n")
cat("- results/tables/direct_univariate_vs_multivariate.csv (comparación directa)\n")
cat("- results/tables/significant_all_coefficients.csv (solo significativos)\n")
cat("- results/tables/comprehensive_summary_stats.csv (estadísticos resumen)\n")
cat("\nGRÁFICOS:\n")
cat("- results/figures/univariate_plots_average.png\n")
cat("- results/figures/univariate_plots_natal.png\n")
cat("- results/figures/univariate_plots_breeding.png\n")
cat("- results/figures/univariate_comparison_grid.png\n")

cat("\n=== ANÁLISIS COMPLETADO ===\n")

# NOTAS PARA EL USUARIO ===============================================

cat("\n=== NOTAS IMPORTANTES ===\n")
cat("1. Si no tienes modelos multivariantes, solo se analizarán los univariantes\n")
cat("2. Los modelos multivariantes deben estar guardados con nombres específicos:\n")
cat("   - results/weibull/[type]/multivariate_median_model.RData\n")
cat("   - results/weibull/[type]/multivariate_long_model.RData\n")
cat("3. El objeto del modelo debe llamarse 'model', 'final_model' o 'selected_model'\n")
cat("4. Revisa la tabla 'direct_univariate_vs_multivariate.csv' para comparaciones directas\n")
cat("5. Variables que pierden significancia en multivariante sugieren confounding\n")