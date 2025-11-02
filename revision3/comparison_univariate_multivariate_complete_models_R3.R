# ==============================================================================
# COMPARACIÓN DE MODELOS UNIVARIANTES vs MULTIVARIANTES CON LOO
# ==============================================================================
#
# DESCRIPCIÓN:
# Script separado para comparar objetivamente modelos univariantes vs 
# multivariantes usando Leave-One-Out Cross-validation (LOO-CV)
#
# ESTRUCTURA DE MODELOS:
# Univariantes: models_average, models_natal, models_breeding (cada uno con múltiples variables)
# Multivariantes: all_short_models (tibble con modelos por dispersal type)
#
# MÉTRICAS:
# - ELPD_LOO: Expected Log Pointwise Predictive Density (mayor = mejor)
# - LOOIC: LOO Information Criterion (menor = mejor) 
# - SE: Standard Error del LOOIC
# - Δ: Diferencia relativa al mejor modelo
#
# INTERPRETACIÓN:
# - Δ ≤ 2: No diferencia sustancial
# - 2 < Δ ≤ 10: Diferencia sustancial  
# - Δ > 10: Evidencia fuerte a favor del mejor modelo
#
# ==============================================================================

# Cargar librerías
suppressPackageStartupMessages({
  library(loo)
  library(brms)
  library(dplyr)
  library(tidyr)
})

# CARGAR TODOS LOS MODELOS ===============================================

cat("=== CARGANDO MODELOS ===\n")

# Modelos univariantes
cat("Cargando modelos univariantes...\n")
load("results/weibull/average/univariate_models.RData")
models_average <- models
load("results/weibull/natal/univariate_models.RData")
models_natal <- models
load("results/weibull/breeding/univariate_models.RData")
models_breeding <- models

# Modelos multivariantes
cat("Cargando modelos multivariantes...\n")
load("results/combined/complete_models_combined.RData")

cat("✓ Todos los modelos cargados\n")

# FUNCIÓN PARA CALCULAR LOO DE MODELOS UNIVARIANTES =====================

calculate_univariate_loo <- function(models_list, dispersal_type) {
  
  cat("Calculando LOO para", dispersal_type, "univariantes...\n")
  
  results <- data.frame()
  
  for (model_name in names(models_list)) {
    model <- models_list[[model_name]]
    
    tryCatch({
      # Extraer información del nombre
      parts <- strsplit(model_name, "_")[[1]]
      distance_type <- parts[length(parts)]  # "median" o "long"
      variable <- paste(parts[-length(parts)], collapse = "_")
      
      # Calcular LOO
      loo_result <- loo(model, refresh = 0)
      
      # Extraer métricas
      result_row <- data.frame(
        model_type = "univariate",
        dispersal_type = dispersal_type,
        distance_type = distance_type,
        variable = variable,
        model_id = paste(dispersal_type, distance_type, variable, "uni", sep = "_"),
        elpd_loo = round(loo_result$estimates["elpd_loo", "Estimate"], 2),
        looic = round(loo_result$estimates["looic", "Estimate"], 2),
        se_looic = round(loo_result$estimates["looic", "SE"], 2),
        p_loo = round(loo_result$estimates["p_loo", "Estimate"], 2),
        n_params = length(fixef(model)[,1]) + 1,  # +1 para sigma
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, result_row)
      
    }, error = function(e) {
      cat("  ✗ Error con", model_name, ":", e$message, "\n")
    })
  }
  
  cat("  ✓ Completado:", nrow(results), "modelos\n")
  return(results)
}

# FUNCIÓN PARA CALCULAR LOO DE MODELOS MULTIVARIANTES ===================

calculate_multivariate_loo <- function() {
  
  cat("Calculando LOO para multivariantes...\n")
  
  results <- data.frame()
  
  for (i in 1:nrow(all_complete_models)) {
    model_info <- all_complete_models[i, ]
    model <- model_info$model[[1]]
    
    tryCatch({
      # Calcular LOO
      loo_result <- loo(model, refresh = 0)
      
      # Extraer métricas
      result_row <- data.frame(
        model_type = "multivariate",
        dispersal_type = model_info$age,
        distance_type = model_info$dispersal_mode,
        variable = "full_model",
        model_id = paste(model_info$age, model_info$dispersal_mode, "multi", sep = "_"),
        elpd_loo = round(loo_result$estimates["elpd_loo", "Estimate"], 2),
        looic = round(loo_result$estimates["looic", "Estimate"], 2),
        se_looic = round(loo_result$estimates["looic", "SE"], 2),
        p_loo = round(loo_result$estimates["p_loo", "Estimate"], 2),
        n_params = length(fixef(model)[,1]) + 1,
        stringsAsFactors = FALSE
      )
      
      results <- rbind(results, result_row)
      
    }, error = function(e) {
      cat("  ✗ Error con modelo", i, ":", e$message, "\n")
    })
  }
  
  cat("  ✓ Completado:", nrow(results), "modelos\n")
  return(results)
}

# CALCULAR LOO PARA TODOS LOS MODELOS ===================================

cat("\n=== CALCULANDO LOO ===\n")

# Univariantes
average_loo <- calculate_univariate_loo(models_average, "average")
natal_loo <- calculate_univariate_loo(models_natal, "natal") 
breeding_loo <- calculate_univariate_loo(models_breeding, "breeding")

# Combinar univariantes
all_univariate_loo <- rbind(average_loo, natal_loo, breeding_loo)

# Multivariantes
multivariate_loo <- calculate_multivariate_loo()

# Combinar todos
all_loo_results <- rbind(all_univariate_loo, multivariate_loo)

cat("\n✓ LOO calculado para", nrow(all_loo_results), "modelos totales\n")

# CALCULAR COMPARACIONES DELTA ===========================================

cat("=== CALCULANDO COMPARACIONES ===\n")

# Comparaciones por grupo (dispersal_type + distance_type)
loo_comparisons <- all_loo_results %>%
  group_by(dispersal_type, distance_type) %>%
  mutate(
    # Encontrar el mejor LOOIC (más bajo)
    best_looic = min(looic, na.rm = TRUE),
    best_elpd = max(elpd_loo, na.rm = TRUE),
    
    # Calcular deltas
    delta_looic = round(looic - best_looic, 2),
    delta_elpd = round(elpd_loo - best_elpd, 2),
    
    # Rankings
    looic_rank = rank(looic),
    elpd_rank = rank(-elpd_loo),  # Negativo porque mayor ELPD es mejor
    
    # Interpretación
    performance = case_when(
      delta_looic == 0 ~ "Best model",
      delta_looic <= 2 ~ "Competitive", 
      delta_looic <= 10 ~ "Worse",
      TRUE ~ "Much worse"
    ),
    
    # Diferencia en SE units
    se_units = round(delta_looic / se_looic, 2)
  ) %>%
  ungroup() %>%
  arrange(dispersal_type, distance_type, looic_rank)

# CREAR RESÚMENES ========================================================

cat("=== CREANDO RESÚMENES ===\n")

# 1. Resumen por tipo de modelo
model_type_summary <- loo_comparisons %>%
  group_by(model_type) %>%
  summarise(
    n_models = n(),
    mean_looic = round(mean(looic), 2),
    mean_elpd = round(mean(elpd_loo), 2),
    n_best = sum(delta_looic == 0),
    n_competitive = sum(delta_looic <= 2),
    prop_best = round(n_best / n(), 3),
    prop_competitive = round(n_competitive / n(), 3),
    .groups = "drop"
  )

# 2. Mejores modelos por grupo (sin importar si es uni o multi)
best_models <- loo_comparisons %>%
  group_by(dispersal_type, distance_type) %>%
  filter(looic_rank == 1) %>%
  select(dispersal_type, distance_type, model_type, variable, looic, delta_looic, performance) %>%
  ungroup()

# 3. Modelos competitivos por grupo (Δ ≤ 2)
competitive_models <- loo_comparisons %>%
  filter(delta_looic <= 2) %>%
  arrange(dispersal_type, distance_type, looic_rank) %>%
  select(dispersal_type, distance_type, model_type, variable, looic, delta_looic, looic_rank, performance)

# 4. Resumen de modelos competitivos por grupo
competitive_summary <- competitive_models %>%
  group_by(dispersal_type, distance_type) %>%
  summarise(
    n_competitive = n(),
    best_model_type = first(model_type),
    best_variable = first(variable),
    best_looic = first(looic),
    competitive_range = paste0("Δ = 0 - ", max(delta_looic)),
    model_types_included = paste(unique(model_type), collapse = ", "),
    .groups = "drop"
  )

# 5. Análisis de incertidumbre bayesiana
uncertainty_analysis <- loo_comparisons %>%
  group_by(dispersal_type, distance_type) %>%
  arrange(looic) %>%
  mutate(
    # Diferencia en unidades de SE
    se_difference = round(delta_looic / se_looic, 2),
    # Interpretación bayesiana
    bayesian_interpretation = case_when(
      delta_looic == 0 ~ "Best model",
      delta_looic <= 2 ~ "Competitive (within 2 LOOIC units)",
      delta_looic <= 4 & se_difference < 2 ~ "Possibly competitive (within 2 SE)",
      delta_looic <= 10 ~ "Substantially worse",
      TRUE ~ "Much worse"
    ),
    # Flag para modelos equivalentes
    equivalent_to_best = ifelse(delta_looic <= 2, "Yes", "No")
  ) %>%
  ungroup()

# 6. Comparación de efectividad por tipo de modelo  
model_effectiveness <- loo_comparisons %>%
  group_by(model_type) %>%
  summarise(
    n_models = n(),
    n_wins = sum(delta_looic == 0),
    n_competitive = sum(delta_looic <= 2),
    prop_wins = round(n_wins / n(), 3),
    prop_competitive = round(n_competitive / n(), 3),
    mean_looic = round(mean(looic), 2),
    median_delta = round(median(delta_looic), 2),
    .groups = "drop"
  )

# GUARDAR RESULTADOS =====================================================

cat("=== GUARDANDO RESULTADOS ===\n")

# Crear directorio si no existe
dir.create("results/tables/complete", recursive = TRUE, showWarnings = FALSE)

# Guardar todas las tablas
write.csv(loo_comparisons, "results/tables/loo_model_comparisons_complete.csv", row.names = FALSE)
write.csv(model_type_summary, "results/tables/loo_model_type_summary_complete.csv", row.names = FALSE)
write.csv(best_models, "results/tables/loo_best_models_complete.csv", row.names = FALSE)
write.csv(competitive_models, "results/tables/loo_competitive_models_complete.csv", row.names = FALSE)
write.csv(competitive_summary, "results/tables/loo_competitive_summary_complete.csv", row.names = FALSE)
write.csv(uncertainty_analysis, "results/tables/loo_uncertainty_analysis_complete.csv", row.names = FALSE)
write.csv(model_effectiveness, "results/tables/loo_model_effectiveness_complete.csv", row.names = FALSE)

# IMPRIMIR RESULTADOS ====================================================

cat("\n=== RESULTADOS DE COMPARACIÓN LOO ===\n")

cat("\n1. RESUMEN POR TIPO DE MODELO:\n")
print(model_type_summary)

cat("\n2. MEJORES MODELOS POR GRUPO (cualquier tipo):\n")
print(best_models)

cat("\n3. MODELOS COMPETITIVOS POR GRUPO (Δ ≤ 2):\n")
print(competitive_summary)

cat("\n4. DETALLES DE MODELOS COMPETITIVOS:\n")
print(head(competitive_models, 15))

cat("\n5. EFECTIVIDAD POR TIPO DE MODELO:\n")
print(model_effectiveness)

cat("\n6. INCERTIDUMBRE BAYESIANA - Casos con múltiples modelos equivalentes:\n")
uncertain_cases <- uncertainty_analysis %>%
  group_by(dispersal_type, distance_type) %>%
  filter(sum(equivalent_to_best == "Yes") > 1) %>%
  ungroup() %>%
  filter(equivalent_to_best == "Yes") %>%
  select(dispersal_type, distance_type, model_type, variable, looic, delta_looic, se_difference, bayesian_interpretation)

if (nrow(uncertain_cases) > 0) {
  print(uncertain_cases)
} else {
  cat("No hay casos con múltiples modelos equivalentes\n")
}

cat("\n5. DIFERENCIAS SUSTANCIALES (Δ LOOIC > 2):\n")
substantial_diffs <- loo_comparisons %>%
  filter(delta_looic > 2) %>%
  select(dispersal_type, distance_type, model_type, variable, delta_looic, performance) %>%
  arrange(desc(delta_looic))

if (nrow(substantial_diffs) > 0) {
  print(head(substantial_diffs, 10))
} else {
  cat("No hay diferencias sustanciales\n")
}

# INTERPRETACIÓN FINAL ===================================================

cat("\n=== INTERPRETACIÓN ===\n")
cat("📊 UMBRALES BAYESIANOS:\n")
cat("   • Δ ≤ 2: Modelos competitivos/equivalentes\n")
cat("   • 2 < Δ ≤ 4: Posiblemente diferentes (revisar SE)\n")
cat("   • 4 < Δ ≤ 10: Sustancialmente diferentes\n")
cat("   • Δ > 10: Evidencia fuerte\n\n")

cat("🔬 CRITERIOS DE INCERTIDUMBRE:\n")
cat("   • SE difference < 2: Diferencia no significativa\n")
cat("   • SE difference ≥ 2: Diferencia estadísticamente clara\n\n")

cat("🏆 RECOMENDACIÓN POR GRUPO:\n")
recommendations <- best_models %>%
  left_join(
    competitive_summary %>% select(dispersal_type, distance_type, n_competitive),
    by = c("dispersal_type", "distance_type")
  ) %>%
  mutate(
    recommendation = case_when(
      n_competitive == 1 ~ paste("Usar", model_type, ifelse(model_type == "univariate", paste("(", variable, ")"), "(full model)"), "- ganador claro"),
      n_competitive > 1 ~ paste("Múltiples modelos equivalentes (n =", n_competitive, ") - considerar parsimonia"),
      TRUE ~ "Revisar"
    )
  )

for (i in 1:nrow(recommendations)) {
  cat(sprintf("   • %s %s: %s\n", 
              recommendations$dispersal_type[i], 
              recommendations$distance_type[i],
              recommendations$recommendation[i]))
}

cat("\n📈 PARA TU PAPER:\n")
total_groups <- nrow(best_models)
multi_wins <- sum(best_models$model_type == "multivariate", na.rm = TRUE)
uni_wins <- sum(best_models$model_type == "univariate", na.rm = TRUE)
competitive_groups <- sum(competitive_summary$n_competitive > 1, na.rm = TRUE)

cat(sprintf("   • Multivariantes óptimos: %d/%d grupos (%.1f%%)\n", 
            multi_wins, total_groups, multi_wins/total_groups*100))
cat(sprintf("   • Univariantes óptimos: %d/%d grupos (%.1f%%)\n", 
            uni_wins, total_groups, uni_wins/total_groups*100))
cat(sprintf("   • Grupos con modelos equivalentes: %d/%d (%.1f%%)\n", 
            competitive_groups, total_groups, competitive_groups/total_groups*100))

cat("\n=== ARCHIVOS GENERADOS ===\n")
cat("📁 results/tables/:\n")
cat("   • loo_best_models.csv (MEJOR modelo por grupo)\n")
cat("   • loo_competitive_models.csv (modelos equivalentes Δ ≤ 2)\n")
cat("   • loo_competitive_summary.csv (resumen de equivalencias)\n")
cat("   • loo_uncertainty_analysis.csv (análisis de incertidumbre)\n")
cat("   • loo_model_effectiveness.csv (efectividad por tipo)\n")
cat("   • loo_model_comparisons.csv (comparaciones completas)\n")
cat("   • loo_model_type_summary.csv (resumen por tipo)\n")

cat("\n✅ ANÁLISIS COMPLETADO\n")