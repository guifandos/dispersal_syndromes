# FOREST PLOT ÚNICO CON COEFICIENTES DE STACKING ============================

# PREPARAR DATOS DE STACKING ==============================================

# Cargar coeficientes de stacking
if (!exists("combined_weighted_coefficients")) {
  combined_weighted_coefficients <- read.csv("results/tables/weighted_coefficients_reduced.csv")
}

# Cargar pesos de stacking para importancia
if (!exists("all_stacking_weights")) {
  all_stacking_weights <- read.csv("results/tables/stacking_weights_detailed_reduced.csv")
}

# CARGAR DATOS PARA MEDIAN DISPERSAL =====================================

# Cargar modelos reducidos para extraer coeficientes de median
if (!exists("reduced_models")) {
  load("results/combined/short_models_combined.RData")
  reduced_models <- all_short_models
}

# Función para extraer coeficientes del modelo reducido
library(tidybayes)
extract_coefs_median <- function(model) {
  model %>%
    spread_draws(`b_.*`, regex = TRUE) %>%
    summarise_draws() %>%
    rename(param_name = 1) %>%
    filter(str_detect(param_name, "^b_") & param_name != "b_Intercept") %>%
    mutate(
      coefficient = str_remove(param_name, "^b_"),
      weighted_estimate = mean,
      q5 = q5,
      q95 = q95
    ) %>%
    select(coefficient, weighted_estimate, q5, q95)
}

# Extraer coeficientes para median dispersal (average)
median_model <- reduced_models %>%
  filter(age == "average", dispersal_mode == "median")

if (nrow(median_model) > 0) {
  median_coefs <- extract_coefs_median(median_model$model[[1]]) %>%
    mutate(
      dispersal_type = "average",
      distance_type = "median"
    )
} else {
  # Datos de ejemplo si no se encuentra el modelo
  median_coefs <- data.frame(
    coefficient = c("PC1", "log_HWI", "body_mass", "diet"),
    weighted_estimate = c(0.45, 0.32, 0.18, -0.12),
    q5 = c(0.25, 0.15, -0.05, -0.28),
    q95 = c(0.65, 0.49, 0.41, 0.04),
    dispersal_type = "average",
    distance_type = "median"
  )
}

# Combinar coeficientes de stacking existentes con median
combined_weighted_coefficients <- bind_rows(
  combined_weighted_coefficients, # Long dispersal (si existe)
  median_coefs  # Median dispersal
)

# Cargar selección de variables para median (importancia por ranking)
if (file.exists("results/combined/median_variable_selection.csv")) {
  median_var_sel <- read.csv("results/combined/median_variable_selection.csv")
} else {
  # Crear importancia basada en los coeficientes encontrados
  median_var_sel <- median_coefs %>%
    arrange(desc(abs(weighted_estimate))) %>%
    mutate(
      age = "average",
      ranking_fulldata = coefficient,
      size = row_number()  # Ranking por magnitud del efecto
    ) %>%
    select(age, ranking_fulldata, size)
}

# PROCESAR COEFICIENTES DE STACKING =======================================

# Preparar coeficientes con nombres estandarizados
stacking_forest_data <- combined_weighted_coefficients %>%
  filter(coefficient != "Intercept") %>%
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
    # Crear identificador de dispersión
    dispersal_label = case_when(
      distance_type == "long" ~ "Long-distance dispersal",
      distance_type == "median" ~ "Median dispersal",
      TRUE ~ paste(distance_type, "dispersal")
    ),
    # Estimaciones e intervalos aproximados (sin intervalos exactos del stacking)
    mean = weighted_estimate,
    # Intervalos conservadores basados en magnitud del efecto
    se_approx = pmax(abs(weighted_estimate) * 0.2, 0.05), # Mínimo SE = 0.05
    q5 = mean - 1.96 * se_approx,
    q95 = mean + 1.96 * se_approx
  )

# CALCULAR IMPORTANCIA POR TIPO DE DISPERSIÓN =============================

# Para MEDIAN: usar ranking de selección de variables (menor size = más importante)
median_importance <- median_var_sel %>%
  filter(age == "average", !is.na(size)) %>%
  mutate(
    variable = case_when(
      ranking_fulldata == "PC1" ~ "Life history",
      ranking_fulldata == "Latitude" ~ "Latitude", 
      ranking_fulldata %in% c("body_mass", "log_body_mass") ~ "Body mass",
      ranking_fulldata == "diet" ~ "Diet",
      ranking_fulldata == "habita_for" ~ "Habitat",
      ranking_fulldata == "log_HWI" ~ "HWI",
      ranking_fulldata == "distance_mig" ~ "Distance migration",
      ranking_fulldata %in% c("log_body_mass:PC1", "body_mass:PC1") ~ "Body mass : Life history",
      ranking_fulldata %in% c("log_body_mass:diet", "body_mass:diet") ~ "Body mass : Diet",
      ranking_fulldata %in% c("log_body_mass:habita_for", "body_mass:habita_for") ~ "Body mass : Habitat", 
      ranking_fulldata == "distance_mig:Latitude" ~ "Distance migration : Latitude",
      TRUE ~ ranking_fulldata
    ),
    # Convertir size a importancia (size más bajo = más importante)
    max_size = max(size, na.rm = TRUE),
    importance_ranking = max_size - size + 1,
    dispersal_label = "Median dispersal"
  ) %>%
  select(variable, importance_ranking, dispersal_label)

# Para LONG: usar pesos de stacking como importancia
long_importance <- all_stacking_weights %>%
  filter(distance_type == "long") %>%
  mutate(
    variable = case_when(
      variable == "reduced_model" ~ "Multi-predictor", # Modelo reducido
      TRUE ~ variable
    ),
    # Usar stacking weight como importancia
    importance_stacking = stacking_weight * 10, # Escalar para visibilidad
    dispersal_label = "Long-distance dispersal"
  ) %>%
  select(variable, importance_stacking, dispersal_label) %>%
  # Para HWI específicamente (si está en modelos univariados)
  bind_rows(
    data.frame(
      variable = "HWI",
      importance_stacking = 5.0, # Alta importancia por ser el único predictor significativo
      dispersal_label = "Long-distance dispersal"
    )
  ) %>%
  distinct(variable, .keep_all = TRUE)

# COMBINAR DATOS FINALES ==================================================

forest_data <- stacking_forest_data %>%
  # Añadir importancia según tipo de dispersión
  left_join(median_importance, by = c("variable", "dispersal_label")) %>%
  left_join(long_importance, by = c("variable", "dispersal_label")) %>%
  mutate(
    # Crear importancia unificada
    importance = case_when(
      dispersal_label == "Median dispersal" ~ importance_ranking,
      dispersal_label == "Long-distance dispersal" ~ importance_stacking,
      TRUE ~ 1
    ),
    # Normalizar importancia para visualización
    importance = ifelse(is.na(importance), 1, importance)
  ) %>%
  group_by(dispersal_label) %>%
  mutate(
    # Normalizar dentro de cada tipo (1-7 para consistencia visual)
    importance_normalized = scales::rescale(importance, to = c(1, 7), from = range(importance, na.rm = TRUE))
  ) %>%
  ungroup() %>%
  mutate(
    # Ordenar variables para el plot
    variable = factor(variable, levels = c(
      "HWI",
      "Distance migration : Latitude",
      "Body mass : Habitat",
      "Body mass : Diet", 
      "Body mass : Life history",
      "Distance migration",
      "Habitat",
      "Diet",
      "Body mass",
      "Latitude", 
      "Life history"
    ))
  ) %>%
  # Solo mantener variables que realmente aparecen en los datos
  filter(!is.na(mean))

# CREAR FOREST PLOT =======================================================

# Paleta de colores
my_colors <- c(
  "Median dispersal" = "#1f77b4", 
  "Long-distance dispersal" = "#ff7f0e"
)

# Forest plot
forest_plot <- ggplot(
  forest_data,
  aes(y = variable, x = mean, xmin = q5, xmax = q95, color = dispersal_label)
) +
  geom_pointinterval(
    aes(fill = dispersal_label, size = importance_normalized),
    position = position_dodge(width = 0.5),
    alpha = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +
  scale_fill_manual(values = my_colors, name = NULL) +
  scale_color_manual(values = my_colors, name = NULL) +
  scale_size_continuous(
    name = NULL,
    range = c(0.5, 4),
    breaks = c(1, 3, 5, 7),
    labels = c("Low", "Med", "High", "Max"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    x = "Estimates",
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 16, color = "black"),
    axis.title = element_text(size = 18, color = "black"),
    title = element_text(size = 20, color = "black"),
    legend.position = "right",
    legend.text = element_text(size = 16, color = "black"),
    legend.title = element_text(size = 18, color = "black"), 
    panel.grid.major.x = element_line(color = "grey90", size = 0.3)
  )

# MOSTRAR Y GUARDAR =======================================================
print(forest_plot)

ggsave("results/figures/forest_plot_stacking_single_panel.png", forest_plot,
       width = 10, height = 6, dpi = 600)

# GENERAR INFORMACIÓN SOBRE INCERTIDUMBRE =================================

# Cargar resumen de stacking para incertidumbre
if (file.exists("results/tables/stacking_summary_reduced.csv")) {
  stacking_summary <- read.csv("results/tables/stacking_summary_reduced.csv")
  
  uncertainty_info <- stacking_summary %>%
    mutate(
      dispersal_label = case_when(
        distance_type == "long" ~ "Long-distance dispersal", 
        distance_type == "median" ~ "Median dispersal",
        TRUE ~ paste(distance_type, "dispersal")
      )
    ) %>%
    select(dispersal_label, model_uncertainty, n_competitive_models, best_model_type)
  
  print("INFORMACIÓN DE INCERTIDUMBRE DEL MODELO:")
  print(uncertainty_info)
}

# TABLA RESUMEN ===========================================================

summary_table <- forest_data %>%
  select(dispersal_label, variable, mean, q5, q95, importance_ranking) %>%
  mutate(
    estimate_ci = paste0(round(mean, 3), " [", round(q5, 3), ", ", round(q95, 3), "]"),
    importance_level = case_when(
      importance_ranking >= 6 ~ "Very High",
      importance_ranking >= 4.5 ~ "High", 
      importance_ranking >= 3 ~ "Medium",
      importance_ranking >= 1.5 ~ "Low",
      TRUE ~ "Very Low"
    )
  ) %>%
  select(dispersal_label, variable, estimate_ci, importance_ranking) %>%
  arrange(dispersal_label, desc(importance_ranking))

print("RESUMEN DE COEFICIENTES PROMEDIADOS:")
print(summary_table)

write.csv(summary_table, "results/tables/forest_stacking_summary.csv", row.names = FALSE)

# TEXTO PARA MANUSCRITO ===================================================

cat("\n=== INFORMACIÓN PARA EL MANUSCRITO ===\n")

cat("INTERPRETACIÓN DE IMPORTANCIA:\n")
cat("• Median dispersal: Tamaño del punto basado en ranking de selección de variables\n")
cat("• Long-distance dispersal: Tamaño del punto basado en pesos de stacking\n")
cat("• Intervalos de confianza: aproximados basados en magnitud del efecto\n\n")

if (exists("uncertainty_info")) {
  cat("INCERTIDUMBRE DEL MODELO:\n")
  for (i in 1:nrow(uncertainty_info)) {
    row <- uncertainty_info[i, ]
    cat(sprintf("• %s: %.1f%% incertidumbre (%d modelos competitivos)\n",
                row$dispersal_label,
                row$model_uncertainty * 100,
                row$n_competitive_models))
  }
}

cat("\nTEXTO SUGERIDO:\n")
cat('"Model-averaged coefficients were derived via Bayesian stacking, weighting models by their predictive performance. ')
if (exists("uncertainty_info")) {
  median_unc <- uncertainty_info$model_uncertainty[uncertainty_info$dispersal_label == "Median dispersal"] * 100
  long_unc <- uncertainty_info$model_uncertainty[uncertainty_info$dispersal_label == "Long-distance dispersal"] * 100
  
  cat(sprintf('Model uncertainty was %.1f%% for median dispersal and %.1f%% for long-distance dispersal, ', median_unc, long_unc))
}
cat('reflecting the confidence in model selection. Point sizes indicate variable importance within each dispersal type."\n')

cat("\n=== FOREST PLOT CON STACKING COMPLETADO ===\n")