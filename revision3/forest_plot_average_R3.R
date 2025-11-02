# CARGAR LIBRERÍAS Y DATOS NECESARIOS ====================================
library(tidyverse)
library(tidybayes)
library(rcartocolor)
library(brms)
## Functionality ----------------------------------------------------------
`%nin%` <- Negate(`%in%`) # a function for negation of %in% function

myspread <- function(df, key, value) {
  # quote key
  keyq <- rlang::enquo(key)
  # break value vector into quotes
  valueq <- rlang::enquo(value)
  s <- rlang::quos(!!valueq)
  df %>% gather(variable, value, !!!s) %>%
    unite(temp, !!keyq, variable) %>%
    spread(temp, value)
}

# Load univariate models (only if not already loaded)
if (!exists("univariate_models") || !all(names(univariate_models) %in% paste0(
  c("body_mass", "PC1", "log_HWI", "diet", "habita_for", "distance_mig", "Latitude"),
  "_", rep(distance_types, each = 7)
))) {
  load(paste0("results/weibull/", disp_type, "/univariate_models.RData"))
  univariate_models <- models
}

# Load reduced models only once
if (!exists("reduced_models")) {
  load("results/combined/short_models_combined.RData")
  reduced_models <- all_short_models
}

# EXTRAER COEFICIENTES ===================================================

# Función para extraer coeficientes
extract_coefs <- function(model, model_type = "reduced") {
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
      model_type = model_type
    ) %>%
    select(variable, mean, q5, q95, model_type)
  
  return(coefs)
}

# Extraer coeficientes para median dispersal (modelos reducidos)
median_coefs <- reduced_models %>%
  filter(age == "average", dispersal_mode == "median") %>%
  rowwise() %>%
  mutate(coefs = list(extract_coefs(model, "reduced"))) %>%
  unnest(coefs) %>%
  select(age, function_id, variable, mean, q5, q95, model_type) %>%
  mutate(dispersal_type = "Median dispersal")

# Extraer coeficientes para long dispersal (solo log_HWI)
# Buscar el modelo correcto en la lista
hwi_model_name <- names(univariate_models)[str_detect(names(univariate_models), "log_HWI_long")]

if(length(hwi_model_name) > 0) {
  long_coefs <- extract_coefs(univariate_models[[hwi_model_name[1]]], "univariate") %>%
    filter(variable == "HWI") %>%
    mutate(
      age = "average",
      function_id = NA,
      dispersal_type = "Long-distance dispersal"
    ) %>%
    select(age, function_id, variable, mean, q5, q95, model_type, dispersal_type)
} else {
  # Si no encuentra el modelo, usar log_HWI_long directamente
  long_coefs <- extract_coefs(univariate_models$log_HWI_long, "univariate") %>%
    filter(variable == "HWI") %>%
    mutate(
      age = "average",
      function_id = NA,
      dispersal_type = "Long-distance dispersal"
    ) %>%
    select(age, function_id, variable, mean, q5, q95, model_type, dispersal_type)
}

# CARGAR SELECCIÓN DE VARIABLES PARA IMPORTANCIA =========================
median_var_sel <- read_csv("results/combined/median_variable_selection.csv")
long_var_sel <- read_csv("results/combined/long_variable_selection.csv")

# Función para procesar importancia
process_importance <- function(var_sel_data) {
  var_sel_data %>%
    filter(!is.na(size) & size > 0) %>%
    group_by(age) %>%
    mutate(
      max_size = max(size, na.rm = TRUE),
      importance = max_size - size + 1  # Sin escalar
    ) %>%
    ungroup() %>%
    mutate(
      variable = recode_factor(
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
        "distance_mig:Latitude" = "Distance migration : Latitude"
      )
    ) %>%
    select(age, variable, importance, size)
}

# Procesar importancia para median dispersal
median_importance <- process_importance(median_var_sel) %>%
  filter(age == "average") %>%
  mutate(dispersal_type = "Median dispersal")

# Para long dispersal, HWI tiene importancia máxima (es la única variable)
long_importance <- tibble(
  age = "average",
  variable = "HWI",
  importance = 1,  # Máxima importancia
  size = 1,
  dispersal_type = "Long-distance dispersal"
)

# COMBINAR DATOS ==========================================================
all_importance <- bind_rows(median_importance, long_importance)

forest_data <- bind_rows(median_coefs, long_coefs) %>%
  left_join(all_importance, by = c("age", "variable", "dispersal_type")) %>%
  mutate(
    # Definir orden de variables para el plot
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
    )),
    # Usar importancia original (sin escalar)
    importance = ifelse(is.na(importance), 1, importance)
  ) %>%
  # Normalizar importancia dentro de cada tipo de dispersión
  group_by(dispersal_type) %>%
  mutate(
    importance_normalized = scales::rescale(importance, to = c(1, 7))
  ) %>%
  ungroup()

# CREAR GRÁFICO ===========================================================

# Paleta de colores
my_colors <- c("Median dispersal" = "#1f77b4", "Long-distance dispersal" = "#ff7f0e")

# Forest plot
forest_plot <- ggplot(
  forest_data, 
  aes(y = variable, x = mean, xmin = q5, xmax = q95, color = dispersal_type)
) +
  geom_pointinterval(
    aes(fill = dispersal_type, size = importance_normalized),
    position = position_dodge(width = 0.5),
    alpha = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +
  scale_fill_manual(values = my_colors, name = "Dispersal Type") + 
  scale_color_manual(values = my_colors, name = "Dispersal Type") +
  scale_size_continuous(
    name = "Variable\nImportance\n(within type)",
    range = c(0.5, 4),
    breaks = c(1, 3, 5, 7),
    labels = c("Low", "Med", "High", "Max"),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    x = "Estimates",
    y = "Variable",
    #title = "Dispersal Effects: Reduced Model (Median) vs HWI (Long-distance)",
    #subtitle = "Average dispersal only • Point size = Variable importance (normalized within each dispersal type)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    title = element_text(size = 14),
    legend.position = "right",
    panel.grid.major.x = element_line(color = "grey90", size = 0.3)
  )

# MOSTRAR Y GUARDAR =======================================================
print(forest_plot)


# Guardar gráfico
ggsave("results/combined/forest_plot_reduced_vs_hwi.png", forest_plot, 
       width = 10, height = 6, dpi = 300)

# TABLA RESUMEN ==========================================================
summary_table <- forest_data %>%
  select(dispersal_type, variable, mean, q5, q95, importance, model_type) %>%
  mutate(
    estimate_ci = paste0(round(mean, 3), " [", round(q5, 3), ", ", round(q95, 3), "]")
  ) %>%
  select(dispersal_type, variable, estimate_ci, importance, model_type)

print("Resumen de efectos:")
print(summary_table)

# Guardar tabla
write_csv(summary_table, "results/combined/forest_summary_reduced_vs_hwi.csv")



# FIN DEL SCRIPT ==========================================================
