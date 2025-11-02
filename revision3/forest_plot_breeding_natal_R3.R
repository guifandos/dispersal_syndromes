# CARGAR LIBRERÍAS Y DATOS NECESARIOS ====================================
library(tidyverse)
library(tidybayes)
library(rcartocolor)
library(brms)

# Load reduced models only once
if (!exists("reduced_models")) {
  load("results/combined/short_models_combined.RData")
  reduced_models <- all_short_models
}


# Cargar selección de variables para importancia
median_var_sel <- read_csv("results/combined/median_variable_selection.csv")
long_var_sel <- read_csv("results/combined/long_variable_selection.csv")

# FUNCIÓN PARA EXTRAER COEFICIENTES ======================================
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

# FUNCIÓN PARA PROCESAR IMPORTANCIA ======================================
process_importance <- function(var_sel_data) {
  var_sel_data %>%
    filter(!is.na(size) & size > 0) %>%
    group_by(age) %>%
    mutate(
      max_size = max(size, na.rm = TRUE),
      importance = max_size - size + 1
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

# EXTRAER COEFICIENTES PARA BREEDING Y NATAL =============================

# Extraer coeficientes para breeding dispersal (median y long)
breeding_median_coefs <- reduced_models %>%
  filter(age == "breeding", dispersal_mode == "median") %>%
  rowwise() %>%
  mutate(coefs = list(extract_coefs(model, "reduced"))) %>%
  unnest(coefs) %>%
  select(age, function_id, variable, mean, q5, q95, model_type) %>%
  mutate(dispersal_type = "Median dispersal")

breeding_long_coefs <- reduced_models %>%
  filter(age == "breeding", dispersal_mode == "long") %>%
  rowwise() %>%
  mutate(coefs = list(extract_coefs(model, "reduced"))) %>%
  unnest(coefs) %>%
  select(age, function_id, variable, mean, q5, q95, model_type) %>%
  mutate(dispersal_type = "Long-distance dispersal")

# Extraer coeficientes para natal dispersal (median y long)
natal_median_coefs <- reduced_models %>%
  filter(age == "natal", dispersal_mode == "median") %>%
  rowwise() %>%
  mutate(coefs = list(extract_coefs(model, "reduced"))) %>%
  unnest(coefs) %>%
  select(age, function_id, variable, mean, q5, q95, model_type) %>%
  mutate(dispersal_type = "Median dispersal")

natal_long_coefs <- reduced_models %>%
  filter(age == "natal", dispersal_mode == "long") %>%
  rowwise() %>%
  mutate(coefs = list(extract_coefs(model, "reduced"))) %>%
  unnest(coefs) %>%
  select(age, function_id, variable, mean, q5, q95, model_type) %>%
  mutate(dispersal_type = "Long-distance dispersal")

# PROCESAR IMPORTANCIA ===================================================
median_importance <- process_importance(median_var_sel) %>%
  filter(age %in% c("breeding", "natal")) %>%
  mutate(dispersal_type = "Median dispersal")

long_importance <- process_importance(long_var_sel) %>%
  filter(age %in% c("breeding", "natal")) %>%
  mutate(dispersal_type = "Long-distance dispersal")

# COMBINAR TODOS LOS DATOS ===============================================
all_coefs <- bind_rows(
  breeding_median_coefs, breeding_long_coefs,
  natal_median_coefs, natal_long_coefs
)

all_importance <- bind_rows(median_importance, long_importance)

forest_data <- all_coefs %>%
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
    importance = ifelse(is.na(importance), 1, importance),
    # Crear etiquetas para facets
    age_label = case_when(
      age == "breeding" ~ "Breeding Dispersal",
      age == "natal" ~ "Natal Dispersal"
    )
  ) %>%
  # Normalizar importancia dentro de cada combinación de age y dispersal_type
  group_by(age, dispersal_type) %>%
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
  facet_wrap(~age_label, ncol = 2, scales = "free_y") +
  labs(
    x = "Estimates",
    y = "Variable",
    #title = "Dispersal Effects: Breeding and Natal Dispersal",
    #subtitle = "Reduced multivariate models • Point size = Variable importance (normalized within each type)"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    title = element_text(size = 13),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_line(color = "grey90", size = 0.3)
  )

# MOSTRAR Y GUARDAR =======================================================
print(forest_plot)

# Guardar gráfico
ggsave("results/combined/forest_plot_breeding_natal.png", forest_plot, 
       width = 14, height = 8, dpi = 300)

# TABLA RESUMEN ==========================================================
summary_table <- forest_data %>%
  select(age, dispersal_type, variable, mean, q5, q95, importance, importance_normalized, model_type) %>%
  mutate(
    estimate_ci = paste0(round(mean, 3), " [", round(q5, 3), ", ", round(q95, 3), "]")
  ) %>%
  select(age, dispersal_type, variable, estimate_ci, importance, importance_normalized, model_type) %>%
  arrange(age, dispersal_type, desc(importance_normalized))

print("Resumen de efectos para breeding y natal:")
print(summary_table)

# Guardar tabla
write_csv(summary_table, "results/combined/forest_summary_breeding_natal.csv")
