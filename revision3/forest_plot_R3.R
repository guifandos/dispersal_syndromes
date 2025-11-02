library(tidyverse)
library(tidybayes)
library(rcartocolor)

# CARGAR MODELOS Y DATOS ==================================================

# Cargar modelos combinados
load("results/combined/short_models_combined.RData")

# Cargar selección de variables
median_var_sel <- read_csv("results/combined/median_variable_selection.csv")
long_var_sel <- read_csv("results/combined/long_variable_selection.csv")

# PROCESAR SELECCIÓN DE VARIABLES ========================================

# Función para procesar selección de variables
process_var_selection <- function(var_sel_data) {
  var_sel_data %>%
    filter(!is.na(size) & size > 0) %>%
    group_by(age) %>%
    mutate(
      # Usar ranking invertido como medida de importancia
      # size = 1 (más importante) -> importance alta
      # size = 7 (menos importante) -> importance baja
      max_size = max(size, na.rm = TRUE),
      importance = max_size - size + 1
    ) %>%
    ungroup() %>%
    mutate(
      # Limpiar nombres de variables
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
    select(age, function_id, dispersal_mode, variable, size, importance)
}

# Procesar datos de selección
median_var_processed <- process_var_selection(median_var_sel) %>%
  mutate(descriptor = "Median dispersal")

long_var_processed <- process_var_selection(long_var_sel) %>%
  mutate(descriptor = "Long-distance dispersal")

# EXTRAER COEFICIENTES DE MODELOS ========================================

# Función para extraer coeficientes de un modelo brms
extract_model_coefs <- function(model, age_name, descriptor_name) {
  model %>%
    spread_draws(`b_.*`, regex = TRUE) %>%
    summarise_draws() %>%
    filter(str_detect(variable, "^b_") & variable != "b_Intercept") %>%
    mutate(
      age = age_name,
      descriptor = descriptor_name,
      # Limpiar nombres de variables (quitar "b_")
      variable_raw = str_remove(variable, "^b_"),
      variable = recode_factor(
        variable_raw,
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
    select(variable, mean, q5, q95, age, descriptor)
}

# Extraer coeficientes para todos los modelos
model_coefs <- all_short_models %>%
  rowwise() %>%
  mutate(
    coefs = list(extract_model_coefs(model, age, dispersal_mode))
  ) %>%
  unnest(coefs) %>%
  select(-model, -variable_selection) %>%
  mutate(
    descriptor = case_when(
      dispersal_mode == "median" ~ "Median dispersal",
      dispersal_mode == "long" ~ "Long-distance dispersal"
    )
  )

# COMBINAR COEFICIENTES CON IMPORTANCIA ==================================

# Combinar todos los datos
all_var_processed <- bind_rows(median_var_processed, long_var_processed)

# Unir coeficientes con importancia
forest_data <- model_coefs %>%
  left_join(
    all_var_processed,
    by = c("age", "variable", "descriptor")
  ) %>%
  filter(!is.na(importance)) %>%
  mutate(
    # Escalar importancia para el tamaño de puntos
    size_scaled = scales::rescale(importance, to = c(1, 5))
  )

# CREAR GRÁFICOS =========================================================

# Paleta de colores
my_pal <- rcartocolor::carto_pal(n = 8, name = "Bold")[c(1, 3, 7, 2)]

# 1. Gráfico separado por tipo de dispersión (median vs long)
plot_by_dispersal <- ggplot(
  forest_data, 
  aes(y = variable, x = mean, xmin = q5, xmax = q95)
) +
  geom_pointinterval(
    aes(colour = age, fill = age, size = size_scaled),
    position = position_dodge(width = 0.5),
    alpha = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +
  scale_fill_manual(values = my_pal, name = "Dispersal Age") + 
  scale_color_manual(values = my_pal, name = "Dispersal Age") +
  scale_size_continuous(
    name = "Variable\nImportance",
    range = c(1, 5),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  facet_wrap(~descriptor) +
  labs(
    x = "Estimates",
    y = "Variable",
    title = "Dispersal Syndromes: Variable Effects by Dispersal Type"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "right"
  )

# 2. Gráfico solo para average dispersal, comparando median vs long
plot_average_only <- ggplot(
  forest_data %>% filter(age == "average"), 
  aes(y = variable, x = mean, xmin = q5, xmax = q95)
) +
  geom_pointinterval(
    aes(colour = descriptor, fill = descriptor, size = size_scaled),
    position = position_dodge(width = 0.5),
    alpha = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +
  scale_fill_manual(
    values = c("Median dispersal" = "#1f77b4", "Long-distance dispersal" = "#ff7f0e"),
    name = "Dispersal Type"
  ) +
  scale_color_manual(
    values = c("Median dispersal" = "#1f77b4", "Long-distance dispersal" = "#ff7f0e"),
    name = "Dispersal Type"
  ) +
  scale_size_continuous(
    name = "Variable\nImportance",
    range = c(1, 5),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  labs(
    x = "Estimates",
    y = "Variable",
    title = "Average Dispersal: Median vs Long-distance Effects"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    title = element_text(size = 14),
    legend.position = "right"
  )

# 3. Gráfico por edad, comparando median vs long
plot_by_age <- ggplot(
  forest_data, 
  aes(y = variable, x = mean, xmin = q5, xmax = q95)
) +
  geom_pointinterval(
    aes(colour = descriptor, fill = descriptor, size = size_scaled),
    position = position_dodge(width = 0.5),
    alpha = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, size = 0.5, col = "grey") +
  scale_fill_manual(
    values = c("Median dispersal" = "#1f77b4", "Long-distance dispersal" = "#ff7f0e"),
    name = "Dispersal Type"
  ) +
  scale_color_manual(
    values = c("Median dispersal" = "#1f77b4", "Long-distance dispersal" = "#ff7f0e"),
    name = "Dispersal Type"
  ) +
  scale_size_continuous(
    name = "Variable\nImportance",
    range = c(1, 5),
    guide = guide_legend(override.aes = list(alpha = 1))
  ) +
  facet_wrap(~age, labeller = labeller(age = c(
    "average" = "Average",
    "breeding" = "Breeding", 
    "natal" = "Natal"
  ))) +
  labs(
    x = "Estimates",
    y = "Variable",
    title = "Dispersal Syndromes by Age Class"
  ) +
  theme_classic() +
  theme(
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    strip.text = element_text(size = 11, face = "bold"),
    legend.position = "right"
  )

# MOSTRAR Y GUARDAR GRÁFICOS =============================================

print(plot_by_dispersal)
print(plot_average_only) 
print(plot_by_age)

# Guardar gráficos
ggsave("results/combined/forest_plot_by_dispersal.png", plot_by_dispersal, 
       width = 12, height = 8, dpi = 300)
ggsave("results/combined/forest_plot_average_only.png", plot_average_only, 
       width = 10, height = 6, dpi = 300)
ggsave("results/combined/forest_plot_by_age.png", plot_by_age, 
       width = 14, height = 8, dpi = 300)

# TABLA RESUMEN ==========================================================

# Crear tabla resumen de importancia
importance_summary <- forest_data %>%
  group_by(variable, age, descriptor) %>%
  summarise(
    mean_effect = round(mean, 3),
    importance = unique(importance),
    ranking = unique(size),
    .groups = "drop"
  ) %>%
  arrange(age, descriptor, desc(importance))

print("Resumen de importancia de variables:")
print(importance_summary)

# Guardar tabla
write_csv(importance_summary, "results/combined/variable_importance_summary.csv")