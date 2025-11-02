library(tidyverse)
library(tidybayes)
library(rcartocolor)
library(conflicted)
library(dplyr)
library(MASS) # Un ejemplo de un paquete que también tiene select

# Establecer la preferencia para select()
conflict_prefer("select", "dplyr")
conflicts_prefer(dplyr::filter)

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
    dplyr::select(age, function_id, dispersal_mode, variable, size, importance)
}

# Procesar datos de selección
median_var_processed <- process_var_selection(median_var_sel) %>%
  mutate(descriptor = "Median dispersal")

long_var_processed <- process_var_selection(long_var_sel) %>%
  mutate(descriptor = "Long-distance dispersal")

# EXTRAER COEFICIENTES DE MODELOS ========================================


# Función para extraer coeficientes de un modelo brms
extract_model_coefs <- function(model) {
  model %>%
    spread_draws(`b_.*`, regex = TRUE) %>%
    summarise_draws() %>%
    rename(param_name = 1) %>%  # Renombrar primera columna
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
  select(-model, -variable_selection) %>%
  # Seleccionar y renombrar
  select(age, function_id, dispersal_mode,
         variable = coefs_variable, # Renombrar aquí
         mean = coefs_mean,         # Renombrar aquí
         q5 = coefs_q5,             # Renombrar aquí
         q95 = coefs_q95            # Renombrar aquí
  ) %>%
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


# Paleta de colores para publicación (colorblind-friendly)
#pub_colors <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7")

# Paleta de colores estilo Nature (más sobria y profesional)
# Grises y azules con un acento naranja - muy usada en Nature
pub_colors <- c("#1B4F72", "#DC7633", "#196F3D", "#7D3C98")  # Azul oscuro, gris azulado, gris, naranja

# Alternativamente, paleta más clásica Nature:
# pub_colors <- c("#1B4F72", "#5D6D7E", "#85929E", "#DC7633")  # Azul marino, grises, naranja tierra

# O paleta más vibrante pero profesional:
# pub_colors <- c("#154360", "#5DADE2", "#48C9B0", "#F39C12")  # Azul oscuro, azul claro, turquesa, ámbar

# Tema base para publicación
theme_publication <- theme_classic() +
  theme(
    # Texto general más grande
    text = element_text(size = 14),
    # Ejes
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    # Strips para facets
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_blank(),
    # Leyenda
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(1.2, "cm"),
    # Panel
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

# Función para escalar importancia más dramáticamente
scale_importance_dramatic <- function(importance) {
  # Escalar de manera exponencial para acentuar diferencias
  scaled <- scales::rescale(importance, to = c(0.5, 3))
  # Aplicar transformación exponencial
  exp_scaled <- scaled^2
  # Escalar al rango final deseado
  scales::rescale(exp_scaled, to = c(1, 6))
}

# Aplicar escalado dramático
forest_data <- forest_data %>%
  mutate(
    size_dramatic = scale_importance_dramatic(importance)
  )

# 1. Gráfico separado por tipo de dispersión (median vs long)
plot_by_dispersal <- ggplot(
  forest_data, 
  aes(y = variable, x = mean, xmin = q5, xmax = q95)
) +
  geom_pointinterval(
    aes(colour = age, fill = age, size = size_dramatic),
    position = position_dodge(width = 0.6),
    alpha = 0.8,
    stroke = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 0.8, col = "black", alpha = 0.7) +
  scale_fill_manual(values = pub_colors, name = "Dispersal age") + 
  scale_color_manual(values = pub_colors, name = "Dispersal age") +
  scale_size_continuous(
    name = "Variable\nimportance",
    range = c(1.5, 6),
    breaks = c(2, 4, 6),
    labels = c("Low", "Medium", "High"),
    guide = guide_legend(
      override.aes = list(alpha = 1, stroke = 1),
      title.position = "top"
    )
  ) +
  facet_wrap(~descriptor, scales = "free_x") +
  labs(
    x = "Effect size",
    y = ""
  ) +
  theme_publication

# 2. Gráfico solo para average dispersal, comparando median vs long
plot_average_only <- ggplot(
  forest_data %>% filter(age == "average"), 
  aes(y = variable, x = mean, xmin = q5, xmax = q95)
) +
  geom_pointinterval(
    aes(colour = descriptor, fill = descriptor, size = size_dramatic),
    position = position_dodge(width = 0.6),
    alpha = 0.8,
    stroke = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 0.8, col = "black", alpha = 0.7) +
  scale_fill_manual(
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#E67E22"),
    name = "Dispersal type"
  ) +
  scale_color_manual(
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#E67E22"),
    name = "Dispersal type"
  ) +
  scale_size_continuous(
    name = "Variable\nimportance",
    range = c(1.5, 6),
    breaks = c(2, 4, 6),
    labels = c("Low", "Medium", "High"),
    guide = guide_legend(
      override.aes = list(alpha = 1, stroke = 1),
      title.position = "top"
    )
  ) +
  labs(
    x = "Effect size",
    y = ""
  ) +
  theme_publication

# 3. Gráfico por edad, comparando median vs long
plot_by_age <- ggplot(
  forest_data, 
  aes(y = variable, x = mean, xmin = q5, xmax = q95)
) +
  geom_pointinterval(
    aes(colour = descriptor, fill = descriptor, size = size_dramatic),
    position = position_dodge(width = 0.6),
    alpha = 0.8,
    stroke = 0.8
  ) +
  geom_vline(xintercept = 0, lty = 2, linewidth = 0.8, col = "black", alpha = 0.7) +
  scale_fill_manual(
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#E67E22"),
    name = "Dispersal type"
  ) +
  scale_color_manual(
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#E67E22"),
    name = "Dispersal type"
  ) +
  scale_size_continuous(
    name = "Variable\nimportance",
    range = c(1.5, 6),
    breaks = c(2, 4, 6),
    labels = c("Low", "Medium", "High"),
    guide = guide_legend(
      override.aes = list(alpha = 1, stroke = 1),
      title.position = "top"
    )
  ) +
  facet_wrap(~age, labeller = labeller(age = c(
    "average" = "Average",
    "breeding" = "Breeding", 
    "natal" = "Natal"
  ))) +
  labs(
    x = "Effect size",
    y = ""
  ) +
  theme_publication

# MOSTRAR Y GUARDAR GRÁFICOS =============================================

print(plot_by_dispersal)
print(plot_average_only) 
print(plot_by_age)

# Guardar gráficos en alta resolución (600 DPI)
ggsave("results/combined/short_forest_plot_by_dispersal.png", plot_by_dispersal, 
       width = 12, height = 8, dpi = 600, bg = "white")
ggsave("results/combined/short_forest_plot_average_only.png", plot_average_only, 
       width = 10, height = 7, dpi = 600, bg = "white")
ggsave("results/combined/short_forest_plot_by_age.png", plot_by_age, 
       width = 14, height = 8, dpi = 600, bg = "white")

# También guardar en formato vectorial para máxima calidad
ggsave("results/combined/short_forest_plot_by_dispersal.pdf", plot_by_dispersal, 
       width = 12, height = 8, bg = "white")
ggsave("results/combined/short_forest_plot_average_only.pdf", plot_average_only, 
       width = 10, height = 7, bg = "white")
ggsave("results/combined/short_forest_plot_by_age.pdf", plot_by_age, 
       width = 14, height = 8, bg = "white")

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
write_csv(importance_summary, "results/combined/short_variable_importance_summary.csv")