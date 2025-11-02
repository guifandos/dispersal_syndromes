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

# Cargar modelos combinados COMPLETOS
load("results/combined/complete_models_combined.RData")

# Cargar selección de variables
median_var_sel <- read_csv("results/combined/median_variable_selection.csv")
long_var_sel <- read_csv("results/combined/long_variable_selection.csv")

# PROCESAR SELECCIÓN DE VARIABLES ========================================

# Función para procesar selección de variables (igual que antes)
process_var_selection <- function(var_sel_data) {
  var_sel_data %>%
    filter(!is.na(size) & size > 0) %>%
    group_by(age) %>%
    mutate(
      # Usar ranking invertido como medida de importancia
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

# EXTRAER COEFICIENTES DE MODELOS COMPLETOS ==============================


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
model_coefs_complete <- all_complete_models %>%
  rowwise() %>%
  mutate(
    coefs = list(extract_model_coefs(model))
  ) %>%
  unnest(coefs, names_sep = "_") %>%
  select(-model) %>%
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

# Combinar todos los datos de selección
all_var_processed <- bind_rows(median_var_processed, long_var_processed)

# Unir coeficientes con importancia
forest_data_complete <- model_coefs_complete %>%
  left_join(
    all_var_processed,
    by = c("age", "variable", "descriptor")
  ) %>%
  # Mantener todas las variables de los modelos completos, incluso sin importancia
  mutate(
    # Para variables sin ranking, asignar importancia mínima
    importance = ifelse(is.na(importance), 1, importance)
  )

# CREAR GRÁFICOS =========================================================

# Paleta con mayor contraste para 3 categorías
pub_colors <- c("#1B4F72", "#DC7633", "#196F3D", "#7D3C98")

# Tema base para publicación
theme_publication <- theme_classic() +
  theme(
    text = element_text(size = 14),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15, face = "bold"),
    axis.ticks = element_line(linewidth = 0.8, color = "black"),
    axis.ticks.length = unit(0.3, "cm"),
    axis.line = element_line(linewidth = 0.8, color = "black"),
    strip.text = element_text(size = 14, face = "bold", color = "black"),
    strip.background = element_blank(),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(1.2, "cm"),
    panel.background = element_blank(),
    panel.grid = element_blank()
  )

# Función para escalar importancia dramáticamente
scale_importance_dramatic <- function(importance) {
  scaled <- scales::rescale(importance, to = c(0.5, 3))
  exp_scaled <- scaled^2
  scales::rescale(exp_scaled, to = c(1, 6))
}

# Aplicar escalado dramático
forest_data_complete <- forest_data_complete %>%
  mutate(
    size_dramatic = scale_importance_dramatic(importance)
  )

# 1. Gráfico separado por tipo de dispersión (median vs long) - MODELOS COMPLETOS
plot_complete_by_dispersal <- ggplot(
  forest_data_complete, 
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

# 2. Gráfico solo para average dispersal - MODELOS COMPLETOS
plot_complete_average_only <- ggplot(
  forest_data_complete %>% filter(age == "average"), 
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
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#DC7633"),
    name = "Dispersal type"
  ) +
  scale_color_manual(
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#DC7633"),
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

# 3. Gráfico por edad - MODELOS COMPLETOS
plot_complete_by_age <- ggplot(
  forest_data_complete, 
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
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#DC7633"),
    name = "Dispersal type"
  ) +
  scale_color_manual(
    values = c("Median dispersal" = "#1B4F72", "Long-distance dispersal" = "#DC7633"),
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

# MOSTRAR Y GUARDAR GRÁFICOS COMPLETOS ===================================

print(plot_complete_by_dispersal)
print(plot_complete_average_only) 
print(plot_complete_by_age)

# Guardar gráficos en alta resolución (600 DPI)
ggsave("results/combined/forest_plot_COMPLETE_by_dispersal.png", plot_complete_by_dispersal, 
       width = 12, height = 10, dpi = 600, bg = "white")
ggsave("results/combined/forest_plot_COMPLETE_average_only.png", plot_complete_average_only, 
       width = 10, height = 8, dpi = 600, bg = "white")
ggsave("results/combined/forest_plot_COMPLETE_by_age.png", plot_complete_by_age, 
       width = 14, height = 10, dpi = 600, bg = "white")

# También guardar en formato vectorial
ggsave("results/combined/forest_plot_COMPLETE_by_dispersal.pdf", plot_complete_by_dispersal, 
       width = 12, height = 10, bg = "white")
ggsave("results/combined/forest_plot_COMPLETE_average_only.pdf", plot_complete_average_only, 
       width = 10, height = 8, bg = "white")
ggsave("results/combined/forest_plot_COMPLETE_by_age.pdf", plot_complete_by_age, 
       width = 14, height = 10, bg = "white")

# TABLA RESUMEN MODELOS COMPLETOS ========================================

# Crear tabla resumen de importancia para modelos completos
importance_summary_complete <- forest_data_complete %>%
  group_by(variable, age, descriptor) %>%
  summarise(
    mean_effect = round(mean, 3),
    importance = unique(importance),
    ranking = unique(ifelse(is.na(size), "Not selected", size)),
    .groups = "drop"
  ) %>%
  arrange(age, descriptor, desc(importance))

print("Resumen de importancia de variables - MODELOS COMPLETOS:")
print(importance_summary_complete)

# Guardar tabla
write_csv(importance_summary_complete, "results/combined/variable_importance_summary_COMPLETE.csv")

cat("\nArchivos de modelos COMPLETOS creados:\n")
cat("- forest_plot_COMPLETE_by_dispersal.png/pdf\n")
cat("- forest_plot_COMPLETE_average_only.png/pdf\n") 
cat("- forest_plot_COMPLETE_by_age.png/pdf\n")
cat("- variable_importance_summary_COMPLETE.csv\n")

