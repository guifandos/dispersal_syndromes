# CARGAR LIBRERÍAS Y DATOS NECESARIOS ====================================
library(tidyverse)
library(tidybayes)
library(brms)
library(performance)

# Cargar modelos reducidos
load("results/combined/short_models_combined.RData")
reduced_models <- all_short_models

# Cargar modelos univariados para average long (HWI)
load(paste0("results/weibull/", disp_type, "/univariate_models.RData"))
univariate_models <- models

# Cargar selección de variables para importancia
median_var_sel <- read_csv("results/combined/median_variable_selection.csv")
long_var_sel <- read_csv("results/combined/long_variable_selection.csv")

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
    
    # Verificar qué columnas existen
    r2_names <- names(r2_result)
    cat("  - Columnas R² disponibles:", paste(r2_names, collapse = ", "), "\n")
    
    # Extraer valores con nombres correctos
    if("R2_Bayes_marginal" %in% r2_names) {
      r2_marginal <- round(r2_result$R2_Bayes_marginal, 3)
    }
    
    if("R2_Bayes_conditional" %in% r2_names) {
      r2_conditional <- round(r2_result$R2_Bayes_conditional, 3)
    } else if("R2_Bayes" %in% r2_names) {
      r2_conditional <- round(r2_result$R2_Bayes, 3)
      cat("  - Usando R2_Bayes como condicional\n")
    }
    
    cat("  - R² marginal:", r2_marginal, "| R² condicional:", r2_conditional, "\n")
    
  }, error = function(e) {
    cat("Error en R² para", model_name, ":", e$message, "\n")
    
    # Intentar método alternativo con brms::bayes_R2()
    tryCatch({
      cat("  - Intentando método alternativo brms::bayes_R2()...\n")
      r2_alt <- bayes_R2(model)
      r2_marginal <<- round(mean(r2_alt), 3)
      r2_conditional <<- "Alt method"
      cat("  - Método alternativo exitoso\n")
    }, error = function(e2) {
      cat("  - También falló método alternativo:", e2$message, "\n")
    })
  })
  
  # 3. CALCULAR SEÑAL FILOGENÉTICA
  phylo_signal <- "NA"
  
  tryCatch({
    cat("  - Verificando efectos aleatorios...\n")
    
    # Verificar que hay efectos aleatorios
    ranef_summary <- ranef(model)
    if(length(ranef_summary) == 0) {
      stop("No hay efectos aleatorios en el modelo")
    }
    
    cat("  - Efectos aleatorios encontrados:", names(ranef_summary), "\n")
    
    hyp <- paste(
      "sd_label__Intercept^2 /",
      "(sd_label__Intercept^2 + sigma^2) = 0"
    )
    
    hyp_result <- hypothesis(model, hyp, class = NULL)
    cat("  - Hipótesis evaluada exitosamente\n")
    
    lambda_estimate <- round(hyp_result$hypothesis$Estimate, 3)
    lambda_lower <- round(hyp_result$hypothesis$CI.Lower, 3)
    lambda_upper <- round(hyp_result$hypothesis$CI.Upper, 3)
    phylo_signal <- paste0(lambda_estimate, " [", lambda_lower, ", ", lambda_upper, "]")
    
  }, error = function(e) {
    cat("Error en filogenia para", model_name, ":", e$message, "\n")
    
    # Intentar diagnóstico de componentes
    tryCatch({
      cat("  - Intentando diagnóstico manual...\n")
      posterior_samples <- as_draws_df(model)
      sd_cols <- names(posterior_samples)[str_detect(names(posterior_samples), "sd_")]
      sigma_cols <- names(posterior_samples)[str_detect(names(posterior_samples), "^sigma$")]
      cat("  - Columnas SD encontradas:", paste(sd_cols, collapse = ", "), "\n")
      cat("  - Columnas sigma encontradas:", paste(sigma_cols, collapse = ", "), "\n")
    }, error = function(e2) {
      cat("  - También falló diagnóstico manual:", e2$message, "\n")
    })
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
    `Coefficients` = coefs_summary
  ))
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
hwi_model_names <- names(univariate_models)[str_detect(names(univariate_models), "log_HWI.*long.*average")]
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

# COMBINAR RESULTADOS ====================================================

if(length(model_results) > 0) {
  final_table <- bind_rows(model_results) %>%
    arrange(factor(Age, levels = c("average", "breeding", "natal")), 
            factor(`Dispersal Type`, levels = c("Median dispersal", "Long-distance dispersal")))
  
  # AÑADIR RANKINGS DE IMPORTANCIA
  cat("=== AÑADIENDO RANKINGS ===\n")
  
  final_table <- final_table %>%
    rowwise() %>%
    mutate(
      `Variable Ranking` = case_when(
        Age == "average" & `Dispersal Type` == "Median dispersal" ~ 
          get_importance_ranking(median_var_sel, "average", "Median dispersal"),
        Age == "average" & `Dispersal Type` == "Long-distance dispersal" ~ 
          "HWI (1)",  # Solo HWI para average long
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
  
  # MOSTRAR Y GUARDAR RESULTADOS
  cat("=== TABLA FINAL ===\n")
  print(final_table)
  
  # Guardar tabla
  write_csv(final_table, "results/combined/comprehensive_model_summary.csv")
  
  # Mostrar versión legible
  cat("\n=== VERSIÓN LEGIBLE ===\n")
  for(i in 1:nrow(final_table)) {
    cat(paste(rep("=", 80), collapse = ""), "\n")
    cat("MODEL:", final_table$Model[i], "\n")
    cat("Age:", final_table$Age[i], "| Type:", final_table$`Dispersal Type`[i], "\n")
    cat("Variables:", final_table$`N Variables`[i], "\n")
    cat("R² Marginal:", final_table$`R2 Marginal`[i], "| R² Conditional:", final_table$`R2 Conditional`[i], "\n")
    cat("Phylogenetic Signal:", final_table$`Phylogenetic Signal`[i], "\n")
    cat("Variable Ranking:", str_wrap(final_table$`Variable Ranking`[i], width = 70), "\n")
    cat("Coefficients:", str_wrap(final_table$Coefficients[i], width = 70), "\n\n")
  }
  
} else {
  cat("No se pudieron procesar los modelos\n")
}

