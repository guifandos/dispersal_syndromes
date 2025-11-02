# CREAR TABLA ESTILO TABLE MODEL (FORMATO LARGO) ========================
library(tidyverse)
library(knitr)
library(kableExtra)

# Asumir que ya tenemos final_table del código anterior
# Si no, ejecutar primero comprehensive_model_table.R

# FUNCIÓN PARA CREAR TABLA ESTILO MODEL ==================================

create_table_model <- function(final_table) {
  
  # Separar coeficientes en columnas individuales
  coef_table <- final_table %>%
    select(Model, Age, `Dispersal Type`, Coefficients) %>%
    separate_rows(Coefficients, sep = "; ") %>%
    filter(Coefficients != "") %>%
    extract(Coefficients, into = c("Variable", "Estimate_CI", "Significance"), 
            regex = "^([^:]+): ([^\\*]+)(\\*?)$") %>%
    mutate(
      Variable = str_trim(Variable),
      Estimate_CI = str_trim(Estimate_CI),
      # Separar estimate del CI
      Estimate = str_extract(Estimate_CI, "^[\\-0-9\\.]+"),
      CI = str_extract(Estimate_CI, "\\[.+\\]"),
      # Limpiar formato
      Estimate = as.numeric(Estimate),
      CI = str_remove_all(CI, "[\\[\\]]"),
      Significance = ifelse(Significance == "*", "Yes", "No")
    ) %>%
    select(-Estimate_CI)
  
  # Crear tabla principal con métricas del modelo
  model_metrics <- final_table %>%
    select(-Coefficients, -`Variable Ranking`) %>%
    mutate(
      # Crear ID único para cada modelo
      Model_ID = paste0(str_sub(Age, 1, 1), "_", 
                        str_extract(`Dispersal Type`, "^[ML]"), "_",
                        row_number()),
      # Limpiar phylogenetic signal
      Phylo_Lambda = str_extract(`Phylogenetic Signal`, "^[0-9\\.]+"),
      Phylo_CI = str_extract(`Phylogenetic Signal`, "\\[.+\\]"),
      Phylo_Lambda = as.numeric(Phylo_Lambda),
      Phylo_CI = str_remove_all(Phylo_CI, "[\\[\\]]")
    )
  
  # Tabla estilo tabmodel - una fila por modelo
  table_model <- model_metrics %>%
    select(
      `Model` = Model,
      `Age Class` = Age,
      `Dispersal Type` = `Dispersal Type`,
      `N Predictors` = `N Variables`,
      `R² Marginal` = `R2 Marginal`,
      `R² Conditional` = `R2 Conditional`,
      `λ (Phylo Signal)` = Phylo_Lambda,
      `λ 95% CI` = Phylo_CI
    ) %>%
    arrange(
      factor(`Age Class`, levels = c("average", "breeding", "natal")),
      factor(`Dispersal Type`, levels = c("Median dispersal", "Long-distance dispersal"))
    )
  
  return(list(
    model_table = table_model,
    coef_details = coef_table
  ))
}

# FUNCIÓN PARA CREAR TABLA DE COEFICIENTES DETALLADA ====================

create_detailed_coef_table <- function(final_table) {
  
  # Extraer y procesar coeficientes
  detailed_coefs <- final_table %>%
    select(Model, Age, `Dispersal Type`, Coefficients) %>%
    separate_rows(Coefficients, sep = "; ") %>%
    filter(Coefficients != "") %>%
    extract(Coefficients, into = c("Variable", "Estimate_CI", "Significance"), 
            regex = "^([^:]+): ([^\\*]+)(\\*?)$") %>%
    mutate(
      Variable = str_trim(Variable),
      Estimate_CI = str_trim(Estimate_CI),
      # Separar estimate del CI
      Estimate = as.numeric(str_extract(Estimate_CI, "^[\\-0-9\\.]+")),
      CI_Lower = as.numeric(str_extract(Estimate_CI, "(?<=\\[)[\\-0-9\\.]+")),
      CI_Upper = as.numeric(str_extract(Estimate_CI, "(?<=, )[\\-0-9\\.]+(?=\\])")),
      Significant = ifelse(Significance == "*", "***", ""),
      # Crear columnas para cada modelo
      Model_Label = paste0(
        case_when(
          Age == "average" ~ "Avg",
          Age == "breeding" ~ "Breed", 
          Age == "natal" ~ "Natal"
        ),
        "_",
        case_when(
          `Dispersal Type` == "Median dispersal" ~ "Med",
          `Dispersal Type` == "Long-distance dispersal" ~ "Long"
        )
      ),
      # Combinar estimate con significancia
      Est_Sig = paste0(sprintf("%.3f", Estimate), Significant)
    ) %>%
    select(Variable, Model_Label, Estimate, CI_Lower, CI_Upper, Est_Sig, Significant)
  
  # Crear tabla wide con modelos como columnas
  coef_wide <- detailed_coefs %>%
    select(Variable, Model_Label, Est_Sig) %>%
    pivot_wider(names_from = Model_Label, values_from = Est_Sig, values_fill = "—") %>%
    arrange(Variable)
  
  return(coef_wide)
}

# EJECUTAR FUNCIONES =====================================================

if(exists("final_table")) {
  
  # Crear tabla principal estilo model
  model_tables <- create_table_model(final_table)
  main_table <- model_tables$model_table
  
  # Crear tabla detallada de coeficientes  
  coef_table_wide <- create_detailed_coef_table(final_table)
  
  # MOSTRAR TABLA PRINCIPAL
  cat("=== TABLA PRINCIPAL - ESTILO TABLE MODEL ===\n")
  print(main_table)
  
  # MOSTRAR TABLA DE COEFICIENTES
  cat("\n=== TABLA DE COEFICIENTES DETALLADA ===\n")
  print(coef_table_wide)
  
  # CREAR VERSIÓN FORMATEADA CON KABLE
  if(require(kableExtra, quietly = TRUE)) {
    
    # Tabla principal formateada
    main_kable <- main_table %>%
      kable(format = "html", 
            caption = "Model Summary: Dispersal Analysis",
            digits = 3) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      row_spec(0, bold = TRUE) %>%
      column_spec(1, bold = TRUE) %>%
      collapse_rows(columns = 1:2, valign = "top")
    
    # Tabla de coeficientes formateada
    coef_kable <- coef_table_wide %>%
      kable(format = "html",
            caption = "Coefficient Estimates (*** p < 0.05)",
            align = c("l", rep("c", ncol(coef_table_wide) - 1))) %>%
      kable_styling(bootstrap_options = c("striped", "hover", "condensed")) %>%
      row_spec(0, bold = TRUE) %>%
      column_spec(1, bold = TRUE)
    
    # Guardar tablas HTML
    cat("\n=== GUARDANDO TABLAS FORMATEADAS ===\n")
    
    # Guardar tabla principal
    main_kable %>%
      save_kable("results/combined/model_summary_table.html")
    
    # Guardar tabla de coeficientes  
    coef_kable %>%
      save_kable("results/combined/coefficients_table.html")
    
    cat("Tablas HTML guardadas en results/combined/\n")
  }
  
  # GUARDAR VERSIONES CSV
  write_csv(main_table, "results/combined/model_summary_table.csv")
  write_csv(coef_table_wide, "results/combined/coefficients_table.csv")
  
  # CREAR VERSIÓN PARA MANUSCRITO
  manuscript_table <- main_table %>%
    mutate(
      Model = case_when(
        str_detect(Model, "Average Median") ~ "Average (Median)",
        str_detect(Model, "Average Long") ~ "Average (Long-distance)", 
        str_detect(Model, "Breeding Median") ~ "Breeding (Median)",
        str_detect(Model, "Breeding Long") ~ "Breeding (Long-distance)",
        str_detect(Model, "Natal Median") ~ "Natal (Median)",
        str_detect(Model, "Natal Long") ~ "Natal (Long-distance)",
        TRUE ~ Model
      ),
      # Combinar R² en una columna
      `R² (Marg. / Cond.)` = paste0(`R² Marginal`, " / ", `R² Conditional`),
      # Combinar λ con CI
      `Phylogenetic Signal λ [95% CI]` = paste0(`λ (Phylo Signal)`, " [", `λ 95% CI`, "]")
    ) %>%
    select(
      Model,
      `N Predictors`,
      `R² (Marg. / Cond.)`,
      `Phylogenetic Signal λ [95% CI]`
    )
  
  cat("\n=== TABLA PARA MANUSCRITO ===\n")
  print(manuscript_table)
  
  write_csv(manuscript_table, "results/combined/manuscript_table.csv")
  
} else {
  cat("ERROR: No existe 'final_table'. Ejecuta primero comprehensive_model_table.R\n")
}
# FIN DEL CÓDIGO =========================================================