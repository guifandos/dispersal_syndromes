# ==============================================================================
# PREPARACIÓN DE DATOS Y MATRICES FILOGENÉTICAS PARA MODELOS DE DISPERSIÓN
# ==============================================================================
#
# DESCRIPCIÓN:
# Este script prepara datasets limpios y matrices filogenéticas para análisis
# de dispersión en aves europeas. Procesa tres tipos de dispersión: average 
# (total), natal, y breeding dispersal.
#
# FUNCIONES PRINCIPALES:
# 1. prepare_dispersal_data(): Función principal que:
#    - Filtra datos por tipo de dispersión
#    - Limpia y estandariza nombres de especies
#    - Crea transformaciones logarítmicas
#    - Elimina NAs y ajusta árboles filogenéticos
#    - Genera matrices de covarianza filogenética
# 2. verify_dataset(): Verifica integridad de datos preparados
#
# INPUTS REQUERIDOS:
# - data/data_philo/dispersal_tree500.nex (árbol filogenético)
# - data/data_process/dispersal_traits_total_PCA_20220424.csv (traits)
# - data/dispersal_distance/Table_S13_ species_dispersal_distances.csv (distancias)
#
# OUTPUTS GENERADOS:
# - data/processed/dispersal_average_complete.RData (data_average, tree_average, A_average)
# - data/processed/dispersal_natal_complete.RData (data_natal, tree_natal, A_natal)
# - data/processed/dispersal_breeding_complete.RData (data_breeding, tree_breeding, A_breeding)
#
# DEPENDENCIAS:
# library(dplyr)      # Manipulación de datos
# library(tidyr)      # Pivoteo y limpieza
# library(ape)        # Análisis filogenéticos
# library(readr)      # Lectura de archivos CSV
# library(purrr)      # Programación funcional (modify_if)
#
# VARIABLES CLAVE GENERADAS:
# - body_mass, log_HWI: Variables morfológicas
# - PC1, PC2: Componentes principales de life-history
# - habita_for, diet: Variables ecológicas categóricas
# - distance_mig, Latitude: Variables geográficas
# - Weibull_median_log, Weibull_upper_distance_log: Variables respuesta
#
# NOTAS:
# - Se eliminan especies con NAs en variables clave
# - Se corrigen nombres problemáticos de especies
# - Los datos se estandarizan (mean=0, sd=1) excepto variables categóricas
# - Las matrices filogenéticas (A) se calculan con vcv.phylo()
#
# AUTOR: Guillermo Fandos
# FECHA: 19/09/2025
# ==============================================================================

# Cargar librerías requeridas
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr) 
  library(ape)
  library(readr)
  library(purrr)
})

# Crear directorios si no existen
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/weibull/natal", recursive = TRUE, showWarnings = FALSE)
dir.create("results/weibull/breeding", recursive = TRUE, showWarnings = FALSE)

# PREPARACIÓN DE DATOS PARA NATAL, BREEDING Y AVERAGE ==================

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

foo <- function(file_names) {
  lapply(file_names, load, environment())
  ls()
}

# Función para preparar datos por tipo de dispersión
prepare_dispersal_data <- function(dispersal_traits_total, distance_total_functions, 
                                   rf_tree, dispersal_type = "average") {
  
  cat("Preparando datos para:", dispersal_type, "\n")
  
  # Filtrar por tipo de dispersión
  dispersal_traits_filtered <- dispersal_traits_total %>% 
    filter(type == dispersal_type) %>% 
    dplyr::rename(label = species) %>% 
    distinct(label, .keep_all = TRUE) 
  
  # Preparar datos de distancia
  distance_filtered <- distance_total_functions %>%
    filter(type == dispersal_type) %>% 
    dplyr::select(species, median, upper_distance, function_id) %>% 
    myspread(function_id, c(median, upper_distance)) %>% 
    dplyr::rename(label = species)
  
  # Limpiar nombres de especies
  dispersal_traits_filtered$label <- gsub(" ", "_", dispersal_traits_filtered$label)
  distance_filtered$label <- gsub(" ", "_", distance_filtered$label)
  
  # Crear columnas log-transformadas
  distance_filtered <- distance_filtered %>% 
    mutate(
      Exponential_median_log = log(Exponential_median + 1),
      Exponential_upper_distance_log = log(Exponential_upper_distance + 1),
      Gamma_median_log = log(Gamma_median + 1),
      Gamma_upper_distance_log = log(Gamma_upper_distance + 1),
      Hcauchy_median_log = log(Hcauchy_median + 1),
      Hcauchy_upper_distance_log = log(Hcauchy_upper_distance + 1),
      Weibull_median_log = log(Weibull_median + 1),
      Weibull_upper_distance_log = log(Weibull_upper_distance + 1)
    )
  
  # Corregir nombres de especies problemáticos
  dispersal_traits_filtered <- dispersal_traits_filtered %>%
    mutate(label = dplyr::recode(label, 
                                 "Apus_melba" = "Tachymarptis_melba",
                                 "Chlidonias_hybridus" = "Chlidonias_hybrida",          
                                 "Delichon_urbica" = "Delichon_urbicum",           
                                 "Mergus_albellus" = "Mergellus_albellus",
                                 "Saxicola_torquata" = "Saxicola_torquatus",
                                 "Tetrao_tetrix" = "Lyrurus_tetrix",             
                                 "Stercorarius_skua" = "Catharacta_skua" 
    ))
  
  row.names(dispersal_traits_filtered) <- dispersal_traits_filtered$label
  
  # Encontrar especies en común entre datos y árbol
  name_dispersal <- unique(dispersal_traits_filtered$label)
  names_tree <- rf_tree$tip.label
  species_common <- base::intersect(name_dispersal, names_tree)
  
  cat("Especies en común:", length(species_common), "\n")
  
  # Filtrar datos por especies en común
  dispersal_traits_filtered <- dispersal_traits_filtered %>% 
    filter(label %in% species_common)
  
  # Crear árbol específico para este tipo
  tree_filtered <- keep.tip(rf_tree, species_common)
  tree_filtered <- compute.brlen(tree_filtered, method = "Grafen")
  
  # Verificar que coincidan
  check_result <- name.check(tree_filtered, dispersal_traits_filtered)
  if (class(check_result) != "character") {
    warning("Hay problemas de coincidencia entre árbol y datos para ", dispersal_type)
    print(check_result)
  }
  
  # Crear dataset de análisis
  dispersal_analysis <- dispersal_traits_filtered %>% 
    dplyr::rename(
      "habita_for" = "Habitat.niche.position.along.forest.open.area.gradient",
      "habitat_niche_breadth" = "Habitat.niche.breadth",
      "humid_grad" = "Position.along.humidity.gradient",
      "human_set" = "Position.along.humidity.gradient",
      "diet" = "Diet.niche.position",
      "LHS" = "Life.history.strategy",
      "climatic_pos" = "Climatic.niche.position...C.",
      "climatic_breadth" = "Climatic.niche.breadth...C.",
      "range_size" = "Breeding.range.size..km2.",
      "distance_mig" = "Migration.distance..km.",
      "body_mass" = "Body.mass..log."
    ) %>% 
    mutate(
      dispersal_distance = log(median + 1),
      lon_dispersal_distance = log(upper_distance + 1) 
    ) %>% 
    dplyr::select(
      dispersal_distance, lon_dispersal_distance, body_mass, HWI, distance_mig, 
      PC1, PC2, Latitude, AnnualTemp, TempRange, AnnualPrecip, PrecipRange
    ) 
  
  # Transformaciones y estandarización
  dispersal_analysis <- dispersal_analysis %>% 
    mutate(
      log_body_mass = log(body_mass + 1),
      log_HWI = log(HWI + 1)
    ) %>% 
    modify_if(., is.character, as.numeric) %>% 
    scale(.) %>% 
    as.data.frame(.) 
  
  # Añadir variables categóricas
  label <- dispersal_traits_filtered$label
  dispersal_analysis <- cbind(
    dispersal_analysis, 
    dispersal_traits_filtered[, c("Territoriality.x", "Migration.1", "Diet.niche.position", 
                                  "Habitat.niche.position.along.forest.open.area.gradient", 
                                  "Habitat.niche.breadth")]
  ) %>% 
    dplyr::rename(
      "habita_for" = "Habitat.niche.position.along.forest.open.area.gradient",
      "habitat_niche_breadth" = "Habitat.niche.breadth",
      "diet" = "Diet.niche.position", 
      "territoriality" = "Territoriality.x"
    )
  
  dispersal_analysis <- cbind(label, dispersal_analysis)
  
  # Eliminar NAs
  dispersal_analysis_clean <- dispersal_analysis %>% 
    drop_na()
  
  cat("Especies después de eliminar NAs:", nrow(dispersal_analysis_clean), "\n")
  
  # Ajustar árbol a datos limpios
  matches <- match(dispersal_analysis_clean$label, tree_filtered$tip.label, nomatch = 0)
  dispersal_analysis_clean <- subset(dispersal_analysis_clean, matches != 0)
  row.names(dispersal_analysis_clean) <- dispersal_analysis_clean$label
  
  # Eliminar especies del árbol que no están en datos
  ff <- name.check(tree_filtered, dispersal_analysis_clean)
  if (class(ff) != "character") {
    to_drop <- ff$tree_not_data
    tree_clean <- drop.tip(tree_filtered, to_drop)
  } else {
    tree_clean <- tree_filtered
  }
  
  # Verificación final
  final_check <- name.check(tree_clean, dispersal_analysis_clean)
  if (class(final_check) != "character") {
    stop("Error: árbol y datos no coinciden para ", dispersal_type)
  }
  
  # Añadir datos de distancia
  dispersal_analysis_clean <- left_join(dispersal_analysis_clean, distance_filtered, by = "label")
  
  # Crear matriz de covarianza filogenética
  A <- ape::vcv.phylo(tree_clean)
  
  cat("Dataset final para", dispersal_type, ":", nrow(dispersal_analysis_clean), "especies\n")
  
  return(list(
    data = dispersal_analysis_clean,
    tree = tree_clean,
    A = A,
    type = dispersal_type
  ))
}

# APLICAR FUNCIÓN PARA CADA TIPO ========================================

# Cargar datos originales (asume que ya están cargados)
rf_tree <- read.nexus("data/data_philo/dispersal_tree500.nex")
dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 

# Preparar datos para cada tipo
cat("=== PREPARANDO DATOS PARA TODOS LOS TIPOS ===\n")

data_average_obj <- prepare_dispersal_data(
  dispersal_traits_total, distance_total_functions, rf_tree, "average"
)

data_natal_obj <- prepare_dispersal_data(
  dispersal_traits_total, distance_total_functions, rf_tree, "natal"
)

data_breeding_obj <- prepare_dispersal_data(
  dispersal_traits_total, distance_total_functions, rf_tree, "breeding"
)

# EXTRAER COMPONENTES ===================================================

# Average
data_average <- data_average_obj$data
tree_average <- data_average_obj$tree
A_average <- data_average_obj$A

# Natal
data_natal <- data_natal_obj$data
tree_natal <- data_natal_obj$tree
A_natal <- data_natal_obj$A

# Breeding
data_breeding <- data_breeding_obj$data
tree_breeding <- data_breeding_obj$tree
A_breeding <- data_breeding_obj$A

# RESUMEN DE DATASETS ===================================================

cat("\n=== RESUMEN DE DATASETS ===\n")
cat("Average dispersal:", nrow(data_average), "especies\n")
cat("Natal dispersal:", nrow(data_natal), "especies\n")
cat("Breeding dispersal:", nrow(data_breeding), "especies\n")

# Verificar columnas disponibles
cat("\nColumnas disponibles:\n")
print(names(data_average))

# Especies en común entre datasets
common_all <- Reduce(intersect, list(
  data_average$label, 
  data_natal$label, 
  data_breeding$label
))
cat("\nEspecies en común entre todos los tipos:", length(common_all), "\n")

# GUARDAR DATASETS PREPARADOS ===========================================

save(data_average, tree_average, A_average, 
     file = "data/processed/dispersal_average_complete.RData")

save(data_natal, tree_natal, A_natal, 
     file = "data/processed/dispersal_natal_complete.RData")

save(data_breeding, tree_breeding, A_breeding, 
     file = "data/processed/dispersal_breeding_complete.RData")

# VERIFICACIONES FINALES ================================================

# Función para verificar dataset
verify_dataset <- function(data, tree, A, type) {
  cat("\n=== VERIFICACIÓN", toupper(type), "===\n")
  cat("Filas en datos:", nrow(data), "\n")
  cat("Tips en árbol:", length(tree$tip.label), "\n")
  cat("Dimensiones de matriz A:", dim(A), "\n")
  
  # Verificar coincidencia
  check <- name.check(tree, data)
  if (class(check) == "character") {
    cat("✓ Árbol y datos coinciden perfectamente\n")
  } else {
    cat("✗ Problemas de coincidencia:\n")
    print(check)
  }
  
  # Verificar NAs en variables clave
  key_vars <- c("body_mass", "log_HWI", "PC1", "habita_for", "diet", "distance_mig", "Latitude")
  nas_count <- sapply(data[key_vars], function(x) sum(is.na(x)))
  if (sum(nas_count) == 0) {
    cat("✓ No hay NAs en variables clave\n")
  } else {
    cat("✗ NAs encontrados:\n")
    print(nas_count[nas_count > 0])
  }
  
  # Verificar columnas de dispersión
  dispersal_cols <- c("Weibull_median_log", "Weibull_upper_distance_log")
  missing_cols <- dispersal_cols[!dispersal_cols %in% names(data)]
  if (length(missing_cols) == 0) {
    cat("✓ Columnas de dispersión presentes\n")
  } else {
    cat("✗ Columnas faltantes:", paste(missing_cols, collapse = ", "), "\n")
  }
}

# Verificar todos los datasets
verify_dataset(data_average, tree_average, A_average, "average")
verify_dataset(data_natal, tree_natal, A_natal, "natal")
verify_dataset(data_breeding, tree_breeding, A_breeding, "breeding")

cat("\n=== PREPARACIÓN COMPLETADA ===\n")
cat("Los datasets están listos para los modelos univariantes.\n")
cat("Archivos guardados en data/processed/\n")
