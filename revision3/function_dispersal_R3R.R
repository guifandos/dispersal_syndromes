# MODULAR DISPERSAL ANALYSIS WORKFLOW =======================================
## Packages ---------------------------------------------------------------
install.load.package <- function(x) {
  if (!require(x, character.only = TRUE))
    install.packages(x, repos='http://cran.us.r-project.org')
  require(x, character.only = TRUE)
}
package_vec <- c(
  "readr","dplyr", "caret","tidyr", "tidyverse","mgcv", "MuMIn", "purrr", "reshape2", "lattice", "car", "ape", "geiger", "phytools", "nlme", "raster", "ggplot2", "sjPlot", "MCMCglmm", "plotMCMC", "tidybayes" , "plotMCMC", "loo", "brms",
  "mice", "projpred","geiger", "caper", "phylolm", "knitr", "ggmice", "picante", "broom", "performance", 
  "DHARMa", "DHARMa.helpers"# names of the packages required placed here as character objects
)

sapply(package_vec, install.load.package)

# Libraries
library(tidyverse)
library(ape)
library(brms)
library(projpred)
library(parallel)
library(doParallel)
library(foreach)
library(knitr)

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

# MAIN WORKFLOW FUNCTION =====================================================

run_dispersal_analysis <- function(dispersal_type = "average", 
                                   save_results = TRUE,
                                   verbose = TRUE) {
  
  if (verbose) cat("Starting dispersal analysis for type:", dispersal_type, "\n")
  
  # LOADING DATA ============================================================
  
  if (verbose) cat("Loading data...\n")
  
  ## Load phylogenetic tree
  rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")
  
  ## Load dispersal and trait data 
  dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
  distance_total_functions <- read_csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 
  
  # DATA MANIPULATION ======================================================
  
  if (verbose) cat("Processing data for type:", dispersal_type, "\n")
  
  ## Prepare main dispersal traits dataset
  dispersal_traits <- dispersal_traits_total %>% 
    filter(type == dispersal_type) %>% 
    dplyr::rename(label = species) %>% 
    distinct(label, .keep_all = TRUE) %>%
    # Standardize species names (replace spaces with underscores)
    mutate(label = str_replace_all(label, " ", "_"))
  
  ## Prepare distance functions dataset
  distance_total_functions_join <- distance_total_functions %>%
    filter(type == dispersal_type) %>% 
    dplyr::select(species, median, upper_distance, function_id) %>% 
    myspread(function_id, c(median, upper_distance)) %>% 
    dplyr::rename(label = species) %>%
    # Standardize species names
    mutate(label = str_replace_all(label, " ", "_")) %>%
    # Add log-transformed variables
    mutate(
      across(c(Exponential_median, Exponential_upper_distance, 
               Gamma_median, Gamma_upper_distance,
               Hcauchy_median, Hcauchy_upper_distance,
               Weibull_median, Weibull_upper_distance), 
             ~ log(.x + 1), 
             .names = "{.col}_log")
    )
  
  # SPECIES NAME CORRECTIONS ===============================================
  
  ## Correct taxonomic mismatches
  dispersal_traits <- dispersal_traits %>%
    mutate(label = case_when(
      label == "Apus_melba" ~ "Tachymarptis_melba",
      label == "Chlidonias_hybridus" ~ "Chlidonias_hybrida",
      label == "Delichon_urbica" ~ "Delichon_urbicum",
      label == "Mergus_albellus" ~ "Mergellus_albellus",
      label == "Saxicola_torquata" ~ "Saxicola_torquatus",
      label == "Tetrao tetrix" ~ "Lyrurus tetrix",
      label == "Stercorarius_skua" ~ "Catharacta_skua",
      TRUE ~ label
    ))
  
  # Set rownames for phylogenetic analysis
  row.names(dispersal_traits) <- dispersal_traits$label
  
  # PHYLOGENETIC DATA MATCHING =============================================
  
  if (verbose) cat("Matching species with phylogeny...\n")
  
  ## Find species present in both datasets and phylogeny
  name_dispersal <- unique(dispersal_traits$label)
  names_tree <- rf.tree$tip.label
  species_name <- intersect(name_dispersal, names_tree)
  
  ## Filter datasets to matched species
  dispersal_traits <- dispersal_traits %>% 
    filter(label %in% species_name)
  
  ## Prepare phylogenetic tree
  dispersal_tree <- keep.tip(rf.tree, species_name)
  dispersal_tree <- compute.brlen(dispersal_tree, method = "Grafen")
  
  ## Check name matching
  name_check_result <- name.check(dispersal_tree, dispersal_traits)
  if (!is.character(name_check_result)) {
    if (verbose) cat("All species names match between tree and data\n")
  } else {
    if (verbose) cat("Warning: Name mismatches detected\n")
  }
  
  # VARIABLE PREPARATION ===================================================
  
  if (verbose) cat("Preparing variables...\n")
  
  ## Create analysis dataset with renamed variables
  dispersal_analysis <- dispersal_traits %>% 
    # Rename variables for clarity
    dplyr::rename(
      habita_for = Habitat.niche.position.along.forest.open.area.gradient,
      habitat_niche_breadth = Habitat.niche.breadth,
      humid_grad = Position.along.humidity.gradient,
      human_set = Position.along.humidity.gradient,
      diet = Diet.niche.position,
      LHS = Life.history.strategy,
      climatic_pos = Climatic.niche.position...C.,
      climatic_breadth = Climatic.niche.breadth...C.,
      range_size = Breeding.range.size..km2.,
      distance_mig = Migration.distance..km.,
      body_mass = Body.mass..log.
    ) %>% 
    # Create response variables
    mutate(
      dispersal_distance = log(median + 1),
      lon_dispersal_distance = log(upper_distance + 1)
    ) %>% 
    # Select relevant variables for analysis
    dplyr::select(
      dispersal_distance, lon_dispersal_distance, body_mass, HWI, 
      distance_mig, PC1, PC2, Latitude, AnnualTemp, TempRange, 
      AnnualPrecip, PrecipRange
    ) %>%
    # Add log-transformed variables
    mutate(
      log_body_mass = log(body_mass + 1),
      log_HWI = log(HWI + 1)
    ) %>% 
    # Convert character columns to numeric
    modify_if(is.character, as.numeric) %>% 
    # Scale all variables
    scale() %>% 
    as.data.frame()
  
  ## Add categorical variables (not scaled)
  categorical_vars <- dispersal_traits %>%
    dplyr::select(
      Territoriality.x, Migration.1, Diet.niche.position,
      Habitat.niche.position.along.forest.open.area.gradient,
      Habitat.niche.breadth
    ) %>%
    dplyr::rename(
      territoriality = Territoriality.x,
      migration_status = Migration.1,
      diet = Diet.niche.position,
      habita_for = Habitat.niche.position.along.forest.open.area.gradient,
      habitat_niche_breadth = Habitat.niche.breadth
    )
  
  ## Combine scaled and categorical variables
  dispersal_analysis <- cbind(
    label = dispersal_traits$label,
    dispersal_analysis,
    categorical_vars
  )
  
  # HANDLE MISSING DATA ====================================================
  
  if (verbose) cat("Handling missing data...\n")
  
  ## Create complete cases dataset
  dispersal_analysis_partial <- dispersal_analysis %>% 
    drop_na()
  
  ## Update phylogenetic tree to match complete dataset
  name_check_partial <- name.check(dispersal_tree, dispersal_analysis_partial)
  
  if (length(name_check_partial$tree_not_data) > 0) {
    dispersal_tree_partial <- drop.tip(dispersal_tree, name_check_partial$tree_not_data)
  } else {
    dispersal_tree_partial <- dispersal_tree
  }
  
  ## Ensure data and tree match perfectly
  matches <- match(dispersal_analysis_partial$label, dispersal_tree_partial$tip.label, nomatch = 0)
  dispersal_analysis_partial <- dispersal_analysis_partial[matches != 0, ]
  row.names(dispersal_analysis_partial) <- dispersal_analysis_partial$label
  
  ## Final check
  final_check <- name.check(dispersal_tree_partial, dispersal_analysis_partial)
  if (!is.character(final_check)) {
    if (verbose) cat("Final dataset ready with", nrow(dispersal_analysis_partial), "species\n")
  }
  
  # MERGE WITH DISTANCE FUNCTIONS ==========================================
  
  ## Add distance function estimates
  dispersal_analysis_partial <- left_join(
    dispersal_analysis_partial, 
    distance_total_functions_join, 
    by = "label"
  )
  
  # MODEL SETUP =============================================================
  
  if (verbose) cat("Setting up models...\n")
  
  ## Define predictor sets

  st <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")
  
  ## Prepare data and phylogeny for modeling
  data <- dispersal_analysis_partial
  A <- ape::vcv.phylo(dispersal_tree_partial)
  
  # RUN MODELS ==============================================================
  
  results <- run_bayesian_models(data, A, st, st_full_log, dispersal_type, verbose)
  
  # SAVE RESULTS ============================================================
  
  if (save_results) {
    save_model_results(results, dispersal_type, verbose)
  }
  
  if (verbose) cat("Analysis completed for type:", dispersal_type, "\n")
  
  return(results)
}

# BAYESIAN MODELING FUNCTION =================================================

run_bayesian_models <- function(data, A, st, st_full_log, dispersal_type, verbose = TRUE) {
  
  # MEDIAN DISPERSAL MODELS ================================================
  
  if (verbose) cat("Fitting median dispersal model (without phylogeny)...\n")
  
  ## Base model without phylogeny
  model_dispersal_gauss_log_complete <- brm(
    paste0("Weibull_median_log ~ ", paste(st, collapse = " + ")),
    data = data, 
    family = gaussian(), 
    data2 = list(A = A),
    chains = 2, 
    cores = 2, 
    iter = 4000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE)
  )
  
  ## Add model comparison criteria
  model_dispersal_gauss_log_complete <- add_criterion(
    model_dispersal_gauss_log_complete, 
    "loo", 
    moment_match = TRUE
  )
  
  # VARIABLE SELECTION FOR MEDIAN DISPERSAL ===============================
  
  if (verbose) cat("Running variable selection for median dispersal...\n")
  
  ## K-fold cross-validation setup
  cv_fits <- run_cvfun(model_dispersal_gauss_log_complete, K = 10)
  
  ## Parallel setup for cross-validation
  ncores <- parallel::detectCores(logical = FALSE)
  doParallel::registerDoParallel(ncores)
  
  ## Cross-validated variable selection
  cvvs <- cv_varsel(
    model_dispersal_gauss_log_complete,
    cv_method = "kfold",
    cvfits = cv_fits,
    method = "forward",
    nclusters_pred = 20,
    nterms_max = 11,
    parallel = TRUE,
    verbose = verbose
  )
  
  ## Cleanup parallel processing
  doParallel::stopImplicitCluster()
  foreach::registerDoSEQ()
  
  ## Variable selection results
  size_decided <- suggest_size(cvvs, stat = "mlpd")
  rk <- ranking(cvvs)
  predictors_final <- head(rk[["fulldata"]], size_decided)
  
  if (verbose) cat("Selected predictors for median dispersal:", paste(predictors_final, collapse = ", "), "\n")
  
  ## Alternative variable selection without CV
  vs1 <- varsel(model_dispersal_gauss_log_complete, 
                method = "forward", 
                cv_method = "loo", 
                refit_prj = TRUE)
  
  modsize_decided <- suggest_size(vs1)
  soltrms <- solution_terms(vs1)
  soltrms_final <- head(soltrms, modsize_decided) 
  
  # FINAL MEDIAN DISPERSAL MODELS ==========================================
  
  if (verbose) cat("Fitting final median dispersal models with phylogeny...\n")
  
  ## Complete model with phylogeny
  model_dispersal_gauss_log_phylo_complete <- brm(
    paste0("Weibull_median_log ~ ", paste(st, collapse = " + "), " + (1|gr(label, cov = A))"),
    data = data, 
    family = gaussian(), 
    data2 = list(A = A),
    chains = 2, 
    cores = 2, 
    iter = 4000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE)
  )
  
  model_dispersal_gauss_log_phylo_complete <- add_criterion(
    model_dispersal_gauss_log_phylo_complete, 
    "loo", 
    moment_match = TRUE, 
    reloo = TRUE
  )
  
  ## Reduced model with selected predictors
  model_dispersal_gauss_log_phylo_short <- brm(
    paste0("Weibull_median_log ~ ", paste(predictors_final, collapse = " + "), " + (1|gr(label, cov = A))"),
    data = data, 
    family = gaussian(), 
    data2 = list(A = A),
    chains = 2, 
    cores = 2, 
    iter = 4000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE)
  )
  
  model_dispersal_gauss_log_phylo_short <- add_criterion(
    model_dispersal_gauss_log_phylo_short, 
    "loo", 
    moment_match = TRUE, 
    reloo = TRUE
  )
  
  # LONG DISTANCE DISPERSAL MODELS =========================================
  
  if (verbose) cat("Fitting long distance dispersal model...\n")
  
  ## Base model for long distance dispersal
  model_long_dispersal_gauss_log_complete <- brm(
    paste0("Weibull_upper_distance_log ~ ", paste(st, collapse = " + ")),
    data = data, 
    family = gaussian(), 
    data2 = list(A = A),
    chains = 2, 
    cores = 2, 
    iter = 4000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE)
  )
  
  model_long_dispersal_gauss_log_complete <- add_criterion(
    model_long_dispersal_gauss_log_complete, 
    "loo", 
    moment_match = TRUE
  )
  
  # VARIABLE SELECTION FOR LONG DISTANCE DISPERSAL ========================
  
  if (verbose) cat("Running variable selection for long distance dispersal...\n")
  
  ## K-fold cross-validation for long distance model
  cv_fits_long <- run_cvfun(model_long_dispersal_gauss_log_complete, K = 10)
  
  ## Parallel setup
  doParallel::registerDoParallel(ncores)
  
  ## Cross-validated variable selection for long distance
  cvvs_long <- cv_varsel(
    model_long_dispersal_gauss_log_complete,
    cv_method = "kfold",
    cvfits = cv_fits_long,
    method = "forward",
    nclusters_pred = 20,
    nterms_max = 11,
    parallel = TRUE,
    verbose = verbose
  )
  
  ## Cleanup
  doParallel::stopImplicitCluster()
  foreach::registerDoSEQ()
  
  ## Results for long distance dispersal
  size_decided_long <- suggest_size(cvvs_long, stat = "mlpd")
  rk_long <- ranking(cvvs_long)
  predictors_final_long <- head(rk_long[["fulldata"]], size_decided_long)
  
  if (verbose) cat("Selected predictors for long distance dispersal:", paste(predictors_final_long, collapse = ", "), "\n")
  
  ## Alternative variable selection without CV
  vs1_long <- varsel(model_long_dispersal_gauss_log_complete, 
                     method = "forward", 
                     cv_method = "loo", 
                     refit_prj = TRUE)
  
  modsize_decided_long <- suggest_size(vs1_long)
  soltrms_long <- solution_terms(vs1_long)
  soltrms_final_long <- head(soltrms_long, modsize_decided_long) 
  
  # FINAL LONG DISTANCE DISPERSAL MODELS ==================================
  
  if (verbose) cat("Fitting final long distance dispersal models with phylogeny...\n")
  
  ## Complete model with phylogeny
  model_dispersal_long_phylo_complete <- brm(
    paste0("Weibull_upper_distance_log ~ ", paste(st, collapse = " + "), " + (1|gr(label, cov = A))"),
    data = data, 
    family = gaussian(), 
    data2 = list(A = A),
    chains = 2, 
    cores = 2, 
    iter = 4000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE)
  )
  
  model_dispersal_long_phylo_complete <- add_criterion(
    model_dispersal_long_phylo_complete, 
    "loo", 
    moment_match = TRUE, 
    reloo = TRUE
  )
  
  ## Reduced model with selected predictors
  model_dispersal_long_phylo_short <- brm(
    paste0("Weibull_upper_distance_log ~ ", paste(predictors_final_long, collapse = " + "), " + (1|gr(label, cov = A))"),
    data = data, 
    family = gaussian(), 
    data2 = list(A = A),
    chains = 2, 
    cores = 2, 
    iter = 4000,
    control = list(adapt_delta = 0.95),
    save_pars = save_pars(all = TRUE)
  )
  
  model_dispersal_long_phylo_short <- add_criterion(
    model_dispersal_long_phylo_short, 
    "loo", 
    moment_match = TRUE, 
    reloo = TRUE
  )
  
  # COMPILE RESULTS =========================================================
  
  results <- list(
    # Models
    median_complete = model_dispersal_gauss_log_phylo_complete,
    median_short = model_dispersal_gauss_log_phylo_short,
    long_complete = model_dispersal_long_phylo_complete,
    long_short = model_dispersal_long_phylo_short,
    
    # Variable selection
    cvvs_median = cvvs,
    cvvs_long = cvvs_long,
    vs1_median = vs1,
    vs1_long = vs1_long,
    
    # Selected predictors
    predictors_median = predictors_final,
    predictors_long = predictors_final_long,
    
    # Data
    data = data,
    phylogeny = list(A = A)
  )
  
  return(results)
}

# SAVE RESULTS FUNCTION ======================================================

save_model_results <- function(results, dispersal_type, verbose = TRUE) {
  
  if (verbose) cat("Saving model results for type:", dispersal_type, "\n")
  
  ## Create results directory if it doesn't exist
  results_dir <- paste0("results/revision3/", dispersal_type)
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  ## Save short models
  results_short <- tibble(
    model = list(results$median_short, results$long_short), 
    variable_selection = list(results$cvvs_median, results$cvvs_long), 
    variable_selection_old = list(results$vs1_median, results$vs1_long),
    age = dispersal_type, 
    function_t = "short_model", 
    type = c("median", "long")
  )
  
  save(results_short, file = file.path(results_dir, paste0("model_", dispersal_type, "_short_R3.RData")))
  
  ## Save complete models
  results_complete <- tibble(
    model = list(results$median_complete, results$long_complete), 
    age = dispersal_type, 
    function_t = "complete_model", 
    type = c("median", "long")
  )
  
  save(results_complete, file = file.path(results_dir, paste0("model_", dispersal_type, "_complete_R3.RData")))
  
  if (verbose) cat("Results saved in", results_dir, "\n")
}

# MAIN EXECUTION ==============================================================

# Run analysis for average dispersal
cat("=== RUNNING ANALYSIS FOR AVERAGE DISPERSAL ===\n")
results_average <- run_dispersal_analysis("average", save_results = TRUE, verbose = TRUE)

# Run analysis for breeding dispersal
cat("\n=== RUNNING ANALYSIS FOR BREEDING DISPERSAL ===\n")
results_breeding <- run_dispersal_analysis("breeding", save_results = TRUE, verbose = TRUE)

# Run analysis for natal dispersal
cat("\n=== RUNNING ANALYSIS FOR BREEDING DISPERSAL ===\n")
results_breeding <- run_dispersal_analysis("natal", save_results = TRUE, verbose = TRUE)


# Optional: Run for other types if they exist
# results_natal <- run_dispersal_analysis("natal", save_results = TRUE, verbose = TRUE)

cat("\n=== ALL ANALYSES COMPLETED ===\n")
cat("Results saved for both average and breeding dispersal types.\n")
