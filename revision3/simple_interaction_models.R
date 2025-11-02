# INDEPENDENT INTERACTION ANALYSIS SCRIPT ===================================
# Simple script to fit models with interactions - no variable selection
# ===========================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(brms)
  library(loo)
  library(dplyr)
  library(tidyr)
  library(ape)
  library(stringr)
})

# CONFIGURATION ==============================================================

# Define priority interactions to test
priority_interactions <- c(
  "body_mass * PC1",           # Size × life strategy
  "body_mass * diet",          # Size × trophic requirements  
  "log_HWI * distance_mig",    # Flight efficiency × migration
  "diet * habita_for",         # Diet × habitat specialization
  "PC1 * diet",               # Life strategy × diet
  "Latitude * distance_mig"    # Latitude × migration
)

# Base predictors
base_predictors <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")

# Model settings
model_settings <- list(
  chains = 2,
  cores = 2, 
  iter = 4000,
  control = list(adapt_delta = 0.98, max_treedepth = 12),  # More conservative
  save_pars = save_pars(all = TRUE)
)

# DATA LOADING FUNCTION ======================================================

load_dispersal_data <- function(dispersal_type = "average") {
  
  cat("=== LOADING DATA FOR", toupper(dispersal_type), "DISPERSAL ===\n")
  
  # Load phylogenetic tree
  rf.tree <- read.nexus("data/data_philo/dispersal_tree500.nex")
  cat("✓ Phylogenetic tree loaded\n")
  
  # Load dispersal and trait data 
  dispersal_traits_total <- read.csv("data/data_process/dispersal_traits_total_PCA_20220424.csv")
  distance_total_functions <- read.csv("data/dispersal_distance/Table_S13_ species_dispersal_distances.csv") 
  cat("✓ Data files loaded\n")
  
  # Prepare main dispersal traits dataset
  dispersal_traits <- dispersal_traits_total %>% 
    filter(type == dispersal_type) %>% 
    rename(label = species) %>% 
    distinct(label, .keep_all = TRUE) %>%
    mutate(label = str_replace_all(label, " ", "_"))
  
  # Prepare distance functions dataset - check column names first
  cat("Checking distance data columns...\n")
  print(head(colnames(distance_total_functions)))
  
  distance_functions <- distance_total_functions %>%
    filter(type == dispersal_type) %>%
    # Check what columns are available
    {
      cols <- colnames(.)
      cat("Available columns:", paste(cols, collapse = ", "), "\n")
      
      # Find the right column names (flexible approach)
      species_col <- cols[grepl("species", cols, ignore.case = TRUE)][1]
      median_col <- cols[grepl("median", cols, ignore.case = TRUE)][1]
      upper_col <- cols[grepl("upper", cols, ignore.case = TRUE)][1]
      function_col <- cols[grepl("function", cols, ignore.case = TRUE)][1]
      
      if (is.na(function_col)) {
        function_col <- cols[grepl("distribution", cols, ignore.case = TRUE)][1]
      }
      
      cat("Found columns:\n")
      cat("  Species:", species_col, "\n")
      cat("  Median:", median_col, "\n") 
      cat("  Upper:", upper_col, "\n")
      cat("  Function:", function_col, "\n")
      
      # Select columns that exist
      cols_to_select <- c(species_col, median_col, upper_col, function_col)
      cols_to_select <- cols_to_select[!is.na(cols_to_select)]
      
      if (length(cols_to_select) < 3) {
        cat("Available columns:\n")
        print(cols)
        stop("Cannot find required columns (species, median, upper distance)")
      }
      
      dplyr::select(., all_of(cols_to_select))
    } %>%
    # Rename to standard names
    rename_with(~ case_when(
      grepl("species", .x, ignore.case = TRUE) ~ "species",
      grepl("median", .x, ignore.case = TRUE) ~ "median",
      grepl("upper", .x, ignore.case = TRUE) ~ "upper_distance", 
      grepl("function|distribution", .x, ignore.case = TRUE) ~ "function_id",
      TRUE ~ .x
    )) %>%
    # Filter for Weibull if function column exists
    {
      if ("function_id" %in% colnames(.)) {
        filter(., function_id == "Weibull" | function_id == "weibull")
      } else {
        .  # Keep all if no function column
      }
    } %>%
    # Clean up
    rename(label = species) %>%
    mutate(
      label = str_replace_all(label, " ", "_"),
      Weibull_median_log = log(median + 1),
      Weibull_upper_distance_log = log(upper_distance + 1)
    ) %>%
    dplyr::select(label, Weibull_median_log, Weibull_upper_distance_log)
  
  # Species name corrections
  dispersal_traits <- dispersal_traits %>%
    mutate(label = case_when(
      label == "Apus_melba" ~ "Tachymarptis_melba",
      label == "Chlidonias_hybridus" ~ "Chlidonias_hybrida",
      label == "Delichon_urbica" ~ "Delichon_urbicum",
      label == "Mergus_albellus" ~ "Mergellus_albellus",
      label == "Saxicola_torquata" ~ "Saxicola_torquatus",
      label == "Tetrao_tetrix" ~ "Lyrurus_tetrix",
      label == "Stercorarius_skua" ~ "Catharacta_skua",
      TRUE ~ label
    ))
  
  # Match species with phylogeny
  names_tree <- rf.tree$tip.label
  species_name <- intersect(dispersal_traits$label, names_tree)
  
  # Filter datasets to matched species
  dispersal_traits <- dispersal_traits %>% 
    filter(label %in% species_name)
  
  # Prepare phylogenetic tree
  dispersal_tree <- keep.tip(rf.tree, species_name)
  dispersal_tree <- compute.brlen(dispersal_tree, method = "Grafen")
  
  # Create analysis dataset - flexible column matching
  cat("=== DEBUGGING COLUMN MATCHING ===\n")
  
  analysis_data <- dispersal_traits %>% 
    inner_join(distance_functions, by = "label") %>%
    # Check available columns
    {
      cols <- colnames(.)
      cat("Total columns after join:", length(cols), "\n")
      cat("First 10 columns:", paste(head(cols, 10), collapse = ", "), "\n")
      cat("Columns containing 'diet':", paste(cols[grepl("diet", cols, ignore.case = TRUE)], collapse = ", "), "\n")
      cat("Columns containing 'niche':", paste(cols[grepl("niche", cols, ignore.case = TRUE)], collapse = ", "), "\n")
      cat("Columns containing 'position':", paste(cols[grepl("position", cols, ignore.case = TRUE)], collapse = ", "), "\n")
      
      # Find columns flexibly
      habitat_col <- cols[grepl("habitat.*forest|forest.*open|forest.*area", cols, ignore.case = TRUE)][1]
      diet_col <- cols[grepl("diet", cols, ignore.case = TRUE)][1]  # Simplified search
      pc1_col <- cols[grepl("^PC1$|life.*history", cols, ignore.case = TRUE)][1]
      mass_col <- cols[grepl("body.*mass|mass.*log", cols, ignore.case = TRUE)][1]
      hwi_col <- cols[grepl("^HWI$", cols, ignore.case = TRUE)][1]
      migration_col <- cols[grepl("migration.*distance", cols, ignore.case = TRUE)][1]
      latitude_col <- cols[grepl("^Latitude$|^latitude$", cols, ignore.case = TRUE)][1]
      
      cat("\nColumn matching results:\n")
      cat("  Habitat:", habitat_col, ifelse(is.na(habitat_col), " (NOT FOUND)", " (FOUND)"), "\n")
      cat("  Diet:", diet_col, ifelse(is.na(diet_col), " (NOT FOUND)", " (FOUND)"), "\n")
      cat("  PC1:", pc1_col, ifelse(is.na(pc1_col), " (NOT FOUND)", " (FOUND)"), "\n") 
      cat("  Body mass:", mass_col, ifelse(is.na(mass_col), " (NOT FOUND)", " (FOUND)"), "\n")
      cat("  HWI:", hwi_col, ifelse(is.na(hwi_col), " (NOT FOUND)", " (FOUND)"), "\n")
      cat("  Migration:", migration_col, ifelse(is.na(migration_col), " (NOT FOUND)", " (FOUND)"), "\n")
      cat("  Latitude:", latitude_col, ifelse(is.na(latitude_col), " (NOT FOUND)", " (FOUND)"), "\n")
      
      # If diet_col is still not found, search more broadly
      if (is.na(diet_col)) {
        cat("\nSearching more broadly for diet-related columns:\n")
        diet_candidates <- cols[grepl("diet|niche|food|trophic", cols, ignore.case = TRUE)]
        cat("  Candidates:", paste(diet_candidates, collapse = ", "), "\n")
        if (length(diet_candidates) > 0) {
          diet_col <- diet_candidates[1]
          cat("  Using:", diet_col, "\n")
        }
      }
      
      # Select available columns
      base_cols <- c("label", "Weibull_median_log", "Weibull_upper_distance_log")
      trait_cols <- c(habitat_col, diet_col, pc1_col, mass_col, hwi_col, migration_col, latitude_col)
      trait_cols <- trait_cols[!is.na(trait_cols)]
      
      cat("\nColumns to select:\n")
      cat("  Base:", paste(base_cols, collapse = ", "), "\n")
      cat("  Traits:", paste(trait_cols, collapse = ", "), "\n")
      
      if (length(trait_cols) < 3) {
        cat("ERROR: Not enough trait columns found\n")
        cat("All available columns:\n")
        print(cols)
        stop("Insufficient trait columns for analysis")
      }
      
      dplyr::select(., all_of(c(base_cols, trait_cols)))
    } %>%
    # Rename to standard names
    {
      cat("\nColumns before renaming:", paste(colnames(.), collapse = ", "), "\n")
      
      rename_with(., ~ case_when(
        grepl("habitat.*forest|forest.*open|forest.*area", .x, ignore.case = TRUE) ~ "habita_for",
        grepl("diet", .x, ignore.case = TRUE) ~ "diet",
        grepl("^PC1$|life.*history", .x, ignore.case = TRUE) ~ "PC1",
        grepl("body.*mass|mass.*log", .x, ignore.case = TRUE) ~ "body_mass",
        grepl("^HWI$", .x, ignore.case = TRUE) ~ "HWI", 
        grepl("migration.*distance", .x, ignore.case = TRUE) ~ "distance_mig",
        grepl("^Latitude$|^latitude$", .x, ignore.case = TRUE) ~ "Latitude",
        TRUE ~ .x
      ))
    } %>%
    {
      cat("Columns after renaming:", paste(colnames(.), collapse = ", "), "\n")
      .
    } %>%
    # Create log_HWI if HWI exists
    {
      if ("HWI" %in% colnames(.)) {
        mutate(., log_HWI = log(HWI + 1))
      } else {
        mutate(., log_HWI = NA)
      }
    } %>%
    # Convert to numeric and scale available continuous variables
    {
      continuous_vars <- intersect(c("body_mass", "log_HWI", "PC1", "distance_mig", "Latitude"), 
                                   colnames(.))
      cat("Scaling continuous variables:", paste(continuous_vars, collapse = ", "), "\n")
      
      if (length(continuous_vars) > 0) {
        mutate(., across(all_of(continuous_vars), ~ as.numeric(scale(.x))))
      } else {
        .
      }
    } %>%
    # Convert categorical variables to numeric
    {
      categorical_vars <- intersect(c("diet", "habita_for"), colnames(.))
      cat("Converting categorical variables:", paste(categorical_vars, collapse = ", "), "\n")
      if (length(categorical_vars) > 0) {
        mutate(., across(all_of(categorical_vars), ~ as.numeric(as.character(.x))))
      } else {
        .
      }
    } %>%
    # Remove rows with missing response variables
    filter(!is.na(Weibull_median_log), !is.na(Weibull_upper_distance_log)) %>%
    as.data.frame()
  
  cat("=== END DEBUGGING ===\n")
  
  # Set rownames and match with phylogeny
  rownames(analysis_data) <- analysis_data$label
  
  # Final species matching
  final_species <- intersect(analysis_data$label, dispersal_tree$tip.label)
  analysis_data <- analysis_data[analysis_data$label %in% final_species, ]
  dispersal_tree_final <- keep.tip(dispersal_tree, final_species)
  
  # Create phylogenetic covariance matrix
  A <- vcv.phylo(dispersal_tree_final)
  
  # Check which variables actually made it to the final dataset
  cat("Final dataset columns:", paste(colnames(analysis_data), collapse = ", "), "\n")
  
  # Update base_predictors based on what's actually available
  original_predictors <- c("log_HWI", "body_mass", "habita_for", "PC1", "diet", "distance_mig", "Latitude")
  available_predictors <- intersect(original_predictors, colnames(analysis_data))
  
  if (length(available_predictors) < length(original_predictors)) {
    cat("WARNING: Some predictors not available in final dataset\n")
    cat("Expected:", paste(original_predictors, collapse = ", "), "\n")
    cat("Available:", paste(available_predictors, collapse = ", "), "\n")
    
    # Remove missing predictors
    missing_predictors <- setdiff(original_predictors, available_predictors)
    cat("Missing:", paste(missing_predictors, collapse = ", "), "\n")
  }
  
  # Update the base_predictors list globally
  base_predictors <<- available_predictors
  
  # Update priority_interactions based on available predictors
  original_interactions <- c(
    "body_mass * PC1",
    "body_mass * diet", 
    "log_HWI * distance_mig",
    "diet * habita_for",
    "PC1 * diet",
    "Latitude * distance_mig"
  )
  
  valid_interactions <- c()
  for (interaction in original_interactions) {
    terms <- unlist(strsplit(interaction, " \\* "))
    if (all(terms %in% available_predictors)) {
      valid_interactions <- c(valid_interactions, interaction)
    } else {
      cat("Skipping interaction", interaction, "- missing terms:", 
          paste(setdiff(terms, available_predictors), collapse = ", "), "\n")
    }
  }
  
  # Update the global priority_interactions
  priority_interactions <<- valid_interactions
  
  cat("✓ Final dataset:", nrow(analysis_data), "species\n")
  cat("✓ Available predictors:", paste(available_predictors, collapse = ", "), "\n")
  cat("✓ Valid interactions:", length(valid_interactions), "\n")
  if (length(valid_interactions) > 0) {
    cat("  -", paste(valid_interactions, collapse = "\n  - "), "\n")
  }
  cat("\n")
  
  return(list(
    data = analysis_data,
    A = A,
    tree = dispersal_tree_final
  ))
}

# MODEL FITTING FUNCTION =====================================================

fit_interaction_models <- function(data_results, dispersal_type, verbose = TRUE) {
  
  if (verbose) cat("=== FITTING INTERACTION MODELS FOR", toupper(dispersal_type), "===\n")
  
  data <- data_results$data
  A <- data_results$A
  
  # Storage for results
  results <- list(
    dispersal_type = dispersal_type,
    models = list(),
    comparisons = list(),
    data_info = list(
      n_species = nrow(data),
      base_predictors = base_predictors,  # Use the updated predictors
      interactions_tested = priority_interactions  # Use the updated interactions
    )
  )
  
  # BASE MODELS =============================================================
  
  if (verbose) cat("\n--- FITTING BASE MODELS ---\n")
  
  # Base formula
  base_formula_median <- paste0("Weibull_median_log ~ ", paste(base_predictors, collapse = " + "), " + (1|gr(label, cov = A))")
  base_formula_long <- paste0("Weibull_upper_distance_log ~ ", paste(base_predictors, collapse = " + "), " + (1|gr(label, cov = A))")
  
  if (verbose) cat("Fitting base model for median dispersal...\n")
  base_model_median <- brm(
    base_formula_median,
    data = data,
    family = gaussian(),
    data2 = list(A = A),
    chains = model_settings$chains,
    cores = model_settings$cores,
    iter = model_settings$iter,
    control = model_settings$control,
    save_pars = model_settings$save_pars,
    refresh = 0
  )
  base_model_median <- add_criterion(base_model_median, "loo", moment_match = TRUE, reloo = TRUE)
  
  if (verbose) cat("Fitting base model for long-distance dispersal...\n")
  base_model_long <- brm(
    base_formula_long,
    data = data,
    family = gaussian(),
    data2 = list(A = A),
    chains = model_settings$chains,
    cores = model_settings$cores,
    iter = model_settings$iter,
    control = model_settings$control,
    save_pars = model_settings$save_pars,
    refresh = 0
  )
  base_model_long <- add_criterion(base_model_long, "loo", moment_match = TRUE, reloo = TRUE)
  
  # Store base models
  results$models$base_median <- base_model_median
  results$models$base_long <- base_model_long
  
  if (verbose) cat("✓ Base models fitted successfully\n")
  
  # INTERACTION MODELS ======================================================
  
  if (verbose) cat("\n--- FITTING INTERACTION MODELS ---\n")
  
  for (i in seq_along(priority_interactions)) {
    
    interaction <- priority_interactions[i]
    interaction_name <- gsub(" \\* ", "_x_", interaction)
    
    # Check if interaction terms exist in data
    interaction_terms <- unlist(strsplit(interaction, " \\* "))
    if (!all(interaction_terms %in% base_predictors)) {
      if (verbose) cat("Skipping", interaction, "- terms not available\n")
      next
    }
    
    if (verbose) cat("Fitting models with interaction:", interaction, "\n")
    
    # Median dispersal with interaction
    tryCatch({
      formula_median_int <- paste0(base_formula_median, " + ", interaction)
      
      model_median_int <- brm(
        formula_median_int,
        data = data,
        family = gaussian(),
        data2 = list(A = A),
        chains = model_settings$chains,
        cores = model_settings$cores,
        iter = model_settings$iter,
        control = model_settings$control,
        save_pars = model_settings$save_pars,
        refresh = 0
      )
      model_median_int <- add_criterion(model_median_int, "loo", moment_match = TRUE, reloo = TRUE)
      
      # Store model
      results$models[[paste0("median_", interaction_name)]] <- model_median_int
      
      # Compare with base model - with error handling
      tryCatch({
        comparison_median <- loo_compare(loo(base_model_median), loo(model_median_int))
        results$comparisons[[paste0("median_", interaction_name)]] <- list(
          interaction = interaction,
          delta_looic = comparison_median[2, "elpd_diff"],
          se_diff = comparison_median[2, "se_diff"],
          comparison = comparison_median
        )
      }, error = function(e2) {
        if (verbose) cat("    Warning: LOO comparison failed for", interaction, "- using WAIC instead\n")
        # Fallback to WAIC
        model_median_int <- add_criterion(model_median_int, "waic")
        base_model_median_waic <- add_criterion(base_model_median, "waic")
        
        comparison_median <- loo_compare(waic(base_model_median_waic), waic(model_median_int))
        results$comparisons[[paste0("median_", interaction_name)]] <- list(
          interaction = interaction,
          delta_looic = comparison_median[2, "elpd_diff"],
          se_diff = comparison_median[2, "se_diff"],
          comparison = comparison_median,
          method = "waic"
        )
      })
      
      if (verbose) cat("  ✓ Median model with", interaction, "completed\n")
      
    }, error = function(e) {
      if (verbose) cat("  ✗ Error with median model for", interaction, ":", e$message, "\n")
    })
    
    # Long-distance dispersal with interaction
    tryCatch({
      formula_long_int <- paste0(base_formula_long, " + ", interaction)
      
      model_long_int <- brm(
        formula_long_int,
        data = data,
        family = gaussian(),
        data2 = list(A = A),
        chains = model_settings$chains,
        cores = model_settings$cores,
        iter = model_settings$iter,
        control = model_settings$control,
        save_pars = model_settings$save_pars,
        refresh = 0
      )
      model_long_int <- add_criterion(model_long_int, "loo", moment_match = TRUE, reloo = TRUE)
      
      # Store model
      results$models[[paste0("long_", interaction_name)]] <- model_long_int
      
      # Compare with base model - with error handling
      tryCatch({
        comparison_long <- loo_compare(loo(base_model_long), loo(model_long_int))
        results$comparisons[[paste0("long_", interaction_name)]] <- list(
          interaction = interaction,
          delta_looic = comparison_long[2, "elpd_diff"], 
          se_diff = comparison_long[2, "se_diff"],
          comparison = comparison_long
        )
      }, error = function(e2) {
        if (verbose) cat("    Warning: LOO comparison failed for", interaction, "- using WAIC instead\n")
        # Fallback to WAIC
        model_long_int <- add_criterion(model_long_int, "waic")
        base_model_long_waic <- add_criterion(base_model_long, "waic")
        
        comparison_long <- loo_compare(waic(base_model_long_waic), waic(model_long_int))
        results$comparisons[[paste0("long_", interaction_name)]] <- list(
          interaction = interaction,
          delta_looic = comparison_long[2, "elpd_diff"],
          se_diff = comparison_long[2, "se_diff"], 
          comparison = comparison_long,
          method = "waic"
        )
      })
      
      if (verbose) cat("  ✓ Long-distance model with", interaction, "completed\n")
      
    }, error = function(e) {
      if (verbose) cat("  ✗ Error with long-distance model for", interaction, ":", e$message, "\n")
    })
  }
  
  if (verbose) cat("\n✓ All interaction models fitted\n")
  
  return(results)
}

# RESULTS SUMMARY FUNCTION ===================================================

summarize_results <- function(results, verbose = TRUE) {
  
  if (verbose) cat("\n=== RESULTS SUMMARY FOR", toupper(results$dispersal_type), "===\n")
  
  # Create summary table
  summary_table <- data.frame(
    interaction = character(),
    response = character(),
    delta_looic = numeric(),
    se_diff = numeric(),
    interpretation = character(),
    stringsAsFactors = FALSE
  )
  
  # Process comparisons
  for (comp_name in names(results$comparisons)) {
    comp <- results$comparisons[[comp_name]]
    
    response_type <- ifelse(grepl("median", comp_name), "median", "long_distance")
    
    interpretation <- case_when(
      comp$delta_looic > 2 * comp$se_diff ~ "Strong improvement",
      comp$delta_looic > comp$se_diff ~ "Moderate improvement",
      abs(comp$delta_looic) <= comp$se_diff ~ "No difference",
      comp$delta_looic < -comp$se_diff ~ "Worse performance",
      TRUE ~ "Inconclusive"
    )
    
    summary_table <- rbind(summary_table, data.frame(
      interaction = comp$interaction,
      response = response_type,
      delta_looic = round(comp$delta_looic, 2),
      se_diff = round(comp$se_diff, 2),
      interpretation = interpretation,
      stringsAsFactors = FALSE
    ))
  }
  
  # Print summary
  if (verbose) {
    cat("INTERACTION EFFECTS SUMMARY:\n")
    print(summary_table)
    
    # Highlight significant interactions
    significant <- summary_table[abs(summary_table$delta_looic) > 2, ]
    if (nrow(significant) > 0) {
      cat("\nSIGNIFICANT INTERACTIONS (|Δ LOOIC| > 2):\n")
      for (i in 1:nrow(significant)) {
        cat("•", significant$interaction[i], "on", significant$response[i], 
            "dispersal: Δ =", significant$delta_looic[i], "\n")
      }
    } else {
      cat("\nNo interactions showed strong evidence (|Δ LOOIC| > 2)\n")
    }
  }
  
  return(summary_table)
}

# SAVE RESULTS FUNCTION ======================================================

save_results <- function(results, summary_table, dispersal_type) {
  
  # Create results directory
  results_dir <- paste0("results/interactions_simple/", dispersal_type)
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Save models
  save(results, file = file.path(results_dir, paste0("interaction_models_", dispersal_type, ".RData")))
  
  # Save summary table
  write.csv(summary_table, file.path(results_dir, paste0("interaction_summary_", dispersal_type, ".csv")), 
            row.names = FALSE)
  
  cat("Results saved in:", results_dir, "\n")
}

# MAIN ANALYSIS FUNCTION =====================================================

run_interaction_analysis <- function(dispersal_type = "average", save_results_flag = TRUE) {
  
  cat("=== STARTING INTERACTION ANALYSIS ===\n")
  cat("Dispersal type:", dispersal_type, "\n")
  cat("Interactions to test:", length(priority_interactions), "\n")
  cat("Expected models: ~", length(priority_interactions) * 2 + 2, "(", length(priority_interactions), "interactions × 2 responses + 2 base models )\n\n")
  
  # Load data
  data_results <- load_dispersal_data(dispersal_type)
  
  # Fit models
  results <- fit_interaction_models(data_results, dispersal_type)
  
  # Summarize results
  summary_table <- summarize_results(results)
  
  # Save results
  if (save_results_flag) {
    save_results(results, summary_table, dispersal_type)
  }
  
  cat("\n=== ANALYSIS COMPLETED ===\n")
  
  return(list(
    results = results,
    summary = summary_table
  ))
}

# RUN ANALYSIS FOR ALL DISPERSAL TYPES ======================================

cat("=== INTERACTION ANALYSIS PIPELINE ===\n")
cat("This script fits models with 6 priority interactions\n")
cat("For each dispersal type: ~14 models (2 base + 6×2 interactions)\n\n")

# Store all results
all_results <- list()

# Run for average dispersal
cat("1/3: AVERAGE DISPERSAL\n")
all_results$average <- run_interaction_analysis("average")

cat("\n" + paste(rep("=", 60), collapse = "") + "\n")

# Run for breeding dispersal
cat("2/3: BREEDING DISPERSAL\n") 
all_results$breeding <- run_interaction_analysis("breeding")

cat("\n" + paste(rep("=", 60), collapse = "") + "\n")

# Run for natal dispersal
cat("3/3: NATAL DISPERSAL\n")
all_results$natal <- run_interaction_analysis("natal")

# FINAL SUMMARY ==============================================================

cat("\n" + paste(rep("=", 80), collapse = "") + "\n")
cat("=== COMPREHENSIVE SUMMARY ===\n")

# Combine all summaries
combined_summary <- do.call(rbind, lapply(names(all_results), function(x) {
  summary <- all_results[[x]]$summary
  summary$dispersal_type <- x
  return(summary)
}))

# Save combined results
write.csv(combined_summary, "results/interactions_simple/combined_interaction_summary.csv", row.names = FALSE)

# Print overall summary
cat("INTERACTIONS ACROSS ALL DISPERSAL TYPES:\n")
print(combined_summary)

# Count significant interactions
significant_total <- combined_summary[abs(combined_summary$delta_looic) > 2, ]
cat("\nOVERALL SIGNIFICANT INTERACTIONS:", nrow(significant_total), "/", nrow(combined_summary), "\n")

if (nrow(significant_total) > 0) {
  cat("Most consistent interactions:\n")
  interaction_counts <- table(significant_total$interaction)
  for (interaction in names(sort(interaction_counts, decreasing = TRUE))) {
    cat("•", interaction, ":", interaction_counts[interaction], "contexts\n")
  }
}

cat("\nTotal models fitted:", nrow(combined_summary) + length(all_results) * 2, "\n")
cat("Results saved in: results/interactions_simple/\n")

cat("\n=== INTERACTION ANALYSIS COMPLETE ===\n")