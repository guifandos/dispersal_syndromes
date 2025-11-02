# QUICK MODEL COMPARISON ====================================================
# Simple, direct comparison without functions

library(brms)
library(loo)
library(dplyr)

# LOAD MODELS =================================================================

cat("=== LOADING MODELS ===\n")

# Univariate models
load("results/weibull/average/univariate_models.RData")
univariate_average <- models

load("results/weibull/natal/univariate_models.RData") 
univariate_natal <- models

load("results/weibull/breeding/univariate_models.RData")
univariate_breeding <- models

# Multivariate models
load("results/combined/complete_models_combined.RData")
multivariate_complete <- all_complete_models

load("results/combined/short_models_combined.RData")
multivariate_reduced <- all_short_models

cat("✓ All models loaded\n")
cat("Univariate average models:", length(univariate_average), "\n")
cat("Univariate natal models:", length(univariate_natal), "\n") 
cat("Univariate breeding models:", length(univariate_breeding), "\n")
cat("Complete multivariate models:", nrow(multivariate_complete), "\n")
cat("Reduced multivariate models:", nrow(multivariate_reduced), "\n\n")

# COMPARISON FUNCTION =========================================================

compare_models_simple <- function(model1, model2, name1, name2) {
  
  cat("Comparing", name1, "vs", name2, "...\n")
  
  # Add LOO if needed
  if (!"loo" %in% names(model1$criteria)) {
    model1 <- add_criterion(model1, "loo", moment_match = TRUE, reloo = TRUE)
  }
  if (!"loo" %in% names(model2$criteria)) {
    model2 <- add_criterion(model2, "loo", moment_match = TRUE, reloo = TRUE)
  }
  
  # Get LOOIC
  looic1 <- loo(model1)$estimates["looic", "Estimate"]
  looic2 <- loo(model2)$estimates["looic", "Estimate"]
  
  # Compare
  delta <- looic1 - looic2  # Positive = model1 worse, negative = model1 better
  
  cat("  ", name1, "LOOIC:", round(looic1, 2), "\n")
  cat("  ", name2, "LOOIC:", round(looic2, 2), "\n")
  cat("  Delta (", name1, "-", name2, "):", round(delta, 2), "\n")
  
  if (delta < -2) {
    winner <- paste(name1, "much better")
  } else if (delta < 0) {
    winner <- paste(name1, "better")
  } else if (delta > 2) {
    winner <- paste(name2, "much better")
  } else if (delta > 0) {
    winner <- paste(name2, "better")
  } else {
    winner <- "Similar"
  }
  
  cat("  Winner:", winner, "\n\n")
  
  return(data.frame(
    model1 = name1,
    model2 = name2,
    looic1 = round(looic1, 2),
    looic2 = round(looic2, 2),
    delta = round(delta, 2),
    winner = winner,
    stringsAsFactors = FALSE
  ))
}

# TEST AVERAGE MEDIAN =========================================================

cat("=== AVERAGE MEDIAN DISPERSAL ===\n")

# Get models
body_mass_model <- univariate_average$body_mass_median
pc1_model <- univariate_average$PC1_median

complete_model <- multivariate_complete %>%
  filter(age == "average", dispersal_mode == "median") %>%
  pull(model) %>% .[[1]]

reduced_model <- multivariate_reduced %>%
  filter(age == "average", dispersal_mode == "median") %>%
  pull(model) %>% .[[1]]

# Comparisons
results <- list()

results[[1]] <- compare_models_simple(body_mass_model, complete_model, "body_mass", "complete")
results[[2]] <- compare_models_simple(pc1_model, complete_model, "PC1", "complete")
results[[3]] <- compare_models_simple(reduced_model, complete_model, "reduced", "complete")
results[[4]] <- compare_models_simple(body_mass_model, reduced_model, "body_mass", "reduced")

# TEST BREEDING MEDIAN ========================================================

cat("=== BREEDING MEDIAN DISPERSAL ===\n")

# Get models
body_mass_breed <- univariate_breeding$body_mass_median
pc1_breed <- univariate_breeding$PC1_median

complete_breed <- multivariate_complete %>%
  filter(age == "breeding", dispersal_mode == "median") %>%
  pull(model) %>% .[[1]]

reduced_breed <- multivariate_reduced %>%
  filter(age == "breeding", dispersal_mode == "median") %>%
  pull(model) %>% .[[1]]

# Comparisons
results[[5]] <- compare_models_simple(body_mass_breed, complete_breed, "body_mass_breed", "complete_breed")
results[[6]] <- compare_models_simple(pc1_breed, complete_breed, "PC1_breed", "complete_breed")
results[[7]] <- compare_models_simple(reduced_breed, complete_breed, "reduced_breed", "complete_breed")

# TEST NATAL MEDIAN ===========================================================

cat("=== NATAL MEDIAN DISPERSAL ===\n")

# Get models
body_mass_natal <- univariate_natal$body_mass_median
pc1_natal <- univariate_natal$PC1_median

complete_natal <- multivariate_complete %>%
  filter(age == "natal", dispersal_mode == "median") %>%
  pull(model) %>% .[[1]]

reduced_natal <- multivariate_reduced %>%
  filter(age == "natal", dispersal_mode == "median") %>%
  pull(model) %>% .[[1]]

# Comparisons
results[[8]] <- compare_models_simple(body_mass_natal, complete_natal, "body_mass_natal", "complete_natal")
results[[9]] <- compare_models_simple(pc1_natal, complete_natal, "PC1_natal", "complete_natal")
results[[10]] <- compare_models_simple(reduced_natal, complete_natal, "reduced_natal", "complete_natal")

# COMBINE RESULTS =============================================================

cat("=== SUMMARY TABLE ===\n")
all_results <- do.call(rbind, results)
print(all_results)

# Save results
write.csv(all_results, "quick_comparison_results.csv", row.names = FALSE)

cat("\n=== KEY FINDINGS ===\n")
univariate_wins <- sum(grepl("body_mass|PC1", all_results$winner))
multivariate_wins <- sum(grepl("complete|reduced", all_results$winner))

cat("Univariate wins:", univariate_wins, "/", nrow(all_results), "\n")
cat("Multivariate wins:", multivariate_wins, "/", nrow(all_results), "\n")

cat("\nResults saved in: quick_comparison_results.csv\n")