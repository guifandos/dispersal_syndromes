# Revision 3 - Dispersal Syndromes Analysis

Code for the third revision of the manuscript **"Complex trait syndromes help understand but do not predict dispersal distances in European birds"**.

## Script organization

### 1. Data preparation
- **`prepare_disperal_age_R3.R`**: Data cleaning, preparation of datasets and phylogenetic matrices for three dispersal types (average, natal, breeding).

### 2. Main functions
- **`function_dispersal_R3R.R`**: Complete modular workflow with functions to run dispersal analyses including Bayesian model fitting with `brms` and variable selection with `projpred`.

### 3. Model loading and combination
- **`load_models_R3.R`**: Loads and combines reduced and complete models for the three dispersal types. Extracts and analyzes variable selection results.

### 4. Model comparison
- **`clean_model_comparison_R3.R`**: Comparison of predictive performance across models.
- **`comprehensive_model_comparison_R3.R`**: Comprehensive cross-validation analysis within and between orders.
- **`comparison_univariate_multivariate_complete_models_R3.R`**: Comparison between univariate and multivariate models.

### 5. Visualizations
- **`forest_plot_average_R3_new.R`**: Forest plots for average dispersal models.
- **`complete_forest_plot_R3.R`**: Complete forest plots for all dispersal types.
- **`5c_conditional_effect_plots_review3_only_average.R`**: Conditional effects for average dispersal.
- **`5c_prediction_plots_final_R3.R`**: Final prediction plots.

### 6. Results tables
- **`tabla_comprenhesive_results_R3.R`**: Comprehensive table of all results.
- **`univariate_multivariate_complete_table_plots_R3.R`**: Tables and plots comparing univariate vs multivariate models.

## Recommended workflow

1. `prepare_disperal_age_R3.R` - Prepare data
2. `function_dispersal_R3R.R` - Run main analyses
3. `load_models_R3.R` - Load and combine models
4. Comparison and visualization scripts as needed

## Main dependencies

- `brms`: Bayesian models
- `projpred`: Variable selection
- `ape`, `phytools`: Phylogenetic analyses
- `tidyverse`: Data manipulation
- `ggplot2`: Visualizations

## Notes

- Scripts assume directory structure: `data/`, `results/`, `figures/`
- All analyses use Weibull distribution for dispersal kernels
- Phylogenetic relationships controlled in all models

## Author

Guillermo Fandos  
Revision 3 - September 2025
