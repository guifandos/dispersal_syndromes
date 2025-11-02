# Technical Summary of Scripts - Revision 3

## Main changes from previous versions

### General improvements in R3
- Implementation of modular and reproducible workflow
- Comprehensive documentation of functions and parameters
- Robust error handling and data verification
- Separate analysis for three dispersal types: average, natal, breeding
- Within-order and between-order cross-validation
- Standardized variable nomenclature

---

## Preparation scripts

### `prepare_disperal_age_R3.R`
**Purpose**: Preparation and cleaning of phylogenetic data

**Main functions**:
- `prepare_dispersal_data()`: Filtering by dispersal type, log transformations, phylogenetic tree adjustments
- `verify_dataset()`: Data integrity verification

**Outputs**:
- `dispersal_average_complete.RData`
- `dispersal_natal_complete.RData`  
- `dispersal_breeding_complete.RData`

**Key variables**:
- Morphological: `body_mass`, `log_HWI`
- Life-history: `PC1`, `PC2`
- Ecological: `habita_for`, `diet`
- Geographic: `distance_mig`, `Latitude`
- Response: `Weibull_median_log`, `Weibull_upper_distance_log`

---

## Analysis scripts

### `function_dispersal_R3R.R`
**Purpose**: Complete modularized workflow

**Main function**:
```r
run_dispersal_analysis(dispersal_type = "average", 
                       save_results = TRUE,
                       verbose = TRUE)
```

**Features**:
- Automatic package installation
- Bayesian models with `brms`
- Variable selection with `projpred`
- Automatic diagnostics (DHARMa)
- Phylogenetic signal calculation (λ)

---

### `load_models_R3.R`
**Purpose**: Model combination and analysis

**Functionalities**:
- Loading reduced and complete models for all 3 types
- Combination into unified dataframes
- Variable selection extraction
- Variable cleaning and renaming
- MLPD scaling by age
- Summary tables and plots generation

**Outputs**:
- `short_models_combined.RData`
- `complete_models_combined.RData`
- Variable selection CSVs (raw and scaled)
- Variable selection plots

---

## Comparison scripts

### `clean_model_comparison_R3.R`
**Purpose**: Predictive performance comparison

**Metrics**:
- R² 
- RMSE
- Within and between orders comparison

---

### `comprehensive_model_comparison_R3.R`
**Purpose**: Comprehensive cross-validation analysis

**Methodology**:
- K-fold cross-validation within orders
- Leave-one-order-out between orders
- Univariate vs multivariate (syndrome) model comparison

---

### `comparison_univariate_multivariate_complete_models_R3.R`
**Purpose**: Specific univariate vs syndrome comparison

**Analysis**:
- Predictive performance by individual trait
- Comparison with multi-trait model (syndrome)
- Identification of most informative traits

---

## Visualization scripts

### `complete_forest_plot_R3.R`
**Purpose**: Forest plots for all dispersal types

**Elements**:
- Standardized coefficients
- 95% credible intervals
- Variable importance
- Separation by median/long distance

---

### `forest_plot_average_R3_new.R`
**Purpose**: Specific forest plot for average dispersal

**Improvements**:
- Optimized visualization
- Variable importance ranking

---

### `5c_conditional_effect_plots_review3_only_average.R`
**Purpose**: Conditional effects of traits on average dispersal

**Visualizations**:
- Effects of HWI, body mass, life history
- Confidence intervals
- Univariate models controlling for phylogeny

---

### `5c_prediction_plots_final_R3.R`
**Purpose**: Within/between orders prediction plots

**Content**:
- Observed vs predicted scatter plots
- R² by model and order
- Visual performance comparison

---

## Table scripts

### `tabla_comprenhesive_results_R3.R`
**Purpose**: Comprehensive results table

**Content**:
- Coefficients from all models
- Fit statistics
- Phylogenetic signal
- Selected variables

---

### `univariate_multivariate_complete_table_plots_R3.R`
**Purpose**: Univariate/multivariate comparison tables and plots

**Outputs**:
- R² table by model
- Comparative plots
- Trait ranking

---

## Expected data structure

```
project/
├── data/
│   ├── data_philo/
│   │   └── dispersal_tree500.nex
│   ├── data_process/
│   │   └── dispersal_traits_total_PCA_20220424.csv
│   ├── dispersal_distance/
│   │   └── Table_S13_species_dispersal_distances.csv
│   └── processed/
│       ├── dispersal_average_complete.RData
│       ├── dispersal_natal_complete.RData
│       └── dispersal_breeding_complete.RData
├── results/
│   ├── weibull/
│   │   ├── average/
│   │   ├── natal/
│   │   └── breeding/
│   └── combined/
└── figures/
```

---

## Naming conventions

### Variables
- `log_*`: Log-transformed variable
- `*_log`: Transformed response variable
- `PC1`, `PC2`: Life-history principal components

### Dispersal types
- `average`: Total dispersal (breeding + natal)
- `natal`: Natal dispersal
- `breeding`: Breeding dispersal

### Distance types
- `median`: Median of Weibull kernel
- `long`: 95th percentile of Weibull kernel

### Models
- `short`: Reduced models after variable selection
- `complete`: Complete models with all variables

---

## Key methodological notes

1. **Distribution**: Weibull for dispersal kernels (best fit for long-distance events)
2. **Phylogeny**: Always controlled via phylogenetic covariance matrix
3. **Transformations**: Continuous variables standardized (mean=0, sd=1)
4. **Variable selection**: `projpred` with LOO cross-validation
5. **Diagnostics**: DHARMa for residuals, VIF for multicollinearity
6. **Validation**: Within-order (5-fold CV), Between-order (leave-one-out)

---

**Date**: September 2025  
**Author**: Guillermo Fandos
