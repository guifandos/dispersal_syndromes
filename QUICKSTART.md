# Quick Start Guide

Instructions to get started with the dispersal analysis code.

## 1. Clone repository

```bash
git clone https://github.com/username/dispersal-syndromes-birds.git
cd dispersal-syndromes-birds
```

## 2. Install dependencies

### Option A: Automatic installation (recommended)
```r
# Open R/RStudio in the project directory
source("revision3/function_dispersal_R3R.R")
# Packages will be installed automatically
```

### Option B: Manual installation
```r
install.packages(c(
  "brms", "projpred", "ape", "phytools", "geiger",
  "tidyverse", "ggplot2", "DHARMa", "performance",
  "mice", "caper", "phylolm", "picante"
))
```

### Option C: Use renv (exact reproducibility)
```r
install.packages("renv")
renv::restore()  # Installs exact package versions
```

## 3. Download data

Data are NOT in the GitHub repository due to size.

**Option 1**: Download from Zenodo/Dryad
```bash
# Download and extract to data/ folder
wget [ZENODO_URL]
unzip data.zip -d data/
```

**Option 2**: Example data (testing only)
```r
# Example data are in data/example/
# Allows code testing but with fewer species
```

## 4. Verify folder structure

```bash
# Create necessary folders
mkdir -p data/{raw,processed}
mkdir -p results/{weibull,combined}
mkdir -p figures
```

Or from R:
```r
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/weibull", recursive = TRUE, showWarnings = FALSE)
dir.create("results/combined", recursive = TRUE, showWarnings = FALSE)
dir.create("figures", recursive = TRUE, showWarnings = FALSE)
```

## 5. Run basic analysis

### Prepare data
```r
source("revision3/prepare_disperal_age_R3.R")

# Verify generated files
list.files("data/processed")
# Should show: dispersal_average_complete.RData, etc.
```

### Run analysis for one dispersal type
```r
source("revision3/function_dispersal_R3R.R")

# Complete analysis for average dispersal
results_average <- run_dispersal_analysis(
  dispersal_type = "average",
  save_results = TRUE,
  verbose = TRUE
)
```

⚠️ **Note**: Bayesian models may take 30-60 minutes per dispersal type.

## 6. Generate main figures

```r
# Load and combine models
source("revision3/load_models_R3.R")

# Forest plots
source("revision3/complete_forest_plot_R3.R")

# Predictions
source("revision3/5c_prediction_plots_final_R3.R")

# Figures saved in figures/
```

## 7. Verify outputs

```r
# View model summaries
load("results/combined/short_models_combined.RData")
load("results/combined/complete_models_combined.RData")

# View structure
str(all_short_models)
str(all_complete_models)

# View generated figures
list.files("figures", recursive = TRUE)
```

---

## Complete workflow

```r
# 1. Prepare data
source("revision3/prepare_disperal_age_R3.R")

# 2. Run analyses for all 3 types
source("revision3/function_dispersal_R3R.R")
run_dispersal_analysis(dispersal_type = "average")
run_dispersal_analysis(dispersal_type = "natal")
run_dispersal_analysis(dispersal_type = "breeding")

# 3. Combine results
source("revision3/load_models_R3.R")

# 4. Compare models
source("revision3/comprehensive_model_comparison_R3.R")

# 5. Generate all figures
source("revision3/complete_forest_plot_R3.R")
source("revision3/5c_conditional_effect_plots_review3_only_average.R")
source("revision3/5c_prediction_plots_final_R3.R")

# 6. Generate tables
source("revision3/tabla_comprenhesive_results_R3.R")
source("revision3/univariate_multivariate_complete_table_plots_R3.R")
```

---

## Common troubleshooting

### Error: "cannot open file"
```r
# Check working directory
getwd()
setwd("path/to/dispersal-syndromes-birds")
```

### Error: "package not found"
```r
# Reinstall packages
source("revision3/function_dispersal_R3R.R")
```

### Error: "phylogenetic tree and data don't match"
```r
# Verify data preparation
source("revision3/prepare_disperal_age_R3.R")
```

### Models too slow
```r
# Reduce iterations for testing (NOT for final analysis)
# In function_dispersal_R3R.R, change:
# iter = 1000 (instead of 2000)
# chains = 2 (instead of 4)
```

### Insufficient memory
```r
# Increase available memory
# In .Renviron:
R_MAX_VSIZE = 100Gb
```

---

## Additional resources

- **Complete documentation**: See `docs/`
- **Technical details**: See `docs/TECHNICAL_SUMMARY.md`
- **Version changes**: See `CHANGELOG.md`
- **Repository structure**: See `docs/REPOSITORY_STRUCTURE.md`

---

## Next steps

After running basic analysis:

1. Explore individual models in `results/weibull/`
2. Review variable selection in `results/combined/`
3. Customize figures by modifying scripts in `revision3/`
4. Compare with your own data

---

## Support

If you encounter problems:

1. Check documentation in `docs/`
2. Verify session info: `sessionInfo()`
3. Contact author: gfandos@ucm.es

---

**Ready to start!** 🚀
