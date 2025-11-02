# Changelog

Change history for the European bird dispersal syndromes project.

## [Revision 3] - 2025-09

### Added
- Complete modular workflow in `function_dispersal_R3R.R`
- Separate analysis for average, natal and breeding dispersal
- Comprehensive within-order and between-order cross-validation
- Specific scripts for univariate vs multivariate comparison
- Comprehensive documentation in all scripts
- Automatic data integrity verification

### Changed
- More robust data preparation with error handling
- Standardized variable nomenclature
- Improved forest plots with importance ranking
- Optimized prediction plots
- Enhanced phylogenetic signal (λ) calculation

### Analysis
- Emphasis on Weibull distribution for dispersal kernels
- Bayesian models with `brms` and selection with `projpred`
- Comprehensive comparison between single-trait and dispersal syndrome models
- Optimized conditional effects analysis

### Main scripts
1. `prepare_disperal_age_R3.R` - Data preparation and phylogenetic matrices
2. `function_dispersal_R3R.R` - Main modular workflow
3. `load_models_R3.R` - Model combination and analysis
4. `*_model_comparison_R3.R` - Model comparison
5. `*_plot_R3*.R` - Visualizations
6. `*_table_R3.R` - Results tables

### Key results
- Single-trait models (body mass, life history) outperform dispersal syndromes in between-order prediction
- HWI (flight efficiency) important for natal and breeding dispersal
- Significant phylogenetic signal in all models
- Different drivers for median vs long-distance dispersal

---

## [Revision 2] - [Previous date]

### Changes
- [List main changes from R2]

---

## [Revision 1] - [Previous date]

### First version
- Initial dispersal syndromes analysis
- Basic models with morphological and life-history traits
- First phylogenetic comparisons

---

## Versioning format

This project follows [Semantic Versioning](https://semver.org/):
- **Major**: Large changes in methodology or structure
- **Minor**: Addition of new analyses or features
- **Patch**: Bug fixes or minor improvements

## Contributions

Contributions should be documented in this file following the format:
- **Added**: New features
- **Changed**: Changes to existing features
- **Removed**: Removed features
- **Fixed**: Bug fixes
- **Security**: Vulnerability fixes

---

**Maintainer**: Guillermo Fandos  
**Last update**: November 2025
