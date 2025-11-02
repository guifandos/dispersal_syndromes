# Recommended Repository Structure

```
dispersal-syndromes-birds/
│
├── README.md                          # General project description
├── .gitignore                         # Files to ignore in git
├── LICENSE                            # Code license
│
├── revision3/                         # Revision 3 scripts
│   ├── README.md                      # R3-specific documentation
│   ├── 01_prepare_disperal_age_R3.R
│   ├── 02_function_dispersal_R3R.R
│   ├── 03_load_models_R3.R
│   ├── 04_model_comparison_*.R
│   └── 05_visualization_*.R
│
├── data/                              # Data (don't upload large files to GitHub)
│   ├── raw/                           # Original unmodified data
│   │   ├── dispersal_tree500.nex
│   │   ├── dispersal_traits_total_PCA.csv
│   │   └── species_dispersal_distances.csv
│   ├── processed/                     # Processed data
│   │   ├── dispersal_average_complete.RData
│   │   ├── dispersal_natal_complete.RData
│   │   └── dispersal_breeding_complete.RData
│   └── .gitkeep                       # To maintain empty folders in git
│
├── results/                           # Model results
│   ├── weibull/
│   │   ├── average/
│   │   ├── natal/
│   │   └── breeding/
│   └── combined/
│       ├── short_models_combined.RData
│       └── complete_models_combined.RData
│
├── figures/                           # Generated figures
│   ├── forest_plots/
│   ├── predictions/
│   └── comparisons/
│
├── manuscript/                        # Manuscript and supplementary materials
│   ├── main_text.pdf
│   ├── supplementary_material.pdf
│   └── figures/
│
└── docs/                              # Additional documentation
    ├── TECHNICAL_SUMMARY.md
    ├── METHODOLOGY.md
    └── CHANGELOG.md
```

## Essential files for GitHub

### Include:
- ✅ All scripts (.R)
- ✅ README.md and documentation
- ✅ Small example files
- ✅ Folder structure (with .gitkeep)
- ✅ LICENSE

### Don't include (use .gitignore):
- ❌ Large .RData files
- ❌ Complete raw data (if >50MB)
- ❌ Intermediate results
- ❌ High-resolution figures
- ❌ Temporary files

## Alternatives for large data

1. **Zenodo/Figshare**: Upload complete datasets
2. **GitHub Releases**: For specific data versions
3. **Git LFS**: For large files (if necessary)
4. **README with links**: Links to external data repositories

## Suggested commits

```bash
git add revision3/*.R
git commit -m "Add revision 3 analysis scripts"

git add README.md .gitignore
git commit -m "Add documentation and gitignore"

git add data/raw/.gitkeep data/processed/.gitkeep
git commit -m "Add data folder structure"
```

## Suggested main README

The repository's main README.md should include:

1. **Title and authors**
2. **Project description**
3. **Paper citation**
4. **Repository structure**
5. **Requirements** (R version, packages)
6. **Usage instructions**
7. **Data links** (Zenodo, etc.)
8. **License**
9. **Contact**

## Recommended license

For research code: **MIT License** or **GPL-3.0**

To include in repository:
```bash
# Add LICENSE file
# Option 1: MIT (more permissive)
# Option 2: GPL-3.0 (copyleft)
```
