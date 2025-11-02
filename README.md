# Simple mechanistic traits outperform complex syndromes in predicting avian dispersal distances

**Guillermo Fandos<sup>1,2</sup>,Rob Robinson <sup>3</sup>, Damaris Zurell <sup>1</sup>**

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10713958.svg)](https://zenodo.org/records/10713958)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Corresponding author**: Guillermo Fandos  
**Email**: gfandos@ucm.es

1. Institute for Biochemistry and Biology, University of Potsdam, D-14469 Potsdam, Germany 
2. Department of Biodiversity ecology and evolution, Faculty of Biology, Complutense University, 28040 Madrid, Spain
3.  British Trust for Ornithology, The Nunnery, Thetford, IP24 2PU, United Kingdom


Code (preliminar): [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10713958.svg)](https://zenodo.org/records/10713958)

---

## Overview

This repository contains all R code and analysis pipelines for the manuscript:


### Key findings

- Multi-trait dispersal syndromes help understand trait covariation within related species
- Single-trait models (body mass, life history) better predict dispersal across distant orders
- Flight efficiency (HWI) is important for natal and breeding dispersal
- Different drivers for median vs. long-distance dispersal

---

## Repository structure

### 1. `data/`
All dispersal estimates come from [Fandos et al. 2023](https://doi.org/10.1111/1365-2656.13836) (*Journal of Animal Ecology*).

Contains:
- **`raw/`**: Original phylogenetic trees, trait data, and dispersal distances
- **`processed/`**: Clean datasets with phylogenetic covariance matrices for analysis

### 2. `scripts/` or `revision3/`
All executed code, including pipelines, analysis scripts, and figure generation files.

**The pipeline and workflow for analyzing dispersal syndromes:**

1. **Join dispersal estimates with biological and ecological traits** (`prepare_disperal_age_R3.R`)
2. **Check correlations between variables** (included in main workflow)
3. **Fit Bayesian phylogenetic mixed models** (`function_dispersal_R3R.R`)
4. **Extract results and create plots from dispersal syndromes** (`load_models_R3.R`, plotting scripts)
5. **Evaluate and cross-predict dispersal estimates** (model comparison scripts)

See [`revision3/README.md`](revision3/README.md) for detailed script descriptions.

### 3. `results/`
Model outputs and combined results:
- **`weibull/`**: Individual models by dispersal type (average, natal, breeding)
- **`combined/`**: Merged models and variable selection summaries

### 4. `figures/`
Generated plots:
- Forest plots
- Conditional effects
- Prediction plots
- Model comparisons

### 5. `manuscript/`
Manuscript and supplementary materials (PDFs)

### 6. `docs/`
Additional documentation and technical details

---

## Quick start

### 1. Install dependencies
```r
source("revision3/function_dispersal_R3R.R")
# Packages will be installed automatically
```

### 2. Prepare data
```r
source("revision3/prepare_disperal_age_R3.R")
```

### 3. Run analyses
```r
source("revision3/function_dispersal_R3R.R")

# Run for each dispersal type
run_dispersal_analysis(dispersal_type = "average")
run_dispersal_analysis(dispersal_type = "natal")
run_dispersal_analysis(dispersal_type = "breeding")
```

### 4. Generate figures
```r
source("revision3/load_models_R3.R")
source("revision3/complete_forest_plot_R3.R")
source("revision3/5c_prediction_plots_final_R3.R")
```

For detailed instructions, see [`QUICKSTART.md`](QUICKSTART.md).

---

## Software requirements

- R >= 4.2.0
- Key packages: `brms`, `projpred`, `ape`, `phytools`, `tidyverse`, `ggplot2`, `DHARMa`, `performance`

---

## Data availability

### Dispersal data
Standardized dispersal kernels from [Fandos et al. 2023](https://doi.org/10.1111/1365-2656.13836):
- 234 species (total dispersal)
- 113 species (breeding dispersal)
- 122 species (natal dispersal)

### Trait data
- **Morphology**: [Sheard et al. 2020](https://doi.org/10.1038/s41467-020-16313-6) (*Nature Communications*)
- **Life history**: [Storchová & Hořák 2018](https://doi.org/10.1111/geb.12709) (*Global Ecology and Biogeography*)
- **Ecology**: [Reif et al. 2016](https://doi.org/10.1111/oik.02276) (*Oikos*)
- **Phylogeny**: [Jetz et al. 2012](https://doi.org/10.1038/nature11631) (*Nature*)

Complete datasets available at: [Zenodo DOI - pending]

---

## Citation

If you use this code, please cite:

```bibtex
@article{fandos2024dispersal,
  title={Dispersal syndromes allow understanding but not predicting dispersal ability across the tree of life},
  author={Fandos, Guillermo and Robinson, Robert A and Zurell, Damaris},
  journal={bioRxiv},
  pages={2024.04.01.587575},
  year={2024},
  publisher={Cold Spring Harbor Laboratory},
  doi={10.1101/2024.04.01.587575}
}
```

And the dispersal data:
```bibtex
@article{fandos2023standardised,
  title={Standardised empirical dispersal kernels emphasise the pervasiveness of long-distance dispersal in European birds},
  author={Fandos, Guillermo and Claramunt, Santiago and Martin, Kathy and Derryberry, Elizabeth P and Iris, Dainson and Lisovski, Simeon},
  journal={Journal of Animal Ecology},
  volume={92},
  number={1},
  pages={158--170},
  year={2023},
  doi={10.1111/1365-2656.13836}
}
```

---

## License

This project is licensed under the MIT License - see the [`LICENSE`](LICENSE) file for details.

Data are subject to their original licenses.

---

## Funding

This work was funded by the German Science Foundation DFG (grant no. ZU 361/1-1).

---

## Contact

**Guillermo Fandos**  
- Email: gfandos@ucm.es
- GitHub: [Update with your GitHub]
- ORCID: [Update with your ORCID]

---

## Acknowledgments

Bird ringing data provided by the EURING Data Bank. We thank all ringers and coordinators who have contributed data over decades.

---

**Last updated**: November 2025
