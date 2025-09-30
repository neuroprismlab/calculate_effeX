# Calculate_EffeX
[![DOI](https://zenodo.org/badge/758229305.svg)](https://doi.org/10.5281/zenodo.17238826)

A pipeline for calculating neuroimaging effect sizes for the [BrainEffeX](neuroprismlab.shinyapps.io/braineffex) platform.

## Overview

This repository provides pipelines to process neuroimaging data from subject-level statistical maps through to effect size results, meta-analyses results, and visualizations of effects. 

## Pipeline Structure

The analysis pipeline consists of three main stages:

### 1. Group-level Analysis - `group_level/` 
**Input:** Subject-level statistical maps from data contributors  
**Output:** Group-level statistical maps and test results for each study

- **Main script:** `do_group_level.m`
- **Demo script:** `group_level_example.mlx`
- **Analysis options:**
  - Univariate vs multivariate statistical tests
  - Parcel-wise vs network-level spatial pooling
  - Motion correction strategies (none, regression, thresholding)
- **Supported tests:** Correlation (r), one-sample t-test (t), two-sample t-test (t2), and multivariate versions of each
- **Language:** MATLAB

### 2. Effect Size Calculation + Meta-Analysis - `effect_size/`
**Input:** Group-level statistical maps from stage 1  
**Output:** Effect sizes for each study and meta-analysis results

- **Main script:** `master.R`
- **Example notebook:** `effect_size_example.ipynb`
- **Key functions:** Effect size calculation, meta-analysis
- **Language:** R

### 3. Visualization - `create_figures/`
**Input:** Effect size data from stage 2 (this data can also be accessed directly from [OSF](https://osf.io/cwnjd/files/osfstorage))  
**Output:** Figures and visualizations for the BrainEffeX Shiny app or related manuscripts

- **Main script:** `generate_figures.R`
- **Example notebook:** `example_generate_figures.ipynb`
- **Output:** Figures as png files
- **Language:** R

---
## Repository Structure
Relavent files:
```
calculate_effex/
├── README.md                           # information and instructions for the repo
├── LICENSE                             # License file
└── group_level                         # Compute group-level statistical maps from subject-level statistical maps
  ├── group_level_example.mlx           # interactive matlab file with example of how to run group level analysis
  ├── do_group_level.m                  # main script to run group level analysis
  └── helper_scripts                    # scripts referenced in do_group_level.m
└── effect_size                         # Compute effect maps from group-level statistical maps, run meta-analyses
  ├── master.R                          # main script to run effect size analysis
  ├── effect_size_example.ipynb         # notebook with example of how to run effect size analysis
  └── scripts                           # scripts referenced in master.R
└── create_figures                      # Create figures for BrainEffeX Shiny app
  ├── example_generate_figures.ipynb    # notebook with example of how to create figures
  ├── generate_figures.R                # main script to generate figures for BrainEffeX Shiny app  
```

---

## License

See [LICENSE](LICENSE) file for details.

## Citation

If you use this pipeline in your research, please cite the [BrainEffeX project](github.com/neuroprismlab/braineffex).

---

**Authors:** Hallee Shearer, Stephanie Noble

