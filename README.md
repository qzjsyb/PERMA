# Metabolomics Analysis Pipeline for Soil Metabolites as Early Warning Indicators of Permafrost Thaw Transitions


This repository contains R scripts and associated datasets for analyzing metabolomics data, including mosiac diversity, logistic regression modeling, hierarchical clustering, and temporal analyses.

## 📁 Repository Contents

### 🔬 R Scripts

- **`0-ordination.R`**: Performs ordination analyses as presented in Supplementary Figures 1 and 2, and Figures 2A, 2C, 2D, Supplementary Figures 4 and 5.

- **`1-mu_diversity.R`**: Calculates mu-diversity metrics and generates plots corresponding to Figures 2B and 2E.

- **`2-logistic_model.R`**: Builds logistic regression models for classification tasks, supporting Figure 3E.

- **`2.1-predictive_feature_analysis.R`**: Analyzes and visualizes predictive features, supporting Figures 3B, 3C, 3D.

- **`2.2-logistic_by_depth.R`**: Constructs logistic regression models stratified by depth and compares them to the main model, generating Supplementary Figure 3.

- **`3-clustering.R`**: Implements Weighted Correlation Network Analysis (WGCNA) and produces Figure 4.

- **`4-temporal_analysis.R`**: Conducts linear regression analyses on predictive features over time, resulting in Figure 5.

- **`optional_IS_Frag_Finder_new.R`**: Optional script to identify insource fragments. Mgf file will be needed. Insource fragment listed are provided.
- **`optional_core_cleaning.R`**:

### 📊 Data Files

- **`compound_link.csv`**: Contains raw feature properties.

- **`HNEG.csv`** and **`RPPOS.csv`**: Datasheets linked to `compound_link.csv` via the `compound_ID` column.

- **`metadata.all.csv`**: Stores all sample metadata.

- **`in_source_HN.csv`** and **`in_source_RP.csv`**: List in-source fragments present in raw data that should be removed during downstream analysis.
