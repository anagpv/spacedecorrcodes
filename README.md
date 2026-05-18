# SpaceDecorr Analysis

Reproducibility code for Vasconcelos, Ana Gabriela, et al. ["Accounting for Spatial Structure in Network Analysis of Spatial Transcriptomics Data."](https://www.biorxiv.org/content/10.1101/2025.07.23.666450v1) bioRxiv (2025): 2025-07.
This repository contains all R scripts to reproduce the analyses and figures.

## Data

Intermediate data files are available on Zenodo: [https://doi.org/10.5281/zenodo.20275986 ](https://doi.org/10.5281/zenodo.20276257) 

Raw Lupus nephritis data: https://doi.org/10.6084/m9.figshare.c.7373860  

Raw NSCLC CosMx data: https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/

Download all data files and set the paths at the top of each script before running.

## Repository structure

    R/                        # Helper functions sourced by analysis scripts
    └── FunctionsSmooth.R
    
    scripts/                  # Analysis scripts
    ├── 00_preprocess_lung.R          # Preprocessing: cleaning lung NSCLC data
    ├── 01_get_spacedecorr_lung.R     # SpaceDecorr: lung NSCLC data
    ├── 02_results_real_lung.R        # Figures: lung real data analysis
    ├── 03_sim_pairwise.R             # Simulation: pairwise correlation
    ├── 04_sim_pairwise_robustcor.R   # Simulation: robust correlation
    ├── 05_results_pairwise.R         # Figures: pairwise simulation results
    ├── 06_sim_wgcna.R                # Simulation: module recovery (WGCNA)
    ├── 07_results_wgcna.R            # Figures: WGCNA simulation results
    ├── 08_sim_diffnet_niche.R        # Simulation: differential network (niche)
    ├── 09_results_diffnet_niche.R    # Figures: diffnet niche results
    ├── 10_sim_diffnet_samples.R      # Simulation: differential network (samples)
    ├── 11_results_diffnet_samples.R  # Figures: diffnet samples results
    ├── 12_get_spacedecorr_lupus.R    # Spacedecorr: lupus nephritis data
    └── 13_results_real_lupus.R       # Figures: lupus real data analysis



