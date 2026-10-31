# SpaceDecorr analysis

Code to reproduce the analyses, figures and tables of the revised manuscript

> Vasconcelos AG, Li W, Danaher P, Wakefield J, Shojaie A. Accounting for Spatial Structure in Network Analysis of Spatial Transcriptomics Data. *Under review.*

The SpaceDecorr R package is available at [github.com/anagpv/spacedecorr](https://github.com/anagpv/spacedecorr) and can be installed with

```r
remotes::install_github("anagpv/spacedecorr")
```

## Repository structure

    R/
    ├── Modified_cor_test.R                 # modified t-test (effective sample size), used in the NSCLC and liver analyses
    └── Verify_spacedecorr_numerically.R    # numerical check of the SpaceDecorr implementation against mgcv
    scripts/
    ├── nsclc/                              # CosMx NSCLC analysis (Figure 2)
    ├── liver/                              # MERFISH liver analysis (Figure 4)
    └── simulations/                        # simulation study (Figure 3), see scripts/simulations/README.md

## Data

**Raw data** (download from the original sources)

- CosMx NSCLC: [Bruker Spatial Biology](https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/) (He et al., *Nat Biotechnol* 2022). We used section Lung6.
- MERFISH human liver: [Dryad, doi:10.5061/dryad.37pvmcvsg](https://doi.org/10.5061/dryad.37pvmcvsg) (Watson et al., *Nat Commun* 2025). 

**Reference databases for the NSCLC analysis** (known interactions and negative-control pairs)

- BioGRID v4.4.245 (`BIOGRID-ALL-4.4.245.tab3.zip`), STRING v12, MSigDB
- TRRUST (`trrust_rawdata.human.tsv`), DoRothEA (confidence A–C, from the R package), ReMap (`ReMap_ChIP-seq.gmt`)

**Intermediate data**: fitted models, correlations and simulation results are on Zenodo ([doi:10.5281/zenodo.23841454](https://doi.org/10.5281/zenodo.23841454)), so figures can be reproduced without refitting.

## Running the code

Each script sets its directories at the top. Edit them before running:

```r
dir.proj <- "/path/to/repo/"              # this repository
dir.data <- "/path/to/data/"              # raw data, databases and the Zenodo download
dir.out  <- "/path/to/output/"            # where results and figures are written
```

Run the scripts in each folder in numerical order.

### NSCLC (`scripts/nsclc/`)

| Script | Does | Paper |
|---|---|---|
| `00_preprocess_nsclc.R` | selects tumor cells of Lung6 and filters genes | |
| `FindIndependentGenes.R` | builds known-interaction, high-confidence and negative-control pairs from the databases | |
| `01_fit_models_nsclc.R` | fits SpaceDecorr (K = 100, 200, 1000) and the comparison methods | |
| `02_correlations_nsclc.R` | correlations, modified t-test, modules, enrichment and module shares | |
| `03_plot_figure2_supplement_nsclc.R` | figures and tables | Figure 2, Supp. Figures S1–S3, SN6–SN7, Supp. Table S2 |

### Liver (`scripts/liver/`)

| Script | Does | Paper |
|---|---|---|
| `01_fit_models_liver.R` | fits SpaceDecorr (K = 100, 200, 1000), VST and spMOCA per sample | |
| `02_correlations_liver.R` | correlations, replicated networks, modules, comparison between conditions, leave-one-out, modified t-test | |
| `03_plot_figure4_supplement_liver.R` | figures and tables | Figure 4, Supp. Figures S5–S6, SN5, SN9, Supp. Table S3 |
| `plot_ess.R` | effective sample size and modified t-test (NSCLC, liver and simulations) | Supp. Figure SN10 |

### Simulations (`scripts/simulations/`)

See [`scripts/simulations/README.md`](scripts/simulations/README.md). They produce Figure 3 and Supplementary Figures S4, SN1–SN4 and SN8.

## Packages

R 4.6.0. SpaceDecorr: spacedecorr (0.1.0), mgcv (1.9.4). Comparison methods: sctransform (0.4.3, the VST-normalized method), Giotto (4.2.2), spMOCA (1.0). Modules and networks: clusterProfiler (4.20.0), igraph (2.3.2). Figures: tidyverse (2.0.0), cowplot (1.2.0). Full session information is in `sessionInfo.txt`, and the simulations use additional packages listed in `scripts/simulations/README.md`.

## Citation

Until the paper is published, please cite the preprint:

> Vasconcelos AG, Danaher P, Wakefield J, Shojaie A. Accounting for Spatial Structure in Network Analysis of Spatial Transcriptomics Data. *bioRxiv* (2025). [doi:10.1101/2025.07.23.666450](https://doi.org/10.1101/2025.07.23.666450)