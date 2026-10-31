# Simulation study: code

Every script sets three directories at the top; edit them before running.

```r
dir.proj <- "/path/to/repo/"              # this folder
dir.data <- "/path/to/zenodo/download/"   # the data deposit (results and the two real-data inputs)
dir.out  <- "/path/to/output/"            # where a run writes its results and a figure script its PDFs
```

A `run_*.R` script simulates and applies the methods; it writes one CSV to `dir.out` and can be stopped and restarted (settings already written are skipped). A `figure_*.R` script reads the CSV of the deposit from `dir.data` and writes PDFs to `dir.out/figures`. The helpers are in `util/`.

| run script | writes | figure script | figure |
|---|---|---|---|
| `run_pairwise.R` | `pairwise_results.csv` | `figure_pairwise.R` | main figure panels A (type 1 error) and B (power); supplementary pairwise figure (SpaceDecorr basis dimension) |
| `run_niche.R` | `niche_results.csv` | `figure_niche.R` | main figure panel C (niche map, niche power); supplementary niche-sensitivity figure; supplementary scenario II error figure |
| `run_block.R` | `block_results.csv` | `figure_block.R` | supplementary block co-expression figure |
| `run_diffnet_niche.R` | `diffnet_niche_results.csv` | `figure_diffnet.R` | supplementary differential-network figure, panels B and C |
| `run_diffnet_niche_spacex.R` | `diffnet_niche_spacex_results.csv` | `figure_diffnet.R` | panels B and C (SpaceX) |
| `run_diffnet_samples.R` | `diffnet_samples_results.csv` | `figure_diffnet.R` | panels D and E |

Environment variables of the run scripts: `N_CORES` (workers; default all cores), `MAX_TASKS` (run the first n tasks only), `DRY_RUN=1` (list the tasks and stop). `run_diffnet_niche_spacex.R` takes `SPACEX_CORES` and `N_TASKS` instead of `N_CORES`.

## Packages

R >= 4.3. Methods: spacedecorr, spMOCA, MERINGUE, Giotto (>= 4), GraphR, sctransform (the VST-normalized method); SpaceX, MSFA and PQLseq for `run_diffnet_niche_spacex.R` only (PQLseq is in `PQLseq/`: `R CMD INSTALL PQLseq`). Simulation: fields, MASS, igraph. Results: filelock, data.table, R.utils (reads the gzipped block file). Figures: tidyverse, cowplot, patchwork, ggh4x, mclust, RColorBrewer; the PDFs are written with `cairo_pdf` in Arial.
