#----------------------#
#- Getting Lupus data -#
#----------------------#

library(readxl)
library(spacedecorr)
library(tidyverse)
library(sctransform)
library(Giotto)

# Raw Lupus data must be downloaded from: https://doi.org/10.6084/m9.figshare.c.7373860
# Place downloaded files in dir.data before running
dir.data <- "/path/to/figshare/download/"   # where cleaneddata.RData lives
dir.out  <- "/path/to/output/"              # where processed outputs will be saved
dir.proj <- "/path/to/repo/"               # root of this GitHub repo

source(file.path(dir.proj, "R", "FunctionsSmooth.R"))

nCores <- 10
args <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

#### load data --------------------
load(file.path(dir.data, "cleaneddata.RData"))
rownames(customlocs) = annot$cell_ID

# code to get contamination metrics in the end
contam <- read.csv(file.path(dir.out, "ContaminationMetricsGenes_Lupus.csv"))

genes_PCT <- contam %>% filter(clusts == "PCT", ratio <  1) %>% pull(target)
genes_Podocyte <- contam %>% filter(clusts == "Podocyte", ratio < 1) %>% pull(target)
rownames(annot) <- annot$cell_ID

# Parameters -----
tissues <- sort(unique(annot$tissuename), decreasing = T)
tissues <- tissues[!tissues %in% c("SLE8.2","SLE8.3")]

params <- expand_grid(tissues = tissues,
                      method = c("vst","giotto", "svg", "spacedecorr","vst+spacedecorr"),
                      k = c(1:3)) %>%
  filter(!(method %in% c("vst","giotto","svg") & k %in% c(2)))
tissue_i <- params$tissues[args]
method <-params$method[args]
kid <- params$k[args]
ksplines <- c(200, 1000)[kid]
basis <- "ts"
dir.create(file.path(dir.out, tissue_i), recursive = TRUE, showWarnings = FALSE)

# PCT ---------------------------------------------------------------------
id_tissue <- annot$cell_ID[which(annot$tissuename == tissue_i & annot$clusts == "PCT")]

raw_tissue <- t(as.matrix(raw)[,id_tissue])
locs_tissue <- as.data.frame(customlocs[id_tissue,])
annot_tissue <- annot[id_tissue,]
annot_tissue <- cbind(annot_tissue, locs_tissue)
rownames(annot_tissue) <- annot_tissue$cell_ID

# vst
if(method == "vst"){
  time_init <- Sys.time()
  data_vst <- vst(umi=t(raw_tissue[,genes_PCT]),
                  cell_attr = annot_tissue,
                  latent_var = "totalcounts"
  )
  time_end <- Sys.time()
  time_vst <- as.numeric(time_end - time_init,units = "secs")
  save(data_vst, time_vst, file = file.path(dir.out, tissue_i, paste0("Data_PCT_",tissue_i,"_vst.Rdata")))
}


# spacedecorr
if(method == "spacedecorr"){
  time_spacedecorr <- system.time(data_spacedecorr <-
                                    spacedecorr(
                                      assay_matrix = raw_tissue[,genes_PCT],
                                      metadata =  annot_tissue,
                                      libsize_col = "totalcounts",
                                      k=ksplines,
                                      basis = basis,
                                      nCores = nCores-1,verbose = T
                                    ))
  save(data_spacedecorr, time_spacedecorr, file = file.path(dir.out, tissue_i, paste0("Data_PCT_",tissue_i,"_spacedecorr_k",ksplines,"_basis",basis,".Rdata")))
}

# Giotto
if(method == "svg"){
  locs_tissue$cell_ID <- rownames(locs_tissue)
  gem_pct <- createGiottoObject(raw_exprs = t(raw_tissue),
                                spatial_locs =  locs_tissue,
                                cell_metadata = annot_tissue
  )
  gem_pct <- createSpatialNetwork(gobject = gem_pct, minimum_k = 2)
  gem_pct <- normalizeGiotto(gem_pct, scalefactor = 10e4)
  km_spat_pct <- binSpect(gem_pct, bin_method = 'kmeans')
  save(km_spat_pct, file = file.path(dir.out, tissue_i, paste0("km_spat_SVG_pct_",tissue_i,".Rdata")))
}

if(method == "giotto"){
  time_giotto <- system.time(data_giotto <- get_smoothed(locs_tissue,raw_tissue[,genes_PCT], 0.1 ))
  save(time_giotto, data_giotto, file = file.path(dir.out, tissue_i, paste0("Data_PCT_",tissue_i,"_Giotto.Rdata")))
}


# Get contamination metrics -----------------------------------------------
rownames(customlocs) = annot$cell_ID
rownames(annot) <- annot$cell_ID

neighbors <- InSituCor:::nearestNeighborGraph(x = customlocs[, 1], y = customlocs[, 2], N = 50, subset = annot$tissuename)
neighborssymm <- 1 * ((neighbors + Matrix::t(neighbors)) != 0)
rownames(neighborssymm) <- colnames(neighborssymm) <- colnames(raw)
contam <- smiDE:::contamination_ratio_metric(assay_matrix = raw,
                                             metadata = data.frame(sdimx = customlocs[, 1],
                                                                   sdimy = customlocs[, 2],
                                                                   cell_ID = annot$cell_ID,
                                                                   tissue = annot$tissuename,
                                                                   clusts = annot$clusts),
                                             adjacency_matrix = neighborssymm,
                                             cluster_col = c("clusts"),
                                             cellid_col = "cell_ID",
                                             grouping_col = tissue,
                                             sdimx_col = "sdimx",
                                             sdimy_col = "sdimy",
                                             verbose = TRUE)

write.csv(contam, file = file.path(dir.out, "ContaminationMetricsGenes_Lupus.csv"))
