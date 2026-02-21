#---------------------------------
#- Get co-expression lung data
#---------------------------------

dir.data <- "/path/to/zenodo/download/"   # folder where Zenodo data was downloaded
dir.out  <- "/path/to/output/"            # folder where outputs will be saved
dir.proj <- "/path/to/repo/"              # root of this GitHub repo

library(tidyverse)
library(Giotto)
library(Matrix)
library(spacedecorr)
library(GraphR)
library(sctransform)
library(SPARK)
library(MERINGUE)
library(data.table)
source(file.path(dir.proj, "R", "FunctionsSmooth.R"))

# Designed to run as a SLURM array job. To run locally, set args <- 1
nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
load(file.path(dir.data, "Lung_Giotto_Filtered_Tumor.RData"))

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)

params <- expand_grid(tissues = "Lung6", #c("Lung5-4","Lung12","Lung13","Lung6"),
                      method = c("smooth","giotto","spacedecorr","vst"),
                      k = c(1:4)
) %>%
  filter(!(method %in% c("vst","giotto","svg") & k != 1))
tissue <- params$tissues[args]
method <-params$method[args]
kid <- params$k[args]
basis <- "ts"

dir.create(file.path(dir.out, tissue), recursive = TRUE, showWarnings = FALSE)

# Specific tissue ----
print(tissue)
id_tissue <- gem@cell_metadata$cell_ID[gem@cell_metadata$tissue == tissue]
gem_tissue <- subsetGiotto(gem, cell_ids = id_tissue)
rm(list = "gem")

ksplines <- c(100, 200, 1000, round(0.1*dim(gem_tissue@raw_exprs)[2]))[kid]

# Contamination
contamination0 <- read.csv(paste0(dir.out,"Manuscript/RealData/Lung/ContaminationMetricsGenes_Lung.csv"))
contamination <- contamination0 %>% filter(avg_cluster > 0.1, ratio < 1, tissue == "Lung6", cell_type == "tumor")
filtered_genes <- contamination$target

gem_tissue <- subsetGiotto(gem_tissue, gene_ids = filtered_genes)
rownames(gem_tissue@cell_metadata) <- gem_tissue@cell_metadata$cell_ID


# Giotto ---------------------------------------------------------------
# Smoothed counts
if(method == "smooth"){
  time_smooth <- system.time(cor_smooth <-  detectSpatialCorGenes(gem_tissue,network_smoothing=1,
                                                                  method = 'network',
                                                                  spatial_network_name = 'Delaunay_network'))
  save(cor_smooth, time_smooth, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_smooth.Rdata")))
  }
# Default
if(method == "giotto"){
  print("Giotto")
  # All
  km_spatialgenes <- binSpect(gem_tissue, bin_method = 'kmeans')
  save(km_spatialgenes, file = file.path(dir.out, tissue, paste0("km_spatialgenes_",tissue,".RData")))

  time_giotto <- system.time(cor_giotto <-  detectSpatialCorGenes(gem_tissue,
                                                                  method = 'network',
                                                                  spatial_network_name = 'Delaunay_network'))
  save(cor_giotto, time_giotto, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_giotto.Rdata")))
  }

# SpaceDecorr ---------------------------------------------------------------
if(method == "spacedecorr"){
  print("SpaceDecorr")
  metadata <- as.data.frame(cbind(gem_tissue@spatial_locs[,c("sdimx","sdimy")],
                                  gem_tissue@cell_metadata))
  rownames(metadata) <- metadata$cell_ID

  time_spacedecorr <- system.time(data_spacedecorr <-
                                    spacedecorr(
                                      assay_matrix = t(as.matrix(gem_tissue@raw_exprs)),
                                      metadata = metadata,
                                      libsize_col = "totalcounts",
                                      k = ksplines,
                                      nCores = nCores,verbose = T
                                    ))
  save(data_spacedecorr, time_spacedecorr, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_spacedecorr_k",ksplines,".Rdata")))
  }

# VST -------------------------------------------------------------------
if(method == "vst"){
  print("vst")
  time_init <- Sys.time()
  data_vst <- vst(umi=as.matrix(gem_tissue@raw_exprs),
                  cell_attr = gem_tissue@cell_metadata,
                  latent_var = "totalcounts"
  )
  time_end <- Sys.time()
  time_vst <- as.numeric(time_end - time_init,units = "secs")
  save(data_vst, time_vst, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_vst.Rdata")))
  }
