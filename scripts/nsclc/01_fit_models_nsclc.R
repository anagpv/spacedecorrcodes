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
library(data.table)
library(spMOCA)
source(file.path(dir.proj, "R", "FunctionsSmooth.R"))

# Designed to run as a SLURM array job. To run locally, set args <- 1
nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
load(file.path(dir.data, "Lung_Giotto_Filtered_Tumor_GiottoSuites.RData"))

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)

params <- expand_grid(tissues = "Lung6", #c("Lung5-4","Lung12","Lung13","Lung6"),
                      method = c("giotto","spacedecorr","vst","spmoca"),
                      k = c(1:3)
) %>% filter(!(method %in% c("vst","giotto","svg") & k != 1))
tissue <- params$tissues[args]
method <-params$method[args]
kid <- params$k[args]
basis <- "ts"

dir.create(file.path(dir.out, tissue), recursive = TRUE, showWarnings = FALSE)

# Specific tissue ----
print(tissue)
id_tissue <- pDataDT(gem, spat_unit = "cell")$cell_ID[pDataDT(gem, spat_unit = "cell")$tissue == tissue]
gem_tissue <- subsetGiotto(gem, cell_ids = id_tissue)
rm(list = "gem")

ksplines <- c(100, 200, 1000)[kid]

# Contamination
contamination0 <- read.csv(paste0(dir.out,"ContaminationMetricsGenes_Lung.csv"))
contamination <- contamination0 %>% filter(avg_cluster > 0.1, ratio < 1, tissue == "Lung6", cell_type == "tumor")
filtered_genes <- contamination$target

gem_tissue <- subsetGiotto(gem_tissue, feat_ids = filtered_genes)
locs <- as.data.frame(gem_tissue@spatial_locs$cell$raw@coordinates)
metadata <- as.data.frame(cbind(locs,
                                niche = pDataDT(gem_tissue, spat_unit = "cell")$niche,
                                totalcounts = pDataDT(gem_tissue, spat_unit = "cell")$totalcounts))
rownames(metadata) <- metadata$cell_ID

# Giotto ---------------------------------------------------------------
f(method == "giotto"){
  print("Giotto")
  # All
  km_spatialgenes <- binSpect(gem_tissue, bin_method = 'kmeans')
  save(km_spatialgenes, file = file.path(dir.out, tissue, paste0("km_spatialgenes_",tissue,".RData")))
  
  time_giotto <- system.time(cor_giotto <-  detectSpatialCorFeats(gem_tissue,
                                                                  method = 'network',
                                                                  spatial_network_name = 'Delaunay_network'))
  save(cor_giotto, time_giotto, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_giotto.Rdata")))
}

# SpaceDecorr ---------------------------------------------------------------
if(method == "spacedecorr"){
  print("SpaceDecorr")
  time_spacedecorr <- system.time(data_spacedecorr <-
                                    spacedecorr(
                                      assay_matrix = t(as.matrix(gem_tissue@expression$cell$rna$normalized@exprMat)),
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
  data_vst <- vst(umi=as.matrix(gem_tissue@expression$cell$rna$normalized@exprMat),
                  cell_attr = pDataDT(gem_tissue, spat_unit = "cell")$niche,
                  latent_var = "totalcounts"
  )
  time_end <- Sys.time()
  time_vst <- as.numeric(time_end - time_init,units = "secs")
  save(data_vst, time_vst, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_vst.Rdata")))
}

# Aggregate into spots ----------------------------------------------------
if(method == "spmoca"){
  counts <- getExpression(gem_tissue,spat_unit = "cell",feat_type = "rna",values = "raw")@exprMat
  norm <- getExpression(gem_tissue,spat_unit = "cell",feat_type = "rna",values = "normalized")@exprMat
  locs <- getSpatialLocations(gem_tissue,spat_unit = "cell")@coordinates
  rownames(locs) <- locs$cell_ID
  colnames(locs)[1:2] <- c("x", "y")
  rm(ls = "gem_tissue")
  
  sum_cols_by_group <- function(mat, grp) t(rowsum(t(mat), grp))
  aggregate_pixels <- function(coords, counts, meta, n_pixels = 2500) {
    
    xr <- range(coords$x)
    yr <- range(coords$y)
    
    # number of bins
    nx <- round(sqrt(n_pixels))
    ny <- round(sqrt(n_pixels))
    
    # bin size
    pixel_size_x <- diff(xr) / nx
    pixel_size_y <- diff(yr) / ny
    
    # assign each cell to a pixel
    ix <- pmin(floor((coords$x - xr[1]) / pixel_size_x) + 1L, nx)
    iy <- pmin(floor((coords$y - yr[1]) / pixel_size_y) + 1L, ny)
    
    # unique pixel id
    pid <- (iy - 1L) * nx + ix
    
    fac <- as.integer(factor(pid))
    K <- max(fac)
    
    # aggregate counts
    agg <- sum_cols_by_group(counts, fac)
    
    spot_ids <- paste0("px_", seq_len(K))
    colnames(agg) <- spot_ids
    
    # recover x/y indices of occupied pixels
    lev <- as.integer(levels(factor(pid)))
    px <- ((lev - 1L) %% nx) + 1L
    py <- ((lev - 1L) %/% nx) + 1L
    
    spot_positions <- data.frame(
      spot_id = spot_ids,
      x = xr[1] + (px - 0.5) * pixel_size_x,
      y = yr[1] + (py - 0.5) * pixel_size_y,
      totalcounts = as.numeric(rowsum(meta$totalcounts, fac))
    )
    rownames(spot_positions) <- spot_positions$spot_ids
    cell_to_spot <- factor(spot_ids[fac], levels = spot_ids)
    
    list(
      counts = agg,
      cell_to_spot = cell_to_spot,
      spot_positions = spot_positions,
      qc = list(
        nx = nx,
        ny = ny,
        n_spots = K
      )
    )
  }
  
  spots <- aggregate_pixels(locs, as.matrix(counts), metadata, n_pixels = 100^2)
  
  cbind(t(spots$counts),spots$spot_positions) %>%
    ggplot(aes(x = x, y = y, col = VEGFA)) +
    geom_point()
  
  ## Run spMOCA ------------------------------------------------------------------
  t0 <- proc.time()
  # Create object
  object <- createSPMOCAObject(spatial_countMat = as.matrix(spots$counts),
                               spatial_location = spots$spot_positions[,-1],
                               feature_list = rownames(counts))
  # Normalization
  object <- normalizeSPMOCA(object)
  
  # Feature Selection
  # Contamination
  object <- featureSelectionSPMOCA(object,
                                   feature.selection = NULL,
                                   feature.list = contamination$target)
  
  # Construct Spatial Kernel
  object <- createSpatialKernel(object)
  
  # Estimate Gene-Gene Correlation
  object <- estimatespMOCA(object)
  cor_spmoca <- object@corr_est
  
  time_spmoca <- unname((proc.time() - t0)["elapsed"])
  save(cor_spmoca, time_spmoca, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_spmocaspots.Rdata")))
  
  
  # SpaceDecorr spot based ------------------------------------------------------
  ksplines <- 200
  print("SpaceDecorr")
  
  rownames(spots$spot_positions) <- spots$spot_positions$spot_id
  colnames(spots$spot_positions) <- c("spot_id","sdimx","sdimy")
  time_spacedecorr <- system.time(data_spacedecorr <-
                                    spacedecorr(
                                      assay_matrix = t(as.matrix(spots$counts)),
                                      metadata = spots$spot_positions,
                                      libsize_col = "totalcounts",
                                      k = ksplines,
                                      nCores = 1,verbose = T
                                    ))
  
  save(data_spacedecorr,time_spacedecorr, file = file.path(dir.out, tissue, paste0("Data_",tissue,"_spacedecorrspots_k",ksplines,".Rdata")))
}