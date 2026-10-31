#---------------------------------
#- Get co-expression lung data
#---------------------------------

dir.data <- "/path/to/zenodo/download/"   # folder where Zenodo data was downloaded
dir.out  <- "/path/to/output/"            # folder where outputs will be saved
dir.proj <- "/path/to/repo/"              # root of this GitHub repo

library(spMOCA)
library(tidyverse)
library(Matrix)
library(spacedecorr)
library(sctransform)

# Designed to run as a SLURM array job. To run locally, set args <- 1
nCores <- 5
args   <- commandArgs(trailingOnly = TRUE)
args <- as.integer(args[1])

params <- expand_grid(tissues = c("AM048", "AM061", "AM042", "AM031", "AM066", "AM062", "AM072"),
                      method = c("spacedecorr","vst","spmoca"),
                      k = c(1:3)
) %>% filter(!(method %in% c("vst","spmoca") & k != 1))
tissue <- params$tissues[args]
method <-params$method[args]
kid <- params$k[args]
basis <- "ts"

dir.create(file.path(dir.out, tissue), recursive = TRUE, showWarnings = FALSE)


# Read data ----------------------------------
meta_data <- read_csv(file.path(dir.data,"cell_properties_healthy_diseased_merfish.csv"))
raw <- as.data.frame(read_csv(file.path(dir.data,"X_raw_healthy_diseased_merfish.csv"), col_names = FALSE))
norm <- as.data.frame(read_csv(file.path(dir.data,"X_healthy_diseased_merfish.csv"), col_names = FALSE))
gene_names<- read_csv(file.path(dir.data,"gene_names_healthy_diseased_merfish.csv"), col_names = FALSE)

colnames(raw) <- colnames(norm) <- gene_names$X1
rownames(raw) <- rownames(norm) <- meta_data$index
tissues <- unique(meta_data$sample_id)
mean_expr <- sapply(tissues, function(s) {
  id <- meta_data$sample_id == s & meta_data$Cell_Type_final %in% c("Central Hep", "Portal Hep", "Fibrotic Hep 1", "Fibrotic Hep 2")
  colMeans(raw[id, ])
})
genes_keep <- rownames(mean_expr)[apply(mean_expr, 1, min) >= 0.1]

raw <- raw[,genes_keep]
norm <- norm[,genes_keep]

# Specific tissue ----
print(tissue)
id_tissue <- (meta_data$sample_id == tissue) 
meta_tissue <- meta_data[id_tissue,] %>%
  dplyr::rename(cell_ID = index, sdimx = x, sdimy = y, libsize = n_counts) %>% as.data.frame
rownames(meta_tissue) <- meta_tissue$cell_ID
raw_tissue <- raw[id_tissue,]
norm_tissue <- norm[id_tissue,]
rm(list = c("meta_data","raw","norm"))


# Hepatocytes -------
id_hep <- (meta_tissue$Cell_Type_final %in% c("Central Hep","Portal Hep","Fibrotic Hep 1", "Fibrotic Hep 2"))
meta_hep <- meta_tissue[id_hep,]

raw_hep <- raw_tissue[id_hep,]
ids_hep <- which(!is.na(rowSums(raw_hep)))
count_hep <- as.matrix(t(raw_hep[ids_hep,]))
meta_hep <- meta_hep[ids_hep,]
meta_hep$cell_type_cov <- meta_hep$Cell_Type_final
meta_hep$cell_type_cov[meta_hep$cell_type_cov %in% c("Central Hep","Portal Hep")] <- "Zonation"
ksplines_hep <- c(100, 200, 1000, round(0.1*nrow(meta_hep)))[kid]

# Zonation Hepatocytes -------
id_zonhep <- (meta_tissue$Cell_Type_final %in% c("Central Hep","Portal Hep"))
meta_zonhep <- meta_tissue[id_zonhep,]

raw_zonhep <- raw_tissue[id_zonhep,]
ids_zonhep <- which(!is.na(rowSums(raw_zonhep)))
count_zonhep <- as.matrix(t(raw_zonhep[ids_zonhep,]))
meta_zonhep <- meta_zonhep[ids_zonhep,]



# SpaceDecorr ---------------------------------------------------------------
if(method == "spacedecorr"){
  print("SpaceDecorr")
  # Hepatocytes with cell state covariate
  time_spacedecorr_hepcov <- system.time(data_spacedecorr_hepcov <-
                                           spacedecorr(
                                             assay_matrix = t(count_hep),
                                             metadata = meta_hep,
                                             libsize_col = "libsize",
                                             covariates = "cell_type_cov",
                                             k = ksplines_hep,
                                             nCores = nCores,verbose = T
                                           ))
  save(data_spacedecorr_hepcov, time_spacedecorr_hepcov, 
       file = file.path(dir.out, paste0("Data_",tissue,"_spacedecorr_k",ksplines,"_hepfibcov.Rdata")))
}

# VST -------------------------------------------------------------------
if(method == "vst"){
  print("vst")
  # Hepatocytes
  time_init <- Sys.time()
  data_vst_hep <- vst(umi = count_hep,
                      cell_attr = meta_hep,
                      latent_var = "libsize"
  )
  time_end <- Sys.time()
  time_vst_hep <- as.numeric(time_end - time_init,units = "secs")
  save(data_vst_hep, time_vst_hep, file = file.path(dir.out, paste0("Data_",tissue,"_vst_hepfib.Rdata")))
  
}


# spMOCA ------------------------------------------------------------------
if(method == "spmoca"){
  print("Hepatocytes")
  ## Zonation Hepatocytes ---
  t0 <- proc.time()
  # Create object
  obj_hep <- createSPMOCAObject(spatial_countMat = count_zonhep,
                                spatial_location = meta_zonhep[,c("sdimx","sdimy")],
                                feature_list = rownames(count_zonhep))
  # Normalization
  obj_hep <- normalizeSPMOCA(obj_hep)
  
  # Feature Selection
  obj_hep <- featureSelectionSPMOCA(obj_hep,
                                    feature.selection = NULL,
                                    feature.list = rownames(count_zonhep))
  
  # Construct Spatial Kernel
  obj_hep <- createSpatialKernel(obj_hep)
  
  # Estimate Gene-Gene Correlation
  obj_hep <- estimatespMOCA(obj_hep)
  cor_spmoca_hep <- obj_hep@corr_est
  
  time_spmoca_hep <- unname((proc.time() - t0)["elapsed"])
  save(cor_spmoca_hep, time_spmoca_hep, file = file.path(dir.out, paste0("Data_",tissue,"_spmoca_hep.Rdata")))
}
