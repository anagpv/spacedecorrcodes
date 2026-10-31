#---------------------------#
#- Preprocessing Lung data -#
#---------------------------#

library(tidyverse)
library(Giotto)
library(smiDE)

dir.data <- "/path/to/zenodo/download/"   # where Lung_Giotto.RData lives
dir.out  <- "/path/to/output/"            # where processed outputs will be saved

# Raw data: Lung_Giotto.RData must be downloaded from:
# https://brukerspatialbiology.com/products/cosmx-spatial-molecular-imager/ffpe-dataset/nsclc-ffpe-dataset/
# Place the downloaded file in dir.data before running this script.

# Change to be adequate to newer versions
load(paste0(dir.data, "Lung_Giotto.RData"))
new_gem <- createGiottoObject(
  expression    = list(
    raw        = gem@expression$rna$raw
  ),
  spatial_locs = gem@spatial_locs$raw,
  cell_metadata = gem@cell_metadata$rna
)
new_gem <- normalizeGiotto(new_gem,
                           spat_unit = "cell", feat_type = "rna",
                           expression_values = "raw",
                           norm_methods = "standard",
                           scalefactor = 6000, log_norm = TRUE,
                           scale_feats = FALSE, scale_cells = FALSE)

# Filtering Giotto ------------------------------------------------------------
new_gem@cell_metadata$cell$rna$niche <- ifelse(pDataDT(new_gem, spat_unit = "cell")$niche %in%
                                                 c("stroma","myeloid-enriched stroma","plasmablast-enriched stroma"),"stroma",
                                               pDataDT(new_gem, spat_unit = "cell")$niche)

# Getting subset of cells
selected_cells <- pDataDT(new_gem, spat_unit = "cell") %>%
  as.data.frame %>%
  filter(tissue %in% c("Lung12","Lung13","Lung5-4","Lung6"),
         niche %in% c("tumor interior","tumor-stroma boundary","stroma"),
         cell_type %in% c("tumor 12","tumor 13","tumor 5","tumor 6","tumor 9")
  ) %>%
  pull(cell_ID)

gem <- subsetGiotto(new_gem, cell_ids = selected_cells)

# create network
gem <- createSpatialNetwork(gobject = gem, minimum_k = 2)
save(gem, file = file.path(dir.data, "Lung_Giotto_Filtered_Tumor.RData"))

# Get cell contamination metric -------------------------------------------
selected_cells <- pDataDT(new_gem, spat_unit = "cell") %>%
  filter(tissue %in% c("Lung12","Lung13","Lung5-4","Lung6")) %>%
  pull(cell_ID)
new_gem <- subsetGiotto(new_gem, cell_ids = selected_cells)

new_gem@cell_metadata$cell$rna$cell_type <-ifelse(grepl("tumor",pDataDT(new_gem, spat_unit = "cell")$cell_type),"tumor",
                                                  pDataDT(new_gem, spat_unit = "cell")$cell_type)

overlap_metrics <-
  smiDE::overlap_ratio_metric(assay_matrix = getExpression(new_gem,spat_unit = "cell",
                                                           feat_type = "rna",values = "normalized")@exprMat
                              ,metadata = merge(pDataDT(new_gem, spat_unit = "cell"),
                                                getSpatialLocations(new_gem,spat_unit = "cell")@coordinates)
                              ,cellid_col = "cell_ID"
                              ,cluster_col = "cell_type"
                              ,sdimx_col = "sdimx"
                              ,sdimy_col = "sdimy"
                              ,grouping_col = "tissue"
                              ,radius = 0.05
                              ,verbose = T
  )
write.csv(overlap_metrics, file = file.path(dir.out, "Lung/ContaminationMetricsGenes_Lung.csv"))

pre_de_obj <-
  pre_de(adjacencies_only = F,
         ,counts = getExpression(new_gem,spat_unit = "cell",
                                 feat_type = "rna",values = "raw")@exprMat
         ,ref_celltype = "tumor"
         ,metadata = merge(pDataDT(new_gem, spat_unit = "cell"),
                           getSpatialLocations(new_gem,spat_unit = "cell")@coordinates)
         ,cell_type_metadata_colname = "cell_type"
         ,split_neighbors_by_colname = "tissue"
         ,mm_radius = 0.05
         ,sdimx_colname = "sdimx"
         ,sdimy_colname = "sdimy"
         ,verbose=TRUE
  )


# SVGs --------------------------------------------------------------------
for(tissue in c("Lung5-4","Lung12","Lung13","Lung6")){
  id_tissue <- pDataDT(new_gem, spat_unit = "cell")$cell_ID[gem@cell_metadata$tissue == tissue]
  gem_tissue <- subsetGiotto(gem, cell_ids = id_tissue)
  km_spat_pct <- binSpect(gem_tissue, bin_method = 'kmeans')
  save(km_spat_pct, file = file.path(dir.out, tissue, paste0("km_spatialgenes_",tissue,".Rdata")))
}


