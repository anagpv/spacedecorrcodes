#----------------------------------
#---- Real data Lung WGCNA
#----------------------------------

dir.data <- "/path/to/zenodo/download/"   # Zenodo data folder
dir.out  <- "/path/to/output/"            # where figures/tables will be saved

library(tidyverse)
library(Giotto)
library(igraph)
library(RANN)
library(GSEABase)
library(dorothea)

# STRINGdb and Pathway distance scores -----
stringscores <- read_csv(paste0(dir.data,"KnownInteractionsDatabases_Lung.csv"))[,-1] %>%
  mutate(prop_pathways = 1-pathway_distance,
         prop_pathways = ifelse(is.na(prop_pathways),1,prop_pathways)) %>%
  filter(gene_to != gene_from)
stringscores$biogrid <- ifelse(is.na(stringscores$biogrid), 0, 1)

indepgenes <- stringscores %>%
  filter(pathway_distance == 1,
         stringdp_score == 0,
         biogrid == 0,
         is.na(trrust), is.na(dorothea), is.na(remap), is.na(n_shared_tf))
pair_indepgenes <- paste(indepgenes$gene_from, indepgenes$gene_to, sep = "-")


# Lung 6 - Correlation from all models -----
load(file.path(dir.data, "Lung_Giotto_Filtered_Tumor.RData"))
id_tissue <- pDataDT(gem)$cell_ID[pDataDT(gem)$tissue == "Lung6"]
gem_tissue <- subsetGiotto(gem, cell_ids = id_tissue)

locs <- gem_tissue@spatial_locs$cell$raw@coordinates
rownames(locs) <- locs$cell_ID
normexpr <- as.data.frame(as.matrix(t(gem_tissue@expression$cell$rna$normalized@exprMat)))
normexpr$cell_ID <- rownames(normexpr)
normexpr <- merge(normexpr,locs, by = "cell_ID" )


# SVG
load(paste0(dir.out,"Lung6/km_spatialgenesLung6.RData"))
svgs <- km_spatialgenes %>% filter(adj.p.value < 0.05)


# Contamination
contamination0 <- read.csv(paste0(dir.out,"ContaminationMetricsGenes_Lung.csv"))
contamination0 <- contamination0 %>% filter(tissue == "Lung6", cell_type == "tumor")
contamination <- contamination0 %>% filter(avg_cluster > 0.1, ratio < 1, tissue == "Lung6", cell_type == "tumor")

# Giotto
load(paste0(dir.out,"Lung6/Data_Lung6_giotto.Rdata"))
cor_giotto_long <- cor_giotto$cor_DT %>%
  dplyr::select(feat_ID, variable, spat_cor) %>%
  dplyr::rename(Var1 = feat_ID, Var2 = variable, Giotto=spat_cor)


# Spacedecorr
# k 100
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k100.Rdata"))
edf100 <- data_spacedecorr$edf
summary(as.numeric(edf100))
spacedecorr_k100 <- cbind(data_spacedecorr$residuals,pDataDT(gem_tissue, spat_unit = "cell"))
spacedecorr_spatial_k100 <- cbind(data_spacedecorr$splines,pDataDT(gem_tissue, spat_unit = "cell"),gem_tissue@spatial_locs$cell$raw@coordinates)
cor_spacedecorr_long_k100 <- reshape2::melt(cor(data_spacedecorr$residuals)) %>%
  dplyr::rename(spacedecorr100 = value)
time_spacedecorr_k100 <- time_spacedecorr
# 648.630s - 10.81min


# k 200
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k200.Rdata"))
edf200 <- data_spacedecorr$edf
summary(as.numeric(edf200))

spacedecorr_k200 <- cbind(data_spacedecorr$residuals,pDataDT(gem_tissue, spat_unit = "cell"))
spacedecorr_spatial_k200 <- cbind(data_spacedecorr$splines,pDataDT(gem_tissue, spat_unit = "cell"),gem_tissue@spatial_locs$cell$raw@coordinates)
cor_spacedecorr_long_k200 <- reshape2::melt(cor(data_spacedecorr$residuals)) %>%
  dplyr::rename(spacedecorr200 = value)
time_spacedecorr_k200 <- time_spacedecorr
# 1232.404s - 20.54min


# k 1000
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k1000.Rdata"))
edf1000 <- data_spacedecorr$edf
spacedecorr_k1000 <- cbind(data_spacedecorr$residuals,pDataDT(gem_tissue, spat_unit = "cell"))
spacedecorr_spatial_k1000 <- cbind(data_spacedecorr$splines,pDataDT(gem_tissue, spat_unit = "cell"),gem_tissue@spatial_locs$cell$raw@coordinates)

cor_spacedecorr_long_k1000 <- reshape2::melt(cor(data_spacedecorr$residuals)) %>%
  dplyr::rename(spacedecorr1000 = value)
time_spacedecorr_k1000 <- time_spacedecorr
# 73880.47s - 1231.34min - 20.52h

# k 200 spots
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorrspots_k200.Rdata"))
cor_spacedecorr_long_k200spots <- reshape2::melt(cor(data_spacedecorr$residuals)) %>%
  dplyr::rename(spacedecorr200spots = value)
time_spacedecorr_k200spots <- time_spacedecorr

edf_long <- bind_rows(
  data.frame(gene = colnames(edf100),  edf = as.numeric(edf100),  K = 100),
  data.frame(gene = colnames(edf200),  edf = as.numeric(edf200),  K = 200),
  data.frame(gene = colnames(edf1000), edf = as.numeric(edf1000), K = 1000)
)

write.csv(edf_long, file = file.path(dir.out, "edf_lung.csv"))

#spMOCA
load(paste0(dir.out,"Lung6/Data_Lung6_spmocaspots.Rdata"))
cor_spmoca_long <- reshape2::melt(cor_spmoca) %>%
  dplyr::rename(spmoca = value)

# Pearson
cor_pearson_long <- reshape2::melt(cor(t(as.matrix(gem_tissue@expression$cell$rna$raw@exprMat)))) %>%dplyr::rename(pearson = value)
Ylog <- apply(t(as.matrix(gem_tissue@expression$cell$rna$raw@exprMat)),2, function(x)log(1e6 * x/pDataDT(gem_tissue, spat_unit = "cell")$totalcounts +1))
cor_logpearson_long <- reshape2::melt(cor(Ylog)) %>%dplyr::rename(logpearson = value)
cor_spearman_long <- reshape2::melt(cor(t(as.matrix(gem_tissue@expression$cell$rna$raw@exprMat)), method = "spearman")) %>%dplyr::rename(spearman = value)


# vst
load(paste0(dir.out,"Lung6/Data_Lung6_vst.Rdata"))
cor_vst_long <- reshape2::melt(cor(t(data_vst$y))) %>%
  dplyr::rename(vst = value)

# combining all
cor_all_long <- cbind(cor_spacedecorr_long_k200) %>%
  left_join(cor_spacedecorr_long_k1000) %>%
  left_join(cor_spacedecorr_long_k100) %>%
  left_join(cor_spacedecorr_long_k200spots) %>%
  left_join(cor_spmoca_long) %>%
  left_join(cor_giotto_long) %>%
  left_join(cor_pearson_long) %>% 
  left_join(cor_logpearson_long) %>% 
  left_join(cor_spearman_long) %>%
  left_join(cor_vst_long) %>%
  merge(.,stringscores, by.x = c("Var1","Var2"), by.y = c("gene_from","gene_to"))  %>%
  filter(Var1 != Var2) %>%
  merge(., contamination0[,-1], by.x = "Var1", by.y = "target",all.x = T ) %>%
  merge(., contamination0[,-1], by.x = "Var2", by.y = "target",all.x = T, suffixes = c(".v1",".v2") ) %>%
  #svg
  merge(., svgs[,c(1,5,6)], by.x = "Var1", by.y = "genes",all.x = T ) %>%
  merge(., svgs[,c(1,5,6)], by.x = "Var2", by.y = "genes",all.x = T, suffixes = c(".v1",".v2") ) %>%
  filter(avg_cluster.v1 > .1, avg_cluster.v2 > .1) %>%
  filter(Var1 %in% contamination$target, Var2 %in% contamination$target) %>%
  mutate(independent = ifelse(pathway_distance == 1 &
                                stringdp_score < 100 &
                                biogrid == 0 &
                                is.na(trrust) & is.na(dorothea) & is.na(remap) & is.na(n_shared_tf) , T, F))
write.csv(cor_all_long, file = file.path(dir.out, "cor_all_long_lung.csv"))


# Modified t-test ---------------------------------------------------------
source(file.path(dir.proj, "modified_cor_test.R"))

# Cells, coordinates, genes ---------------------------------------------------
cell_ids <- pDataDT(gem_tissue, spat_unit = "cell")$cell_ID       # row order of the residuals
locs     <- as.data.frame(gem_tissue@spatial_locs$cell$raw@coordinates)
coords   <- as.matrix(locs[match(cell_ids, locs$cell_ID), c("sdimx", "sdimy")])
coords <- coords * 1000 # to um

# coordinate units: typical distance to the nearest cell should be ~10 µm
nn1 <- RANN::nn2(coords, k = 2)$nn.dists[, 2]

genes <- contamination$target                                      # same gene filter as cor_all_long

# Cell x gene matrices for the three methods ---------------------------------
load(paste0(dir.out, "Lung6/Data_Lung6_vst.Rdata"))
R_vst <- t(as.matrix(data_vst$y))                                  # cells x genes
R_vst <- R_vst[cell_ids, ]

R_log <- Ylog; rownames(R_log) <- cell_ids                         # log-normalized, from the loading code

# residuals as cells x genes matrices, rows named by cell ID
R_k <- list()
for (k in c(100, 200, 1000)) {
  d <- as.data.frame(get(paste0("spacedecorr_k", k)))
  M <- as.matrix(d[, genes])
  rownames(M) <- d$cell_ID
  R_k[[paste0("spacedecorr", k)]] <- M
}
mats <- c(R_k, list(vst = R_vst, lognorm = R_log))

# modified t-test per method (same settings as before)
mt_k <- list()
for (m in names(mats)) {
  print(c(m, format(Sys.time())))
  mt_k[[m]] <- modified_cor_test(mats[[m]][cell_ids, genes], coords,
                                 short_max = 200, n_query = 10000, k_max = 3000)
}

res_k <- bind_rows(lapply(mt_k, function(x) x$results), .id = "method") %>%
  group_by(method) %>%
  mutate(fdr_naive    = p.adjust(p_naive, "BH"),
         fdr_modified = p.adjust(p_modified, "BH")) %>%
  ungroup()
save(mt_k, res_k, file = file.path(dir.out, "res_modified_ttest_Lung6.Rdata"))



# WGCNA -------------------------------------------------------------------
## Aux functions ----------
clusterGenes <- function (cor_matrix,k = 10, hclust_method = "ward.D"){
  cor_dist <- stats::as.dist(1 - cor_matrix)
  cor_h <- stats::hclust(d = cor_dist, method = hclust_method)
  cor_clus <- stats::cutree(cor_h, k = k)
  return(list(hclust = cor_h, clusters = cor_clus))
}

getEnrichClust <- function(clusters) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
  universe <- mapIds(
    org.Hs.eg.db,
    keys = names(clusters),
    column = "ENTREZID",
    keytype = "SYMBOL"
  )
  hallmark <- msigdbr(species = "Homo sapiens", category = "H") |>
    dplyr::select(gs_name, gene_symbol)
  
  c8 <- msigdbr(species = "Homo sapiens", category = "C8")
  
  do.call(rbind, lapply(1:max(clusters), function(clustID){
    genes <- universe[clusters == clustID]
    
    # GO
    go_res <- enrichGO(gene          = genes,
                       universe     = universe,
                       OrgDb        = org.Hs.eg.db,
                       ont          = "BP",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
    go_res <- tryCatch(clusterProfiler::simplify(go_res), error = function(e) go_res)
    go_df  <- if (!is.null(go_res) && nrow(go_res@result) > 0) {
      d <- go_res@result
      d$source  <- "GO"
      d$clustID <- clustID
      d
    } else NULL
    
    # Hallmark
    hm_res <- enricher(gene          = genes,
                       universe     = universe,
                       TERM2GENE    = hallmark,
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)
    hm_df  <- if (!is.null(hm_res) && nrow(hm_res@result) > 0) {
      d <- hm_res@result
      d$source  <- "Hallmark"
      d$clustID <- clustID
      d
    } else NULL
    
    
    dplyr::bind_rows(go_df, hm_df)
  }))
}

## Get modules and enrichment for varying k ----
for(k in c(5, 10, 20, 30)){
  print(k)
  #spacedecorr
  print("spacedecorr")
  # k = 100
  load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k100.Rdata"))
  cormat_spacedecorr_k100 <- cor(data_spacedecorr$residuals[, colnames(data_spacedecorr$residuals) %in% contamination$target])
  clust_spacedecorr_k100 <- clusterGenes(cormat_spacedecorr_k100, k = k)
  enrich_clust_spacedecorr_k100 <- getEnrichClust(clust_spacedecorr_k100$clusters)
  enrich_clust_spacedecorr_k100$model <- "spacedecorr k=100"
  
  
  # k = 200
  load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k200.Rdata"))
  cormat_spacedecorr_k200 <- cor(data_spacedecorr$residuals[, colnames(data_spacedecorr$residuals) %in% contamination$target])
  clust_spacedecorr_k200 <- clusterGenes(cormat_spacedecorr_k200, k = k)
  enrich_clust_spacedecorr_k200 <- getEnrichClust(clust_spacedecorr_k200$clusters)
  enrich_clust_spacedecorr_k200$model <- "spacedecorr k=200"
  
  # k = 1000
  load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k1000.Rdata"))
  cormat_spacedecorr_k1000 <- cor(data_spacedecorr$residuals[, colnames(data_spacedecorr$residuals) %in% contamination$target])
  clust_spacedecorr_k1000 <- clusterGenes(cormat_spacedecorr_k1000, k = k)
  enrich_clust_spacedecorr_k1000 <- getEnrichClust(clust_spacedecorr_k1000$clusters)
  enrich_clust_spacedecorr_k1000$model <- "spacedecorr k=1000"
  
  # k = 200 spots
  load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorrspots_k200.Rdata"))
  cormat_spacedecorr_k200spots <- cor(data_spacedecorr$residuals[, colnames(data_spacedecorr$residuals) %in% contamination$target])
  clust_spacedecorr_k200spots <- clusterGenes(cormat_spacedecorr_k200spots, k = k)
  enrich_clust_spacedecorr_k200spots <- getEnrichClust(clust_spacedecorr_k200spots$clusters)
  enrich_clust_spacedecorr_k200spots$model <- "spacedecorr k=200 (spots)"
  
  # spMOCA
  print("spMOCA")
  load(paste0(dir.out,"Lung6/Data_Lung6_spmocaspots.Rdata"))
  cormat_spmoca <- cor_spmoca[rownames(cor_spmoca) %in% contamination$target,
                              colnames(cor_spmoca) %in% contamination$target]
  clust_spmoca <- clusterGenes(cormat_spmoca, k = k)
  enrich_clust_spmoca <- getEnrichClust(clust_spmoca$clusters)
  enrich_clust_spmoca$model <- "spmoca"
  
  # Giotto
  print("Giotto")
  cormat_giotto <- cor_giotto$cor_DT %>%
    filter(feat_ID %in% contamination$target,
           variable %in% contamination$target) %>%
    graph_from_data_frame(., directed = F) %>%
    get.adjacency(., attr="spat_cor", sparse=FALSE)
  
  clust_giotto <- clusterGenes(cormat_giotto, k = k)
  enrich_clust_giotto <- getEnrichClust(clust_giotto$clusters)
  enrich_clust_giotto$model <- "giotto"
  
  # vst
  print("vst")
  cormat_vst <- cor(t(data_vst$y[rownames(data_vst$y)%in% contamination$target,]))
  clust_vst <- clusterGenes(cormat_vst, k = k)
  enrich_clust_vst <- getEnrichClust(clust_vst$clusters)
  enrich_clust_vst$model <- "vst"
  
  
  # log CPM
  print("log CPM")
  cormat_logpearson <- cor(Ylog[,colnames(Ylog) %in% contamination$target])
  clust_logpearson <- clusterGenes(cormat_logpearson, k = k)
  enrich_clust_logpearson <- getEnrichClust(clust_logpearson$clusters)
  enrich_clust_logpearson$model <- "logPearson"
  
  save(cormat_spacedecorr_k100,clust_spacedecorr_k100,enrich_clust_spacedecorr_k100,
       cormat_spacedecorr_k200,clust_spacedecorr_k200,enrich_clust_spacedecorr_k200,
       cormat_spacedecorr_k1000,clust_spacedecorr_k1000,enrich_clust_spacedecorr_k1000,
       cormat_spacedecorr_k200spots,clust_spacedecorr_k200spots,enrich_clust_spacedecorr_k200spots,
       cormat_giotto,clust_giotto,enrich_clust_giotto,
       cormat_vst,clust_vst,enrich_clust_vst,
       cormat_logpearson,clust_logpearson,enrich_clust_logpearson,
       cormat_spmoca,clust_spmoca,enrich_clust_spmoca,
       file = file.path(dir.out, paste0("ExampleWGCNA_k",k,".Rdata")))
}
load(paste0(dir.out,"ExampleWGCNA_k20.Rdata"))


# Biological coherence ----------------------------------------------------
cor_all_long <- cor_all_long %>%
  mutate(physical = ifelse(biogrid == 1 | stringdp_score > 700, 1, 0),
         regulatory = ifelse(!is.na(trrust) | (!is.na(dorothea) & confidence %in% c("A","B","C")) | 
                               !is.na(remap) | !is.na(shared_tf), 1,0),
         any = ifelse(physical == 1 | regulatory == 1,1,0))


bio_coherence_modules <- do.call(rbind, lapply(c(5,10,20,30), function(k){
  load(file.path(dir.out, paste0("ExampleWGCNA_k",k,".Rdata")))
  do.call(rbind, lapply(1:k, function(clustid){
    df_sd_k100 <- cor_all_long %>% filter(Var1 %in% names(clust_spacedecorr_k100$clusters[clust_spacedecorr_k100$clusters == clustid]),
                                          Var2 %in% names(clust_spacedecorr_k100$clusters[clust_spacedecorr_k100$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)
      ) %>%
      mutate(cluster = clustid,
             model = "spacedecorrk100")
    df_sd_k200 <- cor_all_long %>% filter(Var1 %in% names(clust_spacedecorr_k200$clusters[clust_spacedecorr_k200$clusters == clustid]),
                                          Var2 %in% names(clust_spacedecorr_k200$clusters[clust_spacedecorr_k200$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)
      ) %>%
      mutate(cluster = clustid,
             model = "spacedecorrk200")
    
    df_sd_k1000 <- cor_all_long %>% filter(Var1 %in% names(clust_spacedecorr_k1000$clusters[clust_spacedecorr_k1000$clusters == clustid]),
                                           Var2 %in% names(clust_spacedecorr_k1000$clusters[clust_spacedecorr_k1000$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)
      ) %>%
      mutate(cluster = clustid,
             model = "spacedecorrk1000")
    
    df_sd_k200spots <- cor_all_long %>% filter(Var1 %in% names(clust_spacedecorr_k200spots$clusters[clust_spacedecorr_k200spots$clusters == clustid]),
                                               Var2 %in% names(clust_spacedecorr_k200spots$clusters[clust_spacedecorr_k200spots$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)
      ) %>%
      mutate(cluster = clustid,
             model = "spacedecorrk200spots")
    
    df_spmoca <- cor_all_long %>% filter(Var1 %in% names(clust_spmoca$clusters[clust_spmoca$clusters == clustid]),
                                         Var2 %in% names(clust_spmoca$clusters[clust_spmoca$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)) %>%
      mutate(cluster = clustid,
             model = "spMOCA")
    
    df_giotto <- cor_all_long %>% filter(Var1 %in% names(clust_giotto$clusters[clust_giotto$clusters == clustid]),
                                         Var2 %in% names(clust_giotto$clusters[clust_giotto$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)) %>%
      mutate(cluster = clustid,
             model = "Giotto")
    
    df_vst <- cor_all_long %>% filter(Var1 %in% names(clust_vst$clusters[clust_vst$clusters == clustid]),
                                      Var2 %in% names(clust_vst$clusters[clust_vst$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)) %>%
      mutate(cluster = clustid,
             model = "SPARK")
    
    df_logpearson <- cor_all_long %>% filter(Var1 %in% names(clust_logpearson$clusters[clust_logpearson$clusters == clustid]),
                                             Var2 %in% names(clust_logpearson$clusters[clust_logpearson$clusters == clustid])) %>%
      summarise(pathways = mean(prop_pathways),
                biogrid = mean(biogrid),
                physical = mean(physical),
                regulatory = mean(regulatory),
                any = mean(any)) %>%
      mutate(cluster = clustid,
             model = "logPearson")
    
    rbind(df_sd_k100,df_sd_k200,df_sd_k1000, df_sd_k200spots,
          df_spmoca,df_giotto,df_vst,df_logpearson) %>% 
      mutate(k = k)
  }))
}))
write.csv(bio_coherence_modules, file = file.path(dir.out, "bio_coherence_modules.csv"))


# Regulator enrichment ----------------------------------------------------
load(paste0(dir.out,"ExampleWGCNA_k20.Rdata"))

# Modules of each gene 
modules_df <- bind_rows(
  data.frame(model = "spacedecorr", gene = names(clust_spacedecorr_k200$clusters),
             module = clust_spacedecorr_k200$clusters),
  data.frame(model = "giotto",  gene = names(clust_giotto$clusters),
             module = clust_giotto$clusters),
  data.frame(model = "spark",   gene = names(clust_vst$clusters),
             module = clust_vst$clusters),
  data.frame(model = "spmoca",  gene = names(clust_spmoca$clusters),
             module = clust_spmoca$clusters),
  data.frame(model = "log-CPM", gene = names(clust_logpearson$clusters),
             module = clust_logpearson$clusters)
)
panel_genes <- unique(modules_df$gene)
num_panel <- length(panel_genes)

## TF databases ----
trrust <- read.delim(file.path(dir.data, "trrust_rawdata.human.tsv"), header = FALSE,
                     col.names = c("tf", "target", "mode", "pmid"))

doro <- dorothea::dorothea_hs
doro <- doro[doro$confidence %in% c("A", "B", "C"), ]

remap_gmt <- getGmt(file.path(dir.data, "ReMap_ChIP-seq.gmt"))
remap <- bind_rows(lapply(remap_gmt, function(s)
  data.frame(tf = setName(s), target = geneIds(s))))

regulons <- bind_rows(
  data.frame(tf = trrust$tf, target = trrust$target, source = "TRRUST"),
  data.frame(tf = doro$tf,   target = doro$target,   source = "DoRothEA"),
  data.frame(tf = remap$tf,  target = remap$target,  source = "ReMap")
)

## At least 3 targets on the panel; at most half the panel
regulons <- regulons %>%
  filter(target %in% panel_genes) %>%
  distinct() %>%
  group_by(source, tf) %>%
  filter(n() >= 3, n() <= 0.5 * num_panel) %>%
  ungroup()

count(regulons, source, name = "n_pairs")
regulons %>% distinct(source, tf) %>% count(source, name = "n_tfs")

## Enrichment test -------------------------
modules <- split(modules_df$gene, paste(modules_df$model, modules_df$module, sep = "|"))
modules <- modules[lengths(modules) >= 3]
targets <- split(regulons$target, paste(regulons$source, regulons$tf, sep = "|"))

tf_enrich <- bind_rows(lapply(names(modules), function(m) {
  genes <- modules[[m]]
  data.frame(
    module_id         = m,
    reg_id            = names(targets),
    num_module        = length(genes),
    num_panel_targets = lengths(targets),
    num_mod_targets   = sapply(targets, function(t) sum(genes %in% t))
  )
})) %>%
  filter(num_mod_targets > 0) %>%
  mutate(pvalue = phyper(num_mod_targets - 1, 
                         num_panel_targets, 
                         num_panel - num_panel_targets, num_module, lower.tail = FALSE)) %>%
  tidyr::separate(module_id, c("model", "module"), sep = "\\|") %>%
  tidyr::separate(reg_id, c("source", "tf"), sep = "\\|") %>%
  group_by(model, module, source) %>%
  mutate(fdr = p.adjust(pvalue, method = "BH")) %>% # fdr by module by model
  ungroup() %>%
  arrange(model, module, source, pvalue)


module_summary <- tf_enrich %>%
  group_by(model, module, source, size = num_module) %>%
  summarise(n_sig_tf = sum(fdr < 0.05),
            top_tf   = tf[which.min(pvalue)],
            top_fdr  = min(fdr),
            .groups  = "drop") %>%
  mutate(any_sig = n_sig_tf > 0)

module_summary %>% filter(any_sig) %>% arrange(model, module)

tf_enrich %>% filter(fdr < 0.05)

write.csv(tf_enrich,  file.path(dir.out,"regulator_enrichment.csv"), row.names = FALSE)


# Module category ---------------------------------------------------------
# All cells of Lung6 (every cell type), with coordinates ---
load(file.path(dir.data, "Lung_Giotto.RData"))
new_gem <- createGiottoObject(
  expression    = list(raw = gem@expression$rna$raw),
  spatial_locs  = gem@spatial_locs$raw,
  cell_metadata = gem@cell_metadata$rna
)

meta <- as.data.frame(pDataDT(new_gem, spat_unit = "cell"))
locs <- as.data.frame(new_gem@spatial_locs$cell$raw@coordinates)
cells_all <- merge(meta, locs[, c("cell_ID", "sdimx", "sdimy")], by = "cell_ID") %>%
  filter(tissue == "Lung6")

group_map <- c(
  "tumor 5" = "tumor", "tumor 6" = "tumor", "tumor 9" = "tumor",
  "tumor 12" = "tumor", "tumor 13" = "tumor",
  "epithelial" = "epithelial", "fibroblast" = "fibroblast",
  "endothelial" = "endothelial",
  "macrophage" = "myeloid", "monocyte" = "myeloid", "neutrophil" = "myeloid",
  "mast" = "myeloid", "mDC" = "myeloid", "pDC" = "myeloid",
  "T CD4 naive" = "T/NK", "T CD4 memory" = "T/NK", "T CD8 naive" = "T/NK",
  "T CD8 memory" = "T/NK", "Treg" = "T/NK", "NK" = "T/NK",
  "B-cell" = "B/plasma", "plasmablast" = "B/plasma"
)
cells_all$group <- unname(group_map[cells_all$cell_type])
table(cells_all$group)

## SpaceDecorr residuals ---
load(file.path(dir.out, "Lung6", "Data_Lung6_spacedecorr_k200.Rdata"))
resid <- t(data_spacedecorr$residuals)
colnames(resid) <- rownames(data_spacedecorr$metadata)

## keep tumor cells that have coordinates, in the same order everywhere
tumor_ids <- intersect(colnames(resid), cells_all$cell_ID)
resid <- resid[, tumor_ids]
length(tumor_ids)

## Neighbourhood composition of each tumor cell --------------------------------
## For each tumor cell: proportion of each cell-type group among its k nearest
kneigh <- 30
xy <- as.matrix(cells_all[, c("sdimx", "sdimy")])
focal <- match(tumor_ids, cells_all$cell_ID)
nn <- nn2(xy, xy[focal, ], k = kneigh + 1)
nb_groups <- matrix(cells_all$group[nn$nn.idx[, -1]], nrow = length(focal))

groups <- setdiff(sort(unique(cells_all$group)), "tumor")   # tumor = reference
nb <- sapply(groups, function(g) rowMeans(nb_groups == g))
colnames(nb) <- make.names(groups)

## Log-normalized expression for the same genes and tumor cells
load(file.path(dir.data, "Lung_Giotto_Filtered_Tumor.RData"))
logexpr <- as.matrix(gem@expression$cell$rna$normalized@exprMat[rownames(resid), tumor_ids])
Z <- t(scale(t(logexpr)))
Z[is.na(Z)] <- 0

## Niche and FOV per tumor cell
cell_meta <- cells_all[match(tumor_ids, cells_all$cell_ID), c("niche", "fov")]
cell_meta$fov <- factor(cell_meta$fov)

# Gene-gene correlations before and after adjustment ---
load(file.path(dir.out, "Lung6", "Data_Lung6_vst.Rdata"))
cor_unadj <- cor(t(data_vst$y))                    # VST, no spatial adjustment
cor_adj   <- cor(data_spacedecorr$residuals)       # SpaceDecorr

## Results ---
mod_cat <- data.frame()
for (m in unique(modules_df$model)) {
  for (mod in unique(modules_df$module[modules_df$model == m])) {
    
    genes <- modules_df$gene[modules_df$model == m & modules_df$module == mod]
    genes <- intersect(genes, rownames(resid))
    if (length(genes) < 3) next
    
    ## mean pairwise correlation before / after adjustment
    C_un <- cor_unadj[genes, genes]
    C_ad <- cor_adj[genes, genes]
    c_unadj <- mean(C_un[upper.tri(C_un)])
    c_adj   <- mean(C_ad[upper.tri(C_ad)])
    
    ## neighbourhood share
    R    <- resid[genes, ]
    X    <- cbind(1, nb)
    R_nb <- t(lm.fit(X, t(R))$residuals)          # residuals after removing neighbourhood composition
    c_res <- mean(cor(t(R))[upper.tri(diag(nrow(R)))])
    c_nb  <- mean(cor(t(R_nb))[upper.tri(diag(nrow(R)))])
    nb_share <- ifelse(c_res > 0, 1 - c_nb / c_res, NA)
    
    
    ## FOV share: variation explained by FOV beyond the niche,
    ## using log expression (keeps FOV effects in)
    score_log <- colMeans(Z[genes, ])
    r2_int  <- summary(lm(score_log ~ 1, data = cell_meta))$r.squared
    r2_fov    <- summary(lm(score_log ~ 1 + fov, data = cell_meta))$r.squared
    fov_share <- r2_fov - r2_int
    
    mod_cat <- rbind(mod_cat, data.frame(
      model = m, module = mod, size = length(genes),
      c_unadj = c_unadj, c_adj = c_adj,
      broad_share = ifelse(c_unadj > 0, 1 - c_adj / c_unadj, NA),
      nb_share = nb_share,
      fov_share = fov_share))
  }
}

write.csv(mod_cat, file.path(dir.out, "module_decomposition.csv"), row.names = FALSE)


# Supplementary table -----------------------------------------------------
library(tidyverse)
dir.out <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real data/output/Lung/"

load(paste0(dir.out, "ExampleWGCNA_k20.Rdata"))
n_modules <- 20

# Labels of each method in the different objects and files
# (model is the name shown in the table, as in the paper)
methods <- tribble(
  ~model,           ~clust,                   ~cormat,                   ~go_lab,             ~reg_lab,
  "SpaceDecorr",    "clust_spacedecorr_k200", "cormat_spacedecorr_k200", "spacedecorr k=200", "spacedecorr",
  "Giotto",         "clust_giotto",           "cormat_giotto",           "giotto",            "giotto",
  "VST-normalized", "clust_vst",              "cormat_vst",              "vst",               "spark",
  "spMOCA",         "clust_spmoca",           "cormat_spmoca",           "spmoca",            "spmoca",
  "log-normalized", "clust_logpearson",       "cormat_logpearson",       "logPearson",        "log-CPM"
)

# Results shared by all methods
reg <- read.csv(file.path(dir.out, "regulator_enrichment.csv"))
dec_all <- read.csv(file.path(dir.out, "module_decomposition.csv"))

enrich_df <- do.call(rbind, lapply(c(5,10,20,30), function(k){
  load(file.path(dir.out, paste0("ExampleWGCNA_k",k,".Rdata")))
  rbind(enrich_clust_spacedecorr_k100,
        enrich_clust_spacedecorr_k200,
        enrich_clust_spacedecorr_k1000,
        enrich_clust_spacedecorr_k200spots,
        enrich_clust_spmoca,
        enrich_clust_giotto,
        enrich_clust_logpearson,
        enrich_clust_vst
  ) %>% mutate(k=k)
}))

# -----------------------------------------------------------------------------
# Build the table for each method
# -----------------------------------------------------------------------------
module_table <- data.frame()
for (i in seq_len(nrow(methods))) {
  
  clust  <- get(methods$clust[i])
  cormat <- get(methods$cormat[i])
  
  # 1. Size, hub genes (highest intramodular connectivity), full gene list
  modules <- data.frame(gene = names(clust$clusters),
                        module = as.character(clust$clusters))
  hubs <- modules %>%
    group_by(module) %>%
    summarise(size = n(),
              hub_genes = {
                k_in <- rowSums(cormat[gene, gene, drop = FALSE]) - 1
                paste(head(gene[order(k_in, decreasing = TRUE)], 5), collapse = ", ")
              },
              all_genes = paste(sort(gene), collapse = ", "),
              .groups = "drop")
  
  # 2. Top GO terms (up to ten, BH-adjusted p < 0.05, after simplify)
  go_top <- enrich_df %>%
    filter(k == n_modules, model == methods$go_lab[i], source == "GO",
           p.adjust < 0.05) %>%
    mutate(module = as.character(clustID)) %>%
    arrange(module, p.adjust) %>%
    group_by(module) %>%
    summarise(n_go_terms   = n(),
              top_go_terms = paste(head(Description, 10), collapse = "; "),
              top_go_fdr   = min(p.adjust),
              .groups = "drop")
  
  # 3. Enriched regulators (FDR < 0.05); T = TRRUST, D = DoRothEA, R = ReMap
  reg_top <- reg %>%
    filter(model == methods$reg_lab[i], fdr < 0.05) %>%
    mutate(module = as.character(module),
           db = recode(source, TRRUST = "T", DoRothEA = "D", ReMap = "R")) %>%
    group_by(module, tf) %>%
    summarise(fdr = min(fdr), db = paste(sort(unique(db)), collapse = ","),
              .groups = "drop") %>%
    arrange(module, fdr) %>%
    group_by(module) %>%
    summarise(regulators = paste0(head(tf, 5), " (", head(db, 5), ")", collapse = ", "),
              .groups = "drop")
  
  # 4. Decomposition values (modules with < 3 genes are blank)
  dec <- dec_all %>%
    filter(model == methods$reg_lab[i]) %>%
    mutate(module = as.character(module)) %>%
    dplyr::select(module, c_adj, broad_share, nb_share)
  
  # 5. Combine
  tab <- hubs %>%
    left_join(go_top,  by = "module") %>%
    left_join(reg_top, by = "module") %>%
    left_join(dec,     by = "module") %>%
    mutate(model = methods$model[i],
           module = as.integer(module),
           n_go_terms = replace_na(n_go_terms, 0L)) %>%
    arrange(module)
  
  module_table <- bind_rows(module_table, tab)
}

module_table <- module_table %>%
  mutate(model = factor(model, levels = methods$model),
         top_go_fdr = signif(top_go_fdr, 2),
         across(c(c_adj, broad_share, nb_share), ~ round(.x, 2))) %>%
  arrange(model, module) %>%
  relocate(model)

module_table %>% count(model)

# Full table (with every gene) as CSV, NA saved as blank
write.csv(module_table, file.path(dir.out, "Supp_module_annotation_all.csv"),
          row.names = FALSE, na = "")