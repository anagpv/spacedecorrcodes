#-------------------------------------#
#- Differential corr Lupus Real data -#
#-------------------------------------#
# https://htmlpreview.github.io/?https://github.com/andymckenzie/DGCA/blob/master/inst/doc/DGCA.html

# Lupus data must be downloaded from: https://doi.org/10.6084/m9.figshare.c.7373860
dir.data <- "/path/to/figshare/download/"   # where cleaneddata.RData lives
dir.out  <- "/path/to/output/"              # where intermediate .Rdata files were saved by get_spacedecorr_lupus.R
dir.proj <- "/path/to/repo/"               # root of this GitHub repo


library(DGCA)
library(tidyverse)
library(tidygraph)
library(GOstats)
library(HGNChelper)
library(org.Hs.eg.db)
library(clusterProfiler)
library(igraph)
library(ggraph)
library(ggnet)
library(viridis)
library(ggrepel)
library(cowplot)

# Reading data
load(file.path(dir.data, "cleaneddata.RData"))
rownames(customlocs) = annot$cell_ID
contam <- read.csv(file.path(dir.out, "ContaminationMetricsGenes_Lupus.csv"))
# filtering average count 0.1 per cell
id_PCT <- annot$cell_ID[which(annot$clusts == "PCT")]
raw_PCT <- as.matrix(raw)[,id_PCT]
genes_PCT <- contam %>% filter(clusts == "PCT", ratio <  1,
                               target %in% rownames(raw_PCT)[(rowMeans(raw_PCT) > 0.1)]) %>% pull(target)


tissues <- sort(unique(annot$tissuename), decreasing = T)
tissues <- tissues[!tissues %in% c("SLE8.2","SLE8.3")]
rm(list = c("raw","annot","raw_PCT"))


# STRING and BioGRID ------------------------------------------------------
load(file.path(dir.out, "Lung_Giotto_Filtered_Tumor.RData"))
stringscores <- read_csv(file.path(dir.data, "KnownInteractionsDatabases_Lupus.csv"))[,-1] %>%
  mutate(prop_pathways = 1-pathway_distance,
         prop_pathways = ifelse(is.na(prop_pathways),1,prop_pathways)) %>%
  filter(gene_to != gene_from)
stringscores$stringdp_factor <- cut(stringscores$stringdp_score,
                                    breaks = c(0,1, 400, 700,900, 1000),
                                    labels = c("No known interaction","Low ", "Medium", "High", "Very high"),
                                    include.lowest = TRUE)
stringscores$proppath_factor <- cut(stringscores$prop_pathways,
                                    breaks = c(0,0.0001, 0.01, 0.1, 0.25, 1),
                                    labels = c("None","Less than 1%","1% to 10%", "10% to 25%", "More than 25%"),
                                    include.lowest = TRUE)
stringscores$biogrid <- ifelse(is.na(stringscores$biogrid), 0, 1)
stringscores <- stringscores %>%
  mutate(indep = ifelse(pathway_distance == 1 &
                          stringdp_score < 100 &
                          biogrid == 0 &
                          is.na(trrust) & is.na(dorothea) & is.na(remap), T, F))

indepgenes <- stringscores %>%
  filter(pathway_distance == 1,
         stringdp_score < 100,
         biogrid == 0,
         is.na(trrust), is.na(dorothea), is.na(remap))

indepgenes <- indepgenes %>% rowwise %>%
  mutate(pair = paste(min(gene_from, gene_to),
                      max(gene_from, gene_to),
                      sep = "-"))
pair_indepgenes <- paste(indepgenes$gene_from, indepgenes$gene_to, sep = "-")



# SVG ---------------------------------------------------------------------
svg_pct0 <- do.call(plyr::rbind.fill, lapply(tissues, function(tissue){
  load(file.path(dir.out, tissue, paste0("km_spat_SVG_pct_",tissue,".Rdata")))

  df <- km_spat_pct %>% as.data.frame %>% mutate(tissue = tissue)
  df
}))

svg_pct <- svg_pct0 %>%
  group_by(genes) %>%
  summarise(minpvalue = min(adj.p.value),
            cctpvalue = cct(adj.p.value),
            score = mean(score),
            av_expr = mean(av_expr))

svg_podocyte %>% arrange(genes %in% c("IFIH1","S100A9")) %>%
  ggplot(aes(x = score, y = av_expr, col = (genes %in% c("IFIH1","S100A9")))) + geom_point()

# PCT ----
## Getting results ----
### giotto ----
giotto_lupus_pct <- do.call(plyr::rbind.fill, lapply(tissues, function(tissue){
  load(file.path(dir.out, tissue, paste0("Data_PCT_",tissue,"_Giotto.Rdata")))

  df <- data_giotto %>% as.data.frame %>% mutate(tissue = tissue)
  df$cell_ID <- rownames(df)
  df
}))

giotto_lupus_pct <- giotto_lupus_pct %>%
  arrange(cell_ID) %>%
  mutate(SLE = ifelse(grepl("SLE", tissue),"1","0" ))
rownames(giotto_lupus_pct) <- giotto_lupus_pct$cell_ID
giotto_lupus_pct<-giotto_lupus_pct[,colnames(giotto_lupus_pct) %in% c(genes_PCT,"cell_ID","SLE")]

expmat_pct_giotto <- as.data.frame(t(giotto_lupus_pct[,!colnames(giotto_lupus_pct)%in%c("tissue","cell_ID","SLE")]))
expmat_pct_giotto <- apply(expmat_pct_giotto,2,as.numeric)
design_pct_giotto <- model.matrix(~SLE-1, data = giotto_lupus_pct)
rownames(expmat_pct_giotto) <- colnames(giotto_lupus_pct)[!colnames(giotto_lupus_pct)%in%c("tissue","cell_ID","SLE")]

# Get differential expressed genes --
ddcor_pct_pearson_giotto = ddcorAll(inputMat = as.data.frame(expmat_pct_giotto), # p x n, gene x subject
                                    design = design_pct_giotto, # indicator matrix of groups/condition/celltype
                                    compare = c("SLE0", "SLE1"), # two columns of design to compare
                                    #adjust = "BH",
                                    heatmapPlot = F,
                                    #corrType = "spearman", #slower
                                    nPerm = 10, # number of permutations for test
                                    nPairs = "all", # number of top pairs to show
                                    classify = T
)
rm(list = c("giotto_lupus_pct","expmat_pct_giotto","design_pct_giotto"))

ddcor_pct_pearson_giotto <- ddcor_pct_pearson_giotto %>%
  mutate(SLE0_pValadj = p.adjust(SLE0_pVal, method = "BH"),
         SLE1_pValadj = p.adjust(SLE0_pVal, method = "BH"),
         classSLE0 = as.numeric(SLE0_pValadj < 0.05 & abs(SLE0_cor) > 0.1) * sign(SLE0_cor),
         classSLE0 = case_when(classSLE0 == -1 ~ "-",
                               classSLE0 == 0 ~ "0",
                               classSLE0 == 1 ~ "+"),
         classSLE1 = as.numeric(SLE1_pValadj < 0.05 & abs(SLE1_cor) > 0.1) * sign(SLE1_cor),
         classSLE1 = case_when(classSLE1 == -1 ~ "-",
                               classSLE1 == 0 ~ "0",
                               classSLE1 == 1 ~ "+"),
         class = paste(classSLE0,classSLE1, sep = "/"))

table(ddcor_pct_pearson_giotto$pValDiff_adj<0.05)
ddcor_pct_pearson_giotto %>% filter(pValDiff_adj<0.05) %>%
  mutate(diffCor = SLE1_cor - SLE0_cor) %>%
  filter(abs(diffCor) > 0.1) %>% nrow
ddcor_pct_pearson_giotto %>% filter(pValDiff_adj<0.05) %>%
  mutate(diffCor = SLE1_cor - SLE0_cor) %>%
  filter(abs(diffCor) > 0.1) %>% count(Classes)


### vst ----
vst_lupus_pct <- do.call(plyr::rbind.fill, lapply(tissues, function(tissue){
  load(file.path(dir.out, tissue, paste0("Data_PCT_",tissue,"_vst.Rdata")))

  df <- t(data_vst$y) %>% as.data.frame %>% mutate(tissue = tissue)
  df$cell_ID <- rownames(df)
  df

}))

vst_lupus_pct <- vst_lupus_pct %>%
  arrange(cell_ID) %>%
  mutate(SLE = ifelse(grepl("SLE", tissue),"1","0" ))
rownames(vst_lupus_pct) <- vst_lupus_pct$cell_ID
vst_lupus_pct<-vst_lupus_pct[,colnames(vst_lupus_pct) %in% c(genes_PCT,"cell_ID","SLE")]

expmat_pct_vst <- as.data.frame(t(vst_lupus_pct[,!colnames(vst_lupus_pct)%in%c("tissue","cell_ID","SLE")]))
expmat_pct_vst <- apply(expmat_pct_vst,2,as.numeric)
design_pct_vst <- model.matrix(~SLE-1, data = vst_lupus_pct)
rownames(expmat_pct_vst) <- colnames(vst_lupus_pct)[!colnames(vst_lupus_pct)%in%c("tissue","cell_ID","SLE")]

# Get differential expressed genes --
ddcor_pct_pearson_vst = ddcorAll(inputMat = as.data.frame(expmat_pct_vst), # p x n, gene x subject
                                 design = design_pct_vst, # indicator matrix of groups/condition/celltype
                                 compare = c("SLE0", "SLE1"), # two columns of design to compare
                                 #adjust = "BH",
                                 heatmapPlot = F,
                                 #corrType = "spearman", #slower
                                 nPerm = 10, # number of permutations for test
                                 nPairs = "all", # number of top pairs to show
                                 classify = T
)
rm(list = c("vst_lupus_pct","expmat_pct_vst","design_pct_vst"))

ddcor_pct_pearson_vst <- ddcor_pct_pearson_vst %>%
  mutate(SLE0_pValadj = p.adjust(SLE0_pVal, method = "BH"),
         SLE1_pValadj = p.adjust(SLE0_pVal, method = "BH"),
         classSLE0 = as.numeric(SLE0_pValadj < 0.05 & abs(SLE0_cor) > 0.1) * sign(SLE0_cor),
         classSLE0 = case_when(classSLE0 == -1 ~ "-",
                               classSLE0 == 0 ~ "0",
                               classSLE0 == 1 ~ "+"),
         classSLE1 = as.numeric(SLE1_pValadj < 0.05 & abs(SLE1_cor) > 0.1) * sign(SLE1_cor),
         classSLE1 = case_when(classSLE1 == -1 ~ "-",
                               classSLE1 == 0 ~ "0",
                               classSLE1 == 1 ~ "+"),
         class = paste(classSLE0,classSLE1, sep = "/"))

table(ddcor_pct_pearson_vst$pValDiff_adj<0.05)
table(ddcor_pct_pearson_vst$Classes)
ddcor_pct_pearson_vst %>% filter(pValDiff_adj<0.05) %>%
  mutate(diffCor = SLE1_cor - SLE0_cor) %>%
  filter(abs(diffCor) > 0.1) %>% nrow
ddcor_pct_pearson_vst %>% filter(pValDiff_adj<0.05) %>%
  mutate(diffCor = SLE1_cor - SLE0_cor) %>%
  filter(abs(diffCor) > 0.1) %>% count(Classes)


### Spacedecorr ----
residuals_lupus_pct <- do.call(rbind, lapply(c("SLE8.1","SLE7","SLE3","SLE2","SLE1","Control4","Control3", "Control2", "Control1"), function(tissue){
  load(file.path(dir.out, tissue, paste0("Data_PCT_",tissue,"_spacedecorr_k200_basists.Rdata")))

  df <- data_spacedecorr$Residuals %>% as.data.frame %>% mutate(tissue = tissue)
  df$cell_ID <- rownames(data_spacedecorr$metadata)
  df
}))


ids_pct <- residuals_lupus_pct$cell_ID[!residuals_lupus_pct$tissue %in% c("SLE8.2","SLE8.3")]
residuals_lupus_pct <- residuals_lupus_pct %>%
  filter(cell_ID %in% ids_pct) %>%
  arrange(cell_ID) %>%
  mutate(SLE = ifelse(grepl("SLE", tissue),"1","0" ))
rownames(residuals_lupus_pct) <- residuals_lupus_pct$cell_ID
residuals_lupus_pct<-residuals_lupus_pct[,colnames(residuals_lupus_pct) %in% c(genes_PCT,"cell_ID","SLE")]

expmat_pct <- as.data.frame(t(residuals_lupus_pct[,!colnames(residuals_lupus_pct)%in%c("tissue","cell_ID","SLE")]))
design_pct <- model.matrix(~SLE-1, data = residuals_lupus_pct)

ddcor_pct_pearson = ddcorAll(inputMat = expmat_pct, # p x n, gene x subject
                             design = design_pct, # indicator matrix of groups/condition/celltype
                             compare = c("SLE0", "SLE1"), # two columns of design to compare
                             heatmapPlot = F,
                             nPerm = 10, # number of permutations for test
                             nPairs = "all", # number of top pairs to show
                             classify = T
)
rm(list = c("residuals_lupus_pct","expmat_pct","design_pct"))

ddcor_pct_pearson <- ddcor_pct_pearson %>%
  mutate(SLE0_pValadj = p.adjust(SLE0_pVal, method = "BH"),
         SLE1_pValadj = p.adjust(SLE0_pVal, method = "BH"),
         classSLE0 = as.numeric(SLE0_pValadj < 0.05 & abs(SLE0_cor) > 0.1) * sign(SLE0_cor),
         classSLE0 = case_when(classSLE0 == -1 ~ "-",
                               classSLE0 == 0 ~ "0",
                               classSLE0 == 1 ~ "+"),
         classSLE1 = as.numeric(SLE1_pValadj < 0.05 & abs(SLE1_cor) > 0.1) * sign(SLE1_cor),
         classSLE1 = case_when(classSLE1 == -1 ~ "-",
                               classSLE1 == 0 ~ "0",
                               classSLE1 == 1 ~ "+"),
         class = paste(classSLE0,classSLE1, sep = "/"))

table(ddcor_pct_pearson$pValDiff_adj<0.05)
table(ddcor_pct_pearson$Classes)
ddcor_pct_pearson %>% filter(pValDiff_adj<0.05) %>%
  mutate(diffCor = SLE1_cor - SLE0_cor) %>%
  filter(abs(diffCor) > 0.1) %>% nrow
ddcor_pct_pearson %>% filter(pValDiff_adj<0.05) %>%
  mutate(diffCor = SLE1_cor - SLE0_cor) %>%
  filter(abs(diffCor) > 0.1) %>% count(diffCor>0)


## Scatterplot estimates pvalue ----
### vst ----
all_pearson_pct <- full_join(ddcor_pct_pearson,ddcor_pct_pearson_vst, by = c("Gene1","Gene2"), suffix = c(".sd",".vst") )
all_pearson_pct <- all_pearson_pct %>% merge(.,svg_pct, by.x = "Gene1", by.y = "genes", all.x = T) %>%
  merge(.,svg_pct, by.x = "Gene2", by.y = "genes", all.x = T) %>%
  mutate(
    minpvalue = pmin(minpvalue.x, minpvalue.y, na.rm = TRUE),
    score = rowMeans(cbind(as.numeric(score.x), as.numeric(score.y)), na.rm = TRUE),
    av_expr = rowMeans(cbind(as.numeric(av_expr.x), as.numeric(av_expr.y)), na.rm = TRUE),
    cctpvalue = pmap_dbl(list(cctpvalue.x, cctpvalue.y), ~ cct(c(...)))
  )
all_pearson_pct <- all_pearson_pct %>%
  rowwise %>%
  mutate(pair = paste(min(Gene1, Gene2),
                      max(Gene1, Gene2),
                      sep = "-"))
all_pearson_pct <- all_pearson_pct %>% left_join(indepgenes)

comp_pvalue_pct <- all_pearson_pct %>%arrange(score) %>%
  ggplot(aes(x = -log10(pValDiff_adj.sd), y = -log10(pValDiff_adj.vst),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = -log10(0.05)) +
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "-log10(pvalue SpaceDecorr)",
       y = "-log10(pvalue SPARK)",
       col = "Average spatial autocorrelation score")+
  theme_bw()+
  theme(legend.position = "top",legend.box = "horizontal",
        legend.key.width = unit(1, "cm"),
        text = element_text(size = 15),
        legend.title = element_text(margin = margin(r = 20, b=10)),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.spacing = unit(1, "lines"),              # spacing between keys
        legend.margin = margin(t = 10, b = 10),          # margin around the whole legend
        legend.key.size = unit(1.5, "lines"),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())


comp_corSLE0_pct <- all_pearson_pct %>% arrange(score) %>%
  filter(indep) %>%
  ggplot(aes(x = (SLE0_cor.sd), y = (SLE0_cor.vst),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline() +
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "Correlation for controls\nwith SpaceDecorr",
       y = "Correlation for controls\nwith SPARK",
       col = "Average spatial autocorrelation score")+
  theme_bw()+
  theme(legend.position = "top",#legend.box = "horizontal",
        text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

comp_corSLE1_pct <-all_pearson_pct %>% arrange(score) %>%
  filter(indep) %>%
  ggplot(aes(x = (SLE1_cor.sd), y = (SLE1_cor.vst), label = paste(Gene1, Gene2),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline()+
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "Correlation for SLE\nwith SpaceDecorr",
       y = "Correlation for SLE\nwith SPARK",
       col = "Average spatial autocorrelation score")+
  theme_bw()+
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

comp_diffcor_pct <- all_pearson_pct %>% arrange(score) %>%
  filter(indep) %>%
  ggplot(aes(x = (SLE1_cor.sd - SLE0_cor.sd),
             y = (SLE1_cor.vst-SLE0_cor.vst), label = paste(Gene1, Gene2),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline()+
  geom_hline(yintercept = c(-0.1,0.1), lty = 2, col = "darkgrey") +
  geom_vline(xintercept = c(-0.1,0.1), lty = 2, col = "darkgrey")+
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "Correlation Difference\nSpaceDecorr",
       y = "Correlation Difference\nSPARK",
       col = "Average spatial autocorrelation score")+
  theme_bw()+
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = 15),
        strip.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())


cowplot::plot_grid(legend = get_plot_component(comp_pvalue_pct, 'guide-box-top', return_all = TRUE),
                   plot_grid(comp_pvalue_pct + theme(legend.position = "none"),
                             comp_corSLE0_pct+ theme(legend.position = "none"),
                             comp_corSLE1_pct+ theme(legend.position = "none"),
                             comp_diffcor_pct+ theme(legend.position = "none"),
                             nrow = 1, rel_widths = c(1,1,1) ,labels = c("A","B","C","D")
                   ),
                   nrow = 2,rel_heights = c(0.15,0.85))

### With labels ----
svgs_forplot_pct <- svg_pct0 %>%
  mutate(treat = ifelse(grepl("SLE",tissue),"SLE","Control")) %>%
  group_by(genes, treat) %>%
    score = mean(score),
    av_expr = mean(av_expr)
  ) %>%
  pivot_wider(names_from = treat, values_from = c(score, av_expr))

plot_svgs_score_pct <- svgs_forplot_pct%>%
  ggplot(aes(x = score_SLE, y = score_Control)) +
  geom_point(cex = 1) +
  geom_label_repel(data = svgs_forplot_pct %>% filter(genes %in% c("GPX3","WIF1")),
                   aes(label = genes),col = "red",
                   box.padding = 1, min.segment.length = 0.01, size = 2.5)+
  theme_bw()+
  labs(x = "Average spatial autocorrelation score in SLEs",
       y = "Average spatial autocorrelation score in Controls") +
  theme(text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_svgs_av_expr_pct <- svgs_forplot_pct%>%
  ggplot(aes(x = av_expr_SLE, y = av_expr_Control)) +
  geom_point(cex = 1) +
  geom_label_repel(data = svgs_forplot_pct %>% filter(genes %in% c("GPX3","WIF1")),
                   aes(label = genes),col="red",
                   box.padding = 1, min.segment.length = 0.01, size = 2.5)+
  theme_bw()+
  labs(x = "Average Expression in SLEs",
       y = "Average Expression in Controls") +
  theme(text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())

### Giotto -----
all_pearson_pct_giotto <- full_join(ddcor_pct_pearson,ddcor_pct_pearson_giotto, by = c("Gene1","Gene2"), suffix = c(".sd",".gt") )
all_pearson_pct_giotto<-all_pearson_pct_giotto %>% merge(.,svg_pct, by.x = "Gene1", by.y = "genes", all.x = T) %>%
  merge(.,svg_pct, by.x = "Gene2", by.y = "genes", all.x = T) %>%
  mutate(
    score = rowMeans(cbind(as.numeric(score.x), as.numeric(score.y)), na.rm = TRUE)
  )



comp_pvalue_pct_giotto <- all_pearson_pct_giotto %>%arrange(score) %>%
  ggplot(aes(x = -log10(pValDiff_adj.sd), y = -log10(pValDiff_adj.gt),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = -log10(0.05)) +
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "-log10(pvalue SpaceDecorr)",
       y = "-log10(pvalue Giotto)",
       col = "Average spatial autocorrelation score")+
  theme_bw()+
  theme(legend.position = "top",legend.box = "horizontal",
        legend.key.width = unit(1, "cm"),
        text = element_text(size = 10),
        legend.title = element_text(margin = margin(r = 20, b=10)),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.spacing = unit(1, "lines"),
        legend.margin = margin(t = 10, b = 10),
        legend.key.size = unit(1.5, "lines"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())


comp_corSLE0_pct_giotto <- all_pearson_pct_giotto %>%arrange(score) %>%
  ggplot(aes(x = (SLE0_cor.sd), y = (SLE0_cor.gt),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline() +
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "Correlation for controls with SpaceDecorr",
       y = "Correlation for controls with Giotto",
       col = "Average spatial autocorrelation score")+
  #geom_hline(yintercept = c(-0.1,0.1), lty = 2, col ="darkgrey") +
  #geom_vline(xintercept = c(-0.1,0.1), lty = 2, col ="darkgrey") +
  theme_bw()+
  theme(legend.position = "top",#legend.box = "horizontal",
        text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())

comp_corSLE1_pct_giotto <-all_pearson_pct_giotto %>% arrange(score) %>%
  ggplot(aes(x = (SLE1_cor.sd), y = (SLE1_cor.gt), label = paste(Gene1, Gene2),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_abline()+
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "Correlation for SLE with SpaceDecorr",
       y = "Correlation for SLE with Giotto",
       col = "Average spatial autocorrelation score")+
  #geom_hline(yintercept = c(-0.1,0.1), lty = 2, col ="darkgrey") +
  #geom_vline(xintercept = c(-0.1,0.1), lty = 2, col ="darkgrey") +
  theme_bw()+
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())

comp_diffcor_pct_giotto <- all_pearson_pct_giotto %>% arrange(score) %>%
  ggplot(aes(x = (SLE1_cor.sd - SLE0_cor.sd),
             y = (SLE1_cor.gt-SLE0_cor.gt), label = paste(Gene1, Gene2),
             col = score)) + geom_point(cex=0.7) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = c(-0.1,0.1), lty = 2, col ="darkgrey") +
  geom_vline(xintercept = c(-0.1,0.1), lty = 2, col ="darkgrey") +
  geom_abline()+
  scale_color_viridis(option = "C", direction = -1) +
  labs(x = "Correlation Difference SpaceDecorr",
       y = "Correlation Difference Giotto",
       col = "Average spatial autocorrelation score")+
  theme_bw()+
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = 10),
        strip.text = element_text(size = 10),
        legend.text = element_text(size = 10),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        strip.background = element_blank())


scatter_pct_giotto <- cowplot::plot_grid(legend = get_plot_component(comp_pvalue_pct_giotto, 'guide-box-top', return_all = TRUE),
                                         plot_grid(comp_pvalue_pct_giotto + theme(legend.position = "none"),
                                                   comp_corSLE0_pct_giotto+ theme(legend.position = "none"),
                                                   comp_corSLE1_pct_giotto+ theme(legend.position = "none"),
                                                   comp_diffcor_pct_giotto+ theme(legend.position = "none"),
                                                   nrow = 1, rel_widths = c(1,1,1) ,labels = c("A","B","C","D")
                                         ),
                                         nrow = 2,rel_heights = c(0.15,0.85))
scatter_pct_giotto

## Graph - Hub genes ----
### giotto ----
sig_edges_pct_giotto <- ddcor_pct_pearson_giotto %>%
  filter(pValDiff < 0.05,# class == "0/+",
         abs(SLE1_cor - SLE0_cor) > 0.3,
         zScoreDiff > 0)%>%
  rowwise() %>%
  mutate(
    delta_rho = SLE1_cor - SLE0_cor,
    sign = ifelse(delta_rho > 0, "SLE > Control", "Control > SLE"),
    pair = paste(sort(c(Gene1, Gene2)), collapse = "__"),
    known = ifelse(pair %in% stringscores$pair, "known","not known")
  )

graph_pct_giotto <- graph_from_data_frame(
  d = sig_edges_pct_giotto %>% dplyr::select(Gene1, Gene2, delta_rho, sign,known,Classes),
  directed = FALSE
)

# compute modules
set.seed(124)
clust_pct_giotto <- cluster_louvain(graph_pct_giotto)
V(graph_pct_giotto)$module_louvain <- membership(clust_pct_giotto)#factor(membership(clust_pct_giotto), levels = c(1,3,2), labels = c(1:3))

# get hubs per module
V(graph_pct_giotto)$degree <- igraph::degree(graph_pct_giotto)
hub_genes_pct_giotto <- names(V(graph_pct_giotto))[V(graph_pct_giotto)$degree >= quantile(V(graph_pct_giotto)$degree, 0.90)]
node_df_pct_giotto <- tibble(
  gene = V(graph_pct_giotto)$name,
  module_louvain = V(graph_pct_giotto)$module_louvain,
  degree = V(graph_pct_giotto)$degree
)

hubs_by_module_pct_giotto <- node_df_pct_giotto %>%
  group_by(module_louvain) %>%
  filter(degree >= quantile(degree, 0.95)) %>%
  arrange(desc(degree)) %>%
  ungroup()

V(graph_pct_giotto)$is_hub <- V(graph_pct_giotto)$name %in% hubs_by_module_pct_giotto$gene
V(graph_pct_giotto)$label <- ifelse(V(graph_pct_giotto)$is_hub, V(graph_pct_giotto)$name, NA)


plot_graph_pct_giotto <- ggraph(as_tbl_graph(graph_pct_giotto), layout = "fr") +  # "fr" = Fruchterman-Reingold layout
  geom_edge_link( alpha = 0.7) +
  geom_node_point(aes(color = as.factor(module_louvain), pch = is_hub)) +
  geom_node_label(aes(label = label, col =  as.factor(module_louvain)), force = 3, show.legend = F,
                  repel = TRUE, size = 4,box.padding = 1, min.segment.length = 0.01, max.overlaps = 50,
                  arrow = arrow(length = unit(0.015, "npc")),
                  point.padding = 0.2) +
  scale_edge_color_manual(values = c("SLE > Control" = "grey30", "Control > SLE" = "grey")) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8","3" = "#4DAF4A","4" = "#984EA3","5" = "#FF7F00","6" = "#B79F00")) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_size(range = c(2, 6)) +
  theme_void() +
  labs(edge_color = "Correlation", color = "Module", pch = "Hub gene") +
  #ggtitle("PCT (Giotto)") +
  theme(legend.position = "top",
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        #panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

### vst ----
sig_edges_pct_vst <- ddcor_pct_pearson_vst %>%
  filter(pValDiff < 0.05, #class == "0/+",
         abs(SLE1_cor - SLE0_cor) > 0.1,
         zScoreDiff > 0)%>%
  rowwise() %>%
  mutate(
    delta_rho = SLE1_cor - SLE0_cor,
    sign = ifelse(delta_rho > 0, "SLE > Control", "Control > SLE"),
    pair = paste(sort(c(Gene1, Gene2)), collapse = "__"),
    known = ifelse(pair %in% stringscores$pair, "known","not known")
  )

graph_pct_vst <- graph_from_data_frame(
  d = sig_edges_pct_vst %>% dplyr::select(Gene1, Gene2, delta_rho, sign,known,Classes),
  directed = FALSE
)

# compute modules
set.seed(124)
clust_pct_vst <- cluster_louvain(graph_pct_vst)
V(graph_pct_vst)$module_louvain <- membership(clust_pct_vst)#factor(membership(clust_pct_vst), levels = c(1,3,2), labels = c(1:3))

# get hubs per module
V(graph_pct_vst)$degree <- igraph::degree(graph_pct_vst)
hub_genes_pct_vst <- names(V(graph_pct_vst))[V(graph_pct_vst)$degree >= quantile(V(graph_pct_vst)$degree, 0.90)]
node_df_pct_vst <- tibble(
  gene = V(graph_pct_vst)$name,
  module_louvain = V(graph_pct_vst)$module_louvain,
  degree = V(graph_pct_vst)$degree
)

hubs_by_module_pct_vst <- node_df_pct_vst %>%
  group_by(module_louvain) %>%
  filter(degree >= quantile(degree, 0.95)) %>%
  arrange(desc(degree)) %>%
  ungroup()

V(graph_pct_vst)$is_hub <- V(graph_pct_vst)$name %in% hubs_by_module_pct_vst$gene
V(graph_pct_vst)$label <- ifelse(V(graph_pct_vst)$is_hub, V(graph_pct_vst)$name, NA)


plot_graph_pct_vst <- ggraph(as_tbl_graph(graph_pct_vst), layout = "fr") +  # "fr" = Fruchterman-Reingold layout
  geom_edge_link( alpha = 0.7) +
  geom_node_point(aes(color = as.factor(module_louvain), pch = is_hub)) +
  geom_node_label(aes(label = label, col =  as.factor(module_louvain)), force = 3, show.legend = F,
                  repel = TRUE, size = 4,box.padding = 1, min.segment.length = 0.01, max.overlaps = 50,
                  arrow = arrow(length = unit(0.015, "npc")),
                  point.padding = 0.2) +
  scale_edge_color_manual(values = c("SLE > Control" = "grey30", "Control > SLE" = "grey")) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8","3" = "#4DAF4A","4" = "#984EA3","5" = "#FF7F00","6" = "#B79F00")) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_size(range = c(2, 6)) +
  theme_void() +
  labs(edge_color = "Correlation", color = "Module", pch = "Hub gene") +
  #ggtitle("SPARK") +
  theme(legend.position = "top",
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        #panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )

# Table
genesModules_pct_vst <- data.frame(gene = names(membership(clust_pct_vst)),
                                   module_louvain = as.vector(membership(clust_pct_vst))) %>%
  group_by(module_louvain) %>%
  summarise(genes = paste(gene, collapse = ", ")) %>%
  left_join(hubs_by_module_pct_vst %>%
              group_by(module_louvain) %>%
              summarise(Hubs = paste(gene, collapse = ", ")) %>%
              mutate(module_louvain = as.numeric(module_louvain))) %>%
  mutate(method = "SPARK")


# Modules --
graph_tbl_add_pct_vst <- as_tbl_graph(graph_pct_vst)

# Split by module
modules_add_pct_vst <- graph_tbl_add_pct_vst %N>% dplyr::pull(module_louvain) %>% unique()

# Create a list of subgraphs
subgraphs_add_pct_vst <- lapply(modules_add_pct_vst, function(m_add_pct_vst) {
  induced_subgraph(graph_tbl_add_pct_vst, V(graph_tbl_add_pct_vst)$module_louvain == m_add_pct_vst)
})
names(subgraphs_add_pct_vst) <- modules_add_pct_vst

# Make individual ggraph plots
make_plot_add_pct_vst <- function(graph_sub, module_id) {
  ggraph(graph_sub, layout = "fr") +
    geom_edge_link(alpha = 1) +
    geom_node_point(aes(color = as.factor(module_louvain), pch = is_hub), size = 10) +
    geom_node_text(aes(label = name), size = 1.5, repel = FALSE) +
    scale_edge_color_manual(values = c("SLE > Control" = "grey30", "Control > SLE" = "grey")) +
    scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
                                  "4" = "#984EA3", "5" = "#FF7F00", "6" = "#B79F00")) +
    guides(
      color = "none",
      shape = guide_legend(override.aes = list(size = 2.5), title = "Hub gene"),
      edge_color = guide_legend(title = "")
    ) +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste("Module", module_id))
}

# Generate plots without legends
plots_add_pct_vst <- lapply(names(subgraphs_add_pct_vst), function(m_add_pct_vst) {
  make_plot_add_pct_vst(subgraphs_add_pct_vst[[m_add_pct_vst]], m_add_pct_vst)
})

legend_plot_add <- ggraph(subgraphs_add_pct_vst[[1]], layout = "fr") +
  geom_edge_link( alpha = 1) +
  geom_node_point(aes(color = as.factor(module_louvain), pch = is_hub), size = 10) +
  scale_edge_color_manual(values = c("SLE > Control" = "grey30", "Control > SLE" = "grey")) +
  guides(
    color = "none",
    shape = guide_legend(override.aes = list(size = 3), title = "Hub gene"),
    edge_color = guide_legend(title = "")
  ) +
  theme_void() +
  theme(legend.position = "top", text = element_text(size = 15))+
  theme(
    legend.position = "top",
    legend.box = "horizontal",                        # stack guides vertically
    legend.spacing = unit(1, "lines"),              # spacing between keys
    legend.margin = margin(t = 10, b = 10),          # margin around the whole legend
    legend.key.size = unit(1.5, "lines"),           # size of the legend key boxes
    #text = element_text(size = 15)
  )


# Combine legend and plots
Graph_PCT_VST_Modules <- cowplot::plot_grid(
  get_plot_component(legend_plot_add, 'guide-box-top', return_all = TRUE),
  cowplot::plot_grid(plotlist = plots_add_pct_vst, ncol = 4, align = "hv"),
  ncol = 1,
  rel_heights = c(0.1, 1)
)


### sd ----
sig_edges_pct <- ddcor_pct_pearson %>%
  filter(pValDiff < 0.05, #class == "0/+",
         abs(SLE1_cor - SLE0_cor) > 0.1,
         zScoreDiff >0
  )%>%
  rowwise %>%
  mutate(
    delta_rho = SLE1_cor - SLE0_cor,
    sign = ifelse(delta_rho > 0, "SLE > Control", "Control > SLE"),
    pair = paste(sort(c(Gene1, Gene2)), collapse = "__"),
    known = ifelse(pair %in% stringscores$pair, "known","not known")
  )

graph_pct <- graph_from_data_frame(
  d = sig_edges_pct %>% dplyr::select(Gene1, Gene2, delta_rho, sign, known,Classes),
  directed = FALSE
)

# compute modules
set.seed(124)
clust_pct <- cluster_louvain(graph_pct)
V(graph_pct)$module_louvain <- membership(clust_pct)

names(V(graph_pct))[!names(V(graph_pct)) %in% names(V(graph_pct_vst))]

# get hubs per module
V(graph_pct)$degree <- igraph::degree(graph_pct)
hub_genes_pct <- names(V(graph_pct))[V(graph_pct)$degree >= quantile(V(graph_pct)$degree, 0.90)]
node_df_pct <- tibble(
  gene = V(graph_pct)$name,
  module_louvain = V(graph_pct)$module_louvain,
  degree = V(graph_pct)$degree
)

hubs_by_module_pct <- node_df_pct %>%
  group_by(module_louvain) %>%
  filter(degree >= quantile(degree, 0.95)) %>%
  arrange(desc(degree)) %>%
  ungroup()

V(graph_pct)$is_hub <- V(graph_pct)$name %in% hubs_by_module_pct$gene
V(graph_pct)$label <- ifelse(V(graph_pct)$is_hub, V(graph_pct)$name, NA)


plot_graph_pct_sd <- ggraph(as_tbl_graph(graph_pct), layout = "fr") +  # "fr" = Fruchterman-Reingold layout
  geom_edge_link( alpha = 0.7) +
  geom_node_point(aes(color = as.factor(module_louvain), pch = is_hub)) +
  geom_node_label(aes(label = label, col =  as.factor(module_louvain)), force = 3, show.legend = F,
                  repel = TRUE, size = 4,box.padding = 1, min.segment.length = 0.01, max.overlaps = 50,
                  arrow = arrow(length = unit(0.015, "npc")),
                  point.padding = 0.2) +
  scale_edge_color_manual(values = c("SLE > Control" = "grey30", "Control > SLE" = "grey")) +
  scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8","3" = "#4DAF4A","4" = "#984EA3","5" = "#FF7F00","6" = "#B79F00",
                                "7" = "#A65628","8" = "#F781BF")) +
  scale_edge_width(range = c(0.5, 2)) +
  scale_size(range = c(2, 6)) +
  theme_void() +
  labs(edge_color = "Correlation", color = "Module", pch = "Hub gene") +
  #ggtitle("SpaceDecorr") +
  theme(legend.position = "top",
        text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5),
        #panel.border = element_rect(color = "black", fill = NA, size = 0.5)
  )


### Table ----
genesModules_pct_sd <- data.frame(gene = names(membership(clust_pct)),
                                  module_louvain = as.vector(membership(clust_pct))) %>%
  group_by(module_louvain) %>%
  summarise(genes = paste(gene, collapse = ", ")) %>%
  left_join(hubs_by_module_pct %>%
              group_by(module_louvain) %>%
              summarise(Hubs = paste(gene, collapse = ", ")) %>%
              mutate(module_louvain = as.numeric(module_louvain))) %>%
  mutate(method = "SpaceDecorr")

dir.results <- file.path(dir.out, "results", "lupus")
dir.create(dir.results, recursive = TRUE, showWarnings = FALSE)

rbind(genesModules_pct_vst, genesModules_pct_sd) %>%
  write.csv(file.path(dir.results, "Tables_membership_pct.csv"))

# Modules --
graph_tbl_add_pct <- as_tbl_graph(graph_pct)

# Split by module
modules_add_pct <- graph_tbl_add_pct %N>% dplyr::pull(module_louvain) %>% unique()

# Create a list of subgraphs
subgraphs_add_pct <- lapply(modules_add_pct, function(m_add_pct) {
  induced_subgraph(graph_tbl_add_pct, V(graph_tbl_add_pct)$module_louvain == m_add_pct)
})
names(subgraphs_add_pct) <- modules_add_pct

# Make individual ggraph plots
make_plot_add_pct <- function(graph_sub, module_id) {
  ggraph(graph_sub, layout = "fr") +
    geom_edge_link( alpha = 1) +
    geom_node_point(aes(color = as.factor(module_louvain), pch = is_hub), size = 10) +
    geom_node_text(aes(label = name), size = 2, repel = FALSE) +
    scale_edge_color_manual(values = c("SLE > Control" = "grey30", "Control > SLE" = "grey")) +
    scale_color_manual(values = c("1" = "#E41A1C", "2" = "#377EB8","3" = "#4DAF4A","4" = "#984EA3","5" = "#FF7F00","6" = "#B79F00",
                                  "7" = "#A65628","8" = "#F781BF")) +
    guides(
      color = "none",
      shape = guide_legend(override.aes = list(size = 3), title = "Hub gene"),
      edge_color = guide_legend(title = "")
    ) +
    theme_void() +
    theme(legend.position = "none",
          text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5)) +
    ggtitle(paste("Module", module_id))
}

# Generate plots without legends
plots_add_pct <- lapply(names(subgraphs_add_pct), function(m_add_pct) {
  make_plot_add_pct(subgraphs_add_pct[[m_add_pct]], m_add_pct)
})


# Combine legend and plots
Graph_PCT_SD_Modules <-cowplot::plot_grid(
  get_plot_component(legend_plot_add, 'guide-box-top', return_all = TRUE),
  cowplot::plot_grid(plotlist = plots_add_pct, ncol = 4, align = "hv"),
  ncol = 1,
  rel_heights = c(0.1, 1)
)
