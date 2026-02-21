#----------------------------------
#---- Real data Lung WGCNA
#----------------------------------

dir.data <- "/path/to/zenodo/download/"   # Zenodo data folder
dir.out  <- "/path/to/output/"            # where figures/tables will be saved


library(tidyverse)
library(Giotto)
library(grid)
library(gridExtra)
library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(RColorBrewer)
library(igraph)
library(scico)

# STRINGdb and Pathway distance scores -----
load(paste0(dir.out,"Lung_Giotto_Filtered_Tumor.RData"))
stringscores <- read_csv(paste0(dir.data,"KnownInteractionsDatabases.csv"))[,-1] %>%
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

indepgenes <- stringscores %>%
  filter(pathway_distance == 1,
         stringdp_score < 100,
         biogrid == 0,
         is.na(trrust), is.na(dorothea), is.na(remap))
pair_indepgenes <- paste(indepgenes$gene_from, indepgenes$gene_to, sep = "-")


# Lung 6 - Correlation from all models -----
id_tissue_Lung6 <- pDataDT(gem)$cell_ID[pDataDT(gem)$tissue == "Lung6"]
gem_tissue_Lung6 <- subsetGiotto(gem, cell_ids = id_tissue_Lung6)

locs_lung6 <- gem_tissue_Lung6@spatial_locs
rownames(locs_lung6) <- locs_lung6$cell_ID
normexpr_lung6 <- as.data.frame(as.matrix(t(gem_tissue_Lung6@norm_expr)))
normexpr_lung6$cell_ID <- rownames(normexpr_lung6)
normexpr_lung6 <- merge(normexpr_lung6,locs_lung6, by = "cell_ID" )
normexpr_lung6 <- merge(normexpr_lung6,gem_tissue_Lung6@cell_metadata, by = "cell_ID" )


# SVG
load(paste0(dir.data,"km_spatialgenesLung6.RData"))
svgs_Lung6 <- km_spatialgenes %>% filter(adj.p.value < 0.05)


# Contamination
contamination0 <- read.csv(paste0(dir.out,"ContaminationMetricsGenes_Lung.csv"))
contamination0 <- contamination0 %>% filter(tissue == "Lung6", cell_type == "tumor")
contamination <- contamination0 %>% filter(avg_cluster > 0.1, ratio < 1, tissue == "Lung6", cell_type == "tumor")

# Giotto
load(paste0(dir.data,"Data_Lung6_giotto.Rdata"))
cor_giotto_Lung6 <- cor_giotto
cor_giotto_long_Lung6 <- cor_giotto_Lung6$cor_DT %>%
  dplyr::select(gene_ID, variable, spat_cor) %>%
  dplyr::rename(Var1 = gene_ID, Var2 = variable, Giotto=spat_cor)


# Spacedecorr
# k 100
load(paste0(dir.data,"Data_Lung6_spacedecorr_k100.Rdata"))
edf100 <- data_spacedecorr$edf
summary(as.numeric(edf100))
spacedecorr_Lung6_k100 <- cbind(data_spacedecorr$Residuals,gem_tissue_Lung6@cell_metadata)
spacedecorr_Lung6_spatial_k100 <- cbind(data_spacedecorr$Splines,gem_tissue_Lung6@cell_metadata,gem_tissue_Lung6@spatial_locs)
cor_spacedecorr_long_Lung6_k100 <- reshape2::melt(cor(data_spacedecorr$Residuals)) %>%
  dplyr::rename(spacedecorr100 = value)
time_spacedecorr_k100 <- time_spacedecorr
# 648.630s - 10.81min


# k 200
load(paste0(dir.data,"Data_Lung6_spacedecorr_k200.Rdata"))
edf200 <- data_spacedecorr$edf
summary(as.numeric(edf200))

spacedecorr_Lung6_k200 <- cbind(data_spacedecorr$Residuals,gem_tissue_Lung6@cell_metadata)
spacedecorr_Lung6_spatial_k200 <- cbind(data_spacedecorr$Splines,gem_tissue_Lung6@cell_metadata,gem_tissue_Lung6@spatial_locs)
cor_spacedecorr_long_Lung6_k200 <- reshape2::melt(cor(data_spacedecorr$Residuals)) %>%
  dplyr::rename(spacedecorr200 = value)
time_spacedecorr_k200 <- time_spacedecorr
# 1232.404s - 20.54min

plot(spacedecorr_Lung6_k100$TYK2,spacedecorr_Lung6_k200$TYK2, )
plot(normexpr_lung6$TYK2,spacedecorr_Lung6_k100$TYK2 )

# k 1000
load(paste0(dir.data,"Data_Lung6_spacedecorr_k1000.Rdata"))
spacedecorr_Lung6_k1000 <- cbind(data_spacedecorr$Residuals,gem_tissue_Lung6@cell_metadata)
spacedecorr_Lung6_spatial_k1000 <- cbind(data_spacedecorr$Splines,gem_tissue_Lung6@cell_metadata,gem_tissue_Lung6@spatial_locs)

cor_spacedecorr_long_Lung6_k1000 <- reshape2::melt(cor(data_spacedecorr$Residuals)) %>%
  dplyr::rename(spacedecorr1000 = value)
time_spacedecorr_k1000 <- time_spacedecorr
# 73880.47s - 1231.34min - 20.52h


# Pearson
cor_pearson_long_Lung6 <- reshape2::melt(cor(t(as.matrix(gem_tissue_Lung6@raw_exprs)))) %>%dplyr::rename(pearson = value)
Ylog <- apply(t(as.matrix(gem_tissue_Lung6@raw_exprs)),2, function(x)log(10e4 * x/gem_tissue_Lung6@cell_metadata$totalcounts +1))
cor_logpearson_long_Lung6 <- reshape2::melt(cor(Ylog)) %>%dplyr::rename(logpearson = value)
cor_spearman_long_Lung6 <- reshape2::melt(cor(t(as.matrix(gem_tissue_Lung6@raw_exprs)), method = "spearman")) %>%dplyr::rename(spearman = value)


# vst
load(paste0(dir.data,"Data_Lung6_vst.Rdata"))
cor_vst_long_Lung6 <- reshape2::melt(cor(t(data_vst$y))) %>%
  dplyr::rename(vst = value)

# combining all
cor_all_long_Lung6 <- cbind(cor_spacedecorr_long_Lung6_k200) %>%
  left_join(cor_spacedecorr_long_Lung6_k1000) %>%
  left_join(cor_spacedecorr_long_Lung6_k100) %>%
  left_join(cor_giotto_long_Lung6) %>%
  left_join(cor_pearson_long_Lung6) %>% left_join(cor_logpearson_long_Lung6) %>% left_join(cor_spearman_long_Lung6) %>%
  left_join(cor_vst_long_Lung6) %>%
  merge(.,stringscores, by.x = c("Var1","Var2"), by.y = c("gene_from","gene_to"))  %>%
  filter(Var1 != Var2) %>%
  merge(., contamination0[,-1], by.x = "Var1", by.y = "target",all.x = T ) %>%
  merge(., contamination0[,-1], by.x = "Var2", by.y = "target",all.x = T, suffixes = c(".v1",".v2") ) %>%
  #svg
  merge(., svgs_Lung6[,c(1,5,6)], by.x = "Var1", by.y = "genes",all.x = T ) %>%
  merge(., svgs_Lung6[,c(1,5,6)], by.x = "Var2", by.y = "genes",all.x = T, suffixes = c(".v1",".v2") ) %>%
  filter(avg_cluster.v1 > .1, avg_cluster.v2 > .1) %>%
  filter(Var1 %in% contamination$target, Var2 %in% contamination$target) %>%
  mutate(independent = ifelse(stringdp_score == 0 & pathway_distance == 1 & biogrid == 0 , T, F))
write.csv(cor_all_long_Lung6, file = file.path(dir.out, "cor_all_long_Lung6.csv"))

pvalue_all_Lung6 <- cor_all_long_Lung6 %>% pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  mutate(pvalue = test1cor(cor,66261)) %>%
  group_by(model) %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))


# Plots -----
## Correlation ----
textsize=12
cor_map <- function(gene1label, gene2label, k=500, lower=-1, upper=1){
  rdbu <- brewer.pal(9, "RdBu")
  red_side  <- rdbu[1:5]
  blue_side <- rdbu[5:9]

  flat_red  <- red_side[1]
  flat_blue <- blue_side[5]

  info_clust <- kmeans(normexpr_lung6[,c('sdimx','sdimy')], k, iter.max = 50)
  normexpr_lung6 <- normexpr_lung6 %>% dplyr::mutate(cluster = info_clust$cluster)

  normexpr_lung6$gene1 <- normexpr_lung6[,gene1label]
  normexpr_lung6$gene2 <- normexpr_lung6[,gene2label]
  corrbyclust <- normexpr_lung6 %>% group_by(cluster) %>%
    summarise(cor = cor(gene1, gene2))
  normexpr_lung6 %>% left_join(corrbyclust, by = "cluster") %>%
    ggplot(aes(x = sdimx, y = sdimy, col = cor)) +
    geom_point(cex = 0.5) +
    scale_color_gradient2(mid = "lightgrey", low = flat_blue, high =flat_red, limits = c(lower, upper),
                          breaks = c(-0.3,0,0.3))+
    labs(col = "Correlation") +
    facet_wrap(paste(gene1label,gene2label, sep = " and ")~.) +
    theme_bw() + #coord_fixed() +
    guides(col = guide_colorbar(title.vjust = 0.7,
                                barheight = 0.5,
                                barwidth = 8))+
    theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
          legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
          legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
          panel.grid.minor = element_blank()
    )
}

plot_cormap <- cor_map("MZT2A","DDR1", k =50, lower = -0.3, upper = 0.3);plot_cormap

# correlation for other models
pvalue_all_Lung6 %>%
  filter(independent) %>%
  group_by(model) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  filter(Var2 == "MZT2A", Var1 == "DDR1") %>%
  mutate(sign = as.numeric(padj < 0.05)) %>%
  dplyr::select(Var1,Var2,stringdp_score,pathway_distance, model, cor, pvalue, padj,sign )

## Random effects ----
load(paste0(dir.out,"RealData/Lung6/spacedecorr_Lung6_spatial.Rdata"))

plot_spatialeffect <- spacedecorr_Lung6_spatial_k200 %>%
  pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
  filter(!is.na(Expression))%>%
  arrange(Expression) %>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=.8) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_wrap(gene~.) +
  theme_bw() + labs(col = "Spatial Effect")+
  guides(col = guide_colorbar(title.vjust = 0.7,
                              barheight = 0.5,
                              barwidth = 8))+
  theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
        legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
        legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
        panel.grid.minor = element_blank()
  )

plot_spatialeffect_diffk <- spacedecorr_Lung6_spatial_k1000 %>%
  dplyr::select(sdimx, sdimy,MZT2A,DDR1 ) %>%
  pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
  mutate(k = "k=1000") %>%
  rbind(spacedecorr_Lung6_spatial_k200 %>%
          dplyr::select(sdimx, sdimy,MZT2A,DDR1 ) %>%
          pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
          mutate(k = "k=200"))%>%
  rbind(spacedecorr_Lung6_spatial_k100 %>%
          dplyr::select(sdimx, sdimy,MZT2A,DDR1 ) %>%
          pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
          mutate(k = "k=100"))%>%
  mutate(k = factor(k, levels = c("k=100","k=200","k=1000"))) %>%
  filter(!is.na(Expression))%>%
  arrange(Expression) %>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=.8) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_grid(k~gene) +
  theme_bw() + labs(col = "Spatial Effect")+
  guides(col = guide_colorbar(title.vjust = 0.7,
                              barheight = 0.5,
                              barwidth = 8))+
  theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
        legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
        legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
        panel.grid.minor = element_blank(),
  );plot_spatialeffect_diffk


## Density plots ----
mycol <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"  )

# independent
pvalue_all_Lung6 %>%
  filter(independent) %>%
  group_by(model) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  mutate(t1err = as.numeric(padj < 0.05)) %>%
  group_by(model) %>%
  summarise(t1err = sprintf("%.1f%%",100*mean(t1err)))

cor_all_long_Lung6 %>%
  filter(independent) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("logpearson","vst","Giotto","graphR","spacedecorr100","spacedecorr200","spacedecorr1000")) %>%
  mutate(model = factor(model,
                        levels = c("pearson","spearman","logpearson",
                                   "vst","Giotto","graphR","spacedecorr100","spacedecorr200","spacedecorr1000"),
                        labels = c("Pearson","Spearman","log-CPM",
                                   "SPARK","Giotto","GraphR","SpaceDecorr k=100","SpaceDecorr k=200","SpaceDecorr k=1000"))) %>%
  mutate(t1err = as.numeric(abs(cor) > 0.1)) %>%
  group_by(model) %>%
  summarise(t1err = sprintf("%.1f%%",100*mean(t1err)))


plot_density_indep <- cor_all_long_Lung6 %>%
  filter(independent) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("logpearson","vst","Giotto","graphR","spacedecorr200")) %>%
  mutate(model = factor(model,
                        levels = c("pearson","spearman","logpearson",
                                   "vst","Giotto","graphR","spacedecorr200"),
                        labels = c("Pearson","Spearman","log-CPM",
                                   "SPARK","Giotto","GraphR","SpaceDecorr"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",size = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("Negative Control set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 2) +
  xlim(-0.25,0.25) +
  scale_color_manual(values = mycol) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        legend.text = element_text(size = textsize+1),
        strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_density_indep_diffk <- cor_all_long_Lung6 %>%
  filter(independent) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("spacedecorr100","spacedecorr200","spacedecorr1000")) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorr100","spacedecorr200","spacedecorr1000"),
                        labels = c("k=100","k=200","k=1000"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",size = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("Negative Control set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 1,alpha = 0.7) +
  xlim(-0.05,0.05) +
  scale_color_manual(values = c("#FF7F00" ,"#984EA3","#00A9CE")) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        legend.text = element_text(size = textsize+1),
        strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

# dependent
pvalue_all_Lung6 %>%
  filter(proppath_factor == "More than 25%", stringdp_factor %in% c("High","Very high")) %>%
  mutate(model = factor(model,
                        levels = c("pearson","spearman","logpearson",
                                   "vst","Giotto","graphR","spacedecorr100"),
                        labels = c("Pearson","Spearman","Logpearson",
                                   "SPARK","Giotto","GraphR","SpaceDecorr"))) %>%
  mutate(t1err = as.numeric(padj < 0.05) ) %>%
  group_by(model) %>%
  summarise(t1err = sprintf("%.1f%%",100*mean(t1err)))

plot_density_dep <- cor_all_long_Lung6 %>%
  filter(proppath_factor == "More than 25%", stringdp_factor %in% c("High","Very high")) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("logpearson","vst","Giotto","graphR","spacedecorr200")) %>%
  mutate(model = factor(model,
                        levels = c("pearson","spearman","logpearson",
                                   "vst","Giotto","graphR","spacedecorr200"),
                        labels = c("Pearson","Spearman","log-CPM",
                                   "SPARK","Giotto","GraphR","SpaceDecorr"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",size = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("High confidence set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 2) +
  xlim(-1,1) +
  scale_color_manual(values = mycol) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_density_dep_diffk <- cor_all_long_Lung6 %>%
  filter(proppath_factor == "More than 25%", stringdp_factor %in% c("High","Very high")) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("spacedecorr100","spacedecorr200","spacedecorr1000")) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorr100","spacedecorr200","spacedecorr1000"),
                        labels = c("k=100","k=200","k=1000"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",size = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("High confidence set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 1,alpha=0.7) +
  xlim(-1,1) +
  scale_color_manual(values = c("#FF7F00" ,"#984EA3","#00A9CE")) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        strip.text = element_text(size = textsize+1),

        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_densities <- plot_grid(legend = get_plot_component(plot_density_indep, 'guide-box-top', return_all = TRUE),
                            plot_grid(plot_density_indep + theme(legend.position = "none"),
                                      plot_density_dep+ theme(legend.position = "none"),
                                      labels = c("","")),
                            ncol = 1, rel_heights =  c(0.15,0.85))
plot_densities_diffk <- plot_grid(legend = get_plot_component(plot_density_indep_diffk, 'guide-box-top', return_all = TRUE),
                                  plot_grid(plot_density_indep_diffk + theme(legend.position = "none"),
                                            plot_density_dep_diffk+ theme(legend.position = "none"),
                                            labels = c("","")),
                                  ncol = 1, rel_heights =  c(0.15,0.85))
plot_densities_diffk



# FOVs --------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(purrr)

fov_boxes <- tribble(
  ~fov, ~xmin, ~xmax, ~ymin, ~ymax,
  1,  28.6, 29.6, -19.6, -19.0,
  2,  29.6, 30.6, -19.6, -19.0,
  3,  30.6, 31.6, -19.6, -19.0,
  4,  31.6, 32.5, -19.6, -19.0,
  5,  32.5, 33.5, -19.6, -19.0,
  6,  28.6, 29.6, -19.0, -18.3,
  7,  29.6, 30.6, -19.0, -18.3,
  8,  30.6, 31.6, -19.0, -18.3,
  9,  31.6, 32.5, -19.0, -18.3,
  10,  32.5, 33.5, -19.0, -18.3,
  11,  28.6, 29.6, -18.3, -17.7,
  12,  29.6, 30.6, -18.3, -17.7,
  13,  30.6, 31.6, -18.3, -17.7,
  14,  31.6, 32.5, -18.3, -17.7,
  15,  32.5, 33.5, -18.3, -17.7,
  16,  28.6, 29.6, -17.7, -17.0,
  17,  29.6, 30.6, -17.7, -17.0,
  18,  30.6, 31.6, -17.7, -17.0,
  19,  31.6, 32.5, -17.7, -17.0,
  20,  32.5, 33.5, -17.7, -17.0,
  21,  28.6, 29.6, -17.0, -16.3,
  22,  29.6, 30.6, -17.0, -16.3,
  23,  30.6, 31.6, -17.0, -16.3,
  24,  31.6, 32.5, -17.0, -16.3,
  25,  32.5, 33.5, -17.0, -16.3,
  26,  28.6, 29.6, -16.3, -15.7,
  27,  29.6, 30.6, -16.3, -15.7,
  28,  30.6, 31.6, -16.3, -15.7,
  29,  31.6, 32.5, -16.3, -15.7,
  30,  32.5, 33.5, -16.3, -15.7
)
fov_boxes <- fov_boxes %>%
  mutate(
    center_x = (xmin + xmax) / 2,
    center_y = (ymin + ymax) / 2
  )
# 2. Create the plot
plot_fov <- normexpr_lung6 %>%
  ggplot(aes(x = sdimx, y = sdimy)) +  # Replace with the actual color variable
  geom_point(cex = 0.5, col = "gray") +

  # FoV borders
  geom_rect(data = fov_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = NA, color = "black", size = 0.3) +

  # FoV labels
  geom_text(data = fov_boxes, aes(x = center_x, y = center_y, label = fov),
            inherit.aes = FALSE,
            size = 3.5, color = "black", fontface = "bold") +
  facet_wrap("Field of Views (FOV)"~.)+
  theme_bw() +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "white", colour = "white"),
    legend.box.spacing = unit(0, "pt"),
    text = element_text(size = 15),
    legend.title = element_text(margin = margin(r = 20, b = 10))
  )

# WGCNA -------------------------------------------------------------------
# Aux functions ----------
clusterGenes <- function (cor_matrix,k = 10, hclust_method = "ward.D")
{
  cor_dist <- stats::as.dist(1 - cor_matrix)
  cor_h <- stats::hclust(d = cor_dist, method = hclust_method)
  cor_clus <- stats::cutree(cor_h, k = k)
  return(list(hclust = cor_h, clusters = cor_clus))
}


heatmapclust = function (clustobj,cor_matrix,use_clus_name = NULL, show_cluster_annot = TRUE,
                         show_row_dend = FALSE, show_column_dend = FALSE, show_row_names = FALSE,
                         show_column_names = FALSE, show_plot = NULL, return_plot = NULL,
                         show_clust_legend = TRUE, label_font_size = 10,
                         save_plot = NULL, save_param = list(), default_save_name = "heatmSpatialCorFeats",
                         ...)
{
  ha <- NULL
  hclust_part <- clustobj[["hclust"]]
  clusters_part <- clustobj[["clusters"]]
  clusters_part <- factor(clusters_part, levels = sort(unique(clusters_part)))
  uniq_clusters <- unique(clusters_part)
  uniq_clusters <- factor(uniq_clusters, levels = c(1:length(uniq_clusters)))
  mycolors <- getDistinctColors(length(uniq_clusters))
  names(mycolors) <- uniq_clusters
  ha <- ComplexHeatmap::HeatmapAnnotation(bar = (clusters_part),
                                          col = list(bar = mycolors),
                                          annotation_name_gp = gpar(fontsize = 0),
                                          annotation_legend_param = list(title = NULL),
                                          show_legend = show_clust_legend)
  # heatmap_colors <- colorRamp2(
  #   seq(min(as.matrix(cor_matrix)), max(as.matrix(cor_matrix)), length.out = 100),
  #   colorRampPalette(brewer.pal(9, "RdBu"))(100)
  # )

  heatm <- ComplexHeatmap::Heatmap(matrix = as.matrix(cor_matrix),
                                   cluster_rows = hclust_part, cluster_columns = hclust_part,
                                   show_row_dend = show_row_dend, show_column_dend = show_column_dend,
                                   show_row_names = show_row_names, show_column_names = show_column_names,
                                   top_annotation = ha,#col=heatmap_colors,
                                   heatmap_legend_param = list(
                                     title = NULL,
                                     legend_height = unit(3, "cm"),
                                     labels_gp = gpar(fontsize = label_font_size)  # increase label font size
                                   ),
                                   ...)
  return(heatm)
}

getEnrichClust <- function(clusters) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  universe <- mapIds(
    org.Hs.eg.db,
    keys = names(clusters),
    column = "ENTREZID",
    keytype = "SYMBOL"
  )
  do.call(rbind, lapply(1:max(clusters), function(clustID){
    enrichdat <- enrichGO(gene   =  universe[clusters == clustID],
                          universe      =  universe,
                          OrgDb         = org.Hs.eg.db,
                          ont           = "ALL",        # You can also use "CC" or "MF"
                          pAdjustMethod = "fdr",
                          pvalueCutoff  = 0.05,
                          qvalueCutoff  = 0.05)
    if(nrow(enrichdat@result) > 0){
      enrichdat@result$clustID <- clustID
      return(enrichdat@result)
    }
  }))
}

#-------------------------
k <- 20
#spacedecorr
# k = 100
load(paste0(dir.data,"Data_Lung6_spacedecorr_k100.Rdata"))
cormat_spacedecorr_Lung6_k100 <- cor(data_spacedecorr$Residuals[, colnames(data_spacedecorr$Residuals) %in% contamination$target])
clust_spacedecorr_Lung6_k100 <- clusterGenes(cormat_spacedecorr_Lung6_k100, k = k)
enrich_clust_spacedecorr_Lung6_k100 <- getEnrichClust(clust_spacedecorr_Lung6_k100$clusters)
enrich_clust_spacedecorr_Lung6_k100$model <- "spacedecorr k=100"


# k = 200
load(paste0(dir.data,"Data_Lung6_spacedecorr_k200.Rdata"))
cormat_spacedecorr_Lung6_k200 <- cor(data_spacedecorr$Residuals[, colnames(data_spacedecorr$Residuals) %in% contamination$target])
clust_spacedecorr_Lung6_k200 <- clusterGenes(cormat_spacedecorr_Lung6_k200, k = k)
enrich_clust_spacedecorr_Lung6_k200 <- getEnrichClust(clust_spacedecorr_Lung6_k200$clusters)
enrich_clust_spacedecorr_Lung6_k200$model <- "spacedecorr k=200"

# k = 1000
load(paste0(dir.data,"Data_Lung6_spacedecorr_k1000.Rdata"))
cormat_spacedecorr_Lung6_k1000 <- cor(data_spacedecorr$Residuals[, colnames(data_spacedecorr$Residuals) %in% contamination$target])
clust_spacedecorr_Lung6_k1000 <- clusterGenes(cormat_spacedecorr_Lung6_k1000, k = k)
enrich_clust_spacedecorr_Lung6_k1000 <- getEnrichClust(clust_spacedecorr_Lung6_k1000$clusters)
enrich_clust_spacedecorr_Lung6_k1000$model <- "spacedecorr k=1000"

# Giotto
cormat_giotto_Lung6 <- cor_giotto_Lung6$cor_DT %>%
  filter(gene_ID %in% contamination$target,
         variable %in% contamination$target) %>%
  graph_from_data_frame(., directed = F) %>%
  get.adjacency(., attr="spat_cor", sparse=FALSE)

clust_giotto_Lung6 <- clusterGenes(cormat_giotto_Lung6, k = k)
enrich_clust_giotto_Lung6 <- getEnrichClust(clust_giotto_Lung6$clusters)
enrich_clust_giotto_Lung6$model <- "giotto"

# vst
cormat_vst_Lung6 <- cor(t(data_vst$y[rownames(data_vst$y)%in% contamination$target,]))
clust_vst_Lung6 <- clusterGenes(cormat_vst_Lung6, k = k)
enrich_clust_vst_Lung6 <- getEnrichClust(clust_vst_Lung6$clusters)
enrich_clust_vst_Lung6$model <- "vst"


# log CPM
cormat_logpearson_Lung6 <- cor(Ylog[,colnames(Ylog) %in% contamination$target])
clust_logpearson_Lung6 <- clusterGenes(cormat_logpearson_Lung6, k = k)
enrich_clust_logpearson_Lung6 <- getEnrichClust(clust_logpearson_Lung6$clusters)
enrich_clust_logpearson_Lung6$model <- "logPearson"

save(cormat_spacedecorr_Lung6_k100,clust_spacedecorr_Lung6_k100,enrich_clust_spacedecorr_Lung6_k100,
     cormat_spacedecorr_Lung6_k200,clust_spacedecorr_Lung6_k200,enrich_clust_spacedecorr_Lung6_k200,
     cormat_spacedecorr_Lung6_k1000,clust_spacedecorr_Lung6_k1000,enrich_clust_spacedecorr_Lung6_k1000,
     cormat_giotto_Lung6,clust_giotto_Lung6,enrich_clust_giotto_Lung6,
     cormat_vst_Lung6,clust_vst_Lung6,enrich_clust_vst_Lung6,
     cormat_logpearson_Lung6,clust_logpearson_Lung6,enrich_clust_logpearson_Lung6,
     file = file.path(dir.out, "ExampleWGCNA_k20.Rdata.csv"))
load(paste0(dir.out,"Manuscript/RealData/Lung/ExampleWGCNA_k20.Rdata"))


# Extract RdBu colors
rdbu <- brewer.pal(9, "RdBu")
red_side  <- rdbu[1:5]  # red to white
blue_side <- rdbu[5:9]  # white to blue

# Choose flat red and blue shades
flat_red  <- red_side[1]
flat_blue <- blue_side[5]

# Define custom color scale
heatmap_colors_sd100 <- colorRamp2(
  c(min(cormat_spacedecorr_Lung6_k100),  0, max(upper(cormat_spacedecorr_Lung6_k100))),
  c(flat_blue,  "white", flat_red)
)
plot_wgcna_spacedecorr_k100 <- heatmapclust(clust_spacedecorr_Lung6_k100,cormat_spacedecorr_Lung6_k100,
                                            column_title = "k=100",col=heatmap_colors_sd100,
                                            show_heatmap_legend = F,show_clust_legend = F)

heatmap_colors_sd200 <- colorRamp2(
  c(min(cormat_spacedecorr_Lung6_k200),  0, max(upper(cormat_spacedecorr_Lung6_k200))),
  c(flat_blue,  "white", flat_red)
)
plot_wgcna_spacedecorr_k200 <- heatmapclust(clust_spacedecorr_Lung6_k200,cormat_spacedecorr_Lung6_k200,
                                            column_title = "k=200",col=heatmap_colors_sd200,
                                            show_heatmap_legend = F,show_clust_legend = F)

heatmap_colors_sd1k <- colorRamp2(
  c(min(cormat_spacedecorr_Lung6_k1000),  0, max(upper(cormat_spacedecorr_Lung6_k200))),
  c(flat_blue,  "white", flat_red)
)
plot_wgcna_spacedecorr_k1000 <- heatmapclust(clust_spacedecorr_Lung6_k1000,cormat_spacedecorr_Lung6_k1000,
                                             column_title = "k=1000",col=heatmap_colors_sd1k,
                                             show_heatmap_legend = T,show_clust_legend = F)

heatmap_colors_giotto <- colorRamp2(
  c(min(cormat_giotto_Lung6),  0, max(upper(cormat_giotto_Lung6))),
  c(flat_blue,  "white", flat_red)
)

plot_wgcna_giotto <- heatmapclust(clust_giotto_Lung6,cormat_giotto_Lung6,
                                  column_title = "Giotto",col=heatmap_colors_giotto,
                                  show_heatmap_legend = T,
                                  show_clust_legend = F);plot_wgcna_giotto

cluster_legend <- Legend(
  labels = sort(unique(clust_spacedecorr_Lung6_k1000[["clusters"]])),
  title = "Cluster",
  legend_gp = gpar(fill = getDistinctColors(length(unique(clust_spacedecorr_Lung6_k1000[["clusters"]])))),
  ncol = 2
)

## plot WGCNA ----
plot_clusters_wgcna_diffk <- plot_grid(
  grid.grabExpr(draw(plot_wgcna_spacedecorr_k100)),
  NULL,
  grid.grabExpr(draw(plot_wgcna_spacedecorr_k200)),
  NULL,
  grid.grabExpr(draw(plot_wgcna_spacedecorr_k1000,
                     annotation_legend_list = list(cluster_legend))),
  ncol = 5,rel_widths = c(1,0.1,1,0.1,1.21)
);plot_clusters_wgcna_diffk

plot_clusters_wgcna_diffk <- ggdraw() +
  draw_plot(plot_clusters_wgcna_diffk,
            x = 0.05, y = 0.05, width = 0.9, height = 0.9)

## Biological coherence ----
cor_all_long_Lung6 <- cor_all_long_Lung6 %>%
  mutate(physical = ifelse(biogrid == 1 | stringdp_score > 700, 1, 0),
         regulatory = ifelse(!is.na(trrust) | (!is.na(dorothea) & confidence %in% c("A","B")) | !is.na(remap), 1,0),
         any = ifelse(physical == 1 | regulatory == 1,1,0))

string_path_wgcna <- do.call(rbind, lapply(1:k, function(clustid){
  df_sd_k100 <- cor_all_long_Lung6 %>% filter(Var1 %in% names(clust_spacedecorr_Lung6_k100$clusters[clust_spacedecorr_Lung6_k100$clusters == clustid]),
                                              Var2 %in% names(clust_spacedecorr_Lung6_k100$clusters[clust_spacedecorr_Lung6_k100$clusters == clustid])) %>%
    summarise(string = mean(stringdp_score),
              pathways = mean(prop_pathways),
              biogrid = mean(biogrid),
              physical = mean(physical),
              regulatory = mean(regulatory),
    ) %>%
    mutate(cluster = clustid,
           model = "spacedecorrk100")
  df_sd_k200 <- cor_all_long_Lung6 %>% filter(Var1 %in% names(clust_spacedecorr_Lung6_k200$clusters[clust_spacedecorr_Lung6_k200$clusters == clustid]),
                                              Var2 %in% names(clust_spacedecorr_Lung6_k200$clusters[clust_spacedecorr_Lung6_k200$clusters == clustid])) %>%
    summarise(string = mean(stringdp_score),
              pathways = mean(prop_pathways),
              biogrid = mean(biogrid),
              physical = mean(physical),
              regulatory = mean(regulatory),
    ) %>%
    mutate(cluster = clustid,
           model = "spacedecorrk200")

  df_sd_k1000 <- cor_all_long_Lung6 %>% filter(Var1 %in% names(clust_spacedecorr_Lung6_k1000$clusters[clust_spacedecorr_Lung6_k1000$clusters == clustid]),
                                               Var2 %in% names(clust_spacedecorr_Lung6_k1000$clusters[clust_spacedecorr_Lung6_k1000$clusters == clustid])) %>%
    summarise(string = mean(stringdp_score),
              pathways = mean(prop_pathways),
              biogrid = mean(biogrid),
              physical = mean(physical),
              regulatory = mean(regulatory),
    ) %>%
    mutate(cluster = clustid,
           model = "spacedecorrk1000")

  df_giotto <- cor_all_long_Lung6 %>% filter(Var1 %in% names(clust_giotto_Lung6$clusters[clust_giotto_Lung6$clusters == clustid]),
                                             Var2 %in% names(clust_giotto_Lung6$clusters[clust_giotto_Lung6$clusters == clustid])) %>%
    summarise(string = mean(stringdp_score),
              pathways = mean(prop_pathways),
              biogrid = mean(biogrid),
              physical = mean(physical),
              regulatory = mean(regulatory)) %>%
    mutate(cluster = clustid,
           model = "Giotto")

  df_vst <- cor_all_long_Lung6 %>% filter(Var1 %in% names(clust_vst_Lung6$clusters[clust_vst_Lung6$clusters == clustid]),
                                          Var2 %in% names(clust_vst_Lung6$clusters[clust_vst_Lung6$clusters == clustid])) %>%
    summarise(string = mean(stringdp_score),
              pathways = mean(prop_pathways),
              biogrid = mean(biogrid),
              physical = mean(physical),
              regulatory = mean(regulatory)) %>%
    mutate(cluster = clustid,
           model = "SPARK")

  df_logpearson <- cor_all_long_Lung6 %>% filter(Var1 %in% names(clust_logpearson_Lung6$clusters[clust_logpearson_Lung6$clusters == clustid]),
                                                 Var2 %in% names(clust_logpearson_Lung6$clusters[clust_logpearson_Lung6$clusters == clustid])) %>%
    summarise(string = mean(stringdp_score),
              pathways = mean(prop_pathways),
              biogrid = mean(biogrid),
              physical = mean(physical),
              regulatory = mean(regulatory)) %>%
    mutate(cluster = clustid,
           model = "logPearson")

  rbind(df_sd_k100,df_sd_k200,df_sd_k1000, df_giotto,df_vst,df_logpearson)
}))
mycol <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3"  )

# Metrics
string_path_wgcna %>% pivot_longer(string:regulatory, names_to = "metric", values_to = "values") %>% group_by(metric, model) %>%
  filter(metric %in% c("pathways","physical")) %>%
  summarise(mean = mean(values),
            median = median(values),
            thirdquartile  = quantile(values,0.75),
            min = min(values),
            max = max(values))

# enrichment
enrich_df_Lung6 <-  rbind(enrich_clust_spacedecorr_Lung6_k100,
                          enrich_clust_spacedecorr_Lung6_k200,
                          enrich_clust_spacedecorr_Lung6_k1000,
                          enrich_clust_giotto_Lung6,
                          enrich_clust_logpearson_Lung6,
                          enrich_clust_vst_Lung6
)

# number of total pathways
enrich_df_Lung6 %>%
  filter(p.adjust < 0.05) %>%
  group_by(model) %>%
  summarise(n_terms = n())

# uniquely enriched pathways
#Shows how many potentially meaningful biological processes are captured.
enrich_df_Lung6 %>%
  filter(p.adjust < 0.05) %>%
  group_by(model) %>%
  summarise(n_terms = n_distinct(ID))

# Plot
plot_modules_coherence <- string_path_wgcna %>%
  mutate(model = factor(model,
                        levels = c("Pearson","Spearman","logPearson",
                                   "SPARK","Giotto","spacedecorrk200"),
                        labels = c("Pearson","Spearman","log-CPM",
                                   "SPARK","Giotto","SpaceDecorr"))) %>%
  pivot_longer(string:regulatory, names_to = "database", values_to = "prop") %>%
  filter(database %in% c("physical"), !is.na(model)) %>%
  ggplot(aes(x = model, y = prop, col = model)) +
  geom_boxplot() +
  facet_wrap(""~.)+
  scale_color_manual(values = mycol) +
  theme_bw() +
  labs(x = "", y = "Proportion of interactions known", col = "")+
  theme(legend.position = "none", strip.background = element_blank(),
        #axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size = textsize),
        axis.text.y = element_text(size = textsize - 1)
  );plot_modules_coherence

plot_modules_coherence_diffk <- string_path_wgcna %>%
  mutate(model = factor(model,
                        levels = c("spacedecorrk100","spacedecorrk200","spacedecorrk1000"),
                        labels = c("k=100","k=200","k=1000"))) %>%
  pivot_longer(string:regulatory, names_to = "database", values_to = "prop") %>%
  filter(database %in% c("physical"), !is.na(model)) %>%
  ggplot(aes(x = model, y = prop, col = model)) +
  geom_boxplot() +
  facet_wrap(""~.)+
  scale_color_manual(values = c("#FF7F00", "#984EA3", "#00A9CE")) +
  theme_bw() +
  labs(x = "", y = "Proportion of interactions known", col = "")+
  theme(legend.position = "none", strip.background = element_blank(),
        #axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size = textsize),
        axis.text.y = element_text(size = textsize - 1)
  );plot_modules_coherence_diffk

### Plot metagene ----
compute_metagenes <- function(expr_mat, gene_clusters) {
  # Ensure matching gene names
  common_genes <- intersect(colnames(expr_mat), names(gene_clusters))
  expr_mat <- expr_mat[, common_genes, drop = FALSE]
  gene_clusters <- gene_clusters[common_genes]

  # Turn clusters into factors
  gene_clusters <- factor(gene_clusters)

  # For each cluster, compute mean expression across its genes
  metagene_matrix <- sapply(levels(gene_clusters), function(cl) {
    genes_in_cluster <- names(gene_clusters)[gene_clusters == cl]
    rowMeans(expr_mat[, genes_in_cluster, drop = FALSE])
  })

  # Ensure output is matrix with proper dimension and names
  if (is.vector(metagene_matrix)) {
    metagene_matrix <- matrix(metagene_matrix, ncol = 1)
  }
  rownames(metagene_matrix) <- rownames(expr_mat)
  colnames(metagene_matrix) <- paste0("cluster", seq_len(ncol(metagene_matrix)))

  return(metagene_matrix)
}

spacedecorr_Lung6_spatial <- as.data.frame(spacedecorr_Lung6_spatial_k200)

metagenes_spacedecorr <- cbind(compute_metagenes(spacedecorr_Lung6_spatial[,colnames(spacedecorr_Lung6_spatial)%in% contamination$target],
                                                 clust_spacedecorr_Lung6_k200$clusters),locs_lung6)
metagenes_giotto <- cbind(compute_metagenes(spacedecorr_Lung6_spatial[,colnames(spacedecorr_Lung6_spatial)%in% contamination$target],
                                            clust_giotto_Lung6$clusters),locs_lung6)
metagenes_logpearson <- cbind(compute_metagenes(spacedecorr_Lung6_spatial[,colnames(spacedecorr_Lung6_spatial)%in% contamination$target],
                                                clust_logpearson_Lung6$clusters),locs_lung6)
metagenes_vst <- cbind(compute_metagenes(spacedecorr_Lung6_spatial[,colnames(spacedecorr_Lung6_spatial)%in% contamination$target],
                                         clust_vst_Lung6$clusters),locs_lung6)

clust_logpearson_Lung6$clusters
# Giotto cluster 16 / SD cluster 14
plot_metagiotto_16 <- metagenes_giotto %>% pivot_longer(cluster1:cluster20, names_to = "cluster", values_to = "Expression") %>%
  mutate(cluster = factor(cluster, levels = paste0("cluster",c(1:20)),
                          labels = paste("Cluster",c(1:20),"(Giotto)"))) %>%
  filter(cluster  %in% paste("Cluster",c(16),"(Giotto)")) %>%
  arrange(Expression) %>%
  filter(!is.na(Expression))%>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=0.5) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_wrap(cluster~.) +#coord_fixed() +
  theme_bw() + labs(col = "Average\nSpatial Effect") +
  guides(col = guide_colorbar(title.vjust = 0.7,
                              barheight = 0.5,
                              barwidth = 8))+
  theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
        legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
        legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
        panel.grid.minor = element_blank()
  );plot_metagiotto_16

# SpaceDecorr clusters 9 and 70
plot_metasd_7_8 <- metagenes_spacedecorr %>% pivot_longer(cluster1:cluster20, names_to = "cluster", values_to = "Expression") %>%
  mutate(cluster = factor(cluster, levels = paste0("cluster",c(1:20)),
                          labels = paste("Cluster",c(1:20),"(SpaceDecorr)"))) %>%
  filter(cluster  %in% paste("Cluster",c(9,10),"(SpaceDecorr)")) %>%
  arrange(Expression) %>%
  filter(!is.na(Expression))%>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=0.5) +
  xlim(31.5,33) + ylim(-17.5, -15.7) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_wrap(cluster~.) +
  theme_bw() +labs(col = "Average\nSpatial Effect") +
  guides(col = guide_colorbar(title.vjust = 0.7,
                              barheight = 0.5,
                              barwidth = 8))+
  theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
        legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
        legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(1, "lines")
  );plot_metasd_7_8

# Giotto cluster 15
# Run FOVs boxes before
plot_metagiotto_15 <- metagenes_giotto %>% pivot_longer(cluster1:cluster20, names_to = "cluster", values_to = "Expression") %>%
  mutate(cluster = factor(cluster, levels = paste0("cluster",c(1:20)),
                          labels = paste("Cluster",c(1:20),"(Giotto)"))) %>%
  filter(cluster  %in% paste("Cluster",c(15),"(Giotto)")) %>%
  arrange(Expression) %>%
  filter(!is.na(Expression))%>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=0.5) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  geom_rect(data = fov_boxes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = NA, color = "black", size = 0.3) +
  geom_text(data = fov_boxes, aes(x = center_x, y = center_y, label = fov),
            inherit.aes = FALSE,
            size = 3.5, color = "black", fontface = "bold") +
  facet_wrap(cluster~.) +
  theme_bw() + labs(col = "Average\nSpatial Effect") +
  guides(col = guide_colorbar(title.vjust = 0.7,
                              barheight = 0.5,
                              barwidth = 8))+
  theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
        legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
        legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
        panel.grid.minor = element_blank()
  );plot_metagiotto_15

# all Giotto clusters
plot_metagiotto_all <- metagenes_giotto %>% pivot_longer(cluster1:cluster20, names_to = "cluster", values_to = "Expression") %>%
  mutate(cluster = factor(cluster, levels = paste0("cluster",c(1:20)),
                          labels = paste("Cluster",c(1:20)))) %>%
  arrange(Expression) %>%
  filter(!is.na(Expression))%>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=0.5) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_wrap(cluster~.) +
  theme_bw() + labs(col = "Average\nSpatial Effect") +
  theme(strip.background = element_rect(fill="white", colour = "white"),
        text = element_text(size = 15)) + coord_fixed()

# all SpaceDecorr clusters
plot_metasd_all <- metagenes_spacedecorr %>% pivot_longer(cluster1:cluster20, names_to = "cluster", values_to = "Expression") %>%
  mutate(cluster = factor(cluster, levels = paste0("cluster",c(1:20)),
                          labels = paste("Cluster",c(1:20)))) %>%
  arrange(Expression) %>%
  filter(!is.na(Expression))%>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=0.5) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_wrap(cluster~.) +
  theme_bw() + labs(col = "Average\nSpatial Effect") +
  theme(strip.background = element_rect(fill="white", colour = "white"),
        text = element_text(size = 15)) + coord_fixed()

# Batch effect
cl_sd <- clust_spacedecorr_Lung6_k200$clusters
cl_gt <- clust_giotto_Lung6$clusters


### Graphs technical artifact ---
# Technical artifact
clust_gt <- 15
graph_gt15 <- cor_all_long_Lung6 %>%
  filter(Var1 %in% names(cl_gt[cl_gt == clust_gt] ),
         Var2 %in% names(cl_gt[cl_gt == clust_gt] )) %>%
  dplyr::select(Var1, Var2, independent,physical,prop_pathways,regulatory) %>%
  filter(physical == 1 | prop_pathways > 0.01|regulatory == 1) %>%
  graph_from_data_frame(., directed = F, vertices =names(cl_gt[cl_gt == clust_gt] ) )

graph_sd17 <- cor_all_long_Lung6 %>%
  filter(Var1 %in% names(cl_sd[cl_sd == 17] ),
         Var2 %in% names(cl_sd[cl_sd == 17] )) %>%
  dplyr::select(Var1, Var2, independent,physical,prop_pathways,regulatory) %>%
  filter(physical == 1 | prop_pathways > 0.01) %>%
  graph_from_data_frame(., directed = F, vertices =names(cl_sd[cl_sd == 17] ) )

library(ggraph)
library(tidygraph)
plot_graph_gt15 <- ggraph(graph_gt15, layout = "fr") +
  geom_edge_link(alpha = 0.4, color = "gray50") +
  geom_node_point(size = 2, aes(color = name == "MZT2A"), show.legend = FALSE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  theme_void() +
  scale_color_manual(values = c("black","#B2182B"))

plot_graph_sd17 <- ggraph(graph_sd17, layout = "fr") +
  geom_edge_link(alpha = 0.4, color = "gray50") +
  geom_node_point(size = 2, aes(color = name == "MZT2A"), show.legend = FALSE) +
  geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  theme_void()+
  scale_color_manual(values = c("black","#B2182B"))

plot_grid(plot_graph_gt15,plot_graph_sd17 )

### Compare clusters k ----
cl_sd100 <- clust_spacedecorr_Lung6_k100$clusters
cl_sd200 <- clust_spacedecorr_Lung6_k200$clusters
cl_sd1000 <- clust_spacedecorr_Lung6_k1000$clusters

plot_comp_clusters_diffk <- list(`100`=cl_sd100, `200`=cl_sd200, `1000`=cl_sd1000) %>%
  imap_dfr(~tibble(gene=names(.x), otherk=.y, cl=.x)) %>%
  pivot_wider(id_cols=gene, names_from=otherk, values_from=cl) %>%
  pivot_longer(cols=c(`100`,`1000`), names_to="otherk", values_to="cl_x") %>%
  count(otherk, `200`, cl_x, name="n") %>%
  group_by(`200`) %>% mutate(p = n / sum(n)) %>% ungroup() %>%
  complete(otherk, `200`, cl_x, fill=list(n=0, p=0)) %>%
  mutate(otherk = factor(otherk, levels = c(100,1000), labels = c("k=100","k=1000"))) %>%
  ggplot(aes(factor(cl_x), factor(`200`), fill=p)) +
  geom_tile() +
  scale_fill_gradient(limits=c(0,1)) +
  facet_wrap(~otherk) +
  scale_fill_gradient(low = "white", high = "black") +
  theme_bw() +
  labs(x = "Clusters (k=100 or 1000)", y = "Clusters k=200",
       fill = "Proportion in common\nwith k=200") +
  theme(strip.background = element_blank())

### Save plot ----
plot_densities <- plot_grid(legend = get_plot_component(plot_density_indep, 'guide-box-top', return_all = TRUE),
                            plot_grid(plot_density_indep + theme(legend.position = "none"),
                                      plot_density_dep+ theme(legend.position = "none"),
                                      plot_modules_coherence+coord_flip(),
                                      nrow=1, align = "v", axis = "tb",
                                      labels = c("C","D","E")),
                            ncol = 1, rel_heights =  c(0.15,0.85));plot_densities

plot_grid(
  plot_grid(plot_cormap, plot_spatialeffect, ncol = 2, rel_widths = c(1, 2),
            labels = c("A","B")),
  plot_densities,
  plot_grid(plot_metagiotto_16,
            ggplot() + theme_void(),
            plot_grid(plot_graph_gt15,plot_graph_sd17,ncol=1,labels = c("H","J") ),
            plot_metagiotto_15, nrow = 1, rel_widths = c(1,.5,.5,1),
            labels = c("F","G","","I")),
  ncol = 1
)

plot_densities_diffk <- plot_grid(legend = get_plot_component(plot_density_indep_diffk, 'guide-box-top', return_all = TRUE),
                                  plot_grid(plot_density_indep_diffk + theme(legend.position = "none"),
                                            plot_density_dep_diffk+ theme(legend.position = "none"),
                                            plot_modules_coherence_diffk+coord_flip(),
                                            nrow=1, align = "v", axis = "tb",
                                            labels = c("A","B","C")),
                                  ncol = 1, rel_heights =  c(0.15,0.85))
### Save tables ----
modules_wgcna <- merge(as.data.frame(clust_logpearson_Lung6$clusters),
                       as.data.frame(clust_vst_Lung6$clusters),
                       by = "row.names"
) %>%
  merge(merge(as.data.frame(clust_giotto_Lung6$clusters),
              as.data.frame(clust_spacedecorr_Lung6_k200$clusters),
              by = "row.names"
  ),  by = "Row.names")
