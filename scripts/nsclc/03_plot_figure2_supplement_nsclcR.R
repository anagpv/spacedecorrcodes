#--------------------------------------
#---- Plot figure 2 and supplement
#--------------------------------------

dir.data <- "/path/to/zenodo/download/"  
dir.out  <- "/path/to/output/"          

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
textsize <- 12


# Read data -----
# Raw and normalized data
load(file.path(dir.data, "Lung_Giotto_Filtered_Tumor.RData"))
id_tissue <- pDataDT(gem)$cell_ID[pDataDT(gem)$tissue == "Lung6"]
gem_tissue <- subsetGiotto(gem, cell_ids = id_tissue)

locs <- gem_tissue@spatial_locs$cell$raw@coordinates
rownames(locs) <- locs$cell_ID
normexpr <- as.data.frame(as.matrix(t(gem_tissue@expression$cell$rna$normalized@exprMat)))
normexpr$cell_ID <- rownames(normexpr)
normexpr <- merge(normexpr,locs, by = "cell_ID" )

# Correlations
cor_all_long <- read.csv(file = file.path(dir.out, "cor_all_long_lung.csv"))

test1cor <- function(r,n){
  fisher <- r / sqrt(1 - r^2) * sqrt(n - 2)
  pvalue <- 2*pt(-abs(fisher), n - 2)
  return(pvalue)
}

pvalue_all <- cor_all_long %>%
  filter(as.character(Var1) < as.character(Var2)) %>%   # each pair once
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  mutate(pvalue = test1cor(cor, 66267)) %>%
  group_by(model) %>%
  mutate(padj = p.adjust(pvalue, method = "BH"))


# MZT2A and DDR1 ---------------------------------------------------
pvalue_all %>%
  filter(independent) %>%
  group_by(model) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  filter(Var2 == "MZT2A", Var1 == "DDR1") %>%
  mutate(sign = as.numeric(padj < 0.05)) %>%
  dplyr::select(Var1,Var2,stringdp_score,pathway_distance, model, cor, pvalue, padj,sign ) %>%
  mutate_at(6, ~round(.x,3))



# Plot Random effects (1a) ------------------------------------------------
# Raw counts and metadata (Lung6 tumor cells)
raw  <- gem_tissue@expression$cell$rna$raw@exprMat
meta <- as.data.frame(pDataDT(gem_tissue, spat_unit = "cell"))

# k 200
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k200.Rdata"))
spacedecorr_spatial_k200 <- cbind(data_spacedecorr$splines,pDataDT(gem_tissue, spat_unit = "cell"),gem_tissue@spatial_locs$cell$raw@coordinates)

plot_spatialeffect <- spacedecorr_spatial_k200 %>%
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

# Different basis dimension
# k 100
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k100.Rdata"))
spacedecorr_spatial_k100 <- cbind(data_spacedecorr$splines,pDataDT(gem_tissue, spat_unit = "cell"),gem_tissue@spatial_locs$cell$raw@coordinates)

# k 1000
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorr_k1000.Rdata"))
spacedecorr_spatial_k1000 <- cbind(data_spacedecorr$splines,pDataDT(gem_tissue, spat_unit = "cell"),gem_tissue@spatial_locs$cell$raw@coordinates)

# k 200 spots
load(paste0(dir.out,"Lung6/Data_Lung6_spacedecorrspots_k200.Rdata"))
spacedecorr_spatial_k200spots <- cbind(data_spacedecorr$splines,data_spacedecorr$metadata)

plot_spatialeffect_diffk <- spacedecorr_spatial_k1000 %>%
  dplyr::select(sdimx, sdimy,MZT2A,DDR1 ) %>%
  pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
  mutate(k = "k=1000") %>%
  rbind(spacedecorr_spatial_k200 %>%
          dplyr::select(sdimx, sdimy,MZT2A,DDR1 ) %>%
          pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
          mutate(k = "k=200"))%>%
  rbind(spacedecorr_spatial_k100 %>%
          dplyr::select(sdimx, sdimy,MZT2A,DDR1 ) %>%
          pivot_longer(all_of(c("MZT2A","DDR1")), names_to = "gene", values_to = "Expression") %>%
          mutate(k = "k=100"))%>%
  mutate(k = factor(k, levels = c("k=100","k=200","k=200 (spots)","k=1000"))) %>%
  filter(!is.na(Expression))%>%
  arrange(Expression) %>%
  ggplot(aes(x = sdimx, y=sdimy, col = Expression)) +
  geom_point(cex=.5) +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1) +
  facet_grid(gene~k) +
  theme_bw() + labs(col = "Spatial Effect")+
  guides(col = guide_colorbar(title.vjust = 0.7,
                              barheight = 0.5,
                              barwidth = 8))+
  theme(legend.position = "top",strip.background = element_rect(fill="white", colour = "white"),
        legend.box.spacing = unit(0, "pt"),text = element_text(size = textsize),
        legend.title = element_text(margin = margin(r = 20, b=10),size = textsize-2),
        panel.grid.minor = element_blank(),
  );plot_spatialeffect_diffk
ggsave(paste0(dir.out,"Images/Random_effects_diffk.jpeg"), height = 150, width = 250, units = "mm")


# Plot densities (1 b-c) --------------------------------------------------
mycol <- c("SpaceDecorr"         = "#8C510A",
           "spMOCA"              = "#FF7F00",
           "log-normalized"      = "#000000",
           "VST-normalized"      = "#7F7F7F",
           "Giotto"              = "#6BAED6",
           "MERINGUE"            = "#08519C",
           "GraphR"              = "#1B7837",
           "SpaceDecorr (spots)" = "#e377c2")

# independent
pvalue_all %>% ungroup %>%
  mutate(independent = ifelse(pathway_distance == 1 &
                                stringdp_score < 100 &
                                biogrid == 0 &
                                is.na(trrust) & is.na(dorothea) & is.na(remap) & is.na(n_shared_tf) , T, F)) %>%
  filter(independent) %>%
  group_by(model) %>%
  mutate(t1err = as.numeric(padj < 0.05)) %>%
  group_by(model) %>%
  summarise(t1err = sprintf("%.1f%%",100*mean(t1err)))

plot_density_indep <- cor_all_long %>%
  filter(independent) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("logpearson","vst","Giotto","spacedecorr200", "spmoca")) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorr200", "spmoca", "logpearson","vst","Giotto"),
                        labels = c("SpaceDecorr","spMOCA","log-normalized","VST-normalized","Giotto"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",linewidth = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("Negative Control set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 2) +
  xlim(-0.25,0.25) +
  scale_color_manual(values = mycol) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        # legend.text = element_text(size = textsize+1),
        # strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

# dependent
plot_density_dep <- cor_all_long %>%
  filter(prop_pathways > 0.15, stringdp_score > 700, biogrid == 1) %>% 
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("logpearson","vst","Giotto","spacedecorr200", "spmoca")) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorr200", "spmoca", "logpearson","vst","Giotto"),
                        labels = c("SpaceDecorr","spMOCA","log-normalized","VST-normalized","Giotto"))) %>%
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
        # strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank());plot_density_dep

# Biological coherence (1d) -----------------------------------------------
bio_coherence_modules <- read.csv(file = file.path(dir.out, "bio_coherence_modules.csv"))

# Metrics
bio_coherence_modules %>% pivot_longer(pathways:any, names_to = "metric", values_to = "values") %>% group_by(metric, model) %>%
  filter(metric %in% c("regulatory"), k == 20) %>%
  summarise(mean = mean(values),
            median = median(values),
            thirdquartile  = quantile(values,0.75),
            min = min(values),
            max = max(values)) %>% View

# Plot
# By module number
plot_coherence_k <- (bio_coherence_modules %>%
                       mutate(model = factor(model,
                                             levels = c("spacedecorrk200","spMOCA", "logPearson","SPARK","Giotto"),
                                             labels = c("SpaceDecorr","spMOCA","log-normalized","VST-normalized","Giotto"))) %>%
                       filter(!is.na(model)) %>%
                       pivot_longer(pathways:any, names_to = "database", values_to = "prop") %>%
                       filter(database %in% c("physical","regulatory","any") ) %>%
                       mutate(database = factor(database, levels = c("physical","regulatory","any"),
                                                labels = c("Physical or Functional","Regulatory","Any"))) %>%
                       ggplot(aes(x = k, y = prop, col = model)) +
                       stat_summary(geom = "point")+
                       stat_summary(geom = "line")+
                       scale_color_manual(values = mycol) +
                       scale_x_continuous(breaks=c(5,10,20,30))+
                       theme_bw() +
                       facet_wrap(~database, scales = "free") +
                       labs(x = "Number of modules", y = "Average proportion of interactions known", col = "")+
                       theme(legend.position = "top", strip.background = element_blank(),
                             text = element_text(size = textsize),
                             axis.text.y = element_text(size = textsize - 1)
                       ))


# Physical or functional
plot_modules_coherence <- bio_coherence_modules %>%
  filter( k == 20 ) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorrk200","spMOCA", "logPearson","SPARK","Giotto"),
                        labels = c("SpaceDecorr","spMOCA","log-normalized","VST-normalized","Giotto"))) %>%
  pivot_longer(pathways:any, names_to = "database", values_to = "prop") %>%
  filter(database %in% c("physical"), !is.na(model)) %>%
  ggplot(aes(x = model, y = prop, col = model)) +
  geom_boxplot() +
  facet_wrap(""~.)+
  facet_wrap(~database, scales = "free") +
  scale_color_manual(values = mycol) +
  theme_bw() +
  labs(x = "", y = "Proportion of physical/functional\ninteractions known", col = "")+
  theme(legend.position = "none", strip.background = element_blank(),
        #axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size = textsize),
        axis.text.y = element_text(size = textsize - 1)
  ) + coord_flip();plot_modules_coherence

# Regulatory
bio_coherence_modules %>%
  filter( k == 20 ) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorrk200","spMOCA", "logPearson","SPARK","Giotto"),
                        labels = c("SpaceDecorr","spMOCA","log-normalized","VST-normalized","Giotto"))) %>%
  pivot_longer(pathways:any, names_to = "database", values_to = "prop") %>%
  filter(database %in% c("regulatory"), !is.na(model)) %>%
  ggplot(aes(x = model, y = prop, col = model)) +
  geom_boxplot() +
  facet_wrap(""~.)+
  # facet_wrap(~database, scales = "free") +
  scale_color_manual(values = mycol) +
  theme_bw() +
  labs(x = "", y = "Proportion of regulatory/shared TF\nininteractions known", col = "")+
  theme(legend.position = "none", strip.background = element_blank(),
        #axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size = textsize),
        axis.text.y = element_text(size = textsize - 1)
  ) 
ggsave(paste0(dir.out,"Images/Biological_coherence_regulatory.jpeg"), height = 150, width = 200, units = "mm")


## Pathway enrichment ----
# enrichment
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

# Number of modules with at least one pathway enriched
plot_numodenrich_k <- enrich_df %>% filter(Count > 2) %>%
  filter(p.adjust < 0.05) %>%
  dplyr::count(model, clustID, k) %>%
  dplyr::count(model, k) %>%
  mutate(model = factor(model,
                        levels = c("logPearson","vst","giotto","spmoca","spacedecorr k=200"),
                        labels = c("log-normalized","VST-normalized","Giotto","spMOCA","SpaceDecorr"))) %>%
  filter(!is.na(model)) %>%
  ggplot(aes(x = k, y = n/k, col = model)) +
  stat_summary(geom = "point")+
  stat_summary(geom = "line")+
  scale_color_manual(values = mycol) +
  scale_x_continuous(breaks=c(5,10,20,30))+
  theme_bw() +
  facet_wrap(~"", scales = "free") +
  labs(x = "Number of modules", y = "Proportion of modules with at least one\npathway enriched", col = "")+
  theme(legend.position = "top", strip.background = element_blank(),
        #axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size = textsize),
        axis.text.y = element_text(size = textsize - 1)
  )
plot_grid(legend = get_plot_component(plot_coherence_k, 'guide-box-top', return_all = TRUE),
          plot_grid(plot_coherence_k + theme(legend.position = "none"),
                    plot_numodenrich_k + theme(legend.position = "none"),
                    labels = c("",""), rel_widths = c(3,1)),
          ncol = 1, rel_heights =  c(0.15,0.85))
ggsave(paste0(dir.out,"Images/Biological_coherence_varykmodules.jpeg"), height = 150, width = 300, units = "mm")


# Module decomposition ----------------------------------------------------
mod_cat <- read.csv(file.path(dir.out, "module_decomposition.csv"))

plot_broad_scale <-mod_cat %>% filter(model %in% c("spacedecorr","giotto", "spark")) %>%
  mutate(label = paste0(ifelse(model == "spacedecorr", "S", ifelse(model == "spark","V","G")), module),
         model = factor(model, levels = c("spacedecorr","spark","giotto"),
                        labels = c("SpaceDecorr","VST-normalized","Giotto"))) %>%
  ggplot(aes(broad_share, nb_share, colour = model)) +
  geom_point() +
  geom_text_repel(aes(label = module), size = 3, show.legend = FALSE, min.segment.length = 0,max.overlaps = Inf) +
  labs(x = "Broad-scale share of co-expression",
       y = "Neighbourhood share\nof fine-scale variation",
       colour = NULL) +
  facet_wrap(model~.) +
  theme_bw() +
  scale_color_manual(values = mycol)+
  theme(strip.background = element_blank(),
        legend.position = "none")

# Alluvial ----------------------------------------------------------------
library(ggalluvial)
load(paste0(dir.out,"ExampleWGCNA_k20.Rdata"))

# VST modules to show, ordered from preserved to changed
vst_lab <- c("10" = "Mitosis\n(M10)",
             "11" = "S phase\n(M11)",
             "6"  = "Squamous\nkeratins\n(M6)",
             "17" = "Hypoxia\n(M17)",
             "8"  = "Antigen\npresentation\n+ epithelial\n(M8)",
             "14" = "Interferon\n(M14)")

# SpaceDecorr modules, same order; regulators in parentheses
sd_lab <- c("9"  = "Mitosis\nFOXM1, MYBL2\n(M9)",
            "10" = "S phase\nE2F1\n(M10)",
            "14" = "Squamous\nkeratins\n(M14)",
            "6"  = "Hypoxic\nbasal\n(M6)",
            "7"  = "Hypoxia/lipid\n(M7)",
            "8"  = "MHC class I\n(M8)",
            "11" = "Epithelial\nadhesion\n(M11)",
            "3"  = "Growth-factor\nsignaling\n(M3)",
            "4"  = "Large mixed\nmodule\n(M4)")

# colours shared with the K alluvial (interferon added)
vst_cols <- setNames(c("#2166AC", "#67A9CF", "#1B7837", "#E08214", "#B2182B", "#7B3294"), vst_lab)

# one row per VST module x SpaceDecorr module
flow <- data.frame(gene = names(clust_vst$clusters),
                   vst  = as.character(clust_vst$clusters)) %>%
  filter(vst %in% names(vst_lab)) %>%
  mutate(sd  = as.character(clust_spacedecorr_k200$clusters[gene]),
         vst = factor(vst_lab[vst], levels = vst_lab),
         sd  = ifelse(sd %in% names(sd_lab), sd_lab[sd], "Other"),
         sd  = factor(sd, levels = c(sd_lab, "Other"))) %>%
  dplyr::count(vst, sd) %>%
  mutate(other = sd == "Other")

w <- 0.3   # same width for strata and flows
plot_flow_spark <- ggplot(flow, aes(axis1 = vst, axis2 = sd, y = n)) +
  geom_alluvium(aes(fill = vst, alpha = other), width = w, knot.pos = 0.35) +
  geom_stratum(width = w, fill = "white", colour = "grey40", linewidth = 0.2) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, lineheight = 0.85) +
  scale_x_discrete(limits = c("VST", "SpaceDecorr"), expand = c(0.15, 0.15)) +
  scale_fill_manual(values = vst_cols) +
  scale_alpha_manual(values = c("FALSE" = 0.85, "TRUE" = 0.3)) +
  labs(x = NULL, y = "Number of genes") +
  theme_bw() +
  theme(legend.position = "none", panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = textsize, face = "bold"),
        text = element_text(size = textsize))
plot_flow_spark

### Save plot ----
plot_densities <- plot_grid(legend = get_plot_component(plot_density_indep, 'guide-box-top', return_all = TRUE),
                            plot_grid(plot_density_indep + theme(legend.position = "none"),
                                      plot_density_dep+ theme(legend.position = "none"),
                                      nrow=1, align = "h", axis = "tb",
                                      labels = c("B","C","D")),
                            ncol = 1, rel_heights =  c(0.15,0.85));plot_densities
plot_grid(
  plot_grid(
    plot_grid(plot_spatialeffect, plot_densities, ncol = 2, rel_widths = c(1, 1.5),
              labels = c("A","")),
    plot_grid(plot_modules_coherence, 
              plot_broad_scale,
              nrow = 1, rel_widths = c(2,3),align = "h", axis = "tb",
              labels = c("D","E")),
    ncol = 1, rel_heights =  c(0.5,0.6)
  ), plot_flow_spark, rel_widths = c(3,1.5), labels = c("","F"))
ggsave(paste0(dir.out,"Images/PlotsExampleRealDataWGCNA_review.pdf"), height = 230, width = 450, units = "mm", dpi = 400)



# Supplement -------------------------------------------------------------
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

# Modified t-test ---------------------------------------------------------
load(file = file.path(dir.out, "res_modified_ttest_Lung6.Rdata"))

# ESS and significant pairs per method
res_k %>%
  group_by(method) %>%
  mutate(ess = ess/66267) %>%
  summarise(n_cells = 66267, n_pairs = n(),
            ess_median = median(ess), ess_q10 = quantile(ess, 0.1), ess_q90 = quantile(ess, 0.9),
            sig_naive = sum(fdr_naive < 0.05), sig_modified = sum(fdr_modified < 0.05),
            .groups = "drop")

mycol_ess <- c("SpaceDecorr" = "#8C510A", "VST-normalized" = "#7F7F7F", "log-normalized" = "#000000")
plot_df <- res_nsclc %>%
  mutate(method = factor(method, levels = c("spacedecorr","vst","lognorm"),
                         labels = c("SpaceDecorr","VST-normalized","log-normalized")),
         ess_prop = ess / 66267,
         lp_naive    = -log10(pmax(fdr_naive, 1e-300)),      # avoid Inf
         lp_modified = -log10(pmax(fdr_modified, 1e-300)))

# (A) effective sample size as a proportion of cells
p_ess <- plot_df %>%
  ggplot(aes(method, ess_prop, col = method)) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_boxplot(outlier.size = 0.3) +
  scale_color_manual(values = mycol) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Effective sample size (% of cells)") +
  theme_bw() +
  theme(legend.position = "none")

# (B) standard vs modified t-test, -log10 FDR
p_scatter <- plot_df %>%
  ggplot(aes(lp_naive, lp_modified, col = method)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_abline(lty = 2, col = "grey40") +
  geom_hline(yintercept = -log10(0.05), lty = 3, col = "grey60") +
  geom_vline(xintercept = -log10(0.05), lty = 3, col = "grey60") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_color_manual(values = mycol) +
  facet_wrap(~ method) +
  labs(x = expression(-log[10]~"FDR, standard t-test"),
       y = expression(-log[10]~"FDR, modified t-test")) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank())

plot_grid(p_ess, p_scatter, labels = c("A", "B"), nrow = 1, rel_widths = c(1.5, 2.5), align = "h", axis = "tb")
ggsave(paste0(dir.out,"Images/Modified_ESS_lung.jpeg"), height = 120, width = 300, units = "mm")


# negative controls significant, standard vs modified test ------------------------
neg_pairs <- cor_all_long %>%
  filter(independent) %>%
  transmute(gene1 = as.character(Var1), gene2 = as.character(Var2)) %>%
  distinct()

res_k %>%
  semi_join(neg_pairs, by = c("gene1", "gene2")) %>%
  group_by(method) %>%
  summarise(n = n(),
            sig_naive    = 100 * mean(fdr_naive < 0.05),
            sig_modified = 100 * mean(fdr_modified < 0.05), .groups = "drop")

# Plot correlation matrix and modules (supp) ------------------------------
load(paste0(dir.out,"ExampleWGCNA_k20.Rdata"))
heatmapclust <- function (clustobj,cor_matrix,use_clus_name = NULL, show_cluster_annot = TRUE,
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
upper <- function(cormat)cormat[upper.tri(cormat)]

# Extract RdBu colors
rdbu <- brewer.pal(9, "RdBu")
red_side  <- rdbu[1:5]  # red to white
blue_side <- rdbu[5:9]  # white to blue

# Choose flat red and blue shades
flat_red  <- red_side[1]
flat_blue <- blue_side[5]

# Define custom color scale
heatmap_colors_sd100 <- colorRamp2(
  c(min(cormat_spacedecorr_k100),  0, max(upper(cormat_spacedecorr_k100))),
  c(flat_blue,  "white", flat_red)
)
plot_wgcna_spacedecorr_k100 <- heatmapclust(clust_spacedecorr_k100,cormat_spacedecorr_k100,
                                            column_title = "k=100",col=heatmap_colors_sd100,
                                            show_heatmap_legend = F,show_clust_legend = F)

heatmap_colors_sd200 <- colorRamp2(
  c(min(cormat_spacedecorr_k200),  0, max(upper(cormat_spacedecorr_k200))),
  c(flat_blue,  "white", flat_red)
)
plot_wgcna_spacedecorr_k200 <- heatmapclust(clust_spacedecorr_k200,cormat_spacedecorr_k200,
                                            column_title = "k=200",col=heatmap_colors_sd200,
                                            show_heatmap_legend = F,show_clust_legend = F)

heatmap_colors_sd1k <- colorRamp2(
  c(min(cormat_spacedecorr_k1000),  0, max(upper(cormat_spacedecorr_k200))),
  c(flat_blue,  "white", flat_red)
)
plot_wgcna_spacedecorr_k1000 <- heatmapclust(clust_spacedecorr_k1000,cormat_spacedecorr_k1000,
                                             column_title = "k=1000",col=heatmap_colors_sd1k,
                                             show_heatmap_legend = T,show_clust_legend = F)

heatmap_colors_giotto <- colorRamp2(
  c(min(cormat_giotto),  0, max(upper(cormat_giotto))),
  c(flat_blue,  "white", flat_red)
)

plot_wgcna_giotto <- heatmapclust(clust_giotto,cormat_giotto,
                                  column_title = "Giotto",col=heatmap_colors_giotto,
                                  show_heatmap_legend = T,
                                  show_clust_legend = F);plot_wgcna_giotto

cluster_legend <- Legend(
  labels = sort(unique(clust_spacedecorr_k1000[["clusters"]])),
  title = "Cluster",
  legend_gp = gpar(fill = getDistinctColors(length(unique(clust_spacedecorr_k1000[["clusters"]])))),
  ncol = 2
)

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


# Pseudopots sensitivity -------------------------------------------------------
plot_density_indep_spots <- cor_all_long %>%
  filter(independent) %>%
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("spacedecorr200spots","spacedecorr200","spmoca")) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorr200","spacedecorr200spots","spmoca"),
                        labels = c("SpaceDecorr","SpaceDecorr (spots)","spMOCA"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",size = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("Negative Control set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 1,alpha = 0.7) +
  xlim(-0.25,0.25) +
  scale_color_manual(values = my_col) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        legend.text = element_text(size = textsize+1),
        strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

plot_density_dep_spots <- cor_all_long %>%
  filter(prop_pathways > 0.15, stringdp_score > 800, biogrid == 1) %>% 
  pivot_longer(spacedecorr200:vst, names_to = "model", values_to = "cor") %>%
  filter(model %in% c("spacedecorr200spots","spacedecorr200","spmoca")) %>%
  mutate(model = factor(model,
                        levels = c("spacedecorr200","spacedecorr200spots","spmoca"),
                        labels = c("SpaceDecorr","SpaceDecorr (spots)","spMOCA"))) %>%
  ggplot(aes(x = cor, col = model)) +
  stat_density(geom="line",position="identity",size = 0.8,alpha = 1,adjust = 3) +
  facet_wrap("High confidence set"~.)+
  geom_vline(xintercept = 0, col = "grey", lty = 1,alpha=0.7) +
  xlim(-1,1) +
  scale_color_manual(values = mycol) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        strip.text = element_text(size = textsize+1),
        
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())


plot_modules_coherence_spots <- bio_coherence_modules %>%
  mutate(model = factor(model,
                        levels = c("spacedecorrk200","spacedecorrk200spots","spMOCA"),
                        labels = c("SpaceDecorr","SpaceDecorr\n(spots)","spMOCA"))) %>%
  pivot_longer(pathways:any, names_to = "database", values_to = "prop") %>%
  filter(database %in% c("physical"), !is.na(model), k == 20) %>%
  ggplot(aes(x = model, y = prop, col = model)) +
  geom_boxplot(aes(group = (model))) +
  facet_wrap(""~.)+
  scale_color_manual(values = c("#8C510A", "#e377c2","#FF7F00" )) +
  theme_bw() +
  labs(x = "", y = "Proportion of physical/functional\ninteractions known", col = "")+
  theme(legend.position = "none", strip.background = element_blank(),
        #axis.text.x = element_text(angle=45, vjust=1, hjust=1),
        text = element_text(size = textsize),
        # axis.text.y = element_text(size = textsize - 1)
  )+coord_flip();plot_modules_coherence_spots

plot_grid(legend = get_plot_component(plot_density_indep_spots, 'guide-box-top', return_all = TRUE),
          plot_grid(plot_density_indep_spots + theme(legend.position = "none"),
                    plot_density_dep_spots+ theme(legend.position = "none"),
                    plot_modules_coherence_spots,
                    labels = c("A","B","C"), nrow =1, align = "h", axis = "tb"),
          ncol = 1, rel_heights =  c(0.15,0.85))
ggsave(paste0(dir.out,"Images/Densities_spots.jpeg"), height = 120, width = 300, units = "mm")


# Basis dimension sensitivity ------------------------------------------------------
# Negative control genes
plot_density_indep_diffk <- cor_all_long %>%
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
  scale_color_manual(values = c("#FF7F00" ,"#984EA3","#00A9CE","pink")) +
  labs(x = "Estimated correlation", y = "Density", col = "") +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top",legend.box = "horizontal",
        text = element_text(size = textsize),
        legend.text = element_text(size = textsize+1),
        strip.text = element_text(size = textsize+1),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())

# High confidence genes
plot_density_dep_diffk <- cor_all_long %>%
  filter(prop_pathways > 0.15, stringdp_score > 700, biogrid == 1) %>%
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

# module support by K, physical and regulatory, with VST and log-normalized as reference
plot_modules_coherence_diffk <- bio_coherence_modules %>%
  mutate(model = factor(model,
                        levels = c("spacedecorrk100","spacedecorrk200","spacedecorrk1000"),
                        labels = c("k=100","k=200","k=1000"))) %>%
  pivot_longer(pathways:any, names_to = "database", values_to = "prop") %>%
  filter(database %in% c("physical","regulatory"), !is.na(model), k == 20) %>%
  mutate(database = factor(database, c("physical","regulatory"), c("Physical or functional","Regulatory"))) %>%
  ggplot(aes(x = model, y = prop, col = model)) +
  geom_boxplot(aes(group = interaction(model, k))) +
  facet_wrap(database~.)+
  scale_color_manual(values = c("#FF7F00", "#984EA3", "#00A9CE")) +
  theme_bw() +
  labs(x = "", y = "Proportion of interactions known", col = "")+
  theme(legend.position = "none", strip.background = element_blank(),
        text = element_text(size = textsize),
        axis.text.y = element_text(size = textsize - 1)
  ) + coord_flip();plot_modules_coherence_diffk


# alluvial ----
load(file.path(dir.out, paste0("ExampleWGCNA_k",20,".Rdata")))
cl <- list(k100  = clust_spacedecorr_k100$clusters,
           k200  = clust_spacedecorr_k200$clusters,
           k1000 = clust_spacedecorr_k1000$clusters,
           vst   = clust_vst$clusters,
           log   = clust_logpearson$clusters)
g <- Reduce(intersect, lapply(cl, names))
hl_ord <- c(8, 9, 10, 14, 6)
hl_lab <- c("MHC class I (M8)", "Mitotic (M9)", "S-phase (M10)", "Keratin (M14)", "Hypoxia-basal (M6)")
epi    <- c("EPCAM", "KRT8", "KRT18", "KRT19")   # grouped with MHC class I at k=100

al <- data.frame(gene  = g,
                 k100  = cl$k100[g],
                 k200  = cl$k200[g],
                 k1000 = cl$k1000[g]) %>%
  mutate(module = case_when(k200 %in% hl_ord ~ hl_lab[match(k200, hl_ord)],
                            gene %in% epi ~ "Epithelial markers",
                            TRUE ~ "Other"),
         module = factor(module, c(hl_lab, "Epithelial markers", "Other")))

# keep the highlighted genes and every gene that shares a K = 100 or K = 1000 module with them
hl <- !al$module %in% c("Other", "Epithelial markers")
al <- al %>% filter(hl | k100 %in% k100[hl] | k1000 %in% k1000[hl])


# order strata to reduce crossings: K = 200 with the highlighted modules on top,
# K = 100 and K = 1000 by the median K = 200 position of their genes
lev200   <- c(hl_ord, setdiff(sort(unique(al$k200)), hl_ord))
al$k200  <- factor(al$k200, levels = lev200)
pos      <- as.integer(al$k200)
al$k100  <- factor(al$k100,  levels = names(sort(tapply(pos, al$k100,  median))))
al$k1000 <- factor(al$k1000, levels = names(sort(tapply(pos, al$k1000, median))))
al <- al %>% arrange(desc(module))           # highlighted flows drawn on top

mod_cols <- c("MHC class I (M8)" = "#B2182B", "Mitotic (M9)" = "#2166AC", "S-phase (M10)" = "#67A9CF",
              "Keratin (M14)" = "#1B7837", "Hypoxia-basal (M6)" = "#E08214",
              "Epithelial markers" = "#4D4D4D", "Other" = "grey70")

plot_alluvial_k <- ggplot(al, aes(axis1 = k100, axis2 = k200, axis3 = k1000)) +
  geom_alluvium(aes(fill = module, alpha = module == "Other"), width = 1/5, knot.pos = 0.35) +
  geom_stratum(width = 1/5, fill = "white", colour = "grey40", linewidth = 0.2) +
  geom_text(stat = "stratum", aes(label = paste0("M", after_stat(stratum))), size = 2.3) +
  scale_x_discrete(limits = c("K = 100", "K = 200", "K = 1000"), expand = c(0.1, 0.1)) +
  scale_fill_manual(values = mod_cols, name = NULL) +
  scale_alpha_manual(values = c("TRUE" = 0.35, "FALSE" = 0.9), guide = "none") +
  labs(x = NULL, y = "Genes") +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_blank(),
        axis.ticks.x = element_blank(), axis.text.x = element_text(size = 11),
        legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 3))
plot_alluvial_k


# broad-scale share vs stability across K ----
m200 <- split(names(cl$k200), cl$k200)
best_match <- data.frame(module = as.integer(names(m200)), size = lengths(m200))
for (kk in c("k100", "k1000")) {
  mk <- split(names(cl[[kk]]), cl[[kk]])
  best_match[[kk]] <- sapply(m200, function(a)
    max(sapply(mk, function(b) length(intersect(a, b)) / length(union(a, b)))))
}
bs <- mod_cat %>% filter(model == "spacedecorr") %>% dplyr::select(module, broad_share)
st <- best_match %>% left_join(bs, by = "module") %>% filter(size >= 3)

lab <- c("8" = "MHC I", "9" = "Mitotic", "10" = "S-phase", "14" = "Keratin", "6" = "Hypoxia-basal")
plot_share_k <- st %>%
  pivot_longer(c(k100, k1000), names_to = "K", values_to = "jaccard") %>%
  mutate(K = factor(K, c("k100", "k1000"), c("k=100", "k=1000")),
         label = ifelse(as.character(module) %in% names(lab), lab[as.character(module)], NA)) %>%
  ggplot(aes(broad_share, jaccard, col = K)) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.5) +
  geom_point(aes(size = size), alpha = 0.8) +
  ggrepel::geom_text_repel(aes(label = label), size = 2.6, show.legend = FALSE, na.rm = TRUE) +
  scale_color_manual(values = c("#FF7F00", "#00A9CE")) +
  scale_size_continuous(range = c(1, 4), name = "Module size") +
  labs(x = "Broad-scale share (k=200)", y = "Jaccard index with k=200 module", col = "") +
  theme_bw() +
  theme(legend.position = "bottom",
        text = element_text(size = textsize),
        panel.grid.minor = element_blank())
plot_share_k


# combine ----
plot_grid(legend = get_plot_component(plot_density_indep_diffk, 'guide-box-top', return_all = TRUE),
          plot_grid(plot_density_indep_diffk + theme(legend.position = "none"),
                    plot_density_dep_diffk+ theme(legend.position = "none"),
                    plot_modules_coherence_diffk,
                    labels = c("A","B","C"), nrow =1, align = "h", axis = "tb", rel_widths = c(1, 1, 1.4)),
          plot_grid(plot_alluvial_k, plot_share_k, labels = c("D","E"), nrow = 1, rel_widths = c(1.3, 1)),
          ncol = 1, rel_heights =  c(0.06, 0.42, 0.52))
ggsave(paste0(dir.out,"Images/Densities_diffk.jpeg"), height = 260, width = 350, units = "mm")
