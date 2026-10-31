#-------------------------#
#- Liver merfish dataset -#
#-------------------------#
library(tidyverse)
library(ggrepel)
library(cowplot)

dir.data <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real Data/data/Liver"
dir.out  <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real Data/output/Liver"
dir.proj  <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real Data/scripts"

tissues    <- c("AM048", "AM061", "AM042", "AM031", "AM066", "AM062", "AM072")
conditions <- c("healthy", "healthy", "healthy", "diseased", "diseased", "diseased", "diseased")
names(conditions) <- tissues
hep_types <- c("Central Hep", "Portal Hep", "Fibrotic Hep 1", "Fibrotic Hep 2")

textsize <- 12
zone_cols <- c("Pericentral" = "#B07A2E", "Periportal" = "#3F6FA8",
               "Fibrotic 1" = "#723A5A", "Fibrotic 2" = "#5E4F85")

prog_cols <- c("Pericentral zonation" = "#B07A2E", "Periportal zonation" = "#3F6FA8",
               "Injury program" = "#1B7837", "Other" = "grey70")

# Read data ----
meta_data  <- read.csv(file.path(dir.data, "cell_properties_healthy_diseased_merfish.csv"))
raw        <- as.data.frame(read.csv(file.path(dir.data, "X_raw_healthy_diseased_merfish.csv"), header = FALSE))
gene_names <- read.csv(file.path(dir.data, "gene_names_healthy_diseased_merfish.csv"), header = FALSE)
colnames(raw) <- gene_names$V1
rownames(raw) <- meta_data$index

mean_expr <- sapply(tissues, function(s) {
  id <- meta_data$sample_id == s & meta_data$Cell_Type_final %in% hep_types
  colMeans(raw[id, ])
})
genes_keep <- rownames(mean_expr)[apply(mean_expr, 1, min) >= 0.1]
length(genes_keep) # 197

# Gene programs -------
# Watson et al. 2025, Yakubovsky et al. 2026
programs <- list(
  "Pericentral zonation" = c("CYP2E1", "CYP1A2", "CYP3A4", "GLUL"),
  "Periportal zonation"  = c("SDS", "ASS1", "CYP2A6", "PCK1"),
  "Injury program"       = c("SLC7A2", "CPS1", "ATP9A", "DTX1", "POR", "EXT1")
)
programs <- lapply(programs, function(g) intersect(g, genes_keep))   # genes passing the filter
programs
N <- rowSums(raw)


# Splines maps zonation ---------------------------------------------------
tissue <- "AM042"

load(file.path(dir.out, tissue, paste0("Data_", tissue, "_spacedecorr_k200_hepfibcov.Rdata")))
spl_tissue <- cbind(data_spacedecorr_hepcov$metadata, data_spacedecorr_hepcov$splines)   

# (A) cell types ----
plot_celltypes <- meta_data %>%
  filter(Cell_Type_final %in% hep_types, sample_id == tissue) %>%
  mutate(Cell_Type_final = factor(Cell_Type_final, levels = c("Central Hep","Portal Hep",
                                                              "Fibrotic Hep 1","Fibrotic Hep 2"),
                                  labels = c("Pericentral","Periportal","Fibrotic 1","Fibrotic 2"))) %>%
  ggplot(aes(x = x, y = y, col = Cell_Type_final)) +
  geom_point(size = 0.05) +
  coord_equal() +
  scale_color_manual(values = zone_cols) +
  labs(col = "") +
  theme_bw() +
  guides(col = guide_legend(nrow = 2, override.aes = list(size = 3))) +
  theme(legend.position = "top",
        text = element_text(size = textsize),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())

# (B) fitted spatial trends ----
plot_splines <- spl_tissue %>%
  pivot_longer(c(CYP2E1, SDS), names_to = "gene", values_to = "spline") %>%
  mutate(gene = factor(gene, c("CYP2E1", "SDS"),
                       c("italic(CYP2E1)~(pericentral)", "italic(SDS)~(periportal)"))) %>%
  ggplot(aes(sdimx, sdimy, col = spline)) +
  geom_point(size = 0.05) +
  coord_equal() +
  scale_color_scico(palette = "romaO", midpoint = 0, direction = -1, limits = c(-2, 2), oob = scales::squish) +
  facet_wrap(gene ~ .,labeller = label_parsed) +      
  labs(col = "Spatial trend") +
  theme_bw() +
  guides(col = guide_colorbar(barwidth = unit(4, "cm"), title.vjust = 0.8)) +
  theme(legend.position = "top",
        text = element_text(size = textsize),
        strip.text = element_text(size = textsize, face = "italic"),
        strip.background = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank(),
        panel.grid = element_blank())


# spatial R2 per gene vs mean expression (healthy) ------------------------
spatial_r2 <- read.csv(file = file.path(dir.out, "spatial_r2_hepfib.csv")) %>%
  filter(gene %in% genes_keep)

r2_summary <- spatial_r2 %>%
  group_by(condition, gene) %>%
  summarise(r2 = mean(r2), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = r2) %>%
  left_join(data.frame(gene = rownames(mean_expr),
                       mean_healthy = rowMeans(mean_expr[, tissues[conditions == "healthy"]])),
            by = "gene") %>%
  arrange(desc(healthy)) %>%
  mutate(rank_healthy = row_number(),
         marker = case_when(gene %in% programs[[1]] ~ names(programs)[1],
                            gene %in% programs[[2]] ~ names(programs)[2],
                            gene %in% programs[[3]] ~ names(programs)[3],
                            TRUE ~ "Other"))

plot_r2 <- r2_summary %>%
  ggplot(aes(mean_healthy, healthy, col = marker)) +
  geom_point(data = ~ filter(.x, marker == "Other"), size = 0.6) +
  geom_point(data = ~ filter(.x, marker != "Other"), size = 2) +
  geom_text_repel(data = ~ filter(.x, marker != "Other"), aes(label = gene),
                  size = 3, fontface = "italic", max.overlaps = Inf,
                  min.segment.length = 0, box.padding = 0.4, show.legend = FALSE) +
  scale_x_log10() +
  scale_color_manual(values = prog_cols, breaks = names(prog_cols)) +
  labs(x = "Mean counts per hepatocyte (healthy)", y = expression("Variance explained by spatial trend ("*R^2*")"), col = NULL) +
  theme_bw() +
  guides(col = guide_legend(nrow = 2, override.aes = list(size = 2))) +
  theme(legend.position = "top",
        text = element_text(size = textsize),
        panel.grid.minor = element_blank());plot_r2



# Within condition network ------------------------------------------------
fdr_cut <- 0.05
load(file = file.path(dir.out, "pair_cor_liver_hepfib_modttest.Rdata"))

# long format with naive Fisher P-values from r and the number of cells
pc_long <- pair_cor %>% 
  filter(gene1 %in% genes_keep, gene2 %in% genes_keep) %>%
  pivot_longer(vst_r:spacedecorr_pmod, names_to = c("method","metric"), values_to = "value", names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(sample_id, cell_type, method) %>%
  mutate(fdr = p.adjust(p, "BH"),
         fdr_mod = p.adjust(pmod, "BH")) %>%
  ungroup()


# replication within each condition: significant in every donor, same sign
replicated <- pc_long %>%
  filter(cell_type == "Zonated Hep") %>%
  group_by(condition, cell_type, method, gene1, gene2, type) %>%
  summarise(n_donors  = n(),
            mean_r    = mean(r),
            same_sign = all(r > 0) | all(r < 0),
            rep       = all(fdr < 0.05) & same_sign,
            rep_mod   = all(fdr_mod < 0.05) & same_sign,
            .groups = "drop")


# replicated pairs per condition, cell type and method
rep_counts <- replicated %>%
  group_by(condition, cell_type, method) %>%
  summarise(n_pairs = n(), n_rep = sum(rep), n_repmod = sum(rep_mod), .groups = "drop") %>%
  arrange(cell_type, condition, method)
rep_counts

# by pair type, zonated hepatocytes, healthy
replicated %>%
  filter(cell_type == "Zonated Hep", condition == "healthy", rep) %>%
  dplyr::count(method, type) %>%
  pivot_wider(names_from = method, values_from = n, values_fill = 0)

## Network of zonated hep (healthy) ----------------------------------------------------
library(igraph)
library(ggraph)

m           <- "spacedecorr"    
min_module  <- 3              

# replicated pairs in healthy zonated hepatocytes
edges_zh <- replicated %>%
  filter(cell_type == "Zonated Hep", condition == "healthy", method == m, rep) %>%
  transmute(from = gene1, to = gene2, mean_r, weight = abs(mean_r))
nrow(edges_zh)

g_zh <- graph_from_data_frame(edges_zh, directed = FALSE)

# modules (Louvain on |r|)
set.seed(1)
mod_zh <- cluster_louvain(g_zh, weights = E(g_zh)$weight)
modules_zh <- data.frame(name = V(g_zh)$name, module = as.integer(membership(mod_zh)),
                         degree = igraph::degree(g_zh)) %>%
  group_by(module) %>% mutate(n_genes = n()) %>% ungroup() %>%
  filter(n_genes >= min_module) %>%
  mutate(module = paste("Module", as.integer(factor(module))),    
         program = case_when(name %in% programs[[1]] ~ names(programs)[1],
                             name %in% programs[[2]] ~ names(programs)[2],
                             name %in% programs[[3]] ~ names(programs)[3],
                             TRUE ~ "Other"))

# genes per module (for the text)
modules_zh %>% arrange(module, desc(degree)) %>%
  group_by(module) %>%
  summarise(n_genes = n(), genes = paste(name, collapse = ", "), .groups = "drop") %>%
  print(width = Inf)

# keep only edges within modules
edges_within <- edges_zh %>%
  mutate(m1 = modules_zh$module[match(from, modules_zh$name)],
         m2 = modules_zh$module[match(to,   modules_zh$name)]) %>%
  filter(!is.na(m1), m1 == m2)

# sign of edges within each module
edges_within %>% group_by(module = m1) %>%
  summarise(n_edges = n(), frac_positive = mean(mean_r > 0), mean_abs_r = mean(abs(mean_r)), .groups = "drop")

g_within <- graph_from_data_frame(edges_within %>% dplyr::select(from, to, mean_r),
                                  vertices = modules_zh %>% dplyr::select(name, module, program),
                                  directed = FALSE)

# module labels with size, for the facet titles
lab_mod <- modules_zh %>% dplyr::count(module) %>%
  mutate(label = paste0(module, " (", n, " genes)"))
V(g_within)$module_lab <- lab_mod$label[match(V(g_within)$module, lab_mod$module)]

set.seed(1)
plot_zonated_modules <- ggraph(g_within, layout = "stress") +
  geom_edge_link(aes(edge_colour = mean_r > 0, edge_width = abs(mean_r)), alpha = 0.6) +
  geom_node_point(aes(fill = program), shape = 21, size = 4) +
  geom_node_text(aes(label = name), size = 3.5, fontface = "italic", repel = TRUE,
                 max.overlaps = Inf) +
  facet_nodes(~ module_lab, scales = "free") +
  scale_edge_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "#E7298A"),
                           labels = c("TRUE" = "Positive", "FALSE" = "Negative"), name = NULL) +
  scale_edge_width(range = c(0.3, 2), guide = "none") +
  scale_fill_manual(values = prog_cols, breaks = names(prog_cols), name = NULL) +
  guides(fill = guide_legend(override.aes = list(size = 4))) +
  theme_graph() +
  theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        strip.text = element_text(size = 13),
        panel.border = element_rect(fill = NA, colour = "grey70"))
plot_zonated_modules
ggsave(file.path(dir.out,"Images/Modules_spacedecorr_liver.jpeg"), height = 220, width = 300, units = "mm")

# Plot single module
mod_pick <- "Module 5"

nodes_one <- modules_zh %>% filter(module == mod_pick) %>% dplyr::select(name, program)
edges_one <- edges_within %>% filter(m1 == mod_pick) %>% dplyr::select(from, to, mean_r)

g_one <- graph_from_data_frame(edges_one, vertices = nodes_one, directed = FALSE)

set.seed(1)
plot_one_module <- ggraph(g_one, layout = "stress") +    
  geom_edge_link(aes(edge_colour = mean_r > 0, edge_width = abs(mean_r)), alpha = 0.6) +
  geom_node_point(aes(fill = program), shape = 21, size = 3.5) +
  geom_node_text(aes(label = name), size = 3, fontface = "italic", repel = TRUE) +
  scale_edge_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "#E7298A"),
                           labels = c("TRUE" = "Positive", "FALSE" = "Negative"), name = NULL) +
  scale_edge_width(range = c(0.3, 2), guide = "none") +
  scale_fill_manual(values = prog_cols, name = NULL) +
  theme_graph() +       
  guides(fill = guide_legend(nrow = 2), edge_colour = guide_legend(nrow = 2)) +
  theme(legend.position = "top",
        text = element_text(size = textsize),
        plot.margin = margin(5.5, 5.5, 5.5, 5.5),
        panel.border = element_rect(fill = NA, colour = "grey20"))
plot_one_module


## Module stability --------------------------------------------------------
library(igraph)
library(ggraph)
library(cowplot)

rep_h     <- replicated %>% filter(cell_type == "Zonated Hep", condition == "healthy", rep)
e_sd      <- rep_h %>% filter(method == "spacedecorr")
e_vst_fdr <- rep_h %>% filter(method == "vst")
e_vst_top <- e_vst_fdr %>% slice_max(abs(mean_r), n = nrow(e_sd), with_ties = FALSE)
nets <- list("SpaceDecorr" = e_sd, "VST (matched FDR)" = e_vst_fdr, "VST (matched density)" = e_vst_top)

methods_cl <- c("Louvain", "Walktrap")
injury <- programs[["Injury program"]]
plots_full <- list(); plots_inj <- list()

prog_fill <- c("Pericentral zonation" = "#B2182B", "Periportal zonation" = "#2166AC",
               "Injury program" = "#E08214", "Other" = "grey85")

for (nm in names(nets)) {
  g <- graph_from_data_frame(nets[[nm]] %>% transmute(from = gene1, to = gene2,
                                                      weight = abs(mean_r), mean_r = mean_r),
                             directed = FALSE)
  w <- E(g)$weight
  set.seed(1)
  lay <- layout_with_fr(g, weights = w)
  
  # module memberships
  set.seed(1)
  memb <- list(Louvain  = membership(cluster_louvain(g, weights = w)),
               Walktrap = membership(cluster_walktrap(g, weights = w)))
  
  V(g)$program <- case_when(V(g)$name %in% programs[[1]] ~ names(programs)[1],
                            V(g)$name %in% programs[[2]] ~ names(programs)[2],
                            V(g)$name %in% programs[[3]] ~ names(programs)[3],
                            TRUE ~ "Other")
  V(g)$injury <- V(g)$name %in% injury
  
  for (cl in methods_cl) {
    mm       <- memb[[cl]][V(g)$name]
    inj_mod  <- unique(mm[V(g)$injury])
    inj_size <- sum(mm %in% inj_mod)
    sub_title <- paste0(length(unique(mm)), " modules; injury genes in ", length(inj_mod),
                        " module(s), ", inj_size, " genes")
    
    # (a) full network coloured by module
    V(g)$module <- factor(mm)
    plots_full[[paste(nm, cl)]] <- ggraph(g, layout = "manual", x = lay[, 1], y = lay[, 2]) +
      geom_edge_link(colour = "grey80", alpha = ifelse(nm == "VST (matched FDR)", 0.05, 0.4), width = 0.2) +
      geom_node_point(aes(fill = module, size = injury), shape = 21,
                      colour = ifelse(V(g)$injury, "black", "grey40"),
                      stroke = ifelse(V(g)$injury, 0.8, 0.2)) +
      geom_node_text(aes(label = ifelse(injury, name, "")), size = 3.2, fontface = "bold.italic",
                     repel = TRUE, max.overlaps = Inf, min.segment.length = 0) +
      scale_size_manual(values = c("FALSE" = 2, "TRUE" = 4), guide = "none") +
      labs(subtitle = paste0(length(unique(mm)), " modules, injury genes in ", length(inj_mod),
                             ifelse(length(inj_mod) == 1, " module", " modules"),
                             " (", inj_size, " genes)")) +
      theme_graph() +
      theme(legend.position = "none",
            plot.subtitle = element_text(size = 11, hjust = 0.5),
            plot.margin = margin(5, 5, 5, 5),
            panel.border = element_rect(fill = NA, colour = "grey70"))
    
    # (b) only the module(s) containing injury genes
    sub <- induced_subgraph(g, V(g)$name[mm %in% inj_mod])
    V(sub)$module <- factor(mm[V(sub)$name])
    set.seed(1)
    plots_inj[[paste(nm, cl)]] <- ggraph(sub, layout = "fr") +
      geom_edge_link(aes(edge_colour = mean_r > 0, edge_width = abs(mean_r)), alpha = 0.4) +
      geom_node_point(aes(fill = program, colour = module), shape = 21, size = 2.5, stroke = 0.8) +
      geom_node_text(aes(label = name), size = ifelse(vcount(sub) > 60, 1.5, 2.2), repel = TRUE) +
      scale_edge_colour_manual(values = c("TRUE" = "grey40", "FALSE" = "#E7298A"), guide = "none") +
      scale_edge_width(range = c(0.1, 1), guide = "none") +
      scale_fill_manual(values = prog_cols, name = NULL) +
      guides(colour = "none") +
      labs(title = paste0(nm, " - ", cl), subtitle = paste0(vcount(sub), " genes in injury module(s)")) +
      theme_graph() +
      theme(legend.position = "bottom", plot.title = element_text(size = 10), plot.subtitle = element_text(size = 8))
  }
}

# column headers (networks) and row labels (clustering method)
header <- plot_grid(plotlist = lapply(names(nets), function(nm)
  ggdraw() + draw_label(nm, fontface = "bold", size = 13)), nrow = 1)
rowlab <- plot_grid(plotlist = lapply(methods_cl, function(cl)
  ggdraw() + draw_label(cl, fontface = "bold", size = 13, angle = 90)), ncol = 1)

grid_full <- plot_grid(plotlist = plots_full, nrow = length(methods_cl), byrow = F)

plot_grid(NULL, header,
          rowlab, grid_full,
          ncol = 2, rel_widths = c(0.04, 1), rel_heights = c(0.06, 1))
ggsave(file.path(dir.out, "Images/Modules_sensitivity_liver.jpeg"), height = 200, width = 280, units = "mm")


# Per sample per cell type pairwise correlations --------------------------
load(file = file.path(dir.out, "res_liver_welch_by_celltype.Rdata"))
welch_res <- welch_res %>% filter(gene1 %in% genes_keep, gene2 %in% genes_keep,
                                  cell_type == "Zonated Hep") 

# differential pairs (FDR < 0.05) by cell type, method, direction; and how many reach the minimum permutation p
welch_res %>% 
  filter(fdr < 0.10) %>%
  group_by(cell_type, method) %>%
  summarise(n = n(), 
            separated = sum(as.numeric(separated)),
            n_perm_min = sum(p_perm <= min_p_perm + 1e-12), .groups = "drop") %>%
  arrange(cell_type, method)


# Healthy vs fibrotic correlation for every pair, in zonated hepat --------
plot_scatter <- welch_res %>%
  filter(cell_type == "Zonated Hep", method == "spacedecorr") %>%
  mutate(differential = fdr < 0.10) %>%
  arrange(differential) %>%
  ggplot(aes(r_healthy, r_diseased, col = differential, size = differential)) +
  geom_abline(lty = 2, col = "grey50") +
  geom_hline(yintercept = 0, col = "grey85") + geom_vline(xintercept = 0, col = "grey85") +
  geom_point(alpha = 0.7) +
  coord_equal(xlim = c(-0.15, 0.3), ylim = c(-0.15, 0.3)) +   # zooms without dropping points
  scale_color_manual(values = c("FALSE" = "grey75", "TRUE" = "#B2182B"),
                     labels = c("Not differential", "Differential (FDR < 0.10)")) +
  scale_size_manual(values = c("FALSE" = 0.5, "TRUE" = 1.5), guide = "none") +
  labs(x = "Correlation, healthy", y = "Correlation, fibrotic", col = NULL) +
  theme_bw() +
  guides(col = guide_legend(nrow = 2, override.aes = list(size = 2))) +
  theme(legend.position = "top",
        text = element_text(size = textsize),
        panel.grid.minor = element_blank())
plot_scatter


# Save plot ---------------------------------------------------------------

plot_grid(plot_grid(plot_celltypes, plot_splines, nrow = 1, rel_widths = c(1, 2),
                    labels = c("A", "B"), align = "h", axis = "tb"),
          plot_grid(plot_r2, plot_one_module, plot_scatter, nrow = 1,
                    labels = c("C", "D", "E"), rel_widths = c(1.1, 1, 0.9)),
          ncol = 1, rel_heights = c(1, 1))
ggsave(file.path(dir.out, "Images/PlotsLiver.pdf"), height = 220, width = 330, units = "mm", dpi = 400)

# Sensitivity K -----------------------------------------------------------
load(file = file.path(dir.out, "res_liver_K_sensitivity.Rdata"))

# keep the 101 genes, redo BH within K and sample, and redo the pair type
sens_r2 <- sens_r2 %>% filter(gene %in% genes_keep)
sens_pc <- sens_pc %>%
  filter(gene1 %in% genes_keep, gene2 %in% genes_keep) %>%
  group_by(K, sample_id) %>%
  mutate(fdr = p.adjust(p, "BH")) %>%
  ungroup() %>%
  mutate(type = case_when(
    gene1 %in% programs[[1]] & gene2 %in% programs[[1]] ~ names(programs)[1],
    gene1 %in% programs[[2]] & gene2 %in% programs[[2]] ~ names(programs)[2],
    gene1 %in% programs[[3]] & gene2 %in% programs[[3]] ~ names(programs)[3],
    (gene1 %in% programs[[1]] & gene2 %in% programs[[2]]) |
      (gene1 %in% programs[[2]] & gene2 %in% programs[[1]]) ~ "Opposite zones",
    TRUE ~ "Other"))

# (1) spatial R2 of program genes by K: the spline should capture more zonation as K grows
sens_r2 %>%
  mutate(program = case_when(gene %in% programs[[1]] ~ names(programs)[1],
                             gene %in% programs[[2]] ~ names(programs)[2],
                             gene %in% programs[[3]] ~ names(programs)[3],
                             TRUE ~ "Other")) %>%
  group_by(K, program) %>%
  summarise(mean_r2 = mean(r2), .groups = "drop") %>%
  pivot_wider(names_from = K, values_from = mean_r2, names_prefix = "K")

# (2) mean residual correlation by pair type and K (mean over healthy donors)
sens_type <- sens_pc %>%
  group_by(K, sample_id, type) %>%
  summarise(mean_r = mean(r, na.rm = TRUE), .groups = "drop") %>%
  group_by(K, type) %>%
  summarise(mean_r = mean(mean_r), .groups = "drop")
sens_type %>% pivot_wider(names_from = K, values_from = mean_r, names_prefix = "K")

# (3) pairs replicated in all healthy donors (FDR < 0.05, same sign), by K and pair type
sens_rep <- sens_pc %>%
  group_by(K, gene1, gene2, type) %>%
  summarise(rep = all(fdr < fdr_cut) & (all(r > 0) | all(r < 0)), mean_r = mean(r), .groups = "drop")
sens_rep %>% filter(rep) %>% dplyr::count(K, type) %>%
  pivot_wider(names_from = K, values_from = n, values_fill = 0, names_prefix = "K")

# overlap of replicated pairs with K = 200
rep200 <- sens_rep %>% filter(K == 200, rep) %>% mutate(pair = paste(gene1, gene2)) %>% pull(pair)
sens_rep %>% filter(rep) %>% mutate(pair = paste(gene1, gene2)) %>%
  group_by(K) %>% summarise(n_rep = n(), shared_with_K200 = sum(pair %in% rep200), .groups = "drop")
# check: n_rep at K = 200 should be 491, as in the main text

# Injury modules across k
library(igraph)
injury <- programs[["Injury program"]]

inj_k <- list()
for (K in c(100, 200, 1000)) {
  e <- sens_rep %>% filter(K == !!K, rep) %>% transmute(from = gene1, to = gene2, weight = abs(mean_r))
  g <- graph_from_data_frame(e, directed = FALSE)
  set.seed(1)
  mm <- membership(cluster_louvain(g, weights = E(g)$weight))
  inj_k[[as.character(K)]] <- names(mm)[mm %in% unique(mm[intersect(injury, names(mm))])]
  print(c(K = K, n_modules = length(unique(mm)), injury_modules = length(unique(mm[intersect(injury, names(mm))])),
          size = length(inj_k[[as.character(K)]])))
}

# Jaccard with the K = 200 injury module
sapply(inj_k, function(x) length(intersect(x, inj_k[["200"]])) / length(union(x, inj_k[["200"]])))

type_lev  <- c(names(programs), "Opposite zones", "Other")
type_cols <- c("Pericentral zonation" = "#B07A2E", "Periportal zonation" = "#3F6FA8",
               "Injury program" = "#1B7837", "Opposite zones" = "#4D4D4D", "Other" = "grey60")

# (A) spatial R2 by program
plot_r2_k <- sens_r2 %>%
  mutate(program = case_when(gene %in% programs[[1]] ~ names(programs)[1],
                             gene %in% programs[[2]] ~ names(programs)[2],
                             gene %in% programs[[3]] ~ names(programs)[3],
                             TRUE ~ "Other")) %>%
  group_by(K, program) %>%
  summarise(mean_r2 = mean(r2), .groups = "drop") %>%
  mutate(program = factor(program, type_lev)) %>%
  ggplot(aes(factor(K), mean_r2, group = program, col = program)) +
  geom_line() + geom_point() +
  scale_color_manual(values = type_cols, drop = FALSE) +
  labs(x = "Spline basis dimension (K)", y = expression("Mean spatial "*R^2), col = NULL) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank())

# (B) mean residual correlation by pair type
plot_cor_k <- sens_type %>%
  mutate(type = factor(type, type_lev)) %>%
  ggplot(aes(factor(K), mean_r, group = type, col = type)) +
  geom_hline(yintercept = 0, lty = 2, col = "grey60") +
  geom_line() + geom_point() +
  scale_color_manual(values = type_cols, drop = FALSE) +
  labs(x = "Spline basis dimension (K)", y = "Mean residual correlation", col = NULL) +
  theme_bw() +
  guides(col = guide_legend(nrow = 1)) +
  theme(legend.position = "top", panel.grid.minor = element_blank())


# (C) edf by program --------------------------------------------------------
# keep the 101 genes, as for sens_r2
sens_edf <- sens_edf %>% filter(gene %in% genes_keep)
n_distinct(sens_edf$gene)                                  # should be 101

plot_edf_k <- sens_edf %>%
  group_by(K, gene, program) %>%
  summarise(edf = mean(edf), .groups = "drop") %>%         # mean over healthy samples
  group_by(K, program) %>%
  summarise(median_edf = median(edf), .groups = "drop") %>%
  mutate(program = factor(program, type_lev)) %>%
  ggplot(aes(factor(K), median_edf, group = program, col = program)) +
  geom_point(data = data.frame(K = c(100, 200, 1000), median_edf = c(99, 199, 999)),
             aes(factor(K), median_edf), inherit.aes = FALSE,
             shape = 95, size = 8, col = "grey60") +       # maximum edf, K - 1
  geom_line() + geom_point() +
  scale_y_log10() +
  scale_color_manual(values = type_cols, drop = FALSE) +
  labs(x = "Spline basis dimension (K)", y = "Median edf", col = NULL) +
  theme_bw() +
  theme(legend.position = "none", panel.grid.minor = element_blank())

plot_grid(get_plot_component(plot_cor_k, "guide-box-top", return_all = TRUE),
          plot_grid(plot_edf_k, plot_r2_k, plot_cor_k + theme(legend.position = "none"),
                    labels = c("A", "B", "C"), nrow = 1, align = "h", axis = "tb"),
          ncol = 1, rel_heights = c(0.08, 1))
ggsave(file.path(dir.out, "Images/Liver_K_sensitivity.jpeg"), height = 100, width = 330, units = "mm")



# EDF ---------------------------------------------------------------------
edf_gene <- sens_edf %>%
  group_by(K, gene, program) %>%
  summarise(edf = mean(edf), .groups = "drop")

edf_gene %>%
  group_by(K) %>%
  summarise(median = median(edf), q1 = quantile(edf, 0.25), q3 = quantile(edf, 0.75), max = max(edf),
            pct_near_max   = 100 * mean(edf > 0.9 * (K - 1)),
            pct_within_199 = 100 * mean(edf < 199), .groups = "drop")

# by program
edf_gene %>%
  group_by(K, program) %>%
  summarise(median_edf = median(edf), .groups = "drop") %>%
  pivot_wider(names_from = K, values_from = median_edf, names_prefix = "K")



# Modified test -----------------------------------------------------------
# ESS and significant pairs per method
mycol <- c("SpaceDecorr" = "#8C510A", "VST-normalized" = "#7F7F7F", "log-normalized" = "#000000")
plot_df <- pc_long %>% filter(cell_type == "Pooled",
                              gene1 %in% genes_keep, gene2 %in% genes_keep) %>%
  mutate(method = factor(method, levels = c("spacedecorr","vst"),
                         labels = c("SpaceDecorr","VST-normalized")),
         ess_prop = ess / n_cells,
         lp_naive    = -log10(pmax(fdr, 1e-300)),      # avoid Inf
         lp_modified = -log10(pmax(fdr_mod, 1e-300)))

# (A) effective sample size as a proportion of cells
p_ess_liver <- plot_df %>%
  ggplot(aes(method, ess_prop, col = method)) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_boxplot(outlier.size = 0.3) +
  scale_color_manual(values = mycol) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Effective sample size (% of cells)") +
  theme_bw() +
  theme(legend.position = "none")

# (B) standard vs modified t-test, -log10 FDR
p_scatter_liver <- plot_df %>%
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

plot_grid(p_ess_liver, p_scatter_liver, labels = c("A", "B"), nrow = 1, rel_widths = c(1.5, 2.5), align = "h", axis = "tb")
ggsave(file.path(dir.out,"Images/Modified_ESS_liver.jpeg"), height = 120, width = 300, units = "mm")


plot_df %>% filter( lp_modified < 0.0000001, method == "SpaceDecorr", fdr < 0.05) %>%
  count(gene1) %>% arrange(-n)


# spMOCA ------------------------------------------------------------------
# spMOCA could not be run for AM066 (177,538 hepatocytes).
library(igraph)
library(ggraph)

# t-test for a correlation, as in the lung analysis
test1cor <- function(r, n) {
  fisher <- r / sqrt(1 - r^2) * sqrt(n - 2)
  2 * pt(-abs(fisher), n - 2)
}

injury     <- programs[["Injury program"]]
min_module <- 3
type_lev   <- c(names(programs), "Opposite zones", "Other")
type_cols  <- c("Pericentral zonation" = "#B07A2E", "Periportal zonation" = "#3F6FA8",
                "Injury program" = "#1B7837", "Opposite zones" = "#4D4D4D", "Other" = "grey75")
prog_cols  <- c("Pericentral zonation" = "#B07A2E", "Periportal zonation" = "#3F6FA8",
                "Injury program" = "#1B7837", "Other" = "grey85")

# 1. spMOCA correlations per sample, healthy donors ----------------------------
moca_long <- data.frame()
for (tissue in tissues[conditions == "healthy"]) {
  print(tissue)
  load(file.path(dir.out, tissue, paste0("Data_", tissue, "_spmoca_hep.Rdata")))
  C <- as.matrix(cor_spmoca_hep)
  
  genes <- sort(intersect(genes_keep, colnames(C)))
  C  <- C[genes, genes]
  ij <- which(upper.tri(C), arr.ind = TRUE)
  
  # number of zonated hepatocytes spMOCA was fitted on
  n <- sum(meta_data$sample_id == tissue & meta_data$Cell_Type_final %in% c("Central Hep", "Portal Hep"))
  
  moca_long <- bind_rows(moca_long,
                         data.frame(sample_id = tissue, n = n,
                                    gene1 = genes[ij[, 1]], gene2 = genes[ij[, 2]], r = C[ij]) %>%
                           mutate(p   = test1cor(r, n),
                                  fdr = p.adjust(p, "BH")))
}

# same pair orientation and pair types as the SpaceDecorr results
sd_pairs <- replicated %>%
  filter(cell_type == "Zonated Hep", condition == "healthy", method == "spacedecorr") %>%
  dplyr::select(gene1, gene2, type, mean_r_sd = mean_r)

moca_long <- moca_long %>%
  left_join(sd_pairs %>% dplyr::select(gene1, gene2), by = c("gene1", "gene2"), keep = TRUE) %>%
  mutate(flip = is.na(gene1.y)) %>%
  transmute(sample_id, n,
            gene1 = ifelse(flip, gene2.x, gene1.x),
            gene2 = ifelse(flip, gene1.x, gene2.x),
            r, p, fdr)

# 2. Replication in healthy donors: FDR < 0.05 in every donor, same sign --------
moca_rep <- moca_long %>%
  group_by(gene1, gene2) %>%
  summarise(mean_r    = mean(r),
            same_sign = all(r > 0) | all(r < 0),
            rep       = all(fdr < 0.05) & same_sign, .groups = "drop") %>%
  left_join(sd_pairs, by = c("gene1", "gene2")) %>%
  left_join(replicated %>%
              filter(cell_type == "Zonated Hep", condition == "healthy", method == "vst") %>%
              dplyr::select(gene1, gene2, mean_r_vst = mean_r, rep_vst = rep),
            by = c("gene1", "gene2")) %>%
  mutate(type = ifelse(type == "Other pairs", "Other", type))

# replicated pairs
c(pairs = nrow(moca_rep), spmoca = sum(moca_rep$rep))

# mean correlation by pair type 
moca_rep %>%
  group_by(type) %>%
  summarise(n = n(), spmoca = mean(mean_r), spacedecorr = mean(mean_r_sd), vst = mean(mean_r_vst),
            .groups = "drop")

# agreement across all pairs 
c(spmoca_vst = cor(moca_rep$mean_r, moca_rep$mean_r_vst, method = "spearman"),
  spmoca_sd  = cor(moca_rep$mean_r, moca_rep$mean_r_sd,  method = "spearman"),
  sd_vst     = cor(moca_rep$mean_r_sd, moca_rep$mean_r_vst, method = "spearman"))

# 3. Network overlap with VST at matched FDR -------------------
e_moca <- moca_rep %>% filter(rep)
e_vst  <- moca_rep %>% filter(rep_vst)
shared <- nrow(semi_join(e_moca, e_vst, by = c("gene1", "gene2")))
c(n_spmoca = nrow(e_moca), n_vst = nrow(e_vst),
  jaccard = shared / (nrow(e_moca) + nrow(e_vst) - shared))

# 4. Louvain modules and the injury genes (text: two modules, CPS1 with zonation genes) ---
g <- graph_from_data_frame(e_moca %>% transmute(from = gene1, to = gene2, mean_r, weight = abs(mean_r)),
                           directed = FALSE)
set.seed(1)
louv <- membership(cluster_louvain(g, weights = E(g)$weight))
louv[intersect(injury, names(louv))]                        # module of each injury gene
table(louv)                                                 # module sizes

mod_df <- data.frame(name = V(g)$name, module = as.integer(louv[V(g)$name])) %>%
  group_by(module) %>% mutate(n_genes = n()) %>% ungroup() %>%
  filter(n_genes >= min_module) %>%
  mutate(module = paste0("Module ", as.integer(factor(module)), " (", n_genes, " genes)"),
         program = case_when(name %in% programs[[1]] ~ names(programs)[1],
                             name %in% programs[[2]] ~ names(programs)[2],
                             name %in% programs[[3]] ~ names(programs)[3],
                             TRUE ~ "Other"))
e_within <- e_moca %>%
  transmute(from = gene1, to = gene2, mean_r) %>%
  mutate(m1 = mod_df$module[match(from, mod_df$name)],
         m2 = mod_df$module[match(to,   mod_df$name)]) %>%
  filter(!is.na(m1), m1 == m2)
g_within <- graph_from_data_frame(e_within %>% dplyr::select(from, to, mean_r),
                                  vertices = mod_df %>% dplyr::select(name, module, program),
                                  directed = FALSE)


# (A) correlations, spMOCA vs SpaceDecorr and VST ----
plot_moca <- moca_rep %>%
  pivot_longer(c(mean_r_sd, mean_r_vst), names_to = "ref", values_to = "r_ref") %>%
  mutate(ref  = factor(ref, c("mean_r_sd", "mean_r_vst"), c("SpaceDecorr", "VST-normalized")),
         type = factor(type, type_lev)) %>%
  arrange(desc(type)) %>%                      # program pairs drawn on top
  ggplot(aes(r_ref, mean_r, col = type, size = type == "Other", alpha = type == "Other")) +
  geom_abline(lty = 2, col = "grey50") +
  geom_hline(yintercept = 0, col = "grey85") + geom_vline(xintercept = 0, col = "grey85") +
  geom_point() +
  scale_size_manual(values = c("TRUE" = 0.5, "FALSE" = 1.8), guide = "none") +
  scale_alpha_manual(values = c("TRUE" = 0.4, "FALSE" = 0.9), guide = "none") +
  scale_color_manual(values = type_cols) +
  facet_wrap(~ ref) +
  labs(x = "Correlation (mean over healthy donors)",
       y = "spMOCA correlation", col = NULL) +
  theme_bw() +
  guides(col = guide_legend(nrow = 1, override.aes = list(size = 2.5, alpha = 1))) +
  theme(legend.position = "top", strip.background = element_blank(),
        strip.text = element_text(size = 11), panel.grid.minor = element_blank())

# (B) modules for spMOCA, matched FDR and matched density, Louvain and Walktrap ----
n_sd   <- replicated %>% filter(cell_type == "Zonated Hep", condition == "healthy",
                                method == "spacedecorr", rep) %>% nrow()     # 491
e_dens <- moca_rep %>% filter(same_sign) %>% slice_max(abs(mean_r), n = n_sd, with_ties = FALSE)
nets_m <- list("Matched FDR" = e_moca, "Matched density" = e_dens)

plots_b <- list()
for (nm in names(nets_m)) {
  e <- nets_m[[nm]] %>% transmute(from = gene1, to = gene2, mean_r, weight = abs(mean_r))
  g <- graph_from_data_frame(e, directed = FALSE)
  set.seed(1)
  memb <- list(Louvain  = membership(cluster_louvain(g, weights = E(g)$weight)),
               Walktrap = membership(cluster_walktrap(g, weights = E(g)$weight)))
  
  # backbone for drawing: three strongest edges per gene (all edges for the sparse network)
  e_draw <- bind_rows(e %>% mutate(gg = from), e %>% mutate(gg = to)) %>%
    group_by(gg) %>% slice_max(weight, n = 3, with_ties = FALSE) %>% ungroup() %>%
    distinct(from, to, .keep_all = TRUE)
  if (nm == "Matched density") e_draw <- e
  gd <- graph_from_data_frame(e_draw %>% dplyr::select(from, to, mean_r),
                              vertices = data.frame(name = V(g)$name), directed = FALSE)
  V(gd)$program <- case_when(V(gd)$name %in% programs[[1]] ~ names(programs)[1],
                             V(gd)$name %in% programs[[2]] ~ names(programs)[2],
                             V(gd)$name %in% programs[[3]] ~ names(programs)[3],
                             TRUE ~ "Other")
  set.seed(1)
  lay <- graphlayouts::layout_with_stress(gd)  
  
  for (cl in names(memb)) {
    mm <- memb[[cl]][V(gd)$name]
    V(gd)$module <- factor(mm)
    inj_mod <- unique(mm[V(gd)$name %in% injury])
    plots_b[[paste(nm, cl)]] <- ggraph(gd, layout = "manual", x = lay[, 1], y = lay[, 2]) +
      geom_edge_link(colour = "grey70", alpha = 0.4, width = 0.2) +
      geom_node_point(aes(fill = module, size = program != "Other"), shape = 21,
                      colour = "grey30", stroke = 0.3) +
      geom_node_text(aes(label = ifelse(program != "Other", name, ""), colour = program),
                     size = 3.2, fontface = "bold.italic", repel = TRUE,
                     max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) +
      scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 3.5), guide = "none") +
      scale_fill_hue(l = 80, c = 50, guide = "none") +
      scale_colour_manual(values = prog_cols) +
      labs(subtitle = paste0(length(unique(mm)), " modules, injury genes in ", length(inj_mod),
                             ifelse(length(inj_mod) == 1, " module", " modules"))) +
      theme_graph() +
      theme(plot.subtitle = element_text(size = 11, hjust = 0.5),
            plot.margin = margin(5, 5, 5, 5),
            panel.border = element_rect(fill = NA, colour = "grey70"))
  }
}

# grid: columns = networks, rows = clustering
header <- plot_grid(plotlist = lapply(names(nets_m), function(nm)
  ggdraw() + draw_label(nm, fontface = "bold", size = 12)), nrow = 1)
rowlab <- plot_grid(plotlist = lapply(c("Louvain", "Walktrap"), function(cl)
  ggdraw() + draw_label(cl, fontface = "bold", size = 12, angle = 90)), ncol = 1)
grid_b <- plot_grid(plotlist = plots_b, nrow = 2, byrow = FALSE)

# legend for the label colours (programs)
leg_prog <- get_plot_component(
  ggplot(data.frame(p = factor(names(prog_cols)[1:3], names(prog_cols)[1:3])), aes(p, 1, col = p)) +
    geom_point(size = 3) + scale_colour_manual(values = prog_cols, name = "Gene labels") +
    theme_bw() + theme(legend.position = "top"),
  "guide-box-top", return_all = TRUE)

plot_moca_modules <- plot_grid(leg_prog,
                               plot_grid(NULL, header, rowlab, grid_b, ncol = 2,
                                         rel_widths = c(0.04, 1), rel_heights = c(0.06, 1)),
                               ncol = 1, rel_heights = c(0.06, 1))

# combine with (A) ----
plot_grid(plot_moca, plot_moca_modules, ncol = 1, labels = c("A", "B"), rel_heights = c(0.7, 1.3))
ggsave(file.path(dir.out, "Images/spMOCA_liver.jpeg"), height = 300, width = 260, units = "mm")

