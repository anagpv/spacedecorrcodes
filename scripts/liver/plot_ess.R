#------------------------------------------------------------------------------#
#- Supplementary Note: modified t-test and effective sample size (figure)     -#
#------------------------------------------------------------------------------#
library(tidyverse)
library(cowplot)

dir.sim   <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/"
dir.lung  <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real data/output/Lung/"
dir.liver <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real data/output/Liver/"

mycol <- c("SpaceDecorr" = "#8C510A", "VST-normalized" = "#7F7F7F", "log-normalized" = "#000000")
method_labels <- c(spacedecorr = "SpaceDecorr",
                   vst = "VST-normalized", spark = "VST-normalized",
                   lognorm = "log-normalized", logpearson = "log-normalized")
fdr_cut <- 0.05


# 1. Simulations ------------------------------------------------------------------
pairwise_modified_ttest_results <- read_csv(file.path(dir.sim, "pairwise_modified_ttest_results2.csv"))

sim_labels <- c("Additive"   = "Spatial covariance\n(Additive)",
                "MV"         = "Spatial covariance\n(Matrix variate)",
                "FOV_diff"   = "FOV technical factor",
                "subtype_v2" = "Unobserved subtype")
unique(pairwise_modified_ttest_results$sim) 

t1_df <- pairwise_modified_ttest_results %>%
  filter(rho == 0, sim %in% names(sim_labels)) %>%
  mutate(sim    = factor(sim_labels[sim], levels = sim_labels),
         method = factor(method_labels[model], levels = names(mycol)),
         test   = factor(test, levels = c("cor", "modified"),
                         labels = c("Standard t-test", "Modified t-test"))) %>%
  group_by(sim, num_cells, method, test) %>%
  summarise(t1err = mean(pvalue < 0.05, na.rm = TRUE), n_reps = n(), .groups = "drop")

count(t1_df, n_reps)
t1_df %>%
  dplyr::select(-n_reps) %>%
  pivot_wider(names_from = test, values_from = t1err) %>%
  arrange(sim, method, num_cells) %>%
  print(n = Inf)


p_modsim <- t1_df %>%
  ggplot(aes(factor(num_cells), t1err, col = method, lty = test, group = interaction(method, test))) +
  geom_hline(yintercept = 0.05, lty = 2, col = "grey60") +
  geom_line() +
  geom_point(size = 1.2) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_discrete(labels = function(x) scales::comma(as.numeric(x))) +
  scale_color_manual(values = mycol) +
  scale_linetype_manual(values = c("Standard t-test" = "solid", "Modified t-test" = "dashed")) +
  labs(x = "Number of cells", y = "Type I error", col = NULL, lty = NULL) +
  theme_bw() +
  theme(legend.position = "top", strip.background = element_blank(),
        panel.grid.minor = element_blank(), legend.key.width = unit(1.2, "cm"))


# 2. NSCLC (Lung6) ------------------------------------------------------------------
load(file.path(dir.lung, "res_modified_ttest_Lung6.Rdata"))   # res_nsclc
n_cells_lung <- 66267

plot_df <- res_k %>%
  mutate(method = factor(method, levels = c("spacedecorr200","vst","lognorm"),
                         labels = c("SpaceDecorr","VST-normalized","log-normalized")),
         ess_prop    = ess / n_cells_lung,
         lp_naive    = -log10(pmax(fdr_naive, 1e-300)),      # avoid Inf
         lp_modified = -log10(pmax(fdr_modified, 1e-300))) %>%
  filter(!is.na(method))

# ESS as % of cells and significant pairs per method
res_k %>%
  mutate(ess_prop    = ess / n_cells_lung,
         lp_naive    = -log10(pmax(fdr_naive, 1e-300)),      # avoid Inf
         lp_modified = -log10(pmax(fdr_modified, 1e-300))) %>%
  group_by(method) %>%
  summarise(n_pairs = n(),
            ess_median = 100 * median(ess_prop),
            ess_q10 = 100 * quantile(ess_prop, 0.1), ess_q90 = 100 * quantile(ess_prop, 0.9),
            ess_at_min = sum(ess <= 3),
            sig_naive = sum(fdr_naive < fdr_cut), sig_modified = sum(fdr_modified < fdr_cut),
            .groups = "drop")

# (B) effective sample size as a proportion of cells
p_ess <- plot_df %>%
  ggplot(aes(method, ess_prop, col = method)) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_boxplot(outlier.size = 0.3) +
  scale_color_manual(values = mycol) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Effective sample size (% of cells)") +
  theme_bw() +
  theme(legend.position = "none")

# (C) standard vs modified t-test, -log10 FDR
p_scatter <- plot_df %>%
  ggplot(aes(lp_naive, lp_modified, col = method)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_abline(lty = 2, col = "grey40") +
  geom_hline(yintercept = -log10(fdr_cut), lty = 3, col = "grey60") +
  geom_vline(xintercept = -log10(fdr_cut), lty = 3, col = "grey60") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_color_manual(values = mycol) +
  facet_wrap(~ method) +
  labs(x = expression(-log[10]~"FDR, standard t-test"),
       y = expression(-log[10]~"FDR, modified t-test")) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank())


# 3. Liver, zonated hepatocytes -----------------------------------------------------
# per-sample output of liver_modttest_compare.R (array "ours"), copied from the cluster
fdr_cut <- 0.05
load(file = file.path(dir.liver, "pair_cor_liver_hepfib_modttest.Rdata"))

# long format with naive Fisher P-values from r and the number of cells
pc_long <- pair_cor %>% 
  filter(gene1 %in% genes_keep, gene2 %in% genes_keep) %>%
  pivot_longer(vst_r:spacedecorr_pmod, names_to = c("method","metric"), values_to = "value", names_pattern = "(.*)_(.*)") %>%
  pivot_wider(names_from = metric, values_from = value) %>%
  group_by(sample_id, cell_type, method) %>%
  mutate(fdr = p.adjust(p, "BH"),
         fdr_mod = p.adjust(pmod, "BH")) %>%
  ungroup()

plot_df_liver <- pc_long %>% filter(cell_type == "Pooled",
                                    gene1 %in% genes_keep, gene2 %in% genes_keep) %>%
  mutate(method = factor(method, levels = c("spacedecorr","vst"),
                         labels = c("SpaceDecorr","VST-normalized")),
         ess_prop = ess / n_cells,
         lp_naive    = -log10(pmax(fdr, 1e-300)),      # avoid Inf
         lp_modified = -log10(pmax(fdr_mod, 1e-300)))

# numbers for the text, per sample (should be 5,050 pairs per sample and method)
ess_liver <- plot_df_liver %>%
  group_by(sample_id, condition, method) %>%
  summarise(n_cells = dplyr::first(n_cells), n_pairs = n(),
            ess_median = 100 * median(ess_prop),
            ess_at_min = sum(ess <= 3),
            sig_naive = sum(fdr < fdr_cut), sig_modified = sum(fdr_mod < fdr_cut),
            .groups = "drop")
ess_liver %>% arrange(method, condition, sample_id) 

# range of the per-sample median ESS
ess_liver %>%
  group_by(method) %>%
  summarise(ess_median_min = min(ess_median), ess_median_max = max(ess_median), .groups = "drop")

# replicated pairs in the healthy samples 
plot_df_liver %>%
  filter(condition == "healthy") %>%
  group_by(method, gene1, gene2) %>%
  summarise(same_sign = all(r > 0) | all(r < 0),
            rep       = all(fdr < fdr_cut) & same_sign,
            rep_mod   = all(fdr_mod < fdr_cut) & same_sign, .groups = "drop") %>%
  group_by(method) %>%
  summarise(replicated = sum(rep), replicated_modified = sum(rep_mod), .groups = "drop")

# (D) effective sample size as a proportion of cells, all seven samples
p_ess_liver <- plot_df_liver %>%
  ggplot(aes(method, ess_prop, col = method)) +
  geom_hline(yintercept = 1, lty = 2, col = "grey60") +
  geom_boxplot(outlier.size = 0.3) +
  scale_color_manual(values = mycol) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = NULL, y = "Effective sample size (% of cells)") +
  theme_bw() +
  theme(legend.position = "none")

# (E) standard vs modified t-test, -log10 FDR
p_scatter_liver <- plot_df_liver %>%
  ggplot(aes(lp_naive, lp_modified, col = method)) +
  geom_point(size = 0.3, alpha = 0.3) +
  geom_abline(lty = 2, col = "grey40") +
  geom_hline(yintercept = -log10(fdr_cut), lty = 3, col = "grey60") +
  geom_vline(xintercept = -log10(fdr_cut), lty = 3, col = "grey60") +
  scale_x_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_y_continuous(trans = "log1p", breaks = c(0, 1, 3, 10, 30, 100, 300)) +
  scale_color_manual(values = mycol) +
  facet_wrap(~ method) +
  labs(x = expression(-log[10]~"FDR, standard t-test"),
       y = expression(-log[10]~"FDR, modified t-test")) +
  theme_bw() +
  theme(legend.position = "none", strip.background = element_blank())


# 4. One figure ----------------------------------------------------------------------
# shared legend on top, then one row per analysis with a title
leg <- get_plot_component(p_modsim, "guide-box-top", return_all = TRUE)

row_sim   <- plot_grid(p_modsim + theme(legend.position = "none"), labels = "A")
row_lung  <- plot_grid(p_ess, p_scatter, labels = c("B", "C"), nrow = 1,
                       rel_widths = c(1.5, 2.5), align = "h", axis = "tb")
row_liver <- plot_grid(p_ess_liver, p_scatter_liver, labels = c("D", "E"), nrow = 1,
                       rel_widths = c(1.5, 2.5), align = "h", axis = "tb")

p_all <- plot_grid(leg, row_lung, row_liver,row_sim,
                   ncol = 1, rel_heights = c(0.08,  1.1,  1.1,  1))
p_all
ggsave(file.path(dir.sim, "ModTtest_ESS.jpeg"), p_all, height = 250, width = 260, units = "mm")

