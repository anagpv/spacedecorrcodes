#-------------------------#
#- Liver merfish dataset -#
#-------------------------#
library(tidyverse)
library(ggrepel)
library(cowplot)

dir.data <- "/path/to/zenodo/download/" 
dir.out  <- "/path/to/output/"          
dir.proj <- "/path/to/repo/"       

tissues    <- c("AM048", "AM061", "AM042", "AM031", "AM066", "AM062", "AM072")
conditions <- c("healthy", "healthy", "healthy", "diseased", "diseased", "diseased", "diseased")
names(conditions) <- tissues
hep_types <- c("Central Hep", "Portal Hep", "Fibrotic Hep 1", "Fibrotic Hep 2")

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
genes_keep <- rownames(mean_expr)[apply(mean_expr, 1, min) >= 0.10]
length(genes_keep) # 101

# Number of cells in each dataset
n_cells <- meta_data %>%
  filter(Cell_Type_final %in% hep_types) %>%
  group_by(sample_id) %>%
  summarise(n_hep     = n(),
            n_zonated = sum(Cell_Type_final %in% c("Central Hep", "Portal Hep")), .groups = "drop")
n_cells
n_cells %>% summarise(across(c(n_hep, n_zonated), list(mean = mean, min = min, max = max)))

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

# Spatial R2 per gene ----------------------------------------------
spatial_r2 <- do.call(bind_rows, lapply(tissues, function(tissue){
  print(tissue)
  
  # SpaceDecorr
  f <- file.path(dir.out, tissue, paste0("Data_", tissue, "_spacedecorr_k200_hepfibcov.Rdata"))
  while (grepl("dataless", system(paste("ls -lO", shQuote(f)), intern = TRUE))) Sys.sleep(5)
  load(f)
  spl <- as.matrix(data_spacedecorr_hepcov$splines)   
  rownames(spl) <- rownames(data_spacedecorr_hepcov$metadata)
  res <- as.matrix(data_spacedecorr_hepcov$residuals)
  rownames(res) <- rownames(data_spacedecorr_hepcov$metadata)
  
  # VST
  load(file.path(dir.out, tissue, paste0("Data_", tissue, "_vst_hepfib.Rdata")))
  vst <- t(as.matrix(data_vst_hep$y))                 # cells x genes
  
  # hepatocytes and genes present in all objects
  hep   <- meta_data$index[meta_data$sample_id == tissue & meta_data$Cell_Type_final %in% hep_types]
  cells <- Reduce(intersect, list(hep, rownames(res), rownames(vst)))
  genes <- Reduce(intersect, list(genes_keep, colnames(res), colnames(vst)))
  
  lognorm <- as.matrix(log1p(raw[cells, genes] / N[cells] * 1e6))
  
  # spatial R2: squared correlation of log-normalized expression with the fitted spline
  r2 <- diag(cor(lognorm, spl[cells, genes]))^2
  data.frame(sample_id = tissue, 
             condition = unname(conditions[tissue]),
             gene = genes, r2 = r2)
}))
write.csv(spatial_r2, file = file.path(dir.out, "spatial_r2_hepfib.csv"), row.names = FALSE)

# Sensitivity to the spline flexibility K -------------------------------------
Ks        <- c(100, 200, 1000)
tissues_h <- tissues[conditions == "healthy"]
fdr_cut   <- 0.05

sens <- lapply(Ks, function(K) {
  lapply(tissues_h, function(tissue) {
    print(c(K, tissue))
    f <- file.path(dir.out, tissue, paste0("Data_", tissue, "_spacedecorr_k", K, "_hepfibcov.Rdata"))
    load(f)
    obj <- data_spacedecorr_hepcov                      
    spl <- as.matrix(obj$splines);   rownames(spl) <- rownames(obj$metadata)
    res <- as.matrix(obj$residuals); rownames(res) <- rownames(obj$metadata)
    
    cells <- meta_data$index[meta_data$sample_id == tissue &
                               meta_data$Cell_Type_final %in% c("Central Hep", "Portal Hep")]
    cells <- intersect(cells, rownames(res))
    genes <- intersect(genes_keep, colnames(res))
    
    # edf per gene (data.frame, 1 x genes)
    e   <- setNames(as.numeric(obj$edf[1, ]), colnames(obj$edf))
    edf <- data.frame(K = K, sample_id = tissue, gene = genes, edf = e[genes])
    
    # spatial R2 per gene
    lognorm <- as.matrix(log1p(raw[cells, genes] / N[cells] * 1e6))
    r2 <- data.frame(K = K, sample_id = tissue, gene = genes,
                     r2 = diag(cor(lognorm, spl[cells, genes]))^2)
    
    # residual correlations with naive P-values
    C  <- cor(res[cells, genes], use = "pairwise.complete.obs")
    ut <- upper.tri(C); ij <- which(ut, arr.ind = TRUE)
    n  <- length(cells)
    pc <- data.frame(K = K, sample_id = tissue, gene1 = genes[ij[, 1]], gene2 = genes[ij[, 2]],
                     r = C[ut]) %>%
      mutate(p = 2 * pt(-abs(r) * sqrt((n - 2) / (1 - r^2)), n - 2),
             fdr = p.adjust(p, "BH"))
    list(r2 = r2, pc = pc, edf = edf)                      # edf
  })
})
sens_r2  <- bind_rows(lapply(unlist(sens, recursive = FALSE), `[[`, "r2"))
sens_edf <- bind_rows(lapply(unlist(sens, recursive = FALSE), `[[`, "edf")) %>%     # edf
  mutate(program = case_when(gene %in% programs[[1]] ~ names(programs)[1],
                             gene %in% programs[[2]] ~ names(programs)[2],
                             gene %in% programs[[3]] ~ names(programs)[3],
                             TRUE ~ "Other"))
sens_pc <- bind_rows(lapply(unlist(sens, recursive = FALSE), `[[`, "pc")) %>%
  mutate(type = case_when(
    gene1 %in% programs[[1]] & gene2 %in% programs[[1]] ~ names(programs)[1],
    gene1 %in% programs[[2]] & gene2 %in% programs[[2]] ~ names(programs)[2],
    gene1 %in% programs[[3]] & gene2 %in% programs[[3]] ~ names(programs)[3],
    (gene1 %in% programs[[1]] & gene2 %in% programs[[2]]) |
      (gene1 %in% programs[[2]] & gene2 %in% programs[[1]]) ~ "Opposite zones",
    TRUE ~ "Other pairs"))
save(sens_r2, sens_pc, sens_edf, file = file.path(dir.out, "res_liver_K_sensitivity.Rdata"))

# Pairwise correlations per sample, per hepatocyte state and pooled ----------
groups <- list("Zonated Hep"    = c("Central Hep", "Portal Hep"),
               "Fibrotic Hep 1" = "Fibrotic Hep 1",
               "Fibrotic Hep 2" = "Fibrotic Hep 2",
               "Pooled"         = hep_types)
min_cells <- 20
source(file.path(dir.proj, "modified_cor_test.R"))

pair_cor <- do.call(bind_rows, lapply(tissues, function(tissue){
  print(tissue)
  
  # SpaceDecorr
  f <- file.path(dir.out, tissue, paste0("Data_", tissue, "_spacedecorr_k200_hepfibcov.Rdata"))
  load(f)
  res <- as.matrix(data_spacedecorr_hepcov$residuals)
  rownames(res) <- rownames(data_spacedecorr_hepcov$metadata)
  
  # VST
  f <- file.path(dir.out, tissue, paste0("Data_", tissue, "_vst_hepfib.Rdata"))
  load(f)
  vst <- t(as.matrix(data_vst_hep$y))                 # cells x genes
  
  genes <- Reduce(intersect, list(genes_keep, colnames(res), colnames(vst)))
  
  do.call(bind_rows, lapply(names(groups), function(grp){
    print(grp)
    cells <- meta_data$index[meta_data$sample_id == tissue & meta_data$Cell_Type_final %in% groups[[grp]]]
    cells <- Reduce(intersect, list(cells, rownames(res), rownames(vst)))
    if (length(cells) < min_cells) return(NULL)
    coords <- as.matrix(meta_data[match(cells, meta_data$index), c("x", "y")])
    
    mt_vst <- modified_cor_test(vst[cells, genes], coords)$results
    mt_sd  <- modified_cor_test(res[cells, genes], coords)$results
    
    data.frame(sample_id = tissue,
               condition = unname(conditions[tissue]),
               cell_type = grp,
               n_cells   = length(cells),
               gene1 = mt_vst$gene1, gene2 = mt_vst$gene2,
               vst_r            = mt_vst$r,
               vst_ess        = mt_vst$ess,
               vst_p          = mt_vst$p_naive,
               vst_pmod      = mt_vst$p_modified,
               spacedecorr_r       = mt_sd$r,
               spacedecorr_ess   = mt_sd$ess,
               spacedecorr_p     = mt_sd$p_naive,
               spacedecorr_pmod = mt_sd$p_modified)
  }))
}))

# annotate pairs by program
pair_cor <- pair_cor %>%
  mutate(type = case_when(
    gene1 %in% programs[[1]] & gene2 %in% programs[[1]] ~ names(programs)[1],
    gene1 %in% programs[[2]] & gene2 %in% programs[[2]] ~ names(programs)[2],
    gene1 %in% programs[[3]] & gene2 %in% programs[[3]] ~ names(programs)[3],
    (gene1 %in% programs[[1]] & gene2 %in% programs[[2]]) |
      (gene1 %in% programs[[2]] & gene2 %in% programs[[1]]) ~ "Opposite zones",
    TRUE ~ "Other pairs"))
save(pair_cor, file = file.path(dir.out, "pair_cor_liver_hepfib_modttest.Rdata"))

# Welch test between conditions on per-sample Fisher z -------------------------
colnames(pair_cor)[c(7,11)] <- c("vst","spacedecorr")

welch_res <- data.frame()
for (m in c("vst", "spacedecorr")) {
  for (ct in unique(pair_cor$cell_type)) {
    
    # pairs x samples matrix of Fisher z
    Z <- pair_cor %>%
      filter(cell_type == ct) %>%
      mutate(z = atanh(pmin(pmax(.data[[m]], -0.999), 0.999))) %>%
      dplyr::select(gene1, gene2, type, sample_id, z) %>%
      pivot_wider(names_from = sample_id, values_from = z)
    smp <- intersect(tissues, colnames(Z))                  # samples with this cell type
    Zm  <- as.matrix(Z[, smp])
    h   <- conditions[smp] == "healthy"
    n1  <- sum(h); n2 <- sum(!h)
    if (n1 < 2 | n2 < 2) next                               # need 2+ samples per group for a variance
    
    # Welch t (diseased - healthy)
    v1 <- apply(Zm[, h], 1, var) / n1
    v2 <- apply(Zm[, !h], 1, var) / n2
    t_obs <- (rowMeans(Zm[, !h]) - rowMeans(Zm[, h])) / sqrt(v1 + v2)
    df    <- (v1 + v2)^2 / (v1^2 / (n1 - 1) + v2^2 / (n2 - 1))
    p     <- 2 * pt(-abs(t_obs), df)
    
    # exact permutation over all relabelings of the available samples
    perms  <- combn(length(smp), n1)
    t_perm <- matrix(NA, nrow(Zm), ncol(perms))
    for (k in 1:ncol(perms)) {
      hk <- seq_along(smp) %in% perms[, k]
      t_perm[, k] <- (rowMeans(Zm[, !hk]) - rowMeans(Zm[, hk])) /
        sqrt(apply(Zm[, hk], 1, var) / n1 + apply(Zm[, !hk], 1, var) / n2)
    }
    
    separated <- apply(Zm[, !h, drop = FALSE], 1, max) < apply(Zm[, h, drop = FALSE], 1, min)
    
    welch_res <- bind_rows(welch_res, data.frame(
      method = m, cell_type = ct, n_healthy = n1, n_diseased = n2,
      gene1 = Z$gene1, gene2 = Z$gene2, type = Z$type,
      r_healthy = tanh(rowMeans(Zm[, h])), r_diseased = tanh(rowMeans(Zm[, !h])),
      diff_z = rowMeans(Zm[, !h]) - rowMeans(Zm[, h]),
      t = t_obs, df = df, p = p,
      fdr = p.adjust(p, "BH"),
      p_perm = rowMeans(abs(t_perm) >= abs(t_obs) - 1e-12),
      min_p_perm = 1 / ncol(perms),
      separated = separated))
  }
}

welch_res <- welch_res %>%
  mutate(stronger_in = ifelse(abs(r_healthy) > abs(r_diseased), "Healthy", "Diseased"))

save(welch_res, file = file.path(dir.out, "res_liver_welch_by_celltype.Rdata"))


# Leave one out -----------------------------------------------------------
r_col   <- "spacedecorr"         # correlation column in pair_cor (adjusted SpaceDecorr)
p_col   <- "spacedecorr_p"    # P-value column used for the replicated network
fdr_rep <- 0.05                     # within-sample FDR for replication
fdr_dc  <- 0.10                     # FDR for differences between conditions
injury  <- programs[["Injury program"]]

loo_dc <- data.frame()
for (m in c("vst", r_col)) {
  Zfull <- pair_cor %>%
    filter(cell_type == "Zonated Hep") %>%
    mutate(z = atanh(pmin(pmax(.data[[m]], -0.999), 0.999))) %>%
    dplyr::select(gene1, gene2, sample_id, z) %>%
    pivot_wider(names_from = sample_id, values_from = z)
  
  for (drop in c("none", tissues)) {
    smp <- setdiff(intersect(tissues, colnames(Zfull)), drop)
    Zm  <- as.matrix(Zfull[, smp])
    h   <- conditions[smp] == "healthy"
    n1  <- sum(h); n2 <- sum(!h)
    v1  <- apply(Zm[, h], 1, var) / n1;  v2 <- apply(Zm[, !h], 1, var) / n2
    t_obs <- (rowMeans(Zm[, !h]) - rowMeans(Zm[, h])) / sqrt(v1 + v2)
    df    <- (v1 + v2)^2 / (v1^2 / (n1 - 1) + v2^2 / (n2 - 1))
    fdr   <- p.adjust(2 * pt(-abs(t_obs), df), "BH")
    loo_dc <- bind_rows(loo_dc, data.frame(method = m, dropped = drop,
                                           pair = paste(Zfull$gene1, Zfull$gene2), sig = fdr < fdr_dc,
                                           separated = separated,
                                           lower_fibrotic = abs(tanh(rowMeans(Zm[, !h]))) < abs(tanh(rowMeans(Zm[, h])))))
  }
}

# pairs significant with all donors, per method
full_sig <- loo_dc %>%
  filter(dropped == "none", sig) %>%
  dplyr::select(method, pair)

# % of those still significant after removing each donor
loo_dc_summary <- loo_dc %>%
  semi_join(full_sig, by = c("method", "pair")) %>%
  group_by(method, dropped) %>%
  summarise(retained = mean(sig), .groups = "drop")

loo_table <- loo_dc %>%
  group_by(method, dropped) %>%
  summarise(n_sig = sum(sig),
            pct_lower = ifelse(sum(sig) > 0, 100 * mean(lower_fibrotic[sig]), NA), .groups = "drop") %>%
  left_join(loo_dc_summary %>% dplyr::select(method, dropped, retained), by = c("method", "dropped")) %>%
  mutate(retained = ifelse(dropped == "none", NA, 100 * retained),
         method = ifelse(method == "vst", "vst", "sd")) %>%
  pivot_wider(names_from = method, values_from = c(n_sig, retained, pct_lower)) %>%
  mutate(condition = ifelse(dropped == "none", "", unname(conditions[dropped])),
         dropped = ifelse(dropped == "none", "None", dropped),
         dropped = factor(dropped, c("None", tissues[conditions == "healthy"], tissues[conditions == "diseased"]))) %>%
  arrange(dropped) %>%
  dplyr::select(`Donor removed` = dropped, Condition = condition,
                `VST: pairs` = n_sig_vst, `VST: retained (%)` = retained_vst, `VST: lower in fibrotic (%)` = pct_lower_vst,
                `SpaceDecorr: pairs` = n_sig_sd, `SpaceDecorr: retained (%)` = retained_sd,
                `SpaceDecorr: lower in fibrotic (%)` = pct_lower_sd)

write.csv(loo_table, file.path(dir.out, "loo_between_tissue.csv"), row.names = FALSE)


# Time --------------------------------------------------------------------
find_file <- function(tissue, fn) {
  f <- c(file.path(dir.out, tissue, fn), file.path(dir.out, fn))
  f[file.exists(f)][1]
}

times <- data.frame()
for (tissue in tissues) {
  
  # SpaceDecorr, hepatocytes with cell-state covariate (system.time, 5 cores)
  for (K in c(100, 200, 1000)) {
    f <- find_file(tissue, paste0("Data_", tissue, "_spacedecorr_k", K, "_hepfibcov.Rdata"))
    if (is.na(f)) next
    e <- new.env(); load(f, envir = e)
    times <- bind_rows(times, data.frame(tissue = tissue, method = "SpaceDecorr", K = K,
                                         n_cells = nrow(e$data_spacedecorr_hepcov$residuals),
                                         sec = unname(e$time_spacedecorr_hepcov["elapsed"])))
  }
  
  # VST (seconds)
  f <- find_file(tissue, paste0("Data_", tissue, "_vst_hepfib.Rdata"))
  if (!is.na(f)) {
    e <- new.env(); load(f, envir = e)
    times <- bind_rows(times, data.frame(tissue = tissue, method = "VST", K = NA,
                                         n_cells = ncol(e$data_vst_hep$y),
                                         sec = e$time_vst_hep))
  }
  
  # spMOCA (elapsed seconds)
  f <- find_file(tissue, paste0("Data_", tissue, "_spmoca_hep.Rdata"))
  if (!is.na(f)) {
    e <- new.env(); load(f, envir = e)
    times <- bind_rows(times, data.frame(tissue = tissue, method = "spMOCA", K = NA,
                                         n_cells = NA,
                                         sec = e$time_spmoca_hep))
  }
}

# fill n_cells for spMOCA from the SpaceDecorr fit of the same tissue (same cells)
times <- times %>%
  group_by(tissue) %>%
  mutate(n_cells = ifelse(is.na(n_cells), max(n_cells, na.rm = TRUE), n_cells)) %>%
  ungroup() %>%
  mutate(condition = conditions[tissue], min = sec / 60)

# per tissue
times %>%
  mutate(method = ifelse(method == "SpaceDecorr", paste0("SpaceDecorr_k", K), method)) %>%
  dplyr::select(tissue, n_cells, method, min) %>%
  pivot_wider(names_from = method, values_from = min) %>%
  arrange(n_cells) %>%
  print(width = Inf)

# summary by method
times %>%
  group_by(method, K) %>%
  summarise(n_tissues = n(), median_min = median(min), min_min = min(min), max_min = max(min),
            .groups = "drop")

write.csv(times, file.path(dir.out, "times_liver.csv"), row.names = FALSE)

# time vs number of cells
times %>%
  mutate(method = ifelse(method == "SpaceDecorr", paste0("SpaceDecorr (K = ", K, ")"), method)) %>%
  ggplot(aes(n_cells, min, col = method)) +
  geom_point() + geom_line() +
  scale_x_log10(labels = scales::comma) + scale_y_log10() +
  labs(x = "Number of hepatocytes", y = "Time (minutes)", col = NULL) +
  theme_bw() +
  theme(legend.position = "top", panel.grid.minor = element_blank())



