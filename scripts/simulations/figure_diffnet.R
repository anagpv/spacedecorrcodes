# figure_diffnet.R - supplementary differential-network figure -> supp_diffnet.pdf
#   A  blank (filled in the layout software)
#   B  two real niches (diffnet_niche_results.csv, diffnet_niche_spacex_results.csv):
#      precision-recall curves per network
#   C  two real niches: false discovery proportion per seed at BH-adjusted p < 0.05
#   D  two simulated tissues (diffnet_samples_results.csv): precision-recall curves
#   E  two simulated tissues: false discovery proportion per seed

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

library(tidyverse)
library(cowplot)
library(patchwork)

fig_dir <- file.path(dir.out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

ALPHA       <- 0.05
ALPHA_COL   <- "#CB181D"
ALPHA_LABEL <- sprintf("Target FDR = %s", ALPHA)
Y_LIM_FDP   <- c(-0.07, 1)
PVAL_GRID   <- c(1e-20, 1e-10, 1e-6, 1e-4, 1e-3, 0.01, ALPHA, 0.1, 0.2, 0.5, 1, Inf)   # curve thresholds

METHODS <- tibble::tribble(
  ~model,        ~label,           ~col,      ~legend,
  "spacedecorr", "SpaceDecorr",    "#8C510A", "SpaceDecorr",
  "spMOCA",      "spMOCA",         "#FF7F00", "spMOCA",
  "spacex",      "SpaceX",         "#E6AB02", "SpaceX (only in panels B, C)",
  "logpearson",  "log-normalized", "#000000", "log-normalized",
  "spark",       "VST-normalized", "#7F7F7F", "VST-normalized",
  "meringue",    "MERINGUE",       "#08519C", "MERINGUE",
  "giotto",      "Giotto",         "#6BAED6", "Giotto")
method_cols <- setNames(METHODS$col, METHODS$label)
draw_levels <- c(setdiff(METHODS$label, "SpaceDecorr"), "SpaceDecorr")

SIMS <- c(sample_Additive = "Additive", sample_MV = "Matrix Variate")

FONT        <- "Arial"
BASE_SIZE   <- 14
LEGEND_PT   <- BASE_SIZE + 1
FIG_W       <- 14
ROW_HEIGHTS <- c(1, 1.5)
ROW_IN      <- 4.5
DESIGN <- "
AABBBBBBBCCCCC
DDDDDDDDEEEEEE
"
set_null_device(function(width, height) cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height))

# Data: per-pair rows -> per-seed metrics and precision-recall curves ----
read_results <- function(name)
  read.csv(file.path(dir.data, name), stringsAsFactors = FALSE) %>%
    filter(model %in% METHODS$model) %>%
    mutate(model = factor(METHODS$label[match(model, METHODS$model)], levels = METHODS$label))

NETWORK_SPEC <- tibble::tribble(
  ~network, ~true,        ~est,       ~padj,
  "g0",     "true_edge0", "cor_g0",   "padj_g0",
  "g1",     "true_edge1", "cor_g1",   "padj_g1",
  "diff",   "true_diff",  "cor_diff", "padj_diff")

pairs_long <- function(raw)
  bind_rows(lapply(seq_len(nrow(NETWORK_SPEC)), function(i) {
    s <- NETWORK_SPEC[i, ]
    tibble(sim = raw$sim, seed = raw$seed, model = raw$model, network = s$network,
           true = raw[[s$true]], est = raw[[s$est]], padj = raw[[s$padj]])
  }))

# detection at padj < thr; a missing p-value is never a discovery
confusion <- function(d, thr)
  d %>%
    mutate(pred = coalesce(padj, 1) < thr) %>%
    group_by(sim, seed, model, network) %>%
    summarise(tp = sum(pred & true == 1, na.rm = TRUE),
              fp = sum(pred & true == 0, na.rm = TRUE),
              fn = sum(!pred & true == 1, na.rm = TRUE), .groups = "drop")

seed_fdp <- function(d)
  confusion(d, ALPHA) %>% mutate(fdp = ifelse(tp + fp == 0, 0, fp / (tp + fp)))

pr_curves <- function(d, graph_of)
  bind_rows(lapply(PVAL_GRID, function(thr) mutate(confusion(d, thr), threshold = thr))) %>%
    mutate(precision = ifelse(tp + fp == 0, NA_real_, tp / (tp + fp)),
           recall    = ifelse(tp + fn == 0, NA_real_, tp / (tp + fn)),
           graph     = graph_of(network)) %>%
    group_by(sim, model, graph, threshold) %>%
    summarise(mean_precision = mean(precision, na.rm = TRUE),
              mean_recall    = mean(recall,    na.rm = TRUE), .groups = "drop") %>%
    filter(is.finite(mean_precision), is.finite(mean_recall)) %>%
    arrange(sim, model, graph, threshold)

n_pairs <- choose(20, 2)

## niches (B, C)
raw_niche <- bind_rows(read_results("diffnet_niche_results.csv"), read_results("diffnet_niche_spacex_results.csv"))
stopifnot(nrow(raw_niche) %% n_pairs == 0,
          !any(duplicated(raw_niche[, c("sim", "seed", "model", "gene_i", "gene_j")])))
n0 <- unique(raw_niche$n_g0); n1 <- unique(raw_niche$n_g1)
stopifnot(length(n0) == 1, length(n1) == 1)

niche_graph_levels <- c(sprintf("Myeloid-enriched\nstroma (n=%s)", format(n1, big.mark = ",")),
                        sprintf("Stroma\n(n=%s)",                  format(n0, big.mark = ",")),
                        "Difference")
niche_graph <- function(network)
  factor(unname(c(g1 = niche_graph_levels[1], g0 = niche_graph_levels[2], diff = niche_graph_levels[3])[network]),
         levels = niche_graph_levels)

long_niche    <- pairs_long(raw_niche)
metrics_niche <- seed_fdp(long_niche) %>% mutate(graph = niche_graph(network))
curves_niche  <- pr_curves(long_niche, niche_graph)

## simulated tissues (D, E): the two tissues pooled as "Tissues"
raw_samples <- read_results("diffnet_samples_results.csv")
stopifnot(nrow(raw_samples) %% n_pairs == 0, all(raw_samples$sim %in% names(SIMS)),
          !any(duplicated(raw_samples[, c("sim", "seed", "model", "gene_i", "gene_j")])))
samples_graph <- function(network)
  factor(ifelse(network == "diff", "Difference", "Tissues"), levels = c("Tissues", "Difference"))
samples_sim <- function(sim) factor(SIMS[sim], levels = SIMS)

long_samples    <- pairs_long(raw_samples)
metrics_samples <- seed_fdp(long_samples) %>% mutate(graph = samples_graph(network), sim = samples_sim(sim))
curves_samples  <- pr_curves(long_samples, samples_graph) %>% mutate(sim = samples_sim(sim))

# Panels ----
theme_diffnet <- theme_bw(base_size = BASE_SIZE, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text       = element_text(size = BASE_SIZE),
        legend.position  = "none",
        plot.tag         = element_text(face = "bold", size = BASE_SIZE + 2))

fdp_panel <- function(d, tag, facets = NULL) {
  p <- ggplot(d, aes(graph, fdp, colour = model)) +
    geom_hline(yintercept = ALPHA, colour = ALPHA_COL, linetype = "dashed", linewidth = 0.4) +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.45) +
    scale_colour_manual(values = method_cols, breaks = METHODS$label, labels = METHODS$legend, drop = FALSE) +
    labs(x = "Network", y = "FDP", colour = "Methods", tag = tag) +
    coord_cartesian(ylim = Y_LIM_FDP) +
    annotate("text", x = n_distinct(d$graph) + 0.5, y = 0, label = ALPHA_LABEL, family = FONT,
             hjust = 1, vjust = 1.1, size = 12 / .pt, colour = ALPHA_COL) +
    theme_diffnet
  if (!is.null(facets)) p <- p + facet_grid(facets)
  p
}

pr_panel <- function(curves, tag, facets) {
  curves <- curves %>% mutate(model = factor(model, levels = draw_levels)) %>% arrange(sim, model, graph, threshold)
  ggplot(curves, aes(mean_recall, mean_precision, colour = model, group = model)) +
    geom_path(linewidth = 0.6) +
    geom_point(data = filter(curves, threshold == ALPHA), size = 2.5, shape = 16) +
    scale_colour_manual(values = method_cols, breaks = METHODS$label, labels = METHODS$legend, drop = FALSE) +
    scale_x_continuous("Recall",    limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous("Precision", limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    facet_grid(facets) +
    labs(colour = "Methods", tag = tag) +
    theme_diffnet
}

pA <- ggplot() + theme_void(base_size = BASE_SIZE, base_family = FONT) + labs(tag = "A") +
  theme(plot.tag = theme_diffnet$plot.tag)
pB <- pr_panel(curves_niche, "B", . ~ graph)
pC <- fdp_panel(metrics_niche, "C")
pD <- pr_panel(curves_samples, "D", sim ~ graph)
pE <- fdp_panel(metrics_samples, "E", facets = . ~ sim)

# Legend, assemble, save ----
legend_src <- pB +
  guides(colour = guide_legend(nrow = 1)) +
  theme(legend.position = "top",
        legend.title = element_text(size = LEGEND_PT, margin = margin(r = 6)),
        legend.text  = element_text(size = LEGEND_PT),
        legend.key.width = unit(1.2, "line"), legend.spacing.x = unit(0.4, "line"),
        legend.box.spacing = unit(0, "pt"))
legend   <- get_plot_component(legend_src, "guide-box-top", return_all = TRUE)
legend_h <- round(grid::convertHeight(sum(legend$heights), "in", valueOnly = TRUE) + 0.1, 2)
fig_h    <- legend_h + ROW_IN * sum(ROW_HEIGHTS)

body <- pA + pB + pC + pD + pE + plot_layout(design = DESIGN, heights = ROW_HEIGHTS)
fig  <- plot_grid(legend, body, ncol = 1, rel_heights = c(legend_h, fig_h - legend_h))
ggsave(file.path(fig_dir, "supp_diffnet.pdf"), fig, width = FIG_W, height = fig_h, device = cairo_pdf)
message("wrote supp_diffnet.pdf")
