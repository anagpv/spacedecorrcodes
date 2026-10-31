# figure_block.R - supplementary block co-expression figure from block_results.csv.gz
# (run_block.R) -> supp_block.pdf
#   A  the true block correlation matrix of seed 1, ordered by the clustering that
#      defines the block labels
#   B  expression maps of seed 1: genes 1-2 under each model, zero mass 0% and 30%
#   C  precision-recall curves over BH-adjusted p-value thresholds (point: p < 0.05)
#   D  false discovery proportion per seed at BH-adjusted p < 0.05
#   E  adjusted Rand index per seed: hierarchical clustering of the estimated
#      matrix cut at 4 blocks against the true blocks

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

library(tidyverse)
library(cowplot)
library(patchwork)
library(data.table)

suppressPackageStartupMessages(source(file.path(dir.proj, "util", "simulator_mat.R")))   # sim_block()
source(file.path(dir.proj, "util", "util.R"))                                             # mat_ari()
stopifnot(requireNamespace("ggh4x", quietly = TRUE))

fig_dir <- file.path(dir.out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

NUMBLOCKS <- 4
PI_ZERO   <- c(0, 0.3)
NUM_CELLS <- 2000
NUM_GENES <- 100
RANGE     <- 0.25
ILLUSTRATION_SEED <- 1

ALPHA       <- 0.05
ALPHA_COL   <- "#CB181D"
ALPHA_LABEL <- sprintf("Target FDR = %s", ALPHA)
Y_LIM_FDP   <- c(-0.07, 1)
PVAL_GRID   <- c(1e-20, 1e-10, 1e-6, 1e-4, 1e-3, 0.01, ALPHA, 0.1, 0.2, 0.5, 1, Inf)

METHODS <- tibble::tribble(
  ~model,        ~label,            ~col,
  "spacedecorr", "SpaceDecorr",     "#8C510A",
  "spMOCA",      "spMOCA",          "#FF7F00",
  "logpearson",  "log-normalized",  "#000000",
  "spark",       "VST-normalized",  "#7F7F7F",
  "meringue",    "MERINGUE",        "#08519C",
  "giotto",      "Giotto",          "#6BAED6")
method_cols <- setNames(METHODS$col, METHODS$label)
draw_levels <- c(setdiff(METHODS$label, "SpaceDecorr"), "SpaceDecorr")

SIMS <- c(block_Additive = "Additive", block_MV = "Matrix Variate")
block_sim <- function(sim) factor(SIMS[sim], levels = SIMS)
zero_f    <- function(pi_zero) factor(sprintf("Zero mass: %g%%", 100 * pi_zero), levels = sprintf("Zero mass: %g%%", 100 * PI_ZERO))

COR_COLS   <- c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#CB181D")
BLOCK_COLS <- setNames(RColorBrewer::brewer.pal(8, "Dark2")[seq_len(NUMBLOCKS)], seq_len(NUMBLOCKS))
EXPR_MID   <- 5
EXPR_LOW   <- c("#DECDC4", "#D8B2A2", "#DD9A79")
EXPR_HIGH  <- c("#D27953", "#B55231", "#91301C", "#721A1B", "#5D1023", "#4F1530")
expr_scale <- function(expr) {
  hi <- max(expr)
  stopifnot(hi > EXPR_MID)
  at <- c(seq(0, EXPR_MID, length.out = length(EXPR_LOW)),
          seq(EXPR_MID, hi, length.out = length(EXPR_HIGH) + 1)[-1])
  scale_colour_gradientn(colours = c(EXPR_LOW, EXPR_HIGH), values = at / hi, limits = c(0, hi))
}

FONT      <- "Arial"
BASE_SIZE <- 14
FIG_W     <- 14
HEAT_IN   <- 3.15                                  # side of A's matrix = height of B's two rows of maps
MAP_GAP   <- 0.2
MAP_IN    <- (HEAT_IN - MAP_GAP) / 2
PANEL_H   <- 2.5                                   # panels of C, D, E
TREE_W    <- 14                                    # A's dendrogram, bar and gap, in matrix cells
BAR_W     <- 3
BAR_GAP   <- 0.8
A_LEFT    <- TREE_W + BAR_W + BAR_GAP + 0.5
PANE_GAP  <- 0.3

set_null_device(function(width, height) cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height))
with_cairo <- function(expr) {
  cairo_pdf(tempfile(fileext = ".pdf"), width = FIG_W, height = FIG_W)
  on.exit(dev.off())
  force(expr)
}

# Data: per-pair rows -> per-seed metrics and precision-recall curves ----
raw <- fread(file.path(dir.data, "block_results.csv.gz"), showProgress = FALSE) %>%
  as_tibble() %>%
  filter(model %in% METHODS$model, numblocks == NUMBLOCKS, pi_zero %in% PI_ZERO) %>%
  mutate(model = factor(METHODS$label[match(model, METHODS$model)], levels = METHODS$label),
         sim   = block_sim(sim), zero = zero_f(pi_zero))
stopifnot(nrow(raw) %% choose(NUM_GENES, 2) == 0,
          !any(duplicated(raw[, c("sim", "zero", "seed", "model", "gene_i", "gene_j")])))

# detection at padj < thr; a missing p-value is never a discovery
confusion <- function(d, thr)
  d %>%
    mutate(pred = coalesce(padj, 1) < thr) %>%
    group_by(sim, zero, seed, model) %>%
    summarise(tp = sum(pred & true_edge0 == 1, na.rm = TRUE),
              fp = sum(pred & true_edge0 == 0, na.rm = TRUE),
              fn = sum(!pred & true_edge0 == 1, na.rm = TRUE), .groups = "drop")

ari_of_rows <- function(d) {
  genes <- unique(c(d$gene_i, d$gene_j))
  M <- matrix(NA_real_, length(genes), length(genes), dimnames = list(genes, genes))
  M[cbind(d$gene_i, d$gene_j)] <- d$cor
  M[cbind(d$gene_j, d$gene_i)] <- d$cor
  diag(M) <- 1
  lab <- c(d$block_i, d$block_j)[match(genes, c(d$gene_i, d$gene_j))]
  mat_ari(M, lab, length(unique(lab)))
}

metrics <- inner_join(
  confusion(raw, ALPHA) %>% mutate(fdp = ifelse(tp + fp == 0, 0, fp / (tp + fp))),
  raw %>% group_by(sim, zero, seed, model) %>%
    summarise(ari = ari_of_rows(pick(gene_i, gene_j, block_i, block_j, cor)), .groups = "drop"),
  by = c("sim", "zero", "seed", "model"))

curves <- bind_rows(lapply(PVAL_GRID, function(thr) mutate(confusion(raw, thr), threshold = thr))) %>%
  mutate(precision = ifelse(tp + fp == 0, NA_real_, tp / (tp + fp)),
         recall    = ifelse(tp + fn == 0, NA_real_, tp / (tp + fn))) %>%
  group_by(sim, zero, model, threshold) %>%
  summarise(mean_precision = mean(precision, na.rm = TRUE), mean_recall = mean(recall, na.rm = TRUE), .groups = "drop") %>%
  filter(is.finite(mean_precision), is.finite(mean_recall)) %>%
  arrange(sim, zero, model, threshold)

# Panels ----
theme_block <- theme_bw(base_size = BASE_SIZE, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text       = element_text(size = BASE_SIZE),
        legend.position  = "none",
        plot.tag         = element_text(face = "bold", size = BASE_SIZE + 2))
facet_model_zero <- function() ggh4x::facet_nested(cols = vars(sim, zero))

dendro_segments <- function(hc) {
  n <- length(hc$order)
  leaf_pos <- match(seq_len(n), hc$order)
  node_pos <- numeric(n - 1)
  bind_rows(lapply(seq_len(n - 1), function(i) {
    child <- hc$merge[i, ]
    p  <- ifelse(child < 0, leaf_pos[abs(child)], node_pos[abs(child)])
    h0 <- ifelse(child < 0, 0, hc$height[abs(child)])
    node_pos[i] <<- mean(p)
    tibble(pos = c(p, p[1]), pos_end = c(p, p[2]), h = c(h0, hc$height[i]), h_end = hc$height[i])
  }))
}

## A: the true matrix ordered by the truth clustering, block bar and dendrogram to its left
panel_A <- function(truth) {
  hc <- hclust(dist(truth$delta_raw))
  stopifnot(identical(unname(cutree(hc, NUMBLOCKS)), truth$labels))
  n  <- length(hc$order)
  M  <- pmax(truth$delta0, 0)[hc$order, hc$order]
  bar_right <- 0.5 - BAR_GAP
  bar_left  <- bar_right - BAR_W
  d_heat <- tibble(row = rep(seq_len(n), n), col = rep(seq_len(n), each = n), cor = as.vector(M))
  d_bar  <- tibble(row = seq_len(n), block = truth$labels[hc$order])
  d_tree <- dendro_segments(hc) %>% mutate(across(c(h, h_end), ~ bar_left - .x / max(hc$height) * TREE_W))

  ggplot() +
    geom_tile(aes(col, n + 1 - row, fill = cor, colour = cor), d_heat, linewidth = 0.1) +
    geom_tile(aes((bar_left + bar_right) / 2, n + 1 - row), d_bar, width = BAR_W,
              fill = BLOCK_COLS[d_bar$block], colour = BLOCK_COLS[d_bar$block], linewidth = 0.1) +
    geom_segment(aes(h, n + 1 - pos, xend = h_end, yend = n + 1 - pos_end), d_tree, linewidth = 0.3, lineend = "square") +
    scale_fill_gradientn(colours = COR_COLS, limits = c(0, 1), aesthetics = c("fill", "colour")) +
    coord_cartesian(xlim = c(0.5 - A_LEFT, n + 0.5), ylim = c(0.5, n + 0.5), expand = FALSE) +
    labs(tag = "A") +
    theme_void(base_size = BASE_SIZE, base_family = FONT) +
    theme(legend.position = "none", plot.tag = theme_block$plot.tag)
}

## B: expression maps, log1p of the counts scaled to the median library size
panel_B <- function(sims, grid) {
  genes  <- 1:2
  gene_f <- factor(paste("Gene", genes), levels = paste("Gene", genes))
  d <- bind_rows(lapply(seq_along(sims), function(i) {
    Y   <- as.matrix(sims[[i]]$Y[, genes, drop = FALSE])
    lib <- sims[[i]]$data$libsize
    co  <- sims[[i]]$coordinates
    tibble(zero = zero_f(grid$pi_zero[i]), sim = block_sim(grid$sim[i]), gene = rep(gene_f, each = nrow(Y)),
           sdimx = rep(co$sdimx, length(genes)), sdimy = rep(co$sdimy, length(genes)),
           expr  = as.vector(log1p(Y / lib * median(lib))))
  })) %>% arrange(expr)

  ggplot(d, aes(sdimx, sdimy, colour = expr)) +
    geom_point(size = 2.1, alpha = 0.85, stroke = 0) +
    expr_scale(d$expr) +
    ggh4x::facet_nested(rows = vars(zero), cols = vars(sim, gene)) +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    scale_x_continuous(breaks = c(-1, 0, 1), labels = c("-1", "0", "1")) +
    scale_y_continuous(breaks = c(-1, 0, 1), labels = c("-1", "0", "1")) +
    labs(x = "sdimx", y = "sdimy", colour = "Expression", tag = "B") +
    theme_block +
    theme(legend.position = "right", panel.grid = element_blank(),
          legend.title.position = "right",
          legend.title = element_text(angle = 270, hjust = 0.5, margin = margin(l = 4)),
          panel.background = element_rect(fill = EXPR_LOW[1], colour = NA),
          axis.text = element_text(size = BASE_SIZE - 4), panel.spacing = unit(MAP_GAP, "in"))
}

## C: precision-recall curves, the p < ALPHA point marked
panel_C <- function(curves) {
  curves <- curves %>% mutate(model = factor(model, levels = draw_levels)) %>% arrange(model, threshold)
  ggplot(curves, aes(mean_recall, mean_precision, colour = model, group = model)) +
    geom_path(linewidth = 0.6) +
    geom_point(data = filter(curves, threshold == ALPHA), size = 2.5) +
    facet_model_zero() +
    scale_colour_manual(values = method_cols, breaks = METHODS$label, drop = FALSE) +
    scale_x_continuous("Recall",    limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    scale_y_continuous("Precision", limits = c(0, 1), breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    labs(colour = NULL, tag = "C") +
    theme_block +
    theme(panel.spacing.x = unit(PANE_GAP, "in"), panel.heights = unit(PANEL_H, "in"))
}

## D, E: one box per method in each subplot
model_box <- function(d, yvar, ylab, tag, ylim = NULL, href = NULL) {
  p <- ggplot(d, aes(model, .data[[yvar]], colour = model))
  if (!is.null(href))
    p <- p + geom_hline(yintercept = href, colour = ALPHA_COL, linetype = "dashed", linewidth = 0.4)
  p <- p +
    geom_boxplot(outlier.size = 0.7, linewidth = 0.5) +
    facet_model_zero() +
    scale_colour_manual(values = method_cols, breaks = METHODS$label, drop = FALSE) +
    labs(x = NULL, y = ylab, colour = NULL, tag = tag) +
    theme_block +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.ticks.length.x = unit(0, "pt"),
          panel.spacing.x = unit(PANE_GAP, "in"), panel.heights = unit(PANEL_H, "in"))
  if (!is.null(ylim)) p <- p + coord_cartesian(ylim = ylim)
  if (!is.null(href))
    p <- p + annotate("text", x = n_distinct(d$model) + 0.5, y = 0, label = ALPHA_LABEL, family = FONT,
                      hjust = 1, vjust = 1.1, size = 12 / .pt, colour = ALPHA_COL)
  p
}

# Assemble and save ----
grid <- expand_grid(pi_zero = PI_ZERO, sim = names(SIMS))
sims <- Map(function(pz, s)
  sim_block(rept = ILLUSTRATION_SEED, sim = s, num_cells = NUM_CELLS, num_genes = NUM_GENES,
            range = RANGE, numblocks = NUMBLOCKS, pi_zero = pz),
  grid$pi_zero, grid$sim)

pA <- panel_A(sims[[1]]$truth)
pB <- panel_B(sims, grid)
pC <- panel_C(curves)
pD <- model_box(metrics, "fdp", "FDP", "D", ylim = Y_LIM_FDP, href = ALPHA)
pE <- model_box(metrics, "ari", "ARI", "E", ylim = c(0, 1))

legend_src <- pC +
  guides(colour = guide_legend(nrow = 1, title = "Methods")) +
  theme(legend.position = "top",
        legend.title = element_text(size = BASE_SIZE), legend.text = element_text(size = BASE_SIZE),
        legend.key.width = unit(2, "line"), legend.spacing.x = unit(1, "line"))
legend <- get_plot_component(legend_src, "guide-box-top", return_all = TRUE)

a_w    <- HEAT_IN * (NUM_GENES + A_LEFT) / NUM_GENES
b_w    <- 4 * MAP_IN + 3 * MAP_GAP
row_AB <- pA + plot_spacer() + pB +
  plot_layout(widths = grid::unit.c(unit(a_w, "in"), unit(1, "null"), unit(b_w, "in")), heights = unit(HEAT_IN, "in"))

rows <- with_cairo(list(legend, patchworkGrob(row_AB), ggplotGrob(pC), ggplotGrob(pD), ggplotGrob(pE)))
hs   <- with_cairo(vapply(rows, function(g) grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE), numeric(1)))
fig  <- plot_grid(plotlist = rows, ncol = 1, rel_heights = hs)
ggsave(file.path(fig_dir, "supp_block.pdf"), fig, width = FIG_W, height = sum(hs), device = cairo_pdf)
message("wrote supp_block.pdf")
