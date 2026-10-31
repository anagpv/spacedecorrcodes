# figure_pairwise.R - Scenario I figures from pairwise_results.csv (run_pairwise.R):
#   main figure, panel A   type 1 error vs number of cells: spatial covariance
#                          (Additive), FOV technical factor (FOV_diff), unobserved
#                          subtype (subtype_v2, shift 0.05)        -> main_A_type1error.pdf
#   main figure, panel B   power vs gene-gene correlation, Additive and Matrix
#                          Variate, zero mass 0% and 30%            -> main_B_power.pdf
#   supplementary          type 1 error, power and elapsed time with SpaceDecorr
#                          at spline-basis dimensions 10% / 200 / 1000 -> supp_pairwise_basis.pdf
# The main panels carry no legend or axis titles (added in the layout software).

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

library(tidyverse)
library(cowplot)
library(data.table)
library(patchwork)

fig_dir <- file.path(dir.out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Shared ----
ALPHA     <- 0.05
ALPHA_COL <- "#CB181D"
FONT      <- "Arial"
RANGE     <- 0.25
NUM_CELLS_GRID <- c(500, 1000, 2000, 5000)

METHODS <- tibble::tribble(
  ~model,        ~label,            ~col,
  "spacedecorr", "SpaceDecorr",     "#8C510A",
  "spMOCA",      "spMOCA",          "#FF7F00",
  "logpearson",  "log-normalized",  "#000000",
  "spark",       "VST-normalized",  "#7F7F7F",
  "meringue",    "MERINGUE",        "#08519C",
  "giotto",      "Giotto",          "#6BAED6",
  "GraphR",      "GraphR",          "#1B7837")
method_cols <- setNames(METHODS$col, METHODS$label)
draw_levels <- c(setdiff(METHODS$label, "SpaceDecorr"), "SpaceDecorr")   # SpaceDecorr drawn on top

rows <- fread(file.path(dir.data, "pairwise_results.csv"),
              select = c("seed", "sim", "model", "num_cells", "range", "pi_zero", "rho",
                         "subtype_shift", "k", "pvalue", "time", "status"), showProgress = FALSE) %>%
  as_tibble() %>%
  filter(range == RANGE, model %in% METHODS$model, status == "ok") %>%
  mutate(model = factor(METHODS$label[match(model, METHODS$model)], levels = draw_levels))

# Text is measured on the cairo device the PDFs are written with.
with_cairo <- function(expr) {
  scratch <- tempfile(fileext = ".pdf")
  cairo_pdf(scratch, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
  on.exit({ dev.off(); unlink(scratch) })
  force(expr)
}
set_null_device(function(width, height) cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height))

# Main figure, panels A and B ----
# Every method at SpaceDecorr's default basis dimension (k = 200).
main <- rows %>%
  filter(k == 200 | is.na(k)) %>%
  group_by(sim, model, pi_zero, subtype_shift, num_cells, rho) %>%
  summarise(power = mean(pvalue < ALPHA), .groups = "drop")

PANEL_W   <- 3.33                            # in: plot area of one subplot
PANEL_H   <- 2.1
PANEL_GAP <- 0.4                             # cm between subplots
Y_LIM     <- c(-0.07, 1)                     # room below 0 for the alpha label

pane_theme <- theme_bw(base_size = 20, base_family = FONT) +
  theme(legend.position  = "none",
        axis.title       = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(PANEL_GAP, "cm"))

pane_size <- function(n, w = PANEL_W, h = PANEL_H)
  theme(panel.widths  = unit(n * w + (n - 1) * PANEL_GAP / 2.54, "in"),
        panel.heights = unit(h, "in"))

check_pane <- function(d, facet, pane) {
  n <- count(d, .data[[facet]], model, .drop = FALSE)
  if (any(n$n == 0))
    stop(pane, ": no result for ", paste(sprintf("%s in '%s'", n$model[n$n == 0], n[[1]][n$n == 0]), collapse = "; "))
  d
}

line_pane <- function(d, xvar, facet, x_breaks, href = NULL) {
  ref_line <- if (is.null(href)) NULL else
    geom_hline(yintercept = href, colour = ALPHA_COL, linetype = "dashed", linewidth = 0.4)
  ggplot(d, aes(.data[[xvar]], power, colour = model, group = model)) +
    ref_line +
    geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
    facet_wrap(vars(.data[[facet]]), nrow = 1) +
    scale_colour_manual(values = method_cols, breaks = METHODS$label) +
    scale_x_continuous(breaks = x_breaks) +
    coord_cartesian(ylim = Y_LIM) +
    pane_theme + pane_size(n_distinct(d[[facet]]))
}

save_pane <- function(plot, name) {
  gt   <- with_cairo(ggplotGrob(plot))
  size <- with_cairo(c(w = grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE),
                       h = grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)))
  ggsave(file.path(fig_dir, name), gt, width = size[["w"]], height = size[["h"]],
         device = cairo_pdf, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
  message("wrote ", name)
}

## panel A: type 1 error (rho = 0, no zero mass) vs number of cells
type1 <- filter(main, rho == 0, pi_zero == 0, num_cells %in% NUM_CELLS_GRID)
d_A <- bind_rows(
  "Spatial Covariance"   = filter(type1, sim == "Additive"),
  "FOV Technical Factor" = filter(type1, sim == "FOV_diff"),
  "Unobserved Subtype"   = filter(type1, sim == "subtype_v2", subtype_shift == 0.05),
  .id = "scenario") %>%
  mutate(scenario = factor(scenario, levels = c("Spatial Covariance", "FOV Technical Factor", "Unobserved Subtype")))

panel_A <- check_pane(d_A, "scenario", "panel A") %>%
  line_pane("num_cells", "scenario", x_breaks = c(500, 2000, 5000), href = ALPHA) +
  annotate("text", x = max(NUM_CELLS_GRID), y = 0, label = sprintf("alpha == %s", ALPHA), parse = TRUE,
           family = FONT, hjust = 1, vjust = 1.1, size = 12 / .pt, colour = ALPHA_COL)
save_pane(panel_A, "main_A_type1error.pdf")

## panel B: power vs rho at 2000 cells, Additive then Matrix Variate, zero mass 0% | 30%
power_pane <- function(sim_val) {
  main %>%
    filter(sim == sim_val, num_cells == 2000, rho > 0) %>%
    mutate(zero = factor(paste0("Zero mass: ", 100 * pi_zero, "%"), levels = c("Zero mass: 0%", "Zero mass: 30%"))) %>%
    check_pane("zero", paste("panel B", sim_val)) %>%
    line_pane("rho", "zero", x_breaks = seq(0, 0.9, 0.4))
}
pane_add <- power_pane("Additive") + theme(plot.margin = margin_part(r = 0))
pane_mv  <- power_pane("MV") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.margin = margin_part(l = 1, unit = "cm"))
panel_B <- cbind(with_cairo(ggplotGrob(pane_add)), with_cairo(ggplotGrob(pane_mv)), size = "first")
size_B  <- with_cairo(c(w = grid::convertWidth(sum(panel_B$widths), "in", valueOnly = TRUE),
                        h = grid::convertHeight(sum(panel_B$heights), "in", valueOnly = TRUE)))
ggsave(file.path(fig_dir, "main_B_power.pdf"), panel_B, width = size_B[["w"]], height = size_B[["h"]],
       device = cairo_pdf, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
message("wrote main_B_power.pdf")

# Supplementary: SpaceDecorr's spline-basis dimension ----
# Additive and MV without zero mass. A SpaceDecorr row is drawn under every rule
# its k meets (at 2000 cells 10% of the cells is 200); the other methods are one line each.
K_RULES     <- c("10%" = 0.1, "200" = 200, "1000" = 1000)
K_LINETYPES <- c("10%" = "solid", "200" = "dashed", "1000" = "dotdash")
rule_k <- function(rule, n) if (rule < 1) round(rule * n) else rule

supp_rows <- filter(rows, sim %in% c("Additive", "MV"), pi_zero == 0, num_cells %in% NUM_CELLS_GRID)
sd_rows   <- filter(supp_rows, model == "SpaceDecorr")
by_rule   <- lapply(names(K_RULES), function(nm)
  sd_rows %>% filter(k == rule_k(K_RULES[[nm]], num_cells)) %>% mutate(basis = nm))
per_seed <- bind_rows(c(by_rule, list(mutate(filter(supp_rows, model != "SpaceDecorr"), basis = "other")))) %>%
  mutate(basis = factor(basis, levels = c("other", names(K_RULES))),
         sim   = factor(sim, levels = c("Additive", "MV")),
         line  = factor(paste(model, basis),
                        levels = c(paste(setdiff(draw_levels, "SpaceDecorr"), "other"),
                                   paste("SpaceDecorr", rev(names(K_RULES))))))

supp <- per_seed %>%
  group_by(sim, model, basis, line, num_cells, rho) %>%
  summarise(power = mean(pvalue < ALPHA), mean_time = mean(time), .groups = "drop")

d_a <- filter(supp, rho == 0)
d_b <- filter(supp, num_cells == 2000, rho > 0)
d_c <- filter(supp, sim == "Additive", rho == 0)

BASE_SIZE <- 12
FIG_W     <- 15
BODY_H    <- 4.6
theme_fig <- theme_bw(base_size = BASE_SIZE, base_family = FONT) +
  theme(legend.position  = "none",
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text       = element_text(size = BASE_SIZE),
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(PANEL_GAP, "cm"),
        plot.tag         = element_text(face = "bold", size = BASE_SIZE + 2))
scale_methods <- function() scale_colour_manual(values = method_cols, breaks = METHODS$label, drop = FALSE)
scale_basis   <- function() scale_linetype_manual(values = c(K_LINETYPES, other = "solid"), breaks = names(K_RULES), drop = FALSE)
lines_layer   <- function() list(geom_line(aes(linetype = basis), linewidth = 0.7), geom_point(size = 0.9))

p_a <- ggplot(d_a, aes(num_cells, power, colour = model, group = line)) +
  geom_hline(yintercept = ALPHA, colour = ALPHA_COL, linetype = "dashed", linewidth = 0.4) +
  lines_layer() +
  annotate("text", x = max(NUM_CELLS_GRID), y = 0, label = sprintf("alpha == %s", ALPHA), parse = TRUE,
           family = FONT, hjust = 1, vjust = 1.1, size = (BASE_SIZE - 1) / .pt, colour = ALPHA_COL) +
  facet_wrap(vars(sim), nrow = 1) +
  scale_methods() + scale_basis() +
  scale_x_continuous(breaks = c(500, 2000, 5000)) +
  coord_cartesian(ylim = c(-0.06, 1)) +
  labs(x = "Number of cells", y = "Type 1 Error", colour = "Methods", linetype = "Basis dimension (K)", tag = "A") +
  theme_fig

p_b <- ggplot(d_b, aes(rho, power, colour = model, group = line)) +
  lines_layer() +
  facet_wrap(vars(sim), nrow = 1) +
  scale_methods() + scale_basis() +
  scale_x_continuous(breaks = seq(0, 0.8, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = expression("Correlation  " * rho), y = "Power", tag = "B") +
  theme_fig

p_c <- ggplot(d_c, aes(num_cells, mean_time, colour = model, group = line)) +
  lines_layer() +
  scale_methods() + scale_basis() +
  scale_x_continuous(breaks = c(500, 2000, 5000)) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1, 1, 10, 60, 600),
                labels = c("0.001s", "0.01s", "0.1s", "1s", "10s", "1m", "10m")) +
  labs(x = "Number of cells", y = "Elapsed time", tag = "C") +
  theme_fig

legend_src <- p_a +
  guides(colour = guide_legend(nrow = 1, order = 1), linetype = guide_legend(nrow = 1, order = 2)) +
  theme(legend.position = "top", legend.box = "horizontal",
        legend.title = element_text(size = BASE_SIZE - 1, margin = margin(r = 6)),
        legend.text  = element_text(size = BASE_SIZE - 1),
        legend.key.width = unit(1.8, "line"), legend.key.spacing.x = unit((BASE_SIZE - 1) / 2, "pt"),
        legend.spacing.x = unit(1.5, "line"), legend.box.spacing = unit(0, "pt"))
legend   <- get_plot_component(legend_src, "guide-box-top", return_all = TRUE)
legend_h <- round(grid::convertHeight(sum(legend$heights), "in", valueOnly = TRUE) + 0.1, 2)

body <- p_a + p_b + p_c + plot_layout(nrow = 1, widths = c(2, 2, 1))
fig  <- plot_grid(legend, body, ncol = 1, rel_heights = c(legend_h, BODY_H))
ggsave(file.path(fig_dir, "supp_pairwise_basis.pdf"), fig, width = FIG_W, height = legend_h + BODY_H,
       device = cairo_pdf, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
message("wrote supp_pairwise_basis.pdf")
