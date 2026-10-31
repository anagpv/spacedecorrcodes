# figure_niche.R - Scenario II figures from niche_results.csv (run_niche.R):
#   main figure, panel C   map of the two real niches            -> main_C_niche_map.pdf
#                          power vs within-niche correlation on the real tissue
#                          geometry: Tissue / Niche / Differential -> main_C_niche_power.pdf
#   supplementary          the simulated disc niche: A rejection rate vs niche radius,
#                          B expression maps of the four niche settings, C type 1
#                          error of the inert pair per setting and method
#                                                               -> supp_niche_sensitivity.pdf
#   supplementary          type 1 error of the inert pair on the real tissue geometry
#                          vs within-niche correlation            -> supp_scenario2_type1error.pdf
# The main panels carry no legend or axis titles (added in the layout software).

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

library(tidyverse)
library(cowplot)
library(data.table)
library(patchwork)

suppressPackageStartupMessages(source(file.path(dir.proj, "util", "simulator_mat.R")))   # sim_niche(), mat_lung()

fig_dir <- file.path(dir.out, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Shared ----
ALPHA     <- 0.05
ALPHA_COL <- "#CB181D"
FONT      <- "Arial"

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

TESTS <- c("pvalue", "pvalue_g1", "pvalue_diff", "pvalue_null", "pvalue_g1_null")
niche <- fread(file.path(dir.data, "niche_results.csv"),
               select = c("seed", "sim", "model", "range", "pi_zero", "rho", "rho_niche", "radius",
                          "range_niche", "fov_hi", TESTS), showProgress = FALSE) %>%
  as_tibble() %>%
  filter(model %in% METHODS$model, pi_zero == 0, rho == 0) %>%
  mutate(model = factor(METHODS$label[match(model, METHODS$model)], levels = draw_levels))
stopifnot(!any(duplicated(niche[, c("seed", "sim", "model", "range", "rho_niche", "radius", "range_niche", "fov_hi")])))

# rejection rate at ALPHA; NA when fewer than half of the seeds have a p-value
rate <- function(p) if (mean(!is.na(p)) < 0.5) NA_real_ else mean(p < ALPHA, na.rm = TRUE)

check_complete <- function(d, facets, panel, yvar = "rate") {
  n <- d %>% filter(is.finite(.data[[yvar]])) %>% count(across(all_of(facets)), model, .drop = FALSE)
  if (any(n$n == 0)) {
    gap <- n[n$n == 0, ]
    stop(panel, ": no result for ",
         paste(sprintf("%s in '%s'", gap$model, apply(gap[, facets, drop = FALSE], 1, paste, collapse = " / ")), collapse = "; "))
  }
  invisible(d)
}

with_cairo <- function(expr, width = 14, height = 14) {
  scratch <- tempfile(fileext = ".pdf")
  cairo_pdf(scratch, width = width, height = height, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
  on.exit({ dev.off(); unlink(scratch) })
  force(expr)
}
set_null_device(function(width, height) cairo_pdf(tempfile(fileext = ".pdf"), width = width, height = height))
grob_in <- function(g, what = c("widths", "heights")) {
  what <- match.arg(what)
  conv <- if (what == "widths") grid::convertWidth else grid::convertHeight
  with_cairo(conv(sum(g[[what]]), "in", valueOnly = TRUE))
}

# the real-geometry slice of the main panel and the scenario II error figure:
# range 0.25, no FOV effect, one shared spatial field
real <- filter(niche, sim == "niche_real", range == 0.25, fov_hi == 0, is.na(range_niche))
stopifnot(all(is.na(real$radius)))

# Main figure, panel C ----
PANEL_W   <- 3.33                            # in: plot area of one subplot
PANEL_H   <- 4.66
PANEL_GAP <- 0.4                             # cm between subplots
MAP_RATIO <- 1.4                             # map box height : width
Y_LIM     <- c(-0.07, 1)

pane_theme <- theme_bw(base_size = 20, base_family = FONT) +
  theme(legend.position  = "none",
        axis.title       = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        panel.grid.minor = element_blank(),
        panel.spacing    = unit(PANEL_GAP, "cm"))

save_pane <- function(plot, name) {
  gt   <- with_cairo(ggplotGrob(plot))
  size <- with_cairo(c(w = grid::convertWidth(sum(gt$widths), "in", valueOnly = TRUE),
                       h = grid::convertHeight(sum(gt$heights), "in", valueOnly = TRUE)))
  ggsave(file.path(fig_dir, name), gt, width = size[["w"]], height = size[["h"]],
         device = cairo_pdf, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
  message("wrote ", name)
}

## the map: the Lung5-3 macrophages, myeloid-enriched stroma (grey) and stroma (blue)
lung <- mat_lung()
niche_cols <- c("Myeloid-enriched stroma" = "grey60", "Stroma" = "#21A1DC")
d_map <- tibble(lung$coordinates, niche = factor(lung$niche, levels = c(1, 0), labels = names(niche_cols)))
map <- ggplot(d_map, aes(sdimx, sdimy, colour = niche)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_colour_manual(values = niche_cols) +
  labs(x = NULL, y = NULL) +
  pane_theme +
  theme(axis.text = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
        panel.grid = element_blank(),
        panel.widths = unit(PANEL_H / MAP_RATIO, "in"), panel.heights = unit(PANEL_H, "in"))
save_pane(map, "main_C_niche_map.pdf")

## power vs rho_niche: the true pair over the tissue, within the niche, niche vs rest
d_C <- real %>%
  group_by(model, rho_niche) %>%
  summarise(Tissue = rate(pvalue), Niche = rate(pvalue_g1), Differential = rate(pvalue_diff), .groups = "drop") %>%
  pivot_longer(c(Tissue, Niche, Differential), names_to = "level", values_to = "rate") %>%
  mutate(level = factor(level, levels = c("Tissue", "Niche", "Differential"))) %>%
  check_complete("level", "panel C")

power <- ggplot(d_C, aes(rho_niche, rate, colour = model, group = model)) +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  facet_wrap(vars(level), nrow = 1) +
  scale_colour_manual(values = method_cols, breaks = METHODS$label) +
  scale_x_continuous(breaks = seq(0, 0.9, 0.4)) +
  coord_cartesian(ylim = Y_LIM) +
  pane_theme +
  theme(axis.text = element_text(size = 20),
        panel.widths  = unit(3 * PANEL_W + 2 * PANEL_GAP / 2.54, "in"),
        panel.heights = unit(PANEL_H, "in"))
save_pane(power, "main_C_niche_power.pdf")

# Supplementary: niche sensitivity (simulated disc niche) ----
METHODS6 <- filter(METHODS, model != "GraphR")
disc <- filter(niche, sim == "niche", model %in% METHODS6$label) %>% mutate(model = droplevels(model))

A_RHO_NICHE <- 0.5
A_RANGES    <- c(0.25, 0.8)
RADIUS_GRID <- c(0.2, 0.4, 0.6, 0.8, 1)
C_RANGE     <- 0.25
C_RADIUS    <- 0.6
C_RHO_NICHE <- 0.5
ILLUSTRATION_SEED <- 1

# the four niche settings: rows of B and C
SETTINGS <- tibble::tribble(
  ~key,  ~fov_hi, ~range_niche, ~label,
  "s1",  0,       NA_real_,     "Shared field\nno FOV effect",
  "s2",  1,       NA_real_,     "Shared field\nFOV effect",
  "s3",  0,       0.1,          "Split field\nno FOV effect",
  "s4",  1,       0.1,          "Split field\nFOV effect")
setting_f <- function(key) factor(SETTINGS$label[match(key, SETTINGS$key)], levels = SETTINGS$label)
setting_of <- function(fov_hi, range_niche) {
  id <- function(f, r) paste(f, ifelse(is.na(r), "NA", r))
  SETTINGS$key[match(id(fov_hi, range_niche), id(SETTINGS$fov_hi, SETTINGS$range_niche))]
}

A_TESTS <- c("Negative Control - Tissue" = "pvalue_null",
             "Truth - Tissue"            = "pvalue",
             "Truth - Niche"             = "pvalue_g1",
             "Truth - Differential"      = "pvalue_diff")
range_f <- function(range) factor(sprintf("Range: %g", range), levels = sprintf("Range: %g", A_RANGES))

EXPR_MID  <- 5
EXPR_LOW  <- c("#DECDC4", "#D8B2A2", "#DD9A79")
EXPR_HIGH <- c("#D27953", "#B55231", "#91301C", "#721A1B", "#5D1023", "#4F1530")
expr_scale <- function(expr) {
  hi <- max(expr)
  stopifnot(hi > EXPR_MID)
  at <- c(seq(0, EXPR_MID, length.out = length(EXPR_LOW)),
          seq(EXPR_MID, hi, length.out = length(EXPR_HIGH) + 1)[-1])
  scale_colour_gradientn(colours = c(EXPR_LOW, EXPR_HIGH), values = at / hi, limits = c(0, hi))
}
GUIDE_SHARED <- list(colour = "#65000B", linetype = "22")      # the niche circle under a shared field
GUIDE_SPLIT  <- list(colour = "black",   linetype = "solid")   # under a split field
GUIDE_FOV    <- list(colour = "black",   linetype = "solid")   # the FOV square

BASE_SIZE    <- 14
FIG_W        <- 14
A_PANEL_H    <- 2
A_GAP        <- 0.3
A_EXTRA_GAP  <- 0.45
B_SHARE      <- 0.5
MAP_GAP      <- 0.15
TILE_GAP     <- 0.06
BC_GAP       <- 0.4
POINT_PER_IN <- 1.45
HEADER_TILT  <- 45

theme_niche <- theme_bw(base_size = BASE_SIZE, base_family = FONT) +
  theme(panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA),
        strip.text       = element_text(size = BASE_SIZE),
        legend.position  = "none",
        plot.tag         = element_text(face = "bold", size = BASE_SIZE + 2))

## A: rejection rate vs niche radius, setting 1, two ranges
d_A <- disc %>%
  filter(fov_hi == 0, is.na(range_niche), rho_niche == A_RHO_NICHE, range %in% A_RANGES) %>%
  group_by(model, range, radius) %>%
  summarise(across(all_of(unname(A_TESTS)), rate), .groups = "drop") %>%
  pivot_longer(all_of(unname(A_TESTS)), names_to = "test", values_to = "rate") %>%
  mutate(test = factor(names(A_TESTS)[match(test, A_TESTS)], levels = names(A_TESTS)), range_f = range_f(range)) %>%
  check_complete(c("range_f", "test"), "panel A")

## C: type 1 error of the inert pair at range 0.25, radius 0.6, across the four settings
d_Cs <- disc %>%
  filter(range == C_RANGE, radius == C_RADIUS, rho_niche == C_RHO_NICHE) %>%
  mutate(setting = setting_f(setting_of(fov_hi, range_niche)))
stopifnot(!anyNA(d_Cs$setting))
d_Cs <- d_Cs %>%
  group_by(setting, model) %>%
  summarise(rate = rate(pvalue_null), .groups = "drop") %>%
  mutate(model = factor(model, levels = METHODS6$label)) %>%
  check_complete("setting", "panel C")

line_block <- function(d, href) {
  ggplot(d, aes(radius, rate, colour = model, group = model)) +
    geom_hline(yintercept = href, colour = "grey55", linetype = "dotted", linewidth = 0.4) +
    geom_line(linewidth = 0.7) + geom_point(size = 1.2) +
    facet_grid(rows = vars(range_f), cols = vars(test)) +
    scale_colour_manual(values = method_cols, breaks = METHODS6$label, drop = FALSE) +
    scale_x_continuous("Niche radius", breaks = RADIUS_GRID) +
    scale_y_continuous("Power", breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(colour = NULL) +
    theme_niche +
    theme(panel.spacing.x = unit(A_GAP, "in"), panel.spacing.y = unit(A_GAP, "in"))
}

panel_A <- function(d, w) {
  null  <- line_block(filter(d, test == names(A_TESTS)[1]), href = ALPHA) +
    labs(tag = "A") + theme(strip.text.y = element_blank())
  truth <- line_block(filter(d, test %in% names(A_TESTS)[2:4]), href = c(ALPHA, 0.95)) +
    theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          plot.margin = margin_part(l = A_EXTRA_GAP, unit = "in"))
  null + truth +
    plot_layout(widths = unit(c(w, 3 * w + 2 * A_GAP), "in"), heights = unit(2 * A_PANEL_H + A_GAP, "in"),
                axis_titles = "collect_x")
}

## B: one simulated data set per setting, genes 1-3, the niche circle and FOV square drawn on
sims <- lapply(seq_len(nrow(SETTINGS)), function(i) {
  s <- SETTINGS[i, ]
  sim_niche(rept = ILLUSTRATION_SEED, sim = "niche", num_cells = 2000, range = C_RANGE, rho = 0,
            rho_niche = C_RHO_NICHE, m = 0.5, ratio = 0.5, libsize = 200, radius = C_RADIUS,
            range_niche = if (is.na(s$range_niche)) NULL else s$range_niche,
            fov_effect = if (s$fov_hi == 1) c(0, 1) else c(0, 0), elevate_effect = c(0, 1), pi_zero = 0)
})

circle_path <- function(r, n = 256) {
  th <- seq(0, 2 * pi, length.out = n)
  tibble(x = r * cos(th), y = r * sin(th))
}
panel_B <- function(sims, side) {
  gene_f <- factor(paste("Gene", 1:3), levels = paste("Gene", 1:3))
  d <- bind_rows(lapply(seq_along(sims), function(i) {
    Y   <- as.matrix(sims[[i]]$Y[, 1:3])
    lib <- sims[[i]]$data$libsize
    co  <- sims[[i]]$coordinates
    tibble(setting = setting_f(SETTINGS$key[i]), gene = rep(gene_f, each = nrow(Y)),
           sdimx = rep(co$sdimx, 3), sdimy = rep(co$sdimy, 3),
           expr = as.vector(log1p(Y / lib * median(lib))))
  })) %>% arrange(expr)
  paths <- function(keys, path)
    expand_grid(setting = setting_f(keys), gene = gene_f) %>% mutate(path = list(path)) %>% unnest(path)
  guide <- function(d, style)
    geom_path(aes(x, y), d, inherit.aes = FALSE, colour = style$colour, linetype = style$linetype, linewidth = 1)
  circle <- circle_path(C_RADIUS)
  square <- tibble(x = c(0, 0.5, 0.5, 0, 0), y = c(-0.5, -0.5, 0, 0, -0.5))

  ggplot(d, aes(sdimx, sdimy, colour = expr)) +
    geom_point(size = POINT_PER_IN * side, alpha = 0.85, stroke = 0) +
    expr_scale(d$expr) +
    guide(paths(SETTINGS$key[is.na(SETTINGS$range_niche)],  circle), GUIDE_SHARED) +
    guide(paths(SETTINGS$key[!is.na(SETTINGS$range_niche)], circle), GUIDE_SPLIT) +
    guide(paths(SETTINGS$key[SETTINGS$fov_hi == 1],          square), GUIDE_FOV) +
    facet_grid(rows = vars(setting), cols = vars(gene), switch = "y") +
    coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    labs(x = NULL, y = NULL, tag = "B") +
    theme_niche +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = EXPR_LOW[1], colour = NA),
          axis.text = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          strip.placement = "outside",
          strip.text.y.left = element_text(size = 11, angle = 0, hjust = 1),
          panel.spacing = unit(MAP_GAP, "in"))
}

## C: one tile per (setting, method) with the rate x 100
text_in <- function(labels, size)
  with_cairo(max(grid::convertWidth(grid::grobWidth(grid::textGrob(labels, gp = grid::gpar(fontsize = size, fontfamily = FONT))),
                                    "in", valueOnly = TRUE)))
panel_C <- function(d, tile_w = Inf) {
  d <- mutate(d, label = sprintf("%d", as.integer(floor(100 * rate + 0.5))))
  tilt <- text_in(METHODS6$label, BASE_SIZE) + 0.1 > tile_w
  header <- if (tilt) element_text(size = BASE_SIZE, angle = HEADER_TILT, hjust = 0, vjust = 0, margin = margin(b = 4))
            else      element_text(size = BASE_SIZE, margin = margin(b = 4))
  ggplot(d, aes(model, 1)) +
    geom_tile(aes(fill = rate), width = 1, height = 1) +
    geom_text(aes(label = label), size = BASE_SIZE / .pt, family = FONT) +
    facet_grid(rows = vars(setting), cols = vars(model), scales = "free_x") +
    scale_fill_gradientn(colours = c("white", "grey70"), limits = c(0, 1)) +
    scale_x_discrete(position = "top", expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, tag = "C") +
    theme_niche +
    theme(panel.grid = element_blank(), panel.border = element_blank(),
          axis.text.y = element_blank(), axis.ticks = element_blank(), axis.ticks.length = unit(0, "pt"),
          axis.text.x.top = header, strip.text = element_blank(),
          panel.spacing.x = unit(TILE_GAP, "in"), panel.spacing.y = unit(MAP_GAP, "in"),
          plot.margin = margin_part(l = BC_GAP, unit = "in"))
}

b_area_w  <- function(side)   3 * side + 2 * MAP_GAP
c_area_w  <- function(tile_w) nrow(METHODS6) * tile_w + (nrow(METHODS6) - 1) * TILE_GAP
bc_area_h <- function(side)   4 * side + 3 * MAP_GAP
row_BC <- function(side, tile_w)
  panel_B(sims, side) + panel_C(d_Cs, tile_w) +
    plot_layout(widths = unit(c(b_area_w(side), c_area_w(tile_w)), "in"), heights = unit(bc_area_h(side), "in"))

# panel widths that fill FIG_W, from the width everything else takes
fit_width <- function(build, n_panels, gaps, target = FIG_W, nominal = 2) {
  other <- grob_in(with_cairo(patchworkGrob(build(nominal)))) - n_panels * nominal - gaps
  (target - other - gaps) / n_panels
}
a_w    <- fit_width(function(w) panel_A(d_A, w), n_panels = 4, gaps = 2 * A_GAP)
side   <- fit_width(function(s) wrap_plots(panel_B(sims, s), widths = unit(b_area_w(s), "in"), heights = unit(bc_area_h(s), "in")),
                    n_panels = 3, gaps = 2 * MAP_GAP, target = B_SHARE * FIG_W)
tile_w <- fit_width(function(w) wrap_plots(panel_C(d_Cs, w), widths = unit(c_area_w(w), "in"), heights = unit(bc_area_h(side), "in")),
                    n_panels = nrow(METHODS6), gaps = (nrow(METHODS6) - 1) * TILE_GAP, target = (1 - B_SHARE) * FIG_W)

legend_src <- line_block(filter(d_A, test %in% names(A_TESTS)[2:4]), href = ALPHA) +
  guides(colour = guide_legend(nrow = 1, title = "Methods")) +
  theme(legend.position = "top",
        legend.title = element_text(size = BASE_SIZE), legend.text = element_text(size = BASE_SIZE),
        legend.key.width = unit(2, "line"), legend.spacing.x = unit(1, "line"))
legend <- get_plot_component(legend_src, "guide-box-top", return_all = TRUE)

rows <- with_cairo(list(legend, patchworkGrob(panel_A(d_A, a_w)), patchworkGrob(row_BC(side, tile_w))))
hs   <- vapply(rows, grob_in, numeric(1), what = "heights")
fig  <- plot_grid(plotlist = rows, ncol = 1, rel_heights = hs)
ggsave(file.path(fig_dir, "supp_niche_sensitivity.pdf"), fig, width = FIG_W, height = sum(hs), device = cairo_pdf)
message("wrote supp_niche_sensitivity.pdf")

# Supplementary: scenario II error control (real tissue geometry) ----
S2_TESTS <- c("Negative Control - Tissue" = "pvalue_null",
              "Negative Control - Niche"  = "pvalue_g1_null")
d_S2 <- real %>%
  group_by(model, rho_niche) %>%
  summarise(across(all_of(unname(S2_TESTS)), rate), .groups = "drop") %>%
  pivot_longer(all_of(unname(S2_TESTS)), names_to = "test", values_to = "rate") %>%
  mutate(test = factor(names(S2_TESTS)[match(test, S2_TESTS)], levels = names(S2_TESTS))) %>%
  check_complete("test", "scenario II error")

fig_S2 <- ggplot(d_S2, aes(rho_niche, rate, colour = model, group = model)) +
  geom_hline(yintercept = ALPHA, colour = ALPHA_COL, linetype = "dashed", linewidth = 0.4) +
  geom_line(linewidth = 0.7) + geom_point(size = 0.9) +
  annotate("text", x = max(d_S2$rho_niche), y = 0, label = sprintf("alpha == %s", ALPHA), parse = TRUE,
           family = FONT, hjust = 1, vjust = 1.1, size = 5, colour = ALPHA_COL) +
  facet_wrap(vars(test), nrow = 1) +
  scale_colour_manual(values = method_cols, breaks = METHODS$label) +
  scale_x_continuous(breaks = seq(0, 0.9, 0.4)) +
  coord_cartesian(ylim = c(-0.04, 1)) +
  labs(x = expression("Correlation within the niche  " * rho[niche]), y = "Type 1 Error", colour = "Methods") +
  guides(colour = guide_legend(nrow = 2, byrow = TRUE)) +
  theme_bw(base_size = BASE_SIZE, base_family = FONT) +
  theme(legend.position = "top", legend.justification = "centre",
        legend.title = element_text(size = BASE_SIZE), legend.text = element_text(size = BASE_SIZE),
        legend.key.width = unit(2, "line"), legend.spacing.x = unit(1, "line"),
        strip.background = element_rect(fill = NA, colour = NA), strip.text = element_text(size = BASE_SIZE),
        panel.grid.minor = element_blank(), panel.spacing = unit(0.4, "cm"))
ggsave(file.path(fig_dir, "supp_scenario2_type1error.pdf"), fig_S2, width = 9, height = 5.4,
       device = cairo_pdf, symbolfamily = cairoSymbolFont(FONT, usePUA = FALSE))
message("wrote supp_scenario2_type1error.pdf")
