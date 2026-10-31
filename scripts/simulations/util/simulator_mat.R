# simulator_mat.R - the p x p benchmarks: sim_block (block co-expression
# structure), sim_niche_mat (differential network between two real niches) and
# sim_samples (differential network between two simulated tissues). Counts are
# drawn as in simulator.R from per-gene Gamma parameters sampled from fixed menus.
source(file.path(dir.proj, "util", "simulator.R"))
suppressPackageStartupMessages(library(igraph))

MAT_ZIP  <- file.path(dir.data, "NearPDmatrix.zip")       # the stored block correlation matrices
MAT_LUNG <- file.path(dir.data, "Lung5-3data.Rdata")      # the real tissue (`fulldata`)

MAT_MENU <- list(
  m       = list(values = c(0.15, 0.25, 0.5),      prob = c(1/3, 2/3, 1/3)),
  m_niche = list(values = c(0.15, 0.25, 0.5, 1),   prob = c(1/3, 2/3, 2/3, 1/3)),
  ratio   = list(values = c(0.1, 0.25, 0.5, 0.75), prob = c(1/4, 2/4, 2/4, 1/4)),
  range_niche = c(0.5, 1, 1.5, 2))

.draw_menu <- function(menu, num_genes)
  sample(menu$values, num_genes, replace = TRUE, prob = menu$prob)

.expand_gene <- function(vec, n) matrix(vec, nrow = n, ncol = length(vec), byrow = TRUE)

.gene_names <- function(p) paste0("Y", seq_len(p))

# Gene networks ----

# Small-world differential network. delta0: Watts-Strogatz ring, edge weights
# U(0.5, 0.7), shifted to a correlation matrix. delta1: a `perturb_prob` fraction
# of the pairs flipped (an edge removed, a non-edge added with weight U(0.5, 0.9)),
# the added weights shrunk by the largest common factor that keeps the matrix
# positive definite, so the unperturbed entries stay identical to delta0's.
.smallworld_deltas <- function(rept, num_genes, perturb_prob, nei = NULL, prob = NULL,
                               min_eig = 1e-3, rem_margin = 0.05) {
  if (is.null(nei))
    nei  <- switch(as.character(num_genes), "10" = 2, "20" = 4, "50" = 6,
                   stop("nei: no default for num_genes = ", num_genes))
  if (is.null(prob))
    prob <- switch(as.character(num_genes), "10" = 0.1, "20" = 0.1, "50" = 0.05,
                   stop("prob: no default for num_genes = ", num_genes))
  genes <- .gene_names(num_genes)
  min_eigen <- function(M) min(eigen(M, symmetric = TRUE, only.values = TRUE)$values)

  set.seed(68545 * rept)
  g_1 <- sample_smallworld(dim = 1, size = num_genes, nei = nei, p = prob)
  delta_1 <- as.matrix(as_adjacency_matrix(g_1))
  delta_1[delta_1 == 1 & upper.tri(delta_1)] <- runif(sum(delta_1 == 1 & upper.tri(delta_1)), min = 0.5, max = 0.7)
  delta_1[lower.tri(delta_1)] <- t(delta_1)[lower.tri(delta_1)]
  base_raw <- delta_1
  s0       <- -(min(as.numeric(eigen(delta_1)$values)) - 0.01)
  delta_1  <- cov2cor(delta_1 + diag(s0, num_genes))

  set.seed(5331 * rept)
  delta_2 <- delta_1
  edges2 <- which(upper.tri(delta_2), arr.ind = TRUE)
  perturb_idx2 <- sample(1:nrow(edges2), size = floor(perturb_prob * nrow(edges2)))
  perturbed <- added <- removed <- matrix(FALSE, num_genes, num_genes)
  for (idx in perturb_idx2) {
    i <- edges2[idx, 1]
    j <- edges2[idx, 2]
    perturbed[i, j] <- perturbed[j, i] <- TRUE
    if (delta_2[i, j] != 0) {
      delta_2[i, j] <- delta_2[j, i] <- 0
      removed[i, j] <- removed[j, i] <- TRUE
    } else {
      delta_2[i, j] <- delta_2[j, i] <- runif(1, 0.5, 0.9)
      added[i, j] <- added[j, i] <- TRUE
    }
  }
  u_added <- delta_2

  M_rem <- base_raw; M_rem[removed] <- 0
  s_rem <- -min_eigen(M_rem) / (1 - rem_margin)
  s     <- max(s0, s_rem)
  base_scale <- s0 / s
  delta_1 <- cov2cor(base_raw + diag(s, num_genes))
  at_c <- function(cc) { M <- delta_1; M[removed] <- 0; M[added] <- cc * u_added[added]; M }
  stopifnot(min_eigen(at_c(0)) >= min_eig - 1e-9)
  c_shrink <- 1
  if (min_eigen(at_c(1)) < min_eig) {
    lo <- 0; hi <- 1
    for (it in 1:60) { mid <- (lo + hi) / 2; if (min_eigen(at_c(mid)) >= min_eig) lo <- mid else hi <- mid }
    c_shrink <- lo
  }
  delta_2 <- at_c(c_shrink)

  dimnames(delta_1) <- dimnames(delta_2) <- dimnames(perturbed) <- list(genes, genes)
  offd <- !diag(num_genes)
  list(delta0 = delta_1, delta1 = delta_2,
       true_edge0 = delta_1 != 0 & offd, true_edge1 = delta_2 != 0 & offd, true_diff = perturbed,
       rho_true0 = delta_1, rho_true1 = delta_2, labels = NULL, genes = genes,
       nei = nei, prob = prob, c_shrink = c_shrink, base_scale = base_scale)
}

# The stored block matrix of one seed, read from the zip. The generating
# correlation is cov2cor() of it; the block labels (for the ARI) come from
# hierarchical clustering of the raw matrix.
.block_delta <- function(rept, num_genes, numblocks) {
  f <- sprintf("NearPDmatrix/matPD_rept%d_ngenes%d_numblocks%d.txt", rept, num_genes, numblocks)
  delta_raw <- as.matrix(read.delim(unz(MAT_ZIP, f), sep = " ", header = FALSE))
  stopifnot(nrow(delta_raw) == num_genes, ncol(delta_raw) == num_genes)
  genes <- .gene_names(num_genes)
  dimnames(delta_raw) <- list(genes, genes)
  labels <- cutree(hclust(dist(delta_raw)), k = numblocks)
  delta  <- cov2cor(delta_raw)
  offd   <- !diag(num_genes)
  list(delta0 = delta, delta1 = NULL, delta_raw = delta_raw,
       true_edge0 = delta_raw != 0 & offd, true_edge1 = NULL, true_diff = NULL,
       rho_true0 = delta, rho_true1 = NULL, labels = unname(labels), genes = genes)
}

# Real data ----
.mat_cache <- new.env(parent = emptyenv())

# Lung5-3: the library sizes of all cells (`totalcounts`) and the macrophages of
# the "myeloid-enriched stroma" (niche 1) and "stroma" (niche 0) niches, ordered
# niche 0 then niche 1, native coordinates.
mat_lung <- function(path = MAT_LUNG) {
  if (is.null(.mat_cache$lung)) {
    e <- new.env()
    load(path, envir = e)
    fd   <- e$fulldata
    keep <- which(fd$cell_type == "macrophage" & fd$niche %in% c("myeloid-enriched stroma", "stroma"))
    niche <- as.integer(fd$niche[keep] == "myeloid-enriched stroma")
    ord   <- order(niche)
    .mat_cache$lung <- list(
      totalcounts = fd$totalcounts,
      coordinates = data.frame(sdimx = fd$sdimx[keep][ord], sdimy = fd$sdimy[keep][ord]),
      niche       = niche[ord],
      Area        = fd$Area[keep][ord])
    rm(fd, e); invisible(gc())
  }
  .mat_cache$lung
}

# Cached factor of the Matern covariance over all niche cells at one range.
mat_field_factor <- function(field_range) {
  key <- paste0("factor|", field_range)
  if (is.null(.mat_cache[[key]])) {
    g <- mat_lung()
    message(sprintf("niche_mat: factoring the Matern covariance of %d cells (range %s)", nrow(g$coordinates), field_range))
    sigma <- fields::Matern(as.matrix(dist(g$coordinates)), smoothness = 0.5, range = field_range)
    .mat_cache[[key]] <- tryCatch(t(chol(sigma)), error = function(e) .field_factor(sigma))
  }
  .mat_cache[[key]]
}

# Fill the cache before forking so the workers share it.
warm_mat_cache <- function(niche_fields = FALSE) {
  mat_lung()
  if (niche_fields) for (r in MAT_MENU$range_niche) mat_field_factor(r)
  invisible(NULL)
}

.assemble <- function(Y, coordinates, libsize) {
  Y <- as.data.frame(Y)
  colnames(Y) <- .gene_names(ncol(Y))
  rownames(Y) <- seq_len(nrow(Y))
  data <- Y %>% cbind(coordinates) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize
  rownames(data) <- rownames(Y)
  list(Y = Y, data = data)
}

# Block benchmark ----
# One tissue of num_cells cells on [-1, 1]^2; the genes' correlation matrix is the
# stored block matrix of this seed. Library sizes are sampled from the real cells'.
# Returns list(Y, coordinates, data, membership = NULL, meta, ext_info, truth).
sim_block <- function(rept, sim = c("block_MV", "block_Additive"), num_cells = 2000,
                      num_genes = 100, range = 0.25, numblocks = 4, pi_zero = 0, smoothness = 0.5) {
  sim  <- match.arg(sim)
  base <- if (sim == "block_MV") "MV" else "Additive"
  lung <- mat_lung()

  set.seed(5675 * rept)
  m     <- .draw_menu(MAT_MENU$m, num_genes)
  ratio <- .draw_menu(MAT_MENU$ratio, num_genes)
  v     <- m / ratio
  libsize <- sample(lung$totalcounts, num_cells)
  coordinates <- data.frame(sdimx = runif(num_cells, -1, 1), sdimy = runif(num_cells, -1, 1))

  sigma <- fields::Matern(as.matrix(dist(coordinates)), smoothness = smoothness, range = range)
  truth <- .block_delta(rept, num_genes, numblocks)

  Y <- .sim_counts(base, sigma, truth$delta0, .expand_gene(m, num_cells), .expand_gene(v, num_cells),
                   libsize, pi_zero = pi_zero, num_genes = num_genes)
  stopifnot(all(colSums(Y) > 0))
  a <- .assemble(Y, coordinates, libsize)

  meta <- data.frame(seed = rept, sim = sim, num_cells = num_cells, num_genes = num_genes,
                     range = range, numblocks = numblocks, pi_zero = pi_zero, stringsAsFactors = FALSE)
  list(Y = a$Y, coordinates = coordinates, data = a$data, membership = NULL,
       meta = meta, ext_info = list(), truth = truth)
}

# Differential network, two real niches ----
# Every niche macrophage at its real position; niche 0 (stroma) carries delta0,
# niche 1 (myeloid-enriched stroma) delta1. Additive model, each gene with its own
# spatial range from the menu; library size = round(Area * libsize_scale).
sim_niche_mat <- function(rept, num_genes = 20, perturb_prob = 0.25, nei = NULL, prob = NULL,
                          libsize_scale = 0.1, pi_zero = 0) {
  lung        <- mat_lung()
  coordinates <- lung$coordinates
  niche       <- lung$niche
  num_cells   <- nrow(coordinates)
  n0 <- sum(niche == 0); n1 <- sum(niche == 1)

  set.seed(5675 * rept); m     <- .draw_menu(MAT_MENU$m_niche, num_genes)
  set.seed(6786 * rept); ratio <- .draw_menu(MAT_MENU$ratio, num_genes)
  v       <- m / ratio
  libsize <- round(lung$Area * libsize_scale)
  range_gene <- sample(MAT_MENU$range_niche, replace = TRUE, num_genes)

  truth <- .smallworld_deltas(rept, num_genes, perturb_prob, nei, prob)

  set.seed(47289 * rept)
  spatial <- vapply(seq_len(num_genes),
                    function(x) as.numeric(.draw_field(mat_field_factor(range_gene[x]), 1)),
                    numeric(num_cells))
  Y <- .sim_counts("Additive", sigma = NULL, delta = list(truth$delta0, truth$delta1),
                   .expand_gene(m, num_cells), .expand_gene(v, num_cells),
                   libsize, pi_zero = pi_zero, num_genes = num_genes, spatial = spatial, group = niche)
  stopifnot(all(colSums(Y) > 0))
  a <- .assemble(Y, coordinates, libsize)

  meta <- data.frame(seed = rept, sim = "niche_mat", num_cells = num_cells, num_genes = num_genes,
                     perturb_prob = perturb_prob, nei = truth$nei, prob = truth$prob,
                     libsize_scale = libsize_scale, n0 = n0, n1 = n1, stringsAsFactors = FALSE)
  list(Y = a$Y, coordinates = coordinates, data = a$data, membership = niche, meta = meta,
       ext_info = list(range_gene = range_gene, c_shrink = truth$c_shrink, base_scale = truth$base_scale),
       truth = truth)
}

# Differential network, two simulated tissues ----
# Two tissues on [-1, 1]^2 (fixed coordinates), tissue 1 with delta0, tissue 2
# with delta1; library sizes sampled from the real cells'.
# Returns list(sample1, sample2, membership = NULL, meta, ext_info, truth).
sim_samples <- function(rept, sim = c("sample_MV", "sample_Additive"), num_genes = 20,
                        num_cells1 = 2000, num_cells2 = 2000, range = 0.25,
                        perturb_prob = 0.25, nei = NULL, prob = NULL, pi_zero = 0, smoothness = 0.5) {
  sim  <- match.arg(sim)
  base <- if (sim == "sample_MV") "MV" else "Additive"
  lung <- mat_lung()

  set.seed(378); sample_coords1 <- data.frame(sdimx = runif(num_cells1, -1, 1), sdimy = runif(num_cells1, -1, 1))
  set.seed(354); sample_coords2 <- data.frame(sdimx = runif(num_cells2, -1, 1), sdimy = runif(num_cells2, -1, 1))
  set.seed(5675 * rept); m     <- .draw_menu(MAT_MENU$m, num_genes)
  set.seed(6786 * rept); ratio <- .draw_menu(MAT_MENU$ratio, num_genes)
  v <- m / ratio
  set.seed(6786 * rept); libsize1 <- sample(lung$totalcounts, num_cells1)
  set.seed(512 * rept);  libsize2 <- sample(lung$totalcounts, num_cells2)

  truth <- .smallworld_deltas(rept, num_genes, perturb_prob, nei, prob)

  set.seed(47289 * rept)
  draw <- function(coords, delta, libsize) {
    n <- nrow(coords)
    sigma <- fields::Matern(as.matrix(dist(coords)), smoothness = smoothness, range = range)
    Y <- .sim_counts(base, sigma, delta, .expand_gene(m, n), .expand_gene(v, n),
                     libsize, pi_zero = pi_zero, num_genes = num_genes)
    stopifnot(all(colSums(Y) > 0))
    a <- .assemble(Y, coords, libsize)
    list(Y = a$Y, coordinates = coords, data = a$data)
  }
  sample1 <- draw(sample_coords1, truth$delta0, libsize1)
  sample2 <- draw(sample_coords2, truth$delta1, libsize2)

  meta <- data.frame(seed = rept, sim = sim, num_genes = num_genes,
                     num_cells1 = num_cells1, num_cells2 = num_cells2, range = range,
                     perturb_prob = perturb_prob, nei = truth$nei, prob = truth$prob, stringsAsFactors = FALSE)
  list(sample1 = sample1, sample2 = sample2, membership = NULL, meta = meta,
       ext_info = list(c_shrink = truth$c_shrink, base_scale = truth$base_scale), truth = truth)
}
