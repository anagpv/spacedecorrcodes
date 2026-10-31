# compare_mat.R - one wrapper per method for the p x p correlation matrix of all
# genes, on all cells and within each of two groups; mat_pair_rows() turns the
# matrices into one row per gene pair with the tests of util.R.
source(file.path(dir.proj, "util", "compare.R"))   # .log_cpm(), the method packages

# Helpers ----

.align_genes <- function(M, genes) {
  M <- as.matrix(M)
  stopifnot(!is.null(dimnames(M)), all(genes %in% rownames(M)), all(genes %in% colnames(M)))
  M <- M[genes, genes, drop = FALSE]
  stopifnot(identical(rownames(M), genes), identical(colnames(M), genes))
  M
}

# compute(idx) -> p x p matrix on the cells idx (NULL = all). `global = FALSE`
# skips the all-cells matrix (methods refit per group only).
.grouped_mat <- function(compute, membership, global = TRUE) {
  est <- if (isTRUE(global) || is.null(membership)) compute(NULL) else NULL
  if (is.null(membership))
    return(list(cor = est, cor_g0 = NULL, cor_g1 = NULL, n_g0 = NA_integer_, n_g1 = NA_integer_))
  i0 <- which(membership == 0)
  i1 <- which(membership == 1)
  stopifnot(length(i0) > 3, length(i1) > 3)
  list(cor = est, cor_g0 = compute(i0), cor_g1 = compute(i1), n_g0 = length(i0), n_g1 = length(i1))
}

.rows <- function(M, idx) if (is.null(idx)) M else M[idx, , drop = FALSE]

.finish <- function(g, n, tm, resid = NULL)
  c(g, list(n = n, time = unname(tm["elapsed"]), resid = resid))

# Methods ----
# Each takes (counts, coords, data, libsize, membership, ...) and returns
# list(cor, cor_g0, cor_g1, n, n_g0, n_g1, time, resid); all matrices p x p.

# log-normalized: Pearson on log-CPM
compare_mat_logpearson <- function(counts, coords, data, libsize = "libsize", membership = NULL, ...) {
  genes <- colnames(counts)
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    g <- .grouped_mat(function(idx) .align_genes(cor(.rows(logY, idx)), genes), membership)
  })
  .finish(g, nrow(logY), tm)
}

# VST-normalized: Pearson on sctransform residuals (one vst on all cells)
compare_mat_spark <- function(counts, coords, data, libsize = "libsize", membership = NULL, ...) {
  genes <- colnames(counts); Y <- as.matrix(counts)
  stopifnot(all(colSums(Y) > 0))
  tm <- system.time({
    Yvst <- t(sctransform::vst(umi = t(Y), cell_attr = data, latent_var = libsize,
                               min_cells = 1, verbosity = 0)$y)
    stopifnot(ncol(Yvst) == ncol(Y))
    g <- .grouped_mat(function(idx) .align_genes(cor(.rows(Yvst, idx)), genes), membership)
  })
  .finish(g, nrow(Y), tm)
}

# Giotto: Delaunay-network spatial correlation, refit per group
.giotto_cormat <- function(counts, logcounts, coords, minimum_k = 2) {
  genes <- colnames(counts)
  cor_DT <- suppressWarnings(suppressMessages({
    gobject <- Giotto::createGiottoObject(
      expression      = list(raw = t(as.matrix(counts)), normalized = t(logcounts)),
      expression_feat = "rna",
      spatial_locs    = as.data.frame(coords), verbose = FALSE)
    gobject <- Giotto::createSpatialNetwork(gobject = gobject, minimum_k = minimum_k, verbose = FALSE)
    Giotto::detectSpatialCorFeats(gobject, method = "network",
                                  spatial_network_name = "Delaunay_network")$cor_DT
  }))
  out <- matrix(NA_real_, length(genes), length(genes), dimnames = list(genes, genes))
  out[cbind(as.character(cor_DT$feat_ID), as.character(cor_DT$variable))] <- cor_DT$spat_cor
  out
}

compare_mat_giotto <- function(counts, coords, data, libsize = "libsize", membership = NULL,
                               minimum_k = 2, global = TRUE, ...) {
  genes <- colnames(counts)
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    g <- .grouped_mat(function(idx)
      .align_genes(.giotto_cormat(.rows(counts, idx), .rows(logY, idx), .rows(coords, idx), minimum_k), genes),
      membership, global)
  })
  .finish(g, nrow(counts), tm)
}

# MERINGUE: spatial cross-correlation, refit per group
compare_mat_meringue <- function(counts, coords, data, libsize = "libsize", membership = NULL,
                                 filterDist = 0.1, global = TRUE, ...) {
  genes <- colnames(counts)
  coords <- as.data.frame(coords); rownames(coords) <- rownames(counts)
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    g <- .grouped_mat(function(idx) {
      w <- MERINGUE::getSpatialNeighbors(.rows(coords, idx), filterDist = filterDist, verbose = FALSE)
      .align_genes(MERINGUE::spatialCrossCorMatrix(mat = t(.rows(logY, idx)), weight = w), genes)
    }, membership, global)
  })
  .finish(g, nrow(counts), tm)
}

# spMOCA: spatial-kernel correlation estimate on log-CPM, refit per group
.spmoca_cormat <- function(counts, logcounts, coords) {
  genes  <- colnames(counts)
  object <- spMOCA::createSPMOCAObject(spatial_countMat = t(as.matrix(counts)),
                                       spatial_location = as.data.frame(coords),
                                       feature_list     = genes)
  object@spatial_normCount <- t(logcounts)
  object <- spMOCA::createSpatialKernel(object)
  object <- spMOCA::estimatespMOCA(object)
  as.matrix(object@corr_est)
}

compare_mat_spMOCA <- function(counts, coords, data, libsize = "libsize", membership = NULL,
                               global = TRUE, ...) {
  genes <- colnames(counts)
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    g <- .grouped_mat(function(idx)
      .align_genes(.spmoca_cormat(.rows(counts, idx), .rows(logY, idx), .rows(coords, idx)), genes),
      membership, global)
  })
  .finish(g, nrow(counts), tm)
}

# SpaceDecorr: one NB-GAM fit on all cells with the library size as an offset,
# correlation of the residuals per group; k = spline basis dimension
compare_mat_spacedecorr <- function(counts, coords, data, libsize = "libsize", membership = NULL,
                                    k = 200, basis = "ts", family = "nb", ...) {
  genes <- colnames(counts); n <- nrow(counts)
  ksplines <- if (!is.null(k) && k < 1) round(k * n) else k
  md <- as.data.frame(data)
  md$libsize_latent <- md[[libsize]]
  if (!"libsize" %in% names(md)) md$libsize <- md[[libsize]]
  form <- sprintf("y ~ offset(log(libsize_latent)) + s(sdimx, sdimy, bs='%s', k=%d)", basis, ksplines)
  tm <- system.time({
    fit <- spacedecorr::spacedecorr(
      assay_matrix = as.data.frame(counts), metadata = md,
      libsize_col = libsize, covariates = "libsize_latent", formula_override = form,
      k = ksplines, family = family, basis = basis,
      nCores = 1, verbose = FALSE, return_splines = FALSE)
    resid <- as.matrix(fit$residuals)[, genes, drop = FALSE]
    g <- .grouped_mat(function(idx) .align_genes(cor(.rows(resid, idx)), genes), membership)
  })
  .finish(g, n, tm, resid = resid)
}

# SpaceX: shared and cluster-specific networks from one fit. Two runtime fixes:
# patch_spacex() (call it once after sourcing; needs SpaceX and PQLseq installed)
# replaces a `class() ==` test that errors on R >= 4.2, and spacex_networks()
# computes the posterior-mean networks in R instead of the package's post-processing.
# Cells must arrive ordered by membership.
patch_spacex <- function() {
  loadNamespace("PQLseq")
  f <- SpaceX:::PQLseq.AI; bt <- deparse(body(f))
  bt <- sub('class(P) == "try-error"', 'inherits(P, "try-error")', bt, fixed = TRUE)
  body(f) <- parse(text = paste(bt, collapse = "\n"))[[1]]
  assignInNamespace("PQLseq.AI", f, ns = "SpaceX")
}

spacex_networks <- function(fit, o = NULL) {
  Phi <- fit$Posterior_samples$Phi          # G x factors x draws
  Lam <- fit$Posterior_samples$Lambda       # one array per cluster
  ns  <- dim(Phi)[3]; G <- dim(Phi)[1]; L <- length(Lam)
  o   <- if (is.null(o)) ns else min(o, ns)
  cor_of <- function(M) { d <- sqrt(diag(M)); d[d == 0] <- NA; M / tcrossprod(d) }
  avg    <- function(f) Reduce(`+`, lapply(seq_len(o), f)) / o

  shared  <- avg(function(t) cor_of(tcrossprod(Phi[, , t])))
  cluster <- array(0, c(G, G, L), dimnames = list(rownames(Phi), rownames(Phi), names(Lam)))
  for (l in seq_len(L))
    cluster[, , l] <- avg(function(t) cor_of(tcrossprod(Phi[, , t]) + tcrossprod(Lam[[l]][, , t])))
  list(Shared_network = shared, Cluster_network = cluster)
}

compare_mat_spacex <- function(counts, coords, data, libsize = "libsize", membership = NULL,
                               numCore = 1, ...) {
  stopifnot(!is.null(membership), !is.unsorted(membership), length(unique(membership)) == 2)
  genes  <- colnames(counts)
  labels <- c("niche0", "niche1")[membership + 1]
  tm <- system.time({
    fit <- suppressMessages(SpaceX::SpaceX(
      Gene_expression_mat = as.data.frame(as.matrix(counts)),
      Spatial_locations   = as.data.frame(coords)[, c("sdimx", "sdimy")],
      Cluster_annotations = labels,
      sPMM = FALSE, Post_process = FALSE, numCore = numCore))
    nets <- spacex_networks(fit)
  })
  stopifnot(dim(nets$Cluster_network)[3] == 2)
  nm <- function(M) { dimnames(M) <- list(genes, genes); M }
  list(cor    = .align_genes(nm(nets$Shared_network), genes),
       cor_g0 = .align_genes(nm(nets$Cluster_network[, , 1]), genes),
       cor_g1 = .align_genes(nm(nets$Cluster_network[, , 2]), genes),
       n = nrow(counts), n_g0 = sum(membership == 0), n_g1 = sum(membership == 1),
       time = unname(tm["elapsed"]), resid = NULL)
}

compare_mat_methods <- list(
  logpearson  = compare_mat_logpearson,
  spark       = compare_mat_spark,
  giotto      = compare_mat_giotto,
  meringue    = compare_mat_meringue,
  spMOCA      = compare_mat_spMOCA,
  spacedecorr = compare_mat_spacedecorr,
  spacex      = compare_mat_spacex)

# Per-pair rows ----
# One row per upper-triangle gene pair: estimate, raw and BH-adjusted p-value and
# confidence interval for the all-cells matrix (mat_edge_test), each group and
# their difference (diffCorNet), plus the truth.
compare_mat_result_cols <- c(
  "model", "gene_i", "gene_j", "block_i", "block_j",
  "true_edge0", "true_edge1", "true_diff", "rho_true0", "rho_true1",
  "cor", "pvalue", "padj", "ci_low", "ci_high",
  "cor_g0", "pvalue_g0", "padj_g0", "ci_low_g0", "ci_high_g0",
  "cor_g1", "pvalue_g1", "padj_g1", "ci_low_g1", "ci_high_g1",
  "cor_diff", "pvalue_diff", "padj_diff", "ci_low_diff", "ci_high_diff",
  "n_g0", "n_g1", "time")

mat_pair_rows <- function(model, res, truth) {
  genes <- truth$genes; p <- length(genes)
  ut <- which(upper.tri(diag(p)), arr.ind = TRUE)
  gi <- ut[, 1]; gj <- ut[, 2]; npair <- nrow(ut)
  up  <- function(M) if (is.null(M)) rep(NA_real_, npair) else { M <- .align_genes(M, genes); M[upper.tri(M)] }
  upl <- function(M) if (is.null(M)) rep(NA_integer_, npair) else as.integer(.align_genes(M, genes)[upper.tri(M)])

  g <- if (!is.null(res$cor)) mat_edge_test(.align_genes(res$cor, genes), res$n) else NULL
  d <- if (!is.null(res$cor_g0) && !is.null(res$cor_g1))
    diffCorNet(.align_genes(res$cor_g0, genes), .align_genes(res$cor_g1, genes), res$n_g0, res$n_g1) else NULL
  cor_g0 <- up(res$cor_g0); cor_g1 <- up(res$cor_g1)

  out <- data.frame(
    model = model, gene_i = genes[gi], gene_j = genes[gj],
    block_i = if (is.null(truth$labels)) NA_integer_ else truth$labels[gi],
    block_j = if (is.null(truth$labels)) NA_integer_ else truth$labels[gj],
    true_edge0 = upl(truth$true_edge0), true_edge1 = upl(truth$true_edge1), true_diff = upl(truth$true_diff),
    rho_true0 = up(truth$rho_true0), rho_true1 = up(truth$rho_true1),
    cor = up(res$cor), pvalue = up(g$p), padj = up(g$padj), ci_low = up(g$ci_low), ci_high = up(g$ci_high),
    cor_g0 = cor_g0, pvalue_g0 = up(d$p0), padj_g0 = up(d$padj0), ci_low_g0 = up(d$ci_low0), ci_high_g0 = up(d$ci_high0),
    cor_g1 = cor_g1, pvalue_g1 = up(d$p1), padj_g1 = up(d$padj1), ci_low_g1 = up(d$ci_low1), ci_high_g1 = up(d$ci_high1),
    cor_diff = cor_g1 - cor_g0, pvalue_diff = up(d$p_diff), padj_diff = up(d$padj_diff),
    ci_low_diff = up(d$ci_low_diff), ci_high_diff = up(d$ci_high_diff),
    n_g0 = if (is.null(res$n_g0)) NA_integer_ else res$n_g0,
    n_g1 = if (is.null(res$n_g1)) NA_integer_ else res$n_g1,
    time = res$time, stringsAsFactors = FALSE)
  out[, compare_mat_result_cols]
}
