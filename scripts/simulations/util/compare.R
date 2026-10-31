# compare.R - one wrapper per method for the correlation of a gene pair on one
# tissue, globally and within each of two groups; run_compare_method() returns
# the methods' estimates, p-values and confidence intervals as one wide row.
library(tidyverse)
library(spacedecorr)
library(MERINGUE)
library(GraphR)
library(spMOCA)

# GraphR 1.x: pred_ind() refers to a variable under the wrong name.
pred_ind_fixed <- GraphR:::pred_ind
body(pred_ind_fixed) <- as.call(c(
  as.name("{"),
  lapply(as.list(body(pred_ind_fixed))[-1], function(e)
    eval(bquote(substitute(.(e), list(w_diag_ind = quote(omega_diag_ind))))))
))
assignInNamespace("pred_ind", pred_ind_fixed, ns = "GraphR")

# Statistics ----

# Cauchy combination of p-values.
cct <- function(pvec, weight = NULL) {
  pvec <- pvec[is.finite(pvec)]
  if (length(pvec) == 0) return(NA_real_)
  J <- length(pvec)
  if (is.null(weight)) weight <- rep(1 / J, J)
  t0 <- sum(tan((0.5 - pvec) * pi) * weight)
  0.5 - atan(t0 / sum(weight)) / pi
}

# Fisher-z confidence interval of a correlation of n observations.
fisher_ci <- function(r, n, conf_level = 0.95) {
  if (is.na(r) || is.na(n) || n <= 3)
    return(c(ci_low = NA_real_, ci_high = NA_real_))
  r    <- min(max(r, -1 + 1e-6), 1 - 1e-6)
  z    <- atanh(r)
  se   <- 1 / sqrt(n - 3)
  crit <- qnorm(1 - (1 - conf_level) / 2)
  c(ci_low = tanh(z - crit * se), ci_high = tanh(z + crit * se))
}

# Two-sample Fisher-z test of r1 - r0; the interval is on the Fisher-z scale.
.two_cor_test <- function(r0, n0, r1, n1, conf.level = 0.95) {
  if (any(is.na(c(r0, n0, r1, n1))) || n0 <= 3 || n1 <= 3)
    return(list(diff = NA_real_, p.value = NA_real_, ci = c(NA_real_, NA_real_)))
  z0   <- atanh(min(max(r0, -1 + 1e-6), 1 - 1e-6))
  z1   <- atanh(min(max(r1, -1 + 1e-6), 1 - 1e-6))
  se   <- sqrt(1 / (n0 - 3) + 1 / (n1 - 3))
  crit <- qnorm(1 - (1 - conf.level) / 2)
  list(diff    = r1 - r0,
       p.value = 2 * (1 - pnorm(abs((z1 - z0) / se))),
       ci      = c((z1 - z0) - crit * se, (z1 - z0) + crit * se))
}

# A cor.test()-shaped result for a method that returns a point estimate only:
# t-approximation p-value unless `p.value` is given, Fisher-z interval.
.as_htest <- function(estimate, n = NA_real_, p.value = NULL,
                      method = "spatial correlation", conf.level = 0.95, data.name = "Y1 and Y2") {
  r    <- unname(estimate)[1]
  stat <- c(t = unname(r / sqrt(1 - r^2) * sqrt(n - 2)))
  if (is.null(p.value)) p.value <- 2 * pt(-abs(stat), df = n - 2)
  ci <- as.numeric(fisher_ci(r, n, conf.level))
  attr(ci, "conf.level") <- conf.level
  structure(list(statistic = stat, parameter = c(df = n - 2), p.value = unname(p.value)[1],
                 estimate = c(cor = r), null.value = c(correlation = 0),
                 alternative = "two.sided", method = method, data.name = data.name,
                 conf.int = ci), class = "htest")
}

# Run `compute(idx)` on all cells (idx = NULL) and, with a 0/1 `membership`, on
# each group. A group of 3 cells or fewer gives NULL; a failing group fails the method.
.grouped <- function(compute, membership) {
  est <- compute(NULL)
  if (is.null(membership))
    return(list(est = est, est_g0 = NULL, est_g1 = NULL, n_g0 = NA_integer_, n_g1 = NA_integer_))
  i0 <- which(membership == 0)
  i1 <- which(membership == 1)
  per_group <- function(idx) if (length(idx) > 3) compute(idx) else NULL
  list(est = est, est_g0 = per_group(i0), est_g1 = per_group(i1),
       n_g0 = length(i0), n_g1 = length(i1))
}

# Normalisation ----

# log-CPM, log(1e5 * count / libsize + 1); `libsize` names the column of `data`.
.log_cpm <- function(counts, data = NULL, libsize = "libsize") {
  counts <- as.matrix(counts)
  sz <- if (!is.null(libsize) && is.character(libsize) && !is.null(data) && libsize %in% colnames(data))
    data[[libsize]] else rowSums(counts)
  out <- apply(counts, 2, function(x) log(10e4 * x / sz + 1))
  dimnames(out) <- dimnames(counts)
  out
}

# sctransform Pearson residuals, cells x genes.
.vst <- function(counts, data, libsize = "libsize")
  t(sctransform::vst(umi = t(counts), cell_attr = data, latent_var = libsize)$y)

# Methods ----
# Each takes (counts, coords, data, libsize, features, k, membership) and returns
# list(est, est_g0, est_g1, n_g0, n_g1, time).

pearson_pair <- function(x, y, membership)
  .grouped(function(idx) {
    if (is.null(idx)) cor.test(x, y, method = "pearson")
    else              cor.test(x[idx], y[idx], method = "pearson")
  }, membership)

# log-normalized: Pearson on log-CPM
compare_logpearson <- function(counts, coords = NULL, data = NULL, libsize = "libsize",
                               features = c(1, 2), k = 200, membership = NULL, ...) {
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    g <- pearson_pair(logY[, features[1]], logY[, features[2]], membership)
  })
  c(g, time = unname(tm["elapsed"]))
}

# VST-normalized: Pearson on sctransform residuals
compare_spark <- function(counts, coords = NULL, data = NULL, libsize = "libsize",
                          features = c(1, 2), k = 200, membership = NULL, ...) {
  tm <- system.time({
    Yvst <- .vst(counts, data, libsize)
    g <- pearson_pair(Yvst[, features[1]], Yvst[, features[2]], membership)
  })
  c(g, time = unname(tm["elapsed"]))
}

# Giotto: Delaunay-network spatial correlation, refit per group
compare_giotto <- function(counts, coords, data = NULL, libsize = "libsize",
                           features = c(1, 2), k = 200, minimum_k = 2, membership = NULL, ...) {
  f <- colnames(counts)[features]
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    giotto_r <- function(idx) {
      ct <- if (is.null(idx)) counts else counts[idx, , drop = FALSE]
      lg <- if (is.null(idx)) logY   else logY[idx, , drop = FALSE]
      cd <- if (is.null(idx)) coords else coords[idx, , drop = FALSE]
      suppressWarnings({
        gobject <- Giotto::createGiottoObject(
          expression      = list(raw = t(as.matrix(ct)), normalized = t(lg)),
          expression_feat = "rna",
          spatial_locs    = as.data.frame(cd))
        gobject <- Giotto::createSpatialNetwork(gobject = gobject, minimum_k = minimum_k)
        cor_dt  <- Giotto::detectSpatialCorFeats(
          gobject, method = "network", spatial_network_name = "Delaunay_network")$cor_DT
      })
      sel <- (cor_dt$feat_ID == f[1] & cor_dt$variable == f[2]) |
             (cor_dt$feat_ID == f[2] & cor_dt$variable == f[1])
      as.numeric(cor_dt$spat_cor[sel])[1]
    }
    g <- .grouped(function(idx)
      .as_htest(giotto_r(idx), if (is.null(idx)) nrow(counts) else length(idx),
                method = "Giotto detectSpatialCorFeats (network)", data.name = paste(f[1], "and", f[2])),
      membership)
  })
  c(g, time = unname(tm["elapsed"]))
}

# GraphR: spatial partial correlation, one fit on all cells, per-location
# predictions averaged within each group, CCT p-value
compare_GraphR <- function(counts, coords, data = NULL, libsize = "libsize",
                           features = c(1, 2), k = 200, membership = NULL, ...) {
  f <- colnames(counts)[features]
  n <- nrow(counts)
  mk <- function(r, pv, nn)
    .as_htest(r, nn, p.value = pv, method = "GraphR (spatial partial correlation, CCT p-value)",
              data.name = paste(f[1], "and", f[2]))
  summ <- function(df, nn) {
    if (nrow(df) == 0) return(mk(NA_real_, NA_real_, nn))
    mk(mean(df$CorrPIP), cct(df$FDR_p), nn)
  }

  tm <- system.time({
    Yvst <- .vst(counts, data, libsize)
    res  <- tryCatch(GraphR::GraphR_est(features = as.matrix(Yvst), cont_external = as.matrix(coords)),
                     error = function(e) e)
    if (inherits(res, "error")) {
      g <- list(est = mk(NA_real_, NA_real_, n), est_g0 = NULL, est_g1 = NULL,
                n_g0 = NA_integer_, n_g1 = NA_integer_)
    } else {
      pred <- GraphR::GraphR_pred(as.data.frame(coords), res)
      cc   <- as.data.frame(coords)
      pair <- function(d) dplyr::filter(d, feature1 %in% f, feature2 %in% f, feature1 != feature2)
      base <- pred %>%
        dplyr::left_join(cc, by = c("sdimx", "sdimy")) %>%
        dplyr::mutate(CorrPIP = Correlation * Pr_inclusion) %>%
        pair()
      est <- summ(base, n)
      if (is.null(membership)) {
        g <- list(est = est, est_g0 = NULL, est_g1 = NULL, n_g0 = NA_integer_, n_g1 = NA_integer_)
      } else {
        keyed <- pred %>%
          dplyr::left_join(cbind(cc, .mem = membership), by = c("sdimx", "sdimy")) %>%
          dplyr::mutate(CorrPIP = Correlation * Pr_inclusion) %>%
          pair()
        i0 <- which(membership == 0); i1 <- which(membership == 1)
        g <- list(est    = est,
                  est_g0 = summ(dplyr::filter(keyed, .mem == 0), length(i0)),
                  est_g1 = summ(dplyr::filter(keyed, .mem == 1), length(i1)),
                  n_g0 = length(i0), n_g1 = length(i1))
      }
    }
  })
  c(g, time = unname(tm["elapsed"]))
}

# MERINGUE: spatial cross-correlation, refit per group
compare_meringue <- function(counts, coords, data = NULL, libsize = "libsize",
                             features = c(1, 2), k = 200, filterDist = 0.1, membership = NULL, ...) {
  f <- colnames(counts)[features]
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    mer_r <- function(idx) {
      lg <- if (is.null(idx)) logY   else logY[idx, , drop = FALSE]
      cd <- if (is.null(idx)) coords else coords[idx, , drop = FALSE]
      w   <- MERINGUE::getSpatialNeighbors(as.data.frame(cd), filterDist = filterDist)
      ccm <- MERINGUE::spatialCrossCorMatrix(mat = t(lg), weight = w)
      as.numeric(ccm[f[1], f[2]])
    }
    g <- .grouped(function(idx)
      .as_htest(mer_r(idx), if (is.null(idx)) nrow(counts) else length(idx),
                method = "MERINGUE spatialCrossCorMatrix", data.name = paste(f[1], "and", f[2])),
      membership)
  })
  c(g, time = unname(tm["elapsed"]))
}

# SpaceDecorr: Pearson on the residuals of one fit on all cells; k = spline basis
# dimension (k < 1 is read as a fraction of the cells)
compare_spacedecorr <- function(counts, coords = NULL, data, libsize = "libsize",
                                features = c(1, 2), k = 200, basis = "ts", family = "nb",
                                membership = NULL, ...) {
  n        <- nrow(counts)
  ksplines <- if (!is.null(k) && k < 1) round(k * n) else k
  tm <- system.time({
    fit <- spacedecorr::spacedecorr(
      assay_matrix = as.data.frame(counts[, features, drop = FALSE]),
      metadata     = data,
      basis        = basis,
      libsize_col  = libsize,
      family       = family,
      k            = ksplines,
      verbose      = FALSE)
    g <- pearson_pair(as.numeric(fit$residuals[, 1]), as.numeric(fit$residuals[, 2]), membership)
  })
  c(g, time = unname(tm["elapsed"]))
}

# spMOCA: spatial-kernel correlation estimate, refit per group
compare_spMOCA <- function(counts, coords, data = NULL, libsize = "libsize",
                           features = c(1, 2), k = 200, membership = NULL, ...) {
  genes <- colnames(counts)
  f     <- genes[features]
  tm <- system.time({
    logY <- .log_cpm(counts, data, libsize)
    moca_r <- function(idx) {
      ct <- if (is.null(idx)) counts else counts[idx, , drop = FALSE]
      lg <- if (is.null(idx)) logY   else logY[idx, , drop = FALSE]
      cd <- if (is.null(idx)) coords else coords[idx, , drop = FALSE]
      object <- spMOCA::createSPMOCAObject(spatial_countMat = t(as.matrix(ct)),
                                           spatial_location = as.data.frame(cd),
                                           feature_list     = genes)
      if (ncol(ct) < 100) object@spatial_normCount <- t(lg) else object <- normalizeSPMOCA(object)
      object <- spMOCA::createSpatialKernel(object)
      object <- spMOCA::estimatespMOCA(object)
      as.numeric(object@corr_est[features[1], features[2]])
    }
    g <- .grouped(function(idx)
      .as_htest(moca_r(idx), if (is.null(idx)) nrow(counts) else length(idx),
                method = "spMOCA", data.name = paste(f[1], "and", f[2])),
      membership)
  })
  c(g, time = unname(tm["elapsed"]))
}

# Dispatch ----

compare_methods <- list(
  logpearson  = compare_logpearson,
  spark       = compare_spark,
  giotto      = compare_giotto,
  GraphR      = compare_GraphR,
  meringue    = compare_meringue,
  spacedecorr = compare_spacedecorr,
  spMOCA      = compare_spMOCA)

compare_result_cols <- c(
  "model",
  "cor", "pvalue", "ci_low", "ci_high",
  "cor_g0", "pvalue_g0", "ci_low_g0", "ci_high_g0",
  "cor_g1", "pvalue_g1", "ci_low_g1", "ci_high_g1",
  "cor_diff", "pvalue_diff", "ci_low_diff", "ci_high_diff",
  "n_g0", "n_g1", "time", "status", "error_message")

.htest_row <- function(ht) {
  if (is.null(ht))
    return(c(cor = NA_real_, pvalue = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  ci <- tryCatch(as.numeric(ht$conf.int), error = function(e) c(NA_real_, NA_real_))
  if (length(ci) < 2) ci <- c(NA_real_, NA_real_)
  c(cor = unname(ht$estimate)[1], pvalue = unname(ht$p.value)[1], ci_low = ci[1], ci_high = ci[2])
}

.compare_row <- function(model, gl, g0, g1, d, n_g0, n_g1, time, status, error_message) {
  data.frame(
    model        = model,
    cor          = gl[["cor"]],    pvalue       = gl[["pvalue"]],
    ci_low       = gl[["ci_low"]], ci_high      = gl[["ci_high"]],
    cor_g0       = g0[["cor"]],    pvalue_g0    = g0[["pvalue"]],
    ci_low_g0    = g0[["ci_low"]], ci_high_g0   = g0[["ci_high"]],
    cor_g1       = g1[["cor"]],    pvalue_g1    = g1[["pvalue"]],
    ci_low_g1    = g1[["ci_low"]], ci_high_g1   = g1[["ci_high"]],
    cor_diff     = d[["diff"]],    pvalue_diff  = d[["p.value"]],
    ci_low_diff  = d[["ci"]][1],   ci_high_diff = d[["ci"]][2],
    n_g0 = n_g0, n_g1 = n_g1, time = time, status = status, error_message = error_message,
    stringsAsFactors = FALSE)
}

# One wide row per method: the global estimate (cor, pvalue, ci), the estimates
# within group 0 / 1 (*_g0, *_g1) and the Fisher-z test of their difference
# (*_diff, interval on the z scale). Without `membership` the group columns are NA.
run_compare_method <- function(name, counts, coords = NULL, data = NULL, libsize = "libsize",
                               features = c(1, 2), membership = NULL, ...) {
  if (!name %in% names(compare_methods))
    stop(sprintf("Unknown method '%s'. Available: %s", name, paste(names(compare_methods), collapse = ", ")))

  fn  <- compare_methods[[name]]
  res <- tryCatch(
    fn(counts = counts, coords = coords, data = data, libsize = libsize,
       features = features, membership = membership, ...),
    error = function(e) e)

  na4 <- c(cor = NA_real_, pvalue = NA_real_, ci_low = NA_real_, ci_high = NA_real_)
  nad <- list(diff = NA_real_, p.value = NA_real_, ci = c(NA_real_, NA_real_))

  if (inherits(res, "condition"))
    return(.compare_row(name, na4, na4, na4, nad, NA_integer_, NA_integer_,
                        NA_real_, "error", conditionMessage(res)))

  gl <- .htest_row(res$est)
  g0 <- .htest_row(res$est_g0)
  g1 <- .htest_row(res$est_g1)
  d  <- if (!is.null(res$est_g0) && !is.null(res$est_g1))
    .two_cor_test(g0[["cor"]], res$n_g0, g1[["cor"]], res$n_g1) else nad

  .compare_row(name, gl, g0, g1, d, res$n_g0, res$n_g1, unname(res$time)[1], "ok", NA_character_)
}
