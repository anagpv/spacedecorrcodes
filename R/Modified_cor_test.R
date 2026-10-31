#------------------------------------------------------------------------------#
#- Modified t-test for correlations between spatially autocorrelated genes    -#
#- (Clifford, Richardson & Hemon 1989; Dutilleul 1993)                        -#
#-                                                                            -#
#- For every gene pair: Pearson r, standard t-test P-value (cells independent)-#
#- effective sample size (ESS), and modified-t P-value using the ESS.         -#
#-                                                                            -#
#- ESS = 1 + n^2 / (n + sum_k N_k * rho_x(k) * rho_y(k))                      -#
#-   N_k   = number of ordered cell pairs in distance class k                 -#
#-   rho(k)= spatial autocorrelation of a gene in class k                     -#
#-           (mean product of standardized values over pairs in the class)    -#
#- Scalable estimation: classes up to `short_max` use all neighbors of a      -#
#- random subset of `n_query` cells; longer classes use `n_rand` random pairs.-#
#------------------------------------------------------------------------------#
library(RANN)

modified_cor_test <- function(X, coords,
                              dist_breaks = c(0, 25, 50, 100, 200, 400, 800, 1600, Inf),
                              short_max   = 200,     # must be one of dist_breaks
                              n_query     = 20000,   # cells used as centers for short-range pairs
                              k_max       = 1000,    # max neighbors returned per cell
                              n_rand      = 3e6,     # random pairs for long-range classes
                              max_pairs   = 3e5,     # pairs per class used to estimate rho
                              chunk       = 5e4,     # pairs per chunk (memory)
                              seed        = 1) {
  stopifnot(short_max %in% dist_breaks, nrow(X) == nrow(coords))
  set.seed(seed)
  X <- as.matrix(X); coords <- as.matrix(coords)
  n <- as.numeric(nrow(X)); G <- ncol(X); genes <- colnames(X)
  K <- length(dist_breaks) - 1
  
  # standardized expression and correlation
  Z <- scale(X); Z[is.na(Z)] <- 0
  r <- crossprod(Z) / (n - 1)
  
  # short-range pairs: all neighbors (within short_max) of a random subset of cells
  q  <- if (n > n_query) sample.int(n, n_query) else 1:n
  nn <- nn2(coords, coords[q, , drop = FALSE], k = min(k_max, n), searchtype = "radius", radius = short_max)
  if (max(rowSums(nn$nn.idx > 0)) >= min(k_max, n))
    warning("k_max reached for some cells: increase k_max or decrease short_max")
  pi_s <- rep(q, ncol(nn$nn.idx)); pj_s <- as.vector(nn$nn.idx); pd_s <- as.vector(nn$nn.dists)
  keep <- pj_s > 0 & pj_s != pi_s
  pi_s <- pi_s[keep]; pj_s <- pj_s[keep]; pd_s <- pd_s[keep]
  rm(nn)
  
  # long-range pairs: uniform random pairs of cells
  pi_l <- sample.int(n, n_rand, replace = TRUE); pj_l <- sample.int(n, n_rand, replace = TRUE)
  keep <- pi_l != pj_l
  pi_l <- pi_l[keep]; pj_l <- pj_l[keep]
  pd_l <- sqrt(rowSums((coords[pi_l, , drop = FALSE] - coords[pj_l, , drop = FALSE])^2))
  
  # pair counts (N_k) and autocorrelation (rho, K x G) per distance class
  N_k <- rep(NA, K)
  rho <- matrix(NA, K, G, dimnames = list(NULL, genes))
  for (k in 1:K) {
    lo <- dist_breaks[k]; hi <- dist_breaks[k + 1]
    if (hi <= short_max) {
      id <- which(pd_s > lo & pd_s <= hi)
      N_k[k] <- length(id) * n / length(q)                   # scaled to all cells
      a <- pi_s; b <- pj_s
    } else {
      id <- which(pd_l > lo & pd_l <= hi)
      N_k[k] <- length(id) / length(pd_l) * n * (n - 1)      # fraction of all ordered pairs
      a <- pi_l; b <- pj_l
    }
    if (length(id) == 0) { N_k[k] <- 0; rho[k, ] <- 0; next }
    if (length(id) > max_pairs) id <- sample(id, max_pairs)
    s <- 0
    for (start in seq(1, length(id), by = chunk)) {
      idc <- id[start:min(start + chunk - 1, length(id))]
      s <- s + colSums(Z[a[idc], , drop = FALSE] * Z[b[idc], , drop = FALSE])
    }
    rho[k, ] <- s / length(id)
  }
  
  # effective sample size per pair and P-values
  # den <= 0 happens when the two genes' autocorrelations have opposite signs and
  # the cross term exceeds n in absolute value: the true ESS is then above n, so
  # the estimate is set to n (no reduction) instead of being clamped to 3
  den <- n + t(rho) %*% (N_k * rho)
  ess <- ifelse(den > 0, 1 + n^2 / den, n)
  ess <- pmin(pmax(ess, 3), n)
  n_invalid <- sum(den[upper.tri(den)] <= 0)
  if (n_invalid > 0) message(n_invalid, " pairs with non-positive ESS denominator set to ESS = n")
  # standard t-test for a correlation (same as cor.test); the modified test below
  # is the same test with n replaced by the effective sample size
  p_naive    <- 2 * pt(-abs(r) * sqrt((n - 2) / (1 - r^2)), n - 2)
  p_modified <- pf(r^2 * (ess - 2) / (1 - r^2), 1, ess - 2, lower.tail = FALSE)
  
  ut <- upper.tri(r); ij <- which(ut, arr.ind = TRUE)
  list(results = data.frame(gene1 = genes[ij[, 1]], gene2 = genes[ij[, 2]],
                            r = r[ut], ess = ess[ut],
                            p_naive = p_naive[ut], p_modified = p_modified[ut]),
       n = n, n_invalid = n_invalid,
       classes = data.frame(lo = dist_breaks[-(K + 1)], hi = dist_breaks[-1], N_pairs = N_k),
       rho = rho)
}

#-----------------------------------------------------------------------------#
#- Example use                                                                -#
#-----------------------------------------------------------------------------#
# Liver (one section, zonated hepatocytes): coordinates in µm
# out_sd  <- modified_cor_test(res[cells, genes], coords)
# out_vst <- modified_cor_test(vst[cells, genes], coords)
#
# NSCLC: one call per method, same cells and coordinates for all methods;
# set dist_breaks/short_max in the units of the coordinates (e.g. µm or mm)
# out <- lapply(list(spacedecorr = R_sd, vst = R_vst, lognorm = R_log, giotto = R_giotto),
#               function(M) modified_cor_test(M, coords_nsclc))
#
# Summaries
# res <- bind_rows(sd = out_sd$results, vst = out_vst$results, .id = "method") %>%
#   group_by(method) %>%
#   mutate(fdr_naive = p.adjust(p_naive, "BH"), fdr_modified = p.adjust(p_modified, "BH"))
# res %>% summarise(n_cells = out_sd$n, ess_median = median(ess),
#                   ess_q10 = quantile(ess, 0.1), ess_q90 = quantile(ess, 0.9),
#                   sig_naive = sum(fdr_naive < 0.05), sig_modified = sum(fdr_modified < 0.05))