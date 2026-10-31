# util.R - writing and resuming result files; the tests on p x p correlation
# matrices (block and differential-network benchmarks); ARI of a matrix.

# Results ----

# A named list of per-draw realisations -> one "key=value; ..." string.
encode_ext_info <- function(lst) {
  if (length(lst) == 0) return("")
  paste(vapply(names(lst), function(nm) paste0(nm, "=", paste(lst[[nm]], collapse = ",")), character(1)),
        collapse = "; ")
}

# Refuse to append to a file whose header differs from `expected_cols`.
validate_or_init_csv <- function(path, expected_cols) {
  if (!file.exists(path) || file.size(path) == 0) return(invisible())
  header <- names(read.csv(path, nrows = 0, check.names = FALSE))
  if (!identical(header, expected_cols))
    stop(sprintf("'%s' has columns\n  %s\nbut this script writes\n  %s", path,
                 paste(header, collapse = ", "), paste(expected_cols, collapse = ", ")))
}

# One string per row of df: the identity of a setting.
setting_key <- function(df, cols)
  do.call(paste, c(lapply(cols, function(cl) df[[cl]]), sep = "\r"))

# Drop the settings whose every method is already in out_csv. A setting holding
# some but not all methods stops the run (it cannot be resumed safely).
filter_pending_settings <- function(settings, out_csv, setting_cols_csv, setting_cols_grid,
                                    methods, model_col = "model") {
  stopifnot(length(setting_cols_csv) == length(setting_cols_grid))
  if (!file.exists(out_csv) || file.size(out_csv) == 0) return(settings)
  done <- read.csv(out_csv, stringsAsFactors = FALSE)
  if (!nrow(done)) return(settings)
  seen <- unique(data.frame(key = setting_key(done, setting_cols_csv),
                            model = as.character(done[[model_col]]), stringsAsFactors = FALSE))
  methods_by_key <- split(seen$model, seen$key)
  want     <- sort(unique(methods))
  grid_key <- setting_key(settings, setting_cols_grid)
  status   <- vapply(grid_key, function(kk) {
    have <- intersect(want, methods_by_key[[kk]])
    if (length(have) == 0L) "run" else if (setequal(have, want)) "done" else "partial"
  }, character(1), USE.NAMES = FALSE)
  if (any(status == "partial"))
    stop(sprintf("%d setting(s) in %s hold some but not all of the methods; delete their rows to rerun them",
                 sum(status == "partial"), out_csv))
  settings[status == "run", , drop = FALSE]
}

# A function that appends rows to out_csv under a file lock (several workers
# write to one file); the header is written with the first rows.
make_appender <- function(out_csv, result_cols) {
  function(df) {
    df <- df[, result_cols, drop = FALSE]
    lk <- filelock::lock(paste0(out_csv, ".lock"))
    on.exit(filelock::unlock(lk))
    fresh <- !file.exists(out_csv) || file.size(out_csv) == 0
    write.table(df, out_csv, sep = ",", append = !fresh, col.names = fresh,
                row.names = FALSE, qmethod = "double")
  }
}

# Matrix tests ----

# BH-adjust the upper triangle of a p-value matrix, mirror it.
adjust_upper_triangle <- function(pmat, adjust.method = "fdr") {
  pmat <- as.matrix(pmat)
  ut   <- upper.tri(pmat)
  pmat[ut] <- p.adjust(pmat[ut], method = adjust.method)
  pmat[lower.tri(pmat)] <- t(pmat)[lower.tri(pmat)]
  pmat
}

fisher_ci_vec <- function(r, n, conf_level = 0.95) {
  r    <- pmin(pmax(r, -1 + 1e-6), 1 - 1e-6)
  z    <- atanh(r)
  se   <- 1 / sqrt(n - 3)
  crit <- qnorm(1 - (1 - conf_level) / 2)
  ok   <- !is.na(r) & !is.na(n) & n > 3
  list(ci_low  = ifelse(ok, tanh(z - crit * se), NA_real_),
       ci_high = ifelse(ok, tanh(z + crit * se), NA_real_))
}

# Is each correlation of a p x p matrix non-zero: t-test on n - 2 df, BH over
# the upper triangle, Fisher-z interval.
mat_edge_test <- function(cor_mat, n) {
  r  <- as.matrix(cor_mat)
  diag(r) <- NA_real_
  tt <- r * sqrt((n - 2) / (1 - r^2))
  p  <- 2 * pt(-abs(tt), df = n - 2)
  ci <- fisher_ci_vec(r, n)
  shape <- function(x) { x <- matrix(x, nrow(r), ncol(r)); dimnames(x) <- dimnames(r); x }
  list(p = p, padj = adjust_upper_triangle(p), ci_low = shape(ci$ci_low), ci_high = shape(ci$ci_high))
}

# Two correlation matrices (group 0 on n0 cells, group 1 on n1): Fisher-z test
# of every pairwise difference (group 1 - group 0; interval on the z scale) and
# the single-group tests of mat_edge_test().
diffCorNet <- function(cor_mat0, cor_mat1, n0, n1, conf_level = 0.95) {
  cor_mat0 <- as.matrix(cor_mat0); cor_mat1 <- as.matrix(cor_mat1)
  diag(cor_mat0) <- diag(cor_mat1) <- NA_real_
  z0 <- atanh(cor_mat0)
  z1 <- atanh(cor_mat1)
  se_diff <- sqrt(1 / (n0 - 3) + 1 / (n1 - 3))
  z_diff  <- z1 - z0
  p_diff  <- 2 * pnorm(-abs(z_diff / se_diff))
  crit    <- qnorm(1 - (1 - conf_level) / 2)
  t0 <- mat_edge_test(cor_mat0, n0)
  t1 <- mat_edge_test(cor_mat1, n1)
  list(p0 = t0$p, padj0 = t0$padj, ci_low0 = t0$ci_low, ci_high0 = t0$ci_high,
       p1 = t1$p, padj1 = t1$padj, ci_low1 = t1$ci_low, ci_high1 = t1$ci_high,
       p_diff = p_diff, padj_diff = adjust_upper_triangle(p_diff),
       ci_low_diff = z_diff - crit * se_diff, ci_high_diff = z_diff + crit * se_diff)
}

# ARI of hierarchical clustering of a correlation matrix (cut at K) against the
# true block labels; NA correlations count as 0.
mat_ari <- function(cor_mat, true_labels, K) {
  x <- as.matrix(cor_mat)
  x[is.na(x)] <- 0
  mclust::adjustedRandIndex(true_labels, cutree(hclust(dist(x)), k = K))
}
