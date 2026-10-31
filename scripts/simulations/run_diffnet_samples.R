# run_diffnet_samples.R - differential network between two simulated tissues
# (20 genes, small-world networks differing in a quarter of the pairs), every
# method fitted on each tissue. Writes diffnet_samples_results.csv, read by
# figure_diffnet.R (supplementary differential-network figure, panels D and E).
# Resumable: settings already in the output file are skipped.
#   N_CORES=<n>  workers (default: all cores)   DRY_RUN=1  list the tasks and stop
#   MAX_TASKS=<n>  run the first n tasks only

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

# Settings ----
methods <- c("spacedecorr", "spMOCA", "logpearson", "spark", "meringue", "giotto")
seeds   <- 1:50
n_cores <- as.integer(Sys.getenv("N_CORES", parallel::detectCores()))
out_csv <- file.path(dir.out, "diffnet_samples_results.csv")
meringue_filter_dist <- 1                   # neighbours within 1 on the [-1, 1]^2 tissue

source(file.path(dir.proj, "util", "simulator_mat.R"))
source(file.path(dir.proj, "util", "compare_mat.R"))
source(file.path(dir.proj, "util", "util.R"))
dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

# Simulation grid ----
settings <- expand.grid(
  rept         = seeds,
  sim          = c("sample_MV", "sample_Additive"),
  num_genes    = 20,
  num_cells1   = 2000,
  num_cells2   = 2000,
  range        = 0.25,
  perturb_prob = 0.25,                     # fraction of the gene pairs that differ between the tissues
  nei          = 4,                        # small-world network: neighbours and rewiring probability
  prob         = 0.1,
  k            = 200,                      # SpaceDecorr spline-basis dimension
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

method_args <- function(model)
  switch(model,
         spacedecorr = list(k = settings$k[1]),
         meringue    = list(filterDist = meringue_filter_dist),
         list())

# Output ----
meta_cols   <- c("seed", "sim", "num_genes", "num_cells1", "num_cells2", "range",
                 "perturb_prob", "nei", "prob", "k")
result_cols <- c(meta_cols, "ext_info", compare_mat_result_cols)

setting_cols_csv  <- meta_cols
setting_cols_grid <- sub("^seed$", "rept", setting_cols_csv)

n_all    <- nrow(settings)
settings <- filter_pending_settings(settings, out_csv, setting_cols_csv, setting_cols_grid, methods)
append_rows <- make_appender(out_csv, result_cols)

# One task: one pair of tissues, every method on both ----
run_setting <- function(setting) {
  sim <- sim_samples(rept = setting$rept, sim = setting$sim, num_genes = setting$num_genes,
                     num_cells1 = setting$num_cells1, num_cells2 = setting$num_cells2,
                     range = setting$range, perturb_prob = setting$perturb_prob,
                     nei = setting$nei, prob = setting$prob)
  sim$meta$k <- setting$k
  stopifnot(setting_key(sim$meta, setting_cols_csv) == setting_key(setting, setting_cols_grid))
  base <- cbind(sim$meta, ext_info = encode_ext_info(sim$ext_info))

  rows <- lapply(methods, function(model) {
    fit <- function(s) do.call(compare_mat_methods[[model]], c(list(s$Y, s$coordinates, s$data), method_args(model)))
    r1 <- fit(sim$sample1)                          # tissue 1 -> group 0
    r2 <- fit(sim$sample2)                          # tissue 2 -> group 1
    res <- list(cor = NULL, cor_g0 = r1$cor, cor_g1 = r2$cor,
                n = NA_integer_, n_g0 = r1$n, n_g1 = r2$n, time = r1$time + r2$time)
    mat_pair_rows(model, res, sim$truth)
  })
  rows <- do.call(rbind, rows)
  stopifnot(nrow(rows) == length(methods) * choose(setting$num_genes, 2))
  append_rows(cbind(base[rep(1, nrow(rows)), ], rows))
  message(sprintf("done  sim=%s seed=%s  (%d rows)", setting$sim, setting$rept, nrow(rows)))
  TRUE
}

# Run ----
main <- function() {
  validate_or_init_csv(out_csv, result_cols)
  tasks <- settings
  rownames(tasks) <- NULL
  max_tasks <- suppressWarnings(as.integer(Sys.getenv("MAX_TASKS", NA)))
  if (!is.na(max_tasks)) tasks <- head(tasks, max_tasks)

  cat(sprintf("diffnet_samples | %d settings, %d done, %d to run x %d methods | %d cores\n%s\n",
              n_all, n_all - nrow(settings), nrow(tasks), length(methods), n_cores, out_csv))
  if (nrow(tasks) == 0L) return(invisible())
  if (nzchar(Sys.getenv("DRY_RUN"))) { print(table(sim = tasks$sim)); return(invisible()) }

  warm_mat_cache()                                  # the real cells' library sizes, loaded once

  run_setting(tasks[1, ])                           # the first task unforked: a systematic error stops here
  if (nrow(tasks) > 1L) {
    res <- parallel::mclapply(2:nrow(tasks), function(i) run_setting(tasks[i, ]),
                              mc.cores = n_cores, mc.preschedule = FALSE)
    failed <- which(!vapply(res, isTRUE, logical(1)))
    if (length(failed))
      stop(sprintf("%d task(s) failed; first: %s", length(failed),
                   paste(deparse(res[[failed[1]]]), collapse = " ")))
  }
  cat(sprintf("done - %d tasks\n", nrow(tasks)))
}

main()
