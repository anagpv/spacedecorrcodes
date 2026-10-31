# run_block.R - block co-expression structure: 100 genes in 4 correlation blocks
# on one tissue, every method's estimated correlation matrix tested pair by pair.
# Writes block_results.csv, read by figure_block.R (supplementary block figure).
# Resumable: settings already in the output file are skipped.
#   N_CORES=<n>  workers (default: all cores)   DRY_RUN=1  list the tasks and stop
#   MAX_TASKS=<n>  run the first n tasks only

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

# Settings ----
methods <- c("spacedecorr", "spMOCA", "logpearson", "spark", "meringue", "giotto")
seeds   <- 1:100
n_cores <- as.integer(Sys.getenv("N_CORES", parallel::detectCores()))
out_csv <- file.path(dir.out, "block_results.csv")
meringue_filter_dist <- 0.1                 # neighbours within 0.1 on the [-1, 1]^2 tissue

source(file.path(dir.proj, "util", "simulator_mat.R"))
source(file.path(dir.proj, "util", "compare_mat.R"))
source(file.path(dir.proj, "util", "util.R"))
dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

# Simulation grid ----
settings <- expand.grid(
  rept      = seeds,
  sim       = c("block_Additive", "block_MV"),
  num_cells = 2000,
  num_genes = 100,
  range     = 0.25,
  numblocks = 4,
  pi_zero   = c(0, 0.3),                   # zero mass
  k         = 200,                         # SpaceDecorr spline-basis dimension
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

method_args <- function(model)
  switch(model,
         spacedecorr = list(k = settings$k[1]),
         meringue    = list(filterDist = meringue_filter_dist),
         list())

# Output ----
meta_cols   <- c("seed", "sim", "num_cells", "num_genes", "range", "numblocks", "pi_zero", "k")
result_cols <- c(meta_cols, "model", "gene_i", "gene_j", "block_i", "block_j", "true_edge0",
                 "cor", "pvalue", "padj")

setting_cols_csv  <- meta_cols
setting_cols_grid <- sub("^seed$", "rept", setting_cols_csv)

n_all    <- nrow(settings)
settings <- filter_pending_settings(settings, out_csv, setting_cols_csv, setting_cols_grid, methods)
append_rows <- make_appender(out_csv, result_cols)

# One task: one simulated data set, every method ----
run_setting <- function(setting) {
  sim <- sim_block(rept = setting$rept, sim = setting$sim, num_cells = setting$num_cells,
                   num_genes = setting$num_genes, range = setting$range,
                   numblocks = setting$numblocks, pi_zero = setting$pi_zero)
  sim$meta$k <- setting$k
  stopifnot(setting_key(sim$meta, setting_cols_csv) == setting_key(setting, setting_cols_grid))

  rows <- lapply(methods, function(model) {
    res <- do.call(compare_mat_methods[[model]],
                   c(list(sim$Y, sim$coordinates, sim$data), method_args(model)))
    mat_pair_rows(model, res, sim$truth)
  })
  rows <- do.call(rbind, rows)
  stopifnot(nrow(rows) == length(methods) * choose(setting$num_genes, 2))
  append_rows(cbind(sim$meta[rep(1, nrow(rows)), ], rows))
  message(sprintf("done  sim=%s seed=%s pi_zero=%s  (%d rows)", setting$sim, setting$rept, setting$pi_zero, nrow(rows)))
  TRUE
}

# Run ----
main <- function() {
  validate_or_init_csv(out_csv, result_cols)
  tasks <- settings
  rownames(tasks) <- NULL
  max_tasks <- suppressWarnings(as.integer(Sys.getenv("MAX_TASKS", NA)))
  if (!is.na(max_tasks)) tasks <- head(tasks, max_tasks)

  cat(sprintf("block | %d settings, %d done, %d to run x %d methods | %d cores\n%s\n",
              n_all, n_all - nrow(settings), nrow(tasks), length(methods), n_cores, out_csv))
  if (nrow(tasks) == 0L) return(invisible())
  if (nzchar(Sys.getenv("DRY_RUN"))) { print(table(sim = tasks$sim, pi_zero = tasks$pi_zero)); return(invisible()) }

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
