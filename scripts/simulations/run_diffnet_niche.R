# run_diffnet_niche.R - differential network between two real niches: the
# Lung5-3 macrophages of the stroma (niche 0) and the myeloid-enriched stroma
# (niche 1), 20 genes whose small-world networks differ in a quarter of the
# pairs; every method fitted within each niche. Writes diffnet_niche_results.csv,
# read by figure_diffnet.R (supplementary differential-network figure, panels B
# and C). Resumable: settings already in the output file are skipped.
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
out_csv <- file.path(dir.out, "diffnet_niche_results.csv")
meringue_filter_dist <- 1                   # neighbours within 1 mm

source(file.path(dir.proj, "util", "simulator_mat.R"))
source(file.path(dir.proj, "util", "compare_mat.R"))
source(file.path(dir.proj, "util", "util.R"))
dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

# Simulation grid ----
settings <- expand.grid(
  rept          = seeds,
  sim           = "niche_mat",
  num_cells     = nrow(mat_lung()$coordinates),
  num_genes     = 20,
  perturb_prob  = 0.25,                    # fraction of the gene pairs that differ between the niches
  nei           = 4,                       # small-world network: neighbours and rewiring probability
  prob          = 0.1,
  libsize_scale = 0.1,                     # library size = round(cell area x libsize_scale)
  k             = 200,                     # SpaceDecorr spline-basis dimension
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# Methods that build a spatial graph are fitted within each niche only.
method_args <- function(model)
  switch(model,
         spacedecorr = list(k = settings$k[1]),
         meringue    = list(filterDist = meringue_filter_dist, global = FALSE),
         giotto      = list(global = FALSE),
         spMOCA      = list(global = FALSE),
         list())

# Output ----
meta_cols   <- c("seed", "sim", "num_cells", "num_genes", "perturb_prob", "nei", "prob",
                 "libsize_scale", "n0", "n1", "k")
result_cols <- c(meta_cols, "ext_info", compare_mat_result_cols)

setting_cols_csv  <- setdiff(meta_cols, c("n0", "n1"))
setting_cols_grid <- sub("^seed$", "rept", setting_cols_csv)

n_all    <- nrow(settings)
settings <- filter_pending_settings(settings, out_csv, setting_cols_csv, setting_cols_grid, methods)
append_rows <- make_appender(out_csv, result_cols)

# One task: one simulated data set, every method ----
run_setting <- function(setting) {
  sim <- sim_niche_mat(rept = setting$rept, num_genes = setting$num_genes,
                       perturb_prob = setting$perturb_prob, nei = setting$nei, prob = setting$prob,
                       libsize_scale = setting$libsize_scale)
  sim$meta$k <- setting$k
  stopifnot(setting_key(sim$meta, setting_cols_csv) == setting_key(setting, setting_cols_grid))
  base <- cbind(sim$meta, ext_info = encode_ext_info(sim$ext_info))

  rows <- lapply(methods, function(model) {
    res <- do.call(compare_mat_methods[[model]],
                   c(list(sim$Y, sim$coordinates, sim$data, membership = sim$membership), method_args(model)))
    mat_pair_rows(model, res, sim$truth)
  })
  rows <- do.call(rbind, rows)
  stopifnot(nrow(rows) == length(methods) * choose(setting$num_genes, 2))
  append_rows(cbind(base[rep(1, nrow(rows)), ], rows))
  message(sprintf("done  seed=%s  (%d rows)", setting$rept, nrow(rows)))
  TRUE
}

# Run ----
main <- function() {
  validate_or_init_csv(out_csv, result_cols)
  tasks <- settings
  rownames(tasks) <- NULL
  max_tasks <- suppressWarnings(as.integer(Sys.getenv("MAX_TASKS", NA)))
  if (!is.na(max_tasks)) tasks <- head(tasks, max_tasks)

  cat(sprintf("diffnet_niche | %d settings, %d done, %d to run x %d methods | %d cores\n%s\n",
              n_all, n_all - nrow(settings), nrow(tasks), length(methods), n_cores, out_csv))
  if (nrow(tasks) == 0L) return(invisible())
  if (nzchar(Sys.getenv("DRY_RUN"))) { print(table(seed = tasks$rept)); return(invisible()) }

  warm_mat_cache(niche_fields = TRUE)               # the tissue and its Matern factors, loaded once

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
