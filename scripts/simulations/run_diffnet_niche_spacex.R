# run_diffnet_niche_spacex.R - SpaceX on the two-real-niche differential
# network of run_diffnet_niche.R (same simulation, same grid; one long fit per
# seed). Writes diffnet_niche_spacex_results.csv, read by figure_diffnet.R
# (supplementary differential-network figure, panels B and C).
# Needs SpaceX, MSFA and PQLseq (the PQLseq/ folder of this repository installs).
#   SPACEX_CORES=<n>  cores per SpaceX fit (default 20, one per gene)
#   N_TASKS=<n>       seeds fitted at the same time (default 1); uses N_TASKS x SPACEX_CORES cores
#   DRY_RUN=1  list the tasks and stop   MAX_TASKS=<n>  run the first n tasks only

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

# Settings ----
methods      <- "spacex"
seeds        <- 1:50
spacex_cores <- as.integer(Sys.getenv("SPACEX_CORES", "20"))
n_tasks      <- as.integer(Sys.getenv("N_TASKS", "1"))
out_csv      <- file.path(dir.out, "diffnet_niche_spacex_results.csv")

source(file.path(dir.proj, "util", "simulator_mat.R"))
source(file.path(dir.proj, "util", "compare_mat.R"))
source(file.path(dir.proj, "util", "util.R"))
patch_spacex()
dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

# Simulation grid (run_diffnet_niche.R's) ----
settings <- expand.grid(
  rept          = seeds,
  sim           = "niche_mat",
  num_cells     = nrow(mat_lung()$coordinates),
  num_genes     = 20,
  perturb_prob  = 0.25,
  nei           = 4,
  prob          = 0.1,
  libsize_scale = 0.1,
  k             = 200,                     # not read by SpaceX; kept so the files join
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

# Output ----
meta_cols   <- c("seed", "sim", "num_cells", "num_genes", "perturb_prob", "nei", "prob",
                 "libsize_scale", "n0", "n1", "k")
result_cols <- c(meta_cols, "ext_info", compare_mat_result_cols)

setting_cols_csv  <- setdiff(meta_cols, c("n0", "n1"))
setting_cols_grid <- sub("^seed$", "rept", setting_cols_csv)

n_all    <- nrow(settings)
settings <- filter_pending_settings(settings, out_csv, setting_cols_csv, setting_cols_grid, methods)
append_rows <- make_appender(out_csv, result_cols)

# One task: one simulated data set, SpaceX ----
run_setting <- function(setting) {
  sim <- sim_niche_mat(rept = setting$rept, num_genes = setting$num_genes,
                       perturb_prob = setting$perturb_prob, nei = setting$nei, prob = setting$prob,
                       libsize_scale = setting$libsize_scale)
  sim$meta$k <- setting$k
  stopifnot(setting_key(sim$meta, setting_cols_csv) == setting_key(setting, setting_cols_grid))
  base <- cbind(sim$meta, ext_info = encode_ext_info(sim$ext_info))

  res  <- compare_mat_spacex(sim$Y, sim$coordinates, sim$data, membership = sim$membership,
                             numCore = spacex_cores)
  rows <- mat_pair_rows("spacex", res, sim$truth)
  append_rows(cbind(base[rep(1, nrow(rows)), ], rows))
  message(sprintf("done  seed=%s  (%d rows, %.0f s)", setting$rept, nrow(rows), rows$time[1]))
  TRUE
}

# Run ----
main <- function() {
  validate_or_init_csv(out_csv, result_cols)
  tasks <- settings
  rownames(tasks) <- NULL
  max_tasks <- suppressWarnings(as.integer(Sys.getenv("MAX_TASKS", NA)))
  if (!is.na(max_tasks)) tasks <- head(tasks, max_tasks)

  cat(sprintf("diffnet_niche_spacex | %d settings, %d done, %d to run | %d task(s) x %d cores\n%s\n",
              n_all, n_all - nrow(settings), nrow(tasks), n_tasks, spacex_cores, out_csv))
  if (nrow(tasks) == 0L) return(invisible())
  if (nzchar(Sys.getenv("DRY_RUN"))) { print(table(seed = tasks$rept)); return(invisible()) }

  warm_mat_cache(niche_fields = TRUE)               # the tissue and its Matern factors, loaded once

  worker <- function(i) {
    withCallingHandlers(run_setting(tasks[i, ]),
                        error = function(e) message(sprintf("seed %s FAILED: %s", tasks$rept[i], conditionMessage(e))))
  }
  res <- if (n_tasks == 1L) lapply(seq_len(nrow(tasks)), worker)
         else parallel::mclapply(seq_len(nrow(tasks)), worker, mc.cores = n_tasks, mc.preschedule = FALSE)
  failed <- which(!vapply(res, isTRUE, logical(1)))
  if (length(failed)) stop(sprintf("%d task(s) failed; first: seed %s", length(failed), tasks$rept[failed[1]]))
  cat(sprintf("done - %d tasks\n", nrow(tasks)))
}

main()
