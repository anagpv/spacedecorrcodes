# run_pairwise.R - Scenario I: the correlation of one gene pair on one tissue,
# every method, under spatial covariance (Additive / MV), an FOV technical
# factor and an unobserved subtype. Writes pairwise_results.csv, read by
# figure_pairwise.R (main figure panels A and B; supplementary pairwise figure).
# Resumable: settings already in the output file are skipped.
#   N_CORES=<n>  workers (default: all cores)   DRY_RUN=1  list the tasks and stop
#   MAX_TASKS=<n>  run the first n tasks only

# Directories ----
dir.proj <- "/path/to/repo/"
dir.data <- "/path/to/zenodo/download/"
dir.out  <- "/path/to/output/"

# Settings ----
methods <- c("spacedecorr", "spMOCA", "logpearson", "spark", "meringue", "giotto", "GraphR")
seeds   <- 1:200
n_cores <- as.integer(Sys.getenv("N_CORES", parallel::detectCores()))
out_csv <- file.path(dir.out, "pairwise_results.csv")

source(file.path(dir.proj, "util", "simulator.R"))
source(file.path(dir.proj, "util", "compare.R"))
source(file.path(dir.proj, "util", "util.R"))
dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

# Simulation grid ----
rho_grid <- c(0, 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.90)

settings <- expand.grid(
  rept          = seeds,
  sim           = c("Additive", "MV", "FOV_diff", "subtype_v2"),
  num_cells     = c(500, 1000, 2000, 5000),
  range         = 0.25,
  m             = 0.5,
  ratio         = 0.5,
  libsize       = 200,
  rho           = rho_grid,
  pi_zero       = c(0, 0.3),               # zero mass
  subtype_rho   = 0,                       # subtype_v2: gene 1-2 coupling inside the subtype
  subtype_shift = c(0.05, 0.1, 0.2),       # subtype_v2: shift on the Gamma rate
  fov_diff_lo   = 0,                       # FOV_diff: per-gene shift ~ U(fov_diff_lo, fov_diff_hi)
  fov_diff_hi   = 1,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
settings <- subset(settings, num_cells == 2000 | rho == 0)                      # rho swept at 2000 cells
settings <- subset(settings, pi_zero == 0 | (sim %in% c("Additive", "MV") & num_cells == 2000))
settings <- unique(mask_unused_knobs(settings, PAIRWISE_KNOB_USERS))

# SpaceDecorr spline-basis dimension: 200 everywhere; also 10% of the cells and
# 1000 for Additive / MV without zero mass (the supplementary pairwise figure).
basis_dims <- function(setting) {
  k <- 200
  if (setting$sim %in% c("Additive", "MV") && setting$pi_zero == 0)
    k <- c(k, round(0.1 * setting$num_cells), 1000)
  unique(k[k <= setting$num_cells])
}

# Output ----
meta_cols   <- c("seed", "sim", "num_cells", "num_genes", "range", "m", "ratio", "v", "libsize",
                 "pi_zero", "rho", "subtype_rho", "subtype_shift", "fov_diff_lo", "fov_diff_hi")
result_cols <- c(meta_cols, "k", "ext_info", compare_result_cols)

setting_cols_csv  <- setdiff(meta_cols, c("num_genes", "v"))
setting_cols_grid <- sub("^seed$", "rept", setting_cols_csv)

n_all    <- nrow(settings)
settings <- filter_pending_settings(settings, out_csv, setting_cols_csv, setting_cols_grid, methods)
append_rows <- make_appender(out_csv, result_cols)

# One task: one simulated data set, every method ----
run_setting <- function(setting) {
  sim <- sim_pairwise(
    rept = setting$rept, sim = setting$sim, num_cells = setting$num_cells, range = setting$range,
    rho = setting$rho, m = setting$m, ratio = setting$ratio, libsize = setting$libsize,
    fov_diff_range = c(setting$fov_diff_lo, setting$fov_diff_hi),
    subtype_rho = setting$subtype_rho, subtype_shift = setting$subtype_shift,
    pi_zero = setting$pi_zero)
  stopifnot(setting_key(sim$meta, setting_cols_csv) == setting_key(setting, setting_cols_grid))
  base <- cbind(sim$meta, ext_info = encode_ext_info(sim$ext_info))

  rows <- lapply(methods, function(model) {
    ks <- if (model == "spacedecorr") basis_dims(setting) else NA_real_
    do.call(rbind, lapply(ks, function(k) {
      args <- list(model, sim$Y, sim$coordinates, sim$data, features = c(1, 2), membership = sim$membership)
      if (!is.na(k)) args$k <- k
      r <- do.call(run_compare_method, args)
      if (!isTRUE(r$status == "ok"))
        stop(sprintf("[sim=%s seed=%s num_cells=%s rho=%s pi_zero=%s] %s failed: %s",
                     setting$sim, setting$rept, setting$num_cells, setting$rho, setting$pi_zero,
                     model, r$error_message))
      cbind(base, k = k, r)
    }))
  })
  out <- do.call(rbind, rows)
  append_rows(out)
  message(sprintf("done  sim=%s seed=%s num_cells=%s rho=%s pi_zero=%s  (%d rows)",
                  setting$sim, setting$rept, setting$num_cells, setting$rho, setting$pi_zero, nrow(out)))
  TRUE
}

# Run ----
main <- function() {
  validate_or_init_csv(out_csv, result_cols)
  tasks <- settings[order(settings$rept), ]
  rownames(tasks) <- NULL
  max_tasks <- suppressWarnings(as.integer(Sys.getenv("MAX_TASKS", NA)))
  if (!is.na(max_tasks)) tasks <- head(tasks, max_tasks)

  cat(sprintf("pairwise | %d settings, %d done, %d to run x %d methods | %d cores\n%s\n",
              n_all, n_all - nrow(settings), nrow(tasks), length(methods), n_cores, out_csv))
  if (nrow(tasks) == 0L) return(invisible())
  if (nzchar(Sys.getenv("DRY_RUN"))) { print(table(sim = tasks$sim, pi_zero = tasks$pi_zero)); return(invisible()) }

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
