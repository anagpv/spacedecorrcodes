# run_niche.R - Scenario II: a gene pair coupled inside a spatial niche only,
# every method, on a simulated disc niche ("niche") and on the real tissue
# geometry ("niche_real"). Each method is also run on the inert pair (genes 1
# and 3; columns suffixed _null). Writes niche_results.csv, read by
# figure_niche.R (main figure panel C; supplementary niche figures).
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
out_csv <- file.path(dir.out, "niche_results.csv")

source(file.path(dir.proj, "util", "simulator.R"))
source(file.path(dir.proj, "util", "compare.R"))
source(file.path(dir.proj, "util", "util.R"))
dir.create(dir.out, recursive = TRUE, showWarnings = FALSE)

# Simulation grid ----
rho_grid <- c(0, 0.01, 0.025, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.90)

niche_grid <- function(sim, num_cells, radius) expand.grid(
  rept        = seeds,
  sim         = sim,
  num_cells   = num_cells,
  range       = c(0.25, 0.8),
  m           = 0.5,
  ratio       = 0.5,
  libsize     = 200,
  rho         = 0,                         # gene 1-2 coupling outside the niche
  rho_niche   = rho_grid,                  # gene 1-2 coupling inside the niche
  radius      = radius,
  range_niche = c(NA_real_, 0.1),          # NA: one shared spatial field; 0.1: separate field inside
  fov_lo      = 0,                         # FOV per-gene shift ~ U(fov_lo, fov_hi); 0, 0 = none
  fov_hi      = c(0, 1),
  elevate_lo  = 0,                         # niche co-elevation of genes 1-2 ~ U(elevate_lo, elevate_hi)
  elevate_hi  = 1,
  pi_zero     = 0,
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)

settings <- rbind(
  niche_grid("niche",      num_cells = 2000, radius = c(0.2, 0.4, 0.6, 0.8, 1.0)),
  niche_grid("niche_real", num_cells = nrow(niche_real_geometry()$coordinates), radius = NA_real_))
settings <- unique(mask_unused_knobs(settings, NICHE_KNOB_USERS))

# Output ----
meta_cols   <- c("seed", "sim", "num_cells", "num_genes", "range", "m", "ratio", "v", "libsize",
                 "pi_zero", "rho", "rho_niche", "radius", "range_niche", "fov_lo", "fov_hi",
                 "elevate_lo", "elevate_hi")
null_cols   <- paste0(setdiff(compare_result_cols, "model"), "_null")
result_cols <- c(meta_cols, "ext_info", compare_result_cols, null_cols)

setting_cols_csv  <- setdiff(meta_cols, c("num_genes", "v"))
setting_cols_grid <- sub("^seed$", "rept", setting_cols_csv)

n_all    <- nrow(settings)
settings <- filter_pending_settings(settings, out_csv, setting_cols_csv, setting_cols_grid, methods)
append_rows <- make_appender(out_csv, result_cols)

# One task: one simulated data set, every method on both gene pairs ----
run_setting <- function(setting) {
  sim <- sim_niche(
    rept = setting$rept, sim = setting$sim, num_cells = setting$num_cells, range = setting$range,
    rho = setting$rho, rho_niche = setting$rho_niche, m = setting$m, ratio = setting$ratio,
    libsize = setting$libsize, radius = setting$radius,
    range_niche = if (is.na(setting$range_niche)) NULL else setting$range_niche,
    fov_effect = c(setting$fov_lo, setting$fov_hi),
    elevate_effect = c(setting$elevate_lo, setting$elevate_hi), pi_zero = setting$pi_zero)
  stopifnot(setting_key(sim$meta, setting_cols_csv) == setting_key(setting, setting_cols_grid))
  base <- cbind(sim$meta, ext_info = encode_ext_info(sim$ext_info))
  tag  <- sprintf("[sim=%s seed=%s range=%s rho_niche=%s radius=%s range_niche=%s fov_hi=%s]",
                  setting$sim, setting$rept, setting$range, setting$rho_niche, setting$radius,
                  setting$range_niche, setting$fov_hi)

  rows <- lapply(methods, function(model) {
    run_pair <- function(features) {
      args <- list(model, sim$Y, sim$coordinates, sim$data, features = features, membership = sim$membership)
      if (model == "spacedecorr") args$k <- 200
      r <- do.call(run_compare_method, args)
      if (!isTRUE(r$status == "ok"))
        stop(sprintf("%s %s failed on genes %s: %s", tag, model, paste(features, collapse = ","), r$error_message))
      if (anyNA(c(r$cor_g0, r$cor_g1, r$cor_diff)))
        stop(sprintf("%s %s returned NA within-niche estimates on genes %s", tag, model, paste(features, collapse = ",")))
      r
    }
    coupled <- run_pair(c(1, 2))
    inert   <- run_pair(c(1, 3))[setdiff(compare_result_cols, "model")]
    names(inert) <- paste0(names(inert), "_null")
    cbind(base, coupled, inert)
  })
  out <- do.call(rbind, rows)
  append_rows(out)
  message(sprintf("done  %s  (%d rows)", tag, nrow(out)))
  TRUE
}

# Run ----
main <- function() {
  validate_or_init_csv(out_csv, result_cols)
  tasks <- settings[order(settings$sim != "niche_real", settings$rept), ]
  rownames(tasks) <- NULL
  max_tasks <- suppressWarnings(as.integer(Sys.getenv("MAX_TASKS", NA)))
  if (!is.na(max_tasks)) tasks <- head(tasks, max_tasks)

  cat(sprintf("niche | %d settings, %d done, %d to run x %d methods | %d cores\n%s\n",
              n_all, n_all - nrow(settings), nrow(tasks), length(methods), n_cores, out_csv))
  if (nrow(tasks) == 0L) return(invisible())
  if (nzchar(Sys.getenv("DRY_RUN"))) { print(table(sim = tasks$sim, range = tasks$range)); return(invisible()) }

  # the real geometry never changes: factor its Matern covariances once, before forking
  real <- tasks[tasks$sim == "niche_real", ]
  if (nrow(real)) warm_niche_real_cache(real$range, real$range_niche)

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
