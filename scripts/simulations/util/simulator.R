# simulator.R - count simulators of Scenario I (sim_pairwise) and Scenario II
# (sim_niche): Matern spatial field + gene correlation -> Gamma rate -> Poisson.
library(tidyverse)

# Gene correlation matrix: identity with genes 1 and 2 coupled by rho.
.delta_mat <- function(rho, num_genes = 3) {
  d <- diag(num_genes)
  d[1, 2] <- d[2, 1] <- rho
  d
}

# Factor a correlation matrix once (half %*% t(half) == sigma), draw fields from it.
.field_factor <- function(sigma) {
  eig <- eigen(sigma, symmetric = TRUE)
  pos <- eig$values > 0
  eig$vectors[, pos] %*% diag(sqrt(eig$values[pos]))
}

.draw_field <- function(half, n_fields)
  half %*% matrix(rnorm(ncol(half) * n_fields), nrow = ncol(half))

.spatial_field <- function(sigma, n_fields) .draw_field(.field_factor(sigma), n_fields)

# A knob a simulation does not read is recorded as NA (in `meta` and in the grid).
SUBTYPE_SIMS <- "subtype_v2"
PAIRWISE_KNOB_USERS <- list(subtype_rho = SUBTYPE_SIMS, subtype_shift = SUBTYPE_SIMS,
                            fov_diff_lo = "FOV_diff", fov_diff_hi = "FOV_diff")
NICHE_KNOB_USERS    <- list(radius = "niche")

mask_unused_knobs <- function(df, users) {
  for (k in names(users)) df[[k]][!df$sim %in% users[[k]]] <- NA
  df
}

# Latent Gaussian, cells x genes, marginally N(0, 1).
#   "MV":      matrix-normal draw A Z t(B), row covariance sigma, column covariance delta.
#   otherwise: gene-correlated draw (columns ~ delta) plus one spatial field per gene,
#              (gene + spatial) / sqrt(2). `delta` may be list(delta0, delta1) with a
#              0/1 `group` per cell; the spatial field is then one joint draw.
.latent_gaussian <- function(sim, sigma, delta, num_genes = 3, spatial = NULL, group = NULL) {
  n <- if (is.null(sigma)) nrow(spatial) else nrow(sigma)

  if (sim == "MV") {
    if (is.list(delta)) stop("sim = 'MV' takes a single delta matrix")
    eig_sigma <- eigen(sigma, symmetric = TRUE)
    eig_delta <- eigen(delta, symmetric = TRUE)
    r <- sum(eig_sigma$values > 0)
    s <- sum(eig_delta$values > 0)
    A <- eig_sigma$vectors[, 1:r] %*% diag(sqrt(eig_sigma$values[1:r]))
    B <- eig_delta$vectors[, 1:s] %*% diag(sqrt(eig_delta$values[1:s]))
    Z <- matrix(rnorm(r * s), nrow = r, ncol = s)
    return(A %*% Z %*% t(B))
  }

  if (is.list(delta)) {
    stopifnot(!is.null(group), length(group) == n)
    gene_part <- matrix(0, n, num_genes)
    i0 <- which(group == 0)
    i1 <- which(group == 1)
    if (length(i0)) gene_part[i0, ] <- MASS::mvrnorm(length(i0), rep(0, num_genes), delta[[1]])
    if (length(i1)) gene_part[i1, ] <- MASS::mvrnorm(length(i1), rep(0, num_genes), delta[[2]])
  } else {
    gene_part <- MASS::mvrnorm(n, rep(0, num_genes), delta)
  }
  if (is.null(spatial)) spatial <- .spatial_field(sigma, num_genes)
  (gene_part + spatial) / sqrt(2)
}

# Latent Gaussian -> Gamma rate multipliers (+ an optional additive `shift`) ->
# Poisson(libsize * rate) -> a pi_zero fraction set to structural zeros.
.sim_counts <- function(sim, sigma, delta, m, v, libsize, pi_zero,
                        num_genes = 3, spatial = NULL, shift = NULL, group = NULL) {
  latent <- .latent_gaussian(sim, sigma, delta, num_genes, spatial, group)
  gamma  <- qgamma(pnorm(latent), shape = (m^2) / v, rate = m / v)
  if (!is.null(shift))
    gamma <- if (is.matrix(shift)) gamma + shift else sweep(gamma, 1, shift, "+")

  n <- nrow(latent)
  Y <- matrix(rpois(n * num_genes, lambda = libsize * gamma), nrow = n, ncol = num_genes)
  zero_mask <- matrix(rbinom(n * num_genes, size = 1, prob = pi_zero) == 1, nrow = n, ncol = num_genes)
  Y[zero_mask] <- 0
  Y
}

# Scenario I ----
# One tissue of num_cells cells on [-1, 1]^2, three genes, genes 1 and 2 coupled by rho.
#   "MV", "Additive"  the two latent models
#   "FOV_diff"        Additive plus a per-gene constant ~ U(fov_diff_range) added to the
#                     rate inside the square {0 < x < 0.5, -0.5 < y < 0}
#   "subtype_v2"      Additive with an unobserved subpopulation (cells where a smooth
#                     spatial field is below its median): genes 1-2 coupled by
#                     subtype_rho there; gene 1 up by subtype_shift in the subtype,
#                     genes 2-3 up by subtype_shift elsewhere
# Returns list(Y, coordinates, data, membership, meta, ext_info).
sim_pairwise <- function(rept, sim, num_cells, range, rho, m, ratio, libsize,
                         fov_diff_range = c(0, 1), subtype_rho = 0, subtype_shift = 0.05,
                         pi_zero = 0, smoothness = 0.5) {
  set.seed(47289 * rept)
  v <- m / ratio
  num_genes <- 3

  coordinates <- data.frame(sdimx = runif(num_cells, -1, 1),
                            sdimy = runif(num_cells, -1, 1))
  sigma    <- fields::Matern(as.matrix(dist(coordinates)), smoothness = smoothness, range = range)
  delta    <- .delta_mat(rho, num_genes)
  membership <- NULL
  ext_info   <- list()

  if (sim %in% SUBTYPE_SIMS) {
    fields_drawn <- .spatial_field(sigma, num_genes + 1)      # 3 expression fields + 1 membership field
    is_sub <- pnorm(fields_drawn[, num_genes + 1]) < 0.5
    stopifnot(sum(is_sub) > 3, sum(!is_sub) > 3)

    shift <- matrix(0, num_cells, num_genes)
    shift[is_sub,  1]   <- subtype_shift
    shift[!is_sub, 2:3] <- subtype_shift

    Y <- .sim_counts("Additive", sigma, list(delta, .delta_mat(subtype_rho, num_genes)),
                     m, v, libsize, pi_zero, num_genes,
                     spatial = fields_drawn[, seq_len(num_genes), drop = FALSE],
                     group = as.integer(is_sub), shift = shift)
    membership <- as.integer(is_sub)
    ext_info   <- list(prop_membership1 = mean(is_sub))
  } else {
    shift <- NULL
    if (sim == "FOV_diff") {
      in_fov     <- with(coordinates, sdimx > 0 & sdimx < 0.5 & sdimy > -0.5 & sdimy < 0)
      gene_shift <- runif(num_genes, fov_diff_range[1], fov_diff_range[2])
      shift      <- outer(as.numeric(in_fov), gene_shift)
      ext_info   <- list(gene_shift = gene_shift)
    }
    Y <- .sim_counts(sim, sigma, delta, m, v, libsize, pi_zero, num_genes = num_genes, shift = shift)
  }

  Y <- as.data.frame(Y)
  colnames(Y) <- paste0("Y", 1:num_genes)
  rownames(Y) <- 1:num_cells
  data <- Y %>% cbind(coordinates) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize

  meta <- data.frame(seed = rept, sim = sim, num_cells = num_cells, num_genes = num_genes,
                     range = range, m = m, ratio = ratio, v = v, libsize = libsize,
                     pi_zero = pi_zero, rho = rho,
                     subtype_rho = subtype_rho, subtype_shift = subtype_shift,
                     fov_diff_lo = fov_diff_range[1], fov_diff_hi = fov_diff_range[2],
                     stringsAsFactors = FALSE)
  meta <- mask_unused_knobs(meta, PAIRWISE_KNOB_USERS)

  list(Y = Y, coordinates = coordinates, data = data,
       membership = membership, meta = meta, ext_info = ext_info)
}

# Scenario II, real tissue geometry ----
# The Lung5-3 macrophages in the "myeloid-enriched stroma" and "stroma" niches,
# native coordinates; "stroma" is the niche (membership 1). Loaded once; the
# Matern factors of the geometry are cached, since it is the same for every seed.
NICHE_REAL_RDATA  <- file.path(dir.data, "Lung5-3data.Rdata")
.niche_real_cache <- new.env(parent = emptyenv())

niche_real_geometry <- function(path = NICHE_REAL_RDATA) {
  if (is.null(.niche_real_cache$geometry)) {
    e <- new.env()
    load(path, envir = e)
    keep <- which(e$fulldata$cell_type == "macrophage" &
                  e$fulldata$niche %in% c("myeloid-enriched stroma", "stroma"))
    coordinates <- data.frame(sdimx = e$fulldata$sdimx[keep], sdimy = e$fulldata$sdimy[keep])
    niche       <- as.integer(e$fulldata$niche[keep] == "stroma")
    rm(e)

    # the FOV square of the simulated box, mapped onto the tissue's bounding rectangle
    cx <- mean(range(coordinates$sdimx)); hx <- diff(range(coordinates$sdimx)) / 2
    cy <- mean(range(coordinates$sdimy)); hy <- diff(range(coordinates$sdimy)) / 2
    .niche_real_cache$geometry <- list(
      coordinates = coordinates, niche = niche,
      fov_box = c(xlo = cx, xhi = cx + 0.5 * hx, ylo = cy - 0.5 * hy, yhi = cy))
  }
  .niche_real_cache$geometry
}

.niche_real_factor <- function(cells = c("all", "out", "in"), field_range) {
  cells <- match.arg(cells)
  key   <- paste(cells, field_range, sep = "|")
  if (is.null(.niche_real_cache[[key]])) {
    g   <- niche_real_geometry()
    idx <- switch(cells, all = seq_along(g$niche), out = which(g$niche == 0), `in` = which(g$niche == 1))
    message(sprintf("niche_real: factoring the Matern covariance of %d cells (%s, range %s)",
                    length(idx), cells, field_range))
    .niche_real_cache[[key]] <- .field_factor(
      fields::Matern(as.matrix(dist(g$coordinates[idx, ])), smoothness = 0.5, range = field_range))
  }
  .niche_real_cache[[key]]
}

# Fill the cache before forking so the workers share the factors.
warm_niche_real_cache <- function(ranges, range_niches = NA_real_) {
  split_ranges <- unique(range_niches[!is.na(range_niches)])
  for (r in unique(ranges)) {
    if (anyNA(range_niches))  .niche_real_factor("all", r)
    if (length(split_ranges)) .niche_real_factor("out", r)
  }
  for (r in split_ranges) .niche_real_factor("in", r)
  invisible(NULL)
}

# Scenario II ----
# Additive model; genes 1 and 2 coupled by rho outside the niche and rho_niche inside.
#   "niche"       num_cells cells on [-1, 1]^2, the niche a disc of radius `radius`
#   "niche_real"  the real tissue geometry above (num_cells and radius are not read)
#   range_niche   NULL: one spatial field over all cells; a number: independent
#                 inside / outside fields, the inside one at this range
#   elevate_effect  when rho_niche != 0, genes 1-2 are also raised inside the niche
#                   by a gene-specific constant ~ U(elevate_effect)
#   fov_effect    per-gene constant ~ U(fov_effect) on all genes inside the FOV square
# Returns list(Y, coordinates, data, membership, meta, ext_info).
sim_niche <- function(rept, sim = c("niche", "niche_real"), num_cells, range, rho, rho_niche,
                      m, ratio, libsize, radius = 0.35, niche_center = c(0, 0),
                      range_niche = NULL, fov_effect = c(0, 0), elevate_effect = c(0, 1),
                      pi_zero = 0) {
  sim <- match.arg(sim)
  set.seed(47289 * rept)
  v <- m / ratio
  num_genes <- 3

  if (sim == "niche") {
    coordinates <- data.frame(sdimx = runif(num_cells, -1, 1),
                              sdimy = runif(num_cells, -1, 1))
    d       <- sqrt((coordinates$sdimx - niche_center[1])^2 + (coordinates$sdimy - niche_center[2])^2)
    niche   <- as.integer(d <= radius)
    fov_box <- c(xlo = 0, xhi = 0.5, ylo = -0.5, yhi = 0)
  } else {
    geometry    <- niche_real_geometry()
    coordinates <- geometry$coordinates
    niche       <- geometry$niche
    fov_box     <- geometry$fov_box
    num_cells   <- nrow(coordinates)
  }
  id0 <- which(niche == 0)
  id1 <- which(niche == 1)
  stopifnot(length(id0) > 3, length(id1) > 3)

  delta0 <- .delta_mat(rho,       num_genes)
  delta1 <- .delta_mat(rho_niche, num_genes)

  shift_mat <- matrix(0, num_cells, num_genes)
  elev_gene_shift <- c(0, 0)
  if (rho_niche != 0) {
    elev_gene_shift <- runif(2, elevate_effect[1], elevate_effect[2])
    shift_mat[id1, 1] <- shift_mat[id1, 1] + elev_gene_shift[1]
    shift_mat[id1, 2] <- shift_mat[id1, 2] + elev_gene_shift[2]
  }
  fov_gene_shift <- rep(0, num_genes)
  if (fov_effect[2] > fov_effect[1]) {
    in_fov <- with(coordinates, sdimx > fov_box[["xlo"]] & sdimx < fov_box[["xhi"]] &
                     sdimy > fov_box[["ylo"]] & sdimy < fov_box[["yhi"]])
    fov_gene_shift <- runif(num_genes, fov_effect[1], fov_effect[2])
    shift_mat[in_fov, ] <- sweep(shift_mat[in_fov, , drop = FALSE], 2, fov_gene_shift, "+")
  }
  shift_use <- if (any(shift_mat != 0)) shift_mat else NULL

  draw_counts <- function(idx, cells, field_range, delta, group = NULL) {
    sigma <- spatial <- NULL
    if (sim == "niche") {
      sigma <- fields::Matern(as.matrix(dist(coordinates[idx, ])), smoothness = 0.5, range = field_range)
    } else {
      spatial <- .draw_field(.niche_real_factor(cells, field_range), num_genes)
    }
    .sim_counts("Additive", sigma, delta, m, v, libsize, pi_zero, num_genes,
                spatial = spatial, group = group,
                shift = if (is.null(shift_use)) NULL else shift_use[idx, , drop = FALSE])
  }

  if (is.null(range_niche)) {
    Y <- draw_counts(seq_len(num_cells), "all", range, list(delta0, delta1), group = niche)
  } else {
    Y <- matrix(NA_real_, nrow = num_cells, ncol = num_genes)
    Y[id0, ] <- draw_counts(id0, "out", range,       delta0)
    Y[id1, ] <- draw_counts(id1, "in",  range_niche, delta1)
  }

  Y <- as.data.frame(Y)
  colnames(Y) <- paste0("Y", 1:num_genes)
  rownames(Y) <- 1:num_cells
  data <- Y %>% cbind(coordinates) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize

  meta <- data.frame(seed = rept, sim = sim, num_cells = num_cells, num_genes = num_genes,
                     range = range, m = m, ratio = ratio, v = v, libsize = libsize,
                     pi_zero = pi_zero, rho = rho, rho_niche = rho_niche, radius = radius,
                     range_niche = if (is.null(range_niche)) NA_real_ else range_niche,
                     fov_lo = fov_effect[1], fov_hi = fov_effect[2],
                     elevate_lo = elevate_effect[1], elevate_hi = elevate_effect[2],
                     stringsAsFactors = FALSE)
  meta <- mask_unused_knobs(meta, NICHE_KNOB_USERS)
  ext_info <- list(elev_gene_shift = elev_gene_shift, fov_gene_shift = fov_gene_shift,
                   prop_membership1 = mean(niche))

  list(Y = Y, coordinates = coordinates, data = data,
       membership = niche, meta = meta, ext_info = ext_info)
}
