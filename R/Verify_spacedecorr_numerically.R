#-----------------------------------------------------------------------------#
#- Numerical verification of SpaceDecorr against the equations in Methods     -#
#- Data generated as in the pairwise simulations (Additive, NB)               -#
#- Checks both fitting paths: gam (n < 5000) and bam (n >= 5000)              -#
#-----------------------------------------------------------------------------#
# remotes::install_github("anagpv/spacedecorr")
library(spacedecorr)
library(mgcv)
dir.data <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real Data/data/"
dir.out  <- "/Users/anagabipv/Library/Mobile Documents/com~apple~CloudDocs/Documents UW/Dissertation/SpaceDecorr project/Real Data/output/"


# Simulation --------------------------------------------------------------
# simulation parameters 
range   <- 0.25
rho     <- 0.3
libsize <- 200
k_used  <- 200

verify <- data.frame()
for (num_cells in c(2000, 5000)) {
  print(num_cells)
  set.seed(47289)
  
  # Additive framework: gene-level correlation + independent spatial field per gene
  sample_coords <- data.frame(sdimx = runif(num_cells, -1, 1), sdimy = runif(num_cells, -1, 1))
  sigma_r  <- fields::Matern(as.matrix(dist(sample_coords)), smoothness = 0.5, range = range)
  deltamat <- matrix(c(1, rho, rho, 1), ncol = 2)
  tmp <- MASS::mvrnorm(num_cells, rep(0, 2), deltamat)
  gamma <- sapply(1:2, function(x) {
    s_spatial <- as.numeric(scale(MASS::mvrnorm(n = 1, mu = rep(0, num_cells), Sigma = sigma_r)))
    z <- (as.numeric(scale(tmp[, x])) + s_spatial) / sqrt(2)
    qgamma(pnorm(z), shape = 0.25, rate = 2)
  })
  rm(sigma_r); gc()
  
  # NB counts
  Y <- matrix(rnbinom(num_cells * 2, mu = libsize * as.numeric(gamma), size = m * ratio), ncol = 2)
  colnames(Y) <- c("Y1", "Y2")
  rownames(Y) <- paste0("cell", 1:num_cells)
  data <- data.frame(sample_coords, libsize = libsize, row.names = rownames(Y))
  
  # SpaceDecorr (as called in the simulations)
  fit_pkg <- spacedecorr(assay_matrix = Y, metadata = data, basis = "ts",
                         libsize_col = "libsize", family = "nb", k = k_used)
  use_bam <- fit_pkg$config$use_bam
  
  # model as written in the Methods
  form <- as.formula(paste0("y ~ offset(log(libsize)) + s(sdimx, sdimy, bs = 'ts', k = ", k_used, ")"))
  
  for (g in colnames(Y)) {
    data$y <- Y[, g]
    fit <- if (use_bam) {
      bam(form, family = nb(), data = data, method = "fREML", discrete = TRUE, select = TRUE)
    } else {
      gam(form, family = nb(), data = data, method = "REML", select = TRUE)
    }
    mu    <- fitted(fit)
    theta <- fit$family$getTheta(TRUE)
    
    # Pearson residuals from the formula, clipped at +-sqrt(n)
    r_formula <- (data$y - mu) / sqrt(mu + mu^2 / theta)
    r_formula <- pmin(pmax(r_formula, -sqrt(num_cells)), sqrt(num_cells))
    s_formula <- mgcv::predict.gam(fit, data, type = "terms")[, "s(sdimx,sdimy)"]
    
    verify <- rbind(verify, data.frame(
      n = num_cells, fit = ifelse(use_bam, "bam", "gam"), gene = g,
      mean_count = round(mean(data$y), 2), theta_hat = round(theta, 2),
      max_diff_variance = max(abs(mu + mu^2 / theta - fit$family$variance(mu))),
      max_diff_resid    = max(abs(r_formula - fit_pkg$residuals[, g])),
      max_diff_trend    = max(abs(s_formula - fit_pkg$splines[, g])))
  }
}
rownames(verify) <- NULL
verify
write.csv(verify, file.path(dir.out, "verification_spacedecorr_simulation.csv"), row.names = FALSE)

# Liver example -----------------------------------------------------------
tissue      <- "AM048"
genes_check <- c("CYP2E1", "SDS", "CPS1", "ALDOB", "POR")

# raw counts (cells x genes) with cell IDs as rownames
meta_data  <- read.csv(file.path(dir.data, "Liver/cell_properties_healthy_diseased_merfish.csv"))
raw        <- read.csv(file.path(dir.data, "Liver/X_raw_healthy_diseased_merfish.csv"), header = FALSE)
gene_names <- read.csv(file.path(dir.data, "Liver/gene_names_healthy_diseased_merfish.csv"), header = FALSE)
colnames(raw) <- gene_names$V1
rownames(raw) <- meta_data$index

# SpaceDecorr fit: residuals, model data and configuration
f <- file.path(dir.out,"Liver/", tissue, paste0("Data_", tissue, "_spacedecorr_k200_hepfibcov.Rdata"))
load(f)
fit_pkg <- data_spacedecorr_hepcov
cfg     <- fit_pkg$config; cfg
dat     <- fit_pkg$metadata                      # sdimx, sdimy, libsize, covariates (model data)
res_pkg <- as.matrix(fit_pkg$residuals); rownames(res_pkg) <- rownames(dat)
spl_pkg <- as.matrix(fit_pkg$splines);   rownames(spl_pkg) <- rownames(dat)
n <- nrow(dat)
table(dat$cell_type_cov)

# model as written in the Methods
form <- as.formula(paste0("y ~ offset(log(libsize)) + ",
                          paste(c(cfg$covariates, paste0("s(sdimx, sdimy, bs = '", cfg$basis, "', k = ", cfg$k, ")")),
                                collapse = " + ")))
form

verify_liver <- data.frame()
for (g in genes_check) {
  print(g)
  dat$y <- raw[rownames(dat), g]
  
  fit <- if (cfg$use_bam) {
    bam(form, family = nb(), data = dat, method = "fREML", discrete = TRUE, select = TRUE)
  } else {
    gam(form, family = nb(), data = dat, method = "REML", select = TRUE)
  }
  mu    <- fitted(fit)
  theta <- fit$family$getTheta(TRUE)
  
  # NB parameterization: Var = mu + mu^2 / theta
  v_formula  <- mu + mu^2 / theta
  ll_formula <- sum(dnbinom(dat$y, size = theta, mu = mu, log = TRUE))
  
  # Pearson residuals from the formula, clipped at +-sqrt(n)
  r_formula <- (dat$y - mu) / sqrt(v_formula)
  r_formula <- pmin(pmax(r_formula, -sqrt(n)), sqrt(n))
  
  # fitted spatial trend
  s_fit <- mgcv::predict.gam(fit, dat, type = "terms")[, "s(sdimx,sdimy)"]
  
  verify_liver <- bind_rows(verify_liver, data.frame(
    gene = g, mean_count = mean(dat$y), theta = theta, edf = sum(fit$edf) - 1,
    max_diff_variance = max(abs(v_formula - fit$family$variance(mu))),
    diff_loglik       = abs(ll_formula - as.numeric(logLik(fit))),
    max_diff_resid    = max(abs(r_formula - res_pkg[, g])),
    max_diff_spline   = max(abs(s_fit - spl_pkg[, g]))))
}
verify_liver

write.csv(verify_liver, file.path(dir.out, "verification_spacedecorr_liver.csv"), row.names = FALSE)