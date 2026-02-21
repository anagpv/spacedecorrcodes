#-------------------------------------#
#------- Simulation Correlation ------#
#-------------------------------------#

dir.out  <- "/path/to/output/"            # folder where outputs will be saved


library(tidyverse)
library(mgcv)
library(data.table)
library(Rcpp)
library(spacedecorr)
library(MERINGUE)
library(CSCORE)
library(Seurat)
library(GraphR)
library(sctransform)
library(igraph)
library(EstimateGroupNetwork)
suppressPackageStartupMessages(library(Giotto))

# Generating data ---------------------------------------------------------
cct <- function(pvec,weight=NULL){ #cct combination test
  J = length(pvec)
  if(is.null(weight)) weight = rep(1/J,J)
  t0 = sum(tan((0.5-pvec)*pi)*weight)
  pval = 0.5-atan(t0/sum(weight))/pi
  pval
}
simCorrelation <- function(rept,sim, num_cells,range,rho, m, ratio,libsize,kprop){
  set.seed(47289*rept)
  ksplines <- ifelse(kprop < 1, round(kprop * num_cells), kprop)
  v <- m/ratio
  sample_coords <- data.frame(sdimx = runif(num_cells,-1,1), sdimy = runif(num_cells,-1,1))
  num_genes = 3

  # Simulating ----
  sigma_r <- telefit::maternCov(as.matrix(dist(sample_coords)),smoothness = 0.5, range = range )
  deltamat <- matrix(c(1,rho,0,
                       rho,1,0,
                       0,0,1),ncol=3)

  if(sim == "MV"){
    # tmp <- MASS::mvrnorm(n=1, mu = rep(0,num_cells*num_genes),Sigma=sigma_r %x% deltamat)
    # tmp <- as.matrix(matrix(tmp, ncol = num_genes,byrow=T))

    # Eigen-decomposition for low-rank approx of U and V
    eig_sigma <- eigen(sigma_r, symmetric = TRUE)
    eig_delta <- eigen(deltamat, symmetric = TRUE)

    # Keep top k components
    r <- length(which(eig_sigma$values > 0))
    s <- length(which(eig_delta$values > 0))
    A <- eig_sigma$vectors[, 1:r] %*% diag(sqrt(eig_sigma$values[1:r]))
    B <- eig_delta$vectors[, 1:s] %*% diag(sqrt(eig_delta$values[1:s]))

    # Simulate Z ~ N(0, I)
    Z <- matrix(rnorm(r * s), nrow = r, ncol = s)

    # Construct X
    tmp <- A %*% Z %*% t(B)
    tmp <- as.data.frame(scale(tmp)*0.5)

    gamma <- do.call(cbind, lapply(1:num_genes, function(x){
      qgamma(pnorm(tmp[,x]), shape = (m^2)/v, rate = m/v)
    }) )

    Y <- do.call(rbind, lapply(1:num_cells, function(sample){
      Y <- rpois(num_genes,libsize*gamma[sample,])
    }))

  }else if(sim == "Additive"){
    tmp <- MASS::mvrnorm(num_cells,rep(0,num_genes),deltamat)
    tmp <- as.data.frame(scale(tmp)*0.5)

    gamma <- do.call(cbind, lapply(1:num_genes, function(x){
      qgamma(pnorm(tmp[,x] + as.numeric(scale(MASS::mvrnorm(n=1, mu = rep(0,num_cells), Sigma=sigma_r))*0.5)),
             shape = (m^2)/v, rate = v/m)
    }) )

    Y <- do.call(rbind, lapply(1:num_cells, function(sample){
      Y <- rpois(num_genes,libsize*gamma[sample,])
    }))
  }
  colnames(Y) <- paste0("Y",1:num_genes)
  rownames(Y) <- 1:num_cells

  data <- Y %>% cbind(sample_coords) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize

  # Normalizing data
  logY <- apply(Y, 2, function(x)log(10e4 * x/libsize +1))

  time.vst <- system.time(Yvst <- t(vst(umi=t(Y),cell_attr = data,latent_var = "libsize")$y))[3]

  # Get results ----
  # Original or log transformed
  time.pearson <- system.time(cor_pearson <- cor.test(Y[,1],Y[,2]))[3]
  time.spearman <- system.time(cor_spearman <- cor.test(Y[,1],Y[,2], method = "spearman"))[3]
  time.logpearson <- system.time(cor_logpearson <- cor.test(logY[,1],logY[,2]))[3]

  # SPARK
  time.spark <- system.time(cor_spark <- cor.test(Yvst[,1],Yvst[,2]))[3]
  time.spatk <- time.spark + time.vst

  # Giotto
  time.init <- Sys.time()
  suppressWarnings({
    gobject <- createGiottoObject(raw_exprs = t(Y),
                                  norm_expr = t(logY),
                                  spatial_locs = sample_coords)
  })
  gobject <- createSpatialNetwork(gobject = gobject, minimum_k = 2)
  cor_giotto <-  detectSpatialCorGenes(gobject,
                                       method = 'network',
                                       spatial_network_name = 'Delaunay_network')$cor_DT
  cor_giotto <- cor_giotto %>% filter(gene_ID =="Y1", variable == "Y2") %>% pull(spat_cor)
  pgiotto <- 2*pt(-abs(cor_giotto/(sqrt(1-cor_giotto^2)) * sqrt(num_cells - 2)), num_cells-2)
  time.end <- Sys.time()
  time.giotto <- as.numeric(time.end-time.init, units = "secs")

  # GraphR
  time_init <- Sys.time()
  res_graphRvst <- tryCatch({GraphR_est(features = as.matrix(Yvst),
                                        cont_external = as.matrix(data[,c("sdimx","sdimy")]))},
                            error = function(e)e)
  if(inherits(res_graphRvst, "error")){
    pred_graphRvst <- data.frame(partialcor=NA )
  }else{
    pred_graphRvst <- GraphR_pred(data[,c("sdimx","sdimy")], res_graphRvst)
    pred_graphRvst <- pred_graphRvst %>% left_join(data[,c("sdimx","sdimy")], by = c("sdimx","sdimy")) %>%
      mutate(CorrPIP = Correlation * Pr_inclusion) %>%
      group_by(feature1, feature2) %>%
      summarise(partialcor = mean(CorrPIP),
                cct = cct(FDR_p)) %>%
      filter(feature1 %in% c("Y1","Y2"),
             feature2 %in% c("Y1","Y2"))
  }
  time_end <- Sys.time()
  time.graphR <- as.numeric(time_end - time_init,units = "secs")
  time.graphR <- time.graphR + time.vst

  # meringue
  time_init <- Sys.time()
  w <- getSpatialNeighbors(sample_coords, filterDist = 0.1)
  cor_mer <- spatialCrossCorMatrix(mat = t(logY), weight = w)[1,2]
  pmer <- 2*pt(-abs(cor_mer/(sqrt(1-cor_mer^2)) * sqrt(num_cells - 2)), num_cells-2)
  time_end <- Sys.time()
  time.mer <- as.numeric(time_end-time_init, units = "secs")


  # spacedecorr
  print("spacedecorr")
  time.spacedecorr <- system.time(res_spacedecorr <- spacedecorr(assay_matrix = data[,c("Y1","Y2")],
                                                                        metadata = data,
                                                                        basis = "ts",
                                                                        libsize_col = "libsize",
                                                                        family = "nb",
                                                                        k = ksplines,
                                                                        verbose = F))[3]
  cor_res <- cor.test(as.numeric(res_spacedecorr$residuals[,1]),as.numeric(res_spacedecorr$residuals[,2]))

  # Results
  data.frame(seed = rept, sim, rho, range,num_cells,
             m,ratio,v,libsize, ksplines,kprop,
             model = c("pearson","spearman","logpearson","SPARK",
                       "meringue","giotto","spacedecorr","graphR(cct)"
             ),
             cor = c(cor_pearson$estimate,cor_spearman$estimate,
                     cor_logpearson$estimate,cor_spark$estimate,
                     cor_mer, cor_giotto,
                     cor_res$estimate,pred_graphRvst$partialcor
             ),
             pvalue = c(cor_pearson$p.value,cor_spearman$p.value,
                        cor_logpearson$p.value,cor_spark$p.value,
                        pmer, pgiotto,
                        cor_res$p.value,pred_graphRvst$cct
             ),
             times = c(time.pearson,time.spearman,time.logpearson,time.spark,time.mer,
                       time.giotto,  time.spacedecorr,time.graphR
             )
  )
}

# This script is designed to run as a SLURM array job.
# Each array task ID corresponds to one row of the parameter grid below.
# To run locally, set args manually, e.g.: args <- 1
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)
params <- expand.grid(rept = 1:200,
                      range = c(0.25, 0.5),
                      m = c(0.5),
                      rho = c(0, 0.01,0.025, 0.05, 0.10,0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.90),
                      ratio = c(0.5),
                      libsize = c(200),
                      sim = c("MV","Additive"),
                      num_cells = c(500,1000,2000,5000,10000,25000),
                      kprop = c(100, 200, 1000)
) %>%
  mutate(rho = ifelse(!(num_cells %in% c(2000)), 0, rho)) %>%
  distinct

results <- with(params[args,],
                simCorrelation(rept=rept,sim=sim, num_cells = num_cells,
                               range=range,rho=rho, m=m, ratio=ratio,libsize=libsize,kprop=kprop))

out_path <- file.path(dir.out, "output", "simulations", "pairwise")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

write.csv(results, file = file.path(out_path,
                                    paste0("sim_", paste(sapply(params[args,], as.character), collapse = "_"), ".csv")))
