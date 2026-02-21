#------------------------------------------#
#------- Simulation Diff Net Samples ------#
#------------------------------------------#

dir.out  <- "/path/to/output/"            # folder where outputs will be saved
dir.proj <- "/path/to/repo/"              # root of this GitHub repo

# Designed to run as a SLURM array job. To run locally, set args <- 1 and nCores <- 2
nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)

library(Rcpp)
library(data.table)
library(tidyverse)
library(MASS)
library(mgcv)
library(spacedecorr)
library(MERINGUE)
library(CSCORE)
library(Seurat)
library(igraph)
library(GraphR)
library(Giotto)
library(sctransform)
library(WRS2)
library(EstimateGroupNetwork)

est_realdata <- read.csv(file.path(dir.data, "NBestimatesRealData_macrophages_filtered.csv"))
metrics <- c("TPR", "FPR", "FDR","Precision", "F1", "MCC")

# Auxiliar functions -----------------------------
diffCorNet <- function(cor_mat1, cor_mat2, n1, n2, adjust.method = "fdr", alpha = 0.05, pvalues = F){
  # Fisher Z-transform for group comparison
  fisher_z1 <- 0.5 * log((1 + cor_mat1) / (1 - cor_mat1))
  fisher_z2 <- 0.5 * log((1 + cor_mat2) / (1 - cor_mat2))

  se_diff <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))

  z_stat_diff <- (fisher_z1 - fisher_z2) / se_diff
  pval_diff <- 2 * pnorm(-abs(z_stat_diff))
  pval_diff_adj <- adjust_upper_triangle(pval_diff)

  # Single-group tests: Test if correlation is different from zero
  t_stat1 <- cor_mat1 * sqrt((n1 - 2) / (1 - cor_mat1^2))
  pval1 <- 2 * pt(-abs(t_stat1), df = n1 - 2)
  pval1_adj <- adjust_upper_triangle(pval1)

  t_stat2 <- cor_mat2 * sqrt((n2 - 2) / (1 - cor_mat2^2))
  pval2 <- 2 * pt(-abs(t_stat2), df = n2 - 2)
  pval2_adj <- adjust_upper_triangle(pval2)

  if(pvalues){
    df <- list(difference = pval_diff_adj,
               cormat1 = pval1_adj,
               cormat2 = pval2_adj)
  }else{
    df <-list(difference = apply(pval_diff_adj, 2, function(x) as.numeric(x < 0.05)),
              cormat1 = apply(pval1_adj, 2, function(x) as.numeric(x < 0.05)),
              cormat2 = apply(pval2_adj, 2, function(x) as.numeric(x < 0.05)))
  }
  return(df)
}

diffCorNet_twocor <- function(datamat1, datamat2,nCores, adjust.method = "fdr", alpha = 0.05, tr = 0.1, pvalues = F) {
  require(WRS2)
  require(parallel)
  p <- ncol(datamat1)

  pval_diff <- matrix(NA, ncol = p, nrow = p)
  colnames(pval_diff) <- colnames(datamat1)
  rownames(pval_diff) <- colnames(datamat1)

  if(nCores == 1){
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        # Difference in correlation between group 1 and group 2
        twocor_result <- tryCatch({
          twocor(datamat1[,i], datamat1[,j], datamat2[,i], datamat2[,j], beta = tr, nboot = 100,corfun = "pbcor")
        }, error = function(e) {
          NULL
        })

        if (!is.null(twocor_result)) {
          pval_diff[i,j] <- twocor_result$p.value
          pval_diff[j,i] <- twocor_result$p.value
        }
      }
    }
  }else{
    cl <- parallel::makeCluster(getOption("cl.cores", nCores))
    # Create all unique (i,j) pairs (upper triangle only)
    idx_pairs <- which(upper.tri(matrix(1, p, p)), arr.ind = TRUE)

    clusterExport(cl, varlist = c("datamat1", "datamat2", "tr","idx_pairs" ), envir = environment())
    clusterEvalQ(cl, library(WRS2))  # Load WRS2 in each worker

    # Define the worker function
    worker_fun <- function(pair_idx) {
      i <- idx_pairs[pair_idx, 1]
      j <- idx_pairs[pair_idx, 2]

      res <- tryCatch({
        twocor(datamat1[,i], datamat1[,j], datamat2[,i], datamat2[,j], tr = tr)$p.value
      }, error = function(e) {
        NA
      })

      return(list(i = i, j = j, pval = res))
    }

    # Run parallelized
    results_list <- parLapply(cl, 1:nrow(idx_pairs), worker_fun)

    stopCluster(cl)  # Stop the cluster

    for (res in results_list) {
      pval_diff[res$i, res$j] <- res$pval
    }

    pval_diff[lower.tri(pval_diff)] <- t(pval_diff)[lower.tri(pval_diff)]
    diag(pval_diff) <- 0
  }


  # Single-group tests: is correlation different from zero in group 1
  pval1 <- pball(datamat1, beta = tr, nboot = 100)$p.value
  pval2 <- pball(datamat2, beta = tr, nboot = 100)$p.value


  # Adjust p-values
  pval_diff_adj <- adjust_upper_triangle(pval_diff, adjust.method = adjust.method)
  pval1_adj <- adjust_upper_triangle(pval1, adjust.method = adjust.method)
  pval2_adj <- adjust_upper_triangle(pval2, adjust.method = adjust.method)

  if(pvalues){
    df <- list(difference = pval_diff_adj,
               cormat1 = pval1_adj,
               cormat2 = pval2_adj)
  }else{
    df <-list(difference = apply(pval_diff_adj, 2, function(x) as.numeric(x < 0.05)),
              cormat1 = apply(pval1_adj, 2, function(x) as.numeric(x < 0.05)),
              cormat2 = apply(pval2_adj, 2, function(x) as.numeric(x < 0.05)))
  }
  return(df)
}

# Simulate data ------------------------------------------------------
sim_diffnet <- function(rept,sim,num_genes, num_cells1, num_cells2, range = 1, kprop = 0.1,perturb_prob){
    set.seed(378); sample_coords1 <- data.frame(sdimx = runif(num_cells1,-1,1), sdimy = runif(num_cells1,-1,1))
    set.seed(354); sample_coords2 <- data.frame(sdimx = runif(num_cells2,-1,1), sdimy = runif(num_cells2,-1,1))

    set.seed(5675*rept);m <- sample(c(0.15,0.25,0.5), num_genes, replace = T, prob = c(1/3,2/3,1/3))
    set.seed(6786*rept);ratio <- sample(c(0.1,0.25,0.5,0.75), num_genes, replace = T,prob = c(1/4,2/4,2/4,1/4))
    v <- m/ratio
    set.seed(6786*rept);libsize1 <- sample(fulldata$totalcounts,num_cells1)
    set.seed(512*rept);libsize2 <- sample(fulldata$totalcounts,num_cells2)

    # Precision matrix ----
    nei <- case_when(num_genes == 10 ~ 2,
                     num_genes == 20 ~ 4,
                     num_genes == 50 ~ 6)
    prob <- case_when(num_genes == 10 ~ 0.1,
                      num_genes == 20 ~ 0.1,
                      num_genes == 50 ~ 0.05)
    # Tissue 1
    set.seed(68545*rept)
    g_1 <-sample_smallworld(dim = 1, size = num_genes, nei = nei, p = prob)
    adjmat1 <- as.matrix(as_adjacency_matrix(g_1))
    delta_1 <- adjmat1
    delta_1[delta_1 == 1 & upper.tri(delta_1)] <- runif(length(delta_1[delta_1 == 1 & upper.tri(delta_1)]),min = 0.5,max = 0.7)
    delta_1[lower.tri(delta_1)] <- t(delta_1)[lower.tri(delta_1)]
    delta_1 <-  delta_1 - diag(rep(min(as.numeric(eigen(delta_1)$values))-0.01,num_genes))
    delta_1 <- cov2cor(delta_1)

    # Tissue 2
    set.seed(5331*rept)
    delta_2 <- delta_1
    edges2 <- which(upper.tri(delta_2), arr.ind = TRUE)
    perturb_idx2 <- sample(1:nrow(edges2), size = floor(perturb_prob * nrow(edges2)))
    for (idx in perturb_idx2) {
      i <- edges2[idx, 1]
      j <- edges2[idx, 2]

      if(delta_2[i, j] != 0){
        delta_2[i, j] <- delta_2[j, i] <- 0
      }else{
        delta_2[i, j] <- delta_2[j, i] <- runif(1, 0.5, 0.9)
      }
    }
    delta_2 <-  delta_2 - diag(rep(min(as.numeric(eigen(delta_2)$values))-0.01,num_genes))

    # Simulate ----
    sigma_r1 <- telefit::maternCov(as.matrix(dist(sample_coords1)),smoothness = 0.5, range = range )
    sigma_r2 <- telefit::maternCov(as.matrix(dist(sample_coords2)),smoothness = 0.5, range = range )
    if(sim == "MV"){
      # Tissue 1---
      # Eigen-decomposition for low-rank approx of U and V
      eig_sigma1 <- eigen(sigma_r1, symmetric = TRUE)
      eig_delta1 <- eigen(delta_1, symmetric = TRUE)

      # Keep top k components
      r1 <- length(which(eig_sigma1$values > 0))
      s1 <- length(which(eig_delta1$values > 0))
      A1 <- eig_sigma1$vectors[, 1:r1] %*% diag(sqrt(eig_sigma1$values[1:r1]))
      B1 <- eig_delta1$vectors[, 1:s1] %*% diag(sqrt(eig_delta1$values[1:s1]))

      # Simulate Z ~ N(0, I)
      set.seed(689*rept)
      Z1 <- matrix(rnorm(r1 * s1), nrow = r1, ncol = s1)

      # Construct X
      tmp1 <- A1 %*% Z1 %*% t(B1)

      gamma1 <- do.call(cbind, lapply(1:num_genes, function(x){
        set.seed(x+51*rept)
        qgamma(pnorm(tmp1[,x]),
               shape = (m[x]^2)/v[x], rate = m[x]/v[x])
      }) )

      Y1 <- do.call(rbind, lapply(1:num_cells1, function(sample){
        set.seed(7455+rept*sample)
        Y1 <- rpois(num_genes,libsize1[sample]*gamma1[sample,])
      }))

      # Tissue 2---
      # Eigen-decomposition for low-rank approx of U and V
      eig_sigma2 <- eigen(sigma_r2, symmetric = TRUE)
      eig_delta2 <- eigen(delta_2, symmetric = TRUE)

      # Keep top k components
      r2 <- length(which(eig_sigma2$values > 0))
      s2 <- length(which(eig_delta2$values > 0))
      A2 <- eig_sigma2$vectors[, 1:r2] %*% diag(sqrt(eig_sigma2$values[1:r2]))
      B2 <- eig_delta2$vectors[, 1:s2] %*% diag(sqrt(eig_delta2$values[1:s2]))

      # Simulate Z ~ N(0, I)
      set.seed(435*rept)
      Z2 <- matrix(rnorm(r2 * s2), nrow = r2, ncol = s2)

      # Construct X
      tmp2 <- A2 %*% Z2 %*% t(B2)

      gamma2 <- do.call(cbind, lapply(1:num_genes, function(x){
        set.seed(x+64*rept)
        qgamma(pnorm(tmp2[,x]),
               shape = (m[x]^2)/v[x], rate = m[x]/v[x])
      }) )

      Y2 <- do.call(rbind, lapply(1:num_cells2, function(sample){
        set.seed(642+rept*sample)
        Y2 <- rpois(num_genes,libsize2[sample]*gamma2[sample,])
      }))
    }else if(sim == "Additive"){
      # Tissue 1---
      set.seed(121*rept)
      tmp1 <- mvrnorm(num_cells1,rep(0,num_genes),delta_1)
      gamma1 <- do.call(cbind, lapply(1:num_genes, function(x){
        set.seed(x+51*rept)
        qgamma(pnorm(tmp1[,x] + MASS::mvrnorm(n=1, mu = rep(0,num_cells1), Sigma=sigma_r1)),
               shape = (m[x]^2)/v[x], rate = v[x]/m[x])
      }) )

      Y1 <- do.call(rbind, lapply(1:num_cells1, function(sample){
        set.seed(7455+rept*sample)
        Y1 <- rpois(num_genes,libsize1[sample]*gamma1[sample,])
      }))

      # Tissue 2---
      set.seed(43*rept)
      tmp2 <- mvrnorm(num_cells2,rep(0,num_genes),delta_2)
      gamma2 <- do.call(cbind, lapply(1:num_genes, function(x){
        set.seed(x+63*rept)
        qgamma(pnorm(tmp2[,x] + MASS::mvrnorm(n=1, mu = rep(0,num_cells1), Sigma=sigma_r2)),
               shape = (m[x]^2)/v[x], rate = v[x]/m[x])
      }) )

      Y2 <- do.call(rbind, lapply(1:num_cells2, function(sample){
        set.seed(7647+rept*sample)
        Y2 <- rpois(num_genes,libsize2[sample]*gamma2[sample,])
      }))
    }
    colnames(Y1) <- paste0("Y",1:num_genes)
    colnames(Y2) <- paste0("Y",1:num_genes)

    # data
    data1 <- as.data.frame(cbind(Y1,sample_coords1))
    data1$cell_ID <- as.factor(1:num_cells1)
    data1$libsize <- libsize1
    rownames(data1) <- 1:num_cells1

    data2 <- as.data.frame(cbind(Y2,sample_coords2))
    data2$cell_ID <- as.factor(1:num_cells2)
    data2$libsize <- libsize2
    rownames(data2) <- 1:num_cells2

  res_true <- diffCorNet(delta_1,delta_2, num_cells1, num_cells2)

  datalog1 <- data1
  datalog1[,1:num_genes] <- apply(datalog1[,1:num_genes], 2, function(x)log(10e4 * x/data1$libsize +1))

  datalog2 <- data2
  datalog2[,1:num_genes] <- apply(datalog2[,1:num_genes], 2, function(x)log(10e4 * x/data2$libsize +1))

  datavst1 <- data1[,1:num_genes]
  datavst1[,which(colSums(datavst1) != 0)] <- t(vst(umi=t(data1[,which(colSums(datavst1) != 0)]),cell_attr = data1,latent_var = "libsize")$y)
  colnames(datavst1) <- paste0("Y",1:num_genes)

  datavst2 <- data2[,1:num_genes]
  datavst2[,which(colSums(datavst2) != 0)] <- t(vst(umi=t(data2[,which(colSums(datavst2) != 0)]),cell_attr = data2,latent_var = "libsize")$y)
  colnames(datavst2) <- paste0("Y",1:num_genes)


  # Get results ----
  ## Original or log transformed ----
  # Pearson
  res_pearson <- diffCorNet(cor(data1[, 1:num_genes]),cor(data2[, 1:num_genes]), num_cells1, num_cells2, pvalues = T)
  cor_pearson <- data.frame(Model = "Pearson",
                            type = "Fisher",
                            true0 = upper(res_true$cormat1),
                            true1 = upper(res_true$cormat2),
                            trueDiff = upper(res_true$difference),
                            Niche0 = upper(res_pearson$cormat1),
                            Niche1 = upper(res_pearson$cormat2),
                            Diff = upper(res_pearson$difference)
  )
  res_pearson2 <- diffCorNet_twocor(data1[, 1:num_genes],data2[, 1:num_genes], nCores-1, pvalues = T)
  cor_pearson <- rbind(cor_pearson,
                       data.frame(Model = "Pearson",
                                  type = "twocor",
                                  true0 = upper(res_true$cormat1),
                                  true1 = upper(res_true$cormat2),
                                  trueDiff = upper(res_true$difference),
                                  Niche0 = upper(res_pearson2$cormat1),
                                  Niche1 = upper(res_pearson2$cormat2),
                                  Diff = upper(res_pearson2$difference)
                       ))
  # #Spearman
  res_spearman <- diffCorNet(cor(data1[, 1:num_genes], method = "spearman"),cor(data2[, 1:num_genes], method = "spearman"), num_cells1, num_cells2, pvalues = T)
  cor_spearman <- data.frame(Model = "Spearman",
                             type = "Fisher",
                             true0 = upper(res_true$cormat1),
                             true1 = upper(res_true$cormat2),
                             trueDiff = upper(res_true$difference),
                             Niche0 = upper(res_spearman$cormat1),
                             Niche1 = upper(res_spearman$cormat2),
                             Diff = upper(res_spearman$difference)
  )


  # logpearson
  res_logpearson <- diffCorNet(cor(datalog1[, 1:num_genes]),cor(datalog2[, 1:num_genes]), num_cells1, num_cells2, pvalues = T)
  cor_logpearson <- data.frame(Model = "logPearson",
                               type = "Fisher",
                               true0 = upper(res_true$cormat1),
                               true1 = upper(res_true$cormat2),
                               trueDiff = upper(res_true$difference),
                               Niche0 = upper(res_logpearson$cormat1),
                               Niche1 = upper(res_logpearson$cormat2),
                               Diff = upper(res_logpearson$difference)
  )
  res_logpearson2 <- diffCorNet_twocor(datalog1[, 1:num_genes],datalog2[, 1:num_genes], nCores-1, pvalues = T)
  cor_logpearson <- rbind(cor_logpearson,
                          data.frame(Model = "logPearson",
                                     type = "twocor",
                                     true0 = upper(res_true$cormat1),
                                     true1 = upper(res_true$cormat2),
                                     trueDiff = upper(res_true$difference),
                                     Niche0 = upper(res_logpearson2$cormat1),
                                     Niche1 = upper(res_logpearson2$cormat2),
                                     Diff = upper(res_logpearson2$difference)
                          ))
  #SPARK ----
  res_spark <- diffCorNet(cor(datavst1[, 1:num_genes]),cor(datavst2[, 1:num_genes]), num_cells1, num_cells2, pvalues = T)
  cor_spark <- data.frame(Model = "SPARK",type = "Fisher",
                          true0 = upper(res_true$cormat1),
                          true1 = upper(res_true$cormat2),
                          trueDiff = upper(res_true$difference),
                          Niche0 = upper(res_spark$cormat1),
                          Niche1 = upper(res_spark$cormat2),
                          Diff = upper(res_spark$difference)
  )
  res_spark2 <- diffCorNet_twocor(datavst1[, 1:num_genes],datavst2[, 1:num_genes], nCores-1, pvalues = T)
  cor_spark <- rbind(cor_spark,
                     data.frame(Model = "SPARK",
                                type = "twocor",
                                true0 = upper(res_true$cormat1),
                                true1 = upper(res_true$cormat2),
                                trueDiff = upper(res_true$difference),
                                Niche0 = upper(res_spark2$cormat1),
                                Niche1 = upper(res_spark2$cormat2),
                                Diff = upper(res_spark2$difference)
                     ))
  ## Giotto and KNN ----
  #niche0
  gobject0 <- createGiottoObject(raw_exprs = t(data1[,1:num_genes]),
                                 norm_expr = t(datalog1[,1:num_genes]),
                                 spatial_locs = data1[,c("sdimx","sdimy")])
  gobject0 <- createSpatialNetwork(gobject = gobject0, minimum_k = 2)
  data_giotto0 <-  detectSpatialCorGenes(gobject0,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network')$cor_DT%>%
    graph_from_data_frame(., directed = F) %>%
    get.adjacency(., attr="spat_cor", sparse=FALSE)

  #niche1
  gobject1 <- createGiottoObject(raw_exprs = t(data2[,1:num_genes]),
                                 norm_expr = t(datalog2[,1:num_genes]),
                                 spatial_locs = data2[,c("sdimx","sdimy")])
  gobject1 <- createSpatialNetwork(gobject = gobject1, minimum_k = 2)
  data_giotto1 <-  detectSpatialCorGenes(gobject1,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network')$cor_DT%>%
    graph_from_data_frame(., directed = F) %>%
    get.adjacency(., attr="spat_cor", sparse=FALSE)

  res_giotto <- diffCorNet(data_giotto0,data_giotto1, num_cells1, num_cells2, pvalues = T)
  cor_giotto <- data.frame(Model = "Giotto",type = "Fisher",
                           true0 = upper(res_true$cormat1),
                           true1 = upper(res_true$cormat2),
                           trueDiff = upper(res_true$difference),
                           Niche0 = upper(res_giotto$cormat1),
                           Niche1 = upper(res_giotto$cormat2),
                           Diff = upper(res_giotto$difference)
  )


  ## meringue ----
  rownames(datalog1) <- 1:num_cells1
  rownames(datalog2) <- 1:num_cells2
  # Niche 0
  w0 <- getSpatialNeighbors(data1[,c("sdimx","sdimy")], filterDist = 1)
  colnames(w0) <- rownames(w0) <- 1:num_cells1
  data_mer0 <- spatialCrossCorMatrix(mat = t(datalog1[,1:num_genes]), weight = w0)

  # Niche 1
  w1 <- getSpatialNeighbors(data2[,c("sdimx","sdimy")], filterDist = 1)
  colnames(w1) <- rownames(w1) <- 1:num_cells2
  data_mer1 <- spatialCrossCorMatrix(mat = t(datalog2[,1:num_genes]), weight = w1)

  res_meringue <- diffCorNet(data_mer0,data_mer1, num_cells1, num_cells2, pvalues = T)
  cor_meringue <- data.frame(Model = "Meringue",type = "Fisher",
                             true0 = upper(res_true$cormat1),
                             true1 = upper(res_true$cormat2),
                             trueDiff = upper(res_true$difference),
                             Niche0 = upper(res_meringue$cormat1),
                             Niche1 = upper(res_meringue$cormat2),
                             Diff = upper(res_meringue$difference)
  )


  ## spacedecorr ----
  data_spacedecorr1 <- spacedecorr(assay_matrix = data1[,1:num_genes],
                                   metadata = data1[,c("sdimx","sdimy","libsize")],
                                   libsize_col = "libsize",
                                   kprop = kprop,
                                   verbose = F,
                                   nCores = nCores-1)
  data_spacedecorr2 <- spacedecorr(assay_matrix = data2[,1:num_genes],
                                   metadata = data2[,c("sdimx","sdimy","libsize")],
                                   libsize_col = "libsize",
                                   kprop = kprop,
                                   verbose = F,
                                   nCores = nCores-1)


  res_spacedecorr <- diffCorNet(cor(data_spacedecorr1$residuals),cor(data_spacedecorr2$residuals), num_cells1, num_cells2, pvalues = T)
  cor_spacedecorr <- data.frame(Model = "SpaceDecorr",type = "Fisher",
                                true0 = upper(res_true$cormat1),
                                true1 = upper(res_true$cormat2),
                                trueDiff = upper(res_true$difference),
                                Niche0 = upper(res_spacedecorr$cormat1),
                                Niche1 = upper(res_spacedecorr$cormat2),
                                Diff = upper(res_spacedecorr$difference)
  )
  res_spacedecorr2 <- diffCorNet_twocor(data_spacedecorr1$residuals,data_spacedecorr2$residuals, nCores-1, pvalues = T)
  cor_spacedecorr <- rbind(cor_spacedecorr,
                           data.frame(Model = "SpaceDecorr",
                                      type = "twocor",
                                      true0 = upper(res_true$cormat1),
                                      true1 = upper(res_true$cormat2),
                                      trueDiff = upper(res_true$difference),
                                      Niche0 = upper(res_spacedecorr2$cormat1),
                                      Niche1 = upper(res_spacedecorr2$cormat2),
                                      Diff = upper(res_spacedecorr2$difference)
                           ))

  df <- rbind(cor_pearson, cor_spearman, cor_logpearson, cor_spark,
              cor_giotto, cor_meringue,
              cor_spacedecorr
  ) %>%
    mutate(rept = rept,
           sim=sim,
           num_genes=num_genes,
           num_cells1=num_cells1,
           num_cells2=num_cells2,
           kprop = kprop,
           perturb_prob=perturb_prob,
           range=range,
           rev=rev,
           fixedrate=fixedrate)
  return(df)
}

params <- expand.grid(rept = 1:50,
                      sim = c("Additive","MV"),
                      num_genes = c(10,20),
                      num_cells = c(500, 2000, 5000),
                      perturb_prob = c(0.25),
                      kprop = c(0.1,0.3),
                      range = c(0.25)
)

results <- with(params[args,],
                sim_diffnet(rept,sim,num_genes,num_cells,num_cells, range, kprop,perturb_prob))

dir.create(file.path(dir.out, "results", "DNsample"), recursive = TRUE, showWarnings = FALSE)
write.csv(results, file = file.path(dir.out, "results", "DNsample",
                                    paste0("sim_", paste(sapply(params[args,], as.character), collapse = "_"), ".csv")))
