#-----------------------------------------#
#------- Simulation Diff Net Niches ------#
#-----------------------------------------#

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
library(GLASSOO)
library(EstimateGroupNetwork)
source(file.path(dir.proj, "R", "FunctionsSmooth.R")) # modified functions from Giotto to get smoothed counts

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

sim_diffnet <- function(rept,num_genes, scenario, kprop = 0.1,perturb_prob){
  if(scenario == 1){
    dataniche <- fulldata %>% filter(cell_type == "macrophage",
                                     niche %in% c("tumor-stroma boundary","stroma"))
  }else if(scenario == 2){
    dataniche <- fulldata %>%
      filter(cell_type == "macrophage",
             niche %in% c("myeloid-enriched stroma","stroma"))
  }else if(scenario == 3){
    dataniche <- fulldata %>% filter(cell_type == "macrophage",
                                     niche %in% c("tumor-stroma boundary","myeloid-enriched stroma"))
  }
  sample_coords <- dataniche[,c("sdimx","sdimy","niche")]
  niches <- as.factor(unique(dataniche$niche))

  sample_coords$niche <- as.numeric(sample_coords$niche == levels(niches)[1])
  id0 <- which(sample_coords$niche == 0)
  id1 <- which(sample_coords$niche == 1)
  num_cells <- nrow(dataniche)
  n0<- length(id0)
  n1<-length(id1)


  set.seed(5675*rept);m <- sample(c(0.15,0.25,0.5,1), num_genes, replace = T, prob = c(1/3,2/3,2/3,1/3))
  set.seed(6786*rept);ratio <- sample(c(0.1,0.25,0.5,0.75), num_genes, replace = T,prob = c(1/4,2/4,2/4,1/4))
  v <- m/ratio
  libsize <- round(dataniche$Area*0.1)
  range <- sample(c(0.5,1, 1.5,2),replace = T, num_genes)

  # Precision matrix ----
  # Niche 0
  nei <- case_when(num_genes == 10 ~ 2,
                   num_genes == 20 ~ 4,
                   num_genes == 50 ~ 6 )
  prob <- case_when(num_genes == 10 ~ 0.1,
                    num_genes == 20 ~ 0.1,
                    num_genes == 50 ~ 0.05 )
  # Niche 0
  set.seed(68545*rept)
  g_1 <-sample_smallworld(dim = 1, size = num_genes, nei = nei, p = prob)
  adjmat1 <- as.matrix(as_adjacency_matrix(g_1))
  delta_1 <- adjmat1
  delta_1[delta_1 == 1 & upper.tri(delta_1)] <- runif(length(delta_1[delta_1 == 1 & upper.tri(delta_1)]),min = 0.5,max = 0.7)
  delta_1[lower.tri(delta_1)] <- t(delta_1)[lower.tri(delta_1)]
  delta_1 <-  delta_1 - diag(rep(min(as.numeric(eigen(delta_1)$values))-0.01,num_genes))
  delta_1 <- cov2cor(delta_1)

  # Niche 1
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

  delta_0 <- delta_1
  delta_1 <- delta_2

  # Simulate ----
  set.seed(67983*rept)
  tmp_niche0 <- MASS::mvrnorm(n=n0, mu = rep(0,num_genes),Sigma = delta_0)
  set.seed(6364*rept)
  tmp_niche1 <- MASS::mvrnorm(n=n1, mu = rep(0,num_genes),Sigma = delta_1)

  gamma_niche0 <- matrix(ncol = num_genes, nrow = n0)
  gamma_niche1 <- matrix(ncol = num_genes, nrow = n1)

  for(x in 1:num_genes){
    set.seed(x+787*rept)
    sigma_r <- telefit::maternCov(as.matrix(dist(sample_coords)),smoothness = 0.5, range = range[x] )
    tmpgamma <- MASS::mvrnorm(n=1, mu = rep(0,num_cells), Sigma=sigma_r)

    set.seed(x+51*rept)
    gamma_niche0[,x] <- qgamma(pnorm(tmp_niche0[,x] + tmpgamma[id0]),shape = (m[x]^2)/v[x], rate = m[x]/v[x])

    set.seed(x+632*rept)
    gamma_niche1[,x] <- qgamma(pnorm(tmp_niche1[,x] + tmpgamma[id1]),shape = (m[x]^2)/v[x], rate = m[x]/v[x])
  }
  # Niche 0
  Y_niche0 <- do.call(rbind, lapply(1:n0, function(sample){
    set.seed(7455+rept*sample)
    rpois(num_genes,libsize[sample]*gamma_niche0[sample,])
  }))
  Y0 <- cbind(Y_niche0,sample_coords[id0,])
  colnames(Y0)[1:num_genes] <- paste0("Y",1:num_genes)

  # Niche 1
  Y_niche1 <- do.call(rbind, lapply(1:n1, function(sample){
    set.seed(54232+rept*sample)
    rpois(num_genes,libsize[sample]*gamma_niche1[sample,])
  }))
  Y1 <- cbind(Y_niche1,sample_coords[id1,])
  colnames(Y1)[1:num_genes] <- paste0("Y",1:num_genes)

  # All
  data <- as.data.frame(rbind(Y0,Y1))
  data$cell_ID <- as.factor(1:num_cells)
  data$libsize <- libsize
  rownames(data) <- 1:num_cells

  res_true <- diffCorNet(delta_0,delta_1, n0, n1)

  id0 <- which(data$niche == 0)
  id1 <- which(data$niche == 1)

  datalog <- data
  datalog[,1:num_genes] <- apply(datalog[,1:num_genes], 2, function(x)log(10e4 * x/libsize +1))

  datavst <- t(vst(umi=t(data[,1:num_genes]),cell_attr = data,latent_var = "libsize")$y)
  colnames(datavst) <- paste0("Y",1:num_genes)


  # Get results ----
  ## Original or log transformed ----
  # Pearson
  res_pearson <- diffCorNet(cor(data[id0, 1:num_genes]),cor(data[id1, 1:num_genes]), n0, n1, pvalues = T)
  cor_pearson <- data.frame(Model = "Pearson",type = "Fisher",
                            true0 = upper(res_true$cormat1),
                            true1 = upper(res_true$cormat2),
                            trueDiff = upper(res_true$difference),
                            Niche0 = upper(res_pearson$cormat1),
                            Niche1 = upper(res_pearson$cormat2),
                            Diff = upper(res_pearson$difference)
  )
  res_pearson2 <- diffCorNet_twocor(data[id0, 1:num_genes],data[id1, 1:num_genes], nCores-1, pvalues = T)
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

  #Spearman
  res_spearman <- diffCorNet(cor(data[id0, 1:num_genes], method = "spearman"),cor(data[id1, 1:num_genes], method = "spearman"), n0, n1, pvalues = T)
  cor_spearman <- data.frame(Model = "Spearman",type = "Fisher",
                             true0 = upper(res_true$cormat1),
                             true1 = upper(res_true$cormat2),
                             trueDiff = upper(res_true$difference),
                             Niche0 = upper(res_spearman$cormat1),
                             Niche1 = upper(res_spearman$cormat2),
                             Diff = upper(res_spearman$difference)
  )

  # logpearson
  res_logpearson <- diffCorNet(cor(datalog[id0, 1:num_genes]),cor(datalog[id1, 1:num_genes]), n0, n1, pvalues = T)
  cor_logpearson <- data.frame(Model = "logPearson",type = "Fisher",
                               true0 = upper(res_true$cormat1),
                               true1 = upper(res_true$cormat2),
                               trueDiff = upper(res_true$difference),
                               Niche0 = upper(res_logpearson$cormat1),
                               Niche1 = upper(res_logpearson$cormat2),
                               Diff = upper(res_logpearson$difference)
  )
  res_logpearson2 <- diffCorNet_twocor(datalog[id0, 1:num_genes],datalog[id1, 1:num_genes], nCores-1, pvalues = T)
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
  res_spark <- diffCorNet(cor(datavst[id0, 1:num_genes]),cor(datavst[id1, 1:num_genes]), n0, n1, pvalues = T)
  cor_spark <- data.frame(Model = "SPARK",type = "Fisher",
                          true0 = upper(res_true$cormat1),
                          true1 = upper(res_true$cormat2),
                          trueDiff = upper(res_true$difference),
                          Niche0 = upper(res_spark$cormat1),
                          Niche1 = upper(res_spark$cormat2),
                          Diff = upper(res_spark$difference)
  )
  res_spark2 <- diffCorNet_twocor(datavst[id0, 1:num_genes],datavst[id1, 1:num_genes], nCores-1, pvalues = T)
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
  gobject0 <- createGiottoObject(raw_exprs = t(data[id0,1:num_genes]),
                                 norm_expr = t(datalog[id0,1:num_genes]),
                                 spatial_locs = data[id0,c("sdimx","sdimy")])
  gobject0 <- createSpatialNetwork(gobject = gobject0, minimum_k = 2)
  data_giotto0 <-  detectSpatialCorGenes(gobject0,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network')$cor_DT%>%
    graph_from_data_frame(., directed = F) %>%
    get.adjacency(., attr="spat_cor", sparse=FALSE)

  #niche1
  gobject1 <- createGiottoObject(raw_exprs = t(data[id1,1:num_genes]),
                                 norm_expr = t(datalog[id1,1:num_genes]),
                                 spatial_locs = data[id1,c("sdimx","sdimy")])
  gobject1 <- createSpatialNetwork(gobject = gobject1, minimum_k = 2)
  data_giotto1 <-  detectSpatialCorGenes(gobject1,
                                         method = 'network',
                                         spatial_network_name = 'Delaunay_network')$cor_DT%>%
    graph_from_data_frame(., directed = F) %>%
    get.adjacency(., attr="spat_cor", sparse=FALSE)

  res_giotto <- diffCorNet(data_giotto0,data_giotto1, n0, n1, pvalues = T)
  cor_giotto <- data.frame(Model = "Giotto",type = "Fisher",
                           true0 = upper(res_true$cormat1),
                           true1 = upper(res_true$cormat2),
                           trueDiff = upper(res_true$difference),
                           Niche0 = upper(res_giotto$cormat1),
                           Niche1 = upper(res_giotto$cormat2),
                           Diff = upper(res_giotto$difference)
  )


  ## meringue ----
  rownames(datalog) <- 1:num_cells
  # Niche 0
  w0 <- getSpatialNeighbors(data[id0,c("sdimx","sdimy")], filterDist = 1)
  colnames(w0) <- rownames(w0) <- id0
  data_mer0 <- spatialCrossCorMatrix(mat = t(datalog[id0,1:num_genes]), weight = w0)

  # Niche 1
  w1 <- getSpatialNeighbors(data[id1,c("sdimx","sdimy")], filterDist = 1)
  colnames(w1) <- rownames(w1) <- id1
  data_mer1 <- spatialCrossCorMatrix(mat = t(datalog[id1,1:num_genes]), weight = w1)

  res_meringue <- diffCorNet(data_mer0,data_mer1, n0, n1, pvalues = T)
  cor_meringue <- data.frame(Model = "Meringue",type = "Fisher",
                             true0 = upper(res_true$cormat1),
                             true1 = upper(res_true$cormat2),
                             trueDiff = upper(res_true$difference),
                             Niche0 = upper(res_meringue$cormat1),
                             Niche1 = upper(res_meringue$cormat2),
                             Diff = upper(res_meringue$difference)
  )

  ## spacedecorr ----
  data_spacedecorr <- spacedecorr(assay_matrix = data[,1:num_genes],
                                  metadata = data[,c("sdimx","sdimy","libsize")],
                                  libsize_col = "libsize",
                                  kprop = kprop,
                                  verbose = F,
                                  nCores = nCores-1)

  res_spacedecorr <- diffCorNet(cor(data_spacedecorr$residuals[id0,]),cor(data_spacedecorr$residuals[id1,]), n0, n1, pvalues = T)
  cor_spacedecorr <- data.frame(Model = "SpaceDecorr",type = "Fisher",
                                true0 = upper(res_true$cormat1),
                                true1 = upper(res_true$cormat2),
                                trueDiff = upper(res_true$difference),
                                Niche0 = upper(res_spacedecorr$cormat1),
                                Niche1 = upper(res_spacedecorr$cormat2),
                                Diff = upper(res_spacedecorr$difference)
  )
  res_spacedecorr2 <- diffCorNet_twocor(data_spacedecorr$residuals[id0, 1:num_genes],data_spacedecorr$residuals[id1, 1:num_genes], nCores-1, pvalues = T)
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
           num_genes=num_genes,
           scenario = scenario,
           kprop = kprop,
           perturb_prob=perturb_prob,
           nearpd=nearpd,
           n0=n0, n1=n1,
           labelniche0 = levels(niches)[1],
           labelniche1 = levels(niches)[2])
  return(df)
}

params <- expand.grid(rept = 1:50,
                      num_genes = c(20),
                      perturb_prob = c(0.25),
                      scenario = c(1:3),
                      kprop = c(0.1,0.3)
)
results <- with(params[args,],
                sim_diffnet(rept,num_genes, scenario, kprop,perturb_prob))

dir.create(file.path(dir.out, "results", "DNniche"), recursive = TRUE, showWarnings = FALSE)
write.csv(results, file = file.path(dir.out, "results", "DNniche",
                                    paste0("sim_", paste(sapply(params[args,], as.character), collapse = "_"), ".csv")))
