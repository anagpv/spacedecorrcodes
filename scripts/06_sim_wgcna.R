#---------------------------------------#
#------- Simulation co-expression ------#
#---------------- WGCNA ----------------#
#---------------------------------------#

dir.out  <- "/path/to/output/"            # folder where outputs will be saved

# Designed to run as a SLURM array job. To run locally, set args <- 1 and nCores <- 2
nCores <- as.numeric(Sys.getenv('SLURM_CPUS_ON_NODE'))
slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)


library(tidyverse)
library(mgcv)
library(data.table)
library(Rcpp)
library(MASS)
library(randnet)
library(CSCORE)
library(Seurat)
library(sctransform)
library(mclust)
library(igraph)
library(MERINGUE)
library(WRS2)
library(spacedecorr)
library(Giotto)

# Get estimates from real data
# load("Lung5-3data.Rdata")
# data <- fulldata[fulldata$cell_type == "macrophage",]
# est_realdata <- do.call(rbind, lapply(38:1017, function(i){
#   data$genei <- data[,..i]
#
#   # Just intercept
#   mod1 <- glm.nb(genei~1,data=data)
#   mu1 <- exp(coef(mod1)["(Intercept)"])
#   theta1 <- mod1$theta  # Dispersion parameter
#   sigma1 <- mu1 + (mu1^2 / theta1)
#
#
#   # Sequencing depth
#   mod2 <- glm.nb(genei~offset(log(totalcounts)),data=data)
#   mu2 <- exp(coef(mod2)["(Intercept)"])
#   theta2 <- mod2$theta  # Dispersion parameter
#   sigma2 <- mu2 + (mu2^2 / theta2)
#
#   print(sprintf("mu1 %.3f mu2 %.3f sigma1 %.3f sigma2 %.3f ",mu1, mu2,sigma1,sigma2))
#   data.frame(gene = colnames(counts)[i],
#              method = c("Just intercept","SeqDepth"),
#              mu = c(mu1,mu2),
#              sigma = c(sigma1,sigma2))
# }))
# write.csv(est_realdata, file = "NBestimatesRealData_macrophages.csv")
# est_realdata <- read.csv("NBestimatesRealData_macrophages.csv")
#
# counts <- data[,38:1017]
# colnames(counts)[apply(counts,2,var)>0.2]
# est_realdata <- est_realdata %>% filter(gene %in% colnames(counts)[apply(counts,2,var)>0.2] & method == "Just intercept")
# write.csv(est_realdata, file = "NBestimatesRealData_macrophages_filtered.csv")
est_realdata <- read.csv(file.path(dir.data, "NBestimatesRealData_macrophages_filtered.csv"))

# Generating data ---------------------------------------------------------

sim_wgcna <- function(rept,sim, num_cells = 1000,num_genes = 100, range = 1,
                      fixedlibsize = F,size = NA, kprop = 0.1,numblocks=4,
                      returnMat = F, load = F
){
  set.seed(5675*rept)
  if(kprop < 1) ksplines <- kprop * num_cells else ksplines <- kprop
  m <- sample(c(0.15,0.25,0.5), num_genes, replace = T, prob = c(1/3,2/3,1/3))
  ratio <- sample(c(0.1,0.25,0.5,0.75), num_genes, replace = T,prob = c(1/4,2/4,2/4,1/4))
  v <- m/ratio
  if(fixedlibsize){
    libsize <- rep(size,num_cells)
  }else{
    libsize<- sample(fulldata$totalcounts,num_cells)
  }

  sample_coords <- data.frame(sdimx = runif(num_cells,-1,1), sdimy = runif(num_cells,-1,1))

  # Getting adjacent matrix
  distmat <- as.matrix(dist(sample_coords, diag = T))
  distances <- reshape2::melt(distmat) %>% filter(value < 0.1) %>% dplyr::rename(from=Var1, to=Var2)
  sample_coords2 <- sample_coords %>% mutate(from = 1:n(), to = 1:n())
  nb <- merge(distances, sample_coords2[,1:3], by = "from")
  nb <- as.data.table(merge(nb,  sample_coords2[,c(1,2,4)], by = "to", suffixes=c("_from", "_to")))
  nb$to <- as.factor(nb$to)
  nb$from <- as.factor(nb$from)

  # Adjacent matrix
  adjW <- Matrix::sparseMatrix(i = distances$from,
                               j = distances$to,
                               x = 1, dims = c(num_cells, num_cells))
  nb_list <- lapply(seq_len(nrow(adjW)), function(i) {
    neighbors <- which(adjW[i, ] == 1)
    neighbors <- neighbors[neighbors != i]
    neighbors
  })
  names(nb_list) <- 1:num_cells

  #Simulating
  sigma_r <- telefit::maternCov(as.matrix(dist(sample_coords)),smoothness = 0.5, range = range )

  #--- Block matrix generated with:---
  # num_genes = 100
  # numblocks = 8
  # for(rept in 1:100){
  #   print(rept)
  #   set.seed(47289*rept)
  #   deltamat0 <- BlockModel.Gen(20,num_genes,K=numblocks
  #                               ,beta=1/1000,alpha=5, rho=0.9,simple=FALSE,power=TRUE)$A
  #   deltamat <- deltamat0
  #   deltamat[deltamat == 1 & upper.tri(deltamat)] <- runif(length(deltamat[deltamat == 1 & upper.tri(deltamat)]),min = 0.5,max = 0.9)
  #   deltamat[lower.tri(deltamat)] <- t(deltamat)[lower.tri(deltamat)]
  #   writeMat(paste0(dir.out,"SimulationWGCNA/OriginalWGCNAmatrix/mat_rept",rept,"_ngenes",num_genes,"_numblocks",numblocks,".mat"),
  #            A = deltamat)
  # }
  # NearPD matrix estimated with http://www.seas.ucla.edu/~vandenbe/publications/spcones.pdf
  deltamat <- as.matrix(read.delim(file.path(dir.data, "NearPDmatrix", paste0("matPD_rept",rept,"_ngenes",num_genes,"_numblocks",numblocks,".txt")),
                                   sep = " ", header = F))

  if(sim == "MV"){
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
      qgamma(pnorm(tmp[,x]),
             shape = (m[x]^2)/v[x], rate = m[x]/v[x])
    }) )

    Y <- do.call(rbind, lapply(1:num_cells, function(sample){
      Y <- rpois(num_genes,libsize[sample]*gamma[sample,])
    }))
  }else if(sim == "Additive"){
    tmp <- mvrnorm(num_cells,rep(0,num_genes),deltamat)
    tmp <- as.data.frame(scale(tmp)*0.5)
    gamma <- do.call(cbind, lapply(1:num_genes, function(x){
      qgamma(pnorm(tmp[,x] + as.numeric(scale(MASS::mvrnorm(n=1, mu = rep(0,num_cells), Sigma=sigma_r))*0.5)),
             shape = (m[x]^2)/v[x], rate = v[x]/m[x])
    }) )

    Y <- do.call(rbind, lapply(1:num_cells, function(sample){
      Y <- rpois(num_genes,libsize[sample]*gamma[sample,])
    }))
  }
  rownames(Y) <- 1:num_cells
  colnames(Y) <- paste0("Y",1:num_genes)

  data <- Y %>% cbind(sample_coords) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize
  rownames(data) <- 1:num_cells

  Yvst <- t(vst(umi=t(Y),cell_attr = data,latent_var = "libsize", min_cells = 1)$y)

  remove <- which(colSums(Y) == 0)
  if(length(remove) >0){
    Y <- Y[,-remove]
    deltamat <- deltamat[-remove,-remove]
  }
  data <- Y %>% cbind(sample_coords) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize
  rownames(data) <- 1:num_cells

  logY <- apply(Y, 2, function(x)log(10e4 * x/libsize +1))


  # Get results ----
  # Original or log transformed
  time.pearson <- system.time(cor_pearson <- cor(Y))[3]
  time.spearman <- system.time(cor_spearman <- cor(Y, method = "spearman"))[3]
  time.logpearson <- system.time(cor_logpearson <- cor(logY))[3]

  # SPARK
  time.spark <- system.time(cor_spark <- cor(Yvst))[3]

  # Giotto
  time.init <- Sys.time()
  gobject <- createGiottoObject(raw_exprs = t(Y),
                                norm_expr = t(logY),
                                spatial_locs = sample_coords)
  gobject <- createSpatialNetwork(gobject = gobject, minimum_k = 2)
  cor_giotto <-  detectSpatialCorGenes(gobject,
                                       method = 'network',
                                       spatial_network_name = 'Delaunay_network')$cor_DT %>%
    graph_from_data_frame(., directed = F) %>%
    get.adjacency(., attr="spat_cor", sparse=FALSE)
  time.end <- Sys.time()
  time.giotto <- as.numeric(time.end-time.init, units = "secs")


  # meringue
  time_init <- Sys.time()
  w <- getSpatialNeighbors(sample_coords, filterDist = .1)
  cor_mer <- spatialCrossCorMatrix(mat = t(logY), weight = w)
  time_end <- Sys.time()
  time.mer <- as.numeric(time_end-time_init, units = "secs")

  # spacedecorr
  time.spacedecorr <- system.time(
    res_spacedecorr <- spacedecorr(assay_matrix = Y,
                                   metadata = data,
                                   kprop = kprop,
                                   libsize_col = "libsize",
                                   verbose = F,
                                   nCores = nCores-1))[3]
  cor_res <- cor(res_spacedecorr$residuals)
  cor_win <- winall(res_spacedecorr$residuals)

  # Results ----
  corMat <- list(
    pearson = cor_pearson,
    spearman = cor_spearman,
    logpearson = cor_logpearson,
    spark = cor_spark,
    giotto = cor_giotto,
    meringue = cor_mer,
    spacedecorr = cor_res,
    spacedecorr_win = cor_win$cor
  )

  # dist cov
  true_labels <- cutree(hclust(dist(deltamat)), k = numblocks)
  ARI <- sapply(corMat, function(x){
    print("-")
    x[is.na(x)] <- 0
    predicted_labels <- cutree(hclust(dist(x)), k = numblocks)
    adjustedRandIndex(true_labels, predicted_labels)
  })

  if(returnMat){
    list(params = data.frame(rept,sim, num_cells ,num_genes , range,
                             fixedlibsize ,size , kprop ,numblocks),
         corMat = corMat,
         deltamat = deltamat,
         ARI=ARI, true_labels=true_labels,
         times = c(time.pearson,time.spearman,time.logpearson,time.spark,time.mer,
                   time.giotto)
    )
  }else{
    data.frame(rept,sim, num_cells ,num_genes , range,
               fixedlibsize ,size , kprop ,numblocks,
               model = c("pearson","spearman","logpearson","SPARK",
                         "meringue",
                         "giotto","spacedecorr","spacedecorr_win"
               ),
               ARI=ARI,
               times = c(time.pearson,time.spearman,time.logpearson,time.spark,time.mer,
                        time.giotto, time.spacedecorr,time.spacedecorr))

  }
}

params <- expand.grid(rept = 1:100,
                      num_genes = c(100),
                      range=c(0.25,0.5),
                      sim = c("Additive","MV"),
                      num_cells = c(2000),
                      numblocks = c(4,8,10),
                      kprop = c(200, 1000)
)


results <- with(params[args,],
                sim_wgcna(rept=rept,sim=sim, num_cells=num_cells,num_genes=num_genes,
                          range=range, fixedlibsize=F, load = F,
                          size=NA,kprop=kprop,numblocks=numblocks,returnMat=F))

dir.create(file.path(dir.out, "results", "wgcna"), recursive = TRUE, showWarnings = FALSE)
write.csv(results, file = file.path(dir.out, "results", "wgcna",
                                    paste0("sim_", paste(sapply(params[args,], as.character), collapse = "_"), ".csv")))
