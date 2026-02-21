#-------------------------------------#
#------- Simulation Correlation ------#
#-------------------------------------#

dir.out <- "~/SpaceDecorr_Analysis/"

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

# Auxiliar function -------------------------------------------------------
getAllCors <- function(y1,y2){
  require(WRS2)
  time.pearson<- system.time(pearson <- cor.test(y1,y2, method = "pearson"))[3]
  time.spearman<- system.time(spearman <- cor.test(y1,y2, method = "spearman", exact = F))[3]
  time.kendall<- system.time(kendall <- cor.test(y1,y2, method = "kendall"))[3]
  time.winsorized<- system.time(winsorized <- wincor(y1,y2, tr = 0.1))[3]
  time.percentbend<- system.time(percentbend <- pbcor(y1,y2, beta = 0.1))[3]

  time.winsorized2<- system.time(winsorized2 <- wincor(y1,y2, tr = 0.2))[3]
  time.percentbend2<- system.time(percentbend2 <- pbcor(y1,y2, beta = 0.2))[3]

  data.frame(method = c("pearson","spearman","kendall",
                        "winsorized(tr 0.1)","winsorized(tr 0.2)","percentbend(tr 0.1)","percentbend(tr 0.2)"),
             estimate = c(pearson$estimate, spearman$estimate, kendall$estimate,
                          winsorized$cor,winsorized2$cor, percentbend$cor,percentbend2$cor),
             pvalue = c(pearson$p.value, spearman$p.value, kendall$p.value,
                        winsorized$p.value,winsorized2$p.value, percentbend$p.value,percentbend2$p.value),
             timecor = c(time.pearson, time.spearman, time.kendall, time.winsorized,time.winsorized2,
                         time.percentbend,time.percentbend2))
}
# Generating data ---------------------------------------------------------
simCorrelation <- function(rept,sim, num_cells = 1000,range = 1,rho=0, m, ratio,libsize, kprop){
  v <- m/ratio

  set.seed(41); sample_coords <- data.frame(sdimx = runif(num_cells,-1,1), sdimy = runif(num_cells,-1,1))
  num_genes = 2

  # Simulating ----
  set.seed(47289*rept)
  sigma_r <- telefit::maternCov(as.matrix(dist(sample_coords)),smoothness = 0.5, range = range )
  deltamat <- matrix(c(1,rho,
                       rho,1),ncol=2)

  if(sim == "MV"){
    set.seed(67983*rept)
    tmp <- MASS::mvrnorm(n=1, mu = rep(0,num_cells*num_genes),Sigma=sigma_r %x% deltamat)
    tmp <- as.matrix(matrix(tmp, ncol = num_genes,byrow=T))

    gamma <- do.call(cbind, lapply(1:num_genes, function(x){
      set.seed(x*51*rept)
      qgamma(pnorm(tmp[,x]),
             shape = (m^2)/v, rate = m/v)
    }) )

    Y <- do.call(rbind, lapply(1:num_cells, function(sample){
      set.seed(7455+rept*sample)
      Y <- rpois(num_genes,libsize*gamma[sample,])
    }))

  }else if(sim == "Additive"){
    set.seed(43121*rept)
    tmp <- MASS::mvrnorm(num_cells,rep(0,num_genes),deltamat)
    gamma <- do.call(cbind, lapply(1:num_genes, function(x){
      set.seed(x*51*rept)
      qgamma(pnorm(tmp[,x] + MASS::mvrnorm(n=1, mu = rep(0,num_cells), Sigma=sigma_r)),
             shape = (m^2)/v, rate = v/m)
    }) )

    Y <- do.call(rbind, lapply(1:num_cells, function(sample){
      set.seed(7455+rept*sample)
      Y <- rpois(num_genes,libsize*gamma[sample,])
    }))
  }
  colnames(Y) <- paste0("Y",1:num_genes)
  rownames(Y) <- 1:num_cells

  data <- Y %>% cbind(sample_coords) %>% mutate(id = as.factor(1:n()))
  data$libsize <- libsize

  # Normalizing data
  time.log <- system.time(Ylog <- apply(Y, 2, function(x)log(10e4 * x/libsize +1)))[3]
  time.vst <- system.time(Yvst <- t(vst(umi=t(Y),cell_attr = data,latent_var = "libsize")$y))[3]
  time.spacenorm <- system.time(Ysd <- spacedecorr(assay_matrix = data[,c("Y1","Y2")],
                                                                metadata = data[,c("sdimx","sdimy","libsize")],
                                                                libsize_col = "libsize",
                                                                family = "nb",
                                                                kprop = kprop,
                                                                verbose = F))[3]
  cor_raw <- getAllCors(Y[,1],Y[,2])
  cor_raw$normalization <- "raw"
  cor_raw$timenorm <- 0

  cor_log <- getAllCors(Ylog[,1],Ylog[,2])
  cor_log$normalization <- "log"
  cor_log$timenorm <- time.log

  cor_vst <- getAllCors(Yvst[,1],Yvst[,2])
  cor_vst$normalization <- "vst"
  cor_vst$timenorm <- time.vst

  cor_sd <- getAllCors(Ysd$residuals[,1],Ysd$residuals[,2])
  cor_sd$normalization <- "spacedecorr"
  cor_sd$timenorm <- time.spacenorm

  rbind(cor_raw,cor_log,cor_vst,cor_sd) %>%
    mutate(rept=rept, sim=sim, num_cells=num_cells, range=range, rho=rho, m=m, ratio=ratio, libsize=libsize, kprop=kprop)

}

slurm_arrayid <- Sys.getenv('SLURM_ARRAY_TASK_ID')
args <- as.numeric(slurm_arrayid)
params <- expand.grid(range = c(0.25),
                      libsize = c(200),
                      sim = c("MV","Additive"),
                      m = c(0.1,0.5),
                      ratio = c(0.5),
                      num_cells = c(1000,2000,5000),
                      kprop = c(0.1,0.3),
                      rept = 1:500
)

if(params$num_cells[args] == 2000){
  rhovalues <- c(c(0, 0.01,0.025, 0.05, 0.10,0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.7, 0.8, 0.9),
                 -c(0.01,0.025, 0.05, 0.10,0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.7, 0.8, 0.9))
}else{
  rhovalues <- 0
}
results <- do.call(rbind, lapply(rhovalues, function(rholoop){
  with(params[args,],
       simCorrelation(rept,sim, num_cells,range,rholoop, m, ratio,libsize, kprop))
}))

out_path <- file.path(dir.out, "output", "simulations", "pairwise_robust")
dir.create(out_path, recursive = TRUE, showWarnings = FALSE)

write.csv(results, file = file.path(out_path,
                paste(sapply(params[args, ], as.character), collapse = "_") ,".csv",sep="_"))



