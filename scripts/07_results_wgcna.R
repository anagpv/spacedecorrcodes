#-----------------------------------------------#
#------- Results Simulation co-expression ------#
#-----------------------------------------------#

library(tidyverse)

# Set directories  --------------------------------------------------------
dir.data <- "/path/to/zenodo/download/"
mycol <- c("#E41A1C", "#377EB8","#4DAF4A","#FF7F00","#984EA3")

# NOTE: The CSV files read below were pre-generated from parallel simulation outputs.
# The code used to generate them is retained here for reference but is commented out,
# as it requires the full set of partial results files which are not included in the repo.

# Results -----------------------------------------------------------------
# filepath <- paste0(dir.out,"SimulationWGCNA/PartialWGCNAfinal/")
# files <- list.files(filepath,pattern = '\\.csv')
# results_wgcna0 <- do.call(rbind, lapply(files, function(file)data.table::fread(paste0(filepath,file))))[,-1]
# results_wgcna <- results_wgcna0 %>%
#   mutate(model = factor(model, levels = c("logpearson","SPARK","giotto","meringue","spacedecorr"),
#                         labels = c("log-CPM","SPARK","GIOTTO","MERINGUE","SpaceDecorr"))) %>%
#   filter(!is.na(model),
#          range == 0.5, num_cells == 2000, kprop == 0.1, numblocks == 4)
# write.csv(results_wgcna, file = paste0(dir.out,"WGCNA_results.csv"))
results_wgcna <- read.csv(file.path(dir.data,"WGCNA_results.csv"))

results_wgcna %>%
  ggplot(aes(x = model, y = ARI, col = model)) +
  facet_wrap(sim  ~ .)+
  geom_boxplot() +
  scale_color_manual(values = mycol)+
  theme_bw() +
  labs(x = "", y = "ARI") +
  theme(legend.position = "none",
        text = element_text(size = 13),
        strip.text = element_text(size=13),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = NA, colour = NA)) ;p_ari
