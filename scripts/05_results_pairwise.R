#-----------------------------------------------#
#------- Results Simulation co-expression ------#
#-----------------------------------------------#

library(tidyverse)
library(cowplot)
library(ggrepel)

# Set directories  --------------------------------------------------------
dir.data <- "/path/to/zenodo/download/"
mycol <- c("#E41A1C", "#377EB8","#4DAF4A","#FF7F00","#984EA3"  )

# NOTE: The CSV files read below were pre-generated from parallel simulation outputs.
# The code used to generate them is retained here for reference but is commented out,
# as it requires the full set of partial results files which are not included in the repo.

# Results -----------------------------------------------------------------
# Combining partial results ----
# filepath <- paste0(dir.out,"Partial/Pairwise/")
# files <- list.files(filepath,pattern = '\\.csv')
# resultsCorr0 <- do.call(rbind, lapply(files, function(file)data.table::fread(paste0(filepath,file))))[,-1]
# resultsCorr <- resultsCorr0%>%
#   mutate(model = factor(model, levels = c("pearson","logpearson","spearman",
#                                           "SPARK","giotto","meringue","graphR(cct)","spacedecorr"),
#                         labels = c("Pearson","log-CPM","Spearman",
#                                    "SPARK","GIOTTO","MERINGUE","GraphR","SpaceDecorr")))
#
# resultsCorr <- resultsCorr %>% filter(!model %in% c("Pearson","Spearman"),!is.na(model),
#        m == 0.5, num_cells > 100, num_cells < 10000, range == 0.5,
#       kprop %in% c(200, 1000, 0.1),libsize ==200)
# write.csv(resultsCorr, file = paste0(dir.out,"Pairwise_results.csv"))

resultsCorr <- read.csv(file.path(dir.data, "Pairwise_results.csv")) %>%
    mutate(model = factor(model, levels = c("Pearson","log-CPM","Spearman",
                                            "SPARK","GIOTTO","MERINGUE","GraphR","SpaceDecorr")))

# Main figure -------------------------------------------------------------
# Type 1 error
plott1err <- resultsCorr %>%
  filter(model != "GraphR", rho == 0, sim != "MV", kprop == 200) %>%
  mutate(t1err = ifelse(pvalue < 0.05, 1, 0 )) %>%
  group_by(sim,model, range,rho, num_cells, m, ratio, kprop) %>%
  summarise(t1err = mean(t1err,na.rm = T), time = mean(times,na.rm = T)) %>%
  ggplot(aes(x = num_cells, y = t1err, col = model, group = model)) +
  geom_hline(yintercept = 0.05, col = "darkgrey") +
  geom_point() + geom_line(linewidth = 1) +
  scale_x_continuous(breaks = c(500,seq(1000,10000, by = 1000))) +
  scale_y_continuous(breaks = c(0,0.05,0.25, 0.5,0.75,1), limits = c(0,1)) +
  scale_color_manual(values = mycol)+
  theme_bw() +
  labs(x = "Number of cells",y = "Type 1 error",
       col = "Model",lty = "Basis dimension (K)")+
  theme(legend.position = "top",
        text = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(0.5, "cm"),
        legend.title = element_text(margin = margin(r = 20)),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 1));plott1err

# Power
plotpower <- resultsCorr %>%filter() %>%
                filter(model != "GraphR", num_cells == 2000,
                       rho >= 0, kprop == 200, sim!= "MV") %>%
                mutate(power = ifelse(pvalue < 0.05, 1, 0 ),
                       kprop = ifelse(model != "SpaceDecorr", 0.1,kprop)) %>%
                group_by(sim,model, range,rho, num_cells, m, ratio) %>%
                summarise(power = mean(power, na.rm = T),
                          time = mean(times,na.rm = T)) %>%
                ggplot(aes(x = rho, y = power, col = model, group = model)) +
                geom_hline(yintercept = c(0.95), col = "darkgrey") +
                geom_line(linewidth = 1) +
                scale_color_manual(values = mycol)+
                scale_y_continuous(breaks = c(0,0.25, 0.5,0.75,0.95, 1),
                                   limits = c(0,1)) +
                theme_bw() +
                facet_grid(paste0(num_cells," cells")~ .) +
                labs(x = "True Correlation",y = "Power", col = "",lty = "")+
                theme(legend.position = "top",
                      text = element_text(size = 15),
                      axis.text.x = element_text(size = 10),
                      strip.background = element_rect(fill = NA, colour = NA))+
                guides(color = guide_legend(nrow = 1));plotpower

# Time
plottime <- resultsCorr %>%
  filter(model != "GraphR", rho ==0, kprop == 200, sim != "MV") %>%
  mutate(kprop = ifelse(model != "SpaceDecorr", 0.1, kprop))%>%
  ggplot(aes(x = num_cells, y = times/3, col = model, group = paste(model,kprop))) +
  stat_summary() +
  stat_summary(geom = "line",linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = mycol)+
  scale_x_continuous(breaks = c(500,seq(1000,10000, by = 1000))) +
  scale_y_continuous(trans=scales::log_trans(),
                     breaks = c(0, 0.1,1, 10, 60,60*2,60*5,60*10,60*30),
                     labels = c("0s","0.1s","1s","10s","1m","2m","5m","10m","30m")) +
  labs(x = "Number of cells",y = "Time per gene", col = "",lty = "")+
  theme(legend.position = "top",
        text = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 2))


library(cowplot)
cowplot::plot_grid(legend = get_plot_component(plott1err, 'guide-box-top',
                                               return_all = TRUE),
                   plot_grid(plott1err + theme(legend.position = "none"),
                             plotpower+ theme(legend.position = "none"),
                             plottime+ theme(legend.position = "none"),
                             nrow = 1, labels = c("A","B","C")
                   ),
                   nrow = 2,rel_heights = c(0.15,0.85))


# Supplement -------------------------------------------------------------
# Type 1 error
plott1err_sup <- resultsCorr %>%
  filter(model != "GraphR", rho == 0) %>%
  mutate(t1err = ifelse(pvalue < 0.05, 1, 0 ),
         kprop = factor(kprop, levels = c(0.1,200,1000),
                        labels = c("10%","200","1000")
         )) %>%
  group_by(sim,model, range,rho, num_cells, m, ratio, kprop) %>%
  summarise(t1err = mean(t1err,na.rm = T), time = mean(times,na.rm = T)) %>%
  ggplot(aes(x = num_cells, y = t1err, col = model, group = paste(model,kprop),
             lty = as.factor(kprop))) +
  geom_hline(yintercept = 0.05, col = "darkgrey") +
  geom_point() + geom_line(linewidth = 1) +
  facet_grid(.~sim) +
  scale_x_continuous(breaks = c(500,seq(1000,10000, by = 1000))) +
  scale_y_continuous(breaks = c(0,0.05,0.25, 0.5,0.75,1), limits = c(0,1)) +
  scale_color_manual(values = mycol)+
  theme_bw() +
  labs(x = "Number of cells",y = "Type 1 error",
       col = "Model",lty = "Basis dimension (K)")+
  theme(legend.position = "top",
        text = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        panel.spacing = unit(0.5, "cm"),
        legend.title = element_text(margin = margin(r = 20)),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 1));plott1err_sup

# Power
plotpower_sup <- resultsCorr %>%
                filter(model != "GraphR",num_cells == 2000, rho >= 0) %>%
                mutate(power = ifelse(pvalue < 0.05, 1, 0 ),
                       kprop = ifelse(model != "SpaceDecorr", 0.1,kprop),
                       kprop = factor(kprop, levels = c(0.1,200,1000),
                                      labels = c("10%","200","1000"))) %>%
                group_by(sim,model, range,rho, num_cells, m, ratio, kprop) %>%
                summarise(power = mean(power, na.rm = T),
                          time = mean(times,na.rm = T)) %>%
                ggplot(aes(x = rho, y = power, col = model,
                           group = paste(model, kprop),
                           lty = as.factor(kprop))) +
                geom_hline(yintercept = c(0.95), col = "darkgrey") +
                geom_line(linewidth = 1) +
                scale_color_manual(values = mycol)+
                scale_y_continuous(breaks = c(0,0.25, 0.5,0.75,0.95, 1),
                                   limits = c(0,1)) +
                facet_grid(paste0(num_cells," cells")~ sim, ) +
                theme_bw() +
                labs(x = "True Correlation",y = "Power", col = "",lty = "")+
                theme(legend.position = "top",
                      text = element_text(size = 15),
                      axis.text.x = element_text(size = 10),
                      strip.background = element_rect(fill = NA, colour = NA))+
                guides(color = guide_legend(nrow = 1));plotpower_sup

# Time
plottime_sup <- resultsCorr %>%
  filter(model != "GraphR", rho ==0) %>%
  mutate(kprop = ifelse(model != "SpaceDecorr", 0.1, kprop),
         kprop = factor(kprop, levels = c(0.1,200,1000),
                        labels = c("10%","200","1000")))%>%
  ggplot(aes(x = num_cells, y = times/3, col = model,
             group = paste(model,kprop),lty = as.factor(kprop)
  )) +
  stat_summary() +
  stat_summary(geom = "line",linewidth = 1) +
  theme_bw() +
  scale_color_manual(values = mycol)+
  scale_x_continuous(breaks = c(500,seq(1000,10000, by = 1000))) +
  scale_y_continuous(trans=scales::log_trans(),
                     breaks = c(0, 0.1,1, 10, 60,60*2,60*5,60*10,60*30),
                     labels = c("0s","0.1s","1s","10s","1m","2m","5m","10m","30m")) +
  labs(x = "Number of cells",y = "Time per gene", col = "",lty = "")+
  facet_wrap(""~.)+
  theme(legend.position = "top",
        text = element_text(size = 15),
        axis.text.x = element_text(size = 10),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 2))

cowplot::plot_grid(legend = get_plot_component(plott1err_sup, 'guide-box-top',
                                               return_all = TRUE),
                   plot_grid(plott1err_sup + theme(legend.position = "none"),
                             plotpower_sup + theme(legend.position = "none"),
                             plottime_sup + theme(legend.position = "none"),
                             nrow = 1, labels = c("A","B","C")
                   ),
                   nrow = 2,rel_heights = c(0.15,0.85))


# GraphR ------------------------------------------------------------------
# Type 1 error
plott1err_graphR <- resultsCorr %>%
  filter(model %in% c("GraphR","SpaceDecorr"), rho == 0) %>%
  mutate(t1err = ifelse(pvalue < 0.05, 1, 0 ),
         kprop = factor(kprop, levels = c(0.1,200,1000),
                        labels = c("10%","200","1000"))) %>%
  group_by(sim,model, range,rho, num_cells, m, ratio, kprop) %>%
  summarise(t1err = mean(t1err,na.rm = T), time = mean(times,na.rm = T)) %>%
  ggplot(aes(x = num_cells, y = t1err, col = model, group = paste(model,kprop),
             lty = as.factor(kprop))) +
  geom_hline(yintercept = 0.05, col = "darkgrey") +
  geom_point() + geom_line(linewidth = 1) +
  facet_grid(.~sim) +
  scale_x_continuous(breaks = c(500,seq(1000,10000, by = 1000))) +
  theme_bw() +
  labs(x = "Number of cells",y = "Type 1 error",
       col = "Model",lty = "Basis dimension (K)")+
  theme(legend.position = "top",
        text = element_text(size = 10),
        legend.title = element_text(margin = margin(r = 20)),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 1));plott1err_graphR

# Power
plotpower_graphR <- resultsCorr %>%
  filter(model %in% c("GraphR","SpaceDecorr"), num_cells == 2000, rho >= 0) %>%
  mutate(power = ifelse(pvalue < 0.05, 1, 0 ),
         kprop = factor(kprop, levels = c(0.1,200,1000),
                        labels = c("10%","200","1000"))) %>%
  group_by(sim,model, range,rho, num_cells, m, ratio, kprop,) %>%
  summarise(power = mean(power, na.rm = T),
            time = mean(times,na.rm = T)) %>%
  ggplot(aes(x = rho, y = power, col = model,
             group = paste(model, kprop), lty = as.factor(kprop))) +
  #geom_point() +
  geom_hline(yintercept = c(0.95), col = "darkgrey") +
  geom_line(linewidth = 1) +
  #stat_smooth(se = F, span = 0.5) +
  facet_grid(paste0(num_cells," cells")~ sim, ) +
  theme_bw() +
  labs(x = "True Correlation",y = "Power",
       col = "",lty = "")+
  theme(legend.position = "top",
        text = element_text(size = 10),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 1))

# Time
plottime_graphR <- resultsCorr %>%
  filter(model %in% c("GraphR","SpaceDecorr"), rho ==0) %>%
  filter(!(model == "GraphR" & kprop != 0.1)) %>%
  ggplot(aes(x = num_cells, y = times, col = model,
             group = paste(model, kprop), lty = as.factor(kprop))) +
  #geom_boxplot() +
  stat_summary() +
  stat_summary(geom = "line",linewidth = 1) +
  theme_bw() +
  scale_x_continuous(breaks = c(500,seq(1000,10000, by = 1000))) +
  scale_y_continuous(trans=scales::log_trans(),
                     breaks = c(0, 0.1,1, 10, 60,60*2,60*5,60*10,60*30),
                     labels = c("0s","0.1s","1s","10s","1m","2m","5m","10m","30m")) +
  labs(x = "Number of cells",y = "Time", col = "",lty = "")+
  facet_wrap(""~.)+
  theme(legend.position = "top",
        text = element_text(size = 10),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 2))


cowplot::plot_grid(legend = get_plot_component(plott1err_graphR, 'guide-box-top',
                                               return_all = TRUE),
                   plot_grid(plott1err_graphR + theme(legend.position = "none"),
                             plotpower_graphR + theme(legend.position = "none"),
                             plottime_graphR + theme(legend.position = "none"),
                             nrow = 1, labels = c("A","B","C"), rel_widths = c(2,2,1)
                   ),
                   nrow = 2,rel_heights = c(0.15,0.85))

# Other correlations ------------------------------------------------------
# filepath <- paste0(dir.out,"Partial/Pairwise_robust")
# files <- list.files(filepath,pattern = '\\.csv')
# resultsCorrothers0 <- do.call(rbind, lapply(files, function(file)data.table::fread(paste0(filepath,file))))[,-1]
# resultsCorrothers <- resultsCorrothers0 %>%
#   filter(normalization == "spacedecorr",
#          method %in% c("pearson","spearman", "percentbend(tr 0.1)"),
#          num_cells == 2000,kprop == 0.1) %>%
#   mutate(method = factor(method, levels = c("pearson","spearman","percentbend(tr 0.1)"),
#                   labels = c("Pearson","Spearman","Percent bend")) )
# write.csv(resultsCorrothers, file = paste0(dir.out,"Pairwise_results_robust.csv"))
resultsCorrothers <- read.csv(file.path(dir.data, "Pairwise_results_robust.csv"))

# Power
resultsCorrothers %>%
    group_by(method, normalization, sim, num_cells, range, rho, m, ratio, libsize, kprop) %>%
    mutate(power = ifelse(pvalue < 0.05, 1, 0 )) %>%
    summarise(power = mean(power, na.rm = T)) %>%
    ggplot(aes(x = rho, y = power, col = method)) +
    geom_hline(yintercept = c(0.95), col = "darkgrey") +
    geom_line(linewidth = 1) + geom_vline(xintercept = 0) +
    facet_grid(sim~paste("m:",m)) +
    theme_bw() +
    xlim(-0.5,0.5)+
    labs(x = "True Correlation",y = "Power", col = "",lty = "k prop")+
    theme(legend.position = "top",
          text = element_text(size = 10),
          strip.background = element_rect(fill = NA, colour = NA))+
    guides(color = guide_legend(nrow = 1))

