#-------------------------------------------------#
#------- Results Simulation DiffNet Samples ------#
#-------------------------------------------------#

library(tidyverse)
library(cowplot)

# Set directories  --------------------------------------------------------
dir.data <- "/path/to/zenodo/download/"
mycol <- c("#E41A1C", "#377EB8","#4DAF4A","#FF7F00","#984EA3"  )

# NOTE: The CSV files read below were pre-generated from parallel simulation outputs.
# The code used to generate them is retained here for reference but is commented out,
# as it requires the full set of partial results files which are not included in the repo.

# Results ------
# filepath <- paste0(dir.out,"SimulationDiffNet/DNsample/")
# files <- list.files(filepath,pattern = '\\.csv')
# my_colors <- c( "#5F4690", "#1D6996", "#38A6A5", "#0F8554", "#73AF48", "#EDAD08", "#E17C05",
#                 "#CC503E", "#94346E", "#6F4070", "#994E95","#666666")
#
# resultsdiffnetSamp0 <- do.call(rbind, lapply(files, function(file)data.table::fread(paste0(filepath,file))))[,-1] %>%
#   mutate(Model = factor(Model, levels = c("logPearson","SPARK","Giotto","Meringue","SpaceDecorr"),
#                         labels = c("logPearson","SPARK","GIOTTO","MERINGUE","SpaceDecorr")))%>%
#   filter(!is.na(Model),num_genes == 20, num_cells1 == 2000, range == 0.25, kprop == 0.1)
# write.csv(resultsdiffnetSamp0, file = paste0(dir.out,"DiffNet_samples_results.csv"))
resultsdiffnetSamp0 <- read.csv(paste0(dir.data,"DiffNet_samples_results.csv"))

colnames(resultsdiffnetSamp0)[3:8] <- c("trueT1","trueT2","trueDiff",
                                        "pvalueT1","pvalueT2","pvalueDiff")

library(pROC)
resultsdiffnetSamp1 <- resultsdiffnetSamp0 %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = starts_with("true"), names_to = "true_name", values_to = "true") %>%
  pivot_longer(cols = starts_with("pvalue"), names_to = "p_name", values_to = "pvalue") %>%
  filter(gsub("true","",true_name) == gsub("pvalue","",p_name)) %>%
  mutate(graph = gsub("true","",true_name))

## F1 ----
resultsdiffnetSamppMetrics <- resultsdiffnetSamp1 %>%
  mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue),
         pred = ifelse(pvalue < 0.05, 1, 0)) %>%
  group_by(Model, type, rept, num_genes, num_cells1, range,sim, kprop, perturb_prob, graph) %>%
  summarise(
    TP = sum(pred == 1 & true == 1),
    FP = sum(pred == 1 & true == 0),
    FN = sum(pred == 0 & true == 1),
    precision = ifelse((TP + FP) == 0, 0, TP / (TP + FP)),
    recall = ifelse((TP + FN) == 0, 0, TP / (TP + FN)),
    F1 = ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall)),
    .groups = "drop"
  )

plot_f1_samp <- resultsdiffnetSamppMetrics %>%
  mutate(graph = ifelse(graph == "Diff", "Difference","Tissues"),
         graph = factor(graph, levels = c("Tissues","Difference"))) %>%
  filter(!(type =="twocor" & Model != "SpaceDecorr")) %>%
  ggplot(aes(x = graph, y = F1, col =Model, lty = type)) +
  geom_boxplot()+
  #stat_summary() + stat_summary(geom = "line") +
  facet_grid(.~sim) +
  scale_color_manual(values = mycol)+
  scale_fill_manual(values = c("white","gray")) +
  theme_bw() + ylim(0,1)+
  labs(x = "Network", col = "",y="", lty = "Correlation test") +
  theme(text = element_text(size = 15),
        #axis.text.x = element_text(angle = 45,  hjust=1),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11),
        legend.text = element_text(size = 15),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(
    color = guide_legend(nrow = 1,position = "top")
  );plot_f1_samp

# roc ----
pval_grid <- c(10e-20,10e-10,10e-5,10e-2,10e-1,0.05, 0.10, 1)
resultsdiffnetSamppROC <- resultsdiffnetSamp1 %>%
  mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue)) %>%
  group_by(Model, type, sim, num_genes, num_cells1, kprop, perturb_prob,range,graph,rept) %>%
  do({
    df <- .
    map_dfr(pval_grid, function(thresh) {
      pred <- ifelse(df$pvalue < thresh, 1, 0)
      TP <- sum(pred == 1 & df$true == 1)
      FP <- sum(pred == 1 & df$true == 0)
      FN <- sum(pred == 0 & df$true == 1)

      precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
      recall <- if ((TP + FN) == 0) NA else TP / (TP + FN)

      tibble(
        threshold = thresh,
        precision = precision,
        recall = recall
      )
    })
  })


point_data_samp <- resultsdiffnetSamppROC %>%
  filter(threshold == 0.05) %>%
  mutate(graph = ifelse(graph == "Diff", "Difference","Tissues"),
         graph = factor(graph, levels = c("Tissues","Difference"))) %>%
  filter(!(type == "twocor" & Model != "SpaceDecorr")) %>%
  mutate(type = factor(type, levels = c("Fisher", "twocor"),
                       labels = c("Standard", "Robust"))) %>%
  group_by(Model, type, sim, num_genes, num_cells1, kprop, perturb_prob,range,graph) %>%
  summarize(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    lower = quantile(precision, 0.025, na.rm = TRUE),
    upper = quantile(precision, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

plot_roc_samp <- resultsdiffnetSamppROC %>%
  mutate(graph = ifelse(graph == "Diff", "Difference","Tissues"),
         graph = factor(graph, levels = c("Tissues","Difference"))) %>%
  filter(!(type == "twocor" & Model != "SpaceDecorr")) %>%
  mutate(type = factor(type, levels = c("Fisher", "twocor"),
                       labels = c("Standard", "Robust"))) %>%
  group_by(Model, type, sim, num_genes, num_cells1, kprop, perturb_prob,range,graph,threshold) %>%
  summarize(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    lower = quantile(precision, 0.025, na.rm = TRUE),
    upper = quantile(precision, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = mean_recall, y = mean_precision, color = Model, lty = type, pch = type)) +
  geom_line() +
  scale_x_continuous("Recall", limits = c(0, 1)) +
  scale_y_continuous("Precision", limits = c(0, 1)) +
  geom_point(data = point_data_samp,size = 2.5) +
  labs(col = "", lty = "Correlation\ntest", pch = "Correlation\ntest") +
  facet_grid(sim ~ graph) +
  theme_bw() +
  scale_color_manual(values = mycol)+
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "top", strip.background = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 11),
        text = element_text(size = 15),
        legend.key.width = unit(2, "line"));plot_roc_samp

plot_grid(legend = get_plot_component(plot_roc_samp, 'guide-box-top', return_all = TRUE),
          plot_grid(plot_f1_samp + theme(legend.position = "none"),
                    plot_roc_samp+ theme(legend.position = "none"),
                    align = "h", axis = "tb",
                    nrow = 1, rel_widths = c(1,1) #,labels = c("i","j","k")
          ),
          nrow = 2,rel_heights = c(0.15,0.85))
