#------------------------------------------------#
#------- Results Simulation DiffNet Niches ------#
#------------------------------------------------#

library(tidyverse)
library(cowplot)
library(pROC)

# Set directories  --------------------------------------------------------
dir.data <- "/path/to/zenodo/download/"
mycol <- c("#E41A1C", "#377EB8","#4DAF4A","#FF7F00","#984EA3"  )

# NOTE: The CSV files read below were pre-generated from parallel simulation outputs.
# The code used to generate them is retained here for reference but is commented out,
# as it requires the full set of partial results files which are not included in the repo.

# Results -----------------------------------------------------------------
# filepath <- paste0(dir.out,"SimulationDiffNet/DNniche/")
# files <- list.files(filepath,pattern = '\\.csv')
# resultsdiffnetNichep0 <- do.call(rbind, lapply(files, function(file)data.table::fread(paste0(filepath,file))))[,-1]%>%
#   mutate(Model = factor(Model, levels = c("logPearson","SPARK","Giotto","Meringue","SpaceDecorr"),
#                         labels = c("log-CPM","SPARK","GIOTTO","MERINGUE","SpaceDecorr"))) %>%
#   filter(!is.na(Model),num_genes == 20, scenario == 2, kprop == 0.1,perturb_prob == 0.25)
# write.csv(resultsdiffnetNichep0, file = paste0(dir.out,"DiffNet_niches_results.csv"))

resultsdiffnetNichep0 <- read.csv(paste0(dir.data,"DiffNet_niches_results.csv"))
colnames(resultsdiffnetNichep0)[3:8] <- c("trueN0","trueN1","trueDiff",
                                          "pvalueN0","pvalueN1","pvalueDiff")

resultsdiffnetNichep1 <- resultsdiffnetNichep0 %>%
  mutate(row_id = row_number()) %>%
  pivot_longer(cols = starts_with("true"), names_to = "true_name", values_to = "true") %>%
  pivot_longer(cols = starts_with("pvalue"), names_to = "p_name", values_to = "pvalue") %>%
  filter(gsub("true","",true_name) == gsub("pvalue","",p_name)) %>%
  mutate(graph = gsub("true","",true_name)) %>%
  dplyr::select(!c(row_id, true_name,p_name))

## F1 ----
resultsdiffnetNichepMetrics <- resultsdiffnetNichep1 %>%
  mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue),
         pred = ifelse(pvalue < 0.05, 1, 0)) %>%
  group_by(Model, type, rept, num_genes, scenario, kprop, perturb_prob, graph) %>%
  summarise(
    TP = sum(pred == 1 & true == 1),
    FP = sum(pred == 1 & true == 0),
    FN = sum(pred == 0 & true == 1),
    precision = ifelse((TP + FP) == 0, 0, TP / (TP + FP)),
    recall = ifelse((TP + FN) == 0, 0, TP / (TP + FN)),
    F1 = ifelse((precision + recall) == 0, 0, 2 * precision * recall / (precision + recall)),
    .groups = "drop"
  )


plot_f1_niche <- resultsdiffnetNichepMetrics %>%
  mutate(graph = factor(graph, levels = c("N1","N0",  "Diff"),
                        labels = c("Myeloid-enriched\nstroma (n=4,594)","Stroma\n(n=1,281)","Difference"))) %>%
  ggplot(aes(x = graph, y = F1, col = Model, lty = type)) +
  geom_boxplot()+
  scale_fill_manual(values = c("white","gray")) +
  theme_bw() +
  scale_color_manual(values = mycol)+
  labs(x = "Network", col = "",y="F1", fill = "") +
  theme(legend.position = "top",
        text = element_text(size = 15),
        legend.text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 11),
        strip.background = element_rect(fill = NA, colour = NA))+
  guides(color = guide_legend(nrow = 1));plot_f1_niche


## roc ----
pval_grid <- c(10e-20,10e-10,10e-5,10e-2,10e-1,0.05, 0.10, 1)
resultsdiffnetNichepROC <- resultsdiffnetNichep1 %>%
  mutate(pvalue = ifelse(is.na(pvalue), 1, pvalue)) %>%
  group_by(Model, type, rept, num_genes, scenario, kprop, perturb_prob, graph) %>%
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


point_data <- resultsdiffnetNichepROC %>%
  filter(threshold==0.05) %>%
  mutate(graph = factor(graph, levels = c("N1","N0",  "Diff"),
                        labels = c("Myeloid-enriched\nstroma (n=4,594)",
                                   "Stroma\n(n=1,281)",
                                   "Difference"))) %>%
  #filter(!(type == "twocor" & Model != "SpaceDecorr")) %>%
  mutate(type = factor(type, levels = c("Fisher", "twocor"),
                       labels = c("Standard", "Robust"))) %>%
  group_by(Model, type, num_genes, scenario, kprop, perturb_prob, graph) %>%
  summarize(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    lower = quantile(precision, 0.025, na.rm = TRUE),
    upper = quantile(precision, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

plot_roc_niche <- resultsdiffnetNichepROC %>%
   mutate(graph = factor(graph, levels = c("N1","N0",  "Diff"),
                        labels = c("Myeloid-enriched\nstroma (n=4,594)",
                                   "Stroma\n(n=1,281)",
                                   "Difference"))) %>%
  mutate(type = factor(type, levels = c("Fisher", "twocor"),
                       labels = c("Standard", "Robust"))) %>%
  group_by(Model, type, num_genes, scenario, kprop, perturb_prob, graph, threshold) %>%
  summarize(
    mean_precision = mean(precision, na.rm = TRUE),
    mean_recall = mean(recall, na.rm = TRUE),
    lower = quantile(precision, 0.025, na.rm = TRUE),
    upper = quantile(precision, 0.975, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = mean_recall, y = mean_precision, color = Model, lty = type, pch = type)) +
  geom_line() +
  scale_color_manual(values = mycol)+
  scale_x_continuous("Recall", limits = c(0, 1), breaks = c(0,0.25,.5,.75,1),
                     labels = c("0","0.25","0.5","0.75","1")) +
  scale_y_continuous("Precision", limits = c(0, 1),breaks = c(0,0.25,.5,.75,1),
                     labels = c("0","0.25","0.5","0.75","1")) +
  geom_point(data = point_data,
             size = 2.5) +
  labs(col = "", lty = "Correlation test", pch = "Correlation test") +
  facet_grid(. ~ graph) +
  theme_bw() +
  guides(color = guide_legend(nrow = 1)) +
  theme(text = element_text(size = 15),
        panel.grid.minor = element_blank(),
        axis.text = element_text(size = 9),
        legend.position = "top", strip.background = element_blank(),
        legend.title = element_text(margin = margin(r = 20)),
        legend.key.width = unit(2, "line"));plot_roc_niche

plot_grid(legend = get_plot_component(plot_roc_niche, 'guide-box-top', return_all = TRUE),
          plot_grid(plot_f1_niche + theme(legend.position = "none"),
                    plot_roc_niche+ theme(legend.position = "none"),
                    align = "h", axis = "tb",
                    nrow = 1, rel_widths = c(1,1) #,labels = c("i","j","k")
          ),
          nrow = 2,rel_heights = c(0.15,0.85))

