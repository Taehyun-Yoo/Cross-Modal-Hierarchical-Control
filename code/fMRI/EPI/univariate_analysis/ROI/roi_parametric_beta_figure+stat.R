library(readxl)
library("tidyverse")
library("ggpubr")
library("grid")
library("rstatix")
library("ggplot2")
library("MBESS")
library(ggprism)
library(svglite)

# set working directory
dir <- ("/data/Project1/MSHP/analysis/")
result_dir <- paste0(dir, 'results_2nd/ROI_analysis/ROI_2voxels/')

#set roi
roi <- c("R", "F", "D", "D(2)")
#roi_R <- c("[-27, -4, 53]")
#roi_F <- c("[-39, 14, 29]")
#roi_D <- c("[-45, 35, 17]")

#title <- c(roi_R, roi_F, roi_D)

hier <- c("R", "F", "D")

# set subject
sub_num = 42;
sub_number <- c(01, 02, 03, 04, 05, 06, 07, 09, 10, 12, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 48, 49, 51, 52, 53)
sub <- sprintf("%02d", sub_number)

data_raw <- data.frame()
for (i in 1:length(sub)) {
  sub_dir <- paste0(dir, 'sub', sub[i], '/analysis/ROI_analysis/')
  for (h in 1:length(hier)) {
    data_sub <- read.delim(paste0(sub_dir, 'ROI_2voxels/', hier[h], '_parametric/', hier[h], '_parametric_beta_sub', sub[i], '.txt'), sep="")
    data_sub$Subject = sub[i]
    data_sub <- data_sub[c("Subject", "Session", "Condition", "ROI", "beta")]
    data_raw <- rbind(data_raw, data_sub)
  }
}

roi_total <- data_raw %>%
  group_by(Subject, Condition, ROI) %>%
  summarise(beta = mean(beta)) %>%
  ungroup() %>%
  mutate(Hier = substring(Condition, 23, 23), .before = Condition) %>%
  mutate(Hier = factor(Hier, levels = hier)) %>%
  mutate(ROI = factor(ROI, levels = roi))

# roi_total <- roi %>%
#   map(~ sprintf("/home/taehyun/Desktop/MSHP/marsbar/ROI/final/final_roi_analysis_hierarchy/%s", .x)) %>%
#   map(list.files, pattern = "*.xlsx", full.name = TRUE, recursive = TRUE) %>%
#   flatten() %>%
#   map_dfr(readxl::read_xlsx, sheet = 1) %>%
#   select(-`...1`) %>%
#   mutate(ROI = fct_inorder(ROI),
#          Hier = fct_inorder(Hier),
#          Comp = fct_inorder(Comp),
#          Condition = fct_inorder(Condition))

#One_way_ANOVA_ROI (each hierarchy)
ROI_one_way_ANOVA <- function(res, r) {
  one_way_ANOVA <- anova_test(data = subset(res, res$ROI == roi[r]), dv = beta, wid = Subject, within = Hier, effect.size = "pes")
  table <- get_anova_table(one_way_ANOVA)
  effect <- c(table[1,"Effect"])
  F_value <- c(paste0("F(", table[1,"DFn"], ",", table[1,"DFd"], ") = ", table[1,"F"]))
  p_value <- c(table[1,"p"])
  p_value_star <- c(symnum(table[1,"p"], cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))
  pes <- c(table[1,"pes"])
  # lims <- conf.limits.ncf(F.value = table[1, "F"], df.1 = table[1, "DFn"], df.2 = table[1, "DFd"], conf.level = .90)
  # ci <- c(lims$Lower.Limit/(lims$Lower.Limit + table[1, "DFn"] + table[1, "DFd"] + 1), lims$Upper.Limit/(lims$Upper.Limit + table[1, "DFn"] + table[1, "DFd"] + 1))
  # ci <- round(ci, 3)
  # ci <- paste0("[", ci[1], ", ", ci[2], "]")
  one_way_ANOVA_ROI_table <- data.frame(effect, F_value, p_value, p_value_star, pes)
  one_way_ANOVA_ROI_table
}

write_csv(ROI_one_way_ANOVA(roi_total, 1), paste0(result_dir, 'parametric/roi_R/beta_one_way_ANOVA_in_roi_R.csv'))
write_csv(ROI_one_way_ANOVA(roi_total, 2), paste0(result_dir, 'parametric/roi_F/beta_one_way_ANOVA_in_roi_F.csv'))
write_csv(ROI_one_way_ANOVA(roi_total, 3), paste0(result_dir, 'parametric/roi_D/beta_one_way_ANOVA_in_roi_D.csv'))
write_csv(ROI_one_way_ANOVA(roi_total, 4), paste0(result_dir, 'parametric/roi_D(2)/beta_one_way_ANOVA_in_roi_D(2).csv'))

beta_posthoc <- function(res, r) {
  post_t_test <- t_test(beta ~ Hier, p.adjust.method = "holm", paired = TRUE, data = subset(res, res$ROI == roi[r]), detailed = TRUE)
  .y. <- c("beta_mean", "beta_mean", "beta_mean")
  group1 <- c(post_t_test$group1)
  group2 <- c(post_t_test$group2)
  comparisons <- c(paste0(post_t_test$group1[1], "-", post_t_test$group2[1]), paste0(post_t_test$group1[2], "-", post_t_test$group2[2]), paste0(post_t_test$group1[3], "-", post_t_test$group2[3]))
  t_value <- c(paste0("t(", post_t_test$df[1], ") = ", round(post_t_test$statistic[1],3)), paste0("t(", post_t_test$df[2], ") = ", round(post_t_test$statistic[2],3)), paste0("t(", post_t_test$df[3], ") = ", round(post_t_test$statistic[3],3)))
  adjusted_p_value <- c(post_t_test$p.adj[1], post_t_test$p.adj[2], post_t_test$p.adj[3])
  adjusted_p_value_star <- c(symnum(post_t_test$p.adj[1], cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), symnum(post_t_test$p.adj[2], cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")), symnum(post_t_test$p.adj[3], cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))
  d <- c(round(abs(post_t_test$statistic[1]/sqrt(post_t_test$n1[1])),3), round(abs(post_t_test$statistic[2]/sqrt(post_t_test$n1[2])),3), round(abs(post_t_test$statistic[3]/sqrt(post_t_test$n1[3])),3))
  ci <- c(paste0("[", round(post_t_test$conf.low[1],3), ", ", round(post_t_test$conf.high[1],3), "]"), paste0("[", round(post_t_test$conf.low[2],3), ", ", round(post_t_test$conf.high[2],3), "]"), paste0("[", round(post_t_test$conf.low[3],3), ", ", round(post_t_test$conf.high[3],3), "]"))
  y.position <- c(0.09, 0.12, 0.105)
  post_beta_table <- data.frame(.y., group1, group2, comparisons, t_value, adjusted_p_value, adjusted_p_value_star, d, ci, y.position)
  post_beta_table
}

write_csv(beta_posthoc(roi_total, 1), paste0(result_dir, 'parametric/roi_R/beta_posthoc_in_roi_R.csv'))
write_csv(beta_posthoc(roi_total, 2), paste0(result_dir, 'parametric/roi_F/beta_posthoc_in_roi_F.csv'))
write_csv(beta_posthoc(roi_total, 3), paste0(result_dir, 'parametric/roi_D/beta_posthoc_in_roi_D.csv'))
write_csv(beta_posthoc(roi_total, 4), paste0(result_dir, 'parametric/roi_D(2)/beta_posthoc_in_roi_D(2).csv'))

hier_names <- c(
  'R'="Response",
  'F'="Feature",
  'D'="Dimension"
)

plot_beta <- function(res, r) {
  
  roi_beta <- list()
  for (i in 1:length(roi)) {
    roi_beta[[i]] = filter(res, ROI == roi[i])
  }
  
  summary_beta <- roi_beta[[r]] %>%
    group_by(Hier) %>%
    summarise(se = sd(beta)/sqrt(n()),
              beta = mean(beta)) %>%
    ungroup()
  
  my_plot_beta <- ggplot(summary_beta,aes(x=Hier,y=beta)) +
    geom_bar(aes(fill = Hier), stat="identity", position="dodge", width = 0.6)+
    geom_errorbar(aes(ymin=beta-se,ymax=beta+se),
                  position=position_dodge(0.9),width=0.1, alpha = 1)+
    #facet_grid(. ~Hier) +
    scale_fill_manual(values = c("#2B6A6C", "#F29724", "#B80D48")) +
    scale_x_discrete(limits=c("R", "F", "D"))+
    coord_cartesian(ylim=c(0, 0.13)) +
    scale_y_continuous(breaks = seq(0, 0.1, 0.05))+
    labs(y = "Parameter estimates") +
    # theme(panel.grid.major=element_line(colour="Grey", linetype="dotted"), panel.grid.minor=element_blank(), panel.background=element_rect(fill="transparent",colour=NA), 
    #       axis.line = element_line(colour="Black"), 
    #       plot.title=element_text(hjust = 0.5, face="bold", size=20, vjust=2),
    #       #axis.title.x = element_text(size=20),
    #       axis.title.x = element_blank(),
    #       axis.text = element_text(size=18, colour = "Black"),
    #       axis.text.y = element_text(size=15, colour = "Black"),
    #       axis.title.y = element_text(size=20),
    #       strip.text = element_text(size = 9, colour = "white", face = "bold"),
    #       legend.position = "none")
    theme(panel.grid.major=element_line(colour="Grey", linetype="dotted"), panel.grid.minor=element_blank(), panel.background=element_rect(fill="transparent",colour=NA), 
          axis.line = element_line(colour="Black"), 
          #plot.title=element_text(hjust = 0.5, face="bold", size=20, vjust=2),
          #axis.title.x = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text.x = element_text(size=20, colour = "Black"),
          axis.text.y = element_text(size=20, colour = "Black"),
          axis.title.y = element_text(size=25),
          #strip.text = element_text(size = 9, colour = "white", face = "bold"),
          legend.position = "none")
  
  total_K_beta <- function(res, r){
    K <- beta_posthoc(res, r) 
    K <- tibble(K)
    K <- subset(K, K$adjusted_p_value_star != "" & K$adjusted_p_value_star != "ns")
    K
  }
  
  if (nrow(total_K_beta(res, r)) != 0){
    my_plot_beta <- my_plot_beta +
      # add_pvalue(total_K_beta(res, r),
      #            xmin = "group1", xmax = "group2",
      #            label = "adjusted_p_value_star",
      #            y.position = "y.position",
      #            bracket.size = 0.8,
      #            label.size = 7)
      add_pvalue(total_K_beta(res, r),
                 xmin = "group1", xmax = "group2",
                 label = "adjusted_p_value_star",
                 y.position = "y.position",
                 bracket.size = 0.8,
                 label.size = 9)
  } else {
    my_plot_beta <- my_plot_beta
  }
  
  my_plot_beta
  figure_dir <- paste0(result_dir, 'parametric/roi_', roi[r], '/')
  ggsave(paste0(figure_dir,"beta_roi_", roi[r], "_final.svg"), width = 90, height = 90, units = "mm")
}

plot_beta(roi_total, 1)
plot_beta(roi_total, 2)
plot_beta(roi_total, 3)
plot_beta(roi_total, 4)
