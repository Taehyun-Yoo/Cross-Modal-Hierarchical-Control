library(data.table)
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
result_dir <- paste0(dir, 'mvpa_2_results_2nd/mvpa_ROI/')
behavior_result_dir <-("/home/taehyun/Desktop/Behavior_MSHP/data/behavior_fMRI")

#set roi
roi <- c("R", "D")
#roi_R <- c("[-27, -4, 53]")
#roi_F <- c("[-39, 14, 29]")
#roi_D <- c("[-45, 35, 17]")

#title <- c(roi_R, roi_F, roi_D)

hier <- c("R", "D")

# set subject
sub_num = 42;
sub_number <- c(01, 02, 03, 04, 05, 06, 07, 09, 10, 12, 14, 15, 16, 17, 18, 19, 21, 22, 23, 24, 25, 26, 27, 28, 31, 32, 33, 35, 36, 37, 40, 41, 42, 43, 44, 45, 46, 48, 49, 51, 52, 53)
sub <- sprintf("%02d", sub_number)

data_raw <- data.frame()
for (i in 1:length(sub)) {
  sub_dir <- paste0(dir, 'sub', sub[i], '/analysis/multivariate_analysis_2/')
  for (h in 1:length(hier)) {
    data_sub <- read.delim(paste0(sub_dir, hier[h], '/between_alternative/mvpa_ROI/', hier[h], '_decoding_accuracy_sub', sub[i], '.txt'), sep="")
    data_sub$Subject = sub[i]
    data_sub$Hier = hier[h]
    data_sub <- data_sub[c("Subject", "Hier", "ROI", "accuracy")]
    data_raw <- rbind(data_raw, data_sub)
  }
}

roi_total <- data_raw %>%
  group_by(Subject, Hier, ROI) %>%
  summarise(accuracy = mean(accuracy)) %>%
  ungroup() %>%
  mutate(ROI = factor(ROI, levels = roi)) %>%
  mutate(Hier = factor(Hier, levels = hier))

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

accuracy_posthoc <- function(res, r) {
  post_t_test <- t_test(accuracy ~ Hier, p.adjust.method = "holm", paired = TRUE, data = subset(res, res$ROI == roi[r]), detailed = TRUE)
  .y. <- c("accuracy")
  group1 <- c(post_t_test$group1)
  group2 <- c(post_t_test$group2)
  comparisons <- c(paste0(post_t_test$group1[1], "-", post_t_test$group2[1]))
  t_value <- c(paste0("t(", post_t_test$df[1], ") = ", round(post_t_test$statistic[1],3)))
  p_value <- c(post_t_test$p[1])
  p_value_star <- c(symnum(post_t_test$p[1], cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")))
  y.position <- c(34)
  post_beta_table <- data.frame(.y., group1, group2, comparisons, t_value, p_value, p_value_star, y.position)
  post_beta_table
}

write_csv(accuracy_posthoc(roi_total, 1), paste0(result_dir, 'roi_R/accuracy_posthoc_in_roi_R.csv'))
write_csv(accuracy_posthoc(roi_total, 2), paste0(result_dir, 'roi_D/accuracy_posthoc_in_roi_D.csv'))

hier_names <- c(
  'R'="Response",
  'D'="Dimension"
)

plot_accuracy <- function(res, r) {
  
  roi_accuracy <- list()
  for (i in 1:length(roi)) {
    roi_accuracy[[i]] = filter(res, ROI == roi[i])
  }
  
  summary_accuracy <- roi_accuracy[[r]] %>%
    group_by(Hier) %>%
    summarise(se = sd(accuracy)/sqrt(n()),
              accuracy = mean(accuracy)) %>%
    ungroup()
  
  one_sample_t_test <- subset(res, res$ROI == roi[r]) %>%
    group_by(Hier) %>%
    t_test(accuracy ~ 1, mu = 25) %>%
    adjust_pvalue() %>%
    mutate(y.position = 32) %>%
    mutate(p.signif = case_when(0 < p & p < 0.0001 ~ "****",
                                0.0001 < p & p < 0.001 ~ "***",
                                0.001 < p & p < 0.01 ~ "**",
                                0.01 < p & p < 0.05 ~ "*",
                                0.05 < p ~ ""))
  
  my_plot_accuracy <- ggplot(summary_accuracy,aes(x=Hier,y=accuracy)) +
    geom_bar(aes(fill = Hier), stat="identity", position="dodge", width = 0.6)+
    geom_errorbar(aes(ymin=accuracy-se,ymax=accuracy+se),
                  position=position_dodge(0.9),width=0.1, alpha = 1)+
    #facet_grid(. ~Hier) +
    scale_fill_manual(values = c("#2B6A6C", "#B80D48")) +
    scale_x_discrete(limits=c("R", "D"))+
    coord_cartesian(ylim=c(20, 35)) +
    geom_hline(yintercept = 25, col = "black", linewidth=0.2, lty=2) +
    scale_y_continuous(breaks = seq(20, 35, 5))+
    labs(y = "Decoding accuracy") +
    theme(panel.grid.major=element_line(colour="Grey", linetype="dotted"), panel.grid.minor=element_blank(), panel.background=element_rect(fill="transparent",colour=NA), 
          axis.line = element_line(colour="Black"), 
          plot.title=element_text(hjust = 0.5, face="bold", size=20, vjust=2),
          #axis.title.x = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text = element_text(size=18, colour = "Black"),
          axis.text.y = element_text(size=15, colour = "Black"),
          axis.title.y = element_text(size=20),
          strip.text = element_text(size = 9, colour = "white", face = "bold"),
          legend.position = "none") +
          stat_pvalue_manual(one_sample_t_test, label = "p.signif", xmin = "Hier", xmax = NULL, size = 10)
  
  total_K_accuracy <- function(res, r){
    K <- accuracy_posthoc(res, r) 
    K <- tibble(K)
    K <- subset(K, K$p_value_star != "" & K$p_value_star != "ns")
    K
  }
  
  if (nrow(total_K_accuracy(res, r)) != 0){
    my_plot_accuracy <- my_plot_accuracy +
      add_pvalue(total_K_accuracy(res, r),
                 xmin = "group1", xmax = "group2",
                 label = "p_value_star",
                 y.position = "y.position",
                 bracket.size = 0.8,
                 label.size = 10)
  } else {
    my_plot_accuracy <- my_plot_accuracy
  }
  
  my_plot_accuracy
  figure_dir <- paste0(result_dir, 'roi_', roi[r], '/')
  ggsave(paste0(figure_dir,"accuracy_roi_", roi[r], "_final.svg"), width = 120, height = 120, units = "mm")
}

plot_accuracy(roi_total, 1)
plot_accuracy(roi_total, 2)

#######Behavior#######

# behavior_data_total <- data.frame()
# 
# for (s in 1:sub_num){
#   data_R <- read_excel(paste0(behavior_result_dir, '/', sub[s], '/', 'data', sub[s], '_R.xlsx'))
#   data_D <- read_excel(paste0(behavior_result_dir, '/', sub[s], '/', 'data', sub[s], '_D.xlsx'))
#   
#   data_R$sheet = 1
#   data_D$sheet = 2
#   
#   data_RD <- rbind(data_R, data_D)
#   
#   subtmp1=subset(data_RD,sheet==1, select = c(Subject, Onset, Code, Acc))
#   subtmp1 = subtmp1 %>% filter(Code != 99)
#   sub_R <- as.data.table(subtmp1)
#   sub_R[, inter_RT :=   Onset - shift(Onset, fill=last(Onset))]
#   RT=subset(sub_R, select = c(inter_RT))
#   RT <- RT[-1,]
#   sub_R <- head(sub_R, -1)
#   sub_R[, RT :=  RT]
#   
#   t = as.numeric(substr(sub_R$Code, 9,11))
#   sub_R[, trial_number :=  t]
#   
#   for (i in 1:64){
#     sub_R$condition[sub_R$trial_number == i] = "R1"
#   }
#   subset_R1 <- subset(sub_R, sub_R$condition == "R1")
#   subset_R1 <- as.data.frame(subset_R1)
#   subset_R1['Subject'] <- c(rep(sprintf("%s", sub[s]), 64))
#   subset_R1['Hier'] <- c(rep("R", times = 64))
#   subset_R1['Comp'] <- c(rep("1", times = 64))
#   subset_R1['Condition'] <- c(rep("R1", times = 64))
#   
#   for (i in 129:192){
#     sub_R$condition[sub_R$trial_number == i] = "R4"
#   }
#   subset_R4 <- subset(sub_R, sub_R$condition == "R4")
#   subset_R4 <- as.data.frame(subset_R4)
#   subset_R4['Subject'] <- c(rep(sprintf("%s", sub[s]), 64))
#   subset_R4['Hier'] <- c(rep("R", times = 64))
#   subset_R4['Comp'] <- c(rep("4", times = 64))
#   subset_R4['Condition'] <- c(rep("R4", times = 64))
#   
#   subset_R <- rbind(subset_R1, subset_R4)
#   
#   subtmp2=subset(data_RD,sheet==2, select = c(Subject, Onset, Code, Acc))
#   subtmp2 = subtmp2 %>% filter(Code != 99)
#   sub_D <- as.data.table(subtmp2)
#   sub_D[, inter_RT :=   Onset - shift(Onset, fill=last(Onset))]
#   RT=subset(sub_D, select = c(inter_RT))
#   RT <- RT[-1,]
#   sub_D <- head(sub_D, -1)
#   sub_D[, RT :=  RT]
#   
#   t = as.numeric(substr(sub_D$Code, 9,11))
#   sub_D[, trial_number :=  t]
#   
#   for (i in 1:64){
#     sub_D$condition[sub_D$trial_number == i] = "D1"
#   }
#   subset_D1 <- subset(sub_D, sub_D$condition == "D1")
#   subset_D1 <- as.data.frame(subset_D1)
#   subset_D1['Subject'] <- c(rep(sprintf("%s", sub[s]), 64))
#   subset_D1['Hier'] <- c(rep("D", times = 64))
#   subset_D1['Comp'] <- c(rep("1", times = 64))
#   subset_D1['Condition'] <- c(rep("D1", times = 64))
#   
#   for (i in 129:192){
#     sub_D$condition[sub_D$trial_number == i] = "D4"
#   }
#   subset_D4 <- subset(sub_D, sub_D$condition == "D4")
#   subset_D4 <- as.data.frame(subset_D4)
#   subset_D4['Subject'] <- c(rep(sprintf("%s", sub[s]), 64))
#   subset_D4['Hier'] <- c(rep("D", times = 64))
#   subset_D4['Comp'] <- c(rep("4", times = 64))
#   subset_D4['Condition'] <- c(rep("D4", times = 64))
#   
#   subset_D <- rbind(subset_D1, subset_D4)
#   
#   data_sub <- rbind(subset_R, subset_D)
#   data_sub <- subset(data_sub, select = c(Subject, Hier, Comp, Condition, Acc, RT))
#   
#   data_sub <- mutate(data_sub, Acc = if_else(Acc == "hit", 100, 0))
#   
#   for (k in 1:length(data_sub$RT)){
#     data_sub$RT[k] = data_sub$RT[k]/10
#   }
#   
#   behavior_data_total <- rbind(behavior_data_total, data_sub)
#   
# }
# 
# behav_data_RT = 
#   # use data frame
#   behavior_data_total %>%
#   filter(Acc == 100) %>%
#   # to filter RT > mean(RT) + 3*sd(RT)
#   group_by(Subject, Hier, Comp, Condition) %>%
#   summarise(mean_rt = median(RT)) %>% ungroup()
# 
# behav_data_RT$Hier <- factor(behav_data_RT$Hier, levels = c("R", "D"))
# behav_data_RT$Comp <- factor(behav_data_RT$Comp, levels = c("1", "4"))
# 
# data_total_RT <- roi_total %>% left_join(subset(behav_data_RT, Condition == "R4" | Condition == "D4"))
# 
# correlation_plot_R_RT <- data_total_RT %>%
#   filter(Hier == "R" & ROI == "R") %>%
#   ggplot(aes(x = accuracy, y = mean_rt)) +
#   xlab("Decoding accuracy") +
#   ylab("Mean reaction time (ms)") +
#   coord_cartesian(xlim=c(15,45)) +
#   scale_x_continuous(breaks = c(15, 30, 45)) +
#   geom_point(aes(colour = factor(Subject))) +
#   geom_smooth(method = 'lm', formula = y~x, se = TRUE) +
#   ggpubr::stat_cor(label.y.npc = 0.7, label.x.npc = 0.7) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20, face="bold.italic"),
#     axis.title.x = element_text(size=16, face="bold", family = "Arial Narrow"),
#     axis.title.y = element_text(size=16, face="bold", family = "Arial Narrow"),
#     axis.text.x =  element_text(size = 12, face = "bold", colour = "black", family = "Arial Narrow"),
#     axis.text.y =  element_text(size = 12, face = "bold", colour = "black", family = "Arial Narrow"),
#     panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(colour = "black", fill = NA),
#     panel.background = element_rect(fill = "transparent"),
#     strip.text.x = element_text(size = 14, face = "bold"),
#     strip.background = element_rect(colour = "black", fill = "transparent", linetype = "solid")
#   )
# 
# correlation_plot_D_RT <- data_total_RT %>%
#   filter(Hier == "D" & ROI == "D") %>%
#   ggplot(aes(x = accuracy, y = mean_rt)) +
#   xlab("Decoding accuracy") +
#   ylab("Mean reaction time (ms)") +
#   coord_cartesian(xlim=c(15,45)) +
#   scale_x_continuous(breaks = c(15, 30, 45)) +
#   geom_point(aes(colour = factor(Subject))) +
#   geom_smooth(method = 'lm', formula = y~x, se = TRUE) +
#   ggpubr::stat_cor(label.y.npc = 0.7, label.x.npc = 0.7) +
#   theme(
#     legend.position = "none",
#     plot.title = element_text(hjust = 0.5, size=20, face="bold.italic"),
#     axis.title.x = element_text(size=16, face="bold", family = "Arial Narrow"),
#     axis.title.y = element_text(size=16, face="bold", family = "Arial Narrow"),
#     axis.text.x =  element_text(size = 12, face = "bold", colour = "black", family = "Arial Narrow"),
#     axis.text.y =  element_text(size = 12, face = "bold", colour = "black", family = "Arial Narrow"),
#     panel.grid.major = element_line(colour = "grey", linetype = "dotted"),
#     panel.grid.minor = element_blank(),
#     panel.border = element_rect(colour = "black", fill = NA),
#     panel.background = element_rect(fill = "transparent"),
#     strip.text.x = element_text(size = 14, face = "bold"),
#     strip.background = element_rect(colour = "black", fill = "transparent", linetype = "solid")
#   )
