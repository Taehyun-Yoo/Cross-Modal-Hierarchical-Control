library(readxl)
library("tidyverse")
library("ggpubr")
library("grid")
library("rstatix")
library("ggplot2")
library("MBESS")
library(ggprism)
library(svglite)

dir <- ("/data/Project1/MSHP/analysis/")
result_dir <- paste0(dir, 'results_2nd/SVD_analysis_yz/')

H_contrast <- read.table(paste0(result_dir, 'Tcom_total.statdata'), quote="\"", comment.char="")

names(H_contrast)[1] <-c("component1")
names(H_contrast)[2] <-c("component2")
names(H_contrast)[3] <-c("Hierarchy")
names(H_contrast)[4] <-c("Directionality")
names(H_contrast)[5] <-c("Subject")

#out <- anova_test(data = H_contrast, dv = component2, wid = Subject, within = Hierarchy)
#K <- aov(component2 ~ factor(Hierarchy), data = H_contrast)
#Anova_H <- data.frame(p_value = out[[1]]$p, significance = out[[1]]$`p<.05`)

R <- subset(H_contrast$component1, H_contrast$Hierarchy == 1)
F <- subset(H_contrast$component1, H_contrast$Hierarchy == 2)
D <- subset(H_contrast$component1, H_contrast$Hierarchy == 3)

my_data <- H_contrast

for (k in 1:length(my_data$Hierarchy)){
  if (my_data$Hierarchy[k] == 1) {
    my_data$Hierarchy[k] = "R"
  } else if (my_data$Hierarchy[k] == 2) {
    my_data$Hierarchy[k] = "F"
  } else{
    my_data$Hierarchy[k] = "D"
  }
}

my_data$Hierarchy <- factor(my_data$Hierarchy, levels = c("R", "F", "D"))

#for (k in 1:length(my_data_RD$Hierarchy)){
#  if (my_data_RD$Hierarchy[k] == 1) {
#    my_data_RD$Hierarchy[k] = "R"
#  } else{
#    my_data_RD$Hierarchy[k] = "D"
#  }
#}

svd_hier <- my_data %>%
  group_by(Hierarchy) %>%
  summarise(se = sd(component1)/sqrt(n()),
            component1 = mean(component1)) %>%
  ungroup()

svd_hier$Hierarchy <- factor(svd_hier$Hierarchy, levels = c("R", "F", "D"))

#One_way_ANOVA_SVD (each hierarchy)
SVD_one_way_ANOVA <- function(res) {
  one_way_ANOVA <- anova_test(data = res, dv = component1, wid = Subject, within = Hierarchy, effect.size = "pes")
  table <- get_anova_table(one_way_ANOVA)
  effect <- c(table[1,"Effect"])
  F_value <- c(paste0("F(", table[1,"DFn"], ",", table[1,"DFd"], ") = ", table[1,"F"]))
  p_value <- c(table[1,"p"])
  p_value_star <- c(symnum(table[1,"p"], cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "†", "ns")))
  pes <- c(table[1,"pes"])
  lims <- conf.limits.ncf(F.value = table[1, "F"], df.1 = table[1, "DFn"], df.2 = table[1, "DFd"], conf.level = .90)
  ci <- c(lims$Lower.Limit/(lims$Lower.Limit + table[1, "DFn"] + table[1, "DFd"] + 1), lims$Upper.Limit/(lims$Upper.Limit + table[1, "DFn"] + table[1, "DFd"] + 1))
  ci <- round(ci, 3)
  ci <- paste0("[", ci[1], ", ", ci[2], "]")
  one_way_ANOVA_SVD_table <- data.frame(effect, F_value, p_value, p_value_star, pes, ci)
  one_way_ANOVA_SVD_table
}

write_csv(SVD_one_way_ANOVA(my_data), paste0(result_dir, 'svd_one_way_ANOVA.csv'))

svd_posthoc <- function(res) {
  post_t_test <- t_test(component1 ~ Hierarchy, p.adjust.method = "holm", data = res, paired = TRUE, detailed = TRUE)
  .y. <- c("projected_y", "projected_y", "projected_y")
  group1 <- c(post_t_test$group1)
  group2 <- c(post_t_test$group2)
  comparisons <- c(paste0(post_t_test$group1[1], "-", post_t_test$group2[1]), paste0(post_t_test$group1[2], "-", post_t_test$group2[2]), paste0(post_t_test$group1[3], "-", post_t_test$group2[3]))
  t_value <- c(paste0("t(", post_t_test$df[1], ") = ", round(post_t_test$statistic[1],3)), paste0("t(", post_t_test$df[2], ") = ", round(post_t_test$statistic[2],3)), paste0("t(", post_t_test$df[3], ") = ", round(post_t_test$statistic[3],3)))
  adjusted_p_value <- c(post_t_test$p.adj[1], post_t_test$p.adj[2], post_t_test$p.adj[3])
  adjusted_p_value_star <- c(symnum(post_t_test$p.adj[1], cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "†", "ns")), symnum(post_t_test$p.adj[2], cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "†", "")), symnum(post_t_test$p.adj[3], cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("***", "**", "*", "†", "")))
  d <- c(round(abs(post_t_test$statistic[1]/sqrt(post_t_test$n1[1])),3), round(abs(post_t_test$statistic[2]/sqrt(post_t_test$n1[2])),3), round(abs(post_t_test$statistic[3]/sqrt(post_t_test$n1[3])),3))
  ci <- c(paste0("[", round(post_t_test$conf.low[1],3), ", ", round(post_t_test$conf.high[1],3), "]"), paste0("[", round(post_t_test$conf.low[2],3), ", ", round(post_t_test$conf.high[2],3), "]"), paste0("[", round(post_t_test$conf.low[3],3), ", ", round(post_t_test$conf.high[3],3), "]"))
  y.position <- c(13, 19, 25)
  post_svd_table <- data.frame(.y., group1, group2, comparisons, t_value, adjusted_p_value, adjusted_p_value_star, d, ci, y.position)
  post_svd_table
}

write_csv(svd_posthoc(my_data), paste0(result_dir, 'svd_one_way_posthoc.csv'))

plot_svd <- function(res) {
  
  svd_hier <- res %>%
    group_by(Hierarchy) %>%
    summarise(se = sd(component1)/sqrt(n()),
              component1 = mean(component1)) %>%
    ungroup()
  
  my_plot_svd <- ggplot(svd_hier,aes(x=Hierarchy,y=component1)) +
    geom_point(data= res, stat = "identity", alpha = .1, size = .2) +
    geom_bar(aes(fill = Hierarchy), stat="identity", position="dodge", width = 0.6)+
    geom_errorbar(aes(ymin=component1-se,ymax=component1+se),
                  position=position_dodge(0.9), width = 0.1, alpha= 1)+
    scale_fill_manual(values = c("#2B6A6C", "#F29724", "#B80D48"))+
    scale_x_discrete(limits=c("R", "F", "D"))+
    coord_cartesian(ylim=c(-15, 20)) +
    scale_y_continuous(breaks = seq(-10, 10, 10))+
    #labs(x = "", y = "Value") + # y = "Projected Y-coordinate")+
    labs(y = "Projected Value") +
    theme(panel.grid.major=element_line(colour="Grey", linetype="dotted"), panel.grid.minor=element_blank(), panel.background=element_rect(fill="transparent",colour=NA), 
          axis.line = element_line(colour="Black"), 
          plot.title=element_text(hjust = 0.5, face="bold", size=20, vjust=2),
          #axis.title.x = element_text(size=20),
          axis.title.x = element_blank(),
          axis.text = element_text(size=18, colour = "Black"),
          axis.text.y = element_text(size=15, colour = "Black"),
          axis.title.y = element_text(size=20),
          strip.text = element_text(size = 9, colour = "white", face = "bold"),
          legend.position = "none")
  
  total_K_svd <- function(res){
    K <- svd_posthoc(res) 
    K <- tibble(K)
    K <- subset(K, K$adjusted_p_value_star != "" & K$adjusted_p_value_star != "ns")
    K
  }
  
  if (nrow(total_K_svd(res)) != 0){
    my_plot_svd <- my_plot_svd +
      add_pvalue(total_K_svd(res),
                 xmin = "group1", xmax = "group2",
                 label = "adjusted_p_value_star",
                 y.position = "y.position",
                 bracket.size = 0.8,
                 label.size = 7)
  } else {
    my_plot_svd <- my_plot_svd
  }
  
  my_plot_svd
  
  ggsave(paste0(result_dir,"SVD_final.svg"), width = 120, height = 120, units = "mm")
}

plot_svd(my_data)

# example <- read.table("~/Desktop/MSHP/SVD/total_data_left_0_01_k27.txt", quote="\"", comment.char="")
# names(example)[1] <-c("x_co")
# names(example)[2] <-c("y_co")
# names(example)[3] <-c("z_co")
# names(example)[4] <-c("Hierarchy")
# names(example)[5] <-c("Directionality")
# names(example)[6] <-c("Subject")
# 
# k <- c()
# for (i in 1:dim(example)[1]){
#   k[i] = sum(example[i,6] == example[,6]);
# }
# 
# index <- c()
# for (j in 1:length(k)){
#   if (k[j] == 3) {
#     index[j] = 1;
#   } else {
#     index[j] = 0;
#   }
# }
# 
# example <- cbind(example, index);
# example <- subset(example, example$index == 1);
# 
# ggplot(example, aes(y_co, z_co)) +
#   scale_x_reverse()+
#   geom_point(shape = 21, mapping = aes(x = y_co, y = z_co), colour = c("#2B6A6C", "#F29724", "#B80D48")[example$Hierarchy], fill = "white", size = 2, stroke = 2)+
#   coord_cartesian(xlim=c(80, -30), ylim=c(-30, 80)) +
#   xlab("Rostral-Caudal") +
#   ylab("Ventral-Dorsal") +
#   theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.background=element_rect(fill="transparent",colour=NA), 
#         axis.line = element_line(colour="Black"), 
#         plot.title=element_text(hjust = 0.5, face="bold", size=20, vjust=2),
#         #axis.title.x = element_text(size=20),
#         axis.title.x = element_text(size=20),
#         axis.text = element_blank(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_text(size=20),
#         strip.text = element_text(size = 9, colour = "white", face = "bold"),
#         legend.position = "none")
# figure_dir <- ("/home/taehyun/Desktop/MSHP/figure_code/SVD/")
# ggsave(paste0(figure_dir,"0928_original_coordinate.svg"), width = 120, height = 120, units = "mm")
# 
# colour <- c("#2B6A6C", "#F29724", "#B80D48")
# colour <- colour[as.numeric(example$Hierarchy)]
# 
# scatter3D(example[,1], example[,2], example[,3], phi = 1.5, cex = 1.2)
# 
# example_rotated <- example
# example_rotated[,1] <- example[,2]
# example_rotated[,2] <- example[,1]
# 
# scatterplot3d(example[,1:3], angle = 145, tick.marks = TRUE, xlab = "Lateral-Medial", ylab = "Rostral-Caudal", zlab = "Ventral-Dorsal", pch = 1, color = colour, grid = FALSE, box = TRUE)
# 
# scatterplot3d(example[,1:3], angle = 230, tick.marks = FALSE, xlab = "Lateral-Medial", ylab = "Rostral-Caudal", zlab = "Ventral-Dorsal", pch = 1, color = colour, grid = FALSE, box = TRUE)
# ggsave(paste0(figure_dir,"0929_3D.svg"), width = 120, height = 120, units = "mm")
# 
# scatterplot3d(example_rotated[,1:3], angle = 210, tick.marks = FALSE, xlab = "Caudal-Rostral", ylab = "Medial-Lateral", zlab = "Ventral-Dorsal", pch = 1, color = colour, grid = FALSE, box = FALSE)
# ggsave(paste0(figure_dir,"0929_3D_rotated.svg"), width = 120, height = 120, units = "mm")

  
#### lme4 and emmeans ----------------------------------------------------
#library('lme4')
#library('emmeans')

# test.lm <- lmer(component2 ~ Hierarchy + (1 | Subject), data = H_contrast)

# repeated measures using aov
#K <- aov(component2 ~ factor(Hierarchy) + Error(factor(Subject)), data = H_contrast)
#summary(K)
# post-hoc test using emmeans with Tukey's correction for MCP
#pairs(emmeans(K, "Hierarchy", adjust = "tukey"))