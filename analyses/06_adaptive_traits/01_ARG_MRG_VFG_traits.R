# =============================================================================
# Glacier-to-river -- adaptive traits / Adaptive traits of MAGs
# =============================================================================
# Combine ARG, MRG, VFG and POP-degradation gene counts with MAG metadata to
# relate mobile-element load to adaptive traits.
#
# Manuscript      : -
# Main outputs    : per-MAG adaptive-trait tables
# Original script : 20gaoyuan/20ARGMRGVFG/1_Adaptive_traits.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))
ARG <- read.delim2("rgiargs_2.txt",header = F)
colnames(ARG) <- c("ORF","drug_class","resistance_mechanism","AMR_gene_family",
                   "Antibiotic")

library(stringr)
library(dplyr)
ARG$genome <- str_extract(ARG$ORF, ".*bin\\.[0-9]+")

RGI_ARG <- table(ARG$genome) %>% as.data.frame()
colnames(RGI_ARG) <- c("genome","ARG")

MRG <- read.delim2("diamond_MRG_out3_type.csv")
colnames(MRG) <- c("genome","MRG")

VFG <- read.delim2("diamond_VFG_out3_type.csv")
colnames(VFG) <- c("genome","VFG")

POPdeg <- read.delim2("diamond_POPDeg_out3_type.csv")
colnames(POPdeg) <- c("genome","POPdeg")

MAGmetadata <- read.delim2("MAG5455_metadata.tsv")
MAGmetadata <- select(MAGmetadata, genome)

RGI_ARG <- RGI_ARG[RGI_ARG$genome %in% MAGmetadata$genome,]
MRG <- MRG[MRG$genome %in% MAGmetadata$genome,]
VFG <- VFG[VFG$genome %in% MAGmetadata$genome,]
POPdeg <- POPdeg[POPdeg$genome %in% MAGmetadata$genome,]

library(dplyr)
library(tidyr)
library(stringr)
DF_res <- read.delim2("merged_5467MAG_defense_finder_systems.tsv",check.names = F,
                      sep = '\t'
)
DF_defense <- DF_res[DF_res$activity == "Defense",]
DF_defense <- DF_defense %>%
  mutate(
    genome = str_extract(sys_id, ".*bin\\.\\d+")  # 包含"bin.数字"
  )
# 创建计数矩阵
count_matrix <- DF_defense %>%
  count(genome, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()
# 设置行名
rownames(count_matrix) <- count_matrix$genome
count_matrix <- count_matrix[-1]  # 移除genome列
genome_defense_num <- rowSums(count_matrix) %>% as.data.frame() %>% setNames("DS_num")
genome_defense_num$genome <- rownames(genome_defense_num)
genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAGmetadata$genome,]

MAGmetadata <- MAGmetadata %>%
  full_join( RGI_ARG,by = "genome") %>%
  full_join(MRG, by = "genome") %>%
  full_join(VFG, by = "genome") %>%
  full_join(POPdeg, by = "genome") %>%
  full_join(genome_defense_num, by = "genome")
MAGmetadata[is.na(MAGmetadata)] <- 0

#AMP
AMP_res <- read.delim2("macrel.out.prediction")
AMP_res$genome <- str_extract(AMP_res$Access, ".*bin\\.[0-9]+")
# 创建计数矩阵
AMP_count_matrix <- table(AMP_res$genome) %>% as.data.frame()
colnames(AMP_count_matrix) <- c("genome","AMP")

#BGC
BGC_res <- read.delim2("antismash_res_type.tsv",check.names = F,row.names = 1)
BGC_res_num <- data.frame(
  genome = rownames(BGC_res),BGC = rowSums(sapply(BGC_res, as.numeric))
  ) 

MAGmetadata <- MAGmetadata %>%
  full_join(AMP_count_matrix,by = "genome") %>%
  full_join(BGC_res_num, by = "genome")
MAGmetadata[is.na(MAGmetadata)] <- 0


#统计各种基因之间的关系
library(psych)

corr_test <- corr.test(MAGmetadata[,-1], method = "spearman")
corr_test$r  # 相关系数矩阵
corr_test$p  # p值矩阵
library(corrplot)
corrplot(corr_test$r, 
         method = "color",
         type = "upper",
         tl.col = "black",
         addCoef.col = "black",
         number.cex = 0.7)

# 进行线性回归分析
model <- lm(ARG ~ DS_num, data = MAGmetadata)
# 查看回归模型的摘要
summary(model)

library(ggplot2)
# 绘制散点图和回归线
p1 <- ggplot(MAGmetadata, aes(x = DS_num, y = ARG)) +
  geom_point() +  # 添加散点图
  geom_smooth(method = "nls", se = TRUE) +  # 添加线性回归线，不显示置信区间
  labs(title = NULL, x = "DS number (per MAG)", y = "ARGs number (per MAG)")+  # 添加标题和坐标轴标签
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p1

p2 <- ggplot(MAGmetadata, aes(x = ARG, y = MRG)) +
  geom_point() +  # 添加散点图
  geom_smooth(method = "nls", se = TRUE) +  # 添加线性回归线，不显示置信区间
  labs(title = NULL, x = "DS number (per MAG)", y = "MRG number (per MAG)")+  # 添加标题和坐标轴标签
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.8,0.7),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p2


df.count <- MAGmetadata %>%
  mutate(DS_num = factor(DS_num, levels = 0:83)) %>%
  mutate(group = case_when(DS_num == 1 ~ "group1", TRUE ~ "group2"))
df.msd <- MAGmetadata %>%
  group_by(DS_num) %>%
  summarise_at(vars(ARG), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>% 
  as.data.frame() %>%
  mutate(DS_num = as.numeric(DS_num))
df.msd
# Visualize
library(scales)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(cowplot)
set.seed(12345)
p1 <- ggplot() +
  geom_violin(data = df.count, aes(x = DS_num, y = ARG, fill = DS_num),
              position = position_dodge(0.5), linewidth = 0.5,drop = F) +
  geom_jitter(data = df.count, aes(x = DS_num, y = ARG), size = 0.01, color = "grey") +
  # scale_x_discrete(breaks = 0:25, labels = c(0:25, ">=25")) +
  scale_y_continuous(trans = log2_trans(), breaks = c(1,2,4,8,16,32,64), labels = c(1,2,4,8,16,32,"≥64")) +
  scale_fill_viridis_d(option = 'C')+
  # scale_fill_manual(values = c("1" = "#FFEFCC", "2" = "#FECCCC", "3" = "#99CCCC", "4" = "#CCCCFF", "5" = "#3399CC"))+
  labs(title = "", x = "DS count", y = "ARG count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        #  axis.line.x = element_line(linewidth = 0.5),
        #  axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(angle = 90,size = 10,hjust = 0.5,vjust = 0.5),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p1


p2 <- ggscatter(df.msd, x = "DS_num", y = "max", add = "reg.line", add.params = list(color = "blue", fill = NA, linetype = "dashed")) +
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), label.x = 2, label.y = 30) +
  stat_regline_equation(label.x = 2, label.y = 35) +
  labs(title = "", x = "DS count", y = "Max ARG count", fill = "") +
  theme_bw() +
  theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p2

#批量绘制
plot_logistic_fit <- function(MAGmetadata, geneA, geneB) {
  
  # 数据处理
  df.count <- MAGmetadata %>%
    mutate(!!geneA := factor(!!sym(geneA), levels = 0:83)) %>%
    mutate(group = case_when(!!sym(geneA) == 1 ~ "group1", TRUE ~ "group2"))
  
  df.msd <- MAGmetadata %>%
    group_by(!!sym(geneA)) %>%
    summarise_at(vars(!!geneB), list(max = max, mean = mean, sd = sd, median = median, skew = moments::skewness)) %>%
    as.data.frame() %>%
    mutate(!!geneA := as.numeric(!!sym(geneA)))
  
  # 图1: Violin 图
  pall <- ggplot() +
    geom_violin(data = df.count, aes(x = !!sym(geneA), y = !!sym(geneB), fill = !!sym(geneA)),
                position = position_dodge(0.5), linewidth = 0.5, drop = F) +
    geom_jitter(data = df.count, aes(x = !!sym(geneA), y = !!sym(geneB)), size = 0.01, color = "grey") +
    scale_y_continuous(trans = log2_trans(), breaks = c(1, 2, 4, 8, 16, 32, 64), labels = c(1, 2, 4, 8, 16, 32, "≥64")) +
    scale_fill_viridis_d(option = 'C') +
    labs(x = paste(geneA, "count"), y = paste(geneB, "count"), fill = "") +
    theme_bw() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = 0.5),
          panel.grid = element_blank(),
          # axis.line.x = element_line(linewidth = 0.5),
          # axis.line.y = element_line(linewidth = 0.5),
          axis.title.x = element_text(size = 12, face = "plain"),
          axis.title.y = element_text(size = 12, face = "plain"),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.y = element_line(linewidth = 0.5),
          axis.ticks.length.x = unit(0.1, "cm"),
          axis.ticks.length.y = unit(0.1, "cm"),
          axis.text.x = element_text(size = 10,angle = 90),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  
  # 图2: Logistic 拟合曲线
  pall2 <- ggplot(df.msd, aes(x = !!sym(geneA), y = mean)) +
    geom_point() +  # 绘制散点
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),method = "spearman") +
    geom_smooth(method = "lm", se = TRUE) +  # 添加线性回归线，不显示置信区间
    #    stat_regline_equation(label.x = 2, label.y = 35) +
    labs(x = paste(geneA, "count"), y = paste("Mean", geneB, "count")) +
    theme_bw() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = 0.5),
          panel.grid = element_blank(),
          # axis.line.x = element_line(linewidth = 0.5),
          # axis.line.y = element_line(linewidth = 0.5),
          axis.title.x = element_text(size = 12, face = "plain"),
          axis.title.y = element_text(size = 12, face = "plain"),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.y = element_line(linewidth = 0.5),
          axis.ticks.length.x = unit(0.1, "cm"),
          axis.ticks.length.y = unit(0.1, "cm"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  # write.csv(df.count,"figs12_1.csv")
  #  write.csv(df.msd,"figs12_2.csv")
  # 返回两个图
  return(list(violin_plot = pall, logistic_plot = pall2))
}


plot_logistic_fit <- function(MAGmetadata, geneA, geneB) {
  
  # 第一步：统计每个geneA分组的样本数并按geneA排序
  group_stats <- MAGmetadata %>%
    group_by(!!sym(geneA)) %>%
    summarise(n = n(), .groups = "drop") %>%
    arrange(!!sym(geneA))
  
  # 第二步：创建合并规则（只要分组样本数<30就一直合并）
  merged_groups <- data.frame(
    original = group_stats[[geneA]],
    merged = group_stats[[geneA]],  # 初始化为不合并
    stringsAsFactors = FALSE
  )
  
  i <- 1
  while(i <= nrow(group_stats)) {
    if(group_stats$n[i] < 10) {
      # 找到需要合并的起始分组
      merge_start <- i
      merge_value <- group_stats[[geneA]][i]
      
      # 持续合并直到最后一个分组
      while(i <= nrow(group_stats)) {
        merged_groups$merged[i] <- merge_value
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  
  # 第三步：应用合并规则到原始数据
  MAGmetadata_merged <- MAGmetadata %>%
    left_join(merged_groups, by = setNames("original", geneA)) %>%
    mutate(!!geneA := merged) %>%
    select(-merged)
  
  # 第四步：处理合并后的数据
  df.count <- MAGmetadata_merged %>%
    mutate(!!geneA := factor(!!sym(geneA))) %>%
    mutate(group = case_when(!!sym(geneA) == 1 ~ "group1", TRUE ~ "group2"))
  
  df.msd <- MAGmetadata_merged %>%
    group_by(!!sym(geneA)) %>%
    summarise(
      max = max(!!sym(geneB), na.rm = TRUE),
      mean = mean(!!sym(geneB), na.rm = TRUE),
      sd = sd(!!sym(geneB), na.rm = TRUE),
      median = median(!!sym(geneB), na.rm = TRUE),
      skew = moments::skewness(!!sym(geneB), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(!!geneA := as.numeric(as.character(!!sym(geneA))))
  # 图1: Violin 图
  pall <- ggplot() +
    geom_violin(data = df.count, aes(x = !!sym(geneA), y = !!sym(geneB), fill = !!sym(geneA)),
                position = position_dodge(0.5), linewidth = 0.5, drop = F) +
    geom_jitter(data = df.count, aes(x = !!sym(geneA), y = !!sym(geneB)), size = 0.01, color = "grey") +
    scale_y_continuous(trans = log2_trans(), breaks = c(1, 2, 4, 8, 16, 32, 64), labels = c(1, 2, 4, 8, 16, 32, "≥64")) +
    scale_fill_viridis_d(option = 'C') +
    labs(x = paste(geneA, "count"), y = paste(geneB, "count"), fill = "") +
    theme_bw() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = 0.5),
          panel.grid = element_blank(),
          # axis.line.x = element_line(linewidth = 0.5),
          # axis.line.y = element_line(linewidth = 0.5),
          axis.title.x = element_text(size = 12, face = "plain"),
          axis.title.y = element_text(size = 12, face = "plain"),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.y = element_line(linewidth = 0.5),
          axis.ticks.length.x = unit(0.1, "cm"),
          axis.ticks.length.y = unit(0.1, "cm"),
          axis.text.x = element_text(size = 10,angle = 90),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  
  # 图2: Logistic 拟合曲线
  pall2 <- ggplot(df.msd, aes(x = !!sym(geneA), y = mean)) +
    geom_point() +  # 绘制散点
    stat_cor(aes(label = paste(..r.label.., ..p.label.., sep = "~`,`~")),method = "spearman") +
    geom_smooth(method = "lm", se = TRUE) +  # 添加线性回归线，不显示置信区间
    #    stat_regline_equation(label.x = 2, label.y = 35) +
    labs(x = paste(geneA, "count"), y = paste("Mean", geneB, "count")) +
    theme_bw() +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"),
          plot.background = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(linewidth = 0.5),
          panel.grid = element_blank(),
          # axis.line.x = element_line(linewidth = 0.5),
          # axis.line.y = element_line(linewidth = 0.5),
          axis.title.x = element_text(size = 12, face = "plain"),
          axis.title.y = element_text(size = 12, face = "plain"),
          axis.ticks.x = element_line(linewidth = 0.5),
          axis.ticks.y = element_line(linewidth = 0.5),
          axis.ticks.length.x = unit(0.1, "cm"),
          axis.ticks.length.y = unit(0.1, "cm"),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          legend.position = "none")
  # write.csv(df.count,"figs12_1.csv")
  #  write.csv(df.msd,"figs12_2.csv")
  # 返回两个图
  return(list(violin_plot = pall, logistic_plot = pall2))
}








# 调用函数
geneA <- "DS_num"
geneB <- "ARG"
plot_DS_arg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_DS_arg$violin_plot
plot_DS_arg$logistic_plot

geneA <- "DS_num"
geneB <- "MRG"
plot_DS_MRG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_DS_MRG$violin_plot
plot_DS_MRG$logistic_plot

geneA <- "DS_num"
geneB <- "VFG"
plot_DS_VFG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_DS_VFG$violin_plot
plot_DS_VFG$logistic_plot

geneA <- "DS_num"
geneB <- "POPdeg"
plot_DS_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_DS_POPdeg$violin_plot
plot_DS_POPdeg$logistic_plot

geneA <- "ARG"
geneB <- "MRG"
plot_ARG_MRG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_ARG_MRG$violin_plot
plot_ARG_MRG$logistic_plot

geneA <- "ARG"
geneB <- "VFG"
plot_ARG_VFG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_ARG_VFG$violin_plot
plot_ARG_VFG$logistic_plot

geneA <- "ARG"
geneB <- "POPdeg"
plot_ARG_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_ARG_POPdeg$violin_plot
plot_ARG_POPdeg$logistic_plot

geneA <- "MRG"
geneB <- "VFG"
plot_MRG_VFG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_MRG_VFG$violin_plot
plot_MRG_VFG$logistic_plot

geneA <- "MRG"
geneB <- "POPdeg"
plot_MRG_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_MRG_POPdeg$violin_plot
plot_MRG_POPdeg$logistic_plot

geneA <- "VFG"
geneB <- "POPdeg"
plot_VFG_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_VFG_POPdeg$violin_plot
plot_VFG_POPdeg$logistic_plot


library(Cairo)
library(cowplot)
row1 <- plot_grid(plot_DS_arg$violin_plot, plot_DS_VFG$violin_plot,
                  # plot_DS_MRG$violin_plot, plot_DS_POPdeg$violin_plot,
                  # plot_ARG_VFG$violin_plot,
                  labels = letters[1:2], nrow = 1
                  #rel_heights = c(1, 1, 1.5),
                  # label_x = 0.09,            # 标签x坐标（1为最右侧）
                  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
row1

row2 <- plot_grid(#plot_DS_arg$violin_plot, plot_DS_VFG$violin_plot,
                  plot_DS_MRG$violin_plot, plot_DS_POPdeg$violin_plot,
                  #plot_ARG_VFG$violin_plot,
                  labels = letters[3:4], nrow = 1
                  #rel_heights = c(1, 1, 1.5),
                  # label_x = 0.09,            # 标签x坐标（1为最右侧）
                  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
row2

row3 <- plot_grid(plot_DS_arg$logistic_plot, plot_DS_VFG$logistic_plot,
  plot_DS_MRG$logistic_plot, plot_DS_POPdeg$logistic_plot,
  #plot_ARG_VFG$violin_plot,
  labels = letters[5:8], nrow = 1
  #rel_heights = c(1, 1, 1.5),
  # label_x = 0.09,            # 标签x坐标（1为最右侧）
  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
row3

final_plot1 <- plot_grid(row1, row2, row3, nrow = 3
                         # rel_heights = c(1, 1, 1.5)
)
final_plot1

setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))
cairo_pdf(width = 10,height = 8,filename = "1_DS_4genes.pdf")
final_plot1
dev.off()


row4 <- plot_grid(
  plot_ARG_VFG$violin_plot, plot_ARG_MRG$violin_plot, 
  plot_ARG_POPdeg$violin_plot,
  labels = letters[1:3], nrow = 1
  #rel_heights = c(1, 1, 1.5),
  # label_x = 0.09,            # 标签x坐标（1为最右侧）
  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
row4

final_plot2 <- plot_grid(row4, 
                         plot_MRG_VFG$violin_plot, 
                         plot_MRG_POPdeg$violin_plot, 
                         plot_VFG_POPdeg$violin_plot, 
                         nrow = 4
                        # rel_heights = c(1, 1, 1.5)
                         )
final_plot2

setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))
cairo_pdf(width = 8,height = 7,filename = "1_4genes_inter.pdf")
final_plot2
dev.off()

final_plot3 <- plot_grid(plot_ARG_VFG$logistic_plot, 
                         plot_ARG_MRG$logistic_plot, 
                         plot_ARG_POPdeg$logistic_plot,
                         plot_MRG_VFG$logistic_plot, 
                         plot_MRG_POPdeg$logistic_plot, 
                         plot_VFG_POPdeg$logistic_plot, 
                         nrow = 2,ncol = 3,
                         labels = letters[1:6]
                         # rel_heights = c(1, 1, 1.5)
)
final_plot3

setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))
cairo_pdf(width = 7.38,height = 4.46,filename = "1_4genes_inter_sperman_cor.pdf")
final_plot3
dev.off()


#geneA和geneB返过来看，规律一致
geneA <- "POPdeg"
geneB <- "MRG"
plot_POPdeg_MRG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_POPdeg_MRG$violin_plot
plot_POPdeg_MRG$logistic_plot

geneA <- "POPdeg"
geneB <- "ARG"
plot_POPdeg_ARG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_POPdeg_ARG$violin_plot
plot_POPdeg_ARG$logistic_plot

geneA <- "POPdeg"
geneB <- "VFG"
plot_POPdeg_VFG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_POPdeg_VFG$violin_plot
plot_POPdeg_VFG$logistic_plot

#AMP BGC的规律####
setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))

geneA <- "DS_num"
geneB <- "AMP"
plot_DS_AMP <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_DS_AMP$violin_plot
plot_DS_AMP$logistic_plot

geneA <- "AMP"
geneB <- "DS_num"
plot_AMP_DS <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_AMP_DS$violin_plot
plot_AMP_DS$logistic_plot

geneA <- "DS_num"
geneB <- "BGC"
plot_DS_BGC <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_DS_BGC$violin_plot
plot_DS_BGC$logistic_plot

geneA <- "BGC"
geneB <- "DS_num"
plot_BGC_DS <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_DS$violin_plot
plot_BGC_DS$logistic_plot

geneA <- "AMP"
geneB <- "ARG"
plot_AMP_ARG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_AMP_ARG$violin_plot
plot_AMP_ARG$logistic_plot

geneA <- "AMP"
geneB <- "MRG"
plot_AMP_MRG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_AMP_MRG$violin_plot
plot_AMP_MRG$logistic_plot

geneA <- "AMP"
geneB <- "VFG"
plot_AMP_VFG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_AMP_VFG$violin_plot
plot_AMP_VFG$logistic_plot

geneA <- "AMP"
geneB <- "POPdeg"
plot_AMP_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_AMP_POPdeg$violin_plot
plot_AMP_POPdeg$logistic_plot

geneA <- "AMP"
geneB <- "BGC"
plot_AMP_BGC <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_AMP_BGC$violin_plot
plot_AMP_BGC$logistic_plot

geneA <- "BGC"
geneB <- "ARG"
plot_BGC_ARG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_ARG$violin_plot
plot_BGC_ARG$logistic_plot

geneA <- "ARG"
geneB <- "BGC"
plot_ARG_BGC <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_ARG_BGC$violin_plot
plot_ARG_BGC$logistic_plot


geneA <- "BGC"
geneB <- "MRG"
plot_BGC_MRG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_MRG$violin_plot
plot_BGC_MRG$logistic_plot

geneA <- "BGC"
geneB <- "VFG"
plot_BGC_VFG <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_VFG$violin_plot
plot_BGC_VFG$logistic_plot

geneA <- "BGC"
geneB <- "POPdeg"
plot_BGC_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_POPdeg$violin_plot
plot_BGC_POPdeg$logistic_plot

geneA <- "BGC"
geneB <- "POPdeg"
plot_BGC_POPdeg <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_POPdeg$violin_plot
plot_BGC_POPdeg$logistic_plot

final_plot4 <- plot_grid(plot_DS_AMP$logistic_plot,
                         plot_DS_BGC$logistic_plot, 
                         plot_AMP_ARG$logistic_plot,
                         plot_AMP_MRG$logistic_plot, 
                         plot_AMP_VFG$logistic_plot, 
                         plot_AMP_POPdeg$logistic_plot, 
                         plot_AMP_BGC$logistic_plot,
                         plot_BGC_ARG$logistic_plot,
                         plot_BGC_MRG$logistic_plot,
                         plot_BGC_VFG$logistic_plot,
                         plot_BGC_POPdeg$logistic_plot,
                         nrow = 3,ncol = 4,
                         labels = letters[7:18]
                         # rel_heights = c(1, 1, 1.5)
)
final_plot4

setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))
cairo_pdf(width = 7.38,height = 5.5,filename = "1_11genes_inter_sperman_cor.pdf")
final_plot4
dev.off()


geneA <- "ARG"
geneB <- "AMP"
plot_ARG_AMP <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_ARG_AMP$violin_plot
plot_ARG_AMP$logistic_plot

geneA <- "MRG"
geneB <- "AMP"
plot_MRG_AMP <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_MRG_AMP$violin_plot
plot_MRG_AMP$logistic_plot

geneA <- "VFG"
geneB <- "AMP"
plot_VFG_AMP <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_VFG_AMP$violin_plot
plot_VFG_AMP$logistic_plot

geneA <- "POPdeg"
geneB <- "AMP"
plot_POPdeg_AMP <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_POPdeg_AMP$violin_plot
plot_POPdeg_AMP$logistic_plot

geneA <- "BGC"
geneB <- "AMP"
plot_BGC_AMP <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_BGC_AMP$violin_plot
plot_BGC_AMP$logistic_plot

geneA <- "ARG"
geneB <- "DS_num"
plot_ARG_DS <- plot_logistic_fit(MAGmetadata, geneA, geneB)
plot_ARG_DS$violin_plot
plot_ARG_DS$logistic_plot


################################################################################
library(Cairo)
library(cowplot)
g_DS_others <- plot_grid(plot_DS_arg$violin_plot, plot_DS_MRG$violin_plot,
                  plot_DS_VFG$violin_plot, plot_DS_AMP$violin_plot,
                  plot_DS_BGC$violin_plot, plot_DS_POPdeg$violin_plot,
                  # plot_DS_MRG$violin_plot, plot_DS_POPdeg$violin_plot,
                  # plot_ARG_VFG$violin_plot,
                  labels = letters[1:6], nrow = 3
                  #rel_heights = c(1, 1, 1.5),
                  # label_x = 0.09,            # 标签x坐标（1为最右侧）
                  #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g_DS_others

g_positive_others <- 
  plot_grid(plot_ARG_MRG$violin_plot, plot_ARG_VFG$violin_plot,plot_ARG_BGC$violin_plot, 
            plot_ARG_POPdeg$violin_plot, plot_ARG_AMP$violin_plot,plot_ARG_DS$violin_plot,
            plot_AMP_ARG$violin_plot,plot_AMP_MRG$violin_plot,plot_AMP_VFG$violin_plot,
            plot_AMP_POPdeg$violin_plot,plot_AMP_BGC$violin_plot,plot_AMP_DS$violin_plot,
                         # plot_DS_MRG$violin_plot, plot_DS_POPdeg$violin_plot,
                         # plot_ARG_VFG$violin_plot,
                         labels = letters[1:12], nrow = 4
                         #rel_heights = c(1, 1, 1.5),
                         # label_x = 0.09,            # 标签x坐标（1为最右侧）
                         #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g_positive_others

# row3 <- plot_grid(plot_DS_arg$logistic_plot, plot_DS_VFG$logistic_plot,
#                   plot_DS_MRG$logistic_plot, plot_DS_POPdeg$logistic_plot,
#                   #plot_ARG_VFG$violin_plot,
#                   labels = letters[5:8], nrow = 1
#                   #rel_heights = c(1, 1, 1.5),
#                   # label_x = 0.09,            # 标签x坐标（1为最右侧）
#                   #  label_y = 0.9             # 标签y坐标（1为最上方）
# )
# row3
# 
# final_plot1 <- plot_grid(row1, row2, row3, nrow = 3
#                          # rel_heights = c(1, 1, 1.5)
# )
# final_plot1

setwd(file.path(GTR_DATA, "20gaoyuan/20ARGMRGVFG"))

cairo_pdf(width = 8.5,height = 8,filename = "1_DS_others.pdf")
g_DS_others
dev.off()

cairo_pdf(width = 8.5,height = 8,filename = "1_positive_others.pdf")
g_positive_others
dev.off()



