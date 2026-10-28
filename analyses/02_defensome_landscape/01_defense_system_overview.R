# =============================================================================
# Glacier-to-river -- defensome landscape / Defensome overview
# =============================================================================
# Build the defense-system count matrix from the DefenseFinder output and
# derive prevalence, per-MAG DS number/diversity and treemap abundance
# panels.
#
# Manuscript      : Fig. 2, Fig. S2, Fig. S3
# Main outputs    : defense-system matrices and overview figures
# Original script : 20gaoyuan/10Defensefinder/1defensefinder_res.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
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

MAG_metadata <- read.csv("../01MAG_info/5455MAG_metadata.csv")

genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
count_matrix <- count_matrix[rownames(count_matrix) %in% MAG_metadata$genome,]
DF_defense <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]

MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]

DF_DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")
DF_DG <- DF_DG %>%
  mutate(
    genome = str_extract(replicon, ".*bin\\.\\d+")  # 包含"bin.数字"
  )
DF_DG <- DF_DG[DF_DG$activity == "Defense",]
DF_DG <- DF_DG[DF_DG$genome %in% MAG_metadata$genome, ]
#36622

#####1. 各个ecosystem的DS组成#########
DF_defense_2 <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]
DF_defense_2 <- merge(DF_defense,MAG_metadata,by.x = "genome",by.y = "genome")

# 统计每个ecosystem中各个type的出现次数
result <- DF_defense_2 %>%
  count(ecosystem, type) %>%  # 计算组合频数
  complete(ecosystem, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  )

#CJ
CJ_DS_num <- result[result$ecosystem == "CJ",]
library(reshape2)
CJ_DS_num <- melt(CJ_DS_num)
colnames(CJ_DS_num) <- c("ecosystem","DS","num")
CJ_defense <- DF_defense_2[DF_defense_2$ecosystem == "CJ",]
CJ_genome_DS <- count_matrix[rownames(count_matrix) %in% CJ_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(CJ_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(CJ_genome_DS)
})
# 转换为数据框并整理格式
CJ_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
CJ_DS_proportion$Proportion <- CJ_DS_proportion$Proportion*100

#LCJ
LCJ_DS_num <- result[result$ecosystem == "LCJ",]
#library(reshape2)
LCJ_DS_num <- melt(LCJ_DS_num)
colnames(LCJ_DS_num) <- c("ecosystem","DS","num")
LCJ_defense <- DF_defense_2[DF_defense_2$ecosystem == "LCJ",]
LCJ_genome_DS <- count_matrix[rownames(count_matrix) %in% LCJ_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(LCJ_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(LCJ_genome_DS)
})
# 转换为数据框并整理格式
LCJ_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
LCJ_DS_proportion$Proportion <- LCJ_DS_proportion$Proportion*100

#HH
HH_DS_num <- result[result$ecosystem == "HH",]
#library(reshape2)
HH_DS_num <- melt(HH_DS_num)
colnames(HH_DS_num) <- c("ecosystem","DS","num")
HH_defense <- DF_defense_2[DF_defense_2$ecosystem == "HH",]
HH_genome_DS <- count_matrix[rownames(count_matrix) %in% HH_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(HH_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(HH_genome_DS)
})
# 转换为数据框并整理格式
HH_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
HH_DS_proportion$Proportion <- HH_DS_proportion$Proportion*100

#YJ
YJ_DS_num <- result[result$ecosystem == "YJ",]
#library(reshape2)
YJ_DS_num <- melt(YJ_DS_num)
colnames(YJ_DS_num) <- c("ecosystem","DS","num")
YJ_defense <- DF_defense_2[DF_defense_2$ecosystem == "YJ",]
YJ_genome_DS <- count_matrix[rownames(count_matrix) %in% YJ_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(YJ_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(YJ_genome_DS)
})
# 转换为数据框并整理格式
YJ_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
YJ_DS_proportion$Proportion <- YJ_DS_proportion$Proportion*100

#Glacier
Glacier_DS_num <- result[result$ecosystem == "Glacier",]
#library(reshape2)
Glacier_DS_num <- melt(Glacier_DS_num)
colnames(Glacier_DS_num) <- c("ecosystem","DS","num")
Glacier_defense <- DF_defense_2[DF_defense_2$ecosystem == "Glacier",]
Glacier_genome_DS <- count_matrix[rownames(count_matrix) %in% Glacier_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Glacier_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Glacier_genome_DS)
})
# 转换为数据框并整理格式
Glacier_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Glacier_DS_proportion$Proportion <- Glacier_DS_proportion$Proportion*100





library(ggplot2)
# CJ
CJ_proportion_plot <- ggplot(CJ_DS_proportion[CJ_DS_proportion$id<=25,], aes(x = reorder(Defense_System, -Proportion), 
                                         y = Proportion)) +
  geom_bar(stat = "identity", fill = "#BEBCDF", width = 0.7) +
  labs(title = "",
       x = "",
       y = "Proportion (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )

print(CJ_proportion_plot)

# LCJ
LCJ_proportion_plot <- ggplot(LCJ_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                   y = Proportion)) +
  geom_bar(stat = "identity", fill = "#FFE2CE", width = 0.7) +
  labs(title = "",
       x = "",
       y = "Proportion (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )

print(LCJ_proportion_plot)

# HH
HH_proportion_plot <- ggplot(HH_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                      y = Proportion)) +
  geom_bar(stat = "identity", fill = "#FB8C62", width = 0.7) +
  labs(title = "",
       x = "",
       y = "Proportion (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )

print(HH_proportion_plot)

# YJ
YJ_proportion_plot <- ggplot(YJ_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                   y = Proportion)) +
  geom_bar(stat = "identity", fill = "#FCB2AF", width = 0.7) +
  labs(title = "",
       x = "",
       y = "Proportion (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )

print(YJ_proportion_plot)

# Glacier
Glacier_proportion_plot <- ggplot(Glacier_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                   y = Proportion)) +
  geom_bar(stat = "identity", fill = "#9BDFDF", width = 0.7) +
  labs(title = "",
       x = "",
       y = "Proportion (%)") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )

print(Glacier_proportion_plot)


library(Cairo)
library(cowplot)
cairo_pdf(filename = "1DS_proportion.pdf" , width = 9,height = 4.4)
g <- plot_grid(Glacier_proportion_plot,HH_proportion_plot,
               LCJ_proportion_plot,CJ_proportion_plot,
               YJ_proportion_plot,
               labels = letters[1:5], ncol = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

###2.MAG数量随DS数量的变化##############
#CJ
CJ_genome_DSnum <- rowSums(CJ_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
CJ_genome_DSnum$genome <- rownames(CJ_genome_DSnum)
CJ_MAG_metadata <- MAG_metadata[MAG_metadata$ecosystem == "CJ",]
CJ_MAG_metadata <- left_join(CJ_MAG_metadata, 
                             CJ_genome_DSnum, 
                             by = c("genome" = "genome"))
CJ_MAG_metadata$DS_num[is.na(CJ_MAG_metadata$DS_num)] <- 0
CJ_MAG_DSnum <- select(CJ_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
CJ_count_data <- CJ_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#HH
HH_genome_DSnum <- rowSums(HH_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
HH_genome_DSnum$genome <- rownames(HH_genome_DSnum)
HH_MAG_metadata <- MAG_metadata[MAG_metadata$ecosystem == "HH",]
HH_MAG_metadata <- left_join(HH_MAG_metadata, 
                             HH_genome_DSnum, 
                             by = c("genome" = "genome"))
HH_MAG_metadata$DS_num[is.na(HH_MAG_metadata$DS_num)] <- 0
HH_MAG_DSnum <- select(HH_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
HH_count_data <- HH_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#LCJ
LCJ_genome_DSnum <- rowSums(LCJ_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
LCJ_genome_DSnum$genome <- rownames(LCJ_genome_DSnum)
LCJ_MAG_metadata <- MAG_metadata[MAG_metadata$ecosystem == "LCJ",]
LCJ_MAG_metadata <- left_join(LCJ_MAG_metadata, 
                             LCJ_genome_DSnum, 
                             by = c("genome" = "genome"))
LCJ_MAG_metadata$DS_num[is.na(LCJ_MAG_metadata$DS_num)] <- 0
LCJ_MAG_DSnum <- select(LCJ_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
LCJ_count_data <- LCJ_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#YJ
YJ_genome_DSnum <- rowSums(YJ_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
YJ_genome_DSnum$genome <- rownames(YJ_genome_DSnum)
YJ_MAG_metadata <- MAG_metadata[MAG_metadata$ecosystem == "YJ",]
YJ_MAG_metadata <- left_join(YJ_MAG_metadata, 
                             YJ_genome_DSnum, 
                             by = c("genome" = "genome"))
YJ_MAG_metadata$DS_num[is.na(YJ_MAG_metadata$DS_num)] <- 0
YJ_MAG_DSnum <- select(YJ_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
YJ_count_data <- YJ_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Glacier
Glacier_genome_DSnum <- rowSums(Glacier_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Glacier_genome_DSnum$genome <- rownames(Glacier_genome_DSnum)
Glacier_MAG_metadata <- MAG_metadata[MAG_metadata$ecosystem == "Glacier",]
Glacier_MAG_metadata <- left_join(Glacier_MAG_metadata, 
                             Glacier_genome_DSnum, 
                             by = c("genome" = "genome"))
Glacier_MAG_metadata$DS_num[is.na(Glacier_MAG_metadata$DS_num)] <- 0
Glacier_MAG_DSnum <- select(Glacier_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Glacier_count_data <- Glacier_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

library(ggplot2)
#CJ
CJ_magnum_dsnum <- ggplot(CJ_count_data, aes(x = factor(DS_num), y = genome_count, 
                                       fill = "#BEBCDF")) +
  geom_col(width = 0.7, fill = "#BEBCDF", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
      #subtitle = paste("Total", sum(CJ_count_data$genome_count), "genomes"),
       x = "",
       y = "") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
CJ_magnum_dsnum

#HH
HH_magnum_dsnum <- ggplot(HH_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#BEBCDF")) +
  geom_col(width = 0.7, fill = "#FB8C62", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
       #subtitle = paste("Total", sum(HH_count_data$genome_count), "genomes"),
       x = "",
       y = "") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
HH_magnum_dsnum

#LCJ
LCJ_magnum_dsnum <- ggplot(LCJ_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#BEBCDF")) +
  geom_col(width = 0.7, fill = "#FFE2CE", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
       #subtitle = paste("Total", sum(LCJ_count_data$genome_count), "genomes"),
       x = "",
       y = "") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
LCJ_magnum_dsnum

#YJ
YJ_magnum_dsnum <- ggplot(YJ_count_data, aes(x = factor(DS_num), y = genome_count, 
                                               fill = "#BEBCDF")) +
  geom_col(width = 0.7, fill = "#FCB2AF", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
       #subtitle = paste("Total", sum(YJ_count_data$genome_count), "genomes"),
       x = "",
       y = "") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
YJ_magnum_dsnum

#Glacier
Glacier_magnum_dsnum <- ggplot(Glacier_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#BEBCDF")) +
  geom_col(width = 0.7, fill = "#9BDFDF", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
       #subtitle = paste("Total", sum(Glacier_count_data$genome_count), "genomes"),
       x = "",
       y = "") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    panel.grid.major.x = element_blank(),
    legend.position = "none"
  )+
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 12, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 12, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 90, size = 8, hjust = 1, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
Glacier_magnum_dsnum

library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_MAGnum_DSnum.pdf" , width = 4,height = 7)
g <- plot_grid(Glacier_magnum_dsnum,HH_magnum_dsnum,
               CJ_magnum_dsnum,LCJ_magnum_dsnum,
               YJ_magnum_dsnum,
               labels = letters[7:11], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###3.数量与MAG大小关系###############

#CJ
CJ_MAG_quality <- MAG_quality[MAG_quality$Name %in% CJ_MAG_DSnum$genome,]

CJ_MAG_DSdensity <- merge(CJ_MAG_DSnum, CJ_MAG_quality, by.x = "genome",by.y = "Name")
CJ_MAG_DSdensity$MAG_size <- CJ_MAG_DSdensity$Genome_Size/1000000
CJ_MAG_DSdensity$DS_density <- CJ_MAG_DSdensity$DS_num/CJ_MAG_DSdensity$Genome_Size*1000000

CJ_MAG_DSdensity$magsize <- floor(CJ_MAG_DSdensity$MAG_size)

library(dplyr)
CJ_df_summary <- CJ_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
CJ_df_summary2 <- CJ_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(CJ_df_summary$magsize, CJ_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(CJ_df_summary2$magsize, CJ_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

CJ_df_summary$ecosystem <- 'CJ'
CJ_df_summary2$ecosystem <- 'CJ'

#HH
HH_MAG_quality <- MAG_quality[MAG_quality$Name %in% HH_MAG_DSnum$genome,]

HH_MAG_DSdensity <- merge(HH_MAG_DSnum, HH_MAG_quality, by.x = "genome",by.y = "Name")
HH_MAG_DSdensity$MAG_size <- HH_MAG_DSdensity$Genome_Size/1000000
HH_MAG_DSdensity$DS_density <- HH_MAG_DSdensity$DS_num/HH_MAG_DSdensity$Genome_Size*1000000

HH_MAG_DSdensity$magsize <- floor(HH_MAG_DSdensity$MAG_size)

library(dplyr)
HH_df_summary <- HH_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
HH_df_summary2 <- HH_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(HH_df_summary$magsize, HH_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(HH_df_summary2$magsize, HH_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

HH_df_summary$ecosystem <- 'HH'
HH_df_summary2$ecosystem <- 'HH'

#LCJ
LCJ_MAG_quality <- MAG_quality[MAG_quality$Name %in% LCJ_MAG_DSnum$genome,]

LCJ_MAG_DSdensity <- merge(LCJ_MAG_DSnum, LCJ_MAG_quality, by.x = "genome",by.y = "Name")
LCJ_MAG_DSdensity$MAG_size <- LCJ_MAG_DSdensity$Genome_Size/1000000
LCJ_MAG_DSdensity$DS_density <- LCJ_MAG_DSdensity$DS_num/LCJ_MAG_DSdensity$Genome_Size*1000000

LCJ_MAG_DSdensity$magsize <- floor(LCJ_MAG_DSdensity$MAG_size)

library(dplyr)
LCJ_df_summary <- LCJ_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
LCJ_df_summary2 <- LCJ_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(LCJ_df_summary$magsize, LCJ_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(LCJ_df_summary2$magsize, LCJ_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

LCJ_df_summary$ecosystem <- 'LCJ'
LCJ_df_summary2$ecosystem <- 'LCJ'

#YJ
YJ_MAG_quality <- MAG_quality[MAG_quality$Name %in% YJ_MAG_DSnum$genome,]

YJ_MAG_DSdensity <- merge(YJ_MAG_DSnum, YJ_MAG_quality, by.x = "genome",by.y = "Name")
YJ_MAG_DSdensity$MAG_size <- YJ_MAG_DSdensity$Genome_Size/1000000
YJ_MAG_DSdensity$DS_density <- YJ_MAG_DSdensity$DS_num/YJ_MAG_DSdensity$Genome_Size*1000000

YJ_MAG_DSdensity$magsize <- floor(YJ_MAG_DSdensity$MAG_size)

library(dplyr)
YJ_df_summary <- YJ_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
YJ_df_summary2 <- YJ_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(YJ_df_summary$magsize, YJ_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(YJ_df_summary2$magsize, YJ_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

YJ_df_summary$ecosystem <- 'YJ'
YJ_df_summary2$ecosystem <- 'YJ'

#Glacier
Glacier_MAG_quality <- MAG_quality[MAG_quality$Name %in% Glacier_MAG_DSnum$genome,]

Glacier_MAG_DSdensity <- merge(Glacier_MAG_DSnum, Glacier_MAG_quality, by.x = "genome",by.y = "Name")
Glacier_MAG_DSdensity$MAG_size <- Glacier_MAG_DSdensity$Genome_Size/1000000
Glacier_MAG_DSdensity$DS_density <- Glacier_MAG_DSdensity$DS_num/Glacier_MAG_DSdensity$Genome_Size*1000000

Glacier_MAG_DSdensity$magsize <- floor(Glacier_MAG_DSdensity$MAG_size)

library(dplyr)
Glacier_df_summary <- Glacier_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Glacier_df_summary2 <- Glacier_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Glacier_df_summary$magsize, Glacier_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Glacier_df_summary2$magsize, Glacier_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Glacier_df_summary$ecosystem <- 'Glacier'
Glacier_df_summary2$ecosystem <- 'Glacier'


DSnum_allecosystem <- rbind(Glacier_df_summary2,CJ_df_summary2,
                            HH_df_summary2,LCJ_df_summary2,
                            YJ_df_summary2)
DSdensity_allecosystem <- rbind(Glacier_df_summary,CJ_df_summary,
                            HH_df_summary,LCJ_df_summary,
                            YJ_df_summary)



library(ggplot2)
DSnum_allecosystem_plot <- 
ggplot(DSnum_allecosystem, aes(x = magsize, y = mean_value,color = ecosystem)) +
  geom_point(aes(color=ecosystem))+
  geom_line(aes(color=ecosystem))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = ecosystem),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  labs(x='MAG size (Mb)',y='DS number (per MAG)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "right",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
DSnum_allecosystem_plot

DSdensity_allecosystem_plot <- 
  ggplot(DSdensity_allecosystem, aes(x = magsize, y = mean_value,color = ecosystem)) +
  geom_point(aes(color=ecosystem))+
  geom_line(aes(color=ecosystem))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = ecosystem),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  labs(x='MAG size (Mb)',y='DS density (per Mb)')+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "right",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
DSdensity_allecosystem_plot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_DS_num-density_allecosystem.pdf" , width = 4.14,height = 5.25)
g <- plot_grid(DSnum_allecosystem_plot,DSdensity_allecosystem_plot,
               labels = letters[8], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###4.比较不同环境的Defese family多少，DS per MAG###############

CJ_DF_num <- sum(CJ_DS_num$num>0)

HH_DF_num <- sum(HH_DS_num$num>0)

LCJ_DF_num <- sum(LCJ_DS_num$num>0)

YJ_DF_num <- sum(YJ_DS_num$num>0)

Glacier_DF_num <- sum(Glacier_DS_num$num>0)

allecosystem_DF_num <- data.frame(ecosystem = c('Glacier','HH','CJ','LCJ','YJ'),
                                  DF_num = c(Glacier_DF_num, HH_DF_num, CJ_DF_num, LCJ_DF_num, YJ_DF_num),
                                  DS_mean = c(sum(Glacier_MAG_metadata$DS_num)/nrow(Glacier_MAG_metadata),
                                              sum(HH_MAG_metadata$DS_num)/nrow(HH_MAG_metadata),
                                              sum(CJ_MAG_metadata$DS_num)/nrow(CJ_MAG_metadata),
                                              sum(LCJ_MAG_metadata$DS_num)/nrow(LCJ_MAG_metadata),
                                              sum(YJ_MAG_metadata$DS_num)/nrow(YJ_MAG_metadata))
                                  )

mytheme <-   theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

# 绘制横条形图
p_compare_DStype <- ggplot(allecosystem_DF_num, aes(x = reorder(ecosystem,DF_num) ,y = DF_num)) +
  geom_bar(stat = "identity", aes(fill = ecosystem),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DF_num,2)),vjust = 0.5,hjust = 1.1,size = 4)+
  labs(x = "", y = "Number of DS type", title = "")+  # 设置标签和标题
  mytheme
p_compare_DStype

p_compare_DSmean <- ggplot(allecosystem_DF_num, aes(x = reorder(ecosystem,DS_mean) ,y = DS_mean)) +
  geom_bar(stat = "identity", aes(fill = ecosystem),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DS_mean,2)),vjust = 0.5,hjust = 1.1,size = 4)+
  labs(x = "", y = "Number of DS type", title = "")+  # 设置标签和标题
  mytheme
p_compare_DSmean

library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_DStype_DSmean_allecosystem.pdf" , width = 2,height = 5.25)
g <- plot_grid(p_compare_DStype,p_compare_DSmean,
               labels = letters[9], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###5.DS作用机制统计############
library(treemap)
library(RColorBrewer)
library(treemapify)
library(ggplot2)

DF_mechanism <- read.delim2('df_List_Systems_lpw.csv',sep = ',')
DF_mechanism$Effector[DF_mechanism$Effector == ""] <- "Unknown"
DF_mechanism$Mechnism[DF_mechanism$Mechnism == ""] <- "Unknown"

CJ_DS_average <- CJ_DS_num
CJ_DS_average$DS_average <- (CJ_DS_average$num) / nrow(CJ_MAG_metadata) 
CJ_DS_average <- left_join(
  CJ_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
CJ_DS_average[is.na(CJ_DS_average)] <- "Unknown"

HH_DS_average <- HH_DS_num
HH_DS_average$DS_average <- (HH_DS_average$num) / nrow(HH_MAG_metadata) 
HH_DS_average <- left_join(
  HH_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
HH_DS_average[is.na(HH_DS_average)] <- "Unknown"

LCJ_DS_average <- LCJ_DS_num
LCJ_DS_average$DS_average <- (LCJ_DS_average$num) / nrow(LCJ_MAG_metadata) 
LCJ_DS_average <- left_join(
  LCJ_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
LCJ_DS_average[is.na(LCJ_DS_average)] <- "Unknown"

YJ_DS_average <- YJ_DS_num
YJ_DS_average$DS_average <- (YJ_DS_average$num) / nrow(YJ_MAG_metadata) 
YJ_DS_average <- left_join(
  YJ_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
YJ_DS_average[is.na(YJ_DS_average)] <- "Unknown"

Glacier_DS_average <- Glacier_DS_num
Glacier_DS_average$DS_average <- (Glacier_DS_average$num) / nrow(Glacier_MAG_metadata) 
Glacier_DS_average <- left_join(
  Glacier_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Glacier_DS_average[is.na(Glacier_DS_average)] <- "Unknown"


gradient_colors <- c("#FFA5B8", "#8AC7F5", "#FFE4A3", "#9ADFDF","#C4D8E9", "#C0C0C0")

CJ_treeplot <- ggplot(CJ_DS_average, 
       aes(area = DS_average, fill = Mechnism, label = DS, subgroup = Mechnism)) +
  # 基础树状图
  geom_treemap(colour = "white", size = 1.2) +  # 白色边框
  # 子组边框（取消注释启用）
  # geom_treemap_subgroup_border(colour = "white", size = 2) +
  # 文本标签
  geom_treemap_text(
    aes(label = paste(DS, "\n",round(DS_average, 2))),
    colour = "black",           # 白色文字
    place = "centre",           # 居中显示
   # grow = TRUE,                # 允许文字扩展
    reflow = TRUE,              # 自动换行
   # family = "sans",            # 字体
    fontface = "plain",          # 加粗
    min.size = 4                # 最小字体大小
  ) +
  # 使用渐变协调色系
  scale_fill_manual(values = gradient_colors) +
  # 主题美化
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )
CJ_treeplot

HH_treeplot <- ggplot(HH_DS_average, 
                      aes(area = DS_average, fill = Mechnism, label = DS, subgroup = Mechnism)) +
  # 基础树状图
  geom_treemap(colour = "white", size = 1.2) +  # 白色边框
  # 子组边框（取消注释启用）
  # geom_treemap_subgroup_border(colour = "white", size = 2) +
  # 文本标签
  geom_treemap_text(
    aes(label = paste(DS, "\n",round(DS_average, 2))),
    colour = "black",           # 白色文字
    place = "centre",           # 居中显示
    # grow = TRUE,                # 允许文字扩展
    reflow = TRUE,              # 自动换行
    # family = "sans",            # 字体
    fontface = "plain",          # 加粗
    min.size = 4                # 最小字体大小
  ) +
  # 使用渐变协调色系
  scale_fill_manual(values = gradient_colors) +
  # 主题美化
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )
HH_treeplot

LCJ_treeplot <- ggplot(LCJ_DS_average, 
                      aes(area = DS_average, fill = Mechnism, label = DS, subgroup = Mechnism)) +
  # 基础树状图
  geom_treemap(colour = "white", size = 1.2) +  # 白色边框
  # 子组边框（取消注释启用）
  # geom_treemap_subgroup_border(colour = "white", size = 2) +
  # 文本标签
  geom_treemap_text(
    aes(label = paste(DS, "\n",round(DS_average, 2))),
    colour = "black",           # 白色文字
    place = "centre",           # 居中显示
    # grow = TRUE,                # 允许文字扩展
    reflow = TRUE,              # 自动换行
    # family = "sans",            # 字体
    fontface = "plain",          # 加粗
    min.size = 4                # 最小字体大小
  ) +
  # 使用渐变协调色系
  scale_fill_manual(values = gradient_colors) +
  # 主题美化
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )
LCJ_treeplot

YJ_treeplot <- ggplot(YJ_DS_average, 
                      aes(area = DS_average, fill = Mechnism, label = DS, subgroup = Mechnism)) +
  # 基础树状图
  geom_treemap(colour = "white", size = 1.2) +  # 白色边框
  # 子组边框（取消注释启用）
  # geom_treemap_subgroup_border(colour = "white", size = 2) +
  # 文本标签
  geom_treemap_text(
    aes(label = paste(DS, "\n",round(DS_average, 2))),
    colour = "black",           # 白色文字
    place = "centre",           # 居中显示
    # grow = TRUE,                # 允许文字扩展
    reflow = TRUE,              # 自动换行
    # family = "sans",            # 字体
    fontface = "plain",          # 加粗
    min.size = 4                # 最小字体大小
  ) +
  # 使用渐变协调色系
  scale_fill_manual(values = gradient_colors) +
  # 主题美化
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )
YJ_treeplot

Glacier_treeplot <- ggplot(Glacier_DS_average, 
                      aes(area = DS_average, fill = Mechnism, label = DS, subgroup = Mechnism)) +
  # 基础树状图
  geom_treemap(colour = "white", size = 1.2) +  # 白色边框
  # 子组边框（取消注释启用）
  # geom_treemap_subgroup_border(colour = "white", size = 2) +
  # 文本标签
  geom_treemap_text(
    aes(label = paste(DS, "\n",round(DS_average, 2))),
    colour = "black",           # 白色文字
    place = "centre",           # 居中显示
    # grow = TRUE,                # 允许文字扩展
    reflow = TRUE,              # 自动换行
    # family = "sans",            # 字体
    fontface = "plain",          # 加粗
    min.size = 4                # 最小字体大小
  ) +
  # 使用渐变协调色系
  scale_fill_manual(values = gradient_colors) +
  # 主题美化
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    legend.position = "none",
    legend.title = element_text(face = "bold")
  )
Glacier_treeplot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_DS_average_treemap_allecosystem.pdf" , width = 9,height = 6)
g <- plot_grid(Glacier_treeplot,HH_treeplot,CJ_treeplot,
               LCJ_treeplot,YJ_treeplot,
               labels = letters[1:5], ncol = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()









####环形柱状图########
# 筛选比例>1%的系统并计算标签位置
plot_data <- CJ_DS_proportion[CJ_DS_proportion$Proportion > 1, ] %>%
  mutate(
    # 按比例排序
    Defense_System = factor(Defense_System, 
                            levels = Defense_System[order(Proportion)]),
    # 计算每个柱子的中点角度（0-360度）
    angle = 90 - 360 * (as.numeric(Defense_System) - 0.5) / nrow(.),
# 调整标签水平对齐方式
hjust = ifelse(angle < -90, 1, 0),
# 将角度限制在-180到180度之间
angle = ifelse(angle < -90, angle + 180, angle),
# 标签y位置（柱子高度+偏移量）
label_y = Proportion + max(Proportion)*0.1
)

# 创建自定义纵轴
y_axis_breaks <- seq(0, max(plot_data$Proportion), by = 10)
y_axis_labels <- paste0(y_axis_breaks, "%")

# 创建环形条形图
circular_plot <- ggplot(plot_data, 
                        aes(x = Defense_System, 
                            y = Proportion, 
                            fill = Proportion)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +  # 使用柱状图
  
  annotate(
    "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
    fill = "white", color = 'black',linewidth = 0.1
  ) +
  
  
  scale_fill_gradient(low = "#BEBCDF", high = "#BEBCDF") +   # 颜色渐变
  
  # 极坐标转换
  coord_polar(start = -pi/nrow(plot_data)) +  # 调整起始角度使第一个柱子位于顶部
  
  # 设置y轴范围（中心留白）
  ylim(-max(plot_data$Proportion)*0.3, max(plot_data$Proportion)*1.2) +
  
  # 添加标签（自动调整角度）
  geom_text(aes(y = label_y, label = Defense_System, angle = angle, hjust = hjust),
            size = 2, color = "black", fontface = "plain") +
  
  # 标题和主题
  labs(title = "") +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 12, margin = margin(b = 10)),
    legend.position = "none",
    plot.margin = unit(c(1, 1, 1, 1), "cm")
  )

print(circular_plot)
# 创建基础图形

 