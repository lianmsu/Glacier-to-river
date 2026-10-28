# =============================================================================
# Glacier-to-river -- defensome landscape / Season comparison
# =============================================================================
# Compare DS number and density between sampling seasons (companion analysis
# to the region comparison).
#
# Manuscript      : -
# Main outputs    : 8genome_season.pdf, 08_season_comparison.RData
# Original script : 20gaoyuan/10Defensefinder/8defensefinder_res_season_compare.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
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
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'Yes',"Plateau",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'No',"Plain",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'glacier',"Glacier",MAG_metadata$Plateau)
colnames(MAG_metadata)[16] <- 'Region'

MAG_metadata$medium <- ifelse(MAG_metadata$Type == 'Cryoconite',"Glacier",MAG_metadata$Type)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Ice',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Snow',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'particle',"Plankton",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'water',"Plankton",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'sediment',"Benthon",MAG_metadata$medium)


genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
count_matrix <- count_matrix[rownames(count_matrix) %in% MAG_metadata$genome,]


MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]


#不同季节中MAG数量统计#####################3
#所有样品来自的国家的统计
MAG_metadata$Season <- 
  factor(MAG_metadata$Season,levels = c("Spring","Summer","Autumn","Winter"))
sample_info_season <- MAG_metadata %>% count(Season, name = "Count")

custom_colors <- c(
  "#fc8d62",  # 淡蓝色
  "#8da0cb",  # 淡橙色
  "#66c2a5",  # 淡绿色
  "#e78ac3"  # 淡红色
)
# 绘制水平柱状图
library(ggplot2)
sample_season <- 
  ggplot(sample_info_season, aes(x = Season, y = Count, fill = Season)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, color = "black") +  # 在柱子右侧显示数量
  coord_flip() +  # 将柱状图水平显示
  scale_y_continuous(limits = c(0, 4300), expand = c(0, 0)) +  # 设置 y 轴范围和扩展
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                               "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
  )) +  # 手动设置颜色
  labs(title = "",
       x = "",
       y = "Number of genomes") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中，加粗，字体大小为16
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    axis.text.x = element_blank(),  # 去除纵轴刻度标签
    axis.ticks.x = element_blank(),  # 去除纵轴刻度线
    #axis.ticks.y = element_line(color = "black"),  # 保留 x 轴刻度线
    legend.position = "none"  # 不显示图例
  )
sample_season
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
ggsave("8genome_season.pdf",
       sample_season,
       width = 4.5,height = 2.5,
       units = 'in')



#####1. 各个season的DS组成#########
DF_defense_2 <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]
DF_defense_2 <- merge(DF_defense,MAG_metadata,by.x = "genome",by.y = "genome")

# 统计每个medium中各个type的出现次数
result <- DF_defense_2 %>%
  count(Season, type) %>%  # 计算组合频数
  complete(Season, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  )

#Spring
Spring_DS_num <- result[result$Season == "Spring",]
library(reshape2)
Spring_DS_num <- melt(Spring_DS_num)
colnames(Spring_DS_num) <- c("Season","DS","num")
Spring_defense <- DF_defense_2[DF_defense_2$Season == "Spring",]
Spring_genome_DS <- count_matrix[rownames(count_matrix) %in% Spring_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Spring_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Spring_genome_DS)
})
# 转换为数据框并整理格式
Spring_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Spring_DS_proportion$Proportion <- Spring_DS_proportion$Proportion*100

#Summer
Summer_DS_num <- result[result$Season == "Summer",]
#library(reshape2)
Summer_DS_num <- melt(Summer_DS_num)
colnames(Summer_DS_num) <- c("Season","DS","num")
Summer_defense <- DF_defense_2[DF_defense_2$Season == "Summer",]
Summer_genome_DS <- count_matrix[rownames(count_matrix) %in% Summer_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Summer_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Summer_genome_DS)
})
# 转换为数据框并整理格式
Summer_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Summer_DS_proportion$Proportion <- Summer_DS_proportion$Proportion*100

#Autumn
Autumn_DS_num <- result[result$Season == "Autumn",]
#library(reshape2)
Autumn_DS_num <- melt(Autumn_DS_num)
colnames(Autumn_DS_num) <- c("Season","DS","num")
Autumn_defense <- DF_defense_2[DF_defense_2$Season == "Autumn",]
Autumn_genome_DS <- count_matrix[rownames(count_matrix) %in% Autumn_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Autumn_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Autumn_genome_DS)
})
# 转换为数据框并整理格式
Autumn_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Autumn_DS_proportion$Proportion <- Autumn_DS_proportion$Proportion*100


library(ggplot2)
# Spring
Spring_proportion_plot <- ggplot(Spring_DS_proportion[c(1:25),], 
                                 aes(x = reorder(Defense_System, -Proportion), 
                                         y = Proportion)) +
  geom_bar(stat = "identity", fill = "#fc8d62", width = 0.7) +
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

print(Spring_proportion_plot)

# Summer
Summer_proportion_plot <- 
  ggplot(Summer_DS_proportion[c(1:25),], 
         aes(x = reorder(Defense_System, -Proportion), 
                                                                                   y = Proportion)) +
  geom_bar(stat = "identity", fill = "#8da0cb", width = 0.7) +
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

print(Summer_proportion_plot)

# Autumn
Autumn_proportion_plot <- 
  ggplot(Autumn_DS_proportion[c(1:25),], 
         aes(x = reorder(Defense_System, -Proportion), 
             y = Proportion)) +
  geom_bar(stat = "identity", fill = "#66c2a5", width = 0.7) +
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

print(Autumn_proportion_plot)

library(Cairo)
library(cowplot)
cairo_pdf(filename = "8DS_proportion.pdf" , width = 3.05,height = 6.86)
g <- plot_grid(Spring_proportion_plot,Summer_proportion_plot,
               Autumn_proportion_plot,
               labels = letters[1:3], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

###2.MAG数量随DS数量的变化##############
# Spring
Spring_genome_DSnum <- 
  rowSums(Spring_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Spring_genome_DSnum$genome <- 
  rownames(Spring_genome_DSnum)
Spring_MAG_metadata <- MAG_metadata[MAG_metadata$Season == "Spring",]
Spring_MAG_metadata <- left_join(Spring_MAG_metadata, 
                                 Spring_genome_DSnum, 
                             by = c("genome" = "genome"))
Spring_MAG_metadata$DS_num[is.na(Spring_MAG_metadata$DS_num)] <- 0
Spring_MAG_DSnum <- select(Spring_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Spring_count_data <- Spring_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Summer
Summer_genome_DSnum <- rowSums(Summer_genome_DS) %>% 
  as.data.frame() %>% setNames("DS_num")
Summer_genome_DSnum$genome <- rownames(Summer_genome_DSnum)
Summer_MAG_metadata <- MAG_metadata[MAG_metadata$Season == "Summer",]
Summer_MAG_metadata <- left_join(Summer_MAG_metadata, 
                                 Summer_genome_DSnum, 
                             by = c("genome" = "genome"))
Summer_MAG_metadata$DS_num[is.na(Summer_MAG_metadata$DS_num)] <- 0
Summer_MAG_DSnum <- select(Summer_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Summer_count_data <- Summer_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Autumn
Autumn_genome_DSnum <- rowSums(Autumn_genome_DS) %>% 
  as.data.frame() %>% setNames("DS_num")
Autumn_genome_DSnum$genome <- rownames(Autumn_genome_DSnum)
Autumn_MAG_metadata <- MAG_metadata[MAG_metadata$Season == "Autumn",]
Autumn_MAG_metadata <- left_join(Autumn_MAG_metadata, 
                                 Autumn_genome_DSnum, 
                             by = c("genome" = "genome"))
Autumn_MAG_metadata$DS_num[is.na(Autumn_MAG_metadata$DS_num)] <- 0
Autumn_MAG_DSnum <- select(Autumn_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Autumn_count_data <- Autumn_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

library(ggplot2)
#Spring
Spring_magnum_dsnum <- 
  ggplot(Spring_count_data, aes(x = factor(DS_num), y = genome_count, 
                                       fill = "#fc8d62")) +
  geom_col(width = 0.7, fill = "#fc8d62", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
      #subtitle = paste("Total", sum(Plankton_count_data$genome_count), "genomes"),
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
Spring_magnum_dsnum

#Summer
Summer_magnum_dsnum <- 
  ggplot(Summer_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#8da0cb")) +
  geom_col(width = 0.7, fill = "#8da0cb", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
       #subtitle = paste("Total", sum(Plain_count_data$genome_count), "genomes"),
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
Summer_magnum_dsnum


#Autumn
Autumn_magnum_dsnum <- 
  ggplot(Autumn_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#66c2a5")) +
  geom_col(width = 0.7, fill = "#66c2a5", linewidth = 0.3) +
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
Autumn_magnum_dsnum

library(Cairo)
library(cowplot)
cairo_pdf(filename = "8_MAGnum_DSnum.pdf" , width = 4.33,height = 4.5)
g2 <- plot_grid(Spring_magnum_dsnum,Summer_magnum_dsnum,
               Autumn_magnum_dsnum,
               labels = letters[7:11], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###3.数量与MAG大小关系###############

#Spring
Spring_MAG_quality <- 
  MAG_quality[MAG_quality$Name %in% Spring_MAG_DSnum$genome,]

Spring_MAG_DSdensity <- 
  merge(Spring_MAG_DSnum, Spring_MAG_quality, by.x = "genome",by.y = "Name")
Spring_MAG_DSdensity$MAG_size <- Spring_MAG_DSdensity$Genome_Size/1000000
Spring_MAG_DSdensity$DS_density <- 
  Spring_MAG_DSdensity$DS_num/Spring_MAG_DSdensity$Genome_Size*1000000

Spring_MAG_DSdensity$magsize <- floor(Spring_MAG_DSdensity$MAG_size)

library(dplyr)
Spring_df_summary <- Spring_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Spring_df_summary2 <- Spring_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Spring_df_summary$magsize, 
         Spring_df_summary$mean_value, 
         alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Spring_df_summary2$magsize, 
         Spring_df_summary2$mean_value,
         alternative = "two.sided",method = "spearman",conf.level = 0.95)

Spring_df_summary$Season <- 'Spring'
Spring_df_summary2$Season <- 'Spring'

# Summer
Summer_MAG_quality <- MAG_quality[MAG_quality$Name %in% Summer_MAG_DSnum$genome,]

Summer_MAG_DSdensity <- 
  merge(Summer_MAG_DSnum, Summer_MAG_quality, by.x = "genome",by.y = "Name")
Summer_MAG_DSdensity$MAG_size <- Summer_MAG_DSdensity$Genome_Size/1000000
Summer_MAG_DSdensity$DS_density <- 
  Summer_MAG_DSdensity$DS_num/Summer_MAG_DSdensity$Genome_Size*1000000

Summer_MAG_DSdensity$magsize <- floor(Summer_MAG_DSdensity$MAG_size)

library(dplyr)
Summer_df_summary <- Summer_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Summer_df_summary2 <- Summer_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Summer_df_summary$magsize, Summer_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Summer_df_summary2$magsize, Summer_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Summer_df_summary$Season <- 'Summer'
Summer_df_summary2$Season <- 'Summer'

#Autumn
Autumn_MAG_quality <- MAG_quality[MAG_quality$Name %in% Autumn_MAG_DSnum$genome,]

Autumn_MAG_DSdensity <- merge(Autumn_MAG_DSnum, Autumn_MAG_quality, by.x = "genome",by.y = "Name")
Autumn_MAG_DSdensity$MAG_size <- Autumn_MAG_DSdensity$Genome_Size/1000000
Autumn_MAG_DSdensity$DS_density <- 
  Autumn_MAG_DSdensity$DS_num/Autumn_MAG_DSdensity$Genome_Size*1000000

Autumn_MAG_DSdensity$magsize <- floor(Autumn_MAG_DSdensity$MAG_size)

library(dplyr)
Autumn_df_summary <- Autumn_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Autumn_df_summary2 <- Autumn_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Autumn_df_summary$magsize, Autumn_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Autumn_df_summary2$magsize, Autumn_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Autumn_df_summary$Season <- 'Autumn'
Autumn_df_summary2$Season <- 'Autumn'


DSnum_allseason <- rbind(Spring_df_summary2, Summer_df_summary2,
                            Autumn_df_summary2)
DSdensity_allseason <- rbind(Spring_df_summary,Summer_df_summary,
                            Autumn_df_summary)



library(ggplot2)
DSnum_allseason_plot <- 
ggplot(DSnum_allseason, aes(x = magsize, y = mean_value,color = Season)) +
  geom_point(aes(color=Season))+
  geom_line(aes(color=Season))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = Season),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                               "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
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
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
DSnum_allseason_plot

DSdensity_allseason_plot <- 
  ggplot(DSdensity_allseason, 
         aes(x = magsize, y = mean_value,color = Season)) +
  geom_point(aes(color=Season))+
  geom_line(aes(color=Season))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = Season),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                                "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
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
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
DSdensity_allseason_plot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "8_DS_num-density_allseason.pdf" , width = 4.14,height = 5.25)
g3 <- plot_grid(DSnum_allseason_plot,DSdensity_allseason_plot,
               labels = letters[8], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g3
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###4.比较不同环境的Defese family多少，DS per MAG###############

Spring_DF_num <- sum(Spring_DS_num$num>0)

Summer_DF_num <- sum(Summer_DS_num$num>0)

Autumn_DF_num <- sum(Autumn_DS_num$num>0)

allseason_DF_num <- 
  data.frame(Season = c('Spring','Summer','Autumn'),
             DF_num = c(Spring_DF_num, Summer_DF_num, Autumn_DF_num),
            DS_mean = c(sum(Spring_MAG_metadata$DS_num)/nrow(Spring_MAG_metadata),
             sum(Summer_MAG_metadata$DS_num)/nrow(Summer_MAG_metadata),
              sum(Autumn_MAG_metadata$DS_num)/nrow(Autumn_MAG_metadata))
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
p_compare_DStype <- 
  ggplot(allseason_DF_num, aes(x = reorder(Season,DF_num) ,y = DF_num)) +
  geom_bar(stat = "identity", aes(fill = Season),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  # scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +  # 手动设置颜色
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                                "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
  )) +  # 手动设置颜色
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  # geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DF_num,2)),vjust = 0.5,hjust = 1.1,size = 4)+
  labs(x = "", y = "Number of DS type", title = "")+  # 设置标签和标题
  mytheme
p_compare_DStype

p_compare_DSmean <- 
  ggplot(allseason_DF_num, aes(x = reorder(Season,DS_mean) ,y = DS_mean)) +
  geom_bar(stat = "identity", aes(fill = Season),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                               "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
  )) +  # 手动设置颜色
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  # geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  # geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DS_mean,2)),vjust = 0.5,hjust = 1.1,size = 4)+
  labs(x = "", y = "Number of DS type", title = "")+  # 设置标签和标题
  mytheme
p_compare_DSmean

library(Cairo)
library(cowplot)
cairo_pdf(filename = "8_DStype_DSmean_allseason.pdf" , width = 2,height = 4.6)
g4 <- plot_grid(p_compare_DStype,p_compare_DSmean,
               labels = letters[9], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g4
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###5.DS作用机制统计############
library(treemap)
library(RColorBrewer)
library(treemapify)
library(ggplot2)

DF_mechanism <- read.delim2('df_List_Systems_lpw.csv',sep = ',')
DF_mechanism$Effector[DF_mechanism$Effector == ""] <- "Unknown"
DF_mechanism$Mechnism[DF_mechanism$Mechnism == ""] <- "Unknown"

Spring_DS_average <- Spring_DS_num
Spring_DS_average$DS_average <- 
  (Spring_DS_average$num) / nrow(Spring_MAG_metadata) 
Spring_DS_average <- left_join(
  Spring_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Spring_DS_average[is.na(Spring_DS_average)] <- "Unknown"

Summer_DS_average <- Summer_DS_num
Summer_DS_average$DS_average <- 
  (Summer_DS_average$num) / nrow(Summer_MAG_metadata) 
Summer_DS_average <- left_join(
  Summer_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Summer_DS_average[is.na(Summer_DS_average)] <- "Unknown"

Autumn_DS_average <- Autumn_DS_num
Autumn_DS_average$DS_average <- 
  (Autumn_DS_average$num) / nrow(Autumn_MAG_metadata) 
Autumn_DS_average <- left_join(
  Autumn_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Autumn_DS_average[is.na(Autumn_DS_average)] <- "Unknown"


gradient_colors <- c("#FFA5B8", "#8AC7F5", "#FFE4A3", "#9ADFDF","#C4D8E9", "#C0C0C0")

Spring_treeplot <- ggplot(Spring_DS_average, 
       aes(area = DS_average, fill = Mechnism, 
           label = DS, subgroup = Mechnism)) +
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
Spring_treeplot

Summer_treeplot <- ggplot(Summer_DS_average, 
                      aes(area = DS_average, fill = Mechnism, 
                          label = DS, subgroup = Mechnism)) +
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
Summer_treeplot


Autumn_treeplot <- ggplot(Autumn_DS_average, 
                      aes(area = DS_average, 
                          fill = Mechnism, label = DS, 
                          subgroup = Mechnism)) +
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
Autumn_treeplot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "8_DS_average_treemap_allRegion.pdf" , 
          width = 6,height = 6)
g5 <- plot_grid(Spring_treeplot,Summer_treeplot,Autumn_treeplot,
               labels = letters[1:3], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g5
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

# 
# Region_MAG_num <- table(MAG_metadata$Region) %>% as.data.frame()
# colnames(Region_MAG_num) <- c("Region","Count")
# Region_MAG_num$Region <- factor(Region_MAG_num$Region,levels = c("Plain","Plateau","Glacier"))
# # 绘制水平柱状图
# sample_Region <- 
#   ggplot(Region_MAG_num, aes(x = Region, y = Count, fill = Region)) +
#   geom_bar(stat = "identity",width = 0.8) +
#   geom_text(aes(label = Count), hjust = -0.2, size = 3.5, color = "black") +  # 在柱子右侧显示数量
#   coord_flip() +  # 将柱状图水平显示
#   scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +  # 设置 y 轴范围和扩展
#   scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
#                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
#                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
#   )) +  # 手动设置颜色
#   labs(title = "",
#        x = "",
#        y = "Number of genomes") +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中，加粗，字体大小为16
#     panel.grid.major = element_blank(),  # 去除主要网格线
#     panel.grid.minor = element_blank(),  # 去除次要网格线
#     axis.text.x = element_blank(),  # 去除纵轴刻度标签
#     axis.ticks.x = element_blank(),  # 去除纵轴刻度线
#     #axis.ticks.y = element_line(color = "black"),  # 保留 x 轴刻度线
#     legend.position = "none"  # 不显示图例
#   )
# sample_Region
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
# ggsave("2_genome_Region.pdf",sample_Region,width = 3,height = 1.7,units = 'in')


#6.三种season的DS number比较####
genome_defense_num <- 
  genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
MAG_metadata <- merge(MAG_metadata, genome_defense_num, by = "genome", all = T)
MAG_metadata[,63][is.na(MAG_metadata[,63])] <- 0

MAG_DS_metadata_season <- 
  select(MAG_metadata,genome,ecosystem,Type,Region,Season,DS_num)
MAG_DS_metadata_season$Season <- 
  as.character(MAG_DS_metadata_season$Season) 
MAG_DS_metadata_season <- MAG_DS_metadata_season[MAG_DS_metadata_season$Season %in% 
                                                   c("Spring","Summer","Autumn"),]
MAG_DS_metadata_season$Season <- factor(MAG_DS_metadata_season$Season,
                                levels = c("Spring","Summer","Autumn"))
# MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
#                                             MAG_metadata_plateau_plain$Region)
library(rstatix)
md_p_val_DS <- MAG_DS_metadata_season %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DS_num ~ Season) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Season",dodge=0.5)
md_p_val_DS$y.position <- md_p_val_DS$y.position - 42

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)
p_DS_num_season <- 
  ggplot(MAG_DS_metadata_season, 
         aes(x = Season, y = DS_num)) +
  geom_half_boxplot(aes(fill = Season), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color=Season),
                  width = 0.2,size = 0.5, alpha = 0.7)+
  stat_pvalue_manual(md_p_val_DS,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="DS number (per MAG)")+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                               "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
  )) +  # 手动设置颜色
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA","Spring"="#fc8d62",
                               "Summer"="#8da0cb","Autumn"="#66c2a5","Winter"="#e78ac3"
  )) +  # 手动设置颜色
  #修改主题
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12, face = "plain"),  # 标题居中，加粗，字体大小为16
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),  # 设置图形边距
    plot.background = element_blank(),  # 去除图形背景
    panel.background = element_blank(),  # 去除面板背景
    panel.border = element_rect(linewidth = 0.6),  # 设置面板边框
    panel.grid = element_blank(),  # 去除网格线
    axis.title.x = element_text(size = 10, face = "plain",color = 'black'),  # 设置 x 轴标题
    axis.title.y = element_text(size = 10, face = "plain",color = 'black'),  # 设置 y 轴标题
    axis.ticks.x = element_line(linewidth = 0.5),  # 设置 x 轴刻度线
    axis.ticks.y = element_line(linewidth = 0.5),  # 设置 y 轴刻度线
    axis.ticks.length.x = unit(0.1, "cm"),  # 设置 x 轴刻度长度
    axis.ticks.length.y = unit(0.1, "cm"),  # 设置 y 轴刻度长度
    axis.text.x = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  ) +
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)))   +
  coord_cartesian(ylim = c(-1, 20))
p_DS_num_season


library(Cairo)
library(cowplot)
cairo_pdf(filename = "8_DSseasonDStype_DSmean.pdf" , width = 7.87,height = 3)
g6 <- plot_grid(p_DS_num_season, 
               DSnum_allseason_plot,DSdensity_allseason_plot,
               labels = letters[8], ncol = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g6
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()



MAG_DS_type_num <- 
  count_matrix
MAG_DS_type_num$MAG <- rownames(MAG_DS_type_num)
write.table(MAG_DS_type_num, "8_MAG_DS_type_num.tsv",
            sep = "\t", row.names = F)

save.image("./08_season_comparison.RData")
load("./08_season_comparison.RData")

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

 