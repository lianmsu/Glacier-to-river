# =============================================================================
# Glacier-to-river -- defensome landscape / Medium / lifestyle comparison
# =============================================================================
# Compare DS content across snow, ice, cryoconite, water, particle and
# sediment, and between planktonic and benthic lifestyles.
#
# Manuscript      : Fig. 3b, Fig. 3d, Fig. S4
# Main outputs    : 3genome_medium.pdf, 2_genome_Region.pdf
# Original script : 20gaoyuan/10Defensefinder/3defensefinder_res_plankton_benthon.R
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

MAG_metadata <- read.csv(file.path(GTR_DATA, "20gaoyuan/01MAG_info/5455MAG_metadata.csv"))
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


#不同介质中MAG数量统计#####################3
#所有样品来自的国家的统计
MAG_metadata$medium <- factor(MAG_metadata$medium,levels = c("Benthon","Plankton","Glacier"))
sample_info_medium <- MAG_metadata %>% count(medium, name = "Count")

custom_colors <- c(
  "#9ecae1",  # 淡蓝色
  "#fdae6b",  # 淡橙色
  "#a1d99b",  # 淡绿色
  "#fc9272",  # 淡红色
  "#bcbddc",  # 淡紫色
  "#ffeda0",  # 淡黄色
  "#c6dbef",  # 浅蓝色
  "#fdd0a2",  # 浅橙色
  "#c7e9c0",  # 浅绿色
  "#fcbba1",  # 浅红色
  "#dadaeb",  # 浅紫色
  "#ffa69e"   # 浅黄色
)
# 绘制水平柱状图
library(ggplot2)
sample_medium <- 
  ggplot(sample_info_medium, aes(x = medium, y = Count, fill = medium)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, color = "black") +  # 在柱子右侧显示数量
  coord_flip() +  # 将柱状图水平显示
  scale_y_continuous(limits = c(0, 4300), expand = c(0, 0)) +  # 设置 y 轴范围和扩展
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
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
sample_medium
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
ggsave("3genome_medium.pdf",sample_medium,width = 4.5,height = 2.5,units = 'in')










#####1. 各个medium的DS组成#########
DF_defense_2 <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]
DF_defense_2 <- merge(DF_defense,MAG_metadata,by.x = "genome",by.y = "genome")

# 统计每个medium中各个type的出现次数
result <- DF_defense_2 %>%
  count(medium, type) %>%  # 计算组合频数
  complete(medium, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  )

#Glacier
Glacier_DS_num <- result[result$medium == "Glacier",]
library(reshape2)
Glacier_DS_num <- melt(Glacier_DS_num)
colnames(Glacier_DS_num) <- c("medium","DS","num")
Glacier_defense <- DF_defense_2[DF_defense_2$medium == "Glacier",]
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

#Plankton
Plankton_DS_num <- result[result$medium == "Plankton",]
#library(reshape2)
Plankton_DS_num <- melt(Plankton_DS_num)
colnames(Plankton_DS_num) <- c("medium","DS","num")
Plankton_defense <- DF_defense_2[DF_defense_2$medium == "Plankton",]
Plankton_genome_DS <- count_matrix[rownames(count_matrix) %in% Plankton_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Plankton_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Plankton_genome_DS)
})
# 转换为数据框并整理格式
Plankton_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Plankton_DS_proportion$Proportion <- Plankton_DS_proportion$Proportion*100

#Benthon
Benthon_DS_num <- result[result$medium == "Benthon",]
#library(reshape2)
Benthon_DS_num <- melt(Benthon_DS_num)
colnames(Benthon_DS_num) <- c("medium","DS","num")
Benthon_defense <- DF_defense_2[DF_defense_2$medium == "Benthon",]
Benthon_genome_DS <- count_matrix[rownames(count_matrix) %in% Benthon_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Benthon_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Benthon_genome_DS)
})
# 转换为数据框并整理格式
Benthon_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Benthon_DS_proportion$Proportion <- Benthon_DS_proportion$Proportion*100


library(ggplot2)
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

# Plankton
Plankton_proportion_plot <- ggplot(Plankton_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                   y = Proportion)) +
  geom_bar(stat = "identity", fill = "#c7e9c0", width = 0.7) +
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

print(Plankton_proportion_plot)

# Benthon
Benthon_proportion_plot <- ggplot(Benthon_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                      y = Proportion)) +
  geom_bar(stat = "identity", fill = "#c38e70", width = 0.7) +
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

print(Benthon_proportion_plot)

library(Cairo)
library(cowplot)
cairo_pdf(filename = "3DS_proportion.pdf" , width = 3.05,height = 6.86)
g <- plot_grid(Glacier_proportion_plot,Plankton_proportion_plot,
               Benthon_proportion_plot,
               
               labels = letters[1:5], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

###2.MAG数量随DS数量的变化##############
#Plankton
Plankton_genome_DSnum <- rowSums(Plankton_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Plankton_genome_DSnum$genome <- rownames(Plankton_genome_DSnum)
Plankton_MAG_metadata <- MAG_metadata[MAG_metadata$medium == "Plankton",]
Plankton_MAG_metadata <- left_join(Plankton_MAG_metadata, 
                             Plankton_genome_DSnum, 
                             by = c("genome" = "genome"))
Plankton_MAG_metadata$DS_num[is.na(Plankton_MAG_metadata$DS_num)] <- 0
Plankton_MAG_DSnum <- select(Plankton_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Plankton_count_data <- Plankton_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Benthon
Benthon_genome_DSnum <- rowSums(Benthon_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Benthon_genome_DSnum$genome <- rownames(Benthon_genome_DSnum)
Benthon_MAG_metadata <- MAG_metadata[MAG_metadata$medium == "Benthon",]
Benthon_MAG_metadata <- left_join(Benthon_MAG_metadata, 
                             Benthon_genome_DSnum, 
                             by = c("genome" = "genome"))
Benthon_MAG_metadata$DS_num[is.na(Benthon_MAG_metadata$DS_num)] <- 0
Benthon_MAG_DSnum <- select(Benthon_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Benthon_count_data <- Benthon_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Glacier
Glacier_genome_DSnum <- rowSums(Glacier_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Glacier_genome_DSnum$genome <- rownames(Glacier_genome_DSnum)
Glacier_MAG_metadata <- MAG_metadata[MAG_metadata$medium == "Glacier",]
Glacier_MAG_metadata <- left_join(Glacier_MAG_metadata, 
                             Glacier_genome_DSnum, 
                             by = c("genome" = "genome"))
Glacier_MAG_metadata$DS_num[is.na(Glacier_MAG_metadata$DS_num)] <- 0
Glacier_MAG_DSnum <- select(Glacier_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Glacier_count_data <- Glacier_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

library(ggplot2)
#Plankton
Plankton_magnum_dsnum <- ggplot(Plankton_count_data, aes(x = factor(DS_num), y = genome_count, 
                                       fill = "#c7e9c0")) +
  geom_col(width = 0.7, fill = "#c7e9c0", linewidth = 0.3) +
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
Plankton_magnum_dsnum

#Benthon
Benthon_magnum_dsnum <- ggplot(Benthon_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#c38e70")) +
  geom_col(width = 0.7, fill = "#c38e70", linewidth = 0.3) +
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
Benthon_magnum_dsnum


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
cairo_pdf(filename = "3_MAGnum_DSnum.pdf" , width = 4.33,height = 4.5)
g <- plot_grid(Glacier_magnum_dsnum,Plankton_magnum_dsnum,
               Benthon_magnum_dsnum,
               labels = letters[7:11], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###3.数量与MAG大小关系###############

#Plankton
Plankton_MAG_quality <- MAG_quality[MAG_quality$Name %in% Plankton_MAG_DSnum$genome,]

Plankton_MAG_DSdensity <- merge(Plankton_MAG_DSnum, Plankton_MAG_quality, by.x = "genome",by.y = "Name")
Plankton_MAG_DSdensity$MAG_size <- Plankton_MAG_DSdensity$Genome_Size/1000000
Plankton_MAG_DSdensity$DS_density <- Plankton_MAG_DSdensity$DS_num/Plankton_MAG_DSdensity$Genome_Size*1000000

Plankton_MAG_DSdensity$magsize <- floor(Plankton_MAG_DSdensity$MAG_size)

library(dplyr)
Plankton_df_summary <- Plankton_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Plankton_df_summary2 <- Plankton_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Plankton_df_summary$magsize, Plankton_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Plankton_df_summary2$magsize, Plankton_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Plankton_df_summary$medium <- 'Plankton'
Plankton_df_summary2$medium <- 'Plankton'

#Benthon
Benthon_MAG_quality <- MAG_quality[MAG_quality$Name %in% Benthon_MAG_DSnum$genome,]

Benthon_MAG_DSdensity <- merge(Benthon_MAG_DSnum, Benthon_MAG_quality, by.x = "genome",by.y = "Name")
Benthon_MAG_DSdensity$MAG_size <- Benthon_MAG_DSdensity$Genome_Size/1000000
Benthon_MAG_DSdensity$DS_density <- Benthon_MAG_DSdensity$DS_num/Benthon_MAG_DSdensity$Genome_Size*1000000

Benthon_MAG_DSdensity$magsize <- floor(Benthon_MAG_DSdensity$MAG_size)

library(dplyr)
Benthon_df_summary <- Benthon_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Benthon_df_summary2 <- Benthon_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Benthon_df_summary$magsize, Benthon_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Benthon_df_summary2$magsize, Benthon_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Benthon_df_summary$medium <- 'Benthon'
Benthon_df_summary2$medium <- 'Benthon'

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

Glacier_df_summary$medium <- 'Glacier'
Glacier_df_summary2$medium <- 'Glacier'


DSnum_allmedium <- rbind(Glacier_df_summary2,Plankton_df_summary2,
                            Benthon_df_summary2)
DSdensity_allmedium <- rbind(Glacier_df_summary,Plankton_df_summary,
                            Benthon_df_summary)



library(ggplot2)
DSnum_allmedium_plot <- 
ggplot(DSnum_allmedium, aes(x = magsize, y = mean_value,color = medium)) +
  geom_point(aes(color=medium))+
  geom_line(aes(color=medium))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = medium),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
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
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
DSnum_allmedium_plot

DSdensity_allmedium_plot <- 
  ggplot(DSdensity_allmedium, aes(x = magsize, y = mean_value,color = medium)) +
  geom_point(aes(color=medium))+
  geom_line(aes(color=medium))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = medium),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
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
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
DSdensity_allmedium_plot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "3_DS_num-density_allecosystem.pdf" , width = 4.14,height = 5.25)
g <- plot_grid(DSnum_allmedium_plot,DSdensity_allmedium_plot,
               labels = letters[8], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###4.比较不同环境的Defese family多少，DS per MAG###############

Plankton_DF_num <- sum(Plankton_DS_num$num>0)

Benthon_DF_num <- sum(Benthon_DS_num$num>0)

Glacier_DF_num <- sum(Glacier_DS_num$num>0)

allmedium_DF_num <- data.frame(medium = c('Glacier','Plankton','Benthon'),
                                  DF_num = c(Glacier_DF_num, Plankton_DF_num, Benthon_DF_num),
                                  DS_mean = c(sum(Glacier_MAG_metadata$DS_num)/nrow(Glacier_MAG_metadata),
                                              sum(Plankton_MAG_metadata$DS_num)/nrow(Plankton_MAG_metadata),
                                              sum(Benthon_MAG_metadata$DS_num)/nrow(Benthon_MAG_metadata))
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
p_compare_DStype <- ggplot(allmedium_DF_num, aes(x = reorder(medium,DF_num) ,y = DF_num)) +
  geom_bar(stat = "identity", aes(fill = medium),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DF_num,2)),vjust = 0.5,hjust = 1.1,size = 4)+
  labs(x = "", y = "Number of DS type", title = "")+  # 设置标签和标题
  mytheme
p_compare_DStype

p_compare_DSmean <- ggplot(allmedium_DF_num, aes(x = reorder(medium,DS_mean) ,y = DS_mean)) +
  geom_bar(stat = "identity", aes(fill = medium),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
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
cairo_pdf(filename = "3_DStype_DSmean_allRegion.pdf" , width = 2,height = 4.6)
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

Plankton_DS_average <- Plankton_DS_num
Plankton_DS_average$DS_average <- (Plankton_DS_average$num) / nrow(Plankton_MAG_metadata) 
Plankton_DS_average <- left_join(
  Plankton_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Plankton_DS_average[is.na(Plankton_DS_average)] <- "Unknown"

Benthon_DS_average <- Benthon_DS_num
Benthon_DS_average$DS_average <- (Benthon_DS_average$num) / nrow(Benthon_MAG_metadata) 
Benthon_DS_average <- left_join(
  Benthon_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Benthon_DS_average[is.na(Benthon_DS_average)] <- "Unknown"

Glacier_DS_average <- Glacier_DS_num
Glacier_DS_average$DS_average <- (Glacier_DS_average$num) / nrow(Glacier_MAG_metadata) 
Glacier_DS_average <- left_join(
  Glacier_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Glacier_DS_average[is.na(Glacier_DS_average)] <- "Unknown"


gradient_colors <- c("#FFA5B8", "#8AC7F5", "#FFE4A3", "#9ADFDF","#C4D8E9", "#C0C0C0")

Plankton_treeplot <- ggplot(Plankton_DS_average, 
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
Plankton_treeplot

Benthon_treeplot <- ggplot(Benthon_DS_average, 
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
Benthon_treeplot


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
cairo_pdf(filename = "3_DS_average_treemap_allRegion.pdf" , width = 6,height = 6)
g <- plot_grid(Glacier_treeplot,Plankton_treeplot,Benthon_treeplot,
               labels = letters[1:3], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


Region_MAG_num <- table(MAG_metadata$Region) %>% as.data.frame()
colnames(Region_MAG_num) <- c("Region","Count")
Region_MAG_num$Region <- factor(Region_MAG_num$Region,levels = c("Plain","Plateau","Glacier"))
# 绘制水平柱状图
sample_Region <- 
  ggplot(Region_MAG_num, aes(x = Region, y = Count, fill = Region)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, color = "black") +  # 在柱子右侧显示数量
  coord_flip() +  # 将柱状图水平显示
  scale_y_continuous(limits = c(0, 4000), expand = c(0, 0)) +  # 设置 y 轴范围和扩展
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
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
sample_Region
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
ggsave("2_genome_Region.pdf",sample_Region,width = 3,height = 1.7,units = 'in')


#6.三种region的DS number比较####
genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
MAG_metadata <- merge(MAG_metadata, genome_defense_num, by = "genome", all = T)
MAG_metadata[,63][is.na(MAG_metadata[,63])] <- 0

MAG_DS_metadata_glacier_Benthon_plankton <- select(MAG_metadata,genome,ecosystem,Type,Region,medium,DS_num)
MAG_DS_metadata_glacier_Benthon_plankton$medium <- as.character(MAG_DS_metadata_glacier_Benthon_plankton$medium) 
MAG_DS_metadata_glacier_Benthon_plankton$medium <- factor(MAG_DS_metadata_glacier_Benthon_plankton$medium,
                                                       levels = c("Glacier","Benthon","Plankton"))
# MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
#                                             MAG_metadata_plateau_plain$Region)
library(rstatix)
md_p_val_DS <- MAG_DS_metadata_glacier_Benthon_plankton %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DS_num ~ medium) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="medium",dodge=0.5)
md_p_val_DS$y.position <- md_p_val_DS$y.position - 45

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)
p_DS_glacier_Benthon_plankton <- 
  ggplot(MAG_DS_metadata_glacier_Benthon_plankton, 
         aes(x = medium, y = DS_num)) +
  geom_half_boxplot(aes(fill = medium), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_half_point(aes(color=medium),
                  side = "r",size = 0.5,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(md_p_val_DS,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="DS number (per MAG)")+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0","Plateau"= "#C4D8E9",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70","Plain"="#d4e09b",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0","Plateau"= "#C4D8E9",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70","Plain"="#d4e09b",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
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
p_DS_glacier_Benthon_plankton


library(Cairo)
library(cowplot)
cairo_pdf(filename = "3_DSmedium_DStype_DSmean.pdf" , width = 7.87,height = 3)
g <- plot_grid(p_DS_glacier_Benthon_plankton, 
               DSnum_allmedium_plot,DSdensity_allmedium_plot,
               labels = letters[8], ncol = 3
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

 