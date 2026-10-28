# =============================================================================
# Glacier-to-river -- defensome landscape / Region comparison (glacier / plateau / plain)
# =============================================================================
# Compare DS number, density and diversity between the glacier, plateau and
# plain regions, with Wilcoxon/Mann-Whitney tests.
#
# Manuscript      : Fig. 3a, Fig. 3e, Fig. S2
# Main outputs    : 2_genome_Region.pdf, 2_two_barplot.pdf, 02_region_comparison.RData
# Original script : 20gaoyuan/10Defensefinder/2defensefinder_res_Glacier_Plateau_Plain.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). The default
# is the current working directory; set the environment variable to point
# somewhere else instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)

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
colnames(MAG_metadata)[17] <- 'Region'

genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
count_matrix <- count_matrix[rownames(count_matrix) %in% MAG_metadata$genome,]


MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]



#####1. 各个region的DS组成#########
DF_defense_2 <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]
DF_defense_2 <- merge(DF_defense,MAG_metadata,by.x = "genome",by.y = "genome")

# 统计每个region中各个type的出现次数
result <- DF_defense_2 %>%
  count(Region, type) %>%  # 计算组合频数
  complete(Region, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  )

#Glacier
Glacier_DS_num <- result[result$Region == "Glacier",]
library(reshape2)
Glacier_DS_num <- melt(Glacier_DS_num)
colnames(Glacier_DS_num) <- c("Region","DS","num")
Glacier_defense <- DF_defense_2[DF_defense_2$Region == "Glacier",]
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

#Plateau
Plateau_DS_num <- result[result$Region == "Plateau",]
#library(reshape2)
Plateau_DS_num <- melt(Plateau_DS_num)
colnames(Plateau_DS_num) <- c("Region","DS","num")
Plateau_defense <- DF_defense_2[DF_defense_2$Region == "Plateau",]
Plateau_genome_DS <- count_matrix[rownames(count_matrix) %in% Plateau_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Plateau_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Plateau_genome_DS)
})
# 转换为数据框并整理格式
Plateau_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Plateau_DS_proportion$Proportion <- Plateau_DS_proportion$Proportion*100

#Plain
Plain_DS_num <- result[result$Region == "Plain",]
#library(reshape2)
Plain_DS_num <- melt(Plain_DS_num)
colnames(Plain_DS_num) <- c("Region","DS","num")
Plain_defense <- DF_defense_2[DF_defense_2$Region == "Plain",]
Plain_genome_DS <- count_matrix[rownames(count_matrix) %in% Plain_defense$genome,]
# 计算每列非零比例
non_zero_proportions <- apply(Plain_genome_DS, 2, function(x) {
  sum(x != 0) / nrow(Plain_genome_DS)
})
# 转换为数据框并整理格式
Plain_DS_proportion <- data.frame(
  Defense_System = names(non_zero_proportions),
  Proportion = round(non_zero_proportions, 4),
  row.names = NULL
) %>%
  arrange(desc(Proportion))  # 按比例降序排列
Plain_DS_proportion$Proportion <- Plain_DS_proportion$Proportion*100


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

# Plateau
Plateau_proportion_plot <- ggplot(Plateau_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                   y = Proportion)) +
  geom_bar(stat = "identity", fill = "#C4D8E9", width = 0.7) +
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

print(Plateau_proportion_plot)

# Plain
Plain_proportion_plot <- ggplot(Plain_DS_proportion[c(1:25),], aes(x = reorder(Defense_System, -Proportion), 
                                                                                      y = Proportion)) +
  geom_bar(stat = "identity", fill = "#d4e09b", width = 0.7) +
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

print(Plain_proportion_plot)

library(Cairo)
library(cowplot)
cairo_pdf(filename = "2DS_proportion.pdf" , width = 3.05,height = 6.86)
g <- plot_grid(Glacier_proportion_plot,Plateau_proportion_plot,
               Plain_proportion_plot,
               
               labels = letters[1:5], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

###2.MAG数量随DS数量的变化##############
#Plateau
Plateau_genome_DSnum <- rowSums(Plateau_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Plateau_genome_DSnum$genome <- rownames(Plateau_genome_DSnum)
Plateau_MAG_metadata <- MAG_metadata[MAG_metadata$Region == "Plateau",]
Plateau_MAG_metadata <- left_join(Plateau_MAG_metadata, 
                             Plateau_genome_DSnum, 
                             by = c("genome" = "genome"))
Plateau_MAG_metadata$DS_num[is.na(Plateau_MAG_metadata$DS_num)] <- 0
Plateau_MAG_DSnum <- select(Plateau_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Plateau_count_data <- Plateau_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Plain
Plain_genome_DSnum <- rowSums(Plain_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Plain_genome_DSnum$genome <- rownames(Plain_genome_DSnum)
Plain_MAG_metadata <- MAG_metadata[MAG_metadata$Region == "Plain",]
Plain_MAG_metadata <- left_join(Plain_MAG_metadata, 
                             Plain_genome_DSnum, 
                             by = c("genome" = "genome"))
Plain_MAG_metadata$DS_num[is.na(Plain_MAG_metadata$DS_num)] <- 0
Plain_MAG_DSnum <- select(Plain_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Plain_count_data <- Plain_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

#Glacier
Glacier_genome_DSnum <- rowSums(Glacier_genome_DS) %>% as.data.frame() %>% setNames("DS_num")
Glacier_genome_DSnum$genome <- rownames(Glacier_genome_DSnum)
Glacier_MAG_metadata <- MAG_metadata[MAG_metadata$Region == "Glacier",]
Glacier_MAG_metadata <- left_join(Glacier_MAG_metadata, 
                             Glacier_genome_DSnum, 
                             by = c("genome" = "genome"))
Glacier_MAG_metadata$DS_num[is.na(Glacier_MAG_metadata$DS_num)] <- 0
Glacier_MAG_DSnum <- select(Glacier_MAG_metadata,genome, DS_num)
# 计算每个DS_num的genome数量
Glacier_count_data <- Glacier_MAG_DSnum %>%
  count(DS_num, name = "genome_count")

library(ggplot2)
#Plateau
Plateau_magnum_dsnum <- ggplot(Plateau_count_data, aes(x = factor(DS_num), y = genome_count, 
                                       fill = "#C4D8E9")) +
  geom_col(width = 0.7, fill = "#C4D8E9", linewidth = 0.3) +
  #geom_text(aes(label = genome_count), vjust = -0.3, size = 3.2, fontface = "bold") +
  #scale_fill_gradient(low = "#E15759", high = "#76B7B2", name = "Count") +
  labs(title = "",
      #subtitle = paste("Total", sum(Plateau_count_data$genome_count), "genomes"),
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
Plateau_magnum_dsnum

#Plain
Plain_magnum_dsnum <- ggplot(Plain_count_data, aes(x = factor(DS_num), y = genome_count, 
                                             fill = "#d4e09b")) +
  geom_col(width = 0.7, fill = "#d4e09b", linewidth = 0.3) +
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
Plain_magnum_dsnum


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
cairo_pdf(filename = "2_MAGnum_DSnum.pdf" , width = 4.33,height = 4.5)
g <- plot_grid(Glacier_magnum_dsnum,Plateau_magnum_dsnum,
               Plain_magnum_dsnum,
               labels = letters[7:11], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###3.数量与MAG大小关系###############

#Plateau
Plateau_MAG_quality <- MAG_quality[MAG_quality$Name %in% Plateau_MAG_DSnum$genome,]

Plateau_MAG_DSdensity <- merge(Plateau_MAG_DSnum, Plateau_MAG_quality, by.x = "genome",by.y = "Name")
Plateau_MAG_DSdensity$MAG_size <- Plateau_MAG_DSdensity$Genome_Size/1000000
Plateau_MAG_DSdensity$DS_density <- Plateau_MAG_DSdensity$DS_num/Plateau_MAG_DSdensity$Genome_Size*1000000

Plateau_MAG_DSdensity$magsize <- floor(Plateau_MAG_DSdensity$MAG_size)

library(dplyr)
Plateau_df_summary <- Plateau_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Plateau_df_summary2 <- Plateau_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    n = n(),
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Plateau_df_summary$magsize, Plateau_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Plateau_df_summary2$magsize, Plateau_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Plateau_df_summary$Region <- 'Plateau'
Plateau_df_summary2$Region <- 'Plateau'

#Plain
Plain_MAG_quality <- MAG_quality[MAG_quality$Name %in% Plain_MAG_DSnum$genome,]

Plain_MAG_DSdensity <- merge(Plain_MAG_DSnum, Plain_MAG_quality, by.x = "genome",by.y = "Name")
Plain_MAG_DSdensity$MAG_size <- Plain_MAG_DSdensity$Genome_Size/1000000
Plain_MAG_DSdensity$DS_density <- Plain_MAG_DSdensity$DS_num/Plain_MAG_DSdensity$Genome_Size*1000000

Plain_MAG_DSdensity$magsize <- floor(Plain_MAG_DSdensity$MAG_size)

library(dplyr)
Plain_df_summary <- Plain_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    mean_value = mean(DS_density),
    sd_value = sd(DS_density)
  )
Plain_df_summary2 <- Plain_MAG_DSdensity %>%
  group_by(magsize) %>%
  summarise(
    n = n(),
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Plain_df_summary$magsize, Plain_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Plain_df_summary2$magsize, Plain_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Plain_df_summary$Region <- 'Plain'
Plain_df_summary2$Region <- 'Plain'

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
    n = n(),
    mean_value = mean(DS_num),
    sd_value = sd(DS_num)
  )

#计算两列数字的相关性及显著性
cor.test(Glacier_df_summary$magsize, Glacier_df_summary$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)
cor.test(Glacier_df_summary2$magsize, Glacier_df_summary2$mean_value,alternative = "two.sided",method = "spearman",conf.level = 0.95)

Glacier_df_summary$Region <- 'Glacier'
Glacier_df_summary2$Region <- 'Glacier'


DSnum_allRegion <- rbind(Glacier_df_summary2,Plateau_df_summary2,
                            Plain_df_summary2)
DSdensity_allRegion <- rbind(Glacier_df_summary,Plateau_df_summary,
                            Plain_df_summary)



library(ggplot2)
DSnum_allRegion_plot <- 
ggplot(DSnum_allRegion, aes(x = magsize, y = mean_value,color = Region)) +
  geom_point(aes(color=Region))+
  geom_line(aes(color=Region))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = Region),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  geom_text(aes(label = n), vjust = -0.8, size = 2.5)+
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
DSnum_allRegion_plot

DSdensity_allRegion_plot <- 
  ggplot(DSdensity_allRegion, aes(x = magsize, y = mean_value,color = Region)) +
  geom_point(aes(color=Region))+
  geom_line(aes(color=Region))+
  geom_errorbar(aes(ymin = mean_value - sd_value,
                    ymax = mean_value + sd_value,
                    color = Region),
                width = 0.1
                
  )+
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
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
DSdensity_allRegion_plot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "2_DS_num-density_allecosystem.pdf" , width = 4.14,height = 5.25)
g <- plot_grid(DSnum_allRegion_plot,DSdensity_allRegion_plot,
               labels = letters[8], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


###4.比较不同环境的Defese family多少，DS per MAG###############

Plateau_DF_num <- sum(Plateau_DS_num$num>0)

Plain_DF_num <- sum(Plain_DS_num$num>0)

Glacier_DF_num <- sum(Glacier_DS_num$num>0)

allRegion_DF_num <- data.frame(Region = c('Glacier','Plateau','Plain'),
                                  DF_num = c(Glacier_DF_num, Plateau_DF_num, Plain_DF_num),
                                  DS_mean = c(sum(Glacier_MAG_metadata$DS_num)/nrow(Glacier_MAG_metadata),
                                              sum(Plateau_MAG_metadata$DS_num)/nrow(Plateau_MAG_metadata),
                                              sum(Plain_MAG_metadata$DS_num)/nrow(Plain_MAG_metadata))
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
p_compare_DStype <- ggplot(allRegion_DF_num, aes(x = reorder(Region,DF_num) ,y = DF_num)) +
  geom_bar(stat = "identity", aes(fill = Region),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  #scale_fill_viridis_c(option = "inferno") + ggtitle("'inferno'")+
  geom_vline(xintercept = 4.5, linetype = "dashed", color = "black")+
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "black")+
  geom_text(aes(label = round(DF_num,2)),vjust = 0.5,hjust = 1.1,size = 4)+
  labs(x = "", y = "Number of DS type", title = "")+  # 设置标签和标题
  mytheme
p_compare_DStype

p_compare_DSmean <- ggplot(allRegion_DF_num, aes(x = reorder(Region,DS_mean) ,y = DS_mean)) +
  geom_bar(stat = "identity", aes(fill = Region),width = 0.5) +  # 使用stat = "identity"因为我们的值已经是计数或比例
  coord_flip() +  # 翻转坐标轴，使条形图水平
  # scale_fill_gradient(low = "#FFF1EB",
  #                     high = "#ACE0F9")+
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
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
cairo_pdf(filename = "2_DStype_DSmean_allRegion.pdf" , width = 2,height = 4.6)
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

Plateau_DS_average <- Plateau_DS_num
Plateau_DS_average$DS_average <- (Plateau_DS_average$num) / nrow(Plateau_MAG_metadata) 
Plateau_DS_average <- left_join(
  Plateau_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Plateau_DS_average[is.na(Plateau_DS_average)] <- "Unknown"

Plain_DS_average <- Plain_DS_num
Plain_DS_average$DS_average <- (Plain_DS_average$num) / nrow(Plain_MAG_metadata) 
Plain_DS_average <- left_join(
  Plain_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Plain_DS_average[is.na(Plain_DS_average)] <- "Unknown"

Glacier_DS_average <- Glacier_DS_num
Glacier_DS_average$DS_average <- (Glacier_DS_average$num) / nrow(Glacier_MAG_metadata) 
Glacier_DS_average <- left_join(
  Glacier_DS_average, 
  DF_mechanism, 
  by = c("DS" = "title")  # 正确指定列名映射
)
Glacier_DS_average[is.na(Glacier_DS_average)] <- "Unknown"


gradient_colors <- c("#FFA5B8", "#8AC7F5", "#FFE4A3", "#9ADFDF","#C4D8E9", "#C0C0C0")

Plateau_treeplot <- ggplot(Plateau_DS_average, 
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
Plateau_treeplot

Plain_treeplot <- ggplot(Plain_DS_average, 
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
Plain_treeplot


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

sum(Plateau_DS_average[Plateau_DS_average$Mechnism=="Unknown",]$DS_average)/sum(Plateau_DS_average$DS_average)
sum(Plain_DS_average[Plain_DS_average$Mechnism=="Unknown",]$DS_average)/sum(Plain_DS_average$DS_average)
sum(Glacier_DS_average[Glacier_DS_average$Mechnism=="Unknown",]$DS_average)/sum(Glacier_DS_average$DS_average)



library(Cairo)
library(cowplot)
cairo_pdf(filename = "2_DS_average_treemap_allRegion.pdf" , width = 6,height = 6)
g <- plot_grid(Glacier_treeplot,Plateau_treeplot,Plain_treeplot,
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
MAG_metadata[,62][is.na(MAG_metadata[,62])] <- 0

MAG_DS_metadata_glacier_plateau_plain <- select(MAG_metadata,genome,ecosystem,Type,Region,DS_num)
MAG_DS_metadata_glacier_plateau_plain$Region <- as.character(MAG_DS_metadata_glacier_plateau_plain$Region) 
MAG_DS_metadata_glacier_plateau_plain$Region <- factor(MAG_DS_metadata_glacier_plateau_plain$Region,
                                                       levels = c("Glacier","Plateau","Plain"))
# MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
#                                             MAG_metadata_plateau_plain$Region)
library(rstatix)
md_p_val_DS <- MAG_DS_metadata_glacier_plateau_plain %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DS_num ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val_DS$y.position <- md_p_val_DS$y.position - 40

# Exact two-sided Mann-Whitney U p-values (normal approximation), Bonferroni-corrected
pw <- pairwise.wilcox.test(
  MAG_DS_metadata_glacier_plateau_plain$DS_num,
  MAG_DS_metadata_glacier_plateau_plain$Region,
  p.adjust.method = "bonferroni", exact = FALSE
)$p.value
# pw is lower-triangular over factor levels: row = later level, col = earlier level
lvl <- levels(MAG_DS_metadata_glacier_plateau_plain$Region)
get_p <- function(g1, g2) {
  hi <- lvl[max(match(g1, lvl), match(g2, lvl))]
  lo <- lvl[min(match(g1, lvl), match(g2, lvl))]
  pw[hi, lo]
}
md_p_val_DS$p.adj <- mapply(get_p, md_p_val_DS$group1, md_p_val_DS$group2)
md_p_val_DS$p.label <- paste0("p = ", format(md_p_val_DS$p.adj, scientific = TRUE, digits = 2))

# n per region (number of MAGs) for x-axis labels
region_n_labels <- MAG_DS_metadata_glacier_plateau_plain %>%
  count(Region) %>%
  mutate(lab = paste0(Region, "\n(n=", n, " MAGs)")) %>%
  {setNames(.$lab, .$Region)}

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)
p_DS_glacier_plateau_plain <- 
  ggplot(MAG_DS_metadata_glacier_plateau_plain, 
         aes(x = Region, y = DS_num)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_half_point(aes(color=Region),
                  side = "r",size = 0.5,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(md_p_val_DS,label = "p.label",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="DS number (per MAG)",
       caption = "Mann-Whitney U test (two-sided), Bonferroni-corrected")+
  scale_x_discrete(labels = region_n_labels)+
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
p_DS_glacier_plateau_plain


library(Cairo)
library(cowplot)
cairo_pdf(filename = "2_DS_average_treemap_allRegion.pdf" , width = 6,height = 6)
g <- plot_grid(Glacier_proportion_plot,
               Plateau_proportion_plot,
               Plain_proportion_plot,
               labels = letters[1:6], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

g1 <- plot_grid(Glacier_magnum_dsnum,
               Plateau_magnum_dsnum,
               Plain_magnum_dsnum,
               labels = letters[1:6], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1


               
g2 <- plot_grid(p_DS_glacier_plateau_plain,
                DSnum_allRegion_plot,DSdensity_allRegion_plot,
               labels = letters[7:9], ncol = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2
cairo_pdf(filename = "2_DS_comparision_MAGsize.pdf" , width = 7.87,height = 3)
g2
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

g3 <- plot_grid(p_compare_DSmean,p_compare_DStype, 
                 ncol = 2
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g3
ggsave("2_two_barplot.pdf",g3, 
       width = 3.93, height = 1.96, units = "in")



g4 <- plot_grid(g1, g3,
                 ncol = 1,
                rel_heights = c(1.5, 0.5)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g4

g5 <- plot_grid(g, g4,
                ncol = 2,
                rel_heights = c(1.5, 0.5)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g5

#7. defense family diversity (per MAG)####
head(DF_defense_2)
head(Glacier_defense)

Glacier_defense_richness <- Glacier_defense %>%
  group_by(genome) %>%
  summarise(
    n_types    = n_distinct(type),
    n_subtypes = n_distinct(subtype),
    .groups    = "drop"
  )
Glacier_defense_richness$Region <- "Glacier"

Plateau_defense_richness <- Plateau_defense %>%
  group_by(genome) %>%
  summarise(
    n_types    = n_distinct(type),
    n_subtypes = n_distinct(subtype),
    .groups    = "drop"
  )
Plateau_defense_richness$Region <- "Plateau"

Plain_defense_richness <- Plain_defense %>%
  group_by(genome) %>%
  summarise(
    n_types    = n_distinct(type),
    n_subtypes = n_distinct(subtype),
    .groups    = "drop"
  )
Plain_defense_richness$Region <- "Plain"

Region_defense_richness <- rbind(Glacier_defense_richness,
                                 Plateau_defense_richness,
                                 Plain_defense_richness)
head(Region_defense_richness)

Region_defense_richness$Region <- 
  factor(Region_defense_richness$Region,
         levels = c("Glacier","Plateau","Plain"))

library(rstatix)
md_p_val_DS_richness <- Region_defense_richness %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(n_subtypes ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val_DS_richness$y.position <-
  md_p_val_DS_richness$y.position - 7



library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
# library(gghalves)
library(ggsci)
# library(ggbeeswarm)
p_DS_richness_glacier_plateau_plain <- 
  ggplot(Region_defense_richness, 
         aes(x = Region, y = n_subtypes)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  # geom_half_point(aes(color=Region),
  #                 side = "r",size = 0.5,
  #                 transformation_params = list(height = 0,
  #                                              width = 0.001,
  #                                              seed = 2))+
  geom_jitter(
    aes(color = Region),
    width = 0.08,
    size = 0.5
  )+
  stat_pvalue_manual(md_p_val_DS_richness,
                     label = "p.adj.signif",
                     label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="Richness diversity (DS subtype per MAG)")+
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
p_DS_richness_glacier_plateau_plain

# 比较不同region的DS number的平均值有无显著差异
# Combine per-MAG DS numbers with region labels
Glacier_MAG_DSnum$Region <- "Glacier"
Plateau_MAG_DSnum$Region <- "Plateau"  
Plain_MAG_DSnum$Region <- "Plain"

# Combine and standardize column names
colnames(Glacier_MAG_DSnum) <- colnames(Plateau_MAG_DSnum) <- colnames(Plain_MAG_DSnum) <- NULL

DS_combined <- rbind(
  data.frame(DS_num = Glacier_MAG_DSnum[,2], Region = "Glacier"),
  data.frame(DS_num = Plateau_MAG_DSnum[,2], Region = "Plateau"),
  data.frame(DS_num = Plain_MAG_DSnum[,2], Region = "Plain")
)
colnames(DS_combined) <- c("DS_num", "Region")
DS_combined$DS_num <- as.numeric(DS_combined$DS_num)

# Summary stats by region
DS_combined |>
  group_by(Region) |>
  summarise(
    n = n(),
    mean = mean(DS_num, na.rm = TRUE),
    sd = sd(DS_num, na.rm = TRUE),
    median = median(DS_num, na.rm = TRUE)
  )

# Check distribution visually
library(ggplot2)
ggplot(DS_combined, aes(x = Region, y = DS_num, fill = Region)) +
  geom_boxplot() +
  labs(y = "Number of Defense Systems", x = NULL) +
  scale_fill_manual(values = c("Glacier" = "#4A90D9", "Plateau" = "#E67E22", "Plain" = "#27AE60")) +
  theme_bw()
# Pairwise Wilcoxon rank-sum tests with Bonferroni correction
pairwise.wilcox.test(DS_combined$DS_num, DS_combined$Region,
                     p.adjust.method = "bonferroni",
                     exact = FALSE)


save.image("./02_region_comparison.RData")

load("./02_region_comparison.RData")

# ####环形柱状图########
# # 筛选比例>1%的系统并计算标签位置
# plot_data <- CJ_DS_proportion[CJ_DS_proportion$Proportion > 1, ] %>%
#   mutate(
#     # 按比例排序
#     Defense_System = factor(Defense_System, 
#                             levels = Defense_System[order(Proportion)]),
#     # 计算每个柱子的中点角度（0-360度）
#     angle = 90 - 360 * (as.numeric(Defense_System) - 0.5) / nrow(.),
# # 调整标签水平对齐方式
# hjust = ifelse(angle < -90, 1, 0),
# # 将角度限制在-180到180度之间
# angle = ifelse(angle < -90, angle + 180, angle),
# # 标签y位置（柱子高度+偏移量）
# label_y = Proportion + max(Proportion)*0.1
# )

# # 创建自定义纵轴
# y_axis_breaks <- seq(0, max(plot_data$Proportion), by = 10)
# y_axis_labels <- paste0(y_axis_breaks, "%")

# # 创建环形条形图
# circular_plot <- ggplot(plot_data, 
#                         aes(x = Defense_System, 
#                             y = Proportion, 
#                             fill = Proportion)) +
#   geom_col(width = 0.7, color = "white", linewidth = 0.2) +  # 使用柱状图
  
#   annotate(
#     "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0,
#     fill = "white", color = 'black',linewidth = 0.1
#   ) +
  
  
#   scale_fill_gradient(low = "#BEBCDF", high = "#BEBCDF") +   # 颜色渐变
  
#   # 极坐标转换
#   coord_polar(start = -pi/nrow(plot_data)) +  # 调整起始角度使第一个柱子位于顶部
  
#   # 设置y轴范围（中心留白）
#   ylim(-max(plot_data$Proportion)*0.3, max(plot_data$Proportion)*1.2) +
  
#   # 添加标签（自动调整角度）
#   geom_text(aes(y = label_y, label = Defense_System, angle = angle, hjust = hjust),
#             size = 2, color = "black", fontface = "plain") +
  
#   # 标题和主题
#   labs(title = "") +
#   theme_void() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "plain", size = 12, margin = margin(b = 10)),
#     legend.position = "none",
#     plot.margin = unit(c(1, 1, 1, 1), "cm")
#   )

# print(circular_plot)
# # 创建基础图形

 