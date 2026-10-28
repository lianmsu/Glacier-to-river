# =============================================================================
# Glacier-to-river -- defensome landscape / Defense-system composition (PCoA / beta diversity)
# =============================================================================
# Ordination of MAGs by defense-system composition, alpha/beta diversity of
# the defensome and the DS-frequency heatmap.
#
# Manuscript      : Fig. 2e-f
# Main outputs    : 5_heatmap_DSfrequency.pdf, 05_defense_system_composition_pcoa.RData
# Original script : 20gaoyuan/10Defensefinder/5_MAG_DS_pcoa.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

#基本数据读取####
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

DF_mechanism <- read.delim2('df_List_Systems_lpw.csv',sep = ',')
DF_mechanism$Effector[DF_mechanism$Effector == ""] <- "Unknown"
DF_mechanism$Mechnism[DF_mechanism$Mechnism == ""] <- "Unknown"

#每个点代表一种DS####
##Bray NMDS###############

library(tidyverse)
library(vegan)
library(ape)
#library(phyloseq)
library(geomtextpath)
library(patchwork)


# 标准化数据，使其适用于 Bray-Curtis 距离计算
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


count_matrix_standard = decostand(count_matrix_standard,method = "hellinger")  
# NMDS（非度量多维标度分析），k=2 表示降维到 2 维
mds_bray = metaMDS(t(count_matrix_standard),
                   distance = "bray", k = 2,trymax = 100)  


# 提取 NMDS 分析的 stress 值，四舍五入保留 4 位小数
# 计算 stress 值，越小表示 NMDS 结果越可靠
stress1 <- round(mds_bray$stress,4) 

# 提取 NMDS 计算出的坐标，并与样本分组信息合并
data.plot <- mds_bray$points %>% as.data.frame() %>%  # 转换为数据框
  rownames_to_column(var="title") %>%  # 添加样本 ID 列
  left_join(., DF_mechanism, by="title")  # 关联样本分组信息
data.plot[is.na(data.plot)] <- "Unknown"


# 绘制 Bray-Curtis NMDS 结果
p1 <- ggplot(data.plot, aes(MDS1, MDS2))+
  # 按 Mechnism 颜色区分样本点
  geom_point( aes(colour = Mechnism),size = 2, show.legend = F)+  
  # 95% 置信椭圆
  stat_ellipse(aes(label = Mechnism, color=Mechnism),
               geom = "textpath",
               hjust = 0.3, level = 0.95, linewidth=0.5,
               linetype=2, show.legend = F) +
  xlab("MDS1") + ylab("MDS2") +  # 轴标签
  # 添加 stress 值到标题
  labs(title = paste("bray stress =", stress1, sep=" ")) +  
  theme_test() +
  theme(plot.title = element_text(vjust=0.5, hjust = 0.5))  
p1


##PCoA######
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]

library(vegan)
library(ape)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

# 计算 Bray-Curtis 距离并进行 PCoA 分析
pcoa <- vegdist(t(count_matrix_standard), method = "bray") %>% 
    # 进行主坐标分析（PCoA），不进行负特征值校正
  pcoa(correction = "none")

summary(pcoa)

# 提取 PCoA 结果并合并分组信息
pcoadata <- data.frame(pcoa$vectors[,1], pcoa$vectors[,2]) %>%  
  # 获取 PC1 和 PC2 坐标
  set_colnames(c("PC1", "PC2")) %>%  # 重命名列名
  rownames_to_column(var = "title") %>%  # 将行名转换为列 sample
  # 关联 Subtype 信息
  left_join(., DF_mechanism %>% 
              dplyr::select(title, Mechnism),by="title")

pcoadata[is.na(pcoadata)] <- "Unknown"
# 颜色设置
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(6)

#绘制pcoa散点图
plot <- ggplot(pcoadata, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Mechnism, fill = Mechnism), size = 4) +  
  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot

################################################################

#携带DS的MAG的alpha diversity变化####
subtype_count_matrix <- DF_defense %>%
  count(genome, subtype) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, subtype, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = subtype,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

# 设置行名
rownames(subtype_count_matrix) <- subtype_count_matrix$genome
subtype_count_matrix <- subtype_count_matrix[-1]  # 移除genome列
subtype_count_matrix <-
  subtype_count_matrix[rowSums(subtype_count_matrix) >0,]

subtype_count_matrix <- 
  subtype_count_matrix[rownames(subtype_count_matrix) %in% 
                         MAG_metadata$genome,]


##alpha diversity , 用DS subtype####

library(tidyverse)
library(vegan)

# 7. 转置矩阵：vegan期望行为样本，列为物种
r_MAG_c_DSsubtype <- subtype_count_matrix

# 方法1：使用vegan包直接计算
calculate_alpha_diversity <- function(abundance_matrix) {
  # 注意：abundance_matrix应该是行为样本，列为物种
  
  # 1. 转换为相对丰度（确保每行和为100%）
  # 你的数据已经是百分比，需要转换为小数（0-1）
  abundance_decimal <- abundance_matrix / 100
  
  # 2. 计算各种α多样性指数
  alpha_results <- data.frame(
    Sample = rownames(abundance_decimal),
    # 丰富度：观测到的基因组数（非零值）
    Observed_Richness = apply(abundance_decimal, 1, function(x) sum(x > 0)),
    
    # Shannon指数
    Shannon = diversity(abundance_decimal, index = "shannon"),
    
    # Simpson指数
    Simpson = diversity(abundance_decimal, index = "simpson"),
    
    # Inverse Simpson指数
    InvSimpson = diversity(abundance_decimal, index = "invsimpson"),
    
    # Pielou均匀度
    Pielou = apply(abundance_decimal, 1, function(x) {
      p <- x[x > 0]
      S <- length(p)
      if(S <= 1) return(1)
      H <- -sum(p * log(p))
      H / log(S)
    })
  )
  
  # 添加有效物种数
  alpha_results$Effective_Species_Shannon <- exp(alpha_results$Shannon)
  alpha_results$Effective_Species_Simpson <- alpha_results$InvSimpson
  
  return(alpha_results)
}

# 计算α多样性
DS_MAG_alpha_diversity <- 
  calculate_alpha_diversity(r_MAG_c_DSsubtype)

print("α多样性结果:")
head(DS_MAG_alpha_diversity)

# 排序查看
DS_MAG_alpha_diversity_sorted <- DS_MAG_alpha_diversity %>%
  arrange(desc(Shannon))

# print("\n按Shannon多样性排序:")
# print(DS_MAG_alpha_diversity_sorted)

DS_MAG_alpha_diversity_sorted <- 
  DS_MAG_alpha_diversity_sorted %>% 
  left_join(., MAG_metadata, by = c("Sample" = "genome"))

DS_MAG_alpha_diversity_sorted$Region <- 
  factor(DS_MAG_alpha_diversity_sorted$Region,
         levels = c("Glacier","Plateau","Plain"))

# # 从样本名提取深度和位置信息
# alpha_diversity <- alpha_diversity %>%
#   mutate(
#     Depth = as.numeric(str_extract(Sample, "\\d+")),
#     Location = str_extract(Sample, "[A-Z]+$"),
#     Depth_Group = case_when(
#       Depth <= 1000 ~ "浅层(≤1000m)",
#       Depth <= 3000 ~ "中层(1000-3000m)",
#       TRUE ~ "深层(>3000m)"
#     ),
#     Location_Type = case_when(
#       Location %in% c("ONE", "TWO", "THREE") ~ "核心站位",
#       Location %in% c("FAR", "NEAR") ~ "边缘站位",
#       TRUE ~ "其他"
#     )
#   )

# 分组箱线图
p_DS_MAG_pielou_region <- 
  ggplot(DS_MAG_alpha_diversity_sorted, aes(x = Region, y = Pielou)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  # geom_text(aes(label = Sample), 
  #           data = alpha_diversity_sorted %>% filter(Pielou > 0.6),
  #           vjust = -1.5, size = 3, color = "darkred") +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  #  scale_fill_brewer(palette = "Blues", name = "Region (m)") +
  #  scale_color_brewer(palette = "Set2", name = "位置") +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    #axis.title = element_blank(),
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    #axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

p_DS_MAG_pielou_region

p_DS_MAG_shannon_region <- 
  ggplot(DS_MAG_alpha_diversity_sorted, aes(x = Region, y = Shannon)) +
#  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
      geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  # geom_text(aes(label = Sample), 
  #           data = alpha_diversity %>% filter(Shannon > 0.6),
  #           vjust = -1.5, size = 3, color = "darkred") +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  #  scale_fill_brewer(palette = "Blues", name = "Depth (m)") +
  #  scale_color_brewer(palette = "Set2", name = "位置") +
  labs(x="",y="Shannon of DS (per MAG)")+
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
  )

p_DS_MAG_shannon_region

p_DS_MAG_simpson_region <- 
  ggplot(DS_MAG_alpha_diversity_sorted, aes(x = Region, y = Simpson)) +
#  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
      geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  # geom_text(aes(label = Sample), 
  #           data = alpha_diversity %>% filter(Shannon > 0.6),
  #           vjust = -1.5, size = 3, color = "darkred") +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  #  scale_fill_brewer(palette = "Blues", name = "Depth (m)") +
  #  scale_color_brewer(palette = "Set2", name = "位置") +
  labs(x="",y="Simpson of DS (per MAG)")+
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
  )

p_DS_MAG_simpson_region

p_DS_MAG_richness_region <- 
  ggplot(DS_MAG_alpha_diversity_sorted, aes(x = Region, y = Observed_Richness)) +
#  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
    geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  # geom_text(aes(label = Sample), 
  #           data = alpha_diversity %>% filter(Shannon > 0.6),
  #           vjust = -1.5, size = 3, color = "darkred") +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  #  scale_fill_brewer(palette = "Blues", name = "Depth (m)") +
  #  scale_color_brewer(palette = "Set2", name = "位置") +
  labs(x="",y="Richness of DS (per MAG)")+
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
  )

p_DS_MAG_richness_region

library(Cairo)
library(cowplot)

g5_DS_MAG_richness_region <- 
  plot_grid(p_DS_MAG_shannon_region,
            p_DS_MAG_simpson_region,p_DS_MAG_richness_region,
            labels = letters[1:3], ncol = 3
            # align = "v",  # 垂直对齐
            # axis = "lr",  # 左和右对齐
            # rel_widths = c(1, 1,1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g5_DS_MAG_richness_region

ggsave(plot = g5_DS_MAG_richness_region,
       '5_g5_DS_MAG_richness_region_new.pdf',width = 8,height = 2.6,units = 'in')

#统计检验#

library(rstatix)

# Shannon, Simpson, Observed（请根据实际列名修改）
diversity_indices <- c("Shannon", "Simpson", "Observed_Richness")  
multi_test_results <- DS_MAG_alpha_diversity_sorted %>%
  # 1. 将几个多样性列“揉”成一个长表
  pivot_longer(
    cols = all_of(diversity_indices),
    names_to = "index",
    values_to = "value"
  ) %>%
  # 2. 按多样性指数分组，各自做 Wilcoxon 检验
  group_by(index) %>%
  wilcox_test(value ~ Region)

# 3. Bonferroni 按指数内 3 次两两比较校正。
#    adjust_pvalue 不认 group_by 会跨 9 行校正，故用 p.adjust 组内校正；
#    mutate 会丢失 rstatix 的 args 属性，先保存再恢复以保留 add_xy_position。
args_saved <- attr(multi_test_results, "args")
multi_test_results <- multi_test_results %>%
  group_by(index) %>%
  mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
  ungroup()
attr(multi_test_results, "args") <- args_saved

# 4. 显著性标注与绘图位置
multi_test_results <- multi_test_results %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Region", dodge = 0.5)
# md_p_val_DS_MAG_richness <- DS_MAG_alpha_diversity_sorted %>% 
#   #  group_by(ecosystem) %>%
#   wilcox_test(n_subtypes ~ Region) %>%
#   adjust_pvalue(p.col="p",method="bonferroni") %>%
#   add_significance(p.col="p.adj") %>% 
#   add_xy_position(x="Region",dodge=0.5)



save.image("./05_defense_system_composition_pcoa.RData")

write.table(DS_MAG_alpha_diversity_sorted,
    "DS_MAG_alpha_diversity_sorted.tsv",
    sep = '\t',quote = F, row.names = F)


#每个点代表一个MAGs###########
##PCoA######
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]

library(vegan)
library(ape)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

# 计算 Bray-Curtis 距离并进行 PCoA 分析
pcoa <- vegdist(count_matrix_standard, method = "bray") %>% 
  # 进行主坐标分析（PCoA），不进行负特征值校正
  pcoa(correction = "none")

summary(pcoa)

# 提取 PCoA 结果并合并分组信息
pcoadata <- data.frame(pcoa$vectors[,1], pcoa$vectors[,2]) %>%  
  # 获取 PC1 和 PC2 坐标
  set_colnames(c("PC1", "PC2")) %>%  # 重命名列名
  rownames_to_column(var = "genome") %>%  # 将行名转换为列 sample
  # 关联 Subtype 信息
  left_join(., MAG_metadata %>% 
              dplyr::select(genome, Region),by="genome")

pcoadata[is.na(pcoadata)] <- "Unknown"

#write.csv(pcoadata,"pcoadata.csv")

pcoadata2 <- read.csv("pcoadata.csv")

mag_baclineage <- read.delim2("gtdbtk.bac120.summary.tsv")
mag_arclineage <- read.delim2("gtdbtk.ar53.summary.tsv")
mag_lineage <- rbind(mag_baclineage, mag_arclineage) %>% select(user_genome, classification)
library(tidyr)
# 拆分 classification 列
mag_lineage_split <- mag_lineage %>%
  separate(
    classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )
# 查看拆分后的数据
head(mag_lineage_split)
colnames(mag_lineage_split)[1] <- 'genome'

pcoadata2 <- read.csv("pcoadata.csv")
pcoadata2 <- left_join(pcoadata2,mag_lineage_split,suffix = c("genome","genome"))
# 统计 phylum 列中每个门的计数
p_counts <- pcoadata2 %>%
  count(phylum)
# 找出计数小于 1000 的门
low_count_p <- p_counts %>%
  filter(n < 50) %>%
  pull(phylum)
# 将计数小于 1000 的门及其后续分类等级设置为 "others"
pcoadata2 <- pcoadata2 %>%
  mutate(across(phylum:species, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                         ifelse(phylum %in% low_count_p, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))
phylum <- table(pcoadata2$phylum) %>% as.data.frame()

#添加每个MAG的size大小
pcoadata3 <- pcoadata2 %>% left_join(., MAG_metadata, by="genome")
colnames(MAG_quality)[1] <- "genome"
pcoadata3 <- pcoadata3 %>% left_join(., MAG_quality, by="genome")
pcoadata3$MAG_size <- pcoadata3$Genome_Size/1000000
pcoadata3$magsize <- floor(pcoadata3$MAG_size)
pcoadata3$magsize <- factor(pcoadata3$magsize)

###1、根据ecosystem绘图####

# 颜色设置
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(6)

library(ggplot2)
#绘制pcoa散点图
plot_MAG_ecosystem <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = ecosystem, fill = ecosystem), size = 0.5) +  
#  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
       #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
       legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_ecosystem

###2、region着色####


#绘制pcoa散点图
plot_MAG_region <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Region.x, fill = Region.x), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_region

###3、medium着色####
#绘制pcoa散点图
plot_MAG_medium <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = medium, fill = medium), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_medium

###4、物种着色####
phylum <- table(pcoadata2$phylum) %>% as.data.frame()
color_palette_phylum <- colorRampPalette(brewer.pal(12, "Paired"))(10)


#绘制pcoa散点图
plot_MAG_phylum <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = phylum, fill = phylum), size = 0.5) +  
  scale_color_manual(values = color_palette_phylum) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_phylum

###5、用MAG size着色####


color_palette_magsize <- colorRampPalette(brewer.pal(12, "Paired"))(12)
#绘制pcoa散点图
plot_MAG_DS_magsize <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = magsize, fill = magsize), size = 0.5) +  
  scale_color_manual(values = color_palette_magsize) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_magsize


###6、用season着色####


color_palette_season <- colorRampPalette(brewer.pal(12, "Paired"))(4)
#绘制pcoa散点图
plot_MAG_DS_season <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Season, fill = Season), size = 0.5) +  
  scale_color_manual(values = color_palette_season) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_season

###7、用domain着色####

pcoadata3$domain <- factor(pcoadata3$domain, levels = c("d__Archaea","d__Bacteria"))
#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_domain <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = domain, fill = domain), size = 0.5) +  
#  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
                                "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_domain

###8、用Elevation着色####

#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_Elevation <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Elevation, fill = Elevation), size = 0.5) +  
  #  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
  #                               "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
        color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_Elevation


library(Cairo)
library(cowplot)
cairo_pdf(filename = "5_DSpcoa.pdf" , width = 8,height = 8)
g <- plot_grid(plot_MAG_ecosystem,plot_MAG_region,plot_MAG_medium,
               plot_MAG_phylum,plot_MAG_DS_magsize,
               labels = letters[5], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

write.table(pcoadata3, "pcoadata3.tsv", sep = "\t", quote = F, row.names = F)




##NMDS####
# 标准化数据，使其适用于 Bray-Curtis 距离计算
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


count_matrix_standard = decostand(count_matrix_standard,method = "hellinger")  
# NMDS（非度量多维标度分析），k=2 表示降维到 2 维
mds_bray = metaMDS(count_matrix_standard,
                   distance = "bray", k = 2,trymax = 100)  


# 提取 NMDS 分析的 stress 值，四舍五入保留 4 位小数
# 计算 stress 值，越小表示 NMDS 结果越可靠
stress1 <- round(mds_bray$stress,4) 

# 提取 NMDS 计算出的坐标，并与样本分组信息合并
data.plot <- mds_bray$points %>% as.data.frame() %>%  # 转换为数据框
  rownames_to_column(var="genome") %>%  # 添加样本 ID 列
  left_join(., MAG_metadata, by="genome")  # 关联样本分组信息
data.plot[is.na(data.plot)] <- "Unknown"


# 绘制 Bray-Curtis NMDS 结果
p_nmds_ecosystem <- ggplot(data.plot, aes(MDS1, MDS2))+
  # 按 ecosystem 颜色区分样本点
  geom_point( aes(colour = ecosystem),size = 2, show.legend = F)+  
  # 95% 置信椭圆
  stat_ellipse(aes(label = ecosystem, color=ecosystem),
               geom = "textpath",
               hjust = 0.3, level = 0.95, linewidth=0.5,
               linetype=2, show.legend = F) +
  xlab("MDS1") + ylab("MDS2") +  # 轴标签
  # 添加 stress 值到标题
  labs(title = paste("bray stress =", stress1, sep=" ")) +  
  theme_test() +
  theme(plot.title = element_text(vjust=0.5, hjust = 0.5))  
p_nmds_ecosystem


##PCA####
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


pca <- count_matrix_standard %>% 
  prcomp(., scale. = TRUE)  # 进行标准化 PCA 计算
# 计算 PCA 解释的方差比例
var_explained <- pca$sdev^2 / sum(pca$sdev^2)

pca


DF_mechanism_pca <- data.frame(title = colnames(count_matrix_standard))
DF_mechanism_pca <- left_join(DF_mechanism_pca, DF_mechanism, by="title")
DF_mechanism_pca[is.na(DF_mechanism_pca)] <- "Unknown"
DF_mechanism_pca <- select(DF_mechanism_pca,title,Mechnism)


library(factoextra)
# 绘制 PCA 双标图（biplot）
p_pca <- fviz_pca_biplot(pca, axes = c(1, 2),  # 选择主成分 PC1 和 PC2
                      geom.ind = c("point"),  # 样本点
                      geom.var = c("arrow", "text"),  # 变量箭头 + 标签
                      pointshape = 20,  # 设定点的形状
                      pointsize = 4,  # 设定点的大小
                      # 变量标签避免重叠，颜色设为黑色
                      label = "var", repel = TRUE, col.var = "black",  
                      labelsize = 0.5,  # 变量标签的大小
                      addEllipses = TRUE,  # 添加 95% 置信椭圆
                      ellipse.level = 0.95,  # 设置置信区间水平
                      col.ind = DF_mechanism_pca$Mechnism  # 按 `Subtype` 分类颜色
) +
  scale_color_manual(values =
                       colorRampPalette(brewer.pal(12,"Paired"))(6)) +  
  # 设置 X 轴标签
  labs(x = paste0("(PC1: ",round(var_explained[1] * 100,2),"%)"),  
       y = paste0("(PC2: ",round(var_explained[2] * 100,2),"%)")) + 
  # 设置 X 轴标题样式
  theme(panel.background=element_rect(fill='white',colour='black'),  
        axis.title.x = element_text(
          colour = "black",size = 12,
          margin = margin(t = 12)), 
        # 设置 Y 轴标题样式
        axis.title.y = element_text(
          colour = "black", size = 12, margin = margin(r = 12)), 
        axis.text = element_text(color = "black"),  # 设置轴刻度文字颜色
        plot.title = element_blank(),  # 不显示标题
        legend.title = element_blank(),  # 不显示图例标题
        legend.key = element_blank(),  # 使图例背景透明
        legend.text = element_text(color = "black", size = 9),  
        # 设置图例文本之间的水平间距
        legend.spacing.x = unit(0.06, 'cm'),  
        legend.key.width = unit(0.01, 'cm'),  # 设置图例的水平大小
        legend.key.height = unit(0.01, 'cm'),  # 设置图例的垂直大小
        legend.background = element_blank(),  # 使图例背景透明
        legend.position = c(1, 0),  # 设置图例在图的右下角
        legend.justification = c(1, 0))  # 使图例对齐到右下角
p_pca


#根据pcoa结果绘制不同cluster的DS组成特点####

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))

pcoadata2 <- read.csv("pcoadata.csv")
pcoadata2 <- pcoadata2[pcoadata2$genome %in% MAG_metadata$genome,]

pcoadata2$type <- ifelse(pcoadata2$PC1 < 0.125, "cluster12", "cluster34")
pcoadata2_cluster12 <- pcoadata2[pcoadata2$type == "cluster12",]

pcoadata2_cluster2_part1 <- pcoadata2_cluster12[pcoadata2_cluster12$PC1>=0 & pcoadata2_cluster12$PC2 <= 0.0625,]
pcoadata2_cluster2_part2 <- pcoadata2_cluster12[pcoadata2_cluster12$PC1>=-0.1 &pcoadata2_cluster12$PC1<0 & pcoadata2_cluster12$PC2 <= 0,]
pcoadata2_cluster2_part3 <- pcoadata2_cluster12[pcoadata2_cluster12$PC1>=-0.25 & pcoadata2_cluster12$PC1< -0.1 & pcoadata2_cluster12$PC2 <= -0.125,]
pcoadata2_cluster2 <- rbind(pcoadata2_cluster2_part1, pcoadata2_cluster2_part2, pcoadata2_cluster2_part3)
pcoadata2_cluster2$type <- "Cluster2"

pcoadata2_cluster1 <- pcoadata2_cluster12[!pcoadata2_cluster12$genome 
                                          %in% pcoadata2_cluster2$genome,]
pcoadata2_cluster1$type <- "Cluster1"

pcoadata2_cluster34 <- pcoadata2[pcoadata2$type == "cluster34",]
pcoadata2_cluster3 <- pcoadata2_cluster34[pcoadata2_cluster34$PC2>0.125,]
pcoadata2_cluster4 <- pcoadata2_cluster34[pcoadata2_cluster34$PC2<0.125,]
pcoadata2_cluster3$type <- "Cluster3"
pcoadata2_cluster4$type <- "Cluster4"

pcoadata2_cluster <- rbind(pcoadata2_cluster1,pcoadata2_cluster2,pcoadata2_cluster3,pcoadata2_cluster4)

dss <- count_matrix
dss$genome <- rownames(count_matrix)
dss <- dss[dss$genome %in% pcoadata2_cluster$genome,]
dss <- merge(dss, pcoadata2_cluster, by = "genome")
dss <- select(dss, genome, type, everything())
dss <- select(dss, -PC1, -PC2, -Region, -X)
#不同DS的均值/出现率

#均值
library(dplyr)
dss_average_cluster <- dss %>%
  select(-genome) %>%               # 移除 genome 列
  group_by(type) %>%                # 按 type 分组
  summarise(across(everything(),    # 对所有列计算均值
            ~mean(., na.rm = TRUE))) %>% # 使用匿名函数格式
  data.frame()
rownames(dss_average_cluster) <- dss_average_cluster$type
dss_average_cluster <- dss_average_cluster[,-1]

#绘制热图
library(pheatmap)

colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p_heatmap <- pheatmap(t(dss_average_cluster),
                #color = colorRampPalette(colors)(1000),
                #color = colorRampPalette(c('white','lightpink','red'))(1000),
                color = colors,
                cluster_cols = T,
                cluster_rows = F,
                #annotation_col = super_phyla[,-1],
                fontsize_row = 10
)
p_heatmap
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))

#出现率
dss_frequency_cluster <- dss %>%
  select(-genome) %>%               # 移除 genome 列
  group_by(type) %>%                # 按 type 分组
  summarise(across(everything(),    # 对所有列计算
            ~ mean(. != 0) * 100)) %>%   # 非零值比例转换为百分比
  data.frame(check.names = F)
rownames(dss_frequency_cluster) <- dss_frequency_cluster$type
dss_frequency_cluster <- dss_frequency_cluster[,-1]

#绘制热图
library(pheatmap)

colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p_heatmap <- pheatmap(t(dss_frequency_cluster),
                      #color = colorRampPalette(colors)(1000),
                      #color = colorRampPalette(c('white','lightpink','red'))(1000),
                      color = colors,
                      cluster_cols = T,
                      cluster_rows = F,
                      #annotation_col = super_phyla[,-1],
                      fontsize_row = 10
)
p_heatmap
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))

dss_frequency_cluster_filter <- dss_frequency_cluster[,colSums(dss_frequency_cluster)>3]

p_heatmap_65 <- pheatmap(t(dss_frequency_cluster_filter),
                      #color = colorRampPalette(colors)(1000),
                      #color = colorRampPalette(c('white','lightpink','red'))(1000),
                      color = colors,
                      cluster_cols = T,
                      cluster_rows = F,
                      #annotation_col = super_phyla[,-1],
                      fontsize_row = 10
)
p_heatmap_65
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
ggsave("5_heatmap_DSfrequency.pdf",p_heatmap_65, width = 3, height = 8.5, units = "in")

write.table(dss_frequency_cluster_filter, "dss_frequency_cluster_filter.tsv",
  sep = "\t", quote = F)

#删除RM和SoFIC再次进行PCoA####
##PCoA######
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]

count_matrix_standard <- count_matrix_standard[,!colnames(count_matrix_standard) %in% c("RM","SoFIC")]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


library(vegan)
library(ape)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

# 计算 Bray-Curtis 距离并进行 PCoA 分析
pcoa <- vegdist(count_matrix_standard, method = "bray") %>% 
  # 进行主坐标分析（PCoA），不进行负特征值校正
  pcoa(correction = "none")

summary(pcoa)

# 提取 PCoA 结果并合并分组信息
pcoadata <- data.frame(pcoa$vectors[,1], pcoa$vectors[,2]) %>%  
  # 获取 PC1 和 PC2 坐标
  set_colnames(c("PC1", "PC2")) %>%  # 重命名列名
  rownames_to_column(var = "genome") %>%  # 将行名转换为列 sample
  # 关联 Subtype 信息
  left_join(., MAG_metadata %>% 
              dplyr::select(genome, Region),by="genome")

pcoadata[is.na(pcoadata)] <- "Unknown"

#write.csv(pcoadata,"pcoadata_noRMSoFIC.csv")
pcoadata2 <- read.csv("pcoadata_noRMSoFIC.csv")

mag_baclineage <- read.delim2("gtdbtk.bac120.summary.tsv")
mag_arclineage <- read.delim2("gtdbtk.ar53.summary.tsv")
mag_lineage <- rbind(mag_baclineage, mag_arclineage) %>% select(user_genome, classification)
library(tidyr)
# 拆分 classification 列
mag_lineage_split <- mag_lineage %>%
  separate(
    classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )
# 查看拆分后的数据
head(mag_lineage_split)
colnames(mag_lineage_split)[1] <- 'genome'

pcoadata2 <- read.csv("pcoadata_noRMSoFIC.csv")
pcoadata2 <- left_join(pcoadata2,mag_lineage_split,suffix = c("genome","genome"))
# 统计 phylum 列中每个门的计数
p_counts <- pcoadata2 %>%
  count(phylum)
# 找出计数小于 1000 的门
low_count_p <- p_counts %>%
  filter(n < 50) %>%
  pull(phylum)
# 将计数小于 1000 的门及其后续分类等级设置为 "others"
pcoadata2 <- pcoadata2 %>%
  mutate(across(phylum:species, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                         ifelse(phylum %in% low_count_p, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))
phylum <- table(pcoadata2$phylum) %>% as.data.frame()

#添加每个MAG的size大小
pcoadata3 <- pcoadata2 %>% left_join(., MAG_metadata, by="genome")
colnames(MAG_quality)[1] <- "genome"
pcoadata3 <- pcoadata3 %>% left_join(., MAG_quality, by="genome")
pcoadata3$MAG_size <- pcoadata3$Genome_Size/1000000
pcoadata3$magsize <- floor(pcoadata3$MAG_size)
pcoadata3$magsize <- factor(pcoadata3$magsize)

###1、根据ecosystem绘图####

# 颜色设置
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(6)

library(ggplot2)
#绘制pcoa散点图
plot_MAG_ecosystem <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = ecosystem, fill = ecosystem), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_ecosystem

###2、region着色####


#绘制pcoa散点图
plot_MAG_region <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Region.x, fill = Region.x), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_region

###3、medium着色####
#绘制pcoa散点图
plot_MAG_medium <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = medium, fill = medium), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_medium

###4、物种着色####
phylum <- table(pcoadata2$phylum) %>% as.data.frame()
color_palette_phylum <- colorRampPalette(brewer.pal(12, "Paired"))(10)


#绘制pcoa散点图
plot_MAG_phylum <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = phylum, fill = phylum), size = 0.5) +  
  scale_color_manual(values = color_palette_phylum) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_phylum

###5、用MAG size着色####


color_palette_magsize <- colorRampPalette(brewer.pal(12, "Paired"))(12)
#绘制pcoa散点图
plot_MAG_DS_magsize <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = magsize, fill = magsize), size = 0.5) +  
  scale_color_manual(values = color_palette_magsize) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_magsize


###6、用season着色####


color_palette_season <- colorRampPalette(brewer.pal(12, "Paired"))(4)
#绘制pcoa散点图
plot_MAG_DS_season <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Season, fill = Season), size = 0.5) +  
  scale_color_manual(values = color_palette_season) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_season

###7、用domain着色####

pcoadata3$domain <- factor(pcoadata3$domain, levels = c("d__Archaea","d__Bacteria"))
#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_domain <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = domain, fill = domain), size = 0.5) +  
  #  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
                                "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_domain

###8、用Elevation着色####

#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_Elevation <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Elevation, fill = Elevation), size = 0.5) +  
  #  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
  #                               "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_Elevation


library(Cairo)
library(cowplot)
cairo_pdf(filename = "5_DSpcoa_noRMSoFIC.pdf", width = 8,height = 8)
g <- plot_grid(plot_MAG_ecosystem,plot_MAG_region,plot_MAG_medium,
               plot_MAG_phylum,plot_MAG_DS_magsize,
               labels = letters[5], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

write.table(pcoadata3, "pcoadata3_2.tsv", sep = "\t", quote = F, row.names = F)

save.image("./05_defense_system_composition_pcoa.RData")
load("./05_defense_system_composition_pcoa.RData")





