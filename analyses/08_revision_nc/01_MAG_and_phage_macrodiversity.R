# =============================================================================
# Glacier-to-river -- revision nc / MAG and phage macro-diversity
# =============================================================================
# Alpha/beta diversity of MAGs and phages along the glacier-to-river
# gradient (reviewer revision).
#
# Manuscript      : Revision
# Main outputs    : 01_MAG_and_phage_macrodiversity.RData,
#                   01_MAG_and_phage_macrodiversity_new.RData, NMDS/PCoA tables
# Original script : 93_NC_comments/01_diversity_comparison/1_MAG_Phage_macrodiversity.R
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

setwd(file.path(GTR_DATA, "93_NC_comments/01_diversity_comparison"))
library(dplyr)
library(stringr)
library(ggplot2)

#MAG diversity ####

##alpha diversity####

# 相对丰度的本质：组成性数据
# 每个样本的所有物种相对丰度之和 = 100%

# α多样性计算的假设：
# 1. 数据是组成性的（总和为常数）
# 2. 物种间可比较
# 相对丰度完美满足这些条件

MAGs_4903_relative_abundance <- 
  read.delim2("./丰度/Bins/all_abundance_rpkm_bins.tsv")

# mag_data_rel_only <- MAGs_544_relative_abundance %>%
#   select(Genome, contains("Relative.Abundance"))
# 将所有非Genome列转换为数值型
MAGs_4903_relative_abundance <- MAGs_4903_relative_abundance %>%
  mutate(across(-Genome, as.numeric))

library(tidyverse)
library(vegan)

# 1. 提取相对丰度矩阵
# 假设所有以"Relative.Abundance"结尾的列都是样本的相对丰度
rel_abundance_cols <- 
  grep("sorted\\.RPKM", colnames(MAGs_4903_relative_abundance), value = TRUE)

# 2. 创建样本名（简化列名）
sample_names <- 
  gsub("\\.sorted\\.RPKM\\.\\.\\.\\.", "", rel_abundance_cols)
sample_names <- gsub(".sorted.RPKM", "", sample_names)  # 可选：进一步简化
#sample_names <- gsub("_", " ", sample_names)      # 可选：用空格替代下划线

# 3. 提取相对丰度数据
rel_abundance_matrix <- MAGs_4903_relative_abundance[, rel_abundance_cols]

# 4. 设置行名和列名
rownames(rel_abundance_matrix) <- MAGs_4903_relative_abundance$Genome
colnames(rel_abundance_matrix) <- sample_names

# 5. 查看数据结构
cat("数据维度:", dim(rel_abundance_matrix), "\n")
cat("样本数:", ncol(rel_abundance_matrix), "\n")
cat("基因组数:", nrow(rel_abundance_matrix), "\n")
cat("样本名称:", colnames(rel_abundance_matrix), "\n")

# 6. 移除"unmapped"行（如果需要）
rel_abundance_matrix <- 
  rel_abundance_matrix[rownames(rel_abundance_matrix) != "unmapped", ]

# 7. 转置矩阵：vegan期望行为样本，列为物种
rel_abundance_t <- t(rel_abundance_matrix)

# 方法1：使用vegan包直接计算
calculate_alpha_diversity <- function(abundance_matrix) {
  # 注意：abundance_matrix 应为行为样本、列为物种
  # vegan 内部会将每行标准化为相对丰度（和为1），因此无需预先 /100

  alpha_results <- data.frame(
    Sample = rownames(abundance_matrix),
    # 丰富度：观测到的基因组数（非零值）
    Observed_Richness = apply(abundance_matrix, 1, function(x) sum(x > 0)),

    # Shannon指数
    Shannon = diversity(abundance_matrix, index = "shannon"),

    # Simpson指数
    Simpson = diversity(abundance_matrix, index = "simpson"),

    # Inverse Simpson指数
    InvSimpson = diversity(abundance_matrix, index = "invsimpson")
  )

  # Pielou均匀度：Shannon / ln(丰富度)，丰富度 <= 1 时无定义，设为 NA
  alpha_results$Pielou <- ifelse(
    alpha_results$Observed_Richness <= 1,
    NA,
    alpha_results$Shannon / log(alpha_results$Observed_Richness)
  )

  # 空样本（丰富度为0）的多样性指数无意义，设为 NA
  empty <- alpha_results$Observed_Richness == 0
  alpha_results[empty, c("Shannon", "Simpson", "InvSimpson")] <- NA

  # 添加有效物种数
  alpha_results$Effective_Species_Shannon <- exp(alpha_results$Shannon)
  alpha_results$Effective_Species_Simpson <- alpha_results$InvSimpson

  return(alpha_results)
}

# 计算α多样性
alpha_diversity <- calculate_alpha_diversity(rel_abundance_t)

print("α多样性结果:")
head(alpha_diversity)

# 排序查看
alpha_diversity_sorted <- alpha_diversity %>%
  arrange(desc(Shannon))

print("\n按Shannon多样性排序:")
print(alpha_diversity_sorted)

# 样本的环境因子信息
EF <- read.delim2("./info_corrected.txt",
                  sep = "\t")
EF <- type.convert(EF, as.is = TRUE)

# 样本的环境因子信息
EF_info <- read.delim2("./info_corrected.txt",
                       sep = "\t")
EF_info <- type.convert(EF, as.is = TRUE)

EF_info$medium <- ifelse(EF_info$Type == 'Cryoconite',"Glacier",EF_info$Type)
EF_info$medium <- ifelse(EF_info$medium == 'Ice',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'Snow',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'particle',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'water',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'sediment',"Benthon",EF_info$medium)


EF_info$ecosystem <- EF_info$River
EF_info$ecosystem <- str_remove(EF_info$ecosystem, "_.*")



alpha_diversity_sorted$Sample <- 
  sub("^[^.]*\\.", "", alpha_diversity_sorted$Sample)

alpha_diversity_sorted <- alpha_diversity_sorted %>%
  left_join(EF, by = c("Sample" = "SampleID"))

alpha_diversity_sorted$Plateau <- 
  ifelse(alpha_diversity_sorted$Plateau == 'Yes',"Plateau",alpha_diversity_sorted$Plateau)
alpha_diversity_sorted$Plateau <- 
  ifelse(alpha_diversity_sorted$Plateau == 'No',"Plain",alpha_diversity_sorted$Plateau)
alpha_diversity_sorted$Plateau <- 
  ifelse(alpha_diversity_sorted$Plateau == 'glacier',"Glacier",alpha_diversity_sorted$Plateau)

colnames(alpha_diversity_sorted)[14] <- 'Region'

alpha_diversity_sorted$Region <- factor(alpha_diversity_sorted$Region,
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

#统计检验#
library(rstatix)
# Shannon, Simpson, Observed（请根据实际列名修改）
diversity_indices <- c("Pielou","Shannon", 
                       "Simpson", "Observed_Richness")  
multi_test_MAG_richness <- alpha_diversity_sorted %>%
  # 1. 将几个多样性列“揉”成一个长表
  pivot_longer(
    cols = all_of(diversity_indices),
    names_to = "index",
    values_to = "value"
  ) %>%
  # 2. 按多样性指数分组，各自做完整的检验
  group_by(index) %>%
  wilcox_test(value ~ Region) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Region", dodge = 0.5) %>%
  ungroup()


# 分组箱线图
p_pielou_depth <- 
  ggplot(alpha_diversity_sorted, aes(x = Region, y = Pielou)) +
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

p_pielou_depth

p_shannon_depth <- 
  ggplot(alpha_diversity_sorted, aes(x = Region, y = Shannon)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
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

p_shannon_depth

p_simpson_depth <- 
  ggplot(alpha_diversity_sorted, aes(x = Region, y = Simpson)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
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

p_simpson_depth

p_richness_region <- 
  ggplot(alpha_diversity_sorted, aes(x = Region, y = Observed_Richness)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
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

p_richness_region



save.image("./01_MAG_and_phage_macrodiversity.RData")

##beta diversity####

# 1. 检查哪些行全为 0
empty_rows <- which(rowSums(rel_abundance_t) == 0)

# 2. 查看具体是哪些样本（如果有）
if (length(empty_rows) > 0) {
  cat("发现", length(empty_rows), "个空样本：\n")
  print(rownames(rel_abundance_t)[empty_rows])
  
  # 3. 移除这些空行
  rel_abundance_t_filtered <- rel_abundance_t[-empty_rows, ]
  
  cat("过滤后剩余样本数：", nrow(rel_abundance_t_filtered), "\n")
} else {
  rel_abundance_t_filtered <- rel_abundance_t
}

# 4. 使用过滤后的矩阵计算距离
bray_dist <- vegdist(rel_abundance_t_filtered, method = "bray")

# 计算Bray-Curtis距离
# bray_dist <- vegdist(rel_abundance_t, method = "bray")

# 查看距离矩阵
cat("Bray-Curtis距离矩阵:\n")
print(as.matrix(bray_dist)[1:5, 1:5])  # 显示前5个样本

# 转换为数据框便于分析
bray_matrix <- as.matrix(bray_dist)

library(ggplot2)
library(ggrepel)
library(vegan)

###pcoa####

# 1. 执行PCoA分析
pcoa_result <- cmdscale(bray_matrix, k = 3, eig = TRUE)

# 计算解释方差
eigenvalues <- pcoa_result$eig
variance_explained <- eigenvalues / sum(eigenvalues[eigenvalues > 0]) * 100

# 2. 准备绘图数据
pcoa_df <- data.frame(
  Sample = rownames(bray_matrix),
  PC1 = pcoa_result$points[, 1],
  PC2 = pcoa_result$points[, 2],
  PC3 = pcoa_result$points[, 3]
)

pcoa_df$Sample <- 
  sub("^[^.]*\\.", "", pcoa_df$Sample)

pcoa_df <- pcoa_df %>%
  left_join(EF_info, by = c("Sample" = "SampleID"))
pcoa_df$Plateau <- 
  ifelse(pcoa_df$Plateau == 'Yes',"Plateau",pcoa_df$Plateau)
pcoa_df$Plateau <- 
  ifelse(pcoa_df$Plateau == 'No',"Plain",pcoa_df$Plateau)
pcoa_df$Plateau <- 
  ifelse(pcoa_df$Plateau == 'glacier',"Glacier",pcoa_df$Plateau)

colnames(pcoa_df)[10] <- 'Region'

pcoa_df$Region <- factor(pcoa_df$Region,
                                        levels = c("Glacier","Plateau","Plain"))


# # 合并样本信息
# pcoa_df <- merge(pcoa_df, alpha_diversity, by = "Sample")
# pcoa_df$Depth <- factor(pcoa_df$Depth)

# 3. 2D PCoA图
p_pcoa2d <- ggplot(pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Region, fill = Region), size = 1, alpha = 1) +
  # geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
  #                 box.padding = 0.5, point.padding = 0.2) +
  stat_ellipse(aes(group = Region, color = Region), 
               level = 0.68, linetype = "dashed", linewidth = 0.5) +
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
  scale_shape_manual(values = 1:7) +
  labs(
    x = sprintf("PCoA 1 (%.1f%%)", variance_explained[1]),
    y = sprintf("PCoA 2 (%.1f%%)", variance_explained[2])
  ) +
  # theme_minimal(base_size = 12) +
  # theme(
  #   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #   plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
  #   legend.position = "right",
  #   panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  #   panel.grid.minor = element_blank(),
  #   aspect.ratio = 1
  # ) +
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

print(p_pcoa2d)
write.table(pcoa_df, "pcoa_df_mag.tsv",sep = "\t",quote = F, row.names = F)
write.table(alpha_diversity_sorted, 
  "alpha_diversity_sorted_mag.tsv",sep = "\t",quote = F, row.names = F)

# 4. 3D PCoA图（可选）
library(plotly)
if(require(plotly)) {
  p_pcoa3d <- plot_ly(pcoa_df, 
                      x = ~PC1, y = ~PC2, z = ~PC3,
                      color = ~Region,
                      colors = c("#9BDFDF", "#C4D8E9", "#d4e09b"),
                      # text = ~paste("样本:", Sample,
                      #               "<br>深度:", Depth, "m",
                      #               "<br>位置:", Location),
                      hoverinfo = "text",
                      marker = list(size = 1)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = paste0("PCoA 1 (", round(variance_explained[1],1), "%)")),
                        yaxis = list(title = paste0("PCoA 2 (", round(variance_explained[2],1), "%)")),
                        zaxis = list(title = paste0("PCoA 3 (", round(variance_explained[3],1), "%)"))),
           title = "β diversity 3D PCoA")
  
  # 保存为HTML
  htmlwidgets::saveWidget(p_pcoa3d, "beta_diversity_3d_pcoa.html")
}
p_pcoa3d



###NMDS####

# 1. NMDS分析
set.seed(123)
nmds_result <- 
  metaMDS(bray_matrix, k = 2, trymax = 100, distance = "bray")

cat("NMDS应力值:", nmds_result$stress, "\n")
cat("应力值解释:\n")
cat("  <0.05: 极好\n")
cat("  0.05-0.1: 好\n")
cat("  0.1-0.2: 可接受\n")
cat("  >0.2: 差\n")

# 2. 准备绘图数据
nmds_scores <- scores(nmds_result, display = "sites")
nmds_df <- data.frame(
  Sample = rownames(nmds_scores),
  NMDS1 = nmds_scores[, 1],
  NMDS2 = nmds_scores[, 2]
)

nmds_df$Sample <- 
  sub("^[^.]*\\.", "", nmds_df$Sample)

nmds_df <- nmds_df %>%
  left_join(EF_info, by = c("Sample" = "SampleID"))
nmds_df$Plateau <- 
  ifelse(nmds_df$Plateau == 'Yes',"Plateau",nmds_df$Plateau)
nmds_df$Plateau <- 
  ifelse(nmds_df$Plateau == 'No',"Plain",nmds_df$Plateau)
nmds_df$Plateau <- 
  ifelse(nmds_df$Plateau == 'glacier',"Glacier",nmds_df$Plateau)

colnames(nmds_df)[9] <- 'Region'

nmds_df$Region <- factor(nmds_df$Region,
                         levels = c("Glacier","Plateau","Plain"))


# nmds_df <- merge(nmds_df, alpha_diversity, by = "Sample")
# nmds_df$Depth <- factor(nmds_df$Depth)


# 3. 绘制NMDS图
p_nmds <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Region, shape = Region), size = 5, alpha = 0.8) +
#  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
  stat_ellipse(aes(group = Region), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )

print(p_nmds)

library(Cairo)
library(cowplot)

g1_MAG_macrodiversity <- 
  plot_grid(p_pielou_depth,p_shannon_depth,
            p_simpson_depth,p_richness_region,
            p_pcoa2d,
            labels = letters[1:5], ncol = 2,
            align = "v",  # 垂直对齐
            axis = "lr",  # 左和右对齐
            rel_widths = c(1, 1,1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g1_MAG_macrodiversity

ggsave(plot = g1_MAG_macrodiversity,
       '1_g1_MAG_macrodiversity.pdf',width = 6,height = 6,units = 'in')

save.image("./01_MAG_and_phage_macrodiversity.RData")

load("./01_MAG_and_phage_macrodiversity.RData")


#Phage diversity ####

##alpha diversity ####

phage_relative_abundance <- 
  read.delim2("./丰度/Virus/all_abundance_rpkm.tsv")

# mag_data_rel_only <- MAGs_544_relative_abundance %>%
#   select(Genome, contains("Relative.Abundance"))
# 将所有非Genome列转换为数值型
phage_relative_abundance <- phage_relative_abundance %>%
  mutate(across(-Contig, as.numeric))

library(tidyverse)
library(vegan)

# 1. 提取相对丰度矩阵
# 假设所有以"Relative.Abundance"结尾的列都是样本的相对丰度
phage_rel_abundance_cols <- 
  grep("sorted\\.RPKM", colnames(phage_relative_abundance), value = TRUE)

# 2. 创建样本名（简化列名）
sample_names <- 
  gsub("\\.sorted\\.RPKM\\.\\.\\.\\.", "", phage_rel_abundance_cols)
sample_names <- gsub(".sorted.RPKM", "", sample_names)  # 可选：进一步简化
#sample_names <- gsub("_", " ", sample_names)      # 可选：用空格替代下划线

# 3. 提取相对丰度数据
phage_rel_abundance_matrix <- phage_relative_abundance[, phage_rel_abundance_cols]

# 4. 设置行名和列名
rownames(phage_rel_abundance_matrix) <- phage_relative_abundance$Contig
colnames(phage_rel_abundance_matrix) <- sample_names

# 5. 查看数据结构
cat("数据维度:", dim(phage_rel_abundance_matrix), "\n")
cat("样本数:", ncol(phage_rel_abundance_matrix), "\n")
cat("基因组数:", nrow(phage_rel_abundance_matrix), "\n")
cat("样本名称:", colnames(phage_rel_abundance_matrix), "\n")

# 6. 移除"unmapped"行（如果需要）
phage_rel_abundance_matrix <- 
  phage_rel_abundance_matrix[rownames(phage_rel_abundance_matrix) != "unmapped", ]

# 7. 转置矩阵：vegan期望行为样本，列为物种
phage_rel_abundance_t <- t(phage_rel_abundance_matrix)

# 方法1：使用vegan包直接计算
calculate_alpha_diversity <- function(abundance_matrix) {
  # 注意：abundance_matrix 应为行为样本、列为物种
  # vegan 内部会将每行标准化为相对丰度（和为1），因此无需预先 /100

  alpha_results <- data.frame(
    Sample = rownames(abundance_matrix),
    # 丰富度：观测到的基因组数（非零值）
    Observed_Richness = apply(abundance_matrix, 1, function(x) sum(x > 0)),

    # Shannon指数
    Shannon = diversity(abundance_matrix, index = "shannon"),

    # Simpson指数
    Simpson = diversity(abundance_matrix, index = "simpson"),

    # Inverse Simpson指数
    InvSimpson = diversity(abundance_matrix, index = "invsimpson")
  )

  # Pielou均匀度：Shannon / ln(丰富度)，丰富度 <= 1 时无定义，设为 NA
  alpha_results$Pielou <- ifelse(
    alpha_results$Observed_Richness <= 1,
    NA,
    alpha_results$Shannon / log(alpha_results$Observed_Richness)
  )

  # 空样本（丰富度为0）的多样性指数无意义，设为 NA
  empty <- alpha_results$Observed_Richness == 0
  alpha_results[empty, c("Shannon", "Simpson", "InvSimpson")] <- NA

  # 添加有效物种数
  alpha_results$Effective_Species_Shannon <- exp(alpha_results$Shannon)
  alpha_results$Effective_Species_Simpson <- alpha_results$InvSimpson

  return(alpha_results)
}

# 计算α多样性
phage_alpha_diversity <- calculate_alpha_diversity(phage_rel_abundance_t)

print("α多样性结果:")
head(phage_alpha_diversity)

# 排序查看
phage_alpha_diversity_sorted <- phage_alpha_diversity %>%
  arrange(desc(Shannon))

print("\n按Shannon多样性排序:")
print(phage_alpha_diversity_sorted)

# 样本的环境因子信息
EF_info <- read.delim2("./info_corrected.txt",
                  sep = "\t")
EF_info <- type.convert(EF, as.is = TRUE)

EF_info$medium <- ifelse(EF_info$Type == 'Cryoconite',"Glacier",EF_info$Type)
EF_info$medium <- ifelse(EF_info$medium == 'Ice',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'Snow',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'particle',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'water',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'sediment',"Benthon",EF_info$medium)


EF_info$ecosystem <- EF_info$River
EF_info$ecosystem <- str_remove(EF_info$ecosystem, "_.*")


phage_alpha_diversity_sorted$Sample <- 
  sub("^[^.]*\\.", "", phage_alpha_diversity_sorted$Sample)

phage_alpha_diversity_sorted <- phage_alpha_diversity_sorted %>%
  left_join(EF, by = c("Sample" = "SampleID"))

phage_alpha_diversity_sorted$Plateau <- 
  ifelse(phage_alpha_diversity_sorted$Plateau == 'Yes',"Plateau",phage_alpha_diversity_sorted$Plateau)
phage_alpha_diversity_sorted$Plateau <- 
  ifelse(phage_alpha_diversity_sorted$Plateau == 'No',"Plain",phage_alpha_diversity_sorted$Plateau)
phage_alpha_diversity_sorted$Plateau <- 
  ifelse(phage_alpha_diversity_sorted$Plateau == 'glacier',"Glacier",phage_alpha_diversity_sorted$Plateau)

colnames(phage_alpha_diversity_sorted)[14] <- 'Region'

phage_alpha_diversity_sorted$Region <- factor(phage_alpha_diversity_sorted$Region,
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


#统计检验#
library(rstatix)
# Shannon, Simpson, Observed（请根据实际列名修改）
diversity_indices <- c("Pielou","Shannon", 
                       "Simpson", "Observed_Richness")  
multi_test_phage_richness <- phage_alpha_diversity_sorted %>%
  # 1. 将几个多样性列“揉”成一个长表
  pivot_longer(
    cols = all_of(diversity_indices),
    names_to = "index",
    values_to = "value"
  ) %>%
  # 2. 按多样性指数分组，各自做完整的检验
  group_by(index) %>%
  wilcox_test(value ~ Region) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Region", dodge = 0.5) %>%
  ungroup()

save.image("./01_MAG_and_phage_macrodiversity_new.RData")
load("./01_MAG_and_phage_macrodiversity_new.RData")

# 分组箱线图
p_phage_pielou_region <- 
  ggplot(phage_alpha_diversity_sorted, aes(x = Region, y = Pielou)) +
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

p_phage_pielou_region

p_phage_shannon_region <- 
  ggplot(phage_alpha_diversity_sorted, aes(x = Region, y = Shannon)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
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

p_phage_shannon_region

p_phage_simpson_region <- 
  ggplot(phage_alpha_diversity_sorted, aes(x = Region, y = Simpson)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
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

p_phage_simpson_region

p_phage_richness_region <- 
  ggplot(phage_alpha_diversity_sorted, aes(x = Region, y = Observed_Richness)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
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

p_phage_richness_region



save.image("./01_MAG_and_phage_macrodiversity.RData")

##beta diversity####
# 计算Bray-Curtis距离
phage_bray_dist <- vegdist(phage_rel_abundance_t, method = "bray")

# 查看距离矩阵
cat("Bray-Curtis距离矩阵:\n")
print(as.matrix(phage_bray_dist)[1:5, 1:5])  # 显示前5个样本

# 转换为数据框便于分析
phage_bray_matrix <- as.matrix(phage_bray_dist)

library(ggplot2)
library(ggrepel)
library(vegan)

###pcoa####

# 1. 执行PCoA分析
phage_pcoa_result <- cmdscale(phage_bray_matrix, k = 3, eig = TRUE)

# 计算解释方差
phage_eigenvalues <- phage_pcoa_result$eig
phage_variance_explained <- 
  phage_eigenvalues / sum(phage_eigenvalues[phage_eigenvalues > 0]) * 100

# 2. 准备绘图数据
phage_pcoa_df <- data.frame(
  Sample = rownames(phage_bray_matrix),
  PC1 = phage_pcoa_result$points[, 1],
  PC2 = phage_pcoa_result$points[, 2],
  PC3 = phage_pcoa_result$points[, 3]
)

phage_pcoa_df$Sample <- 
  sub("^[^.]*\\.", "", phage_pcoa_df$Sample)

phage_pcoa_df <- phage_pcoa_df %>%
  left_join(EF, by = c("Sample" = "SampleID"))
phage_pcoa_df$Plateau <- 
  ifelse(phage_pcoa_df$Plateau == 'Yes',"Plateau",phage_pcoa_df$Plateau)
phage_pcoa_df$Plateau <- 
  ifelse(phage_pcoa_df$Plateau == 'No',"Plain",phage_pcoa_df$Plateau)
phage_pcoa_df$Plateau <- 
  ifelse(phage_pcoa_df$Plateau == 'glacier',"Glacier",phage_pcoa_df$Plateau)

colnames(phage_pcoa_df)[13] <- 'Region'

phage_pcoa_df$Region <- factor(phage_pcoa_df$Region,
                         levels = c("Glacier","Plateau","Plain"))


# # 合并样本信息
# pcoa_df <- merge(pcoa_df, alpha_diversity, by = "Sample")
# pcoa_df$Depth <- factor(pcoa_df$Depth)

# 3. 2D PCoA图
p_phage_pcoa2d <- ggplot(phage_pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Region, fill = Region), size = 1, alpha = 1) +
  # geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
  #                 box.padding = 0.5, point.padding = 0.2) +
  stat_ellipse(aes(group = Region, color = Region), 
               level = 0.68, linetype = "dashed", linewidth = 0.5) +
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
  scale_shape_manual(values = 1:7) +
  labs(
    x = sprintf("PCoA 1 (%.1f%%)", phage_variance_explained[1]),
    y = sprintf("PCoA 2 (%.1f%%)", phage_variance_explained[2])
  ) +
  # theme_minimal(base_size = 12) +
  # theme(
  #   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #   plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
  #   legend.position = "right",
  #   panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  #   panel.grid.minor = element_blank(),
  #   aspect.ratio = 1
  # ) +
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

print(p_phage_pcoa2d)

# 4. 3D PCoA图（可选）
library(plotly)
if(require(plotly)) {
  p_phage_pcoa3d <- plot_ly(phage_pcoa_df, 
                      x = ~PC1, y = ~PC2, z = ~PC3,
                      color = ~Region,
                      colors = c("#9BDFDF", "#C4D8E9", "#d4e09b"),
                      # text = ~paste("样本:", Sample,
                      #               "<br>深度:", Depth, "m",
                      #               "<br>位置:", Location),
                      hoverinfo = "text",
                      marker = list(size = 1)) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = paste0("PCoA 1 (", round(phage_variance_explained[1],1), "%)")),
                        yaxis = list(title = paste0("PCoA 2 (", round(phage_variance_explained[2],1), "%)")),
                        zaxis = list(title = paste0("PCoA 3 (", round(phage_variance_explained[3],1), "%)"))),
           title = "Phage β diversity 3D PCoA")
  
  # 保存为HTML
  htmlwidgets::saveWidget(p_phage_pcoa3d, "phage_beta_diversity_3d_pcoa.html")
}
p_phage_pcoa3d

write.table(phage_pcoa_df, "phage_pcoa_df.tsv",sep = "\t",quote = F, row.names = F)
write.table(phage_alpha_diversity_sorted, 
  "phage_alpha_diversity_sorted.tsv",sep = "\t",quote = F, row.names = F)


###NMDS####

# 1. NMDS分析
set.seed(123)
phage_nmds_result <- 
  metaMDS(phage_bray_matrix, k = 2, trymax = 100, distance = "bray")

cat("NMDS应力值:", phage_nmds_result$stress, "\n")
cat("应力值解释:\n")
cat("  <0.05: 极好\n")
cat("  0.05-0.1: 好\n")
cat("  0.1-0.2: 可接受\n")
cat("  >0.2: 差\n")

# 2. 准备绘图数据
phage_nmds_scores <- scores(phage_nmds_result, display = "sites")
phage_nmds_df <- data.frame(
  Sample = rownames(phage_nmds_scores),
  NMDS1 = phage_nmds_scores[, 1],
  NMDS2 = phage_nmds_scores[, 2]
)

phage_nmds_df$Sample <- 
  sub("^[^.]*\\.", "", phage_nmds_df$Sample)

phage_nmds_df <- phage_nmds_df %>%
  left_join(EF_info, by = c("Sample" = "SampleID"))
phage_nmds_df$Plateau <- 
  ifelse(phage_nmds_df$Plateau == 'Yes',"Plateau",phage_nmds_df$Plateau)
phage_nmds_df$Plateau <- 
  ifelse(phage_nmds_df$Plateau == 'No',"Plain",phage_nmds_df$Plateau)
phage_nmds_df$Plateau <- 
  ifelse(phage_nmds_df$Plateau == 'glacier',"Glacier",phage_nmds_df$Plateau)

colnames(phage_nmds_df)[9] <- 'Region'

phage_nmds_df$Region <- factor(phage_nmds_df$Region,
                         levels = c("Glacier","Plateau","Plain"))


# nmds_df <- merge(nmds_df, alpha_diversity, by = "Sample")
# nmds_df$Depth <- factor(nmds_df$Depth)


# 3. 绘制NMDS图
p_phage_nmds <- ggplot(phage_nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Region), size = 5, alpha = 0.8) +
  #  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
  stat_ellipse(aes(group = Region), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )+
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

print(p_phage_nmds)


library(Cairo)
library(cowplot)

g2_phage_macrodiversity <- 
  plot_grid(p_phage_pielou_region,p_phage_shannon_region,
            p_phage_simpson_region,p_phage_richness_region,
            p_phage_pcoa2d,p_phage_pcoa3d,
            labels = letters[1:5], ncol = 2,
            align = "v",  # 垂直对齐
            axis = "lr",  # 左和右对齐
            rel_widths = c(1, 1,1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
g2_phage_macrodiversity

ggsave(plot = g2_phage_macrodiversity,
       '1_g2_phage_macrodiversity.pdf',width = 6,height = 6,units = 'in')

save.image("./01_MAG_and_phage_macrodiversity.RData")
load("./01_MAG_and_phage_macrodiversity.RData")





#phage 在不同分类情况下的beta diversity####

# 3. 绘制NMDS图
p_phage_nmds_region <- ggplot(phage_nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Region), size = 1, alpha = 1) +
  #  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
  stat_ellipse(aes(group = Region), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", phage_nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )+
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

p_phage_nmds_region

p_phage_nmds_season <- ggplot(phage_nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Season), size = 1, alpha = 1) +
  #  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
  stat_ellipse(aes(group = Season), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", phage_nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )+
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

p_phage_nmds_season
ggsave(plot=p_phage_nmds_season,
       "p_phage_nmds_season_legend.pdf",
       width = 6, height = 6, units = "in")

p_phage_nmds_ecosystem <- ggplot(phage_nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = ecosystem), size = 1, alpha = 1) +
  #  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
  stat_ellipse(aes(group = ecosystem), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","YR"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", phage_nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )+
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
p_phage_nmds_ecosystem

p_phage_nmds_medium <- ggplot(phage_nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = medium), size = 1, alpha = 1) +
  #  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
  stat_ellipse(aes(group = medium), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","YR"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", phage_nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
    legend.position = "right"
  )+
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
p_phage_nmds_medium

library(Cairo)
library(cowplot)

group_plot3_phage_nmds <- 
  plot_grid(p_phage_nmds_region,p_phage_nmds_medium,
            p_phage_nmds_season,p_phage_nmds_ecosystem,
            labels = letters[1:5], ncol = 2,
            align = "v",  # 垂直对齐
            axis = "lr",  # 左和右对齐
            rel_widths = c(1, 1,1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
group_plot3_phage_nmds

ggsave(plot = group_plot3_phage_nmds,
       '1_group_plot3_phage_nmds.pdf',width = 9,height = 6,units = 'in')

write.table(phage_nmds_df, "phage_nmds_df.tsv", sep = "\t", quote = F, row.names = F)


## ANOSIM ####
# 假设 phage_nmds_df 已包含分组信息，且 Sample 列与距离矩阵的样本名一致
# 距离矩阵的行名就是样本名
sample_names_2 <- rownames(phage_bray_matrix)  # 或 attr(phage_bray_dist, "Labels")
sample_names_2 <- sub("^[^.]*\\.", "", sample_names_2)

# 从 phage_nmds_df 中匹配分组
group_region <- 
  phage_nmds_df$Region[match(sample_names_2, phage_nmds_df$Sample)]
group_season <- 
  phage_nmds_df$Season[match(sample_names_2, phage_nmds_df$Sample)]
group_ecosystem <- 
  phage_nmds_df$ecosystem[match(sample_names_2, phage_nmds_df$Sample)]
group_medium <- 
  phage_nmds_df$medium[match(sample_names_2, phage_nmds_df$Sample)]

library(vegan)

# 对 Region 分组
anosim_region <- anosim(phage_bray_dist, 
                        grouping = group_region, 
                        permutations = 999)
print(anosim_region)
# ANOSIM statistic R: 0.7358 
# Significance: 0.001 
# 对 Season
anosim_season <- 
  anosim(phage_bray_dist, grouping = group_season, permutations = 999)
print(anosim_season)
# ANOSIM statistic R: 0.3066 
# Significance: 0.001

# 对 ecosystem
anosim_ecosystem <- 
  anosim(phage_bray_dist, grouping = group_ecosystem, permutations = 999)
print(anosim_ecosystem)
# ANOSIM statistic R: 0.5579 
# Significance: 0.001
# 对 medium
anosim_medium <- 
  anosim(phage_bray_dist, grouping = group_medium, permutations = 999)
print(anosim_medium)
# ANOSIM statistic R: 0.4477 
# Significance: 0.001

plot(anosim_region, main = "ANOSIM: Region")
plot(anosim_medium, main = "ANOSIM: Medium")
plot(anosim_ecosystem, main = "ANOSIM: Ecosystem")
plot(anosim_season, main = "ANOSIM: Season")

## ADONIS ####
# library(vegan)
# 
# # 1. 匹配样本顺序
# sample_order <- rownames(phage_bray_matrix)
# sample_order <- sub("^[^.]*\\.", "", sample_order)
# 
# df_matched <- 
#   phage_nmds_df[match(sample_order, 
#                       phage_nmds_df$Sample), ]
# 
# # 2. 多因子 ADONIS（边际效应）
# adonis_result <- adonis2(
#   phage_bray_dist ~ Region + Season + ecosystem + medium,
#   data = df_matched,
#   permutations = 999,
#   by = "margin"
# )
# print(adonis_result)
# 
# # 3. 方差齐性检验
# disper_test <- betadisper(phage_bray_dist, df_matched$Region)
# permutest(disper_test)
# 
# # 4. （可选）两两比较，以 Region 为例
# if (!require(pairwiseAdonis)) install.packages("pairwiseAdonis")
# library(pairwiseAdonis)
# pairwise.adonis2(phage_bray_dist ~ Region, data = df_matched)

# MAG 在不同分类情况下的 beta多样性####

p_NMDS2d_region <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Region, fill = Region), size = 1, alpha = 1) +
  # geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
  #                 box.padding = 0.5, point.padding = 0.2) +
  # stat_ellipse(aes(group = Region, color = Region), 
  #              level = 0.68, linetype = "dashed", linewidth = 0.5) +
  stat_ellipse(aes(group = Region), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
  #                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  # )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","YR"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  # theme_minimal(base_size = 12) +
  # theme(
  #   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #   plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
  #   legend.position = "right",
  #   panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  #   panel.grid.minor = element_blank(),
  #   aspect.ratio = 1
  # ) +
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

p_NMDS2d_region


p_NMDS2d_medium <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = medium, fill = medium), size = 1, alpha = 1) +
  # geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
  #                 box.padding = 0.5, point.padding = 0.2) +
  # stat_ellipse(aes(group = medium, color = medium), 
  #              level = 0.68, linetype = "dashed", linewidth = 0.5) +
  stat_ellipse(aes(group = medium), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
  #                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  # )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","YR"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  # theme_minimal(base_size = 12) +
  # theme(
  #   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #   plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
  #   legend.position = "right",
  #   panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  #   panel.grid.minor = element_blank(),
  #   aspect.ratio = 1
  # ) +
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

p_NMDS2d_medium

p_NMDS2d_season <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Season, fill = Season), size = 1, alpha = 1) +
  # geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
  #                 box.padding = 0.5, point.padding = 0.2) +
  # stat_ellipse(aes(group = Season, color = Season), 
  #              level = 0.68, linetype = "dashed", linewidth = 0.5) +
  stat_ellipse(aes(group = Season), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  # theme_minimal(base_size = 12) +
  # theme(
  #   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #   plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
  #   legend.position = "right",
  #   panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  #   panel.grid.minor = element_blank(),
  #   aspect.ratio = 1
  # ) +
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
p_NMDS2d_season

p_NMDS2d_ecosystem <- ggplot(nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = ecosystem, fill = ecosystem), size = 1, alpha = 1) +
  # geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 20,
  #                 box.padding = 0.5, point.padding = 0.2) +
  # stat_ellipse(aes(group = ecosystem, color = ecosystem), 
  #              level = 0.68, linetype = "dashed", linewidth = 0.5) +
  stat_ellipse(aes(group = ecosystem), 
               level = 0.68, linetype = "dashed", linewidth = 0.8) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","YR"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_shape_manual(values = 1:7) +
  labs(
    subtitle = sprintf("stress = %.3f", nmds_result$stress),
    x = "NMDS 1",
    y = "NMDS 2"
  ) +
  # theme_minimal(base_size = 12) +
  # theme(
  #   plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
  #   plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray40"),
  #   legend.position = "right",
  #   panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
  #   panel.grid.minor = element_blank(),
  #   aspect.ratio = 1
  # ) +
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
p_NMDS2d_ecosystem

write.table(nmds_df, "nmds_df.tsv", sep = "\t", quote = F, row.names = F)

## ANOSIM ####
# 假设 phage_nmds_df 已包含分组信息，且 Sample 列与距离矩阵的样本名一致
# 距离矩阵的行名就是样本名
sample_names_3 <- rownames(bray_matrix)  # 或 attr(phage_bray_dist, "Labels")
sample_names_3 <- sub("^[^.]*\\.", "", sample_names_3)

# 从 pcoa_df 中匹配分组
group_region_mag <- 
  pcoa_df$Region[match(sample_names_3, pcoa_df$Sample)]
group_season_mag <- 
  pcoa_df$Season[match(sample_names_3, pcoa_df$Sample)]
group_ecosystem_mag <- 
  pcoa_df$ecosystem[match(sample_names_3, pcoa_df$Sample)]
group_medium_mag <- 
  pcoa_df$medium[match(sample_names_3, pcoa_df$Sample)]

library(vegan)

# 对 Region 分组
anosim_region_mag <- anosim(bray_dist, 
                        grouping = group_region_mag, 
                        permutations = 999)
print(anosim_region_mag)
# ANOSIM statistic R: 0.6287
# Significance: 0.001 
# 对 Season
anosim_season_mag <- 
  anosim(bray_dist, 
         grouping = group_season_mag, permutations = 999)
print(anosim_season_mag)
# ANOSIM statistic R: 0.2453
# Significance: 0.001

# 对 ecosystem
anosim_ecosystem_mag <- 
  anosim(bray_dist, grouping = group_ecosystem_mag, permutations = 999)
print(anosim_ecosystem_mag)
# ANOSIM statistic R: 0.4367
# Significance: 0.001
# 对 medium
anosim_medium_mag <- 
  anosim(bray_dist, grouping = group_medium_mag, permutations = 999)
print(anosim_medium_mag)
# ANOSIM statistic R: 0.7987
# Significance: 0.001

plot(anosim_region_mag, main = "ANOSIM: Region")
plot(anosim_medium_mag, main = "ANOSIM: Medium")
plot(anosim_ecosystem_mag, main = "ANOSIM: Ecosystem")
plot(anosim_season_mag, main = "ANOSIM: Season")

library(Cairo)
library(cowplot)

group_plot4_mag_nmds <- 
  plot_grid(p_NMDS2d_region,p_NMDS2d_medium,
            p_NMDS2d_ecosystem,p_NMDS2d_season,
            labels = letters[1:5], ncol = 2,
            align = "v",  # 垂直对齐
            axis = "lr",  # 左和右对齐
            rel_widths = c(1, 1,1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
group_plot4_mag_nmds

ggsave(plot = group_plot4_mag_nmds,
       '1_group_plot4_mag_nmds.pdf',width = 9,height = 6,units = 'in')











# MAG TPM relative abundance####
# MAGs_4903_tpm <- 
#   read.delim2("./丰度/Bins/all_abundance_tpm_bin.tsv")
# 
# # mag_data_rel_only <- MAGs_544_relative_abundance %>%
# #   select(Genome, contains("Relative.Abundance"))
# # 将所有非Genome列转换为数值型
# MAGs_4903_tpm <- MAGs_4903_tpm %>%
#   mutate(across(-Genome, as.numeric))
# 
# library(tidyverse)
# library(vegan)
# 
# # 1. 提取相对丰度矩阵
# # 假设所有以"Relative.Abundance"结尾的列都是样本的相对丰度
# rel_abundance_cols_tpm <- 
#   grep("sorted\\.TPM", 
#        colnames(MAGs_4903_tpm), value = TRUE)
# 
# # 2. 创建样本名（简化列名）
# sample_names_tpm <- 
#   gsub("\\.sorted\\.TPM\\.\\.\\.\\.", "", rel_abundance_cols_tpm)
# sample_names_tpm <- gsub(".sorted.TPM", "", sample_names_tpm)  # 可选：进一步简化
# #sample_names <- gsub("_", " ", sample_names)      # 可选：用空格替代下划线
# 
# # 3. 提取相对丰度数据
# rel_abundance_matrix_tpm <- 
#   MAGs_4903_tpm[, rel_abundance_cols_tpm]
# 
# # 4. 设置行名和列名
# rownames(rel_abundance_matrix_tpm) <- 
#   MAGs_4903_tpm$Genome
# colnames(rel_abundance_matrix_tpm) <- 
#   sample_names_tpm
# 
# # 5. 查看数据结构
# cat("数据维度:", dim(rel_abundance_matrix_tpm), "\n")
# cat("样本数:", ncol(rel_abundance_matrix_tpm), "\n")
# cat("基因组数:", nrow(rel_abundance_matrix_tpm), "\n")
# cat("样本名称:", colnames(rel_abundance_matrix_tpm), "\n")
# 
# # 6. 移除"unmapped"行（如果需要）
# rel_abundance_matrix_tpm <- 
#   rel_abundance_matrix_tpm[rownames(rel_abundance_matrix_tpm) != "unmapped", ]
# 
# # 7. 转置矩阵：vegan期望行为样本，列为物种
# rel_abundance_tpm_t <- t(rel_abundance_matrix_tpm)
# 
# # 1. 检查哪些行全为 0
# empty_rows <- which(rowSums(rel_abundance_tpm_t) == 0)
# 
# # 2. 查看具体是哪些样本（如果有）
# if (length(empty_rows) > 0) {
#   cat("发现", length(empty_rows), "个空样本：\n")
#   print(rownames(rel_abundance_tpm_t)[empty_rows])
#   
#   # 3. 移除这些空行
#   rel_abundance_tpm_t_filtered <- rel_abundance_tpm_t[-empty_rows, ]
#   
#   cat("过滤后剩余样本数：", nrow(rel_abundance_tpm_t_filtered), "\n")
# } else {
#   rel_abundance_tpm_t_filtered <- rel_abundance_tpm_t
# }
# 
# # 4. 使用过滤后的矩阵计算距离
# bray_dist_tpm <- vegdist(rel_abundance_tpm_t_filtered, method = "bray")
# 
# #bray_dist_tpm <- vegdist(rel_abundance_tpm_t, method = "bray")
# bray_matrix_tpm <- as.matrix(bray_dist_tpm)
# 
# # 1. NMDS分析
# set.seed(123)
# nmds_result_tpm <- 
#   metaMDS(bray_matrix_tpm, k = 2, trymax = 100, distance = "bray")
# 
# # cat("NMDS应力值:", nmds_result$stress, "\n")
# # cat("应力值解释:\n")
# # cat("  <0.05: 极好\n")
# # cat("  0.05-0.1: 好\n")
# # cat("  0.1-0.2: 可接受\n")
# # cat("  >0.2: 差\n")
# 
# # 2. 准备绘图数据
# nmds_scores_tpm <- scores(nmds_result_tpm, display = "sites")
# nmds_df_tpm <- data.frame(
#   Sample = rownames(nmds_scores_tpm),
#   NMDS1 = nmds_scores_tpm[, 1],
#   NMDS2 = nmds_scores_tpm[, 2]
# )
# 
# nmds_df_tpm$Sample <- 
#   sub("^[^.]*\\.", "", nmds_df_tpm$Sample)
# 
# nmds_df_tpm <- nmds_df_tpm %>%
#   left_join(EF_info, by = c("Sample" = "SampleID"))
# nmds_df_tpm$Plateau <- 
#   ifelse(nmds_df_tpm$Plateau == 'Yes',"Plateau",nmds_df_tpm$Plateau)
# nmds_df_tpm$Plateau <- 
#   ifelse(nmds_df_tpm$Plateau == 'No',"Plain",nmds_df_tpm$Plateau)
# nmds_df_tpm$Plateau <- 
#   ifelse(nmds_df_tpm$Plateau == 'glacier',"Glacier",nmds_df_tpm$Plateau)
# 
# colnames(nmds_df_tpm)[9] <- 'Region'
# 
# nmds_df_tpm$Region <- factor(nmds_df_tpm$Region,
#                          levels = c("Glacier","Plateau","Plain"))
# 
# 
# # nmds_df <- merge(nmds_df, alpha_diversity, by = "Sample")
# # nmds_df$Depth <- factor(nmds_df$Depth)
# 
# 
# # 3. 绘制NMDS图
# p_nmds_tpm <- ggplot(nmds_df_tpm, aes(x = NMDS1, y = NMDS2)) +
#   geom_point(aes(color = Region, shape = Region), size = 5, alpha = 0.8) +
#   #  geom_text_repel(aes(label = Sample), size = 3.5, max.overlaps = 15) +
#   stat_ellipse(aes(group = Region), 
#                level = 0.68, linetype = "dashed", linewidth = 0.8) +
#   scale_fill_brewer(palette = "Set2") +
#   scale_color_brewer(palette = "Set2") +
#   scale_shape_manual(values = 1:7) +
#   labs(
#     subtitle = sprintf("stress = %.3f", nmds_result$stress),
#     x = "NMDS 1",
#     y = "NMDS 2"
#   ) +
#   theme_minimal() +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold"),
#     plot.subtitle = element_text(hjust = 0.5, color = "gray40"),
#     legend.position = "right"
#   )
# 
# print(p_nmds_tpm)





save.image("./01_MAG_and_phage_macrodiversity.RData")
load("./01_MAG_and_phage_macrodiversity.RData")
