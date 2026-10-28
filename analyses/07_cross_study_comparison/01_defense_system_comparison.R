# =============================================================================
# Glacier-to-river -- cross study comparison / Rivers vs aquifers: defense systems
# =============================================================================
# Compare the river/glacier defensome with the aquifer defensome of a
# companion study.
#
# Manuscript      : -
# Main outputs    : 1_river_aquifer_p.pdf
# Original script : 20gaoyuan/40plateau_groundwater_comparision/1_defense_comparision.R
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
GTR_GW <- normalizePath(Sys.getenv("GTR_GW", unset = file.path(GTR_DATA, "groundwater")), mustWork = FALSE)

setwd(file.path(GTR_DATA, "20gaoyuan/40plateau_groundwater_comparision"))

#基础数据####
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

genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
count_matrix <- count_matrix[rownames(count_matrix) %in% MAG_metadata$genome,]
DF_defense <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]

MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]


#只挑河流，90 5 ，N50 100kb####
MAG_metadata_river <- MAG_metadata[! MAG_metadata$Type %in% c("ice","Snow",
                                                             "Cryoconite"),]
MAG_quality_100kb <- MAG_quality[MAG_quality$Contig_N50>=100000,]

MAG_metadata_river_100kb <- MAG_metadata_river[MAG_metadata_river$genome %in% 
                                                 MAG_quality_100kb$Name,]

genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% 
                                           MAG_metadata_river_100kb$genome,]

MAG_metadata_river_100kb <- merge(MAG_metadata_river_100kb, genome_defense_num,
                                  by = "genome",all=T)

MAG_metadata_river_100kb[is.na(MAG_metadata_river_100kb[, 62]), 62] <- 0
mean(MAG_metadata_river_100kb$DS_num) #4.19

MAG_quality_100kb <- MAG_quality_100kb[MAG_quality_100kb$Name %in% 
                                         MAG_metadata_river_100kb$genome,]
MAG_metadata_river_100kb <- merge(MAG_metadata_river_100kb, MAG_quality_100kb, 
                                  by.x = "genome", by.y = "Name")
MAG_metadata_river_100kb$DS_density_kb <- 
  MAG_metadata_river_100kb$DS_num / (MAG_metadata_river_100kb$Genome_Size)*1000

MAG_metadata_river_100kb$DS_density_kblog <- 
  log10(MAG_metadata_river_100kb$DS_density_kb+1)*1000

MAG_metadata_river_100kb$river_aquifer <- "River"

count_matrix_90_5_100kb_river <- count_matrix
count_matrix_90_5_100kb_river$genome <- rownames(count_matrix_90_5_100kb_river)
count_matrix_90_5_100kb_river <- count_matrix_90_5_100kb_river[count_matrix_90_5_100kb_river$genome %in% 
                                                                 MAG_metadata_river_100kb$genome,]
count_matrix_90_5_100kb_river <- count_matrix_90_5_100kb_river[,-189]
count_matrix_90_5_100kb_river <- 
  count_matrix_90_5_100kb_river[,colSums(count_matrix_90_5_100kb_river)!=0]
# DS diversity 132

#把count_matrix_90_5_100kb_river中大于0的值改为1。
count_matrix_90_5_100kb_river <- count_matrix_90_5_100kb_river %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, .)))
river_DS_diversity <- rowSums(count_matrix_90_5_100kb_river) %>% 
  as.data.frame() %>% setNames("DS_diversity")
river_DS_diversity$genome <- rownames(river_DS_diversity)
MAG_metadata_river_100kb <- merge(MAG_metadata_river_100kb, river_DS_diversity,
                                  by = "genome", all = T)
MAG_metadata_river_100kb[is.na(MAG_metadata_river_100kb[, 79]), 79] <- 0
MAG_metadata_river_100kb$DS_diversity_densitylog <- 
  log10(MAG_metadata_river_100kb$DS_diversity / MAG_metadata_river_100kb$Genome_Size *1000 +1)*1000

MAG_metadata_river_100kb_select <- select(MAG_metadata_river_100kb,
                                          genome, river_aquifer,DS_num,
                                          DS_density_kblog,DS_diversity,
                                          DS_diversity_densitylog)

#地下水defense####
gw_1626_DF <- read.delim2("1626mag_ds_counts_matrix.csv",sep = ',',check.names = F)
rownames(gw_1626_DF) <- gw_1626_DF$genome
genome_name <- gw_1626_DF$genome
gw_1626_DF <- gw_1626_DF[,-1]
gw_1626_DF <- as.data.frame(lapply(gw_1626_DF, as.numeric))
rownames(gw_1626_DF) <- genome_name
defense_num_1626 <- data.frame(DS_num = rowSums(gw_1626_DF),
                               genome = genome_name)


gw_1626_MAG_quality <- read.delim2("1626MAGs_90_5_N50100kb_score70.txt")

gw_1626_MAG_metadata <- merge(gw_1626_MAG_quality, defense_num_1626, 
                              by = "genome", all = T)
gw_1626_MAG_metadata[is.na(gw_1626_MAG_metadata[, 9]), 9] <- 0
gw_1626_MAG_metadata$DS_density_kb <- 
  gw_1626_MAG_metadata$DS_num / gw_1626_MAG_metadata$length *1000
gw_1626_MAG_metadata$DS_density_kblog <-
  log10(gw_1626_MAG_metadata$DS_density_kb + 1) *1000
gw_1626_MAG_metadata$river_aquifer <- "Aquifer"

gw_1626_DF <- gw_1626_DF[,colSums(gw_1626_DF)!=0]
# DS diversity 207 这可能是由于地下水MAG多导致的！

gw_1626_DF <- gw_1626_DF %>%
  mutate(across(everything(), ~ ifelse(. > 0, 1, .)))
gw_DS_diversity <- rowSums(gw_1626_DF) %>% 
  as.data.frame() %>% setNames("DS_diversity")
gw_DS_diversity$genome <- rownames(gw_DS_diversity)

gw_1626_MAG_metadata <- merge(gw_1626_MAG_metadata, gw_DS_diversity,
                              by = "genome", all = T)
gw_1626_MAG_metadata$DS_diversity_densitylog <- 
  log10(gw_1626_MAG_metadata$DS_diversity / gw_1626_MAG_metadata$length *1000 +1)*1000

gw_1626_MAG_metadata_select <- select(gw_1626_MAG_metadata, genome, river_aquifer,
                                      DS_num, DS_density_kblog, DS_diversity,
                                      DS_diversity_densitylog)

river_aquifer_metadata <- rbind(MAG_metadata_river_100kb_select,gw_1626_MAG_metadata_select)

library(ggplot2)
library(ggpubr)  # 用于添加统计检验
library(tidyr)   # 用于数据长宽转换

# 1. 数据准备（长格式转换）
df_long <- river_aquifer_metadata %>%
  pivot_longer(
    cols = c("DS_num", "DS_density_kblog", "DS_diversity", "DS_diversity_densitylog"),
    names_to = "Variable",
    values_to = "Value"
  )

df_long$Variable <- factor(df_long$Variable, levels = c("DS_num","DS_density_kblog",
                                                        "DS_diversity","DS_diversity_densitylog"))

# 2. 绘制箱线图+小提琴图
p <- ggplot(df_long, aes(x = river_aquifer, y = Value, fill = river_aquifer)) +
  # 先画数据点（底层）
  geom_jitter(
    width = 0.2,          # 避免点完全重叠
    alpha = 0.4,          # 透明度（0=完全透明，1=不透明）
    size = 1.5,           # 点大小
    color = "#ced4da"      # 点颜色（深灰色）
  ) +
  geom_violin(alpha = 0.5, trim = FALSE) +  # 小提琴图（透明度0.5）
  geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = NA) +  # 箱线图（隐藏异常点）
  facet_wrap(~Variable, scales = "free_y", ncol = 2) +  # 分面显示不同变量
  scale_fill_manual(values = c("River" = "#00bbf9", "Aquifer" = "#F28E2B")) +  # 自定义颜色
  labs(x = "Group", y = "Value", title = "Comparison of Defense System Metrics") +
  theme_bw() +
  theme(legend.position = "none")

# 3. 添加显著性检验（Wilcoxon秩和检验）
p <- p + stat_compare_means(
  method = "wilcox.test",  # 非参数检验
  comparisons = list(c("River", "Aquifer")),  # 比较组
  label = "p.signif",      # 显示显著性符号（***, **, *）
  symnum.args = list(      # 自定义符号阈值
    cutpoints = c(0, 0.001, 0.01, 0.05, 1),
    symbols = c("***", "**", "*", "ns")
  )
)

# 4. 显示图形
print(p)

setwd(file.path(GTR_DATA, "20gaoyuan/40plateau_groundwater_comparision"))
ggsave("1_river_aquifer_p.pdf",p,width = 6, height = 5.3, units = 'in')

















