# =============================================================================
# Glacier-to-river -- defense islands / DefenseFinder vs DefensePredictor
# =============================================================================
# Compare the defense genes reported by DefenseFinder and DefensePredictor
# across the 5,455 MAGs (counts, density, overlap).
#
# Manuscript      : Fig. S13
# Main outputs    : FigS6B_C_density.pdf, 1_DP_DF_compare_3.pdf
# Original script : 20gaoyuan/16DefensePredictor/1DF_DP.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

#分析defensepredictor和defensefinder的结果。

base_dir  <- file.path(GTR_DATA, "20gaoyuan")
proj_dir  <- file.path(base_dir, "16DefensePredictor")
meta_path <- file.path(base_dir, "01MAG_info/5455MAG_metadata.csv")

setwd(proj_dir)

library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)

# 5455 个 MAG 的清单（基准，所有结果都限定在这 5455 个 MAG 内比较）
MAG5455_metadata <- read.csv(meta_path)
mag_genomes <- MAG5455_metadata$genome   # 5455 个 genome

# DefenseFinder 基因的 activity 标注（Defense / Antidefense），用于排除 Anti-defense 基因。
# 以 (hit_id, sys_id) 作为唯一键，与 DF_DG_position.tsv 的行对齐。
DF_activity <- read.delim("merged_5467MAG_defense_finder_genes.tsv") %>%
  select(hit_id, sys_id, activity)

#1. 密度图##########################

#1.1
#Defense predictor结果（每个 genome 的 defense gene 数量）
DP <- read.csv("defensepredictor_genecount.csv", sep = ',', header = F)
colnames(DP) <- c("genome", "DP_number")
DP <- DP[DP$genome %in% mag_genomes, ]

#Defense Finder结果（每个 genome 的 DG 数量），排除 Anti-defense 基因并去冗余
DF <- read.csv("DF_DG_position.tsv", sep = '\t')
DF <- DF %>%
  mutate(genome = str_extract(Contig, ".*bin\\.\\d+")) %>%
  left_join(DF_activity, by = c("HitID" = "hit_id", "SysID" = "sys_id")) %>%
  filter(activity == "Defense") %>%
  distinct(HitID, .keep_all = TRUE)   # 同一基因因属于多个系统而重复列出，去冗余
DF <- DF[DF$genome %in% mag_genomes, ]

genome_DG_counts <- as.data.frame(table(DF$genome))
colnames(genome_DG_counts) <- c("genome", "DG_number")

# 以 5455 个 MAG 为基准，未检出的 genome 数量记为 0
DP_DF_gene_num <- data.frame(genome = mag_genomes, stringsAsFactors = FALSE)
DP_DF_gene_num <- merge(DP_DF_gene_num, DP, by = "genome", all.x = TRUE)
DP_DF_gene_num <- merge(DP_DF_gene_num, genome_DG_counts, by = "genome", all.x = TRUE)
DP_DF_gene_num[is.na(DP_DF_gene_num)] <- 0

DP_DF_gene_num <- DP_DF_gene_num %>% arrange(DG_number)

mytheme <-  theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( lineheight = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(lineheight = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( lineheight = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(lineheight = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

# 绘制点线图
p1 <- 
ggplot(DP_DF_gene_num, aes(x = factor(genome, levels = genome))) +
  geom_point(aes(y = DP_number, color = "DP_number"), size = 1) +  # 绘制 DP_number 点
#  geom_line(aes(y = DP_number, color = "DP_number")) +   # 绘制 DP_number 线
  geom_point(aes(y = DG_number, color = "DG_number"), size = 1) +  # 绘制 DG_number 点
#  geom_line(aes(y = DG_number, color = "DG_number")) +   # 绘制 DG_number 线
  labs(x = "MAGs", y = "Number of defense genes", color = "Legend") +  # 设置坐标轴和图例标题
  scale_color_manual(values = c("DP_number" = "#efc3e6", "DG_number" = "#b8bedd")) +  # 设置颜色
  mytheme +
  theme(axis.text.x = element_blank(),  # 隐藏横坐标的具体内容
        axis.ticks.x = element_blank())   # 隐藏横坐标的刻度
p1
write.csv(DP_DF_gene_num, "FigS6B.csv")

#密度图
DP_DF_gene_num2 <- DP_DF_gene_num[, -1]

DP_DF_gene_num2 %>% gather(key = "key", value = "val", 1:2) -> data

p2 <- 
ggplot(data, aes(x = val, fill = key, color = key)) +
  geom_density(alpha = 0.25) +
  scale_color_manual(values = c("DP_number" = "#efc3e6", "DG_number" = "#b8bedd")) +
  scale_fill_manual(values = c("DP_number" = "#efc3e6", "DG_number" = "#b8bedd")) +
  labs(x = "Number of defense genes", y = "Density of MAGs") +  # 设置坐标轴标题
  mytheme
p2
write.csv(data, "FigS6C.csv")

subset_data1 <- data[data$key == 'DG_number', ]
density_data1 <- density(subset_data1$val)
peak_x1 <- density_data1$x[which.max(density_data1$y)]
peak_y1 <- max(density_data1$y)

p2 <- p2 +
  geom_vline(xintercept = peak_x1, linetype = "dashed", color = "#f07167") +
  geom_text(aes(x = peak_x1, y = 0, label = round(peak_x1, 2)), 
            vjust = 0, hjust = 0, color = "black")

p2

subset_data2 <- data[data$key == 'DP_number', ]
density_data2 <- density(subset_data2$val)
peak_x2 <- density_data2$x[which.max(density_data2$y)]
peak_y2 <- max(density_data2$y)

p2 <- p2 +
  geom_vline(xintercept = peak_x2, linetype = "dashed", color = "#1BC6CA") +
  geom_text(aes(x = peak_x2, y = 0, label = round(peak_x2, 2)), 
            vjust = 0, hjust = 0, color = "black")

p2


library(Cairo)
library(cowplot)
cairo_pdf("FigS6B_C_density.pdf", width = 6, height = 8)
g <- plot_grid(p1, p2,
               labels = LETTERS[1:2], nrow = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
dev.off()


#2. 统计DP与DF重合的基因########################################################

DG_coords <- read.delim2("DF_DG_position.tsv", sep = '\t',
                        header = T)
DG_coords <- DG_coords %>%
  mutate(
    genome = str_extract(Contig, ".*bin\\.\\d+")  # 包含"bin.数字"
  ) %>%
  left_join(DF_activity, by = c("HitID" = "hit_id", "SysID" = "sys_id")) %>%
  filter(activity == "Defense") %>%    # 排除 Anti-defense 基因
  distinct(HitID, .keep_all = TRUE)    # 同一基因因属于多个系统而重复列出，去冗余
DG_coords <- DG_coords[DG_coords$genome %in% mag_genomes, ]
colnames(DG_coords)[3:4] <- c("start", "end")

DP_gene_info <- read.csv("5466mag_defensepredictor_processed_filtered.csv",
                           sep = ',',
                           header = T)
DP_gene_info <- DP_gene_info %>%
  mutate(
    genome = str_extract(product_accession, ".*bin\\.\\d+")  # 包含"bin.数字"
  )
DP_gene_info <- DP_gene_info[DP_gene_info$genome %in% mag_genomes, ]


# 统计 genome、start、end 三者均一致的重合基因数。
# 只按 start/end 匹配会跨基因组误配（不同 MAG 内坐标相同的基因被算成重合），
# 这里限定同一 genome 后再按坐标匹配，并以唯一坐标计数，避免多对多放大。
matched_rows <- inner_join(DG_coords, DP_gene_info, 
                           by = c("genome", "start", "end"),
                           relationship = "many-to-many")
num_matched_rows <- nrow(distinct(matched_rows, genome, start, end))

# 输出结果
print(paste("Number of overlapping genes (matching genome, start and end):", num_matched_rows))

# 直接从 5455 个 MAG 过滤后的数据重新计算，不再硬编码
chonghe <- num_matched_rows
chonghe_data <- data.frame(type = c("DF", "DP"),
                           gene_num = c(nrow(DG_coords), nrow(DP_gene_info)),
                           chonghe = c(chonghe, chonghe)
                           )
#基因数量柱状图，重合比例图
# 自定义颜色和主题
p3 <- 
ggplot(chonghe_data, aes(x = type, y = gene_num, fill = type)) +
  geom_bar(stat = "identity", width = 0.6) +
  labs(x = "Software", y = "Number of defense genes", fill = "Type") +
  scale_fill_manual(values = c("DF" = "#b8bedd", "DP" = "#efc3e6")) +  # 自定义颜色
  theme_minimal() +
  mytheme
p3

# 计算 chonghe 占 gene_num 的百分比
chonghe_data$both <- (chonghe_data$chonghe / chonghe_data$gene_num) * 100
chonghe_data$single <- 100 - chonghe_data$both
chonghe_data2 <- chonghe_data[, c(1, 4, 5)]
library(reshape)
chonghe_data2_long <- melt(chonghe_data2)

colorsss <- c(
  "#D1CB7E",  # 橙色
  "#EFEFEF",  # 青色
  # "#F5E320",  # 红色
  # "#CC281C",  # 绿色
  "#D1CB7E",  # 紫色
  "#5AAAB2",  # 棕色
  "#E77F5A",  # 粉色
  "#A0E6C5",  # 灰色
  "#70A47F",  # 橄榄色
  "#1f77b4",  # 蓝色
  "#D68112",  # 浅蓝色
  "#C9A79E",  # 浅橙色
  "#98df8a",  # 浅绿色
  #  "#ff9896",  # 浅红色
  "#6C6B6F",  # 浅紫色
  "#CED2DD",  # 浅棕色
  #  "#f7b6d2",  # 浅粉色
  "#4E616A",  # 浅橄榄色
  #  "#9edae5",  # 浅青色
  "#ED6764",  # 深蓝色
  "#74A0A4"   # 深绿色
)
p4 <-  
  ggplot(chonghe_data2_long, aes(x = type, y = value, fill = variable)) +
  geom_bar(stat = "identity", width = 0.6) +  # 使用实际值绘制柱状图
  scale_y_continuous(limits = c(0, 100)) +  # 设置 y 轴的范围和刻度
  #scale_y_continuous(limits=c(0,100)) +  # 将y轴标签转换为百分比
  scale_fill_manual(values = colorsss) +
  # scale_x_discrete(limits=c('Cd11','Cr11','Cu11','phOH11','PNP11','PAP11','CK11',
  #                           'Cd21','Cr21','Cu21','phOH21','PNP21','PAP21','CK21',
  #                           'Seed'))+
  #coord_flip() +  # 翻转坐标轴
  #scale_y_continuous(limits=c(0,100),expand = c(0,0)) +  # 将y轴标签转换为百分比
  labs(x = "Software", y = "Relative proportion (%)", fill = "") +  # 设置轴标签和图例标题
  
  #  geom_vline(xintercept = c(7.5,14.5),linetype=2)+
  guides(fill = guide_legend(title = "", nrow = 1, position = 'bottom')) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.2),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    #legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p4
write.csv(chonghe_data, "FigS6AD.csv")
library(Cairo)
library(cowplot)
cairo_pdf("1_DP_DF_compare_3.pdf", width = 5.64, height = 5.64)

g <- plot_grid(p3, p2, p1, p4,
               labels = letters[1:4], nrow = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g

setwd(proj_dir)
dev.off()



#3比较Defensefinder的HMM结果与DG结果，或者加上DP，说明DG是合理的##########################

#利用DP来注释防御岛上未被注释为DG的基因，绘制寻找新DG的示例#############
