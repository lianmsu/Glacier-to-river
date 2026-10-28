# =============================================================================
# Glacier-to-river -- defensome landscape / DS density vs genome size and elevation
# =============================================================================
# Relate DS number and DS density to genome size and sample elevation;
# builds the density comparison panels.
#
# Manuscript      : Fig. 2h-i, Fig. S2e, Fig. S4h-i, Fig. S5c-d
# Main outputs    : 04_density_and_genome_size.RData
# Original script : 20gaoyuan/10Defensefinder/4defensefinder_res_density_compare.R
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
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'Yes',"Plateau",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'No',"Plain",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'glacier',"Glacier",MAG_metadata$Plateau)
colnames(MAG_metadata)[17] <- 'Region'

MAG_metadata$medium <- MAG_metadata$Type
# MAG_metadata$medium <- ifelse(MAG_metadata$Type == 'Cryoconite',"Glacier",MAG_metadata$Type)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Ice',"Glacier",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Snow',"Glacier",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'particle',"Plankton",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'water',"Plankton",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'sediment',"Benthon",MAG_metadata$medium)


genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
count_matrix <- count_matrix[rownames(count_matrix) %in% MAG_metadata$genome,]


MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]


#不同介质中MAG数量统计#####################3
#所有样品来自的国家的统计
#MAG_metadata$medium <- factor(MAG_metadata$medium,levels = c("Benthon","Plankton","Glacier"))
#sample_info_medium <- MAG_metadata %>% count(medium, name = "Count")

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

all5455MAG_DSsum <- rowSums(count_matrix) %>% as.data.frame() %>% setNames("DS_num")
all5455MAG_DSsum$genome <- rownames(all5455MAG_DSsum)
all5455MAG_metadata <- left_join(MAG_metadata, 
                                 all5455MAG_DSsum, 
                                 by = c("genome" = "genome"))
all5455MAG_metadata$DS_num[is.na(all5455MAG_metadata$DS_num)] <- 0

all5455MAG_metadata <- merge(all5455MAG_metadata, MAG_quality, by.x = "genome",by.y = "Name")

all5455MAG_metadata$MAG_size <- all5455MAG_metadata$Genome_Size/1000000
all5455MAG_metadata$DS_densityperkb <- all5455MAG_metadata$DS_num/all5455MAG_metadata$Genome_Size*1000

all5455MAG_metadata$DS_densitylog <- log10(all5455MAG_metadata$DS_densityperkb+1)*1000


all5455MAG_metadata_extract <- select(all5455MAG_metadata,genome,ecosystem,medium,Region,DS_densitylog)

all5455MAG_metadata_extract$ecosystem <- factor(all5455MAG_metadata_extract$ecosystem,
                                                levels = c("Glacier","HH","CJ","LCJ","YJ"))
all5455MAG_metadata_extract$medium <- factor(all5455MAG_metadata_extract$medium,
                                                levels = c("Snow","Ice","Cryoconite","water","particle","sediment"))
all5455MAG_metadata_extract$Region <- factor(all5455MAG_metadata_extract$Region,
                                             levels = c("Glacier","Plateau","Plain"))


#5. 绘制箱线+小提琴+散点复合图
library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)

#计算显著性差异
library(rstatix)
md_p_val1 <- all5455MAG_metadata_extract %>% 
  group_by(ecosystem) %>%
  wilcox_test(DS_densitylog ~ medium) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="medium",dodge=0.5,group = "ecosystem")

# ds_density$log10density <- log10(ds_density$density+1)
# md_p_val2 <- ds_density %>% 
#   wilcox_test(log10density ~ region) %>%
#   adjust_pvalue(p.col="p",method="bonferroni") %>%
#   add_significance(p.col="p.adj") %>% 
#   add_xy_position(x="region",dodge=0.5)

#md_p_val1$y.position <- seq(from = 0.6, to = 0.6 + (0.03 * (21 - 1)), by = 0.03)
# 调整y轴位置
#md_p_val1$y.position <- max(all5455MAG_metadata_extract$DS_densitylog, na.rm = TRUE) * 
#  (1 + 0.05 * (1:nrow(md_p_val1)))

md_p_val1$y.position <- md_p_val1$y.position-2.5

library(ggplot2)
library(gghalves)
library(cowplot)

# 定义 6 种 medium 的固定颜色
medium_colors <- c("Snow" = "#4DBBD5",
                   "Ice" = "#E64B35",
                   "Cryoconite" = "#00A087",
                   "water" = "#3C5488",
                   "particle" = "#F39B7F",
                   "sediment" = "#8491B4")

# 各 ecosystem 实际存在的 medium 组合
ecosystem_media <- all5455MAG_metadata_extract |>
  group_by(ecosystem, medium) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 ecosystem 独立绘图
plot_list <- list()
ecosystems <- levels(all5455MAG_metadata_extract$ecosystem)
n_media <- c()

for (i in seq_along(ecosystems)) {
  eco <- ecosystems[i]
  
  existing_media <- ecosystem_media |>
    filter(ecosystem == eco) |>
    pull(medium)
  n_media <- c(n_media, length(existing_media))
  
  df_sub <- all5455MAG_metadata_extract |>
    filter(ecosystem == eco) |>
    mutate(medium = factor(medium, levels = existing_media))
  
  # 重新计算显著性
  if (length(existing_media) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(DS_densitylog ~ medium) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "medium", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 2.5
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = medium, y = DS_densitylog)) +
    geom_half_boxplot(aes(fill = medium), color = "black",
                      side = "l", errorbar.draw = T,
                      outlier.shape = NA, width = 0.8) +
    geom_jitter(aes(color = medium), width = 0.1, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = medium_colors) +
    scale_color_manual(values = medium_colors) +
    labs(x = NULL, y = NULL, title = eco) +
    theme_bw() +
    theme(plot.margin = margin(5, 5, 5, 5),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.line = element_line(color = "black", linewidth = 0.4),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(margin = margin(t = 2), color = "black",
                                     size = 10, angle = 45, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 只在第一张图保留 y 轴
  if (i > 1) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.y = element_blank())
  }
  
  if (!is.null(pval_sub)) {
    p <- p + stat_pvalue_manual(pval_sub, label = "p.adj.signif",
                                label.size = 5, hide.ns = T)
  }
  
  plot_list[[eco]] <- p
}

# 合并：各 panel 宽度与 media 数量成正比
p_ecosystem_medium <- plot_grid(plotlist = plot_list, nrow = 1,
                                 rel_widths = n_media)
p_ecosystem_medium



# 定义 ecosystem 固定颜色（保留原配色方案）
ecosystem_colors <- c("YJ" = "#FCB2AF",
                      "Glacier" = "#9BDFDF",
                      "LCJ" = "#FFE2CE",
                      "CJ" = "#BEBCDF",
                      "HH" = "#FB8C62")

# 各 Region 实际存在的 ecosystem 组合（排除 Glacier）
region_ecosystem_media <- all5455MAG_metadata_extract |>
  filter(ecosystem != "Glacier") |>
  group_by(Region, ecosystem) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 Region 独立绘图
plot_list_region <- list()
regions <- levels(droplevels(all5455MAG_metadata_extract$Region[
  all5455MAG_metadata_extract$Region != "Glacier"]))
n_eco <- c()

for (i in seq_along(regions)) {
  reg <- regions[i]
  
  existing_eco <- region_ecosystem_media |>
    filter(Region == reg) |>
    pull(ecosystem)
  n_eco <- c(n_eco, length(existing_eco))
  
  df_sub <- all5455MAG_metadata_extract |>
    filter(ecosystem != "Glacier", Region == reg) |>
    mutate(ecosystem = factor(ecosystem, levels = existing_eco))
  
  # 重新计算显著性
  if (length(existing_eco) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(DS_densitylog ~ ecosystem) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "ecosystem", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 2.5
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = ecosystem, y = DS_densitylog)) +
    geom_half_boxplot(aes(fill = ecosystem), color = "black",
                      side = "l", errorbar.draw = T,
                      outlier.shape = NA, width = 0.8) +
    geom_jitter(aes(color = ecosystem), width = 0.1, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = ecosystem_colors) +
    scale_color_manual(values = ecosystem_colors) +
    labs(x = NULL, y = NULL, title = reg) +
    theme_bw() +
    theme(plot.margin = margin(5, 5, 5, 5),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.line = element_line(color = "black", linewidth = 0.4),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(margin = margin(t = 2), color = "black",
                                     size = 10, angle = 45, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 只在第一张图保留 y 轴
  if (i > 1) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.y = element_blank())
  }
  
  if (!is.null(pval_sub)) {
    p <- p + stat_pvalue_manual(pval_sub, label = "p.adj.signif",
                                label.size = 5, hide.ns = T)
  }
  
  plot_list_region[[reg]] <- p
}

# 合并：各 panel 宽度与 ecosystem 数量成正比
p_Region_ecosystem <- plot_grid(plotlist = plot_list_region, nrow = 1,
                                 rel_widths = n_eco)
p_Region_ecosystem


###海拔与DS的关系##################

# 创建一个新的分组变量：每500单位一个Elevation组
all5455MAG_metadata_elevation <- all5455MAG_metadata %>%
  mutate(Elevation_Group = cut(Elevation, 
                               breaks = seq(floor(min(Elevation, na.rm = TRUE)/500)*500, 
                                            ceiling(max(Elevation, na.rm = TRUE)/500)*500, 
                                            by = 500),
                               include.lowest = TRUE,
                               right = FALSE))

# 计算每组的平均值和标准差
summary_stats <- all5455MAG_metadata_elevation %>%
  group_by(Elevation_Group) %>%
  summarise(mean_DS_densitylog = mean(DS_densitylog, na.rm = TRUE),
            sd_DS_densitylog = sd(DS_densitylog, na.rm = TRUE),
            .groups = 'drop')


# 进行显著性检验：Elevation与DS_densitylog的相关性
cor_result <- cor.test(all5455MAG_metadata_elevation$Elevation, 
                       all5455MAG_metadata_elevation$DS_densitylog, 
                       alternative = "two.sided", 
                       method = "spearman", 
                       conf.level = 0.95)

# 计算每组的平均值和标准差
# 计算每组的平均值和标准差，并提取Elevation中位值（每组中Elevation的真实中位数）
summary_stats <- all5455MAG_metadata_elevation %>%
  group_by(Elevation_Group) %>%
  summarise(mean_DS_densitylog = mean(DS_densitylog, na.rm = TRUE),
            sd_DS_densitylog = sd(DS_densitylog, na.rm = TRUE),
            Elevation_mid = median(Elevation, na.rm = TRUE),
            .groups = 'drop')

# 绘图并添加线性拟合线
p_elevation_DSdensity <- 
  ggplot(summary_stats, aes(x = Elevation_mid, y = mean_DS_densitylog)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_DS_densitylog - sd_DS_densitylog,
                    ymax = mean_DS_densitylog + sd_DS_densitylog), 
                width = 200) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(x = "Elevation (bin midpoints, 500 unit bins)",
       y = "log10 (DS Density +1) × 1e-3",
       title = "DS Density (log) across Elevation Groups with Linear Fit") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))+
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
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "right",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

print(p_elevation_DSdensity)

write.table(summary_stats, "summary_stats.tsv",sep = "\t",quote = F, row.names = F)

# 进行显著性检验：Elevation与DS_densitylog的相关性（整体数据）
cor_result <- cor.test(all5455MAG_metadata$Elevation, 
                       all5455MAG_metadata$DS_densitylog, 
                       alternative = "two.sided", 
                       method = "spearman", 
                       conf.level = 0.95)

print(cor_result)

# 使用分组中位数和均值进行线性拟合相关性检验
lm_result <- lm(mean_DS_densitylog ~ Elevation_mid, data = summary_stats)
summary(lm_result)



library(Cairo)
library(cowplot)
cairo_pdf(filename = "4_geografical_DSdensity1.pdf" , width = 6.6,height = 6.9)
g1 <- plot_grid(p_Region_ecosystem,
               p_elevation_DSdensity,
               labels = letters[2:3], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

g2 <- plot_grid(p_ecosystem_medium,g1,
                labels = letters[1], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2


# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()




#河流内部 高原与平原段比较########
all5455MAG_metadata_extract <- select(all5455MAG_metadata,genome,ecosystem,medium,Region,DS_densitylog)

all5455MAG_metadata_extract$ecosystem <- factor(all5455MAG_metadata_extract$ecosystem,
                                                levels = c("Glacier","HH","CJ","LCJ","YJ"))
all5455MAG_metadata_extract$medium <- factor(all5455MAG_metadata_extract$medium,
                                             levels = c("Snow","Ice","Cryoconite","water","particle","sediment"))

all5455MAG_metadata_extract$Region <- ifelse(all5455MAG_metadata_extract$Region == "Glacier", "Plateau", all5455MAG_metadata_extract$Region)

all5455MAG_metadata_extract$Region <- factor(all5455MAG_metadata_extract$Region,
                                             levels = c("Plateau","Plain"))

# 定义 Region 固定颜色
region_colors <- c("Plateau" = "#C4D8E9", "Plain" = "#d4e09b")

# 各 ecosystem 实际存在的 Region 组合（排除 Glacier 和 YJ）
eco_region_data <- all5455MAG_metadata_extract |>
  group_by(ecosystem, Region) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 ecosystem 独立绘图
plot_list_eco_region <- list()
ecos_for_region <- eco_region_data |> pull(ecosystem) |> unique()
n_region_per_eco <- c()

for (i in seq_along(ecos_for_region)) {
  eco_val <- ecos_for_region[i]
  
  existing_regions <- eco_region_data |>
    filter(ecosystem == eco_val) |>
    pull(Region)
  n_region_per_eco <- c(n_region_per_eco, length(existing_regions))
  
  df_sub <- all5455MAG_metadata_extract |>
    filter(ecosystem == eco_val) |>
    mutate(Region = factor(Region, levels = existing_regions))
  
  # 重新计算显著性
  if (length(existing_regions) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(DS_densitylog ~ Region) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "Region", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 2.5
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = Region, y = DS_densitylog)) +
    geom_half_boxplot(aes(fill = Region), color = "black",
                      side = "l", errorbar.draw = T,
                      outlier.shape = NA, width = 0.8) +
    geom_jitter(aes(color = Region), width = 0.1, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = region_colors) +
    scale_color_manual(values = region_colors) +
    labs(x = NULL, y = NULL, title = eco_val) +
    theme_bw() +
    theme(plot.margin = margin(5, 5, 5, 5),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.line = element_line(color = "black", linewidth = 0.4),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(margin = margin(t = 2), color = "black",
                                     size = 10, angle = 45, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 只在第一张图保留 y 轴
  if (i > 1) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.y = element_blank())
  }
  
  if (!is.null(pval_sub)) {
    p <- p + stat_pvalue_manual(pval_sub, label = "p.adj.signif",
                                label.size = 5, hide.ns = T)
  }
  
  plot_list_eco_region[[eco_val]] <- p
}

# 合并：各 panel 宽度与 Region 数量成正比
p_ecosystem_Region <- plot_grid(plotlist = plot_list_eco_region, nrow = 1,
                                 rel_widths = n_region_per_eco)
p_ecosystem_Region


library(Cairo)
library(cowplot)
cairo_pdf(filename = "4_geografical_DSdensity2.pdf" , width = 6.6,height = 6.9)
g1 <- plot_grid(p_ecosystem_Region,
                p_elevation_DSdensity,
                labels = letters[2:3], ncol = 2,
                rel_widths = c(3, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

g2 <- plot_grid(p_ecosystem_medium,g1,
                labels = letters[1], ncol = 1
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2


setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

##以Region分组，比较不同medium########
#计算显著性差异

library(rstatix)

# 各 Region 实际存在的 medium 组合
region_media <- all5455MAG_metadata_extract |>
  group_by(Region, medium) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 Region 独立绘图
plot_list_region_med <- list()
regions_med <- region_media |> pull(Region) |> unique()
n_med_per_region <- c()

for (i in seq_along(regions_med)) {
  reg_val <- regions_med[i]
  
  existing_media <- region_media |>
    filter(Region == reg_val) |>
    pull(medium)
  n_med_per_region <- c(n_med_per_region, length(existing_media))
  
  df_sub <- all5455MAG_metadata_extract |>
    filter(Region == reg_val) |>
    mutate(medium = factor(medium, levels = existing_media))
  
  # 重新计算显著性
  if (length(existing_media) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(DS_densitylog ~ medium) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "medium", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 5
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = medium, y = DS_densitylog)) +
    geom_half_boxplot(aes(fill = medium), color = "black",
                      side = "l", errorbar.draw = T,
                      outlier.shape = NA, width = 0.8) +
    geom_jitter(aes(color = medium), width = 0.1, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = medium_colors) +
    scale_color_manual(values = medium_colors) +
    labs(x = NULL, y = NULL, title = reg_val) +
    theme_bw() +
    theme(plot.margin = margin(5, 5, 5, 5),
          panel.background = element_blank(),
          plot.background = element_blank(),
          strip.text = element_text(size = 12),
          axis.line = element_line(color = "black", linewidth = 0.4),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_text(margin = margin(t = 2), color = "black",
                                     size = 10, angle = 45, vjust = 1, hjust = 1),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 只在第一张图保留 y 轴
  if (i > 1) {
    p <- p + theme(axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.line.y = element_blank())
  }
  
  if (!is.null(pval_sub)) {
    p <- p + stat_pvalue_manual(pval_sub, label = "p.adj.signif",
                                label.size = 5, hide.ns = T)
  }
  
  plot_list_region_med[[reg_val]] <- p
}

# 合并：各 panel 宽度与 medium 数量成正比
p_Region_medium <- plot_grid(plotlist = plot_list_region_med, nrow = 1,
                              rel_widths = n_med_per_region)
p_Region_medium


library(Cairo)
library(cowplot)
cairo_pdf(filename = "4_geografical_DSdensity5_new.pdf" , width = 8,height = 8)
g1 <- plot_grid(p_ecosystem_Region,
                p_elevation_DSdensity,
                labels = letters[4:5], ncol = 2,
                rel_widths = c(3, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

g2 <- plot_grid(p_Region_ecosystem,p_Region_medium,
                labels = letters[2:3], ncol = 2,
                align = c("v"), axis = "lr"
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2



g3 <- plot_grid(p_ecosystem_medium,g2,g1,
                labels = letters[1], ncol = 1
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g3



dev.off()




# p_Region_ecosystem2 <- ggplot(all5455MAG_metadata_extract[all5455MAG_metadata_extract$ecosystem != "Glacier",], 
#                              aes(x = ecosystem, y = DS_densitylog)) +
#   geom_half_boxplot(aes(fill=ecosystem), color ="black",
#                     side = "l",
#                     errorbar.draw = T,
#                     outlier.shape = NA,
#                     width=0.8
#   )+
#   geom_jitter(aes(color=ecosystem), width = 0.1, size = 0.5, alpha = 0.6)+
#   stat_pvalue_manual(md_p_val2,label = "p.adj.signif",label.size=5,hide.ns = T)+
#   facet_wrap(.~Region,nrow=1)+ #分⾯操作
#   # 定义颜⾊
#   # scale_fill_npg()+
#   # scale_color_npg()+
#   scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
#                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
#                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
#   )) +  # 手动设置颜色
#   scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
#                                 "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
#                                 "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
#   )) +  # 手动设置颜色
#   labs(x="River",y="log10 (DS Density +1) × 1e-3")+
#   #修改主题
#   theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
#         panel.background = element_blank(),
#         plot.background = element_blank(),
#         strip.text = element_text(size=12),
#         axis.line = element_line(color = "black",size = 0.4),
#         axis.text.y = element_text(color="black",size=10),
#         # axis.text.x = element_text(margin = margin(t = 2),color="black",
#         #                            size=10,angle = 45,vjust = 1,hjust = 1),
#         legend.position = "none",
#         axis.title.x = element_blank(),  # 隐藏x轴标题
#         axis.text.x = element_blank(),   # 隐藏x轴标签（坐标值）
#         axis.ticks.x = element_blank(),   # 可选：隐藏x轴刻度线
#         panel.spacing = unit(0,"lines"))+
#   coord_cartesian()
# p_Region_ecosystem2

# p_Region_medium2 <- ggplot(all5455MAG_metadata_extract, 
#                           aes(x = medium, y = DS_densitylog)) +
#   geom_half_boxplot(aes(fill=medium), color ="black",
#                     side = "l",
#                     errorbar.draw = T,
#                     outlier.shape = NA,
#                     width=0.8
#   )+
#   geom_jitter(aes(color=medium), width = 0.1, size = 0.5, alpha = 0.6)+
#   stat_pvalue_manual(md_p_val4,label = "p.adj.signif",label.size=5,hide.ns = T)+
#   facet_wrap(.~Region,nrow=1)+ #分⾯操作
#   # 定义颜⾊
#   scale_fill_npg()+
#   scale_color_npg()+
#   # scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0","Plateau"= "#C4D8E9",
#   #                              "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70","Plain"="#d4e09b",
#   #                              "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
#   # )) +  # 手动设置颜色
#   # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0","Plateau"= "#C4D8E9",
#   #                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70","Plain"="#d4e09b",
#   #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
#   # )) +  # 手动设置颜色
#   labs(x="medium",y="log10 (DS Density +1) × 1e-3")+
#   #修改主题
#   theme(plot.margin=unit(c(0.5,0.5,0.5,0.5),units=,"cm"),
#         panel.background = element_blank(),
#         plot.background = element_blank(),
#         strip.text = element_text(size=12),
#         axis.line = element_line(color = "black",size = 0.4),
#         axis.text.y = element_text(color="black",size=10),
#         # axis.text.x = element_text(margin = margin(t = 2),color="black",
#         #                            size=10,angle = 45,vjust = 1,hjust = 1),
#         legend.position = "none",
#         axis.title.x = element_blank(),  # 隐藏x轴标题
#         axis.text.x = element_blank(),   # 隐藏x轴标签（坐标值）
#         axis.ticks.x = element_blank(),   # 可选：隐藏x轴刻度线
#         panel.spacing = unit(0,"lines"))+
#   coord_cartesian()
# p_Region_medium2

g2 <- plot_grid(p_Region_ecosystem2,p_Region_medium2,
                labels = letters[2:3], ncol = 2
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g2

cairo_pdf(filename = "4_region2ecosystem_medium.pdf" , width = 7.5,height = 2.4)
g2
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

write.table(df_sub, "df_sub_ds_density.tsv",sep = "\t", quote = F, row.names = F)

#直接比较glacier、plateau、plain；四条河流；Benthon和Plankton ####

genome_defense_num <- genome_defense_num[genome_defense_num$genome %in% MAG_metadata$genome,]
MAG_metadata <- merge(MAG_metadata, genome_defense_num, by = "genome", all = T)
MAG_metadata[,63][is.na(MAG_metadata[,63])] <- 0

MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Cryoconite',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Ice',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Snow',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'particle',"Plankton",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'water',"Plankton",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'sediment',"Benthon",MAG_metadata$medium)
MAG_metadata$Region <- factor(MAG_metadata$Region,
                              levels = c("Glacier","Plateau","Plain"))
MAG_metadata$medium <- factor(MAG_metadata$medium,
                              levels = c("Glacier","Benthon","Plankton"))
MAG_metadata$ecosystem <- factor(MAG_metadata$ecosystem,
                              levels = c("Glacier","YJ","LCJ","HH","CJ"))

MAG_metadata <- merge(MAG_metadata, MAG_quality, by.x = "genome",by.y = "Name", all = T)
MAG_metadata$MAG_size <- MAG_metadata$Genome_Size/1000000
MAG_metadata$DS_densityperkb <- MAG_metadata$DS_num/MAG_metadata$Genome_Size*1000
MAG_metadata$DS_densitylog <- log10(MAG_metadata$DS_densityperkb+1)*1000

library(rstatix)
md_p_val_DS1 <- MAG_metadata %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DS_densitylog ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val_DS1$y.position <- md_p_val_DS1$y.position - 3

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)
p_DSdensity_glacier_plateau_plain <- 
  ggplot(MAG_metadata, 
         aes(x = Region, y = DS_densitylog)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color=Region), width = 0.1, size = 0.5, alpha = 0.6)+
  stat_pvalue_manual(md_p_val_DS1,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="log10 (DS Density +1) × 1e-3")+
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
  coord_cartesian(ylim = c(-0.1, 5))
p_DSdensity_glacier_plateau_plain


md_p_val_DS2 <- MAG_metadata %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DS_densitylog ~ ecosystem) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="ecosystem",dodge=0.5)
md_p_val_DS2$y.position <- md_p_val_DS2$y.position - 4

p_DSdensity_ecossytem <- 
  ggplot(MAG_metadata[MAG_metadata$medium != "Glacier",], 
         aes(x = ecosystem, y = DS_densitylog)) +
  geom_half_boxplot(aes(fill = ecosystem), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color=ecosystem), width = 0.1, size = 0.5, alpha = 0.6)+
  #stat_pvalue_manual(md_p_val_DS2,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="log10 (DS Density +1) × 1e-3")+
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
  coord_cartesian(ylim = c(-0.1, 5))
p_DSdensity_ecossytem


md_p_val_DS3 <- MAG_metadata %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DS_densitylog ~ medium) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="medium",dodge=0.5)
md_p_val_DS3$y.position <- md_p_val_DS3$y.position - 4

p_DSdensity_medium <- 
  ggplot(MAG_metadata[MAG_metadata$medium != "Glacier",], 
         aes(x = medium, y = DS_densitylog)) +
  geom_half_boxplot(aes(fill = medium), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color=medium), width = 0.1, size = 0.5, alpha = 0.6)+
  #stat_pvalue_manual(md_p_val_DS3,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="log10 (DS Density +1) × 1e-3")+
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
  coord_cartesian(ylim = c(-0.1, 5))
p_DSdensity_medium

library(Cairo)
library(cowplot)
g4 <- plot_grid(p_DSdensity_glacier_plateau_plain,p_DSdensity_ecossytem,p_DSdensity_medium,
                labels = letters[2:3], ncol = 3
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g4

cairo_pdf(filename = "4_DSdensity_comparision_new.pdf" , width = 8,height = 2.4)
g4
# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()


save.image("./04_density_and_genome_size.RData")

load("./04_density_and_genome_size.RData")

write.table(MAG_metadata, "MAG_metadata.tsv",
    quote = F, sep = "\t", row.names = F)
