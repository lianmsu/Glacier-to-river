# =============================================================================
# Glacier-to-river -- defensome landscape / Defensome patterns across phyla
# =============================================================================
# Per-phylum DS number, density and diversity, plus defense-system frequency
# for Cyanobacteriota, Myxococcota and Gemmatimonadota.
#
# Manuscript      : Fig. S8, Fig. S9
# Main outputs    : 6_phyla_dss_frequency.pdf, 6_Cyanobacteriota_order_dss_frequency.pdf,
#                   6_Myxococcota_class_DS_frequency.pdf,
#                   6_Gemmatimonadota_phylum_DS_frequency.pdf
# Original script : 20gaoyuan/10Defensefinder/6_specialtax_DS.R
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

DF_res <- read.delim2("merged_5467MAG_defense_finder_systems.tsv",
                      check.names = F,
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


#1、主要门中的主要DS的frequency####
dss <- count_matrix

not_in_dss <- MAG_metadata$genome[!MAG_metadata$genome %in% rownames(dss)]
new_rows <- data.frame(
  genome = not_in_dss,  # 重复 genome 值
  stringsAsFactors = FALSE
)
# 为 new_rows 添加与 dss 相同的列，并将这些列的值填充为 0
col_names <- setdiff(names(dss), "genome")  # 获取 dss 中除了 genome 之外的列名
new_rows[col_names] <- lapply(col_names, function(col_name) rep(0, nrow(dss)))
rownames(new_rows) <- new_rows$genome
new_rows <- new_rows[,-1]
#合并new_rows与dss
dss <- rbind(dss,new_rows)

dss$genome <- rownames(dss)
dss <- select(dss, genome, everything())
dss <- merge(dss, mag_lineage_split, by = "genome")

library(dplyr)
#对tax_dss的col2列进行处理，删除_后面的字母
# 假设 tax_dss 是你的数据框，col2 是你想要处理的列
library(stringr)
# 使用 sub 函数替换掉第三个 _ 及其后面的内容
dss$phylum <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", dss$phylum)

#统计col2列每个元素的数量
element_counts <- dss %>%
  count(phylum)

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(phylum) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE) * 100))
  return(output)
}
result_data <- result(dss) %>% as.data.frame()

#过滤数据
library(dplyr)
#筛选>100的mags
phyla10 <- element_counts[element_counts$n>=10,]

result_data <- result_data[result_data$phylum %in% phyla10$phylum,]

rownames(result_data) <- result_data$phylum
result_data <- result_data[,-1]

rownames(result_data) <- gsub("^p__", "", rownames(result_data))

#过滤列#
#删除全为0的列
result_data <- result_data[, !colSums(result_data) == 0]
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=1,]
result_data <- result_data[,colnames(result_data) %in% ds_persentage$Column]

#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(result_data,"Fig2F.csv",quote = F)
ggplot2::ggsave("6_phyla_dss_frequency.pdf", plot = p, width = 4.5, height = 8, units = "in")


#2、主要门的DS family richness density统计####
result_data <- result(dss) %>% as.data.frame()
#过滤数据
library(dplyr)
#筛选>100的mags
phyla10 <- element_counts[element_counts$n>=10,]

result_data <- result_data[result_data$phylum %in% phyla10$phylum,]

rownames(result_data) <- result_data$phylum
result_data <- result_data[,-1]
rownames(result_data) <- gsub("^p__", "", rownames(result_data))

row_sums_result <- apply(result_data, 1, function(row) sum(row > 0))
row_sums_result_df <- as.data.frame(row_sums_result)
names(row_sums_result_df) <- c("Count")

magdsr <- data.frame(genome=dss$genome,dsr=rowSums(dss[,c(-1,-190:-197)] > 0))
dsr_size <- merge(magdsr,MAG_quality,by.x = 'genome',by.y = 'Name')


dsr_density <- data.frame(genome=dsr_size$genome,
                          density=dsr_size$dsr/dsr_size$Genome_Size*1000000,
                          richness=dsr_size$dsr,
                          length_Mb=dsr_size$Genome_Size/1000000
)
dsr_density <- merge(dsr_density,mag_lineage_split,by.x = 'genome',by.y = 'genome')
dsr_density$phylum <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", dsr_density$phylum)

result3 <- function(df) {
  output <- df %>%
    group_by(phylum) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE),  # 添加标准差的计算
          se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))  # 计算标准误
        )
      )
    )
  return(output)
}
dsr <- result3(dsr_density) %>% as.data.frame()

dsr <- dsr[dsr$phylum %in% phyla10$phylum,]

#绘制简单的柱状图
library(ggplot2)
dsr$phylum <- gsub("^p__", "", dsr$phylum)

dsr_phayla_density <- ggplot(dsr, aes(x = phylum, y = density_mean)) +
  geom_bar(stat = "identity", fill = "#52b69a") +
  geom_errorbar(
    aes(ymin = density_mean - 1.96*density_se, ymax = density_mean + 1.96*density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS diversity (per MAG per Mb)") +
  scale_y_continuous(limits = c(0,2),expand = c(0,0))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    #axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_density

dsr_phayla_richness <- ggplot(dsr, aes(x = phylum, y = richness_mean)) +
  geom_bar(stat = "identity", fill = "#34a0a4") +
  geom_errorbar(
    aes(ymin = richness_mean - 1.96*richness_se, ymax = richness_mean + 1.96*richness_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = richness_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$richness_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS diversity (per MAG)") +
  scale_y_continuous(limits = c(0,8.5),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_richness

dsr_phayla_length <- ggplot(dsr, aes(x = phylum, y = length_Mb_mean)) +
  geom_bar(stat = "identity", fill = "#ffa69e") +
  geom_errorbar(
    aes(ymin = length_Mb_mean - 1.96*length_Mb_se, ymax = length_Mb_mean + 1.96*length_Mb_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = length_Mb_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$length_Mb_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,7.5),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_length

library(cowplot)
g <- plot_grid(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
               labels = letters[1:3], nrow = 3,rel_heights = c(1, 1, 1.7),
               label_x = 0.09,            # 标签x坐标（1为最右侧）
               label_y = 0.9             # 标签y坐标（1为最上方）
)
g


# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g,'6_dss_19phyla.pdf',width = 4.5,height = 5.5,units = 'in')


#主要门中的DS number(per MAG) DS number density (per MAG per Mb)

DS_number <- data.frame(DS_num = rowSums(dss[,c(-1,-190:-197)]),
                        genome = dss$genome)
dsr_density <- merge(dsr_density, DS_number, by = "genome")
dsr_density$DS_density <- dsr_density$DS_num / dsr_density$length_Mb

dsr <- result3(dsr_density) %>% as.data.frame()

dsr <- dsr[dsr$phylum %in% phyla10$phylum,]
dsr$phylum <- gsub("^p__", "", dsr$phylum)
#huitu
dsr_phayla_DSnum <- ggplot(dsr, aes(x = phylum, y = DS_num_mean)) +
  geom_bar(stat = "identity", fill = "#be95c4") +
  geom_errorbar(
    aes(ymin = DS_num_mean - 1.96*DS_num_se, ymax = DS_num_mean + 1.96*DS_num_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = DS_num_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$DS_num_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,17),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_DSnum

dsr_phayla_DSdensity <- ggplot(dsr, aes(x = phylum, y = DS_density_mean)) +
  geom_bar(stat = "identity", fill = "#e0b1cb") +
  geom_errorbar(
    aes(ymin = DS_density_mean - 1.96*DS_density_se, ymax = DS_density_mean + 1.96*DS_density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = DS_density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$DS_density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,4),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_phayla_DSdensity

library(cowplot)
g <- plot_grid(dsr_phayla_length, 
               dsr_phayla_DSnum, dsr_phayla_DSdensity, 
               dsr_phayla_richness,dsr_phayla_density,
               labels = letters[1:5], nrow = 5,rel_heights = c(1, 1,1,1, 1.8),
               label_x = 0.09,            # 标签x坐标（1为最右侧）
               label_y = 0.9)             # 标签y坐标（1为最上方）

g

write.table(dsr, "dsr.tsv", sep = "\t", quote = F, row.names = F)
write.table(result_data, "result_data_heatmap.tsv", sep = "\t", quote = F)


# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g,'6_dss_19phyla_5.pdf',width = 4.5,height = 8,units = 'in')


#3. 蓝细菌门单独分析####

Cyanobacteriota_DS <- dss[dss$phylum == "p__Cyanobacteriota",]
#可以在order水平
Cyanobacteriota_DS <- merge(Cyanobacteriota_DS, MAG_quality, 
                            by.x = "genome", by.y = "Name")
Cyanobacteriota_DS <- merge(Cyanobacteriota_DS, MAG_metadata, 
                            by.x = "genome", by.y = "genome")
Cyanobacteriota_DS$DS_num <- rowSums(Cyanobacteriota_DS[,c(2:189)])

##不同order的DS分析####
###DS frequency####

Cyanobacteriota_DSfrequency <- Cyanobacteriota_DS[,c(1:197)]

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(order) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE) * 100))
  return(output)
}
result_data <- result(Cyanobacteriota_DSfrequency) %>% as.data.frame()

#过滤数据
library(dplyr)

element_counts <- Cyanobacteriota_DSfrequency %>%
  count(order)

#筛选>10的mags
order10 <- element_counts[element_counts$n>=10,]

result_data <- result_data[result_data$order %in% order10$order,]

rownames(result_data) <- result_data$order
result_data <- result_data[,-1]

rownames(result_data) <- gsub("^o__", "", rownames(result_data))

#过滤列#
#删除全为0的列
result_data <- result_data[, !colSums(result_data) == 0]
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=0.5,]
result_data <- result_data[,colnames(result_data) %in% ds_persentage$Column]

#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(result_data,"Fig2F.csv",quote = F)
ggplot2::ggsave("6_Cyanobacteriota_order_dss_frequency.pdf", plot = p, width = 3, height = 8, units = "in")

write.table(result_data, "result_data_heatmap_cyan.tsv", sep = "\t", quote = F)


###不同order的DS数量、density、MAG size#####
Cyanobacteriota_DS2 <- Cyanobacteriota_DS[,-2:-189]
#筛选>10的mags
order10 <- element_counts[element_counts$n>=10,]

Cyanobacteriota_DS2 <- Cyanobacteriota_DS2[Cyanobacteriota_DS2$order %in% order10$order,]


magdsr <- data.frame(genome=Cyanobacteriota_DS2$genome, 
                     DS_num= Cyanobacteriota_DS2$DS_num,
                     Genome_size = Cyanobacteriota_DS2$Genome_Size,
                     order = Cyanobacteriota_DS2$order
                     )
magdsr$DS_density <- magdsr$DS_num / magdsr$Genome_size *1000000
magdsr$length_Mb <- magdsr$Genome_size /1000000

Cyanobacteriota_DSfrequency <- Cyanobacteriota_DS[,c(1:197)]
Cyanobacteriota_DSfrequency <- Cyanobacteriota_DSfrequency[Cyanobacteriota_DSfrequency$genome %in% Cyanobacteriota_DS2$genome,]
ds_richness <- data.frame(genome=Cyanobacteriota_DSfrequency$genome,
                          richness=rowSums(Cyanobacteriota_DSfrequency[,-1] > 0))
magdsr <- merge(magdsr, ds_richness, by = "genome")
magdsr$order <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", magdsr$order)
magdsr$order <- gsub("^o__", "", magdsr$order)
magdsr$density <- magdsr$richness / magdsr$Genome_size *1000000

result3 <- function(df) {
  output <- df %>%
    group_by(order) %>%
    summarise(
      across(
        where(is.numeric),
        list(
          mean = ~mean(., na.rm = TRUE),
          sd = ~sd(., na.rm = TRUE),  # 添加标准差的计算
          se = ~sd(., na.rm = TRUE) / sqrt(sum(!is.na(.)))  # 计算标准误
        )
      )
    )
  return(output)
}
dsr <- result3(magdsr) %>% as.data.frame()


#绘制简单的柱状图
library(ggplot2)

dsr_order_density <- ggplot(dsr, aes(x = order, y = density_mean)) +
  geom_bar(stat = "identity", fill = "#52b69a") +
  geom_errorbar(
    aes(ymin = density_mean - 1.96*density_se, ymax = density_mean + 1.96*density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS diversity (per MAG per Mb)") +
  scale_y_continuous(limits = c(0,5),expand = c(0,0))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
    #axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_order_density

dsr_order_richness <- ggplot(dsr, aes(x = order, y = richness_mean)) +
  geom_bar(stat = "identity", fill = "#34a0a4") +
  geom_errorbar(
    aes(ymin = richness_mean - 1.96*richness_se, ymax = richness_mean + 1.96*richness_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = richness_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$richness_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DS diversity (per MAG)") +
  scale_y_continuous(limits = c(0,20),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_order_richness

dsr_order_length <- ggplot(dsr, aes(x = order, y = length_Mb_mean)) +
  geom_bar(stat = "identity", fill = "#ffa69e") +
  geom_errorbar(
    aes(ymin = length_Mb_mean - 1.96*length_Mb_se, ymax = length_Mb_mean + 1.96*length_Mb_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = length_Mb_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$length_Mb_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,6),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_order_length

library(cowplot)
g <- plot_grid(dsr_order_length, dsr_order_richness,dsr_order_density,
               labels = letters[1:3], nrow = 3,rel_heights = c(1, 1, 1.7),
               label_x = 0.09,            # 标签x坐标（1为最右侧）
               label_y = 0.9             # 标签y坐标（1为最上方）
)
g

# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g,'6_dss_19phyla.pdf',width = 4.5,height = 5.5,units = 'in')


#主要门中的DS number(per MAG) DS number density (per MAG per Mb)

#huitu
dsr_order_DSnum <- ggplot(dsr, aes(x = order, y = DS_num_mean)) +
  geom_bar(stat = "identity", fill = "#be95c4") +
  geom_errorbar(
    aes(ymin = DS_num_mean - 1.96*DS_num_se, ymax = DS_num_mean + 1.96*DS_num_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = DS_num_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$DS_num_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,35),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_order_DSnum

dsr_order_DSdensity <- ggplot(dsr, aes(x = order, y = DS_density_mean)) +
  geom_bar(stat = "identity", fill = "#e0b1cb") +
  geom_errorbar(
    aes(ymin = DS_density_mean - 1.96*DS_density_se, ymax = DS_density_mean + 1.96*DS_density_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = DS_density_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(dsr$DS_density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "MAG size (Mb)") +
  scale_y_continuous(limits = c(0,7),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text(face = "bold", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
dsr_order_DSdensity

library(cowplot)
g <- plot_grid(dsr_order_length, 
               dsr_order_DSnum, dsr_order_DSdensity, 
               dsr_order_richness,dsr_order_density,
               labels = letters[1:5], nrow = 5,rel_heights = c(1, 1,1,1, 1.8),
               label_x = 0.09,            # 标签x坐标（1为最右侧）
               label_y = 0.9)             # 标签y坐标（1为最上方）

g

# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g,'6_Cyanobacteriota_order.pdf',width = 2,height = 8,units = 'in')

write.table(dsr, "dsr_cyan.tsv", 
  sep = "\t", quote = F, row.names = F)


#基因组大小分界####


#4. 黏细菌分析####

Myxococcota_DS <- dss[dss$phylum == "p__Myxococcota",]
#可以在class水平
Myxococcota_DS <- merge(Myxococcota_DS, MAG_quality, 
                            by.x = "genome", by.y = "Name")
Myxococcota_DS <- merge(Myxococcota_DS, MAG_metadata, 
                            by.x = "genome", by.y = "genome")
Myxococcota_DS$DS_num <- rowSums(Myxococcota_DS[,c(2:189)])

##不同class的DS分析####
###DS frequency####

Myxococcota_DSfrequency <- Myxococcota_DS[,c(1:197)]

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(class) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE) * 100))
  return(output)
}
result_data <- result(Myxococcota_DSfrequency) %>% as.data.frame()

#过滤数据
library(dplyr)

element_counts <- Myxococcota_DSfrequency %>%
  count(class)

#筛选>5的mags
class5 <- element_counts[element_counts$n>=5,]

result_data <- result_data[result_data$class %in% class5$class,]

rownames(result_data) <- result_data$class
result_data <- result_data[,-1]

rownames(result_data) <- gsub("^c__", "", rownames(result_data))

#过滤列#
#删除全为0的列
result_data <- result_data[, !colSums(result_data) == 0]
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=0.5,]
result_data <- result_data[,colnames(result_data) %in% ds_persentage$Column]

#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(result_data,"Fig2F.csv",quote = F)
ggplot2::ggsave("6_Myxococcota_class_DS_frequency.pdf", plot = p, width = 2.5, height = 8, units = "in")

write.table(result_data, "result_data_myxoco.tsv",
  sep = "\t", quote = F)

#5. p__Gemmatimonadota分析####

Gemmatimonadota_DS <- dss[dss$phylum == "p__Gemmatimonadota",]
#在phylum水平
Gemmatimonadota_DS <- merge(Gemmatimonadota_DS, MAG_quality, 
                        by.x = "genome", by.y = "Name")
Gemmatimonadota_DS <- merge(Gemmatimonadota_DS, MAG_metadata, 
                        by.x = "genome", by.y = "genome")
Gemmatimonadota_DS$DS_num <- rowSums(Gemmatimonadota_DS[,c(2:189)])

##不同class的DS分析####
###DS frequency####

Gemmatimonadota_DSfrequency <- Gemmatimonadota_DS[,c(1:197)]

#分组统计
result <- function(df) {
  output <- df %>%
    group_by(phylum) %>%
    summarise(across(where(is.numeric), ~mean(. != 0, na.rm = TRUE) * 100))
  return(output)
}
result_data <- result(Gemmatimonadota_DSfrequency) %>% as.data.frame()

#过滤数据
library(dplyr)

element_counts <- Gemmatimonadota_DSfrequency %>%
  count(phylum)

#筛选>5的mags
#phylum5 <- element_counts[element_counts$n>=5,]

#result_data <- result_data[result_data$class %in% class5$class,]

rownames(result_data) <- result_data$phylum
result_data <- result_data[,-1]

rownames(result_data) <- gsub("^p__", "", rownames(result_data))

#过滤列#
#删除全为0的列
#result_data <- result_data[, !colSums(result_data) == 0]
#筛选主要的DSs
# 计算每一列非0数字所占的百分比并输出到新数据框
percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=1,]
result_data <- result_data[,colnames(result_data) %in% ds_persentage$Column]

#绘制热图
library(pheatmap)

# 设置颜色
#colors <- colorRampPalette(c('white', colorRampPalette('#F9B14F', '#8F2731')))(1000)
colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p <- pheatmap(t(result_data),
              #color = colorRampPalette(colors)(1000),
              #color = colorRampPalette(c('white','lightpink','red'))(1000),
              color = colors,
              cluster_cols = F,
              cluster_rows = F,
              angle_col = 90
)
p
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
#write.csv(result_data,"Fig2F.csv",quote = F)
ggplot2::ggsave("6_Gemmatimonadota_phylum_DS_frequency.pdf", plot = p, width = 1.9, height = 8, units = "in")

write.table(result_data, "result_data_gemma.tsv",
  sep = "\t", quote = F)

save.image("./07_phylum_specialization.RData")
load("./07_phylum_specialization.RData")


