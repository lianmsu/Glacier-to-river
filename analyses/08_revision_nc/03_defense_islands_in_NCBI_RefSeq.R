# =============================================================================
# Glacier-to-river -- revision nc / Defense islands in NCBI RefSeq genomes
# =============================================================================
# Repeat the defense-island analysis on 22,942 NCBI RefSeq genomes to test
# whether the observed DI patterns generalize.
#
# Manuscript      : Revision
# Main outputs    : 1_summary_stats_DI_n_normalized.pdf,
#                   03_defense_islands_in_NCBI_RefSeq.RData
# Original script : 93_NC_comments/03_DIs_NCBI_RefSeq/1_DIs_NCBI_refseq_genome.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "93_NC_comments/03_DIs_NCBI_RefSeq"))
library(dplyr)
library(stringr)
library(ggplot2)

#DIs no antidefense gene
DI_no <- read.delim2("22942genome_defense_islands_summary.csv",check.names = F)

DI_no <- DI_no[,c(1:6)]
DI_no <- DI_no %>%
  mutate(
    contig = sub("(-[^_]+)_.*", "\\1", Contig),
    Genome = sub("-[^_]+_.*", "", Contig)
  )

bac_22143_info <- 
  read.delim2("NCBI_refseq_22143_bacterial_genomes_20260305_withtax.tsv")
arc_799_info <-
  read.delim2("NCBI_refseq_799_archaeal_genomes_20260305_withtax.tsv")
refseq_22942_info <- rbind(bac_22143_info,arc_799_info)

library(tidyr)
# 拆分 classification 列
mag_lineage_split <- refseq_22942_info %>%
  separate(
    Standard_tax,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )
# 查看拆分后的数据
head(mag_lineage_split)
colnames(mag_lineage_split)[2] <- 'genome'

DI_no <- merge(DI_no, mag_lineage_split, by.x = "Genome",by.y = "genome")


DI_no$phylum <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", DI_no$phylum)
DI_no$phylum <- gsub("^p__", "", DI_no$phylum)
phylum_count <- table(DI_no$phylum) %>% as.data.frame()

DI_no$DI_length <- DI_no$`Last Record Num` - DI_no$`First Record Num` +1
DI_no$DG_proportion <- DI_no$`Total Gene Count` / DI_no$DI_length *100


#不同物种中DI数量、DI中DG数量、DG占DI总gene的比例; 同时用物种系统发育树进行展示####
summary_stats <- DI_no %>%
  group_by(phylum) %>%
  summarise(
    n = n(),
    across(c(`Total Gene Count`, DG_proportion),
           list(
             mean = ~mean(., na.rm = TRUE),
             sd = ~sd(., na.rm = TRUE),
             se = ~sd(., na.rm = TRUE) / sqrt(n())
           )
    ))
summary_stats <- summary_stats[summary_stats$n >=10,]

## DI数量归一化到物种####
head(summary_stats)
mag_lineage_split2 <- mag_lineage_split
table(mag_lineage_split2$phylum)
mag_lineage_split2$phylum <-
  sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", mag_lineage_split2$phylum)
mag_phyla_num <- 
  table(mag_lineage_split2$phylum) %>% as.data.frame()
colnames(mag_phyla_num) <- c("phyla","MAG_num")

mag_phyla_num$phyla <- gsub("^p__", "", mag_phyla_num$phyla)
head(mag_phyla_num)


summary_stats_2 <- summary_stats %>%
  left_join(mag_phyla_num, by = c("phylum" = "phyla")) %>%
  mutate(n_normalized = n / MAG_num)   # 归一化结果存为新列

summary_stats_DI_n_normalized <- 
  ggplot(summary_stats_2, aes(x = phylum, y = n_normalized)) +
  geom_bar(stat = "identity", fill = "#ff9f1c") +
  # geom_errorbar(
  #   aes(ymin = density_mean - 1.96*density_se, ymax = density_mean + 1.96*density_se),
  #   width = 0.2,  # 误差线宽度
  #   position = position_dodge(width = 0.7),  # 与条形图对齐
  #   color = "grey"  # 误差线颜色
  # ) +
  # geom_point(
  #   aes(y = density_mean),
  #   position = position_dodge(width = 0.7),  # 与条形图对齐
  #   color = "grey",  # 点颜色
  #   size = 1.5  # 点大小
  # )+
  #geom_hline(yintercept = mean(dsr$density_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DI number (normalized to genome number)") +
  scale_y_continuous(limits = c(0,0.3),expand = c(0,0))+
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
summary_stats_DI_n_normalized

ggsave("1_summary_stats_DI_n_normalized.pdf", 
       summary_stats_DI_n_normalized, 
       width = 5, 
       height = 4.5, 
       units = "in")


save.image("./03_defense_islands_in_NCBI_RefSeq.RData")
load("./03_defense_islands_in_NCBI_RefSeq.RData")





