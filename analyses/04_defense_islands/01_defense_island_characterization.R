# =============================================================================
# Glacier-to-river -- defense islands / Defense-island characterization
# =============================================================================
# Number of DIs and DG density across the major lineages, host contribution
# by region/ecosystem/medium, DS composition of DIs and the DI O/E heatmap.
#
# Manuscript      : Fig. 5a-d, Fig. S11
# Main outputs    : 1_summary_stats_DI_n_normalized.pdf, 1_DS_volcano_fisher.pdf,
#                   01_defense_island_characterization.RData
# Original script : 20gaoyuan/14_DIs/1_DIs.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))

#DIs no antidefense gene
DI_no <- read.delim2("5455MAG_defense_islands_summary_no_anti.csv",check.names = F)
library(dplyr)
library(stringr)
#提取MAGs名称
DI_no <- DI_no %>%
  mutate(genome = str_extract(Contig, ".*bin\\.\\d+"))
#提取contig名称
DI_no <- DI_no %>%
  mutate(contig = str_extract(Contig, ".*bin\\.\\d+\\-\\d+"))


# MAG_metadata <- read.csv("../01MAG_info/5455MAG_metadata.csv")
MAG_metadata <- read.delim2("../01MAG_info/MAG5455_metadata_new.tsv")

MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'Yes',"Plateau",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'No',"Plain",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'glacier',"Glacier",MAG_metadata$Plateau)
colnames(MAG_metadata)[13] <- 'Region'

MAG_metadata$medium <- ifelse(MAG_metadata$Type == 'Cryoconite',"Glacier",MAG_metadata$Type)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Ice',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Snow',"Glacier",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'particle',"Plankton",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'water',"Plankton",MAG_metadata$medium)
MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'sediment',"Benthon",MAG_metadata$medium)

MAG_quality <- read.delim2("../01MAG_info/quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]

DF_mechanism <- read.delim2('../10Defensefinder/df_List_Systems_lpw.csv',sep = ',')
DF_mechanism$Effector[DF_mechanism$Effector == ""] <- "Unknown"
DF_mechanism$Mechnism[DF_mechanism$Mechnism == ""] <- "Unknown"

mag_baclineage <- read.delim2("../10Defensefinder/gtdbtk.bac120.summary.tsv")
mag_arclineage <- read.delim2("../10Defensefinder/gtdbtk.ar53.summary.tsv")
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


DI_no <- merge(DI_no, mag_lineage_split, by.x = "genome",by.y = "genome")
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
summary_stats <- summary_stats[summary_stats$n >1,]

library(ggplot2)
summary_stats_DI_n <- ggplot(summary_stats, aes(x = phylum, y = n)) +
  geom_bar(stat = "identity", fill = "#00b4d8") +
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
       y = "DI number") +
  scale_y_continuous(limits = c(0,50),expand = c(0,0))+
  #coord_cartesian(ylim = c(0, 1.1 * max(row_sums_result_df$Count)))+
  #theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),  # 横轴标签旋转
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
summary_stats_DI_n

summary_stats_DI_DGnum <- ggplot(summary_stats, aes(x = phylum, y = `Total Gene Count_mean`)) +
  geom_bar(stat = "identity", fill = "#98c1d9") +
  geom_errorbar(
    aes(ymin = `Total Gene Count_mean` - 1.96*`Total Gene Count_se`, 
        ymax = `Total Gene Count_mean` + 1.96*`Total Gene Count_se`),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = `Total Gene Count_mean`),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(summary_stats$`Total Gene Count_mean`), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DG number (per DI)") +
  scale_y_continuous(limits = c(0,9),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
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
summary_stats_DI_DGnum

summary_stats_DI_DGproportion <- ggplot(summary_stats, aes(x = phylum, y = DG_proportion_mean)) +
  geom_bar(stat = "identity", fill = "#90e0ef") +
  geom_errorbar(
    aes(ymin = DG_proportion_mean - 1.96*DG_proportion_se, 
        ymax = DG_proportion_mean + 1.96*DG_proportion_se),
    width = 0.2,  # 误差线宽度
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey"  # 误差线颜色
  ) +
  geom_point(
    aes(y = DG_proportion_mean),
    position = position_dodge(width = 0.7),  # 与条形图对齐
    color = "grey",  # 点颜色
    size = 1.5  # 点大小
  )+
  geom_hline(yintercept = mean(summary_stats$DG_proportion_mean), linetype = 2)+
  theme_minimal() +
  labs(title = NULL,
       x = NULL,
       y = "DG proportion (per DI, %)") +
  scale_y_continuous(limits = c(0,115),expand = c(0,0), label = function(x) sprintf("%.1f", x) )+
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
summary_stats_DI_DGproportion

library(cowplot)
g <- plot_grid(summary_stats_DI_n, summary_stats_DI_DGnum,summary_stats_DI_DGproportion,
               labels = letters[1:3], nrow = 3,rel_heights = c(1, 1, 1.6),
               label_x = 0.15,            # 标签x坐标（1为最右侧）
               label_y = 0.95             # 标签y坐标（1为最上方）
)
g

write.table(summary_stats, "summary_stats_dis.tsv",
  sep = "\t", quote = F, row.names = F)

# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
#write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g,'1_DIs_9phyla.pdf',width = 3.5,height = 7,units = 'in')

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
  geom_bar(stat = "identity", fill = "#00b4d8") +
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
       y = "DI number (normalized to MAG number)") +
  scale_y_continuous(limits = c(0,0.2),expand = c(0,0))+
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

#这些DI宿主的环境来源####
DI_host <- unique(DI_no$genome) %>% as.data.frame() %>% setNames("genome")
DI_host <- merge(DI_host, mag_lineage_split, by = "genome" ,all = F)
DI_host <- merge(DI_host, MAG_metadata, by = "genome", all = F)
table(DI_host$phylum)

DI_host_env <- select(DI_host, ecosystem, Region, Type)

DI_host_env_ecosystem <- DI_host_env %>%
  group_by(ecosystem) %>%
  summarise(
    n = n()) %>%
  mutate(percentage = n / sum(n) * 100)
DI_host_env_ecosystem$X <- "X"
DI_host_env_ecosystem$ecosystem <- factor(DI_host_env_ecosystem$ecosystem, 
                                          levels = c("Glacier","YJ","LCJ","HH","CJ"))

p_DI_host_env_ecosystem <- 
ggplot(DI_host_env_ecosystem, aes(x = X, y = percentage, fill = ecosystem)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  labs(
    x = "Ecosystem",
    y = "DI host source (%)",
    fill = ""
  ) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_DI_host_env_ecosystem

DI_host_env_Region <- DI_host_env %>%
  group_by(Region) %>%
  summarise(
    n = n()) %>%
  mutate(percentage = n / sum(n) * 100)
DI_host_env_Region$X <- "X"
DI_host_env_Region$Region <- factor(DI_host_env_Region$Region, 
                                    levels = c("Glacier","Plateau","Plain"))

p_DI_host_env_Region <- 
  ggplot(DI_host_env_Region, aes(x = X, y = percentage, fill = Region)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  labs(
    x = "Region",
    y = "DI host source (%)",
    fill = ""
  ) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_DI_host_env_Region

DI_host_env_Type <- DI_host_env %>%
  group_by(Type) %>%
  summarise(
    n = n()) %>%
  mutate(percentage = n / sum(n) * 100)
DI_host_env_Type$X <- "X"
DI_host_env_Type$Type <- factor(DI_host_env_Type$Type,
                                    levels = c("Snow","Ice","Cryoconite",
                                               "particle","sediment","water"))

p_DI_host_env_Type <- 
  ggplot(DI_host_env_Type, aes(x = X, y = percentage, fill = Type)) +
  geom_col(position = "stack", width = 0.6) +
  scale_fill_manual(values = c("Snow"="#E74A35","Ice"="#4DBAD4","Cryoconite"="#00A48A",
                               "water"="#3B5387","particle"= "#F39B7E","sediment"="#8491B4")) +  # 手动设置颜色
  labs(
    x = "Medium",
    y = "DI host source (%)",
    fill = ""
  ) +
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
p_DI_host_env_Type

library(cowplot)
g2 <- plot_grid(p_DI_host_env_ecosystem, p_DI_host_env_Region,p_DI_host_env_Type,
               labels = letters[1:3], nrow = 1 #rel_heights = c(1, 1, 1.6),
              # label_x = 0.15,            # 标签x坐标（1为最右侧）
              # label_y = 0.95             # 标签y坐标（1为最上方）
)
g2

# ggarrange(dsr_phayla_length, dsr_phayla_richness,dsr_phayla_density,
#           ncol=1, nrow=3, common.legend = TRUE, legend="bottom"
#           ) 

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
#write.csv(dsr,"Fig2CDE.csv",quote = F)
ggsave(plot = g2,'1_DI_host_source_proportion.pdf',width = 4.3,height = 3,units = 'in')

write.table(DI_host_env_ecosystem, "DI_host_env_ecosystem.tsv",
  sep = "\t", quote = F, row.names = F)
write.table(DI_host_env_Region, "DI_host_env_Region.tsv",
  sep = "\t", quote = F, row.names = F)
write.table(DI_host_env_Type, "DI_host_env_Type.tsv",
  sep = "\t", quote = F, row.names = F)

#DI中各类DS所占比例####
#计算所有DIs中，各类DSs的数量所占比例
alldss <- colSums(DI_no[, 8:75]) %>% as.data.frame() %>% setNames("dssnumber")
sum(DI_no$`Total Gene Count`) #774
sum(alldss) #366

# 筛选出dssnumber 小于5的行，并将它们合并为一行
alldss_filtered <- alldss %>%
  filter(dssnumber < 9) %>%
  summarise(combined_dssgene = sum(dssnumber)) %>% setNames("dssnumber")

# 将合并后的行添加到原始数据框的末尾
alldss <- rbind( alldss,alldss_filtered)
rownames(alldss)[nrow(alldss)] <- "Others"
alldss$ds <- rownames(alldss)

#绘制圈图
library(ggplot2)
library(ggpmisc)

# 绘制中空饼图
mycolor=c("#FB9A99","#7294D4","#FDBF6F","#A6CEE3","#1F78B4","#CAB2D6",
          "#FF7F00","#0B775E","#33A02C","#6A3D9A","#B2DF8A","#8a817c")
pie <- ggplot(alldss[alldss$dssnumber>=9,], aes(x = 2, y = dssnumber, fill = ds)) +
  geom_bar(stat = "identity", width = 1, color = "white") +  # 使用 geom_bar 模拟饼图
  coord_polar(theta = "y") +  # 将坐标转换为极坐标
  xlim(0.5,2.5)+
  scale_fill_manual(values = mycolor) +
  geom_text( aes(label = dssnumber), 
             position = position_stack(vjust = 0.5),
             color = "white", # 文本颜色
             size = 5) + # 文本大小
  guides(fill = guide_legend(title = "DSs"))+
  theme_void() +  # 移除背景元素
  theme(legend.position = "none",  # 移除默认图例
        axis.text.x = element_blank(),
        axis.ticks = element_blank())  # 移除坐标轴标签和刻度
pie
# 创建自定义图例
legend_plot <- ggplot(alldss[alldss$dssnumber>=9,], aes(x = 2, y = dssnumber, fill = ds)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = mycolor) +
  theme_void() +
  guides(fill = guide_legend(title = ""))

# 提取图例
legend <- get_legend(legend_plot)

# 创建组合图
combined_plot <- ggdraw() +
  draw_plot(pie, 0, 0, 1, 1) +
  draw_grob(legend, 0.4, 0.4, 0.2, 0.25)  # 调整图例的位置和大小
combined_plot

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
#write.csv(alldss[alldss$dssnumber>50,],"fig5e.csv")
ggsave(plot = combined_plot,"1_DI_dss_pie.pdf",width = 3.4, height = 3.1, units = "in")

write.table(alldss, "alldss_dis.tsv",
  sep = "\t", quote = F, row.names = F)

#DI中DG之间的共现关系/携带DS的倾向####

#携带DS的倾向
setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))

DF_res <- read.delim2("../10Defensefinder/merged_5467MAG_defense_finder_systems.tsv",check.names = F,
                      sep = '\t')

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
#count_matrix <- count_matrix[-1]  # 移除genome列

#genome_defense_num <- rowSums(count_matrix) %>% as.data.frame() %>% setNames("DS_num")
#genome_defense_num$genome <- rownames(genome_defense_num)

count_matrix_noDI <- count_matrix[!count_matrix$genome %in% DI_no$genome,]
no_di_ds <- count_matrix_noDI[,-1] %>% colSums() %>% as.data.frame() %>% setNames('non DIs')
di_ds <- DI_no[,8:75] %>% colSums() %>% as.data.frame() %>% setNames('DIs')

di_ds <- merge(no_di_ds,di_ds,by = 0)
di_ds[is.na(di_ds)] <- 0
rownames(di_ds) <- di_ds$Row.names
di_ds <- di_ds[,-1]
di_ds_chi <- chisq.test(di_ds)
di_ds_oe <- di_ds_chi$observed/di_ds_chi$expected %>% as.data.frame() %>%
  setNames(c("non DIs","DIs"))
di_ds_oe <- di_ds_oe[!is.nan(di_ds_oe$DIs),]

#绘制热图
library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数

cairo_pdf(width = 8,height = 1.8,filename = "1_di_ds_oe.pdf")
col_fun <- colorRamp2(c(0,1, 2), c("white","#70d6ff", "#ff6b6b"))  # 粉色到绿色的渐变
p_di_ds_oe <- Heatmap(t(di_ds_oe),
                      col = col_fun,
                      cluster_rows = F,
                      cluster_columns = T,
                      name = 'O/E',
                      row_title = NULL,
                      row_split = colnames(di_ds_oe),
                      row_names_side = "left", # 行标题显示在左侧
                      # row_names_gp = gpar(fontsize = 10),
                      column_names_gp = gpar(fontsize = 8),
                      column_names_rot = 60,
                      rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
                      heatmap_legend_param = list(
                        title = "O/E", at = c(0,1, 2), 
                        labels = c("0","1", "≥2")),
)
p_di_ds_oe
setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
dev.off()

write.table(di_ds_oe, "di_ds_oe.tsv", sep = "\t", quote = F)

#岛内系统 vs 无岛基因组系统：逐类型 Fisher 精确检验 + 火山图####
#上方整体卡方检验在大量零单元格下 p 值失效（NaN），故改为逐类型双侧 Fisher 检验并做多重比较校正

# 1) 将防御系统定位到防御岛（contig + 记录号区间重叠）
DF_defense <- DF_defense |>
  mutate(
    contig  = str_extract(sys_id, ".*bin\\.\\d+\\-\\d+"),
    genome  = str_extract(sys_id, ".*bin\\.\\d+"),
    beg_rec = as.integer(str_extract(sys_beg, "\\d+$")),
    end_rec = as.integer(str_extract(sys_end, "\\d+$"))
  )

islands <- DI_no |>
  mutate(
    contig = str_extract(Contig, ".*bin\\.\\d+\\-\\d+"),
    first  = as.integer(`First Record Num`),
    last   = as.integer(`Last Record Num`)
  ) |>
  select(contig, first, last)

in_island <- rep(FALSE, nrow(DF_defense))
for (ct in unique(DF_defense$contig)) {
  idx <- which(DF_defense$contig == ct)
  isl <- islands[islands$contig == ct, , drop = FALSE]
  if (nrow(isl) == 0) next
  for (i in idx) {
    in_island[i] <- any(DF_defense$beg_rec[i] <= isl$last &
                        DF_defense$end_rec[i] >= isl$first)
  }
}
DF_defense$in_island <- in_island

# 2) 岛内系统 vs 无岛基因组系统（逐类型系统数）
island_type <- DF_defense[DF_defense$in_island, ] |>
  count(type, name = "in_island")
free_type <- DF_defense[!DF_defense$genome %in% unique(DI_no$genome), ] |>
  count(type, name = "free_genome")

tab <- merge(
  data.frame(type = free_type$type, free_genome = free_type$free_genome),
  data.frame(type = island_type$type, in_island = island_type$in_island),
  by = "type", all = TRUE
)
tab[is.na(tab)] <- 0
tot_island <- sum(tab$in_island)
tot_free   <- sum(tab$free_genome)

# 3) 逐类型 2x2 Fisher 精确检验（双侧）+ BH/Bonferroni 校正
fisher_per_type <- lapply(seq_len(nrow(tab)), function(i) {
  a <- tab$in_island[i]
  b <- tab$free_genome[i]
  if (a + b == 0) return(NULL)  # 该类型既不在岛内也不在无岛基因组，不参与检验
  m2 <- matrix(c(a, tot_island - a, b, tot_free - b), nrow = 2)
  ft <- fisher.test(m2, alternative = "two.sided")
  data.frame(type = tab$type[i], island = a, free_genome = b,
             OR = as.numeric(ft$estimate), p_value = ft$p.value)
}) |> bind_rows()

fisher_per_type <- fisher_per_type |>
  mutate(
    p_FDR = p.adjust(p_value, "BH"),
    p_bonf = p.adjust(p_value, "bonferroni")
  ) |>
  arrange(p_value)

write.table(fisher_per_type, "di_ds_fisher_systems_island_vs_free.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 4) 火山图（Haldane 校正避免 0/Inf）
library(ggrepel)
set.seed(1)
fisher_plot <- fisher_per_type |>
  mutate(
    OR_h = ((island + 0.5) * (tot_free - free_genome + 0.5)) /
           ((free_genome + 0.5) * (tot_island - island + 0.5)),
    log2OR = log2(OR_h),
    neglogFDR = -log10(p_FDR),
    group = case_when(
      p_FDR >= 0.05 ~ "n.s.",
      OR_h > 1 ~ "enriched in DIs",
      TRUE ~ "depleted in DIs"
    )
  )

p_volcano <- ggplot(fisher_plot, aes(x = log2OR, y = neglogFDR)) +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60") +
  geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
  geom_point(aes(color = group), size = 1.8) +
  scale_color_manual(
    values = c("enriched in DIs" = "#c0392b",
               "depleted in DIs" = "#2471a3",
               "n.s." = "grey80"),
    name = NULL
  ) +
  geom_text_repel(
    data = subset(fisher_plot, group != "n.s."),
    aes(label = type), size = 3, max.overlaps = 20
  ) +
  labs(x = expression(log[2]~"(odds ratio, island vs island-free)"),
       y = expression(-log[10]~"FDR")) +
  theme_bw() +
  theme(legend.position = "bottom")
p_volcano

ggsave("1_DS_volcano_fisher.pdf", p_volcano, width = 6.5, height = 5.5, units = "in")

#DI中，DS之间的共现关系
#读取所有DSs信息
dss_all_sum <- count_matrix[,-1] %>% colSums() %>% as.data.frame() %>% setNames("DS")
dss_all_sum$DF <- rownames(dss_all_sum)
di_ds <- DI_no[,8:75] %>% colSums() %>% as.data.frame() %>% setNames('DI_DS')
di_ds$DF <- rownames(di_ds)

allmag_di_ds <- merge(dss_all_sum,di_ds,by ="DF" ,all = T)
allmag_di_ds[is.na(allmag_di_ds)] <- 0

allmag_di_ds$`non DIs` <- allmag_di_ds$DS - allmag_di_ds$DI_DS
rownames(allmag_di_ds) <- allmag_di_ds$DF
f_mag_di_ds <- select(allmag_di_ds,DI_DS,`non DIs`)
f_mag_di_ds <- f_mag_di_ds[f_mag_di_ds$DI_DS>0,]
f_mag_di_ds <- f_mag_di_ds[f_mag_di_ds$`non DIs`>0,]

#一次性对所有数据进行fisher，超内存
library(rstatix)
df_fisher <- fisher.test(f_mag_di_ds,alternative = "two.sided",
                         workspace = 1e+9,simulate.p.value = TRUE) #超内存
df_fisher <- fisher.test(f_mag_di_ds[c(4,5),],alternative = "two.sided",
                         workspace = 1e+9)

calculate_fisher_tests <- function(f_mag_di_ds) {
  # 初始化一个空的数据框来存储结果
  results_df <- data.frame(
    Row1 = character(),
    Row2 = character(),
    Odds_Ratio = numeric(),
    P_Value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # 获取f_mag_di_ds的行数
  n_rows <- nrow(f_mag_di_ds)
  
  # 对每两行进行循环
  for (i in 1:(n_rows)) {
    for (j in (i):n_rows) {
      # 提取当前行和下一行的行名
      row_names <- c(rownames(f_mag_di_ds)[i], rownames(f_mag_di_ds)[j])
      
      # 提取当前行和下一行的数据
      data_to_test <- f_mag_di_ds[c(i, j), ]
      
      # 确保data_to_test是一个2x2矩阵
      if (ncol(data_to_test) == 2) {
        # 执行Fisher精确检验
        fisher_result <- fisher.test(data_to_test, alternative = "two.sided", workspace = 1e+9)
        
        # 将行名、odds ratio和p.value添加到结果数据框中
        results_df <- rbind(results_df, data.frame(
          Row1 = row_names[1],
          Row2 = row_names[2],
          Odds_Ratio = fisher_result$estimate[1],
          P_Value = fisher_result$p.value,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  return(results_df)
  # 示例调用
  # 假设f_mag_di_ds是一个两列的数据框
  # f_mag_di_ds <- data.frame(
  #   A = c(1, 2, 3, 4),
  #   B = c(5, 6, 7, 8),
  #   row.names = c("Row1", "Row2", "Row3", "Row4")
  # )
  # results <- calculate_fisher_tests(f_mag_di_ds)
  # print(results)
}

results_df <- calculate_fisher_tests(f_mag_di_ds)

library(tidyr)
tb_od = pivot_wider(results_df[,1:3], names_from = c('Row2'), values_from = 'Odds_Ratio') %>% as.data.frame()
#kable(tb_wide_new)
rownames(tb_od) <- tb_od$Row1
tb_od <- tb_od[,-1]
tb_od[is.na(tb_od)] <- 0

try <- tb_od+t(tb_od)
diag(try) <- 1

library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数
col_fun <- colorRamp2(c(0,1,6), c("#70d6ff", "white", "#9f86c0"))
DS_or <- Heatmap(try,
                 col = col_fun,
                 name = 'P-value',
                 cluster_rows = F,
                 cluster_columns = F,
                 row_names_gp = gpar(fontsize = 6),
                 column_names_gp = gpar(fontsize = 6),
                 heatmap_legend_param = list(
                   title = "P-value", at = c(0,1,2,4,6), 
                   labels = c("0","1","2","4", "≥6")),
                 show_heatmap_legend = F
                 )
DS_or
library(ggplot2)
setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
cairo_pdf(width = 6,height = 6,filename = "1_DS_OR_nolegend.pdf")
DS_or
dev.off()

write.table(try, "try_dis_heatmap.tsv",
  sep = "\t", quote = F)

#######绘制p-value
library(tidyr)
tb_p = pivot_wider(results_df[,c(1,2,4)], names_from = c('Row2'), values_from = 'P_Value') %>% as.data.frame()
#kable(tb_wide_new)
rownames(tb_p) <- tb_p$Row1
tb_p <- tb_p[,-1]
tb_p[is.na(tb_p)] <- 0

p_value <- tb_p+t(tb_p)
diag(p_value) <- 1

library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数
col_fun2 <- colorRamp2(c(0,0.05,0.1), c("#f5cb5c", "white", "#f7f4ea"))
DS_p_value <- Heatmap(p_value,
                      col = col_fun2,
                      name = 'P-value',
                      cluster_rows = F,
                      cluster_columns = F,
                      row_names_gp = gpar(fontsize = 6),
                      column_names_gp = gpar(fontsize = 6),
                      heatmap_legend_param = list(
                        title = "P-value", at = c(0,0.05,0.1),
                        labels = c("0","0.05", "≥0.1")),
                      show_heatmap_legend = F
                      )
DS_p_value

write.table(p_value, "p_value_dis_heatmap.tsv",
  sep = "\t", quote = F)


setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
cairo_pdf(width = 6,height = 6,filename = "1_DS_P-value_nolegend.pdf")
DS_p_value
dev.off()


#多少DI与MGE相关？####

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))







#DI上非DG的基因的功能注释####





#携带Antidefense的代表DI的genetic context####

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))

DI_no <- read.delim2("5455MAG_defense_islands_summary_no_anti.csv",check.names = F)
library(dplyr)
library(stringr)
#提取MAGs名称
DI_no <- DI_no %>%
  mutate(genome = str_extract(Contig, ".*bin\\.\\d+"))
#提取contig名称
DI_no <- DI_no %>%
  mutate(contig = str_extract(Contig, ".*bin\\.\\d+\\-\\d+"))

DI_withanti <- read.delim2("5455MAG_defense_islands_summary_with_anti.csv",check.names = F)
#提取MAGs名称
DI_withanti <- DI_withanti %>%
  mutate(genome = str_extract(Contig, ".*bin\\.\\d+"))
#提取contig名称
DI_withanti <- DI_withanti %>%
  mutate(contig = str_extract(Contig, ".*bin\\.\\d+\\-\\d+"))

# 方法1：使用anti_join() - 更简洁
library(dplyr)
# 只选择需要比较的列
DI_no_subset <- DI_no %>% select(Contig, `Defense Types`)
# 找出DI_withanti中独有的行
DI_withanti_unique <- anti_join(DI_withanti, DI_no_subset, 
                                by = c("Contig", "Defense Types"))

DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")
MAG_metadata <- read.delim2("../01MAG_info/MAG5455_metadata.tsv")
DI_contig <- read.delim2("selected_DI_contig.gff", header = F, sep = '\t')
colnames(DI_contig) <- c("seqid","source","type","start","end",
                                  "score","strand","phase","attributes")
#ORFgff <- ORFgff[c(1:37589),]
library(stringr)
DI_contig <- DI_contig %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(attributes, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(seqid, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列



DI_contig <- left_join(DI_contig,DG, by = "hit_id")

DI_contig <- select(DI_contig,seqid, replicon,strand, hit_id, start, end, gene_name, type.y, subtype, activity)


DI_contig$strand <- ifelse(DI_contig$strand == "-", 0, 1)


DI_contig <- DI_contig %>%
  mutate(orf_number = as.numeric(str_extract(hit_id, "\\d+$"))) %>%
  arrange(seqid, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(seqid) %>%
  mutate(
    is_gene_present = !is.na(gene_name),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = list(which(is_gene_present))
  ) %>%
  mutate(
    # 提取每组的第一个和最后一个基因索引
    first_index = sapply(gene_indices, function(x) x[1]),
    last_index = sapply(gene_indices, function(x) x[length(x)]),
    # 标记要保留的行
    to_keep = (row_number() >= first_index & row_number() <= last_index) | !is.na(gene_name)
  ) %>%
  ungroup() %>%
  filter(to_keep) %>%  # 仅保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列

DI_contig$activity <- ifelse(
  DI_contig$activity %in% c("Defense", "Antidefense"),  # 检查是否是 "Defense" 或 "T4SS"
  DI_contig$activity,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
  )

custom_colors <- c("Defense" = "#7779AE", "Antidefense" = "#E8C07C",
                   "IS" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色

library(ggplot2)
library(gggenes)
library(circlize)
p_final_DI <- ggplot(DI_contig,
                     aes(xmin = start, xmax = end, 
                         y = seqid, fill = activity,
                         forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = DI_contig[!is.na(DI_contig$gene_name),],
    aes(x = start, y = seqid, forward = T,fill = activity)
  ) +
  geom_feature_label(
    data = DI_contig[!is.na(DI_contig$gene_name),],
    aes(x = start, y = seqid, label = gene_name, forward = T,color = activity)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,forward = DI_contig$strand,
                  color = NA,  #控制基因箭头的边框
                  arrowhead_width = unit(1, "mm")) +  #data = 
  scale_fill_manual(values = custom_colors) +  # 使用自定义颜色
  scale_color_manual(values = custom_colors) +  # 添加文字颜色的手动映射
  facet_wrap(~ seqid, scales = "free", ncol = 1) +
  # scale_fill_brewer(palette = "Set3")+  #Others的颜色是"#FDB462"
  # scale_fill_manual(values = c('grey','grey','grey','grey',
  #                              'grey','grey','grey'))+
  theme_genes()+
  theme(axis.title.y = element_blank(),# 隐藏y轴标题
        axis.text.y = element_blank()   # 隐藏y轴刻度标签
  )+
  labs(y = "",x=NULL)+
  guides(fill=guide_legend(title = "Genes",nrow =1,position = 'bottom'))
p_final_DI

write.table(DI_contig, "DI_contig.tsv", sep = "\t", quote = F, row.names = F)

setwd(file.path(GTR_DATA, "20gaoyuan/14_DIs"))
pdf(file = "1_specialDI_genetic_context_DFgff.pdf", width = 6, height = 2.4)
p_final_DI
# 生成绘图代码
dev.off()



save.image("./01_defense_island_characterization.RData")
load("./01_defense_island_characterization.RData")





