# =============================================================================
# Glacier-to-river -- antidefensome / Phage antidefensome
# =============================================================================
# Describe the viral catalog (completeness/length, taxonomy) and the defense
# and anti-defense systems encoded by 59,997 phage genomes.
#
# Manuscript      : Fig. S14, Fig. S15
# Main outputs    : 1_p_phage_quality.pdf, 1_phage_tax_sankey.pdf, 1_phage_tax_legend.pdf,
#                   01_phage_antidefensome.RData
# Original script : 20gaoyuan/30allphage_antidefensome/1_allphage_antidefense.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

library(stringr)
library(ggplot2)
library(dplyr)
setwd(file.path(GTR_DATA, "20gaoyuan/30allphage_antidefensome"))

all_59997_phages_DF <- read.delim2("virus_select_all_medium-10kb_high_complete_defense_finder_systems.tsv", 
                                   header = T, check.names = F)
all_59997_phages_DF_anti <- all_59997_phages_DF[all_59997_phages_DF$activity == "Antidefense",]
all_59997_phages_DF_defense <- all_59997_phages_DF[all_59997_phages_DF$activity == "Defense",]

all_59997_phages_DF_anti$contig_id <- str_replace(all_59997_phages_DF_anti$sys_beg, "_\\d+$", "")
all_59997_phages_DF_defense$contig_id <- str_replace(all_59997_phages_DF_defense$sys_beg, "_\\d+$", "")



all_59997_phages_DF_gene <- read.delim2("virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv", 
                                   header = T, check.names = F)
all_59997_phages_DF_gene_anti <- all_59997_phages_DF_gene[all_59997_phages_DF_gene$activity == "Antidefense",]

genomad_res <- read.delim2("virus_select_all_medium-10kb_high_complete_virus_summary_genomad.tsv",
                           header = T, check.names = F)
#对genomad_res 的seq_name列进行清理，删除|及|后的字符，结果作为新列contig

genomad_res$contig <- str_replace(genomad_res$seq_name, "\\|.*", "")

all_55265_phages_DF_anti <- all_59997_phages_DF_anti[all_59997_phages_DF_anti$contig_id %in% genomad_res$contig,]
all_55265_phages_DF_defense <- all_59997_phages_DF_defense[all_59997_phages_DF_defense$contig_id %in% genomad_res$contig,]


#不同antiDS的数量分布####
library(ggplot2)
library(dplyr)  # 用于数据计数

# 1. 准备数据 - 计算每种类型的数量
type_counts <- all_55265_phages_DF_anti %>%
  count(type, name = "count") %>%
  arrange(desc(count))  # 按数量降序排列

# 2. 绘制柱状图
p_antidefense <-
ggplot(type_counts, aes(x = reorder(type, -count), y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_brewer(palette = "Set3") +  # 使用RColorBrewer的Set3调色板
  labs(title = "Distribution of Anti-defense system types",
       x = "",
       y = "Count") +
  geom_text(aes(label = count), vjust = -0.3, size = 3.5) +  # 在柱子上方显示数量
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
    axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_antidefense

#DS分布####
library(RColorBrewer)
# 设置要显示的具体数量
top_n <- 20
DStype_counts <- all_55265_phages_DF_defense %>%
  count(type, name = "count") %>%
  arrange(desc(count)) %>%
  mutate(type = ifelse(row_number() > top_n, "Other", type)) %>%
  group_by(type) %>%
  summarise(count = sum(count)) %>%
  arrange(desc(count))


# 2. 绘制柱状图
p_DS <-
  ggplot(DStype_counts, aes(x = reorder(type, -count), y = count, fill = type)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_manual(values = c(colorRampPalette(brewer.pal(12, "Set3"))(top_n), "gray")) +
  
  labs(title = "Distribution of defense system types",
       x = "",
       y = "Count") +
  geom_text(aes(label = count), vjust = -0.3, size = 3.5) +  # 在柱子上方显示数量
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
    axis.text.x = element_text(angle = 90, size = 10, hjust = 1, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_DS

library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_phage_DS_antiDS_num.pdf" , width = 8,height = 4)
g <- plot_grid(p_DS,p_antidefense,
               labels = letters[3:4], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd(file.path(GTR_DATA, "20gaoyuan/30allphage_antidefensome"))
dev.off()



#phage序列状况、质量状况、物种组成####
head(genomad_res)
checkv_phage_quality <- read.delim2("quality_summary_checkv.tsv")
checkv_phage_quality <- 
  checkv_phage_quality[checkv_phage_quality$contig_id %in% 
                                               genomad_res$contig,]

library(ggplot2)
library(ggbeeswarm)  # 用于更好的散点排列
# 计算每个类别的数量
quality_counts <- checkv_phage_quality %>%
  count(checkv_quality, name = "count") %>%
  mutate(label = paste0(checkv_quality, " (n=", count, ")"))

# 合并计数回原始数据
checkv_phage_quality <- checkv_phage_quality %>%
  left_join(quality_counts, by = "checkv_quality")

# 确保有序因子（包含计数标签）
checkv_phage_quality$label <- factor(
  checkv_phage_quality$label,
  levels = unique(quality_counts$label[order(quality_counts$checkv_quality)]),
  ordered = TRUE
)
# # 确保checkv_quality是有序因子
# checkv_phage_quality$checkv_quality <- factor(
#   checkv_phage_quality$checkv_quality,
#   levels = c("Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete"),
#   ordered = TRUE
# )
# 创建颜色映射
quality_colors <- viridis::viridis_pal(option = "D", 
                                       begin = 0.2, end = 0.8)(length(unique(checkv_phage_quality$checkv_quality)))

# 绘制横向箱线图
p_phage_quality <- 
ggplot(checkv_phage_quality, 
       aes(x = contig_length, 
           y = label, 
           fill = checkv_quality,
           color = checkv_quality)) +
  
  # 添加箱线图
  geom_boxplot(alpha = 0.7, 
               outlier.shape = 1,          # 实心圆点
               outlier.size = 2,            # 离群点大小
               outlier.alpha = 0.7,         # 离群点透明度
               width = 0.6,
               size = 0.4                   # 边框线粗细
               ) +
  # 设置颜色和填充（使用相同颜色方案）
  scale_fill_manual(values = quality_colors) +
  scale_color_manual(values = quality_colors) +
  # 添加散点表示所有数据点
  # geom_quasirandom(
  #   aes(color = checkv_quality),  # 使用颜色区分
  #   alpha = 0.5,                 # 半透明
  #   size = 1.5,                  # 点大小
  #   groupOnX = FALSE             # 在y轴方向排列
  # ) +
  
  # 使用viridis颜色方案
  # scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  # scale_color_viridis_d(option = "D", begin = 0.2, end = 0.8) +
  
  # 对数变换x轴
  scale_x_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  # 
  # 添加标签和主题
  labs(
    title = "Contig length distribution by CheckV quality",
    x = "Contig length (bp, log10 scale)",
    y = "CheckV quality assessment",
    #caption = "Each point represents one contig"
  ) +
  
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
    axis.text.x = element_text(angle = 0, size = 10, hjust = 0.5, vjust = 0.5,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  )
p_phage_quality

setwd(file.path(GTR_DATA, "20gaoyuan/30allphage_antidefensome"))
ggsave("1_p_phage_quality.pdf",p_phage_quality ,width = 4,height = 1.8,units = 'in')


#phage物种组成,桑基图展示####
library(tidyr)
# 拆分 classification 列
genomad_lineage_split <- genomad_res %>%
  separate(
    taxonomy,
    into = c("Virus", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )


library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)

setwd(file.path(GTR_DATA, "20gaoyuan/30allphage_antidefensome"))

df <- genomad_lineage_split

# 统计 k 列中每个k的计数
realm_counts <- df %>%
  count(Realm)

# 统计 phylum 列中每个门的计数
p_counts <- df %>%
  count(Phylum)
# 统计 class 列中每个门的计数
c_counts <- df %>%
  count(Class)

# 统计 order 列中每个门的计数
o_counts <- df %>%
  count(Order)

# 统计 family 列中每个门的计数
f_counts <- df %>%
  count(Family)

df <- df[,c(12:18)]
df <- df[df$Virus != "Unclassified",]

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)
df1 <- df %>% select(1,2) %>% group_by(Virus,Realm) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Virus", "Realm"),
                      values_from = "n")
df2 <- df %>% select(2,4) %>% group_by(Realm,Phylum) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Realm", "Phylum"),
                      values_from = "n")
df3 <- df %>% select(4,5) %>% group_by(Phylum,Class) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Phylum", "Class"),
                      values_from = "n")

df4 <- df %>% select(5,6) %>% group_by(Class,Order) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Class", "Order"),
                      values_from = "n")
df5 <- df %>% select(6,7) %>% group_by(Order,Family) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Order", "Family"),
                      values_from = "n")

p_phage_sankey_4level <- 
  ggplot(data=df1,aes(x = stage,y =n,group = node,
                      edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3,fontface="plain",color="black")+
  # 绘制第 2，3 层级
  geom_sankeyedge(data=df2,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df2,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  geom_text(data=df2 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3,color="black",fontface="plain")+
  # 绘制第 3，4 层级
  geom_sankeyedge(data=df3,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df3,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df3 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  # 绘制第 4，5 层级
  # geom_sankeyedge(data=df4,aes(fill = node),
  #                 position = position_sankey(order = "ascending",v_space ="auto",
  #                                            width = 0.01))+
  # geom_sankeynode(data=df4,aes(fill=node,color=node),
  #                 position = position_sankey(order = "ascending",v_space ="auto",
  #                                            width = 0.05))+
  # geom_text(data=df4 %>% filter(connector=="to"),
  #           aes(label = node),stat = "sankeynode",angle=0,
  #           position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
  #           hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  # 
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "bottom")+
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "in"),
        axis.text.x=element_text(color="black",face="plain",size=10,
                                 margin = margin(b=0.1,unit = "cm")),
        legend.position="none")

p_phage_sankey_4level

setwd(file.path(GTR_DATA, "20gaoyuan/30allphage_antidefensome"))
ggsave("1_phage_tax_sankey.pdf",p_phage_sankey_4level,width = 8,height = 5.3,units = 'in')

write.table(df1,"df1_sankey.tsv",sep = "\t",quote = F, row.names = F)
write.table(checkv_phage_quality,"checkv_phage_quality.tsv",
  sep = "\t",quote = F, row.names = F)
write.table(DStype_counts,"DStype_counts.tsv",
  sep = "\t",quote = F, row.names = F)
write.table(type_counts,"type_counts_antidefense.tsv",
  sep = "\t",quote = F, row.names = F)


#展示不同region ecosystem的病毒物种组成####
head(genomad_lineage_split)

checkv_phage_quality_2 <- 
  checkv_phage_quality %>%
  left_join(genomad_lineage_split, by = c("contig_id"="seq_name"))

checkv_phage_quality_2 <- checkv_phage_quality_2 %>%
  separate(contig_id, into = c("seq_number", "Sample_ID", "River"), 
           sep = "\\+", remove = FALSE)

EF_info <- read.delim2("info_corrected.txt",sep = "\t")

EF_info$Plateau <- ifelse(EF_info$Plateau == 'Yes',"Plateau",EF_info$Plateau)
EF_info$Plateau <- ifelse(EF_info$Plateau == 'No',"Plain",EF_info$Plateau)
EF_info$Plateau <- ifelse(EF_info$Plateau == 'glacier',"Glacier",EF_info$Plateau)
colnames(EF_info)[7] <- 'Region'

EF_info$medium <- ifelse(EF_info$Type == 'Cryoconite',"Glacier",EF_info$Type)
EF_info$medium <- ifelse(EF_info$medium == 'Ice',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'Snow',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'particle',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'water',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'sediment',"Benthon",EF_info$medium)


EF_info$ecosystem <- EF_info$River
EF_info$ecosystem <- str_remove(EF_info$ecosystem, "_.*")

table(EF_info$ecosystem)
table(EF_info$Region)
table(EF_info$Type)
table(EF_info$medium)

checkv_phage_quality_2 <- 
  checkv_phage_quality_2 %>%
  left_join(EF_info, by = c("Sample_ID"="SampleID"))


phage_lineage_geography <- select(checkv_phage_quality_2, contig_id, 
                                Class, ecosystem, Region, Type, medium)

phage_lineage_geography <- 
  phage_lineage_geography[!is.na(phage_lineage_geography$Region),]


  


phyla_count <- table(phage_lineage_geography$Class) %>% as.data.frame()
phyla_count_merged <- phyla_count %>%
  mutate(Var1 = ifelse(Freq < 10, "Others", as.character(Var1))) %>%
  group_by(Var1) %>%
  summarise(Freq = sum(Freq)) %>%
  arrange(desc(Freq))  # 按频率降序排列
phyla_count_merged$Var1[3] <- "Unassigned"
phyla_count_merged$color <- c(
  "#00a8e8", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#CCCCCC", "#FF9896"
)
# 确保Others是灰色
phyla_count_merged$color[phyla_count_merged$Var1 == "Others"] <- "#CCCCCC"

phage_lineage_geography$Class_2 <- ifelse(
  phage_lineage_geography$Class == "",          # 先判断是否为空白字符串
  "Unassigned",                                  # 空白 → Unassigned
  ifelse(
    phage_lineage_geography$Class %in% phyla_count_merged$Var1,
    phage_lineage_geography$Class,               # 匹配 → 保留原值
    "Others"                                      # 不匹配且非空白 → Others
  )
)

phage_lineage_geography <- 
  phage_lineage_geography[!is.na(phage_lineage_geography$Class_2),]

# phage_lineage_geography$Class_2 <- 
#   ifelse(phage_lineage_geography$Class %in% phyla_count_merged$Var1,
#          phage_lineage_geography$Class, "Others")

phage_lineage_geography <- merge(phage_lineage_geography, phyla_count_merged,
                               by.x = "Class_2", by.y = "Var1")


library(ggplot2)

phage_lineage_geography$ecosystem <- 
  factor(phage_lineage_geography$ecosystem,
         levels = c("Glacier","CJ","YR","LCJ","YJ"))

p_ecosystem <- 
  ggplot(phage_lineage_geography, aes(x = ecosystem, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = phage_lineage_geography$Class_2,
                      breaks = phage_lineage_geography$color) +
  labs(x = "Ecosystem", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_ecosystem

phage_lineage_geography$Region <- factor(phage_lineage_geography$Region, 
                                       levels = c("Glacier","Plateau","Plain"))
p_Region <- 
  ggplot(phage_lineage_geography, aes(x = Region, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = phage_lineage_geography$phyla,
                      breaks = phage_lineage_geography$color) +
  labs(x = "Region", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_Region


phage_lineage_geography$medium <- factor(phage_lineage_geography$medium, 
                                       levels = c("Glacier","Benthon","Plankton"))
p_medium <- 
  ggplot(phage_lineage_geography, aes(x = medium, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = phage_lineage_geography$phyla,
                      breaks = phage_lineage_geography$color) +
  labs(x = "Medium", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_medium

phage_lineage_geography$Type <- factor(phage_lineage_geography$Type, 
                                     levels = c("Snow","Ice","Cryoconite",
                                                "sediment","particle","water"))
p_Type <- 
  ggplot(phage_lineage_geography, aes(x = Type, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = phage_lineage_geography$phyla,
                      breaks = phage_lineage_geography$color) +
  labs(x = "Detailed medium", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_Type


library(Cairo)
library(cowplot)

row2 <- plot_grid(p_Region, 
                  p_medium, 
                  nrow = 1,ncol = 2,
                  labels = letters[2:3]
                  # rel_heights = c(1, 1, 1.5)
)
row2
final_plot <- plot_grid(p_ecosystem, 
                        row2, 
                        p_Type, 
                        nrow = 3,ncol = 1,
                        labels = letters[1:3]
                        # rel_heights = c(1, 1, 1.5)
)
final_plot


cairo_pdf(width = 5.6,height = 7, 
          filename = "1_phage_tax_composition.pdf")
final_plot
dev.off()

p_legend <-   
  ggplot(phage_lineage_geography, aes(x = Type, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = phage_lineage_geography$Class_2,
                      breaks = phage_lineage_geography$color)
p_legend

ggplot2::ggsave("1_phage_tax_legend.pdf",p_legend, 
                width = 5.6, height = 7, units = "in")


save.image("./01_phage_antidefensome.RData")
load("./01_phage_antidefensome.RData")

write.table(phage_lineage_geography, "phage_lineage_geography.tsv",sep = "\t",
  quote = F, row.names = F)

