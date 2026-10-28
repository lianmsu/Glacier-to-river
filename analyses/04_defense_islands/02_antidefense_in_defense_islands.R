# =============================================================================
# Glacier-to-river -- defense islands / Anti-defense genes inside defense islands
# =============================================================================
# Map anti-defense genes onto the defense islands and describe the
# representative DIs that carry both defense and anti-defense systems.
#
# Manuscript      : Fig. 5e
# Main outputs    : anti-defense co-localisation tables
# Original script : 20gaoyuan/15_antidefense/4antidefense_colocalization/1_anti_colocalization.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "20gaoyuan/15_antidefense/4antidefense_colocalization"))

DG <- read.delim2("3contig.id_anti.tsv",header = T)
#phageDG <- read.delim2("virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv")

phage <- read.delim2("3contig.id.gff",check.names = F,header = F)
colnames(phage) <- c("seqid","source","CDS","start","end",
                         "score","strand","phase","attributes")


library(dplyr)
library(stringr)
phage <- phage %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(attributes, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(seqid, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

phage <- left_join(phage, DG, by = "hit_id")

phage$strand <- ifelse(phage$strand == "-", 0, 1)


phage2 <- select(phage, seqid, start, end,
                               strand, hit_id, gene_name, type, activity)

phage2 <- phage2 %>%
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

#special_phage_bind <- special_phage_bind[c(-1:-168,-237:-264),]

phage2$type <- ifelse(
  !is.na(phage2$type) ,  # 检查是否是 "Defense" 或 "T4SS"
  phage2$type,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
)

custom_colors <- c("Anti_CBASS" = "#7779AE", "NADP" = "#c38e70",
                   "Anti_Retron" = "#E8C07C", 
                   "Anti_Thoeris" = "#52b2cf", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色

library(ggplot2)
library(gggenes)
library(circlize)
p_special_phage2 <- ggplot(phage2,
                     aes(xmin = start, xmax = end, 
                         y = seqid, fill = type,
                         forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = phage2[!is.na(phage2$gene_name),],
    aes(x = start, y = seqid, forward = T,fill = type)
  ) +
  geom_feature_label(
    data = phage2[!is.na(phage2$gene_name),],
    aes(x = start, y = seqid, label = gene_name, forward = T,color = type)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1, forward = phage2$strand,
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
p_special_phage2

setwd(file.path(GTR_DATA, "20gaoyuan/15_antidefense/4antidefense_colocalization"))
pdf(file = "1_specialphage_4anti_DFgff.pdf", width = 4, height = 3)
p_special_phage2
# 生成绘图代码
dev.off()





