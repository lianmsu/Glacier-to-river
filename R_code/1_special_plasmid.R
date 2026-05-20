
setwd("D:/pku/20gaoyuan/15_antidefense/1special_plasmid")

DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")

plasmid <- read.delim2("selected_plasmid_contig.gff",check.names = F,header = F)
colnames(plasmid) <- c("seqid","source","CDS","start","end",
                         "score","strand","phase","attributes")

plasmidhost <- read.delim2("selected_plasmidhost_contig.gff",check.names = F,header = F)
colnames(plasmidhost) <- c("seqid","source","CDS","start","end",
                       "score","strand","phase","attributes")

library(dplyr)
library(stringr)
plasmid <- plasmid %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(attributes, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(seqid, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

plasmidhost <- plasmidhost %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(attributes, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(seqid, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

plasmid <- left_join(plasmid,DG, by = "hit_id")
plasmidhost <- left_join(plasmidhost, DG, by = "hit_id")

plasmid$strand <- ifelse(plasmid$strand == "-", 0, 1)
plasmidhost$strand <- ifelse(plasmidhost$strand == "-", 0, 1)

special_plasmid_bind <- rbind(plasmidhost,plasmid)

special_plasmid_bind <- select(special_plasmid_bind, seqid, start, end,
                               strand, hit_id, gene_name, activity)
# 1. 将protospacer转换为列表并明确指定数据类型
protospacer_row <- list(
  seqid = 'HH-15BT_WAP-bin.49-243',  # 保持字符型
  start = as.integer(16752),         # 转换为整数
  end = as.integer(16781),           # 转换为整数
  strand = as.integer(0),            # 转换为整数
  hit_id = 'HH-15BT_WAP-bin.49-243', # 保持字符型
  gene_name = 'protospacer',         # 保持字符型
  activity = 'protospacer'           # 保持字符型
)
# 2. 转换为数据框并与原数据框合并
special_plasmid_bind <- rbind(
  special_plasmid_bind,
  as.data.frame(protospacer_row, stringsAsFactors = FALSE)
)

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

special_plasmid_bind$activity <- ifelse(
  special_plasmid_bind$activity %in% c("Defense", "Antidefense","protospacer"),  # 检查是否是 "Defense" 或 "T4SS"
  special_plasmid_bind$activity,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
)

custom_colors <- c("Defense" = "#7779AE", "Antidefense" = "#c38e70",
                   "IS" = "#E8C07C", 
                   "protospacer" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色

library(ggplot2)
library(gggenes)
library(circlize)
p_special_plasmid_bind <- ggplot(special_plasmid_bind,
                     aes(xmin = start, xmax = end, 
                         y = seqid, fill = activity,
                         forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = special_plasmid_bind[!is.na(special_plasmid_bind$gene_name),],
    aes(x = start, y = seqid, forward = T,fill = activity)
  ) +
  geom_feature_label(
    data = special_plasmid_bind[!is.na(special_plasmid_bind$gene_name),],
    aes(x = start, y = seqid, label = gene_name, forward = T,color = activity)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1, forward = special_plasmid_bind$strand,
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
p_special_plasmid_bind

setwd("D:/pku/20gaoyuan/15_antidefense/1special_plasmid")
pdf(file = "1_specialplasmid_genetic_context_DFgff.pdf", width = 7.5, height = 5)
p_special_plasmid_bind
# 生成绘图代码
dev.off()

