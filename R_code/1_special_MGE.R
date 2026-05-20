
setwd("D:/pku/20gaoyuan/13_special_MGE")

DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")
ORFgff <- read.delim2("DG_gff.gff", header = F, sep = '\t')
MAG_metadata <- read.delim2("MAG5455_metadata.tsv")

ORFgff <- ORFgff[c(1:37492),]
library(stringr)
ORFgff <- ORFgff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

DG <- left_join(DG, ORFgff, by = "hit_id")
DG <- select(DG, replicon, hit_id, V4, V5, gene_name, type, subtype, activity)

DG <- DG %>% 
  mutate(genome = str_extract(replicon, ".*bin\\.\\d+"))
DG <- DG[DG$activity == "Defense",]
DG <- DG[DG$genome %in% MAG_metadata$genome,]

#ICE####
setwd("D:/pku/20gaoyuan/13_special_MGE/1ICE")

contig1 <- read.delim2("CJ_14-30WL_WA-bin.46-17.gff",header = F)
colnames(contig1) <- c("seqid","source","type","start","end",
                       "score","strand","phase","attributes")
contig1_ICE <- read.delim2("CJ_14-30WL_WA-bin.46_contig_17_ICE2_gene.gff",header = F)
contig1_ICE <- contig1_ICE[,c(1:9)]
colnames(contig1_ICE) <- c("seqid","source","type","start","end",
                       "score","strand","phase","attributes")
contig1_ICE$strand <- gsub(",","",contig1_ICE$strand)
contig1_ICE$seqid <- contig1$seqid[1]
# 
contig1_ICE$genename <- gsub(".*Name=([^;]+).*", "\\1", contig1_ICE$attributes)
# 修改class列的提取逻辑，处理Note=后无内容的情况
contig1_ICE$class <- ifelse(
  grepl("Note=", contig1_ICE$attributes),
  gsub(".*Note=([^;]*).*", "\\1", contig1_ICE$attributes),
  NA
)
# 将空字符串转换为NA
contig1_ICE$class[contig1_ICE$class == ""] <- NA



contig2 <- read.delim2("Glacier-1611DDS-bin.12-3.gff",header = F)
colnames(contig2) <- c("seqid","source","type","start","end",
                       "score","strand","phase","attributes")
contig2_ICE <- read.delim2("Glacier-1611DDS-bin.12_contig_3_ICE2_gene.gff",header = F)
contig2_ICE <- contig2_ICE[,c(1:9)]
colnames(contig2_ICE) <- c("seqid","source","type","start","end",
                           "score","strand","phase","attributes")
contig2_ICE$strand <- gsub(",","",contig2_ICE$strand)
contig2_ICE$seqid <- contig2$seqid[1]
# 
contig2_ICE$genename <- gsub(".*Name=([^;]+).*", "\\1", contig2_ICE$attributes)
# 修改class列的提取逻辑，处理Note=后无内容的情况
contig2_ICE$class <- ifelse(
  grepl("Note=", contig2_ICE$attributes),
  gsub(".*Note=([^;]*).*", "\\1", contig2_ICE$attributes),
  NA
)
# 将空字符串转换为NA
contig2_ICE$class[contig2_ICE$class == ""] <- NA

contig <- rbind(contig1, contig2)
contig_ICE <- rbind(contig1_ICE, contig2_ICE)

#contig <- left_join(contig, contig_ICE, by = c("start","end"))

# 为contig数据框新增genename和class列，初始值为NA
contig$genename <- NA
contig$class <- NA

# 遍历contig中的每一行
for(i in 1:nrow(contig)) {
  # 获取当前CDS的起始和结束位置
  cds_start <- contig$start[i]
  cds_end <- contig$end[i]
  
  # 在contig_ICE中查找包含当前CDS的基因区间
  matching_gene <- which(
    contig_ICE$start <= cds_start & 
      contig_ICE$end >= cds_end &
      contig_ICE$seqid == contig$seqid[i]
  )
  
  # 如果找到匹配的基因区间，则提取genename和class
  if(length(matching_gene) > 0) {
    # 如果有多个匹配，取第一个
    matching_gene <- matching_gene[1]
    contig$genename[i] <- contig_ICE$genename[matching_gene]
    contig$class[i] <- contig_ICE$class[matching_gene]
  }
}

# 遍历contig中的每一行
for(i in 1:nrow(contig)) {
  # 获取当前CDS的起始和结束位置
  cds_start <- contig$start[i]
  cds_end <- contig$end[i]
  
  # 在contig_ICE中查找包含当前CDS的基因区间
  matching_gene <- which(
    DG$V4 <= cds_start & 
      DG$V5 >= cds_end &
      DG$replicon == contig$seqid[i]
  )
  
  # 如果找到匹配的基因区间，则提取genename和class
  if(length(matching_gene) > 0) {
    # 如果有多个匹配，取第一个
    matching_gene <- matching_gene[1]
    contig$genename[i] <- DG$gene_name[matching_gene]
    contig$class[i] <- DG$activity[matching_gene]
  }
}



library(ggplot2)
library(gggenes)
library(circlize)

# 提取V9列中的ID部分
contig <- contig %>%
  mutate(ID_part = str_extract(attributes, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(seqid, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

contig$genename <- ifelse(
  contig$genename == "hypothetical protein",  # 检查是否是 "Defense" 或 "T4SS"
  NA,                            # 如果是，保持原值
  contig$genename                                      # 否则改为 "Others"
)

final_contig <- contig %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(seqid, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(seqid) %>%
  mutate(
    is_gene_present = !is.na(genename),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = list(which(is_gene_present))
  ) %>%
  mutate(
    # 提取每组的第一个和最后一个基因索引
    first_index = sapply(gene_indices, function(x) x[1]),
    last_index = sapply(gene_indices, function(x) x[length(x)]),
    # 标记要保留的行
    to_keep = (row_number() >= first_index & row_number() <= last_index) | !is.na(genename)
  ) %>%
  ungroup() %>%
  filter(to_keep) %>%  # 仅保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列

final_contig$class <- ifelse(
  final_contig$class %in% c("Defense", "T4SS"),  # 检查是否是 "Defense" 或 "T4SS"
  final_contig$class,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
)

final_contig$class <- factor(final_contig$class, levels = c("Defense",
                                                  "T4SS","Others"
))

final_contig$strand <- ifelse(final_contig$strand == "-", 0, 1)

custom_colors <- c("Defense" = "#7779AE", 
                   "T4SS" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色


p_final_contig <- ggplot(final_contig,
                         aes(xmin = start, xmax = end, 
                             y = seqid, fill = class,
                             forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_contig[!is.na(final_contig$genename),],
    aes(x = start, y = seqid, forward = T,fill = class)
  ) +
  geom_feature_label(
    data = final_contig[!is.na(final_contig$genename),],
    aes(x = start, y = seqid, label = genename, forward = T,color = class)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,forward = final_contig$strand,
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
p_final_contig

setwd("D:/pku/20gaoyuan/13_special_MGE/1ICE")
#ggsave("p_final_contig_2ICE.pdf",p_final_contig,width = 8,height = 3,units = 'in',
#      limitsize = FALSE)


#IS####
library(dplyr)
library(stringr)

setwd("D:/pku/20gaoyuan/13_special_MGE")

DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")
#ORFgff <- read.delim2("DG_gff.gff", header = F, sep = '\t')
ORFgff <- read.delim2("DG_gff_DF.gff", header = F, sep = '\t')
MAG_metadata <- read.delim2("MAG5455_metadata.tsv")

ORFgff <- ORFgff[c(1:37589),]
library(stringr)
ORFgff <- ORFgff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

DG <- left_join(DG, ORFgff, by = "hit_id")
DG <- select(DG, replicon, hit_id, V4, V5, gene_name, type, subtype, activity)

DG <- DG %>% 
  mutate(genome = str_extract(replicon, ".*bin\\.\\d+"))
DG <- DG[DG$activity == "Defense",]
DG <- DG[DG$genome %in% MAG_metadata$genome,]


setwd("D:/pku/20gaoyuan/13_special_MGE/2IS")
selected_IS <- read.delim2("IS_DS_filtered_mix.tsv")

selected_IS_contig <- read.delim2("selected_IS_DS_contig.gff",header = F)
colnames(selected_IS_contig) <- c("seqid","source","type","start","end",
                       "score","strand","phase","attributes")

selected_IS_contig <- selected_IS_contig %>%
  mutate(ID_part = str_extract(attributes, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+")) %>%
  mutate(ORF = paste0(seqid, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

selected_IS_contig$strand <- ifelse(selected_IS_contig$strand == "-", 0, 1)

selected_IS_contig$genename <- NA
selected_IS_contig$class <- NA
# 遍历selected_IS_contig中的每一行
for(i in 1:nrow(selected_IS_contig)) {
  # 获取当前CDS的起始和结束位置
  cds_start <- selected_IS_contig$start[i]
  cds_end <- selected_IS_contig$end[i]
  
  # 在contig_ICE中查找包含当前CDS的基因区间
  matching_gene <- which(
    DG$V4 <= cds_start & 
      DG$V5 >= cds_end &
      DG$replicon == selected_IS_contig$seqid[i]
  )
  
  # 如果找到匹配的基因区间，则提取genename和class
  if(length(matching_gene) > 0) {
    # 如果有多个匹配，取第一个
    matching_gene <- matching_gene[1]
    selected_IS_contig$genename[i] <- DG$gene_name[matching_gene]
    selected_IS_contig$class[i] <- DG$activity[matching_gene]
  }
}


selected_IS_2 <- data.frame(seqid = selected_IS$seq_name,
                            source = "selected_IS",
                            type = selected_IS$family,
                            start = selected_IS$isBegin,
                            end = selected_IS$isEnd,
                            score = selected_IS$score,
                            strand = selected_IS$strand,
                            phase = selected_IS$nGaps,
                            attributes = selected_IS$tir,
                            ORF = selected_IS$sys_end,
                            genename = selected_IS$family,
                            class = "IS"
                            )
selected_IS_2$strand <- ifelse(selected_IS_2$strand == "-", 0, 1)

  
final_IS <- rbind(selected_IS_contig, selected_IS_2)


final_IS <- final_IS %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(seqid, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(seqid) %>%
  mutate(
    is_gene_present = !is.na(genename),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = list(which(is_gene_present))
  ) %>%
  mutate(
    # 提取每组的第一个和最后一个基因索引
    first_index = sapply(gene_indices, function(x) x[1]),
    last_index = sapply(gene_indices, function(x) x[length(x)]),
    # 标记要保留的行
    to_keep = (row_number() >= first_index & row_number() <= last_index) | !is.na(genename)
  ) %>%
  ungroup() %>%
  filter(to_keep) %>%  # 仅保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列










final_IS$class <- ifelse(
  final_IS$class %in% c("Defense", "IS"),  # 检查是否是 "Defense" 或 "T4SS"
  final_IS$class,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
  )

custom_colors <- c("Defense" = "#7779AE", 
                   "IS" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色

library(ggplot2)
library(gggenes)
library(circlize)
p_final_IS <- ggplot(final_IS,
                         aes(xmin = start, xmax = end, 
                             y = seqid, fill = class,
                             forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_IS[!is.na(final_IS$genename),],
    aes(x = start, y = seqid, forward = T,fill = class)
  ) +
  geom_feature_label(
    data = final_IS[!is.na(final_IS$genename),],
    aes(x = start, y = seqid, label = genename, forward = T,color = class)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,forward = final_IS$strand,
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
p_final_IS

setwd("D:/pku/20gaoyuan/13_special_MGE/2IS")
#ggsave("p_final_contig_2ICE.pdf",p_final_contig,width = 8,height = 3,units = 'in',
#      limitsize = FALSE)


pdf(file = "IS_plot_DFgff.pdf", width = 8, height = 15)
p_final_IS
# 生成绘图代码
dev.off()


#prophage####
library(dplyr)
library(stringr)

setwd("D:/pku/20gaoyuan/13_special_MGE")

DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")
#ORFgff <- read.delim2("DG_gff.gff", header = F, sep = '\t')
ORFgff <- read.delim2("DG_gff_DF.gff", header = F, sep = '\t')
MAG_metadata <- read.delim2("MAG5455_metadata.tsv")

ORFgff <- ORFgff[c(1:37589),]
library(stringr)
ORFgff <- ORFgff %>%
  mutate(
    # 提取ID=后面的内容，然后取最后一个_后的数字
    gene_num = str_extract(V9, "ID=[^;]+") %>% 
      str_remove("ID=") %>% 
      str_extract("[^_]+$"),
    hit_id = paste0(V1, "_", gene_num)
  ) %>%
  select(-gene_num)  # 移除临时列

DG <- left_join(DG, ORFgff, by = "hit_id")
DG <- select(DG, replicon, hit_id, V4, V5, gene_name, type, subtype, activity)

DG <- DG %>% 
  mutate(genome = str_extract(replicon, ".*bin\\.\\d+"))
DG <- DG[DG$activity == "Defense",]
DG <- DG[DG$genome %in% MAG_metadata$genome,]


setwd("D:/pku/20gaoyuan/13_special_MGE/3prophage")
selected_prophage <- read.delim2("selected_prophage.tsv")

selected_prophage_contig <- read.delim2("selected_prophage_DS_contig.gff",header = F)
colnames(selected_prophage_contig) <- c("seqid","source","type","start","end",
                                  "score","strand","phase","attributes")

selected_prophage_contig <- selected_prophage_contig %>%
  mutate(ID_part = str_extract(attributes, "ID=[^;]+")) %>%
  mutate(ID_suffix = str_extract(ID_part, "_\\d+$")) %>%
  mutate(ORF = paste0(seqid, ID_suffix)) %>%
  select(-ID_part, -ID_suffix)

selected_prophage_contig$strand <- ifelse(selected_prophage_contig$strand == "-", 0, 1)

selected_prophage_contig$genename <- NA
selected_prophage_contig$class <- NA
selected_prophage_contig <- 
  selected_prophage_contig[selected_prophage_contig$seqid %in% selected_prophage$contig,]

# 遍历selected_IS_contig中的每一行
for(i in 1:nrow(selected_prophage_contig)) {
  # 获取当前CDS的起始和结束位置
  cds_start <- selected_prophage_contig$start[i]
  cds_end <- selected_prophage_contig$end[i]
  
  # 在contig_ICE中查找包含当前CDS的基因区间
  matching_gene <- which(
    DG$V4 <= cds_start & 
      DG$V5 >= cds_end &
      DG$replicon == selected_prophage_contig$seqid[i]
  )
  
  # 如果找到匹配的基因区间，则提取genename和class
  if(length(matching_gene) > 0) {
    # 如果有多个匹配，取第一个
    matching_gene <- matching_gene[1]
    selected_prophage_contig$genename[i] <- DG$gene_name[matching_gene]
    selected_prophage_contig$class[i] <- DG$activity[matching_gene]
  }
}


selected_prophage_2 <- data.frame(seqid = selected_prophage$contig,
                            source = "selected_prophage",
                            type = selected_prophage$MGE,
                            start = selected_prophage$begin,
                            end = selected_prophage$end,
                            score = 0,
                            strand = '-',
                            phase = 0,
                            attributes = selected_prophage$position,
                            ORF = selected_prophage$contig,
                            genename = selected_prophage$MGE,
                            class = "MGE"
)
selected_prophage_2$strand <- ifelse(selected_prophage_2$strand == "-", 0, 1)

final_1 <- selected_prophage_contig %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(seqid, orf_number) %>%  # 按 V1 和末尾数字排序
  group_by(seqid) %>%
  mutate(
    is_gene_present = !is.na(genename),
    # 获取 gene 不为 NA 的行的索引
    gene_indices = list(which(is_gene_present))
  ) %>%
  mutate(
    # 提取每组的第一个和最后一个基因索引
    first_index = sapply(gene_indices, function(x) x[1]),
    last_index = sapply(gene_indices, function(x) x[length(x)]),
    # 标记要保留的行
    to_keep = (row_number() >= first_index & row_number() <= last_index) | !is.na(genename)
  ) %>%
  ungroup() %>%
  filter(to_keep) %>%  # 仅保留需要的行
  select(-orf_number, -is_gene_present, -gene_indices, -first_index, -last_index, -to_keep)  # 清理临时列


final_prophage$class <- ifelse(
  final_prophage$class %in% c("Defense", "MGE"),  # 检查是否是 "Defense" 或 "T4SS"
  final_prophage$class,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
)

final_prophage <- rbind(selected_prophage_contig, selected_prophage_2)

custom_colors <- c("Defense" = "#7779AE", 
                   "MGE" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色

library(ggplot2)
library(gggenes)
library(circlize)
p_final_prophage <- ggplot(final_prophage,
                     aes(xmin = start, xmax = end, 
                         y = seqid, fill = class,
                         forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_prophage[!is.na(final_prophage$genename),],
    aes(x = start, y = seqid, forward = T,fill = class)
  ) +
  geom_feature_label(
    data = final_prophage[!is.na(final_prophage$genename),],
    aes(x = start, y = seqid, label = genename, forward = T,color = class)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,forward = final_prophage$strand,
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
p_final_prophage

setwd("D:/pku/20gaoyuan/13_special_MGE/3prophage")
#ggsave("p_final_contig_2ICE.pdf",p_final_contig,width = 8,height = 3,units = 'in',
#      limitsize = FALSE)
pdf(file = "prophage_plot_DFgff.pdf", width = 8, height = 8)
p_final_prophage
# 生成绘图代码
dev.off()


#观察发现HH-36XT_WAP-bin.37-7保留25到90号基因，YJ-RKZ_Ap-bin.30-6保留21到49
selected_prophage_contig <- selected_prophage_contig[c(25:90,148:176),]


final_1 <- selected_prophage_contig %>%
  mutate(orf_number = as.numeric(str_extract(ORF, "\\d+$"))) %>%
  arrange(seqid, orf_number) %>%  # 按 V1 和末尾数字排序
  select(-orf_number)  # 清理临时列
final_prophage <- rbind(selected_prophage_contig, selected_prophage_2)
final_prophage$class <- ifelse(
  final_prophage$class %in% c("Defense", "MGE"),  # 检查是否是 "Defense" 或 "T4SS"
  final_prophage$class,                            # 如果是，保持原值
  "Others"                                       # 否则改为 "Others"
)

custom_colors <- c("Defense" = "#7779AE", 
                   "MGE" = "#E8C07C", 
                   "type3" = "green", 
                   "Others" = "#91C8ED")  # 包含 "others" 的自定义颜色



prophage_maker_annotation <- read.delim2("final_prophage_marker.tsv",sep = "\t")
prophage_maker_annotation <- prophage_maker_annotation[!is.na(prophage_maker_annotation$maker_annotation),]

for(i in 1:nrow(final_prophage)) {
  # 获取当前CDS的起始和结束位置
  cds_start <- final_prophage$start[i]
  cds_end <- final_prophage$end[i]
  
  # 在contig_ICE中查找包含当前CDS的基因区间
  matching_gene <- which(
    prophage_maker_annotation$start == cds_start & 
      prophage_maker_annotation$end == cds_end 
  )
  
  # 如果找到匹配的基因区间，则提取genename和class
  if(length(matching_gene) > 0) {
    # 如果有多个匹配，取第一个
    matching_gene <- matching_gene[1]
    final_prophage$genename[i] <- prophage_maker_annotation$maker_annotation[matching_gene]
    #final_prophage$class[i] <- prophage_maker_annotation$activity[matching_gene]
  }
}


library(ggplot2)
library(gggenes)
library(circlize)
p_final_prophage <- ggplot(final_prophage,
                           aes(xmin = start, xmax = end, 
                               y = seqid, fill = class,
                               forward = strand)) +
  # scale_fill_manual(values = c("grey"),  # 为特定的分类指定颜色
  #                   limits = c("others")  # 设置图例的顺序和分类
  #                   )+
  geom_feature(
    data = final_prophage[!is.na(final_prophage$genename),],
    aes(x = start, y = seqid, forward = T,fill = class)
  ) +
  geom_feature_label(
    data = final_prophage[!is.na(final_prophage$genename),],
    aes(x = start, y = seqid, label = genename, forward = T,color = class)
  ) + 
  geom_gene_arrow(arrowhead_height = unit(3, "mm"),
                  size = 0.1,forward = final_prophage$strand,
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
p_final_prophage


setwd("D:/pku/20gaoyuan/13_special_MGE/3prophage")
#ggsave("p_final_contig_2ICE.pdf",p_final_contig,width = 8,height = 3,units = 'in',
#      limitsize = FALSE)
pdf(file = "prophage_plot_DFgff_checkVfilter.pdf", width = 8, height = 5)
p_final_prophage
# 生成绘图代码
dev.off()














MAG_metadata[MAG_metadata$genome == "CJ_18-PZH_SW-bin.10",] #Plateau
MAG_metadata[MAG_metadata$genome == "YJ-RKZ_Ap-bin.30",] #Plateau
MAG_metadata[MAG_metadata$genome == "HH-36XT_WAP-bin.37",] #Plain
MAG_metadata[MAG_metadata$genome == "HH-TK_w22-bin.45",] #Plateau
MAG_metadata[MAG_metadata$genome == "LCJ-NZD_pA-bin.87",] #Plain



#CRISPR-Cas9####

setwd("D:/pku/20gaoyuan/13_special_MGE/4CRISPR-Cas9")
all_DS <- read.delim2("merged_5467MAG_defense_finder_systems.tsv")
all_DS <- all_DS %>% 
  mutate(genome = str_extract(sys_beg, ".*bin\\.\\d+"),
         contig = str_extract(sys_beg, ".*bin\\.\\d+\\-\\d+"),)
all_DS <- all_DS[all_DS$activity == "Defense",]
all_DS <- all_DS[all_DS$genome %in% MAG_metadata$genome,]

all_DS <- all_DS[all_DS$type == "Cas",]


#统计cas9附近存在共定位蛋白的比例
library(stringr)
all_cas9 <- all_DS %>% 
  filter(str_detect(name_of_profiles_in_sys, regex("cas9", ignore_case = TRUE)))

all_cas9 <- all_cas9 %>% 
    mutate(gene_range = as.numeric(str_extract(all_cas9$sys_end, "\\d+$"))  - 
             as.numeric(str_extract(all_cas9$sys_beg, "\\d+$")) + 1)

sum(all_cas9$genes_count == all_cas9$gene_range) #96
sum(all_cas9$genes_count != all_cas9$gene_range) #40

cas9_statistic <- data.frame(cas9_type = c("No colocalized gene", 
                                           "With colocalized gene"),
                             number = c(96, 40))
#huitu
library(ggplot2)
p_cas9gene_number <- 
  ggplot(cas9_statistic, aes(x = reorder(cas9_type, -number), y = number, fill = cas9_type)) +
  geom_bar(stat = "identity", width = 0.5, color = "white", alpha = 0.8) +  # 柱状图
  geom_text(aes(label = number), vjust = -0.5, size = 3.5) +               # 添加数值标签
  scale_fill_manual(values = c("No colocalized gene" = "#619CFF", "With colocalized gene"="#fa9500",
                               "Integrons"="#e8a598", "ICEs"="#00BA38",
                               "Chromosome" = "#FCDABA",
                               "ISs" = "#9f7e69",
                               "Phage" = "#ED6764")) +                                    # 使用Set2配色
  labs(
   # title = "",
    x = "CRISPR-Cas9 system",
    y = "Count",
    fill = ""
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
    axis.text.x = element_text(angle = 8, size = 10, hjust = 0.5, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))                 # 调整Y轴范围
p_cas9gene_number
setwd("D:/pku/20gaoyuan/13_special_MGE/4CRISPR-Cas9")
ggsave("1_p_cas9gene_number.pdf",p_cas9gene_number,width = 2.6,height = 2.6,units = "in")


#提取cas9附近的共定位蛋白
all_cas9_2 <- all_cas9[all_cas9$genes_count != all_cas9$gene_range,]
library(tidyverse)
# 首先创建一个函数来提取范围内的基因ID
extract_genes <- function(row) {
  # 获取当前行的相关信息
  genome_id <- row$contig
  sys_beg_num <- as.numeric(str_extract(row$sys_beg, "(?<=_)\\d+$"))
  sys_end_num <- as.numeric(str_extract(row$sys_end, "(?<=_)\\d+$"))
  
  # 提取protein_in_syst中的基因ID
  excluded_genes <- str_split(row$protein_in_syst, ",") %>% unlist()
  
  # 生成sys_beg到sys_end之间的所有可能基因ID
  all_possible_genes <- paste0(genome_id, "_", seq(sys_beg_num, sys_end_num))
  
  # 返回不在excluded_genes中的基因
  setdiff(all_possible_genes, excluded_genes)
}

# 应用函数到每一行
all_cas9_2 <- all_cas9_2 %>%
  rowwise() %>%
  mutate(
    genes_in_range = list(extract_genes(cur_data())),
    genes_in_range_count = length(genes_in_range),
    genes_in_range_str = ifelse(length(genes_in_range) > 0, 
                                paste(genes_in_range, collapse = ","), 
                                NA)
  ) %>%
  ungroup()

# 查看结果
head(all_cas9_2 %>% select(sys_id, genes_in_range_str))

# 1. 将所有逗号分隔的字符串合并并拆分为单个基因ID
all_genes <- unlist(strsplit(all_cas9_2$genes_in_range_str, ","))
# 2. 写入文本文件（每行一个基因ID）
writeLines(all_genes, "cas9_genes_in_range.txt")


#蛋白长度分布
protein_length <- read.delim2("sequence_lengths.tsv",header = F)
colnames(protein_length) <- c("protein","length")

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

protein_length_plot <-
  ggplot(protein_length, aes(x=length)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue")+
  labs(
    #    title = "Density Plot of MAG Completeness",
    x = "Protein length",
    y = "Density"
  ) +
  mytheme
protein_length_plot

protein_55_eggnog <- read.delim2("MM_igenas9b.emapper.annotations_2.tsv")
# protein_55_eggnog <- merge(protein_55_eggnog, protein_length, 
#                            by.x = "query",by.y = "protein")
protein_55_eggnog <- protein_55_eggnog[protein_55_eggnog$COG_category != "-",]

#蛋白注释统计
protein_annotation <- merge(protein_length,protein_55_eggnog,
                            by.x = "protein",by.y = "query", all = T)
protein_annotation[is.na(protein_annotation)] <- "Unknown"
protein_annotation_stat <- table(protein_annotation$PFAMs) %>% as.data.frame()
colnames(protein_annotation_stat) <- c("Annotation","Number")

p_protein_annotation <- 
  ggplot(protein_annotation_stat, aes(x = reorder(Annotation, Number), y = Number, fill = Annotation)) +
  geom_bar(stat = "identity", width = 0.5, color = "white", alpha = 0.8) +  # 柱状图
  geom_text(aes(label = Number), vjust = 0.5, hjust =0,size = 3.5) +               # 添加数值标签
  scale_fill_manual(values = c("No colocalized gene" = "#619CFF", "With colocalized gene"="#fa9500",
                               "Integrons"="#e8a598", "ICEs"="#00BA38",
                               "Chromosome" = "#FCDABA",
                               "ISs" = "#9f7e69",
                               "Unknown" = "#ED6764")) +                                    # 使用Set2配色
  labs(
    # title = "",
    x = "",
    y = "Count",
    fill = ""
  ) +
  coord_flip() +  # 关键：转换为横向
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
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))                 # 调整Y轴范围
p_protein_annotation
setwd("D:/pku/20gaoyuan/13_special_MGE/4CRISPR-Cas9")
ggsave("1_p_protein_annotation.pdf",p_protein_annotation,width = 2.6,height = 2.6,units = "in")


library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_cas9_colocationprotein.pdf" , width = 6.8,height = 5.5)
g1 <- plot_grid(p_cas9gene_number,
                p_protein_annotation,
                protein_length_plot,
                labels = letters[1:3], ncol = 2
                #rel_heights = c(1, 1, 1.5),
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1
setwd("D:/pku/20gaoyuan/13_special_MGE/4CRISPR-Cas9")
dev.off()




