# =============================================================================
# Glacier-to-river -- mobile genetic elements / MGEs and the defensome they carry
# =============================================================================
# Count prophages, plasmids, ICE/IME and IS per MAG, compare DS density on
# chromosomes vs MGEs, and compute the O/E localisation ratios.
#
# Manuscript      : Fig. 4a-c, Fig. S5, Fig. S6
# Main outputs    : 1_p_MGE_number.pdf, 1_mge_ds_density_comparison*.pdf, 1_MGEnum_*.pdf,
#                   1_DSnum_glacier_plateau_plain_comprison.pdf,
#                   01_MGE_cargo_and_defense_density.RData
# Original script : 20gaoyuan/12_5455MAG_MGEs/1_MGE_combination.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

#MGE数量统计、生境分布、携带的DS。

#每个MAG上携带的MGE数量、length统计####
setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
library(dplyr)

#phage/prophage
genomad_virus <- read.delim2("merged_5467MAG_virus_summary.tsv")

prophage <- genomad_virus[genomad_virus$topology == "Provirus",]
# 提取 | 前面的部分
prophage$contig <- sub("\\|.*$", "", prophage$seq_name)

phage <- genomad_virus[genomad_virus$topology != "Provirus",]

#plasmid
plasmid <- read.delim2("merged_5467MAG_plasmid_summary.tsv")

#ICEs/IMEs
ice <- read.delim2("processed_ice_sequences_ID.txt",header = F)
colnames(ice) <- c("seq_name","range","length")

#IS
IS <- read.delim2("5455MAG_filtered_completeIS.tsv")

#integron
integrons <- read.delim2("all5455MAG_contig.integrons")
complete_integrons <- integrons[integrons$type == "complete",]
# 处理数据
complete_integrons_clean <- complete_integrons %>%
  group_by(ID_replicon) %>%          # 按ID_replicon分组
  summarise(
    min_pos_beg = min(pos_beg),      # 最小起始位置
    max_pos_end = max(pos_end),      # 最大结束位置
    length = max_pos_end - min_pos_beg +1,
    .groups = "drop"                 # 取消分组
  ) %>%
  arrange(ID_replicon)               # 按ID排序


#查看这些MGE contig有多少重叠####
genomad_MGEid <- c(genomad_virus$seq_name , plasmid$seq_name, prophage$contig)
iceid <- ice$seq_name
ISid <- unique(IS$seqID)

inter1 <- intersect(genomad_MGEid,iceid) #24,23是质粒,1个prophage
inter2 <- intersect(genomad_MGEid,ISid) #617，563属于质粒，17属于phage，37属于prophage
inter3 <- intersect(iceid,ISid) #8  这8个里有7个的位置完全不重叠
inter4 <- intersect(genomad_MGEid, complete_integrons_clean$ID_replicon) #23 6个质粒、17phage
inter5 <- intersect(iceid, complete_integrons_clean$ID_replicon) #0
inter6 <- intersect(ISid,  complete_integrons_clean$ID_replicon) #9 其中4个完全不重叠

inter1_info <- plasmid[plasmid$seq_name %in% inter1,] #23
inter2_info_plasmid <- plasmid[plasmid$seq_name %in% inter2,] #563
inter2_info_virus <- phage[phage$seq_name %in% inter2,] #17
inter2_info_prophage <- prophage[prophage$contig %in% inter2,] #37
inter2_info_IS <- IS[IS$seqID %in% inter2_info_virus$seq_name,]
inter2_info_IS2 <- IS[IS$seqID %in% inter2_info_prophage$contig,]

inter3_info_plasmid <- plasmid[plasmid$seq_name %in% inter3,] #2
inter3_info_ice <- ice[ice$seq_name %in% inter3,]
inter3_info_IS <- IS[IS$seqID %in% inter3,]

inter4_info_plasmid <- plasmid[plasmid$seq_name %in% inter4,] #6
inter4_info_virus <- genomad_virus[genomad_virus$seq_name %in% inter4,] #17

inter6_info_IS <- IS[IS$seqID %in% inter6,]
inter6_info_integron <- complete_integrons_clean[complete_integrons_clean$ID_replicon %in% inter6,]

#基本上都是plasmid/phage与其他MGE存在重叠，即质粒上携带其他类型MGE。这种统一归到质粒/phage。
#删除其他MGE id中存在的plasmid contig id#####
ICE_filter <- ice[!ice$seq_name %in% inter1,]
IS_filter <- IS[!IS$seqID %in% inter2,]
integrons_filter <- complete_integrons_clean[!complete_integrons_clean$ID_replicon %in% inter4,]
head(prophage)
head(phage)
head(plasmid)


ICE_filter$genome <- sub("-[0-9]+$", "", ICE_filter$seq_name)
ICE_filter$MGE <- "ICEs"
ICE_filter <- select(ICE_filter, genome, MGE, length, seq_name)

IS_filter$genome <- sub("-[0-9]+$", "", IS_filter$seqID)
IS_filter$MGE <- "ISs"
IS_filter <- select(IS_filter, genome, MGE, isLen, seqID)
colnames(IS_filter) <- c("genome","MGE","length","seq_name")

integrons_filter$genome <- sub("-[0-9]+$", "", integrons_filter$ID_replicon)
integrons_filter$MGE <- "Integrons"
integrons_filter <- select(integrons_filter, genome, MGE, length, ID_replicon)
colnames(integrons_filter) <- c("genome","MGE","length","seq_name")

prophage$genome <- sub("-[0-9]+$", "", prophage$contig)
prophage$MGE <- "Prophage"
prophage <- select(prophage, genome, MGE, length, seq_name)

phage$genome <- sub("-[0-9]+$", "", phage$seq_name)
phage$MGE <- "Phage"
phage <- select(phage, genome, MGE, length, seq_name)

plasmid$genome <- sub("-[0-9]+$", "", plasmid$seq_name)
plasmid$MGE <- "Plasmid"
plasmid <- select(plasmid, genome, MGE, length, seq_name)


all_MGE <- rbind(ICE_filter, IS_filter, integrons_filter, prophage, phage, plasmid)

library(dplyr)
library(tidyr)

# 转换数据框
genome_MGE_length <- all_MGE %>%
  group_by(genome, MGE) %>%        # 按genome和MGE分组
  summarise(total_length = sum(length),  # 计算每组length总和
            .groups = "drop") %>%      # 取消分组
  pivot_wider(names_from = MGE,        # 将MGE作为列名
              values_from = total_length,  # 值来自total_length
              values_fill = 0)         # 缺失值填充为0
genome_MGE_length$MGE <- rowSums(genome_MGE_length[,-1])

#write.table(genome_MGE_length, "genome_MGE_length.tsv",quote = F,sep = '\t',row.names = F)

# 查看结果
head(genome_MGE_length)

genome_MGE_number <- all_MGE %>% 
  group_by(genome,MGE) %>%
  summarise(number = length(length),
            .groups = "drop") %>%
  pivot_wider(names_from = MGE,
              values_from = number,
              values_fill = 0)
genome_MGE_number$MGE <- rowSums(genome_MGE_number[,-1])

#write.table(genome_MGE_number, "genome_MGE_number.tsv",quote = F,sep = '\t',row.names = F)

#绘图表达不同MGE的数量

MGE_number <- data.frame(MGEtype = colnames(genome_MGE_number)[2:7],
                         MGEnum = colSums(genome_MGE_number[,c(2:7)]))

library(ggplot2)
p_MGE_number <- 
ggplot(MGE_number, aes(x = reorder(MGEtype, -MGEnum), y = MGEnum, fill = MGEtype)) +
  geom_bar(stat = "identity", width = 0.7, color = "white", alpha = 0.8) +  # 柱状图
  geom_text(aes(label = MGEnum), vjust = -0.5, size = 3.5) +               # 添加数值标签
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integrons"="#e8a598", "ICEs"="#00BA38",
                               "Chromosome" = "#FCDABA",
                               "ISs" = "#9f7e69",
                               "Phage" = "#ED6764")) +                                    # 使用Set2配色
  labs(
    title = "",
    x = "MGE Type",
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
    axis.text.x = element_text(angle = 45, size = 10, hjust = 1, vjust = 1,color = 'black'),  # 设置 x 轴文本
    axis.text.y = element_text(size = 10, color = 'black'),  # 设置 y 轴文本
    legend.position = "none"  # 去除图例
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))                 # 调整Y轴范围
p_MGE_number
setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
ggsave("1_p_MGE_number.pdf",p_MGE_number,width = 3,height = 4,units = "in")




#MGE与DS关系##########
DS <- read.delim2("merged_5467MAG_defense_finder_systems.tsv")
DG <- read.delim2("merged_5467MAG_defense_finder_genes.tsv")

#每个MGE都能容易得到它们（prophage，IS，ICEs，Integrons）在contig上的位置信息。

#需要根据gff文件得到每个ORF在contig上的位置，为DS和DG赋予contig位置信息。

#ORFgff <- read.delim2("DG_gff.gff", header = F, sep = '\t')
ORFgff <- read.delim2("DG_gff_DF.gff", header = F, sep = '\t')
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


##统计MAG上所有的DS####
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


MAG_metadata <- read.delim2("MAG5455_metadata.tsv")
count_matrix <- count_matrix[count_matrix$genome %in% MAG_metadata$genome,]
DF_defense <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]

count_matrix <- count_matrix[-1]  # 移除genome列

##统计质粒上的DS####
head(plasmid)
head(DF_defense)

DF_defense <- DF_defense %>%
  mutate(
    seq_name = str_extract(sys_id, ".*bin\\.\\d+\\-\\d+")  # 包含"bin.数字"
  )

plasmid_DS <- DF_defense[DF_defense$seq_name %in% plasmid$seq_name,]

# 创建计数矩阵
plasmid_DS_statistic <- plasmid_DS %>%
  count(genome, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

##统计phage上的DS####
phage_DS <- DF_defense[DF_defense$seq_name %in% phage$seq_name, ]

# 创建计数矩阵
phage_DS_statistic <- phage_DS %>%
  count(genome, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

##统计prophage上的DS####

prophage <- prophage %>%
  mutate(
    contig = str_extract(seq_name, ".*bin\\.\\d+\\-\\d+"),  # 提取包含"bin.数字-数字"的部分
    position = str_extract(seq_name, "\\d+\\_\\d+"),       # 提取"数字_数字"格式的位置信息
    begin = as.numeric(str_extract(position, "^\\d+")),    # 提取_前的数字作为起始位置
    end = as.numeric(str_extract(position, "\\d+$"))       # 提取_后的数字作为结束位置
  )

#根据DG在contig上的位置信息，计算每个DS在contig上的位置信息

DG <- DG %>% 
  mutate(genome = str_extract(replicon, ".*bin\\.\\d+"))
DG <- DG[DG$activity == "Defense",]
DG <- DG[DG$genome %in% MAG_metadata$genome,]

#用DF_defense的sys_beg去匹配DG数据框的hit_id列，匹配到的行的V4列的值作为DF_defense对应行的begin
#用DF_defense的sys_end去匹配DG数据框的hit_id列，匹配到的行的V5列的值作为DF_defense对应行的end

# 首先为DF_defense创建begin和end列（如果不存在）
if(!"begin" %in% names(DF_defense)) DF_defense$begin <- NA
if(!"end" %in% names(DF_defense)) DF_defense$end <- NA

# 创建从hit_id到V4/V5的映射
v4_mapping <- setNames(DG$V4, DG$hit_id)
v5_mapping <- setNames(DG$V5, DG$hit_id)

# 用sys_beg匹配并填充begin列
DF_defense$begin <- ifelse(
  DF_defense$sys_beg %in% names(v4_mapping),
  v4_mapping[as.character(DF_defense$sys_beg)],
  DF_defense$begin
)

# 用sys_end匹配并填充end列
DF_defense$end <- ifelse(
  DF_defense$sys_end %in% names(v5_mapping),
  v5_mapping[as.character(DF_defense$sys_end)],
  DF_defense$end
)

# 验证结果
head(DF_defense[, c("sys_beg", "sys_end", "begin", "end")])

Prophage_DS <- DF_defense[DF_defense$seq_name %in% prophage$contig,]
prophage_2 <- select(prophage, contig, begin, end)
colnames(prophage_2) <- c("seq_name","begin_p","end_p")
Prophage_DS <- left_join(Prophage_DS, prophage_2, by = "seq_name")
#对Prophage_DS进行过滤，只保留begin到end， 与begin_p到end_p，二者的范围存在交叉的行。
Prophage_DS_filtered <- Prophage_DS[
  (Prophage_DS$begin <= Prophage_DS$end_p & Prophage_DS$end >= Prophage_DS$begin_p) |
    (Prophage_DS$begin_p <= Prophage_DS$end & Prophage_DS$end_p >= Prophage_DS$begin),]
Prophage_DS_filtered <- Prophage_DS_filtered[!is.na(Prophage_DS_filtered$sys_id),]

# 创建计数矩阵
prophage_DS_statistic <- Prophage_DS_filtered %>%
  count(genome, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

##Integrons(complete)上的DS####
integron_middle <- complete_integrons_clean
colnames(integron_middle)[1] <- "seq_name"
integrons_filter <- left_join(integrons_filter, integron_middle, by = "seq_name")

integron_DS <- DF_defense[DF_defense$seq_name %in% integrons_filter$seq_name,]

integron_DS <- left_join(integron_DS, integrons_filter, by = "seq_name")

#对integron_DS进行过滤，只保留begin到end， 与begin_p到end_p，二者的范围存在交叉的行。
integron_DS_filtered <- integron_DS[
  (integron_DS$begin <= integron_DS$max_pos_end & integron_DS$end >= integron_DS$min_pos_beg) |
    (integron_DS$min_pos_beg <= integron_DS$end & integron_DS$max_pos_end >= integron_DS$begin),]

# 创建计数矩阵
integron_DS_statistic <- integron_DS_filtered %>%
  count(genome.x, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome.x, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

##ICE上的DS####
ICE_filter_middle <- ICE_filter
ICE_filter_middle <- left_join(ICE_filter_middle, ice, by = "seq_name")
ICE_filter_middle <- ICE_filter_middle %>%
  mutate(ICEbegin = as.numeric(str_extract(range, "^\\d+")),    # 提取_前的数字作为起始位置
         ICEend = as.numeric(str_extract(range, "\\d+$"))  )

ICE_DS <- DF_defense[DF_defense$seq_name %in% ICE_filter$seq_name,]
ICE_DS <- left_join(ICE_DS, ICE_filter_middle, by = "seq_name")
#对ICE_DS进行过滤，只保留begin到end， 与ICEbegin到ICEend，二者的范围存在交叉的行。
ICE_DS_filtered <- ICE_DS[
  (ICE_DS$begin <= ICE_DS$ICEend & ICE_DS$end >= ICE_DS$ICEbeg) |
    (ICE_DS$ICEbeg <= ICE_DS$end & ICE_DS$ICEend >= ICE_DS$begin),]
ICE_DS_filtered <- ICE_DS_filtered[!is.na(ICE_DS_filtered$sys_id),]

# 创建计数矩阵
ICE_DS_statistic <- ICE_DS_filtered %>%
  count(genome.x, type) %>%  # 计算每个genome-type组合的出现次数
  complete(genome.x, type, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()

##IS上的DS####
IS_middle <- IS
colnames(IS_middle)[1] <- "seq_name"
colnames(IS_middle)[6] <- "length"
IS_filter_middle <- left_join(IS_filter, IS_middle, by = c("seq_name","length") )

IS_DS <- DF_defense[DF_defense$seq_name %in% IS_filter_middle$seq_name,]
IS_DS <- left_join(IS_DS, IS_filter_middle, by = c("seq_name"))
#对IS_DS进行过滤，只保留begin到end， 与isBegin到isEnd，二者的范围存在交叉的行。
IS_DS_filtered <- IS_DS[
  (IS_DS$begin <= IS_DS$isEnd & IS_DS$end >= IS_DS$isBegin) |
    (IS_DS$isBegin <= IS_DS$end & IS_DS$isEnd >= IS_DS$begin),]
IS_DS_filtered <- IS_DS_filtered[!is.na(IS_DS_filtered$sys_id),]

# 创建计数矩阵
IS_DS_statistic <- IS_DS_filtered %>%
  count(genome.x, type.x) %>%  # 计算每个genome-type组合的出现次数
  complete(genome.x, type.x, fill = list(n = 0)) %>%  # 填充缺失组合为0
  pivot_wider(
    names_from = type.x,
    values_from = n,
    values_fill = 0
  ) %>%
  as.data.frame()



##合并分析####
MAG_metadata <- read.delim2("MAG5455_metadata.tsv")
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

count_matrix <- count_matrix[rownames(count_matrix) %in% MAG_metadata$genome,]

MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]
colnames(MAG_quality)[1] <- "genome"

MAG_DSsum <- rowSums(count_matrix) %>% as.data.frame() %>% setNames("DSsum")
MAG_DSsum$genome <- rownames(MAG_DSsum)

MAG_quality <- left_join(MAG_quality, MAG_DSsum, by = "genome")

MAG_quality$MAG_size <- MAG_quality$Genome_Size/1000000
#MAG_quality$DS_density <- MAG_quality$DSsum/MAG_quality$Genome_Size*1000000
MAG_quality$magsize <- floor(MAG_quality$MAG_size)
#MAG_quality[is.na(MAG_quality)] <- 0


plasmid_DS_num <- data.frame(genome = plasmid_DS_statistic$genome,
                             plasmid_DS_num = rowSums(plasmid_DS_statistic[,-1])
                             )
phage_DS_num <- data.frame(genome = phage_DS_statistic$genome,
                             phage_DS_num = rowSums(phage_DS_statistic[,-1])
)
prophage_DS_num <- data.frame(genome = prophage_DS_statistic$genome,
                              prophage_DS_num = rowSums(prophage_DS_statistic[,-1])
)
integron_DS_num <- data.frame(genome = integron_DS_statistic$genome,
                              integron_DS_num = rowSums(integron_DS_statistic[,-1])
)
ICE_DS_num <- data.frame(genome = ICE_DS_statistic$genome,
                         ICE_DS_num = rowSums(ICE_DS_statistic[,-1])
)
IS_DS_num <- data.frame(genome = IS_DS_statistic$genome,
                        IS_DS_num = rowSums(IS_DS_statistic[,-1])
)

MGE_length <- read.delim2("genome_MGE_length.tsv")

MAG_quality <- MAG_quality %>%
  full_join(plasmid_DS_num,by = "genome") %>%
  full_join(phage_DS_num,by = "genome") %>%
  full_join(prophage_DS_num,by = "genome") %>%
  full_join(integron_DS_num,by = "genome") %>%
  full_join(ICE_DS_num,by = "genome") %>%
  full_join(IS_DS_num,by = "genome") 

MAG_quality <- left_join(MAG_quality, MGE_length, by = "genome")
MAG_quality[is.na(MAG_quality)] <- 0

MAG_quality$Chromosome_DS_num <- MAG_quality$DSsum - 
  (MAG_quality$plasmid_DS_num+MAG_quality$phage_DS_num+
     MAG_quality$prophage_DS_num+MAG_quality$integron_DS_num+
     MAG_quality$ICE_DS_num+MAG_quality$IS_DS_num)

MAG_quality$Chromosome_length <- MAG_quality$Genome_Size - MAG_quality$MGE


MAG_quality$Chromosome_DS_density <- MAG_quality$Chromosome_DS_num/MAG_quality$Chromosome_length*1000
MAG_quality$plasmid_DS_density <- MAG_quality$plasmid_DS_num/MAG_quality$Plasmid*1000
MAG_quality$phage_DS_density <- MAG_quality$phage_DS_num/MAG_quality$Phage*1000
MAG_quality$prophage_DS_density <- MAG_quality$prophage_DS_num/MAG_quality$Prophage*1000
MAG_quality$Integrons_DS_density <- MAG_quality$integron_DS_num/MAG_quality$Integrons*1000
MAG_quality$ICE_DS_density <- MAG_quality$ICE_DS_num/MAG_quality$ICEs*1000
MAG_quality$IS_DS_density <- MAG_quality$IS_DS_num/MAG_quality$ISs*1000

#prophage
prophage_ds_density <- select(MAG_quality,Prophage,prophage_DS_num,prophage_DS_density)
prophage_ds_density <- prophage_ds_density[prophage_ds_density$Prophage>0,]
prophage_ds_density$location <- "Prophage"
colnames(prophage_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

#phage上
phage_ds_density <- select(MAG_quality,Phage,phage_DS_num,phage_DS_density)
phage_ds_density <- phage_ds_density[phage_ds_density$Phage>0,]
phage_ds_density$location <- "Phage"
colnames(phage_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

#质粒上
plasmid_ds_density <- select(MAG_quality,Plasmid,plasmid_DS_num,plasmid_DS_density)
plasmid_ds_density <- plasmid_ds_density[plasmid_ds_density$Plasmid>0,]
plasmid_ds_density$location <- "Plasmid"
colnames(plasmid_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

#integron
integron_ds_density <- select(MAG_quality,Integrons,integron_DS_num,Integrons_DS_density)
integron_ds_density <- integron_ds_density[integron_ds_density$Integrons>0,]
integron_ds_density$location <- "Integron"
colnames(integron_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

#ICE
ICE_ds_density <- select(MAG_quality,ICEs,ICE_DS_num,ICE_DS_density)
ICE_ds_density <- ICE_ds_density[ICE_ds_density$ICEs>0,]
ICE_ds_density$location <- "ICE/IME"
colnames(ICE_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

#IS
IS_ds_density <- select(MAG_quality,ISs,IS_DS_num,IS_DS_density)
IS_ds_density <- IS_ds_density[IS_ds_density$ISs>0,]
IS_ds_density$location <- "IS"
colnames(IS_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

#染色体上
chromosome_ds_density <- select(MAG_quality,Chromosome_length,Chromosome_DS_num,Chromosome_DS_density)
chromosome_ds_density <- chromosome_ds_density[chromosome_ds_density$Chromosome_length>0,]
chromosome_ds_density$location <- "Chromosome"
colnames(chromosome_ds_density)[1:3] <- c("size","dsnumber","dsdensity")

mge_ds_density <- rbind(prophage_ds_density,phage_ds_density,
                        plasmid_ds_density,integron_ds_density,
                        ICE_ds_density,IS_ds_density,
                        chromosome_ds_density)
mge_ds_density$location <- factor(mge_ds_density$location, levels = c("Chromosome", "Plasmid","Phage",
                                                                      "Prophage","Integron","ICE/IME","IS"))
##DS density comparison绘图####
library(rstatix)
md_p_val <- mge_ds_density[mge_ds_density$dsnumber>0,] %>% 
  wilcox_test(dsdensity ~ location) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="location",dodge=0.5)


library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means

#md_p_val$y.position <- c(0.4,0.44,0.48,0.52,0.56,0.60,0.64,0.68,0.72,0.76,0.8,0.84,0.88,0.92,0.96)
md_p_val$y.position <- seq(from = 0.3, to = 0.6, length.out = 21)
p_mge_ds_density <- 
  ggplot(mge_ds_density[mge_ds_density$dsnumber>0,],aes(x=location,y=dsdensity))+
  geom_jitter(aes(group=1,color = location),width = 0.2, alpha = 0.3)+  # 显示异常值点
  geom_boxplot(aes(fill = location),width = 0.2,outlier.shape = NA, staplewidth = 0.5,linewidth=0.2) +  # 隐藏箱线图中的异常值点
  geom_violin(trim = T, aes(fill = location),show.legend = F,linewidth=0.2) + # 绘制小提琴图
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  stat_pvalue_manual(md_p_val,label = "p.adj.signif",hide.ns=T,
                     tip.length = 0,label.size = 5,color="black")+
  stat_summary(fun=mean,geom="point",col="black",fill="#F98400",
               shape=23,size=1.2,show.legend = F)+
  #stat_summary(fun=mean, geom="line", aes(group=CPR), col="black")+
  scale_fill_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                               "Integron"="#e8a598", "ICE/IME"="#00BA38",
                               "Chromosome" = "#FCDABA",
                               "IS" = "#9f7e69",
                               "Phage" = "#ED6764")) +
  scale_color_manual(values = c("Plasmid" = "#619CFF", "Prophage"="#fa9500",
                                "Integron"="#e8a598", "ICE/IME"="#00BA38",
                                "Chromosome" = "#FCDABA","IS" = "#9f7e69",
                                "Phage" = "#ED6764")) +
  #ylim(0,15)+
  #scale_y_continuous(limits = c(0,10))+
  coord_cartesian(ylim = c(0,0.6))+
  labs(x=NULL,y='DS density (per kb)')+
  mytheme
p_mge_ds_density

setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
ggsave("1_mge_ds_density_comparison.pdf",p_mge_ds_density,width = 5.5,height = 3.5,units = "in")
write.table(mge_ds_density, "mge_ds_density.tsv", sep = "\t", quote = F, row.names = F)
ggsave("1_mge_ds_density_comparison_new.pdf",
p_mge_ds_density,width = 7.87,height = 3,units = "in")


mytheme <-  theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( size = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(size = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
  )


##O/E分析####
ds_family_plasmid <- colSums(plasmid_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("Plasmid")
ds_family_plasmid$dsfamily <- rownames(ds_family_plasmid)

ds_family_phage <- colSums(phage_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("Phage")
ds_family_phage$dsfamily <- rownames(ds_family_phage)

ds_family_prophage <- colSums(prophage_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("Prophage")
ds_family_prophage$dsfamily <- rownames(ds_family_prophage)

ds_family_integron <- colSums(integron_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("Integron")
ds_family_integron$dsfamily <- rownames(ds_family_integron)

ds_family_ICE <- colSums(ICE_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("ICE")
ds_family_ICE$dsfamily <- rownames(ds_family_ICE)

ds_family_IS <- colSums(IS_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("IS")
ds_family_IS$dsfamily <- rownames(ds_family_IS)


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
count_matrix <- count_matrix[count_matrix$genome %in% MAG_metadata$genome,]

ds_family_chromosome <- colSums(count_matrix[,-1]) %>% as.data.frame() %>% setNames("All")
ds_family_chromosome$dsfamily <- rownames(ds_family_chromosome)

ds_family_bind <- ds_family_chromosome %>%
  full_join(ds_family_phage,by = 'dsfamily') %>%
  full_join(ds_family_plasmid,by = 'dsfamily') %>%
  full_join(ds_family_prophage,by = 'dsfamily') %>%
  full_join(ds_family_integron,by = 'dsfamily') %>%
  full_join(ds_family_ICE,by = 'dsfamily') %>%
  full_join(ds_family_IS,by = 'dsfamily')

ds_family_bind[is.na(ds_family_bind)] <- 0

ds_family_bind$Chromosome <- ds_family_bind$All - 
  ds_family_bind$Plasmid - 
  ds_family_bind$Prophage -
  ds_family_bind$Integron-
  ds_family_bind$ICE -
  ds_family_bind$IS -
  ds_family_bind$Phage

ds_family_final <- select(ds_family_bind,dsfamily,Chromosome,Plasmid,Phage,
                          Prophage,Integron,ICE,IS) %>% as.data.frame()

rownames(ds_family_final) <- ds_family_final$dsfamily
ds_family_final <- ds_family_final[,-1]
sum(rowSums(ds_family_final)>20)
ds_family_final <- ds_family_final[rowSums(ds_family_final)>20,]


oe <- chisq.test(ds_family_final,simulate.p.value = T)

oe_bi <- oe$observed/oe$expected %>% as.data.frame()
oe_bi <- na.omit(oe_bi)
#oe_bi <- oe_bi[oe_bi$Chromosome!=0,]
oe_bi <- oe_bi[!rowSums(oe_bi[,-1]) == 0,]

#绘制热图
library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数

# 定义颜色映射函数
col_fun <- colorRamp2(c(0,1, 2), c("white","#B6B8D8", "#FF7F0E"))  # 粉色到绿色的渐变
#setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
cairo_pdf(width = 8.35,height = 3,filename = "1_mag_69DS_oe.pdf")
p_mge_ds_oe <- Heatmap(t(oe_bi),
                       col = col_fun,
                       cluster_rows = F,
                       cluster_columns = T,
                       name = 'O/E',
                       row_title = NULL,
                       row_split = colnames(oe_bi),
                       row_names_side = "left", # 行标题显示在左侧
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 8),
                       column_names_rot = 90,
                       rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
                       heatmap_legend_param = list(
                         title = "O/E", at = c(0,1, 2),#direction = "horizontal",
                         labels = c("0","1", "≥2"))
)
p_mge_ds_oe
dev.off()
#setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
write.table(oe_bi,"oe_bi.tsv",sep = "\t", quote = F)

#筛选MGE的contig范围与DS有交叉的####
##contig筛选####
head(plasmid_DS)
head(phage_DS)

head(integron_DS_filtered) #共6个，DS均被integron包含
head(ICE_DS_filtered)
head(IS_DS_filtered)
head(Prophage_DS_filtered)

ICE_DS_filtered_mix <- ICE_DS_filtered[
  (ICE_DS_filtered$begin <= ICE_DS_filtered$ICEend & ICE_DS_filtered$ICEend <= ICE_DS_filtered$end) |
    (ICE_DS_filtered$ICEbeg <= ICE_DS_filtered$end & ICE_DS_filtered$ICEbegin >= ICE_DS_filtered$begin),]
#0个

integron_DS_filtered_mix <- integron_DS_filtered[
  (integron_DS_filtered$begin <= integron_DS_filtered$max_pos_end & integron_DS_filtered$max_pos_end <= integron_DS_filtered$end) |
    (integron_DS_filtered$min_pos_beg <= integron_DS_filtered$end & integron_DS_filtered$min_pos_beg >= integron_DS_filtered$begin),]
#0

IS_DS_filtered_mix <- IS_DS_filtered[
  (IS_DS_filtered$begin <= IS_DS_filtered$isEnd & IS_DS_filtered$isEnd <= IS_DS_filtered$end) |
    (IS_DS_filtered$isBegin <= IS_DS_filtered$end & IS_DS_filtered$isBegin >= IS_DS_filtered$begin),]
#33个

Prophage_DS_filtered_mix <- Prophage_DS_filtered[
  (Prophage_DS_filtered$begin <= Prophage_DS_filtered$end_p & Prophage_DS_filtered$end_p <= Prophage_DS_filtered$end) |
    (Prophage_DS_filtered$begin_p <= Prophage_DS_filtered$end & Prophage_DS_filtered$begin_p >= Prophage_DS_filtered$begin),]
#11个


#查看筛选到的MGE的详细质量信息
#prophage
selected_prophage <- prophage[prophage$contig %in% Prophage_DS_filtered_mix$seq_name,]

#ISs
selected_IS <- IS_filter_middle[IS_filter_middle$seq_name %in% IS_DS_filtered_mix$seq_name,]

#ICE
selected_ICE <- ICE_DS_filtered_mix

##绘制selected contig的genetic background。####
save.image()

setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
# loading the workspace
load(".RData")
###ICE contig####
selected_ICE$seq_name

setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
library(readr)
write_tsv(IS_DS_filtered_mix, "IS_DS_filtered_mix.tsv")
write_tsv(selected_prophage,"selected_prophage.tsv")
##MGE插入导致的移码突变####


#统计比较不同生境注释到的MGE的数量####
setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
MAG_MGE_num <- read.delim2("genome_MGE_number.tsv")

MAG_metadata <- read.delim2("MAG5455_metadata.tsv")

MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'Yes',"Plateau",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'No',"Plain",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'glacier',"Glacier",MAG_metadata$Plateau)
colnames(MAG_metadata)[16] <- 'Region'

MAG_metadata$medium <- MAG_metadata$Type
# MAG_metadata$medium <- ifelse(MAG_metadata$Type == 'Cryoconite',"Glacier",MAG_metadata$Type)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Ice',"Glacier",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'Snow',"Glacier",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'particle',"Plankton",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'water',"Plankton",MAG_metadata$medium)
# MAG_metadata$medium <- ifelse(MAG_metadata$medium == 'sediment',"Benthon",MAG_metadata$medium)

MAG_MGE_num <- MAG_MGE_num[MAG_MGE_num$genome %in% MAG_metadata$genome,]

MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]
colnames(MAG_quality)[1] <- "genome"



MAG_metadata <- merge(MAG_metadata, MAG_quality, 
                      by.x = "genome", by.y = "genome")
MAG_metadata <- merge(MAG_metadata, MAG_MGE_num, 
                      by.x = "genome", by.y = "genome", all = T)

MAG_metadata[,75:81][is.na(MAG_metadata[,75:81])] <- 0

MAG_metadata$ecosystem <- factor(MAG_metadata$ecosystem,
                              levels = c("Glacier","HH","CJ","LCJ","YJ"))
MAG_metadata$medium <- factor(MAG_metadata$medium,
                              levels = c("Snow","Ice","Cryoconite","water","particle","sediment"))
MAG_metadata$Region <- factor(MAG_metadata$Region,
                              levels = c("Glacier","Plateau","Plain"))

##绘制箱线+小提琴+散点复合图####

MAG_metadata_plateau_plain <- select(MAG_metadata,genome,ecosystem,medium,Region,MGE)
MAG_metadata_plateau_plain$Region <- as.character(MAG_metadata_plateau_plain$Region) 

MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
                                           MAG_metadata_plateau_plain$Region)

md_p_val0 <- MAG_metadata_plateau_plain %>% 
#  group_by(ecosystem) %>%
  wilcox_test(MGE ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val0$y.position <- md_p_val0$y.position - 50

library(gghalves)
library(ggsci)
p_plateau_plain <- 
  ggplot(MAG_metadata_plateau_plain, 
         aes(x = Region, y = MGE)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  # geom_half_point(aes(color=Region),
  #                 side = "r",size = 0.5,
  #                 transformation_params = list(height = 0,width = 0.001,seed = 2))+
  geom_jitter(aes(color=Region), width = 0.1, size = 0.5, alpha = 0.6)+
  stat_pvalue_manual(md_p_val0,label = "p.adj.signif",label.size=5,hide.ns = T)+
 # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  scale_fill_npg()+
  scale_color_npg()+
  labs(x="",y="MGE number")+
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
  coord_cartesian()
p_plateau_plain

# setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
ggsave("1_MGEnum_plateau_plain_comprison.pdf",p_plateau_plain, 
       width = 2.6, height = 2.78, units = "in")

library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)

#计算显著性差异
library(rstatix)
md_p_val1 <- MAG_metadata %>% 
  group_by(ecosystem) %>%
  wilcox_test(MGE ~ medium) %>%
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

md_p_val1$y.position <- md_p_val1$y.position-50

# 定义 6 种 medium 的固定颜色
medium_colors <- c("Snow" = "#4DBBD5",
                   "Ice" = "#E64B35",
                   "Cryoconite" = "#00A087",
                   "water" = "#3C5488",
                   "particle" = "#F39B7F",
                   "sediment" = "#8491B4")

# 各 ecosystem 实际存在的 medium 组合
ecosystem_media <- MAG_metadata |>
  group_by(ecosystem, medium) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 ecosystem 独立绘图
plot_list_eco_med <- list()
ecosystems <- levels(MAG_metadata$ecosystem)
n_media <- c()

for (i in seq_along(ecosystems)) {
  eco <- ecosystems[i]
  
  existing_media <- ecosystem_media |>
    filter(ecosystem == eco) |>
    pull(medium)
  n_media <- c(n_media, length(existing_media))
  
  df_sub <- MAG_metadata |>
    filter(ecosystem == eco) |>
    mutate(medium = factor(medium, levels = existing_media))
  
  # 重新计算显著性
  if (length(existing_media) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(MGE ~ medium) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "medium", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 50
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = medium, y = MGE)) +
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
  
  plot_list_eco_med[[eco]] <- p
}

# 合并：各 panel 宽度与 media 数量成正比
p_ecosystem_medium <- plot_grid(plotlist = plot_list_eco_med, nrow = 1,
                                 rel_widths = n_media)
p_ecosystem_medium



# 定义 ecosystem 固定颜色
ecosystem_colors <- c("YJ" = "#FCB2AF",
                      "Glacier" = "#9BDFDF",
                      "LCJ" = "#FFE2CE",
                      "CJ" = "#BEBCDF",
                      "HH" = "#FB8C62")

# 各 Region 实际存在的 ecosystem 组合（排除 Glacier）
region_ecosystem_data <- MAG_metadata |>
  filter(ecosystem != "Glacier") |>
  group_by(Region, ecosystem) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 Region 独立绘图
plot_list_region_eco <- list()
regions_for_eco <- region_ecosystem_data |> pull(Region) |> as.character() |> unique()
n_eco_per_region <- c()

for (i in seq_along(regions_for_eco)) {
  reg_val <- regions_for_eco[i]
  
  existing_eco <- region_ecosystem_data |>
    filter(Region == reg_val) |>
    pull(ecosystem)
  n_eco_per_region <- c(n_eco_per_region, length(existing_eco))
  
  df_sub <- MAG_metadata |>
    filter(ecosystem != "Glacier", Region == reg_val) |>
    mutate(ecosystem = factor(ecosystem, levels = existing_eco))
  
  # 重新计算显著性
  if (length(existing_eco) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(MGE ~ ecosystem) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "ecosystem", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 50
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = ecosystem, y = MGE)) +
    geom_half_boxplot(aes(fill = ecosystem), color = "black",
                      side = "l", errorbar.draw = T,
                      outlier.shape = NA, width = 0.8) +
    geom_jitter(aes(color = ecosystem), width = 0.1, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = ecosystem_colors) +
    scale_color_manual(values = ecosystem_colors) +
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
  
  plot_list_region_eco[[reg_val]] <- p
}

# 合并：各 panel 宽度与 ecosystem 数量成正比
p_Region_ecosystem <- plot_grid(plotlist = plot_list_region_eco, nrow = 1,
                                 rel_widths = n_eco_per_region)
p_Region_ecosystem

write.table(df_sub, "df_sub_MGE_num.tsv",sep = "\t",quote = F, row.names = F)

##海拔与DS的关系##################

# 创建一个新的分组变量：每500单位一个Elevation组
MAG_metadata$Elevation <- as.numeric(MAG_metadata$Elevation)
MAG_metadata_elevation <- MAG_metadata %>%
  mutate(Elevation_Group = cut(Elevation, 
                               breaks = seq(floor(min(Elevation, na.rm = TRUE)/500)*500, 
                                            ceiling(max(Elevation, na.rm = TRUE)/500)*500, 
                                            by = 500),
                               include.lowest = TRUE,
                               right = FALSE))

# 计算每组的平均值和标准差
summary_stats <- MAG_metadata_elevation %>%
  group_by(Elevation_Group) %>%
  summarise(mean_MGE = mean(MGE, na.rm = TRUE),
            sd_MGE = sd(MGE, na.rm = TRUE),
            .groups = 'drop')


# 进行显著性检验：Elevation与DS_densitylog的相关性
cor_result <- cor.test(MAG_metadata_elevation$Elevation, 
                       MAG_metadata_elevation$MGE, 
                       alternative = "two.sided", 
                       method = "spearman", 
                       conf.level = 0.95)

# 计算每组的平均值和标准差
# 计算每组的平均值和标准差，并提取Elevation中位值（每组中Elevation的真实中位数）
summary_stats <- MAG_metadata_elevation %>%
  group_by(Elevation_Group) %>%
  summarise(mean_MGE = mean(MGE, na.rm = TRUE),
            sd_MGE = sd(MGE, na.rm = TRUE),
            Elevation_mid = median(Elevation, na.rm = TRUE),
            .groups = 'drop')

# 绘图并添加线性拟合线
p_elevation_MGE <- 
  ggplot(summary_stats, aes(x = Elevation_mid, y = mean_MGE)) +
  geom_point(size = 3) +
  geom_line() +
  geom_errorbar(aes(ymin = mean_MGE - sd_MGE,
                    ymax = mean_MGE + sd_MGE), 
                width = 200) +
  geom_smooth(method = "lm", se = TRUE, color = "blue", linetype = "dashed") +
  labs(x = "Elevation (m)",
       y = "MGE number (per MAG)",
       title = "MGE Number vs. Elevation Groups") +
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

print(p_elevation_MGE)

write.table(summary_stats, "summary_stats_mge.tsv",sep = "\t",quote = F,row.names = F)

# 进行显著性检验：Elevation与DS_densitylog的相关性（整体数据）
cor_result <- cor.test(MAG_metadata$Elevation, 
                       MAG_metadata$MGE, 
                       alternative = "two.sided", 
                       method = "spearman", 
                       conf.level = 0.95)

print(cor_result)

# 使用分组中位数和均值进行线性拟合相关性检验
lm_result <- lm(mean_MGE ~ Elevation_mid, data = summary_stats)
summary(lm_result)



library(Cairo)
library(cowplot)
cairo_pdf(filename = "4_geografical_DSdensity1.pdf" , width = 6.6,height = 6.9)
g1 <- plot_grid(p_Region_ecosystem,
                p_elevation_MGE,
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


##河流内部 高原与平原段比较########
library(dplyr)

MAG_metadata_extract <- select(MAG_metadata,genome,ecosystem,medium,Region,MGE)

MAG_metadata_extract$ecosystem <- factor(MAG_metadata_extract$ecosystem,
                                                levels = c("Glacier","HH","CJ","LCJ","YJ"))
MAG_metadata_extract$medium <- factor(MAG_metadata_extract$medium,
                                             levels = c("Snow","Ice","Cryoconite","water","particle","sediment"))
MAG_metadata_extract$Region <- as.character(MAG_metadata_extract$Region) 

MAG_metadata_extract$Region <- ifelse(MAG_metadata_extract$Region == "Glacier", "Plateau", 
                                      MAG_metadata_extract$Region)

MAG_metadata_extract$Region <- factor(MAG_metadata_extract$Region,
                                             levels = c("Plateau","Plain"))

# 定义 Region 固定颜色
region_colors <- c("Plateau" = "#C4D8E9", "Plain" = "#d4e09b")

# 各 ecosystem 实际存在的 Region 组合
eco_region_data <- MAG_metadata_extract |>
  group_by(ecosystem, Region) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 ecosystem 独立绘图
plot_list_eco_region <- list()
ecos_for_region <- eco_region_data |> pull(ecosystem) |> as.character() |> unique()
n_region_per_eco <- c()

for (i in seq_along(ecos_for_region)) {
  eco_val <- ecos_for_region[i]
  
  existing_regions <- eco_region_data |>
    filter(ecosystem == eco_val) |>
    pull(Region)
  n_region_per_eco <- c(n_region_per_eco, length(existing_regions))
  
  df_sub <- MAG_metadata_extract |>
    filter(ecosystem == eco_val) |>
    mutate(Region = factor(Region, levels = existing_regions))
  
  # 重新计算显著性
  if (length(existing_regions) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(MGE ~ Region) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "Region", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 50
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = Region, y = MGE)) +
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


# setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
dev.off()

##以Region分组，比较不同medium########
# 各 Region 实际存在的 medium 组合
region_media <- MAG_metadata_extract |>
  group_by(Region, medium) |>
  summarise(n = n(), .groups = "drop") |>
  filter(n > 0)

# 为每个 Region 独立绘图
plot_list_region_med <- list()
regions_med <- region_media |> pull(Region) |> unique()
n_med_per_region <- c()
# ds_density$log10density <- log10(ds_density$density+1)
# md_p_val2 <- ds_density %>% 
#   wilcox_test(log10density ~ region) %>%
#   adjust_pvalue(p.col="p",method="bonferroni") %>%
#   add_significance(p.col="p.adj") %>% 
#   add_xy_position(x="region",dodge=0.5)

#md_p_val1$y.position <- seq(from = 0.6, to = 0.6 + (0.03 * (21 - 1)), by = 0.03)
# 调整y轴位置
#md_p_val1$y.position <- max(all5455MAG_metadata_extract$DS_densitylog, na.rm = TRUE) * 
#(1 + 0.05 * (1:nrow(md_p_val1)))
# md_p_val4$xmin <- md_p_val4$xmin -1
# md_p_val4$xmax <- md_p_val4$xmax -1
for (i in seq_along(regions_med)) {
  reg_val <- regions_med[i]
  
  existing_media <- region_media |>
    filter(Region == reg_val) |>
    pull(medium)
  n_med_per_region <- c(n_med_per_region, length(existing_media))
  
  df_sub <- MAG_metadata_extract |>
    filter(Region == reg_val) |>
    mutate(medium = factor(medium, levels = existing_media))
  
  # 重新计算显著性
  if (length(existing_media) >= 2) {
    pval_sub <- df_sub |>
      wilcox_test(MGE ~ medium) |>
      adjust_pvalue(p.col = "p", method = "bonferroni") |>
      add_significance(p.col = "p.adj") |>
      add_xy_position(x = "medium", dodge = 0.5)
    pval_sub$y.position <- pval_sub$y.position - 100
  } else {
    pval_sub <- NULL
  }
  
  p <- ggplot(df_sub, aes(x = medium, y = MGE)) +
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
cairo_pdf(filename = "1_geografical_MGE5_new.pdf" , width = 7.5,height = 8)
g1 <- plot_grid(p_ecosystem_Region,
                p_elevation_MGE,
                labels = letters[4:5], ncol = 2,
                rel_widths = c(3, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

g2 <- plot_grid(p_Region_ecosystem,p_Region_medium,
                labels = letters[2:3], ncol = 2
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


# setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
dev.off()


# p_Region_ecosystem <- ggplot(MAG_metadata[MAG_metadata$ecosystem != "Glacier",], 
#                              aes(x = ecosystem, y = MGE)) +
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
#   labs(x="River",y="MGE number")+
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
# p_Region_ecosystem

# p_Region_medium <- ggplot(MAG_metadata_extract, 
#                           aes(x = medium, y = MGE)) +
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
#   labs(x="Medium",y="MGE number")+
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
# p_Region_medium

# g2 <- plot_grid(p_Region_ecosystem,p_Region_medium,
#                 labels = letters[2:3], ncol = 2
#                 #rel_heights = c(1, 1, 1.5),
#                 # label_x = 0.09,            # 标签x坐标（1为最右侧）
#                 #  label_y = 0.9             # 标签y坐标（1为最上方）
# )
# g2

# cairo_pdf(filename = "1_region2ecosystem_medium.pdf" , width = 7.5,height = 2.4)
# g2

# setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
# dev.off()





#比较glacier，plateau，plain的MGE变化，DS变化；DS与MGE的correlation####

##三个region的MGE比较####
library(dplyr)

MAG_metadata_glacier_plateau_plain <- select(MAG_metadata,genome,ecosystem,medium,Region,MGE)
MAG_metadata_glacier_plateau_plain$Region <- as.character(MAG_metadata_glacier_plateau_plain$Region) 
MAG_metadata_glacier_plateau_plain$Region <- factor(MAG_metadata_glacier_plateau_plain$Region,
                                                    levels = c("Glacier","Plateau","Plain"))
# MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
#                                             MAG_metadata_plateau_plain$Region)

md_p_val0 <- MAG_metadata_glacier_plateau_plain %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(MGE ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val0$y.position <- md_p_val0$y.position - 100

p_glacier_plateau_plain <- 
  ggplot(MAG_metadata_glacier_plateau_plain, 
         aes(x = Region, y = MGE)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  # geom_half_point(aes(color=Region),
  #                 side = "r",size = 0.5,
  #                 transformation_params = list(height = 0,width = 0.001,seed = 2))+
  geom_jitter(aes(color=Region), width = 0.1, size = 0.5, alpha = 0.6)+
  stat_pvalue_manual(md_p_val0,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="MGE number (per MAG)")+
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
  coord_cartesian(ylim = c(-1, 80))
p_glacier_plateau_plain

# setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
# ggsave("1_MGEnum_glacier_plateau_plain_comprison.pdf",p_glacier_plateau_plain, 
#        width = 2.6, height = 2.78, units = "in")
ggsave("1_MGEnum_glacier_plateau_plain_comprison_new.pdf",p_glacier_plateau_plain, 
       width = 2.76, height = 2.4, units = "in")

write.table(MAG_metadata_glacier_plateau_plain, "MAG_metadata_glacier_plateau_plain.tsv",
  sep = '\t', quote = F, row.names = F)

##3个介质的MGE比较####
MAG_metadata$medium2 <- ifelse(MAG_metadata$Type == 'Cryoconite',"Glacier",MAG_metadata$Type)
MAG_metadata$medium2 <- ifelse(MAG_metadata$medium2 == 'Ice',"Glacier",MAG_metadata$medium2)
MAG_metadata$medium2 <- ifelse(MAG_metadata$medium2 == 'Snow',"Glacier",MAG_metadata$medium2)
MAG_metadata$medium2 <- ifelse(MAG_metadata$medium2 == 'particle',"Plankton",MAG_metadata$medium2)
MAG_metadata$medium2 <- ifelse(MAG_metadata$medium2 == 'water',"Plankton",MAG_metadata$medium2)
MAG_metadata$medium2 <- ifelse(MAG_metadata$medium2 == 'sediment',"Benthon",MAG_metadata$medium2)

MAG_metadata$medium2 <- as.character(MAG_metadata$medium2) 
MAG_metadata$medium2 <- factor(MAG_metadata$medium2,
                               levels = c("Glacier","Benthon","Plankton"))
# MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
#                                             MAG_metadata_plateau_plain$Region)

md_p_val_medium <- MAG_metadata %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(MGE ~ medium2) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="medium2",dodge=0.5)
md_p_val_medium$y.position <- md_p_val_medium$y.position - 100

p_glacier_Benthon_Plankton <- 
  ggplot(MAG_metadata, 
         aes(x = medium2, y = MGE)) +
  geom_half_boxplot(aes(fill = medium2), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_jitter(aes(color=Region), width = 0.1, size = 0.5, alpha = 0.6)+
  stat_pvalue_manual(md_p_val_medium,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="MGE number (per MAG)")+
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
  coord_cartesian(ylim = c(-1, 80))
p_glacier_Benthon_Plankton

write.table(MAG_metadata, "MAG_metadata_mge.tsv",sep = "\t",quote = F, row.names = F)



##3个region的DS比较####
library(rstatix)

MAG_DSsum <- MAG_DSsum[MAG_DSsum$genome %in% MAG_metadata$genome,]
MAG_metadata <- merge(MAG_metadata, MAG_DSsum, by = "genome", all = T)
MAG_metadata[,82][is.na(MAG_metadata[,82])] <- 0

MAG_DS_metadata_glacier_plateau_plain <- select(MAG_metadata,genome,ecosystem,medium,Region,MGE,DSsum)
MAG_DS_metadata_glacier_plateau_plain$Region <- as.character(MAG_DS_metadata_glacier_plateau_plain$Region) 
MAG_DS_metadata_glacier_plateau_plain$Region <- factor(MAG_DS_metadata_glacier_plateau_plain$Region,
                                                    levels = c("Glacier","Plateau","Plain"))
# MAG_metadata_plateau_plain$Region <- ifelse(MAG_metadata_plateau_plain$Region == "Glacier", "Plateau", 
#                                             MAG_metadata_plateau_plain$Region)

md_p_val_DS <- MAG_DS_metadata_glacier_plateau_plain %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(DSsum ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val_DS$y.position <- md_p_val_DS$y.position - 40

p_DS_glacier_plateau_plain <- 
  ggplot(MAG_DS_metadata_glacier_plateau_plain, 
         aes(x = Region, y = DSsum)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  geom_half_point(aes(color=Region),
                  side = "r",size = 0.5,
                  transformation_params = list(height = 0,width = 0.001,seed = 2))+
  stat_pvalue_manual(md_p_val_DS,label = "p.adj.signif",label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="DS number (per MAG)")+
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
  coord_cartesian(ylim = c(-1, 20))
p_DS_glacier_plateau_plain

setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
ggsave("1_DSnum_glacier_plateau_plain_comprison.pdf",p_DS_glacier_plateau_plain, 
       width = 2.6, height = 2.78, units = "in")

##DS与MGE的相关####

cor.test(MAG_DS_metadata_glacier_plateau_plain$MGE, 
         MAG_DS_metadata_glacier_plateau_plain$DSsum, method = "spearman")
#ρ0.41, p-value < 2.2e-16
p_DSnum_MGE_cor <- 
ggplot(MAG_DS_metadata_glacier_plateau_plain,aes(x = MGE, y = DSsum, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 150, y = 40, 
           label = "Spearman: ρ = 0.41, P < 2.2e-16", vjust = 1.5, hjust = 1)+
  #scale_y_continuous(breaks = seq(0, max(DHR[-c(40,131),]$ds_richness), by = 2))+
  labs(x='MGE number (per MAG)',y='DS number (per MAG)')+
  scale_color_manual(values = c( "#9c89b8",
                                 "indianred1"
  )) +  # 手动设置每个category的颜色
  mytheme
p_DSnum_MGE_cor

MAG_metadata$DS_densityperkb <- MAG_metadata$DSsum/MAG_metadata$Genome_Size*1000
MAG_metadata$DS_densitylog <- log10(MAG_metadata$DS_densityperkb+1)*1000

cor.test(MAG_metadata$MGE, 
         MAG_metadata$DS_densitylog, method = "spearman")
#ρ0.0.34, p-value < 2.2e-16

cor.test(MAG_metadata$MGE, 
         MAG_metadata$DSsum, method = "spearman")

p_DSdensity_MGE_cor <- 
ggplot(MAG_metadata,aes(x = MGE, y = DS_densitylog, color = "#9c89b8"))+
  geom_point(size=2)+
  geom_smooth(method = 'lm',se = T)+
  # stat_poly_eq(use_label(c("eq", "adj.R2", "p.value.label")),
  #              formula = y ~ x,  parse = TRUE,
  #              size = 5 #公式字体大小
  #              # label.x = 0.05,  #位置 ，0-1之间的比例
  #              # label.y = 0.95
  # )+
  annotate("text", x = 150, y = 5, 
           label = "Spearman: ρ = 0.34, P < 2.2e-16", vjust = 1.5, hjust = 1)+
  #scale_y_continuous(breaks = seq(0, max(DHR[-c(40,131),]$ds_richness), by = 2))+
  labs(x='MGE number (per MAG)',y='log10 (DS density +1) × 1e-3')+
  scale_color_manual(values = c( "#89b0ae","#9c89b8"
  )) +  # 手动设置每个category的颜色
  mytheme
p_DSdensity_MGE_cor

write.table(MAG_metadata, "MAG_metadata_mge.tsv", sep = "\t", quote = F, row.names = F)

library(Cairo)
library(cowplot)
cairo_pdf(filename = "1_MGE4_correlation2.pdf" , width = 7,height = 6.5)
g1 <- plot_grid(p_glacier_plateau_plain,p_glacier_Benthon_Plankton,
                p_DSnum_MGE_cor, p_DSdensity_MGE_cor,
                labels = letters[1:4], ncol = 2
                #rel_widths = c(3, 2)
                # label_x = 0.09,            # 标签x坐标（1为最右侧）
                #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g1

# setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
dev.off()

save.image("./01_MGE_cargo_and_defense_density.RData")
load("./01_MGE_cargo_and_defense_density.RData")


