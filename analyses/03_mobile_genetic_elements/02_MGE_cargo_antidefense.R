# =============================================================================
# Glacier-to-river -- mobile genetic elements / Anti-defense genes carried by MGEs
# =============================================================================
# Count anti-defense genes per MGE type, their density and their O/E
# localisation ratios.
#
# Manuscript      : Fig. S15
# Main outputs    : 2_mge_antids_density_comparison.pdf, 02_MGE_cargo_antidefense.RData
# Original script : 20gaoyuan/12_5455MAG_MGEs/2_MGE_antidefense.R
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


#MGE数量统计、生境分布、携带的antiDS。

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


#MGE与antiDS关系##########
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

DF_defense <- DF_res[DF_res$activity == "Antidefense",]

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
DG <- DG[DG$activity == "Antidefense",]
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
integron_DS_num <- data.frame(genome = integron_DS_statistic$genome.x,
                              integron_DS_num = integron_DS_statistic[,2])

ICE_DS_num <- data.frame(genome = ICE_DS_statistic$genome,
                         ICE_DS_num = ICE_DS_statistic[,-1])

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

mge_ds_density <- mge_ds_density[mge_ds_density$location %in% c("Plasmid","Phage","Prophage"),]
mge_ds_density <- mge_ds_density[mge_ds_density$dsnumber>0,]
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
  #geom_jitter(aes(group=1,color = location),width = 0.2, alpha = 0.3)+  # 显示异常值点
  geom_violin(trim = T, aes(fill = location),show.legend = F,linewidth=0.2) + # 绘制小提琴图
  geom_boxplot(aes(fill = location),width = 0.2,outlier.shape = NA, staplewidth = 0.5,linewidth=0.2) +  # 隐藏箱线图中的异常值点
  geom_jitter(aes(group=1,color = 'black'), fill = 'black',width = 0.2, alpha = 0.3)+  # 显示异常值点
  # stat_compare_means(comparisons = list(c("1", "2"),c("SS", "WA"),c("SS", "WS"),
  #                                       c("SA", "WA"),c("SA", "WS"),
  #                                       c("WA", "WS")
  # ),
  # method = "wilcox.test", test = "p.adjust",label = "p.signif")+
  # stat_pvalue_manual(md_p_val,label = "p.adj.signif",hide.ns=T,
  #                    tip.length = 0,label.size = 5,color="black")+
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
  coord_cartesian(ylim = c(0,0.3))+
  labs(x=NULL,y='Antidefense density (per kb)')+
  mytheme
p_mge_ds_density

setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
ggsave("2_mge_antids_density_comparison.pdf",p_mge_ds_density,width = 3.3,height = 3,units = "in")

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

ds_family_integron <- sum(integron_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("Integron")
ds_family_integron$dsfamily <- colnames(integron_DS_statistic)[2]

ds_family_ICE <- sum(ICE_DS_statistic[,-1]) %>% as.data.frame() %>% setNames("ICE")
ds_family_ICE$dsfamily <- colnames(ICE_DS_statistic)[2]

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
  full_join(ds_family_prophage,by = 'dsfamily')


ds_family_bind[is.na(ds_family_bind)] <- 0

ds_family_bind$Chromosome <- ds_family_bind$All - 
  ds_family_bind$Plasmid - 
  ds_family_bind$Prophage -
  ds_family_bind$Phage

ds_family_final <- select(ds_family_bind,dsfamily,Chromosome,Plasmid,Phage,
                          Prophage) %>% as.data.frame()

rownames(ds_family_final) <- ds_family_final$dsfamily
ds_family_final <- ds_family_final[,-1]
sum(rowSums(ds_family_final)>20)
#ds_family_final <- ds_family_final[rowSums(ds_family_final)>20,]


oe <- chisq.test(ds_family_final,simulate.p.value = T)

oe_bi <- oe$observed/oe$expected %>% as.data.frame()
oe_bi <- na.omit(oe_bi)
#oe_bi <- oe_bi[oe_bi$Chromosome!=0,]
#oe_bi <- oe_bi[!rowSums(oe_bi[,-1]) == 0,]

#绘制热图
library(ComplexHeatmap)
library(circlize)  # 用于 colorRamp2 函数

# 定义颜色映射函数
col_fun <- colorRamp2(c(0,1, 2), c("white","#B6B8D8", "#FF7F0E"))  # 粉色到绿色的渐变
setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
cairo_pdf(width = 4.5,height = 2.5,filename = "2_MGE_antiDS_oe.pdf")
p_mge_ds_oe <- Heatmap(t(oe_bi),
                       col = col_fun,
                       cluster_rows = F,
                       cluster_columns = T,
                       name = 'O/E',
                       row_title = NULL,
                       row_split = colnames(oe_bi),
                       row_names_side = "left", # 行标题显示在左侧
                       row_names_gp = gpar(fontsize = 10),
                       column_names_gp = gpar(fontsize = 10),
                       column_names_rot = 90,
                       rect_gp = gpar(col = "white", lwd = 0.5),    ##热图矩阵线设置
                       heatmap_legend_param = list(
                         title = "O/E", at = c(0,1, 2),#direction = "horizontal",
                         labels = c("0","1", "≥2"))
)
p_mge_ds_oe
dev.off()
setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))

#MGE携带antidefense的热图####
library(ComplexHeatmap)
library(circlize)  # 用于颜色映射

# 定义颜色映射（数值越大颜色越深）
col_fun <- colorRamp2(c(0, max(ds_family_final)), c("white", "firebrick3"))

ds_family_bind2 <- ds_family_chromosome %>%
  full_join(ds_family_phage,by = 'dsfamily') %>%
  full_join(ds_family_plasmid,by = 'dsfamily') %>%
  full_join(ds_family_prophage,by = 'dsfamily') %>%
  full_join(ds_family_integron,by = 'dsfamily') %>%
  full_join(ds_family_ICE,by = 'dsfamily') %>% 
  full_join(ds_family_IS,by = 'dsfamily')

ds_family_bind2[is.na(ds_family_bind2)] <- 0

ds_family_bind2$Chromosome <- ds_family_bind2$All - 
  ds_family_bind2$Plasmid - 
  ds_family_bind2$Prophage -
  ds_family_bind2$Phage -
  ds_family_bind2$Integron -
  ds_family_bind2$ICE -
  ds_family_bind2$IS

ds_family_final2 <- select(ds_family_bind2,dsfamily,Chromosome,Plasmid,Phage,
                          Prophage,Integron,ICE) %>% as.data.frame()

rownames(ds_family_final2) <- ds_family_final2$dsfamily
ds_family_final2 <- ds_family_final2[,-1]
sum(rowSums(ds_family_final2)>20)

ds_family_final_heatmap <- ds_family_final2 %>% as.matrix()

p_heatmap_antidefense_MGE <- 
Heatmap(ds_family_final_heatmap,
        name = "Count",                # 图例标题
        col = col_fun,                 # 颜色函数
        rect_gp = gpar(col = "gray"),  # 单元格边框
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(ds_family_final_heatmap[i, j], x, y, gp = gpar(fontsize = 10))  # 显示数值
        },
        row_names_side = "left",
        column_names_rot = 0,
        #row_title = "Defense System Family",
        #column_title = "Genomic Context",
        heatmap_legend_param = list(title_position = "leftcenter-rot"))
p_heatmap_antidefense_MGE

setwd(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs"))
cairo_pdf(width = 5,height = 2.8,filename = "2_heatmap_antidefense_MGE_2.pdf")
p_heatmap_antidefense_MGE
dev.off()

write.table(ds_family_final_heatmap, "ds_family_final_heatmap.tsv",
  sep = "\t", quote = F)

write.table(oe_bi, "oe_bi_anti_mge.tsv",
  sep = "\t", quote = F)

write.table(mge_ds_density, "mge_ds_density.tsv",
  sep = "\t", quote = F)

save.image("./02_MGE_cargo_antidefense.RData")

load("./02_MGE_cargo_antidefense.RData")











