# =============================================================================
# Glacier-to-river -- MAG catalog / Archaeal MAG tree annotation
# =============================================================================
# Same iTOL annotation workflow as step 04, for the archaeal MAG tree.
#
# Manuscript      : Fig. 1b
# Main outputs    : iTOL annotation files (hub)
# Original script : 20gaoyuan/11_42mag_ar_tree/1_arMAG_tree_annotation.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "20gaoyuan/11_42mag_ar_tree"))


# Tong Zhou, Kuidong Xu, Feng Zhao et. al. itol.toolkit accelerates working with 
# iTOL (Interactive Tree Of Life) by an automated generation of annotation files, 
# Bioinformatics, 2023;, btad339, https://doi.org/10.1093/bioinformatics/btad339

#参考学习https://tongzhou2017.github.io/itol.toolkit/articles/Get_Start.html
rm(list = ls())
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")
library(itol.toolkit)
library(data.table) # 用于元数据文件读取
library(ape) # 用于树文件读写
library(stringr) # 用于字符串处理
library(dplyr)


# 假设您的树存储在一个名为"tree"的变量中
tree <- "bin5455ar.unrooted.tree"
phylo <- ape::read.tree(tree)
hub <- create_hub(tree)

#1.计算树的特征####
# node_depths <- node.depth(tree)
# print(node_depths)
# 
# #节点支持度
# boot_result <- boot.trees(tree, stat="max", B=100) # B是自举次数
#树长
tree_length <- sum(phylo$edge.length)

node_depths <- node.depth(phylo)#

#计算每个节点的系统发育深度
phylo_depth <- vcv.phylo(phylo)#消耗内存巨大，算一次就行
phylo_depth_diag <- diag(phylo_depth) %>% as.data.frame() %>% setNames('phylo_depth')
write.csv(phylo_depth_diag,"phylo_depth_diag.csv")

# 进化树有多少末端分类单元？
Ntip(phylo)
# 查看末端分类单元名称， 取前5个
phylo$tip.label[1:5]

#2.准备itol树的注释文件

#2.1 物种注释信息####
mags_tax <- read.delim2('gtdbtk.ar53.summary.tsv')
# 假设 df 是你的数据框，MAGs 是其中的一列

#处理物种信息
taxonoxy <- function(tax){
  tax <- tax[,1:2]
  max_cols <- max(sapply(strsplit(tax$classification, ";"), length))# 确定最大的拆分列数
  col_names <- paste0("col", 1:max_cols)# 生成列名
  library(tidyr)
  taxonomy <- tax %>%
    separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列
}

mags_tax <- taxonoxy(mags_tax)

# 使用 sub 函数替换掉第三个 _ 及其后面的内容
#mags_tax$col2 <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", mags_tax$col2)
# 删除p__前缀
mags_tax$col2 <- gsub("^p__", "", mags_tax$col2)

#按门设置颜色
library(dplyr)
df_phycolor <- mags_tax %>% select(user_genome, col2)
set.seed(666)

#手动为门添加颜色
colorss <- c("#e0afa0","#bbd0ff","#a4c3b2")
element_counts <- df_phycolor %>%
  count(col2)

element_count_sorted <- element_counts[order(-element_counts$n), ]
element_count_sorted$color <- c(colorss)

#
df_phycolor <- merge(df_phycolor,element_count_sorted,by.x = "col2",by.y = "col2")
df_phycolor <- df_phycolor[, c("user_genome", "col2","color")]

#subtype = "range"才能大面积改变颜色
unit_2 <- create_unit(data = df_phycolor, 
                      key = "phylo_color_range", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                      #color = "wesanderson",  #需要的颜色超过了集合中有的颜色
                      # line_type = c(rep("normal",4),"dashed"),
                      # size_factor = 5, 
                      tree = phylo)
write_unit(unit_2)

unit_3 <- create_unit(data = df_phycolor, 
                      key = "phylo_color_branch", 
                      type = "TREE_COLORS", 
                      subtype = "branch", 
                      #color = "wesanderson",  #需要的颜色超过了集合中有的颜色
                      # line_type = c(rep("normal",4),"dashed"),
                      # size_factor = 5, 
                      tree = phylo)
write_unit(unit_3)

#2.2 MAG size####
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

MAG_metadata <- read.csv(file.path(GTR_DATA, "20gaoyuan/01MAG_info/5455MAG_metadata.csv"))
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

MAG_DSsum <- rowSums(count_matrix) %>% as.data.frame() %>% setnames("DSsum")
MAG_DSsum$genome <- rownames(MAG_DSsum)

colnames(MAG_quality)[1] <- "genome"
MAG_quality <- left_join(MAG_quality, MAG_DSsum, by = "genome")

MAG_quality$MAG_size <- MAG_quality$Genome_Size/1000000
MAG_quality$DS_density <- MAG_quality$DSsum/MAG_quality$Genome_Size*1000000
MAG_quality$magsize <- floor(MAG_quality$MAG_size)

MAG_quality[is.na(MAG_quality)] <- 0


unit_7 <- create_unit(data = MAG_quality[,c(1,16)],
                      key = "mag_length", 
                      type = "DATASET_MULTIBAR", 
                      tree = tree)
write_unit(unit_7)

#2.3 添加主要DSs的信息####

## 计算每一列非0数字所占的百分比并输出到新数据框
dss <- MAG_quality[,c(1,2)]
count_matrix$genome <- rownames(count_matrix)
dss <- left_join(dss, count_matrix, by = "genome")
count_matrix <- count_matrix[,-189]

rownames(dss) <- dss$genome
dss <- dss[,-c(1,2)]
dss[is.na(dss)] <- 0

percent_nonzero <- colSums(dss != 0) / nrow(dss) * 100
new_df <- data.frame(Column = names(percent_nonzero), Percentage = percent_nonzero)
# 按照值的大小对数据集进行降序排列
ds_persentage <- new_df[order(-new_df$Percentage), ]
ds_persentage <- ds_persentage[ds_persentage$Percentage>=2,]

dss_33 <- dss[,colnames(dss) %in% ds_persentage$Column]
#将dss_33中大于0的值全部设置为1
dss_33 <- dss_33 %>%
  mutate(across(everything(), ~if_else(. > 0, 1, .)))

dss_33$user_genome <- rownames(dss_33)
dss_33 <- dss_33 %>%
  select(user_genome, everything())

unit_5 <- create_unit(data = dss_33,
                      key = 'dss33_p-a_heatmap',
                      type = 'DATASET_HEATMAP',
                      tree = tree
)
unit_5@specific_themes$heatmap$tree$tree_display <- 0 #不聚类
unit_5@specific_themes$heatmap$color$min <- "#ffffff"
unit_5@specific_themes$heatmap$color$max <- "#8DD4DA"
unit_5@common_themes$basic_theme$margin <- 10

write_unit(unit_5)
#

#2.4添加DS和DG的密度####

##DS密度####

#ds_density <- read.csv("ds_density_tax.csv")
#ds_density$genome <- gsub(".fasta","",ds_density$genome)


ds_density_d <- MAG_quality[,c('genome',"DS_density")]
ds_density_d$ds_density2 <- ds_density_d$DS_density

unit_6 <- create_unit(data = ds_density_d,
                      key = 'ds_density_heatmap',
                      type = 'DATASET_GRADIENT',
                      method = 'mean',
                      tree = tree
)

unit_6@specific_themes$heatmap$color$min <- "#f3f6f4"
unit_6@specific_themes$heatmap$color$mid <- "#74c69d"
unit_6@specific_themes$heatmap$color$max <- "#1b4332"
unit_6@common_themes$basic_theme$margin <- 10

write_unit(unit_6)

##DG density####
#读取含有dss的mags的数据
setwd(file.path(GTR_DATA, "20gaoyuan/11_42mag_ar_tree"))
dgs <- read.delim2('merged_5467MAG_defense_finder_genes.tsv',sep = '\t',
                   quote = "",header = T
)

dgs <- dgs[dgs$activity == "Defense",]
dgs$genome <- sub("-[0-9]+$", "", dgs$replicon)
DG_num <- dgs %>% 
  count(genome, name = "DGnum")  # 简写方式

DG_num <- DG_num[DG_num$genome %in% MAG_metadata$genome,]

MAG_quality <- left_join(MAG_quality, DG_num, by = "genome")
MAG_quality[is.na(MAG_quality)] <- 0
MAG_quality$DG_density <- MAG_quality$DGnum/MAG_quality$Genome_Size*1000000

dg_density_d <- MAG_quality[,c("genome","DG_density")]
dg_density_d$dg_density2 <- dg_density_d$DG_density

setwd(file.path(GTR_DATA, "20gaoyuan/11_42mag_ar_tree"))
unit_7 <- create_unit(data = dg_density_d,
                      key = 'dg_density_heatmap',
                      type = 'DATASET_GRADIENT',
                      method = 'mean',
                      tree = tree)
unit_7@specific_themes$heatmap$color$min <- "#f3f6f4"
unit_7@specific_themes$heatmap$color$mid <- "#8DB3D8"
unit_7@specific_themes$heatmap$color$max <- "#1969AF"
unit_7@common_themes$basic_theme$margin <- 10
write_unit(unit_7)


#2.5 DI数量####
DI <- read.delim2("5455MAG_defense_islands_summary.csv",check.names = F)
library(stringr)
DI$genome <- str_extract(DI$Contig, ".*bin\\.\\d+")
DI <- DI[DI$genome %in% MAG_quality$genome,]
DI_num <- DI %>% 
  count(genome, name = "DInum")  # 简写方式
MAG_quality <- left_join(MAG_quality, DI_num, by = "genome")
MAG_quality[is.na(MAG_quality)] <- 0

DI_num_d <- MAG_quality[,c("genome","DInum")]
DI_num_d$DI_num2 <- DI_num_d$DInum

setwd(file.path(GTR_DATA, "20gaoyuan/11_42mag_ar_tree"))
unit_8 <- create_unit(data = DI_num_d,
                      key = 'DI_num_heatmap',
                      type = 'DATASET_GRADIENT',
                      method = 'mean',
                      tree = tree)
unit_8@specific_themes$heatmap$color$min <- "#f3f6f4"
unit_8@specific_themes$heatmap$color$mid <- "#faac46"
unit_8@specific_themes$heatmap$color$max <- "#ff7b00"
unit_8@common_themes$basic_theme$margin <- 10
write_unit(unit_8)

#2.6 MGE density####
MGE_num <- read.delim2(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs/genome_MGE_number.tsv"))
MAG_quality <- left_join(MAG_quality, MGE_num, by = "genome")
MAG_quality[is.na(MAG_quality)] <- 0
MAG_quality$MGE_density <- MAG_quality$MGE/MAG_quality$Genome_Size*1000000

MGE_density_d <- MAG_quality[,c("genome","MGE_density")]
MGE_density_d$MGE_density2 <- MGE_density_d$MGE_density
setwd(file.path(GTR_DATA, "20gaoyuan/11_42mag_ar_tree"))
unit_9 <- create_unit(data = MGE_density_d,
                      key = 'MGE_density_heatmap',
                      type = 'DATASET_GRADIENT',
                      method = 'mean',
                      tree = tree)
unit_9@specific_themes$heatmap$color$min <- "#f3f6f4"
unit_9@specific_themes$heatmap$color$mid <- "#b185dd"
unit_9@specific_themes$heatmap$color$max <- "#6247aa"
unit_9@common_themes$basic_theme$margin <- 10
write_unit(unit_9)


#3. Ecosystem Region Medium(Type)注释####
##3.1 Ecosystem####

ecosystem_color <- 
  c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
  "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA") %>% as.data.frame() %>% setNames("color")
ecosystem_color$ecosystem <- rownames(ecosystem_color)
ecosystem_df <- MAG_metadata[,c("genome","ecosystem")]
ecosystem_df <- left_join(ecosystem_df, ecosystem_color, by = "ecosystem")

unit_31 <- create_unit(data = ecosystem_df,
                       key = "All_MAG_ecosystem",
                       type = "DATASET_COLORSTRIP",
                       tree = tree)
unit_31@common_themes$basic_theme$margin <- 10
write_unit(unit_31)

##3.2 Region####

region_color <- 
  c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
    "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
    "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  ) %>% as.data.frame() %>% setNames("color")
region_color$Region <- rownames(region_color)
region_df <- MAG_metadata[,c("genome","Region")]
region_df <- left_join(region_df, region_color, by = "Region")

unit_32 <- create_unit(data = region_df,
                       key = "All_MAG_region",
                       type = "DATASET_COLORSTRIP",
                       tree = tree)
unit_32@common_themes$basic_theme$margin <- 10
write_unit(unit_32)

##3.3 Type####

type_color <- 
  c("Snow"="#E64B35FF","Ice"="#4DBBD5FF",
    "Cryoconite"="#00A087FF","water"= "#3C5488FF",
    "particle"= "#F39B7FFF","sediment"="#8491B4FF") %>%
  as.data.frame() %>% setNames("color")
type_color$Type <- rownames(type_color)
type_df <- MAG_metadata[,c("genome","Type")]
type_df <- left_join(type_df, type_color, by = "Type")

unit_33 <- create_unit(data = type_df,
                       key = "All_MAG_type",
                       type = "DATASET_COLORSTRIP",
                       tree = tree)
unit_33@common_themes$basic_theme$margin <- 10
write_unit(unit_33)





#Storing drawing information
hub <- hub +
  unit_3 +
  unit_6 +
  unit_7 +
  unit_5

write_hub(hub,getwd())




