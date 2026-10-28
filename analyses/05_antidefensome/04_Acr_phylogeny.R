# =============================================================================
# Glacier-to-river -- antidefensome / Anti-CRISPR phylogeny annotation
# =============================================================================
# Annotate the maximum-likelihood Acr tree (44 predicted + 339 validated Acr
# proteins) with iTOL colour strips.
#
# Manuscript      : Fig. S16
# Main outputs    : iTOL annotation files (hub)
# Original script : 20gaoyuan/15_antidefense/3Acr_tree/Acr_tree_annotation.R
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

# Tong Zhou, Kuidong Xu, Feng Zhao et. al. itol.toolkit accelerates working with 
# iTOL (Interactive Tree Of Life) by an automated generation of annotation files, 
# Bioinformatics, 2023;, btad339, https://doi.org/10.1093/bioinformatics/btad339

#参考学习https://tongzhou2017.github.io/itol.toolkit/articles/Get_Start.html
rm(list = ls())
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
library(itol.toolkit)
library(data.table) # 用于元数据文件读取
library(ape) # 用于树文件读写
library(stringr) # 用于字符串处理
library(dplyr)

setwd(file.path(GTR_DATA, "20gaoyuan/15_antidefense/3Acr_tree"))
# 假设您的树存储在一个名为"tree"的变量中
tree <- "383acr_aln.trim.afa.treefile"

phylo <- ape::read.tree(tree)

hub <- create_hub(tree)

#1.计算树的特征####


#2. 准备注释文件####
acr_ref <- readxl::read_excel("./Acr_protein.xlsx")
myacr164 <- read.delim2("44_viral_acr.tsv")
#spacertargetacr27 <- read.delim2("spacer_targeted_prophagewithAnti_CRISPR.txt",header = F)

#2.1 区分我们的与参考数据库的

acr_ref_typeinhibited <- select(acr_ref,ID,TypeInhibited)
myacr_type <- select(myacr164,hit_id,type)
myacr_type$hit_id <- gsub("\\+", "_", myacr_type$hit_id)
#myacr_spacertarget <- select(spacertargetacr27,V2,V23)

#myacr_spacertarget$V23 <- "Spacer"
#myacr_type[myacr_type$hit_id%in%myacr_spacertarget$V2,]$type <- "Spacer"

colnames(acr_ref_typeinhibited) <- c('id','type')
colnames(myacr_type) <- c('id','type')
df_idtype <- rbind(acr_ref_typeinhibited,myacr_type)



colorss <- c(
#   "#17becf",  # 青色
#  "#ff7f0e",  # 橙色
  "#F5E320",  # 红色
  "#CC281C",  # 绿色
  "#D1CB7E",  # 紫色
  "#5AAAB2",  # 棕色
  "#E77F5A",  # 粉色
  "#A0E6C5",  # 灰色
  "#70A47F",  # 橄榄色
# "#1f77b4",  # 蓝色
  "#D68112",  # 浅蓝色
  "#C9A79E",  # 浅橙色
#  "#98df8a",  # 浅绿色
#  "#ff9896",  # 浅红色
  "#6C6B6F",  # 浅紫色
  "#CED2DD",  # 浅棕色
#  "#f7b6d2",  # 浅粉色
  "#4E616A",  # 浅橄榄色
#  "#9edae5",  # 浅青色
  "#ED6764",  # 深蓝色
  "#74A0A4"   # 深绿色
)

colorss <- c(
  "#A6E3E9", "#D1CB7E", "#50C878", "#E77F5A", 
  "#0047AB", "#0087BD", "#6A5ACD", "#4682B4",
  "#01796F", "#87CEEB", "#4B0082", "#00A693",
  "#6C6B6F", "#C9A79E"
)

element_counts <- df_idtype %>%
  count(type)

element_count_sorted <- element_counts[order(-element_counts$n), ]
element_count_sorted$color <- c(colorss)

#
df_idtype <- merge(df_idtype,element_count_sorted,by.x = "type",by.y = "type")
df_idtype <- df_idtype[, c("id", "type","color")]




unit_1 <- create_unit(data = df_idtype, 
                      key = "Acrtype_color_range", 
                      type = "TREE_COLORS", 
                      subtype = "range", 
                     # color = colorss,  #需要的颜色超过了集合中有的颜色
                      # line_type = c(rep("normal",4),"dashed"),
                      # size_factor = 5, 
                      tree = phylo)
write_unit(unit_1)

#2.2 对被spacer靶向的病毒上的Acr进行标记，01矩阵
acr_ref <- readxl::read_excel("./Acr_pAcr_Sequence_Files/Acr_protein.xlsx")
myacr164 <- read.delim2("prophage_164Anti_CRISPR.txt")
spacertargetacr27 <- read.delim2("spacer_targeted_prophagewithAnti_CRISPR.txt",header = F)

acr_ref_typeinhibited <- select(acr_ref,ID,TypeInhibited)
myacr_type <- select(myacr164,hit_id,type)
myacr_spacertarget <- select(spacertargetacr27,V2,V23)

myacr_spacertarget$V23 <- "Spacer"
myacr_type[myacr_type$hit_id%in%myacr_spacertarget$V2,]$type <- 1
myacr_type[myacr_type$type!=1,]$type <- 0
acr_ref_typeinhibited$TypeInhibited <- 0

colnames(acr_ref_typeinhibited) <- c('id','type')
colnames(myacr_type) <- c('id','type')
df_spacer <- rbind(acr_ref_typeinhibited,myacr_type)

unit_2 <- create_unit(data = df_spacer,
                      key = 'spacer_27acr_heatmap',
                      type = 'DATASET_EXTERNALSHAPE',
                      tree = tree
)
write_unit(unit_2)

