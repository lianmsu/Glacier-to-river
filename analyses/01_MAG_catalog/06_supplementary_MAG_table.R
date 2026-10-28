# =============================================================================
# Glacier-to-river -- MAG catalog / Supplementary MAG table
# =============================================================================
# Merge MAG metadata with the CheckM quality report and the GTDB taxonomy to
# build the per-MAG supplementary table.
#
# Manuscript      : Supplementary Data
# Main outputs    : supplementary MAG table
# Original script : 20gaoyuan/Supplementary_tables/1MAG_info.R
# -----------------------------------------------------------------------------
# Input data and output files are resolved from GTR_DATA, the root of the local
# analysis directory (its layout mirrors the original working copy). To run the
# script elsewhere, set the environment variable instead of editing this file:
#
#     Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
#
# =============================================================================
GTR_DATA <- Sys.getenv("GTR_DATA", unset = "D:/02_pku/01_research/03_plateau_defensome")

setwd(file.path(GTR_DATA, "20gaoyuan/Supplementary_tables"))

library(dplyr)
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

MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]

MAG_metadata <- merge(MAG_metadata,MAG_quality, 
                      by.x = "genome",by.y = "Name",all = T)

bac_tax <- read.delim2("gtdbtk.bac120.summary.tsv")
ar_tax <- read.delim2("gtdbtk.ar53.summary.tsv")

all_MAG_tax <- rbind(bac_tax, ar_tax)

#处理物种信息
taxonoxy <- function(tax){
  tax <- tax[,1:2]
  max_cols <- max(sapply(strsplit(tax$classification, ";"), length))# 确定最大的拆分列数
  col_names <- paste0("col", 1:max_cols)# 生成列名
  library(tidyr)
  taxonomy <- tax %>%
    separate(classification, into = col_names, sep = ";", fill = "right")# 拆分字符串到列
}

all_MAG_tax <- taxonoxy(all_MAG_tax)
colnames(all_MAG_tax) <- c("genome","Domain","Phylum","Class",
                           "Order","Family","Genus","Species")
all_MAG_tax$Domain <- gsub("^d__", "", all_MAG_tax$Domain)
all_MAG_tax$Phylum <- gsub("^p__", "", all_MAG_tax$Phylum)
all_MAG_tax$Class <- gsub("^c__", "", all_MAG_tax$Class)
all_MAG_tax$Order <- gsub("^o__", "", all_MAG_tax$Order)
all_MAG_tax$Family <- gsub("^f__", "", all_MAG_tax$Family)
all_MAG_tax$Genus <- gsub("^g__", "", all_MAG_tax$Genus)
all_MAG_tax$Species <- gsub("^s__", "", all_MAG_tax$Species)

MAG_metadata <- merge(MAG_metadata,all_MAG_tax, 
                      by.x = "genome",by.y = "genome",all = T)

write.csv(MAG_metadata,"MAG5455_metadata.csv",quote = F)




##########

DF_res <- read.delim2("merged_5467MAG_defense_finder_systems.tsv",
                      check.names = F,
                      sep = '\t')
DF_defense <- DF_res[DF_res$activity == "Defense",]
library(dplyr)
library(tidyr)
library(stringr)
DF_defense <- DF_defense %>%
  mutate(
    genome = str_extract(sys_id, ".*bin\\.\\d+")  # 包含"bin.数字"
  )
DF_defense <- DF_defense[DF_defense$genome %in% MAG_metadata$genome,]
DF_defense <- select(DF_defense, genome, everything())

write.table(DF_defense,"DF_defense.tsv",quote = F,sep = '\t')

###############
DI <- read.delim2("5455MAG_defense_islands_summary_no_anti.csv",check.names = F)
write.csv(DI,"DI.csv",quote = F)

DI_no <- DI_no[DI_no$genome %in% MAG_metadata$genome,]
write.table(DI_no,"DI.tsv",quote = F,sep = "\t")


#############

write.table(all_55265_phages_DF_anti,"all_55265_phages_DF_anti.tsv",quote = F,sep = "\t")


