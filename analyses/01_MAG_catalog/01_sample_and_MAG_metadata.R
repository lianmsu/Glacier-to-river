# =============================================================================
# Glacier-to-river -- MAG catalog / Sample and MAG metadata assembly
# =============================================================================
# Split the MAG identifiers into region / site / bin and merge them with the
# sample metadata table.
#
# Manuscript      : Fig. S1
# Main outputs    : in-memory metadata table
# Original script : 20gaoyuan/01MAG_info/1MAG_info.R
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

setwd(file.path(GTR_DATA, "20gaoyuan/01MAG_info"))

mag_list <- read.table("mags_list.txt")
colnames(mag_list) <- "genome"

library(tidyr)

# 按"-"拆分成多列
mag_list_split <- mag_list %>%
  separate(
    genome,
    into = c("region", "site", "bin"),  # 指定新列名
    sep = "-",                             # 分隔符
    remove = FALSE                         # 保留原始列
  )

mag_list_split <- mag_list_split[,-4]

mag_list_split <- mag_list_split %>%
  separate(
    region,
    into = c("ecosystem", "year"),  # 指定新列名
    sep = "_",                             # 分隔符
    remove = FALSE                         # 保留原始列
  )


mag_list_split <- mag_list_split %>%
  separate(
    site,
    into = c("sampling_site", "season_medium"),  # 指定新列名
    sep = "_",                             # 分隔符
    remove = FALSE                         # 保留原始列
  )

factor_all <- read.csv("factor_all.csv")
factor_all <- factor_all[,c(2,12:54)]
sample_infos <- read.delim2("info_corrected.txt", sep = "\t")
sample_infos <- sample_infos %>%
  left_join(factor_all, by=c("SampleID"))


mag_list_split <- merge(mag_list_split, sample_infos, by.x = "site", by.y = "SampleID")

write.table(mag_list_split,"MAG5455_metadata_new.tsv",quote = F,sep = '\t',row.names = F)


sample_infos <- sample_infos %>%
  separate(
    River,
    into = c("ecosystem", "year"),  # 指定新列名
    sep = "_",                             # 分隔符
    remove = FALSE                         # 保留原始列
  )

library(dplyr)

group_counts <- sample_infos %>%
  group_by(Plateau, Type, ecosystem) %>%
  summarise(Count = n(), .groups = "drop")

sum(group_counts[group_counts$Plateau=="Yes" & group_counts$ecosystem=="Yellow",]$Count)


