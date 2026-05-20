
setwd("D:/pku/20gaoyuan/01MAG_info")

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


mag_list_split <- merge(mag_list_split, factor_all, by.x = "site", by.y = "SampleID")

write.table(mag_list_split,"MAG5455_metadata.tsv",quote = F,sep = '\t',row.names = F)


factor_all <- factor_all %>%
  separate(
    River,
    into = c("ecosystem", "year"),  # 指定新列名
    sep = "_",                             # 分隔符
    remove = FALSE                         # 保留原始列
  )

library(dplyr)

group_counts <- factor_all %>%
  group_by(Plateau, Type, ecosystem) %>%
  summarise(Count = n(), .groups = "drop")

sum(group_counts[group_counts$Plateau=="Yes" & group_counts$ecosystem=="Yellow",]$Count)


