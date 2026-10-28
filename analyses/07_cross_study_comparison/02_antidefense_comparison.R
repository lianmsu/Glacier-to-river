# =============================================================================
# Glacier-to-river -- cross study comparison / Rivers vs aquifers: phage anti-defense
# =============================================================================
# Compare the anti-defense gene content of river phages with groundwater
# phages.
#
# Manuscript      : -
# Main outputs    : 2_river_phage_aquifer_comparision.pdf
# Original script : 20gaoyuan/40plateau_groundwater_comparision/2_antidefense_comparision.R
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
GTR_GW <- normalizePath(Sys.getenv("GTR_GW", unset = file.path(GTR_DATA, "groundwater")), mustWork = FALSE)

setwd(file.path(GTR_DATA, "20gaoyuan/40plateau_groundwater_comparision"))

#河流phage的Anti-defense####

all_59997_phages_DF <- read.delim2("virus_select_all_medium-10kb_high_complete_defense_finder_systems.tsv", 
                                   header = T, check.names = F)
all_59997_phages_DF_anti <- all_59997_phages_DF[all_59997_phages_DF$activity == "Antidefense",]
all_59997_phages_DF_defense <- all_59997_phages_DF[all_59997_phages_DF$activity == "Defense",]

all_59997_phages_DF_anti$contig_id <- str_replace(all_59997_phages_DF_anti$sys_beg, "_\\d+$", "")
all_59997_phages_DF_defense$contig_id <- str_replace(all_59997_phages_DF_defense$sys_beg, "_\\d+$", "")



# all_59997_phages_DF_gene <- read.delim2("virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv", 
#                                         header = T, check.names = F)
# all_59997_phages_DF_gene_anti <- all_59997_phages_DF_gene[all_59997_phages_DF_gene$activity == "Antidefense",]

genomad_res <- read.delim2("virus_select_all_medium-10kb_high_complete_virus_summary_genomad.tsv",
                           header = T, check.names = F)
#对genomad_res 的seq_name列进行清理，删除|及|后的字符，结果作为新列contig

genomad_res$contig <- str_replace(genomad_res$seq_name, "\\|.*", "")

all_55265_phages_DF_anti <- all_59997_phages_DF_anti[all_59997_phages_DF_anti$contig_id %in% genomad_res$contig,]
all_55265_phages_DF_defense <- all_59997_phages_DF_defense[all_59997_phages_DF_defense$contig_id %in% genomad_res$contig,]

#病毒总数太多，含有Anti-defense或者Defense system的contig很少。不能用所有的来绘图

checkv_phage_quality <- read.delim2("quality_summary_checkv.tsv")
checkv_phage_quality <- checkv_phage_quality[checkv_phage_quality$contig_id %in% 
                                               genomad_res$contig,]
checkv_phage_quality <- checkv_phage_quality[checkv_phage_quality$checkv_quality != "Low-quality",]

#去除冰川来源的病毒
library(dplyr)
library(stringr)
# 过滤掉 contig_id 以 "Glacier" 结尾的行
filtered_phage <- checkv_phage_quality %>%
  filter(!str_detect(contig_id, "Glacier$"))  # "$" 表示字符串末尾

all_52795_phage_DF_anti <- 
  all_55265_phages_DF_anti[all_55265_phages_DF_anti$contig_id %in% filtered_phage$contig_id,]
all_52795_phage_DF_defense <- 
  all_55265_phages_DF_defense[all_55265_phages_DF_defense$contig_id %in% filtered_phage$contig_id,]

1804/52795

389/52795

#地下水phage的####
gw_280420_virus_DF <- read.delim2("GWVC.vOTUs_defense_finder_systems.tsv")
gw_280420_virus_DF$contig_id <- str_replace(gw_280420_virus_DF$sys_beg, "_\\d+$", "")

gw_280420_virus_metadata <- read.delim2("GWVC.vOTUs.information.txt")

gw_280420_virus_metadata <- 
  gw_280420_virus_metadata[gw_280420_virus_metadata$checkv_quality %in% c("Complete","High-quality","Medium-quality"),]

gw_27127_virus_DF <- 
  gw_280420_virus_DF[gw_280420_virus_DF$contig_id %in% gw_280420_virus_metadata$contig_id,]

gw_27127_virus_DF_anti <- gw_27127_virus_DF[gw_27127_virus_DF$activity == "Antidefense",]
gw_27127_virus_DF_defense <- gw_27127_virus_DF[gw_27127_virus_DF$activity == "Defense",]

1998/27127
1283/27127


river_aquifer_virus_metadata <- data.frame(
  Anti_defense = c(1804,1998),Defense = c(389,1283),
  phage_num = c(52795,27127),river_aquifer = c("River", "Aquifer")
)


library(dplyr)
library(tidyr)
library(ggplot2)

# 计算比例并转换为长格式
plot_data <- river_aquifer_virus_metadata %>%
  mutate(
    Anti_defense_ratio = Anti_defense / phage_num,
    Defense_ratio = Defense / phage_num
  ) %>%
  select(river_aquifer, Anti_defense_ratio, Defense_ratio) %>%
  pivot_longer(
    cols = c("Anti_defense_ratio", "Defense_ratio"),
    names_to = "System_type",
    values_to = "Ratio"
  ) %>%
  mutate(
    System_type = case_when(
      System_type == "Anti_defense_ratio" ~ "Anti-defense systems",
      System_type == "Defense_ratio" ~ "Defense systems"
    )
  )

plot_data$System_type <- factor(plot_data$System_type, levels = c("Anti-defense systems","Defense systems"))
plot_data$Ratio <- plot_data$Ratio*100

p <- 
ggplot(plot_data, aes(x = river_aquifer, y = Ratio, fill = river_aquifer)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.6) +
  facet_wrap(~System_type, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = c("River" = "#00bbf9", "Aquifer" = "#F28E2B")) +
  geom_text(
    aes(label = sprintf("%.1f%%", Ratio)),  # 显示比例（保留1位小数）
    vjust = 1.5,  # 标签位置（垂直调整）
    size = 4,
    color = "black"
  ) +
  labs(x = "Group", y = "Proportion of phage genomes", 
       title = "Proportion of defense and anti-defense systems") +
  theme_bw() +
  theme(legend.position = "none")
p

setwd(file.path(GTR_DATA, "20gaoyuan/40plateau_groundwater_comparision"))
ggsave("2_river_aquifer_phage_comparision.pdf",p,width = 2.8, height = 5.3, units = 'in')







