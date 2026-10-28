# =============================================================================
# Glacier-to-river -- defensome landscape / Taxonomic composition of MAGs
# =============================================================================
# Aggregate the MAG taxonomy into phylum-level composition across regions
# and ecosystems.
#
# Manuscript      : Fig. S7
# Main outputs    : 7_tax_legend.pdf
# Original script : 20gaoyuan/10Defensefinder/7_tax_composition.R
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

library(dplyr)
setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
mag_baclineage <- read.delim2("gtdbtk.bac120.summary.tsv")
mag_arclineage <- read.delim2("gtdbtk.ar53.summary.tsv")
mag_lineage <- rbind(mag_baclineage, mag_arclineage) %>% select(user_genome, classification)
library(tidyr)
# 拆分 classification 列
mag_lineage_split <- mag_lineage %>%
  separate(
    classification,
    into = c("domain", "phylum", "class", "order", "family", "genus", "species"),
    sep = ";",
    remove = FALSE  # 如果保留原始列，设为 FALSE；否则设为 TRUE
  )
# 查看拆分后的数据
head(mag_lineage_split)
colnames(mag_lineage_split)[1] <- 'genome'

# 使用 sub 函数替换掉第三个 _ 及其后面的内容
mag_lineage_split$phylum <- sub("^(([^_]*_){2}[^_]*)_.*$", "\\1", 
                                mag_lineage_split$phylum)
mag_lineage_split$phylum <- gsub("^p__", "", mag_lineage_split$phylum)
phyla_count <- table(mag_lineage_split$phylum) %>% as.data.frame()

MAG_metadata <- read.csv("../01MAG_info/5455MAG_metadata.csv")
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

table(MAG_metadata$ecosystem)
table(MAG_metadata$Region)
table(MAG_metadata$Type)
table(MAG_metadata$medium)

mag_lineage_geography <- merge(mag_lineage_split, MAG_metadata, by = "genome")
mag_lineage_geography <- select(mag_lineage_geography, genome, 
                                phylum, ecosystem, Region, Type, medium)

phyla_count_merged <- phyla_count %>%
  mutate(Var1 = ifelse(Freq < 10, "Others", as.character(Var1))) %>%
  group_by(Var1) %>%
  summarise(Freq = sum(Freq)) %>%
  arrange(desc(Freq))  # 按频率降序排列
phyla_count_merged$color <- c(
  "#00a8e8", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD",
  "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF",
  "#AEC7E8", "#FFBB78", "#CCCCCC", "#FF9896", "#C5B0D5",
  "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5", "#98DF8A"
)
# 确保Others是灰色
phyla_count_merged$color[phyla_count_merged$Var1 == "Others"] <- "#CCCCCC"

mag_lineage_geography$phyla <- ifelse(mag_lineage_geography$phylum %in% phyla_count_merged$Var1,
                                      mag_lineage_geography$phylum, "Others")

mag_lineage_geography <- merge(mag_lineage_geography, phyla_count_merged,
                               by.x = "phyla", by.y = "Var1")


library(ggplot2)

mag_lineage_geography$ecosystem <- factor(mag_lineage_geography$ecosystem, 
                                       levels = c("Glacier","CJ","HH","LCJ","YJ"))

p_ecosystem <- 
ggplot(mag_lineage_geography, aes(x = ecosystem, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = mag_lineage_geography$phyla,
                      breaks = mag_lineage_geography$color) +
  labs(x = "Ecosystem", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_ecosystem

mag_lineage_geography$Region <- factor(mag_lineage_geography$Region, 
                                       levels = c("Glacier","Plateau","Plain"))
p_Region <- 
  ggplot(mag_lineage_geography, aes(x = Region, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = mag_lineage_geography$phyla,
                      breaks = mag_lineage_geography$color) +
  labs(x = "Region", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_Region


mag_lineage_geography$medium <- factor(mag_lineage_geography$medium, 
                                       levels = c("Glacier","Benthon","Plankton"))
p_medium <- 
  ggplot(mag_lineage_geography, aes(x = medium, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = mag_lineage_geography$phyla,
                      breaks = mag_lineage_geography$color) +
  labs(x = "Medium", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_medium

mag_lineage_geography$Type <- factor(mag_lineage_geography$Type, 
                                       levels = c("Snow","Ice","Cryoconite",
                                                  "sediment","particle","water"))
p_Type <- 
  ggplot(mag_lineage_geography, aes(x = Type, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = mag_lineage_geography$phyla,
                      breaks = mag_lineage_geography$color) +
  labs(x = "Detailed medium", 
       y = "Relative proportion of MAGs",
       fill = "Phylum") +
  theme_bw() +
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
        plot.background = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(linewidth = 0.5),
        panel.grid = element_blank(),
        # axis.line.x = element_line(linewidth = 0.5),
        # axis.line.y = element_line(linewidth = 0.5),
        axis.title.x = element_text(size = 12, face = "plain"),
        axis.title.y = element_text(size = 12, face = "plain"),
        axis.ticks.x = element_line(linewidth = 0.5),
        axis.ticks.y = element_line(linewidth = 0.5),
        axis.ticks.length.x = unit(0.1, "cm"),
        axis.ticks.length.y = unit(0.1, "cm"),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        legend.position = "none")
p_Type

write.table(mag_lineage_geography,"mag_lineage_geography.tsv",sep = "\t",
  quote = F, row.names = F)

library(Cairo)
library(cowplot)

row2 <- plot_grid(p_Region, 
                  p_medium, 
                  nrow = 1,ncol = 2,
                  labels = letters[2:3]
                  # rel_heights = c(1, 1, 1.5)
)
row2
final_plot <- plot_grid(p_ecosystem, 
                         row2, 
                         p_Type, 
                         nrow = 3,ncol = 1,
                         labels = letters[1:3]
                         # rel_heights = c(1, 1, 1.5)
)
final_plot

setwd(file.path(GTR_DATA, "20gaoyuan/10Defensefinder"))
cairo_pdf(width = 5.6,height = 7,filename = "7_tax_composition.pdf")
final_plot
dev.off()

p_legend <-   
  ggplot(mag_lineage_geography, aes(x = Type, fill = color)) +
  geom_bar(position = "fill", width = 0.7) +
  scale_fill_identity(guide = "legend", 
                      labels = mag_lineage_geography$phyla,
                      breaks = mag_lineage_geography$color)
p_legend

ggplot2::ggsave("7_tax_legend.pdf",p_legend, 
                width = 5.6, height = 7, units = "in")







