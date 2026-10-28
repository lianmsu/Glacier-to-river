# =============================================================================
# Glacier-to-river -- revision nc / VHR and prophage variation
# =============================================================================
# Virus-to-host ratio (VHR) and prophage abundance across regions/mediums;
# regenerates the main Figure 2.
#
# Manuscript      : Revision / Fig. 2
# Main outputs    : Figure2.pdf, 02_VHR_and_prophage_abundance.RData, stat_tests_final.tsv
# Original script : 93_NC_comments/02_VHR_prophage_variation/1_VHR_prophage_abundance.R
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

setwd(file.path(GTR_DATA, "93_NC_comments/02_VHR_prophage_variation"))
library(dplyr)
library(stringr)
library(ggplot2)

#VHR变化####
VHR <- readxl::read_excel("./VHR_from_NC_fig2c.xlsx")

EF <- read.delim2("../01_diversity_comparison/info_corrected.txt",
                  sep = "\t")
VHR <- VHR %>% left_join(EF, by = c("SampleID"="SampleID"))
VHR$Plateau <- 
  ifelse(VHR$Plateau == 'Yes',"Plateau",VHR$Plateau)
VHR$Plateau <- 
  ifelse(VHR$Plateau == 'No',"Plain",VHR$Plateau)
VHR$Plateau <- 
  ifelse(VHR$Plateau == 'glacier',"Glacier",VHR$Plateau)

colnames(VHR)[8] <- 'Region'

VHR$Region <- factor(VHR$Region,
                     levels = c("Glacier","Plateau","Plain"))

VHR <- type.convert(VHR, as.is = TRUE)
VHR$Region <- factor(VHR$Region,
                     levels = c("Glacier","Plateau","Plain"))


library(rstatix)
md_p_val_VHR <- VHR %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(VHR ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val_VHR$y.position <-
  log10(md_p_val_VHR$y.position)



library(ggplot2)
library(ggsignif)#geom_signif
library(ggpubr)#stat_pvalue_manual stat_compare_means
library(gghalves)
library(ggsci)
# library(ggbeeswarm)
plot_VHR_region <- 
  ggplot(VHR, 
         aes(x = Region, y = log10(VHR))) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  # geom_half_point(aes(color=Region),
  #                 side = "r",size = 0.5,
  #                 transformation_params = list(height = 0,
  #                                              width = 0.001,
  #                                              seed = 2))+
  geom_jitter(
    aes(color = Region),
    width = 0.08,
    size = 0.5
  )+
  stat_pvalue_manual(md_p_val_VHR,
                     label = "p.adj.signif",
                     label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="Virus−host abundance ratio (VHR,log10)")+
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
  coord_cartesian(ylim = c(0, 4))
plot_VHR_region

#prophage相对丰度变化####
all_phage_relative_abundance <- read.delim2("all_abundance_tpm.tsv")
all_phage_relative_abundance <- 
  type.convert(all_phage_relative_abundance, as.is = TRUE)


all_phage_lifestyle <- read.delim2("lifecycle_virus.csv",sep = ",")

temperate_pahge <- 
  all_phage_lifestyle[all_phage_lifestyle$lifecycle=="Temperate",]
temperate_phage_tpm <- 
  all_phage_relative_abundance[all_phage_relative_abundance$Contig %in% temperate_pahge$virus_name,]

temperate_phage_tpm <- type.convert(temperate_phage_tpm, as.is = TRUE)

prophage_tpm_sum <- 
  colSums(temperate_phage_tpm[,-1]) %>% 
  as.data.frame() %>%
  setNames("TPM_sum")
prophage_tpm_sum$SampleID <- rownames(prophage_tpm_sum)
prophage_tpm_sum$SampleID <- 
  gsub(".sorted.TPM", "", prophage_tpm_sum$SampleID)  # 可选：进一步简化
prophage_tpm_sum$SampleID <- 
  sub("^[^.]*\\.", "", prophage_tpm_sum$SampleID)

prophage_tpm_sum <- prophage_tpm_sum %>%
  left_join(EF, by = c("SampleID" = "SampleID"))

prophage_tpm_sum$Plateau <- 
  ifelse(prophage_tpm_sum$Plateau == 'Yes',"Plateau",prophage_tpm_sum$Plateau)
prophage_tpm_sum$Plateau <- 
  ifelse(prophage_tpm_sum$Plateau == 'No',"Plain",prophage_tpm_sum$Plateau)
prophage_tpm_sum$Plateau <- 
  ifelse(prophage_tpm_sum$Plateau == 'glacier',"Glacier",prophage_tpm_sum$Plateau)

colnames(prophage_tpm_sum)[8] <- 'Region'

prophage_tpm_sum$Region <- factor(prophage_tpm_sum$Region,
                                  levels = c("Glacier","Plateau","Plain"))


md_p_val_prophage_tpm_sum <- 
  prophage_tpm_sum %>% 
  #  group_by(ecosystem) %>%
  wilcox_test(TPM_sum ~ Region) %>%
  adjust_pvalue(p.col="p",method="bonferroni") %>%
  add_significance(p.col="p.adj") %>% 
  add_xy_position(x="Region",dodge=0.5)
md_p_val_prophage_tpm_sum$y.position <-
  log10(md_p_val_prophage_tpm_sum$y.position)



plot_prophage_tpm_sum_region <- 
  ggplot(prophage_tpm_sum, 
         aes(x = Region, y = TPM_sum)) +
  geom_half_boxplot(aes(fill = Region), color ="black",
                    side = "l",
                    errorbar.draw = T,
                    outlier.shape = NA,
                    width=0.8
  )+
  # geom_half_point(aes(color=Region),
  #                 side = "r",size = 0.5,
  #                 transformation_params = list(height = 0,
  #                                              width = 0.001,
  #                                              seed = 2))+
  geom_jitter(
    aes(color = Region),
    width = 0.08,
    size = 0.5
  )+
  stat_pvalue_manual(md_p_val_prophage_tpm_sum,
                     label = "p.adj.signif",
                     label.size=5,hide.ns = T)+
  # facet_wrap(.~ecosystem,nrow=1)+ #分⾯操作
  # 定义颜⾊
  # scale_fill_npg()+
  # scale_color_npg()+
  labs(x="",y="Relative abundance (TPM) of temperate viruses")+
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
  )
  #scale_y_continuous(expand = expansion(mult = c(0, 0.1)))   +
#  coord_cartesian(ylim = c(0, 4))
plot_prophage_tpm_sum_region

write.table(prophage_tpm_sum, "prophage_tpm_sum.tsv",
  sep = '\t', quote = F, row.names = F)

library(Cairo)
library(cowplot)

group_plot1_VHR_TPM <- 
  plot_grid(plot_VHR_region,plot_prophage_tpm_sum_region,
            # labels = letters[1:2], 
            ncol = 2
            # align = "v",  # 垂直对齐
            # axis = "lr",  # 左和右对齐
            # rel_widths = c(1, 1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
group_plot1_VHR_TPM

ggsave(plot = group_plot1_VHR_TPM,
       '1_group_plot1_VHR_TPM_new.pdf',
       width = 5.52, height = 2.4, units = 'in')





# —— 从RData加载单个子图（非复合图），以便后续统一拼图对齐 ——

# 1. 从4defensefinder_res_density_compare.RData加载DS density的3个单独子图
temp_env <- new.env()
load(file.path(GTR_DATA, "20gaoyuan/10Defensefinder/04_density_and_genome_size.RData"), envir = temp_env)
p_DSdensity_region <- temp_env[["p_DSdensity_glacier_plateau_plain"]]
p_DSdensity_ecosystem <- temp_env[["p_DSdensity_ecossytem"]]
p_DSdensity_medium <- temp_env[["p_DSdensity_medium"]]
rm(temp_env)

# 2. 从5_MAG_DS_pcoa.RData加载DS richness的3个单独子图
temp_env <- new.env()
load(file.path(GTR_DATA, "20gaoyuan/10Defensefinder/05_defense_system_composition_pcoa.RData"), envir = temp_env)
p_DS_shannon <- temp_env[["p_DS_MAG_shannon_region"]]
p_DS_simpson <- temp_env[["p_DS_MAG_simpson_region"]]
p_DS_richness <- temp_env[["p_DS_MAG_richness_region"]]
rm(temp_env)

# 3. 从1_MGE_combination.RData加载MGE图
temp_env <- new.env()
load(file.path(GTR_DATA, "20gaoyuan/12_5455MAG_MGEs/01_MGE_cargo_and_defense_density.RData"), envir = temp_env)
p_MGE_region <- temp_env[["p_glacier_plateau_plain"]]
rm(temp_env)

# —— 用单个3×3 plot_grid统一拼合，替代嵌套拼图，确保面板对齐 ——
# 预先对齐所有子图
plots_aligned <- align_plots(
  p_DSdensity_region, p_DSdensity_ecosystem, p_DSdensity_medium,
  p_DS_shannon,       p_DS_simpson,          p_DS_richness,
  p_MGE_region,       plot_VHR_region,       plot_prophage_tpm_sum_region,
  align = "hv",
  axis = "tblr"
)

group_plot_DS_MGE_VHR_region <- plot_grid(
  plots_aligned[[1]],  plots_aligned[[2]],  plots_aligned[[3]],
  plots_aligned[[4]],  plots_aligned[[5]],  plots_aligned[[6]],
  plots_aligned[[7]],  plots_aligned[[8]],  plots_aligned[[9]],
  ncol = 3,
  rel_heights = c(1, 1, 1),
  rel_widths = c(1, 1, 1)
)
group_plot_DS_MGE_VHR_region

ggsave("Figure2.pdf", 
    group_plot_DS_MGE_VHR_region,
  width = 8, height = 7, units = "in")


save.image("./02_VHR_and_prophage_abundance.RData")
load("./02_VHR_and_prophage_abundance.RData")
