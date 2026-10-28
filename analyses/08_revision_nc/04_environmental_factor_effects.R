# =============================================================================
# Glacier-to-river -- revision nc / Environmental drivers of the defensome
# =============================================================================
# Relate DS number/density and anti-defense content to the measured
# environmental factors.
#
# Manuscript      : Revision
# Main outputs    : 04_environmental_factor_effects.RData
# Original script : 93_NC_comments/04_environmental_factors/1_environmental_factor_effects.R
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

setwd(file.path(GTR_DATA, "93_NC_comments/04_environmental_factors"))
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)

EF_info <- read.delim2("info_corrected.txt",sep = '\t')

EF_factor <- read.delim2("factor_all.csv",sep = ",")
EF_factor <- EF_factor[,c(2, 12:54)]

EF_info <- EF_info %>%
  left_join(EF_factor, 
            by = c("SampleID"="SampleID"))
EF_info <- type.convert(EF_info, as.is = TRUE)

EF_info$medium <- ifelse(EF_info$Type == 'Cryoconite',"Glacier",EF_info$Type)
EF_info$medium <- ifelse(EF_info$medium == 'Ice',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'Snow',"Glacier",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'particle',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'water',"Plankton",EF_info$medium)
EF_info$medium <- ifelse(EF_info$medium == 'sediment',"Benthon",EF_info$medium)


EF_info$ecosystem <- EF_info$River
EF_info$ecosystem <- str_remove(EF_info$ecosystem, "_.*")

EF_info$Plateau <- 
  ifelse(EF_info$Plateau == 'Yes',"Plateau",EF_info$Plateau)
EF_info$Plateau <- 
  ifelse(EF_info$Plateau == 'No',"Plain",EF_info$Plateau)
EF_info$Plateau <- 
  ifelse(EF_info$Plateau == 'glacier',"Glacier",EF_info$Plateau)

colnames(EF_info)[7] <- 'Region'

EF_info$Region <- factor(EF_info$Region,
                         levels = c("Glacier","Plateau","Plain"))


# 如果不行，尝试 GBK（常用于 Windows 中文系统导出的文件）
Factor_from_NC <- 
  read.delim2("sourceDataFig1b_fromNC.tsv", 
              sep = '\t', fileEncoding = "GBK")
Factor_from_NC <- Factor_from_NC[,-2]
head(Factor_from_NC)
# 查看是否有空字符串
table(Factor_from_NC$Feature == "")
# 将空字符串替换为 NA，并删除这些行（或保留，但列名无法处理）
Factor_clean <- Factor_from_NC %>%
  mutate(Feature = ifelse(Feature == "", NA, Feature)) %>%
  filter(!is.na(Feature))

Factor_wide <- Factor_clean %>%
  select(-Measurement.units) %>%           # 去掉单位列（若所有特征单位一致可移除）
  pivot_wider(
    names_from = Feature,                  # 新列名来自 Feature 的值
    values_from = Value                    # 新列的值来自 Value
  )

Factor_wide <- type.convert(Factor_wide, as.is = T)

#挑选部分缺失较少的环境因子进行分析
EF_info <- select(EF_info, 
                  SampleID, River, Season, Type, 
                  Distance, Region, Latitude, Longitude,
                  mainstream_or_not, 
                  pH, NO3, DOC_TOC)

EF_info <- EF_info %>%
  left_join(Factor_wide, by = c("SampleID"="SampleID"))


# 不同区域环境因子比较 ####
head(EF_info)

EF_info$Heat[EF_info$Heat < 0] <- NA
EF_info$LE[EF_info$LE < 0] <- NA

EF_info$DOC_TOC_log <- log10(EF_info$DOC_TOC)
EF_info$NH4_log <- log10(EF_info$NH4)
EF_info$NO3_log <- log10(EF_info$NO3)
EF_info$TP_log <- log10(EF_info$TP)


library(tidyverse)
library(rstatix)
library(ggpubr)

# 选定环境因子（数值型）
env_vars1 <- c("Elevation", "Temperature", "DO", "pH",
              "DOC_TOC_log", "NO3_log", "NH4_log", "TP_log",  
              "Heat", "LE", "srad")

# 1. 转为长格式，并逐变量剔除 NA（保留非缺失值）
EF_long <- EF_info %>%
  select(SampleID, Region, all_of(env_vars)) %>%
  pivot_longer(cols = all_of(env_vars), 
               names_to = "Variable", 
               values_to = "Value") %>%
  drop_na(Value)   # 去除该变量中的 NA 行（不影响其他变量）

library(rstatix)
# Shannon, Simpson, Observed（请根据实际列名修改）
diversity_indices <- c("Pielou","Shannon", 
                       "Simpson", "Observed_Richness")  

EF_long_witout_glacier <- 
  EF_long[!EF_long$Region == "Glacier",]

mytheme <-

  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    axis.title = element_text( size = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    #axis.title = element_blank(),
    axis.text = element_text(size = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    #axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
    axis.text.x = element_blank(),
    axis.ticks = element_line(color = "grey20", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = 'none',
    #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
    legend.title = element_text( size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )

plot_elevation <- 
  ggplot(EF_long[EF_long$Variable=="Elevation",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_elevation

plot_temperature <- 
  ggplot(EF_long[EF_long$Variable=="Temperature",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_temperature

plot_DO <- 
  ggplot(EF_long[EF_long$Variable=="DO",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_DO

plot_pH <- 
  ggplot(EF_long[EF_long$Variable=="pH",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_pH

plot_DOC_TOC_log <- 
  ggplot(EF_long[EF_long$Variable=="DOC_TOC_log",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_DOC_TOC_log

plot_NO3_log <- 
  ggplot(EF_long[EF_long$Variable=="NO3_log",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_NO3_log

plot_NH4_log <- 
  ggplot(EF_long[EF_long$Variable=="NH4_log",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_NH4_log

plot_TP_log <- 
  ggplot(EF_long[EF_long$Variable=="TP_log",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_TP_log

plot_Heat <- 
  ggplot(EF_long[EF_long$Variable=="Heat",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_Heat

plot_LE <- 
  ggplot(EF_long[EF_long$Variable=="LE",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_LE

plot_srad <- 
  ggplot(EF_long[EF_long$Variable=="srad",], 
         aes(x = Region, y = Value, fill = Region)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(aes(color = Region), width = 0.2, size = 1.5, alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y", ncol = 4) +
  # scale_fill_brewer(palette = "Set2") +
  # scale_color_brewer(palette = "Set2") +
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_bw() +
  theme(strip.background = element_rect(fill = "grey90"),
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom") +
  labs(x = NULL, y = NULL)+
  mytheme
plot_srad

library(Cairo)
library(cowplot)

group_plot_envs <- 
  plot_grid(plot_elevation,plot_temperature,plot_DO,plot_pH,
            plot_DOC_TOC_log, plot_NO3_log, plot_NH4_log, plot_TP_log,
            plot_Heat, plot_LE, plot_srad,
#            labels = letters[1:5], 
            ncol = 4,
            align = "v",  # 垂直对齐
            axis = "lr",  # 左和右对齐
            rel_widths = c(1, 1,1,1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
group_plot_envs

ggsave(plot = group_plot_envs,
       '1_group_plot_envs.pdf',width = 8.5,height = 5.7,units = 'in')

## 统计 ####
library(tidyverse)
library(rstatix)

# 假设 EF_long 已存在，包含 SampleID, Region, Variable, Value
# 确保 Region 是因子（已设定水平）

# 自定义检验函数（对单个 Variable 子集）
perform_wilcox <- function(data) {
  # 1. 计算每个 Region 组的样本数
  group_counts <- data %>%
    group_by(Region) %>%
    summarise(n = n(), .groups = "drop")
  
  # 2. 保留样本数 >= 2 的组
  valid_groups <- group_counts %>%
    filter(n >= 2) %>%
    pull(Region)
  
  # 3. 如果有效组数 < 2，则无法进行比较，返回空 tibble
  if (length(valid_groups) < 2) {
    return(tibble())
  }
  
  # 4. 过滤数据，仅保留有效组，并移除未使用的因子水平
  data_filtered <- data %>%
    filter(Region %in% valid_groups) %>%
    droplevels()
  
  # 5. 执行 Wilcoxon 检验，校正 p 值，添加位置信息
  res <- tryCatch({
    data_filtered %>%
      wilcox_test(Value ~ Region) %>%
      adjust_pvalue(p.col = "p", method = "bonferroni") %>%
      add_significance(p.col = "p.adj") %>%
      add_xy_position(x = "Region", dodge = 0.5)
  }, error = function(e) {
    # 若出错（例如某组数据仍不足），返回空
    tibble()
  })
  
  return(res)
}

# 对每个 Variable 应用函数，并合并结果
stat_all_envs <- EF_long %>%
  group_by(Variable) %>%
  group_modify(~ perform_wilcox(.x)) %>%
  ungroup()

# 查看成功的统计结果
print(stat_all_envs)


# 1. 过滤掉 Glacier 组（无数据），并去除未使用的因子水平
EF_long_noglacier <- EF_long %>%
  filter(Region != "Glacier") %>%        # 剔除空组
  mutate(Region = droplevels(Region))    # 删除因子中未使用的水平 "Glacier"

# 2. 再次检查样本量（现在应该只有 Plateau 和 Plain）
table(EF_long_noglacier$Region)

# 3. 重新执行检验（现在只有两组，会直接进行两组间比较）
statistic_envs_NO3 <- 
  EF_long_noglacier %>%
  filter(Variable == "NO3_log", Value != -Inf) %>%   # 正确筛选行
  wilcox_test(Value ~ Region) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Region", dodge = 0.5) %>%
  ungroup()

statistic_envs_NH4 <- 
  EF_long_noglacier %>%
  filter(Variable == "NH4_log", Value != -Inf) %>%   # 正确筛选行
  wilcox_test(Value ~ Region) %>%
  adjust_pvalue(p.col = "p", method = "bonferroni") %>%
  add_significance(p.col = "p.adj") %>%
  add_xy_position(x = "Region", dodge = 0.5) %>%
  ungroup()


# 环境因子影响的相关分析及mantel分析 #####################################
################################ 环境因子数据
EF_info$River <- str_remove(EF_info$River, "_.*")

################################ 原核丰度数据
MAGs_4903_relative_abundance <- 
  read.delim2("../01_diversity_comparison/丰度/Bins/all_abundance_rpkm_bins.tsv")
# 将所有非Genome列转换为数值型
MAGs_4903_relative_abundance <- MAGs_4903_relative_abundance %>%
  mutate(across(-Genome, as.numeric))
# 1. 将第一列 'Genome' 设为行名，并移除该列
MAGs_4903_relative_abundance_2 <- MAGs_4903_relative_abundance
rownames(MAGs_4903_relative_abundance_2) <- MAGs_4903_relative_abundance_2$Genome
MAGs_4903_relative_abundance_2 <- MAGs_4903_relative_abundance_2[, -1]  # 删除 Genome 列
# 2. 转置（行名变列名，列名变行名）
t_MAGs_4903_relative_abundance_2 <- as.data.frame(t(MAGs_4903_relative_abundance_2))
# 3. 将行名（原样本名）保存为新列 'Sample'
t_MAGs_4903_relative_abundance_2$Sample <- rownames(t_MAGs_4903_relative_abundance_2)
# 4. 重新排列列，将 Sample 置于第一列
t_MAGs_4903_relative_abundance_2 <- t_MAGs_4903_relative_abundance_2[, c("Sample", setdiff(names(t_MAGs_4903_relative_abundance_2), "Sample"))]
# 5. 清除行名（以免冗余）
rownames(t_MAGs_4903_relative_abundance_2) <- NULL
# 处理样本名
t_MAGs_4903_relative_abundance_2$Sample <- 
  gsub(".sorted.RPKM", "", t_MAGs_4903_relative_abundance_2$Sample)  # 可选：进一步简化
# sample_names <- gsub("_", " ", sample_names)      # 可选：用空格替代下划线
t_MAGs_4903_relative_abundance_2$Sample <- 
  sub("^[^.]*\\.", "", t_MAGs_4903_relative_abundance_2$Sample)

#################病毒丰度数据#
# 读取
phage_relative_abundance <- 
  read.delim2("../01_diversity_comparison/丰度/Virus/all_abundance_rpkm.tsv")
# 将所有非Genome列转换为数值型
phage_relative_abundance <- phage_relative_abundance %>%
  mutate(across(-Contig, as.numeric))
# 1. 将第一列 'Contig' 设为行名，并移除该列
phage_relative_abundance_2 <- phage_relative_abundance
rownames(phage_relative_abundance_2) <- phage_relative_abundance_2$Contig
phage_relative_abundance_2 <- phage_relative_abundance_2[, -1]  # 删除 Contig 列
# 2. 转置（行名变列名，列名变行名）
t_phage_relative_abundance_2 <- as.data.frame(t(phage_relative_abundance_2))
# 3. 将行名（原样本名）保存为新列 'Sample'
t_phage_relative_abundance_2$Sample <- rownames(t_phage_relative_abundance_2)
# 4. 重新排列列，将 Sample 置于第一列
t_phage_relative_abundance_2 <- 
  t_phage_relative_abundance_2[, c("Sample", setdiff(names(t_phage_relative_abundance_2), "Sample"))]
# 5. 清除行名（以免冗余）
rownames(t_phage_relative_abundance_2) <- NULL
# 处理样本名
t_phage_relative_abundance_2$Sample <- 
  gsub(".sorted.RPKM", "", t_phage_relative_abundance_2$Sample)  # 可选：进一步简化
# sample_names <- gsub("_", " ", sample_names)      # 可选：用空格替代下划线
t_phage_relative_abundance_2$Sample <- 
  sub("^[^.]*\\.", "", t_phage_relative_abundance_2$Sample)

########################统计分析

# 安装 linkET（若未安装）
# devtools::install_github("Hy4m/linkET")
library(vegan)
library(cluster)   # 用于 Gower 距离（处理分类变量）
library(tidyverse)
library(linkET)

# 假设 EF_info 已存在，包含 SampleID 和所有变量
env_df <- EF_info %>%
  select(SampleID, pH, NO3, DOC_TOC, Temperature, DO, Elevation, 
         Heat, LE, NH4, srad, TP, River, Season, Type, Region)

# # 确定数值型和分类型列名
num_vars <- c("pH", "NO3", "DOC_TOC", "Temperature", "DO", "Elevation",
              "Heat", "LE", "NH4", "srad", "TP")
cat_vars <- c("River", "Season", "Type", "Region")
# 
# # 为了计算 Pearson 相关（热图），我们将所有变量转换为数值型：
# # 数值型保持不变，分类变量转为因子水平编码（或虚拟变量，但水平编码更简单）
# env_for_cor <- env_df %>%
#   mutate(across(all_of(cat_vars), ~ as.numeric(factor(.)))) %>%
#   select(-SampleID)   # 移除样本列
# 
# # 相关矩阵计算（pairwise 处理缺失）
# cor_matrix <- cor(env_for_cor, use = "pairwise.complete.obs")

#######2. 群落距离矩阵（MAG 和 Phage）
# 假设已有 mag_dist 和 phage_dist（dist 对象）
# 如果还没有，请先计算：
# 去除全为零的行（若有）
mag_abund <- t_MAGs_4903_relative_abundance_2
rownames(mag_abund) <- mag_abund$Sample
mag_abund <- mag_abund[, -1]
# 获取环境数据中的样本名（用于后续筛选）
#samples_env <- env_df$SampleID
#common <- intersect(rownames(mag_abund), samples_env)
#mag_sub <- mag_abund[common, ]
# mag_sub <- mag_abund[rowSums(mag_abund) > 0, ]
# mag_dist <- vegdist(mag_sub, method = "bray")

phage_abund <- t_phage_relative_abundance_2
rownames(phage_abund) <- phage_abund$Sample
phage_abund <- phage_abund[, -1]
#common <- intersect(rownames(mag_abund), env_for_cor)
#mag_sub <- mag_abund[common, ]
# phage_sub <- phage_abund[rowSums(phage_abund) > 0, ]
# phage_dist <- vegdist(phage_sub, method = "bray")

# 获取环境数据中的样本名（用于后续筛选）
# samples_env <- env_df$SampleID

# # 从距离矩阵中提取与环境样本相同的子集（如果距离矩阵的标签与环境样本不完全一致）
# # 注意：距离矩阵的标签是样本名，需要与环境样本名匹配
# common_samples <- intersect(attr(mag_dist, "Labels"), samples_env)
# mag_dist_sub <- as.dist(as.matrix(mag_dist)[common_samples, common_samples])
# phage_dist_sub <- as.dist(as.matrix(phage_dist)[common_samples, common_samples])

#########################3. 计算每个环境变量与群落的 Mantel 相关
# 计算单变量与群落距离的 Mantel 检验
# get_mantel <- function(var_name, env_df, comm_dist, num_vars, cat_vars) {
#   # var_name: 变量名
#   # env_df: 包含 SampleID 和所有变量的数据框
#   # comm_dist: 群落距离矩阵（dist 对象）
#   
#   # 提取该变量的值
#   vec <- env_df[[var_name]]
#   
#   # 获取非缺失样本
#   idx <- which(!is.na(vec))
#   if (length(idx) < 3) return(NA)  # 样本太少
#   
#   # 获取样本名
#   samps <- env_df$SampleID[idx]
#   # 与群落距离的样本取交集
#   common <- intersect(samps, attr(comm_dist, "Labels"))
#   if (length(common) < 3) return(NA)
#   
#   # 子集该变量的值（按 common 顺序）
#   vec_sub <- vec[env_df$SampleID %in% common]
#   names(vec_sub) <- common
#   # 按 common 顺序重新排序
#   vec_sub <- vec_sub[common]
#   
#   # 计算该变量的距离矩阵
#   if (var_name %in% num_vars) {
#     # 数值型：标准化后欧氏距离
#     vec_scaled <- scale(vec_sub)
#     env_dist <- dist(vec_scaled)
#   } else {
#     # 分类变量：Gower 距离
#     df_temp <- data.frame(factor = as.factor(vec_sub))
#     env_dist <- daisy(df_temp, metric = "gower")
#   }
#   
#   # 从群落距离中提取对应样本的子距离矩阵
#   comm_sub <- as.dist(as.matrix(comm_dist)[common, common])
#   
#   # Mantel 检验
#   test <- mantel(comm_sub, env_dist, method = "pearson", permutations = 999)
#   
#   return(data.frame(
#     from = var_name,
#     r = test$statistic,
#     p = test$signif
#   ))
# }

# ######### 批量计算
# all_vars <- c(num_vars, cat_vars)
# mantel_res_list_MAG <- list()
# 
# for (v in all_vars) {
#   res <- get_mantel(v, env_df, mag_dist_sub, num_vars, cat_vars)
#   if (!is.na(res[1])) {
#     mantel_res_list_MAG[[v]] <- res
#   }
# }
# 
# # 合并成数据框（每个变量一行）
# mantel_MAG_df <- bind_rows(mantel_res_list_MAG) %>%
#   mutate(
#     to = "MAG",   # 先处理 MAG，后面再添加 Phage
#     sig = ifelse(p < 0.001, "***",
#                  ifelse(p < 0.01, "**",
#                         ifelse(p < 0.05, "*", "ns")))
#   )
# 
# # 对 Phage 重复，得到另一组结果（改为 to = "Phage"）
# # 可以合并两个数据集，添加列 to 区分
# mantel_res_list_phage <- list()
# 
# for (v in all_vars) {
#   res <- get_mantel(v, env_df, phage_dist_sub, num_vars, cat_vars)
#   if (!is.na(res[1])) {
#     mantel_res_list_phage[[v]] <- res
#   }
# }
# 
# # 合并成数据框（每个变量一行）
# mantel_phage_df <- bind_rows(mantel_res_list_phage) %>%
#   mutate(
#     to = "Phage",   #  Phage
#     sig = ifelse(p < 0.001, "***",
#                  ifelse(p < 0.01, "**",
#                         ifelse(p < 0.05, "*", "ns")))
#   )
# 
# 
# # 假设我们已得到 mantel_MAG 和 mantel_Phage 两个数据框
# # 合并为长数据
# mantel_combined <- bind_rows(
#   mantel_MAG_df %>% mutate(Community = "MAG"),
#   mantel_phage_df %>% mutate(Community = "Phage")
# ) %>%
#   rename(from2 = from, to2 = Community, rd = r, pd = p) %>%
#   mutate(
#     # 线型：实线（p < 0.05）或虚线（p >= 0.05）
#     linetype = ifelse(pd < 0.05, "solid", "dashed"),
#     # 颜色：正相关（红色）或负相关（蓝色）
#     color = ifelse(rd > 0, "positive", "negative"),
#     # 设定线宽（可按比例缩放）
#     size = abs(rd) * 3   # 乘数可调
#   )
# 
# mantel_combined <- mantel_combined %>%
#   rename(from = from2)
# mantel_combined <- mantel_combined %>%
#   select(-to2)
# #########绘图
# 
# # 1. 环境因子相关矩阵（cor_matrix）
# 
# correlation <- correlate(cor_matrix)   # 生成 linkET 的 correlate 对象
# 
# # 2. 绘制热图
# plot1_EF_community_cor <- 
#   qcorrplot(correlation, type = "upper", diag = FALSE) +
#   # 添加 Mantel 连线
#   geom_couple(
#     data = mantel_combined,
#     aes(colour = color, size = size, linetype = linetype),
#     curvature = 0.25
# #    curvature = nice_curvature()
#   ) +
#   # 手动设置颜色
#   scale_colour_manual(values = c(positive = "#E41A1C", 
#                                  negative = "#377EB8"),
#                       name = "Mantel r") +
#   scale_size_continuous(range = c(0.5, 3), guide = "none") +  # 线宽
#   scale_linetype_identity() +   # 直接使用 linetype 列
#   theme(
#     legend.position = "right",
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
# 
# print(plot1_EF_community_cor)


## Mantel 分析###############################################
env_use <- EF_info %>%
  select(
    SampleID,
    pH, NO3, DOC_TOC,
    Temperature,
    DO,
    Elevation,
    Heat,
    LE,
    NH4,
    srad,
    TP,
    River,
    Season,
    Type,
    Region
  ) %>%
  column_to_rownames("SampleID")

mag_abund <- mag_abund[rowSums(mag_abund) > 0, ]

common <- Reduce(intersect,
                 list(rownames(env_use),
                      rownames(mag_abund),
                      rownames(phage_abund)))

env_use   <- env_use[common, ]
mag_use   <- mag_abund[common, ]
phage_use <- phage_abund[common, ]

spec_use <- cbind(
  mag_use,
  phage_use
)

# spec_select <- list(
#   MAG   = 1:ncol(mag_use),
#   Phage = (ncol(mag_use)+1):ncol(spec_use)
# )

# env_select <- as.list(
#   setNames(
#     seq_len(ncol(env_use)),
#     colnames(env_use)
#   )
# )
# install.packages("FD")
library(FD)

mantel_tbl <- mantel_test(
  spec = spec_use,
  env = env_use,
  
  spec_select = list(
    MAG = 1:ncol(mag_use),
    Phage = (ncol(mag_use)+1):ncol(spec_use)
  ),
  
  env_select = as.list(
    setNames(
      seq_len(ncol(env_use)),
      colnames(env_use)
    )
  ),
  
  spec_dist = "bray",
  env_dist = "gower",
  
  na_omit = TRUE,
  
  permutations = 999
)


# mantel_tbl <- mantel_test(
#   spec = spec_use,
#   env = env_use,
#   
#   spec_select = spec_select,
#   env_select = env_select,
#   
#   spec_dist = "bray",
#   env_dist = "gower",
#   
#   permutations = 999
# )

mantel_plot <- mantel_tbl %>%
  rename(
    from = env,
    to   = spec,
    rd   = r,
    pd   = p
  ) %>%
  mutate(
    colour = ifelse(rd > 0,
                    "positive",
                    "negative"),
    
    linetype = ifelse(pd < 0.05,
                      "solid",
                      "dashed")
  ) %>%
  mutate(
    colour2 = case_when(
      abs(rd) < 0.3 | pd > 0.05 ~ "weak",
      rd > 0 ~ "positive",
      rd < 0 ~ "negative"
    )
  )

env_cor <- env_use %>%
  mutate(
    across(
      c(River, Season, Type, Region),
      ~ as.numeric(factor(.))
    )
  )
correlation <- correlate(env_cor,     method = "pearson",
                         use = "pairwise.complete.obs")

plot_mantel_EF_community <- 
qcorrplot(correlation,
          grid_col = "#00468BFF",
          "white","#42B540FF",
          grid_size = 0.2,
          type = "upper",
          diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(colours = c("#00468BFF",
                                   "white","#42B540FF"),
                       limits = c(-1, 1))+
  geom_mark(size = 4,
            only_mark = T,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05,
            colour = 'white')+
  
  geom_couple(
    data = mantel_plot,
    aes(
      colour = colour2,
      size = abs(rd),
      linetype = linetype
    ),
    curvature = nice_curvature()
  ) +
  
  scale_colour_manual(
    values = c(
      positive = "#ff9f1c",
      negative = "#377EB8",
      weak = "grey70"
    )
  ) +
  
  scale_size(range = c(0.5,1)) +
  scale_linetype_identity() 

plot_mantel_EF_community

write.table(as.data.frame(correlation$r), "correlation_community.tsv", sep = "\t", quote = F, row.names = TRUE)
write.table(mantel_plot, "mantel_plot_community.tsv", sep = "\t", quote = F, row.names = F)


save.image("./04_environmental_factor_effects.RData")
load("./04_environmental_factor_effects.RData")

## 环境因子与DS####
#把每个MAG当作样本#
MAG_DS_type_info <- read.delim2("8_MAG_DS_type_num.tsv")
rownames(MAG_DS_type_info) <- MAG_DS_type_info$MAG
MAG_DS_type_info <- select(MAG_DS_type_info, -MAG)

MAG_meta_info <- read.delim2("MAG5455_metadata_new.tsv")
MAG_meta_info_2 <- 
  MAG_meta_info[MAG_meta_info$genome %in% MAG_DS_type_info$MAG,]
MAG_meta_info_2 <- MAG_meta_info_2[,c(1,2)]
MAG_meta_info_2 <- MAG_meta_info_2 %>%
  left_join(EF_info, by = c("site"="SampleID"))
MAG_meta_info_2 <- MAG_meta_info_2 %>%
  select(
    genome,
    pH, NO3, DOC_TOC,
    Temperature,
    DO,
    Elevation,
    Heat,
    LE,
    NH4,
    srad,
    TP,
    River,
    Season,
    Type,
    Region
  ) %>%
  column_to_rownames("genome")

# Mantel #
library(FD)
MAG_DS_mantel_tbl <- mantel_test(
  spec = MAG_DS_type_info,
  env = MAG_meta_info_2,
  
  spec_select = list(
    DS = 1:ncol(MAG_DS_type_info)
  ),
  
  env_select = as.list(
    setNames(
      seq_len(ncol(MAG_meta_info_2)),
      colnames(MAG_meta_info_2)
    )
  ),
  
  spec_dist = "bray",
  env_dist = "gower",
  
  na_omit = TRUE,
  
  permutations = 999
)


# mantel_tbl <- mantel_test(
#   spec = spec_use,
#   env = env_use,
#   
#   spec_select = spec_select,
#   env_select = env_select,
#   
#   spec_dist = "bray",
#   env_dist = "gower",
#   
#   permutations = 999
# )

MAG_DS_mantel_plot <- MAG_DS_mantel_tbl %>%
  rename(
    from = env,
    to   = spec,
    rd   = r,
    pd   = p
  ) %>%
  mutate(
    colour = ifelse(rd > 0,
                    "positive",
                    "negative"),
    
    linetype = ifelse(pd < 0.05,
                      "solid",
                      "dashed")
  ) %>%
  mutate(
    colour2 = case_when(
      abs(rd) < 0.3 | pd > 0.05 ~ "weak",
      rd > 0 ~ "positive",
      rd < 0 ~ "negative"
    )
  )

MAG_meta_info_2_cor <- MAG_meta_info_2 %>%
  mutate(
    across(
      c(River, Season, Type, Region),
      ~ as.numeric(factor(.))
    )
  )
MAG_meta_info_2_correlation <- 
  correlate(MAG_meta_info_2_cor,
            method = "pearson",
            use = "pairwise.complete.obs")

plot_mantel_MAG_DS <- 
  qcorrplot(MAG_meta_info_2_correlation,
            grid_col = "#00468BFF",
            "white","#42B540FF",
            grid_size = 0.2,
            type = "upper",
            diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(colours = c("#00468BFF",
                                   "white","#42B540FF"),
                       limits = c(-1, 1))+
  geom_mark(size = 4,
            only_mark = T,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05,
            colour = 'white')+
  
  geom_couple(
    data = MAG_DS_mantel_plot,
    aes(
      colour = colour2,
      size = abs(rd),
      linetype = linetype
    ),
    curvature = nice_curvature()
  ) +
  
  scale_colour_manual(
    values = c(
      positive = "#ff9f1c",
      negative = "#377EB8",
      weak = "grey70"
    )
  ) +
  
  scale_size(range = c(0.5,1)) +
  scale_linetype_identity()


plot_mantel_MAG_DS

library(Cairo)
library(cowplot)

group_plot_mantel <- 
  plot_grid(plot_mantel_EF_community,
            plot_mantel_MAG_DS,
            labels = letters[1:2], ncol = 2,
            align = "hv"  # 垂直对齐
          #  axis = "lr",  # 左和右对齐
          #  rel_widths = c(1, 1)
            # label_x = 0.09,            # 标签x坐标（1为最右侧）
            #  label_y = 0.9             # 标签y坐标（1为最上方）
  )
group_plot_mantel

ggsave(plot = group_plot_mantel,
       '1_group_plot_mantel.pdf',
       width = 12.4,height = 6.2,units = 'in')



# CCA分析 环境因子对微生物及病毒的影响 ####
#CCA分析需要样本的环境因子不含NA


# 结构方程模型 ####

# 这个方法也要求数据不能有NA #




save.image("./04_environmental_factor_effects.RData")
load("./04_environmental_factor_effects.RData")

write.table(EF_info, "EF_info.tsv",
            quote = F, row.names = F,sep = "\t")
