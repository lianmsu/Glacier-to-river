#基本数据读取####
setwd("D:/pku/20gaoyuan/10Defensefinder")
library(dplyr)
library(tidyr)
library(stringr)

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

MAG_metadata <- read.csv("D:/pku/20gaoyuan/01MAG_info/5455MAG_metadata.csv")
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

DF_mechanism <- read.delim2('df_List_Systems_lpw.csv',sep = ',')
DF_mechanism$Effector[DF_mechanism$Effector == ""] <- "Unknown"
DF_mechanism$Mechnism[DF_mechanism$Mechnism == ""] <- "Unknown"

#每个点代表一种DS####
##Bray NMDS###############

library(tidyverse)
library(vegan)
library(ape)
#library(phyloseq)
library(geomtextpath)
library(patchwork)


# 标准化数据，使其适用于 Bray-Curtis 距离计算
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


count_matrix_standard = decostand(count_matrix_standard,method = "hellinger")  
# NMDS（非度量多维标度分析），k=2 表示降维到 2 维
mds_bray = metaMDS(t(count_matrix_standard),
                   distance = "bray", k = 2,trymax = 100)  


# 提取 NMDS 分析的 stress 值，四舍五入保留 4 位小数
# 计算 stress 值，越小表示 NMDS 结果越可靠
stress1 <- round(mds_bray$stress,4) 

# 提取 NMDS 计算出的坐标，并与样本分组信息合并
data.plot <- mds_bray$points %>% as.data.frame() %>%  # 转换为数据框
  rownames_to_column(var="title") %>%  # 添加样本 ID 列
  left_join(., DF_mechanism, by="title")  # 关联样本分组信息
data.plot[is.na(data.plot)] <- "Unknown"


# 绘制 Bray-Curtis NMDS 结果
p1 <- ggplot(data.plot, aes(MDS1, MDS2))+
  # 按 Mechnism 颜色区分样本点
  geom_point( aes(colour = Mechnism),size = 2, show.legend = F)+  
  # 95% 置信椭圆
  stat_ellipse(aes(label = Mechnism, color=Mechnism),
               geom = "textpath",
               hjust = 0.3, level = 0.95, linewidth=0.5,
               linetype=2, show.legend = F) +
  xlab("MDS1") + ylab("MDS2") +  # 轴标签
  # 添加 stress 值到标题
  labs(title = paste("bray stress =", stress1, sep=" ")) +  
  theme_test() +
  theme(plot.title = element_text(vjust=0.5, hjust = 0.5))  
p1


##PCoA######
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]

library(vegan)
library(ape)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

# 计算 Bray-Curtis 距离并进行 PCoA 分析
pcoa <- vegdist(t(count_matrix_standard), method = "bray") %>% 
    # 进行主坐标分析（PCoA），不进行负特征值校正
  pcoa(correction = "none")

summary(pcoa)

# 提取 PCoA 结果并合并分组信息
pcoadata <- data.frame(pcoa$vectors[,1], pcoa$vectors[,2]) %>%  
  # 获取 PC1 和 PC2 坐标
  set_colnames(c("PC1", "PC2")) %>%  # 重命名列名
  rownames_to_column(var = "title") %>%  # 将行名转换为列 sample
  # 关联 Subtype 信息
  left_join(., DF_mechanism %>% 
              dplyr::select(title, Mechnism),by="title")

pcoadata[is.na(pcoadata)] <- "Unknown"
# 颜色设置
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(6)

#绘制pcoa散点图
plot <- ggplot(pcoadata, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Mechnism, fill = Mechnism), size = 4) +  
  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot


#每个点代表一个MAGs###########
##PCoA######
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]

library(vegan)
library(ape)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

# 计算 Bray-Curtis 距离并进行 PCoA 分析
pcoa <- vegdist(count_matrix_standard, method = "bray") %>% 
  # 进行主坐标分析（PCoA），不进行负特征值校正
  pcoa(correction = "none")

summary(pcoa)

# 提取 PCoA 结果并合并分组信息
pcoadata <- data.frame(pcoa$vectors[,1], pcoa$vectors[,2]) %>%  
  # 获取 PC1 和 PC2 坐标
  set_colnames(c("PC1", "PC2")) %>%  # 重命名列名
  rownames_to_column(var = "genome") %>%  # 将行名转换为列 sample
  # 关联 Subtype 信息
  left_join(., MAG_metadata %>% 
              dplyr::select(genome, Region),by="genome")

pcoadata[is.na(pcoadata)] <- "Unknown"

#write.csv(pcoadata,"pcoadata.csv")

pcoadata2 <- read.csv("pcoadata.csv")

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

pcoadata2 <- read.csv("pcoadata.csv")
pcoadata2 <- left_join(pcoadata2,mag_lineage_split,suffix = c("genome","genome"))
# 统计 phylum 列中每个门的计数
p_counts <- pcoadata2 %>%
  count(phylum)
# 找出计数小于 1000 的门
low_count_p <- p_counts %>%
  filter(n < 50) %>%
  pull(phylum)
# 将计数小于 1000 的门及其后续分类等级设置为 "others"
pcoadata2 <- pcoadata2 %>%
  mutate(across(phylum:species, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                         ifelse(phylum %in% low_count_p, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))
phylum <- table(pcoadata2$phylum) %>% as.data.frame()

#添加每个MAG的size大小
pcoadata3 <- pcoadata2 %>% left_join(., MAG_metadata, by="genome")
colnames(MAG_quality)[1] <- "genome"
pcoadata3 <- pcoadata3 %>% left_join(., MAG_quality, by="genome")
pcoadata3$MAG_size <- pcoadata3$Genome_Size/1000000
pcoadata3$magsize <- floor(pcoadata3$MAG_size)
pcoadata3$magsize <- factor(pcoadata3$magsize)

###1、根据ecosystem绘图####

# 颜色设置
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(6)

library(ggplot2)
#绘制pcoa散点图
plot_MAG_ecosystem <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = ecosystem, fill = ecosystem), size = 0.5) +  
#  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
       #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
       legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_ecosystem

###2、region着色####


#绘制pcoa散点图
plot_MAG_region <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Region.x, fill = Region.x), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_region

###3、medium着色####
#绘制pcoa散点图
plot_MAG_medium <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = medium, fill = medium), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_medium

###4、物种着色####
phylum <- table(pcoadata2$phylum) %>% as.data.frame()
color_palette_phylum <- colorRampPalette(brewer.pal(12, "Paired"))(10)


#绘制pcoa散点图
plot_MAG_phylum <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = phylum, fill = phylum), size = 0.5) +  
  scale_color_manual(values = color_palette_phylum) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_phylum

###5、用MAG size着色####


color_palette_magsize <- colorRampPalette(brewer.pal(12, "Paired"))(12)
#绘制pcoa散点图
plot_MAG_DS_magsize <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = magsize, fill = magsize), size = 0.5) +  
  scale_color_manual(values = color_palette_magsize) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_magsize


###6、用season着色####


color_palette_season <- colorRampPalette(brewer.pal(12, "Paired"))(4)
#绘制pcoa散点图
plot_MAG_DS_season <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Season, fill = Season), size = 0.5) +  
  scale_color_manual(values = color_palette_season) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_season

###7、用domain着色####

pcoadata3$domain <- factor(pcoadata3$domain, levels = c("d__Archaea","d__Bacteria"))
#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_domain <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = domain, fill = domain), size = 0.5) +  
#  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
                                "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_domain

###8、用Elevation着色####

#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_Elevation <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Elevation, fill = Elevation), size = 0.5) +  
  #  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
  #                               "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
        color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_Elevation


library(Cairo)
library(cowplot)
cairo_pdf(filename = "5_DSpcoa.pdf" , width = 8,height = 8)
g <- plot_grid(plot_MAG_ecosystem,plot_MAG_region,plot_MAG_medium,
               plot_MAG_phylum,plot_MAG_DS_magsize,
               labels = letters[5], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd("D:/pku/20gaoyuan/10Defensefinder")
dev.off()






##NMDS####
# 标准化数据，使其适用于 Bray-Curtis 距离计算
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


count_matrix_standard = decostand(count_matrix_standard,method = "hellinger")  
# NMDS（非度量多维标度分析），k=2 表示降维到 2 维
mds_bray = metaMDS(count_matrix_standard,
                   distance = "bray", k = 2,trymax = 100)  


# 提取 NMDS 分析的 stress 值，四舍五入保留 4 位小数
# 计算 stress 值，越小表示 NMDS 结果越可靠
stress1 <- round(mds_bray$stress,4) 

# 提取 NMDS 计算出的坐标，并与样本分组信息合并
data.plot <- mds_bray$points %>% as.data.frame() %>%  # 转换为数据框
  rownames_to_column(var="genome") %>%  # 添加样本 ID 列
  left_join(., MAG_metadata, by="genome")  # 关联样本分组信息
data.plot[is.na(data.plot)] <- "Unknown"


# 绘制 Bray-Curtis NMDS 结果
p_nmds_ecosystem <- ggplot(data.plot, aes(MDS1, MDS2))+
  # 按 ecosystem 颜色区分样本点
  geom_point( aes(colour = ecosystem),size = 2, show.legend = F)+  
  # 95% 置信椭圆
  stat_ellipse(aes(label = ecosystem, color=ecosystem),
               geom = "textpath",
               hjust = 0.3, level = 0.95, linewidth=0.5,
               linetype=2, show.legend = F) +
  xlab("MDS1") + ylab("MDS2") +  # 轴标签
  # 添加 stress 值到标题
  labs(title = paste("bray stress =", stress1, sep=" ")) +  
  theme_test() +
  theme(plot.title = element_text(vjust=0.5, hjust = 0.5))  
p_nmds_ecosystem


##PCA####
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


pca <- count_matrix_standard %>% 
  prcomp(., scale. = TRUE)  # 进行标准化 PCA 计算
# 计算 PCA 解释的方差比例
var_explained <- pca$sdev^2 / sum(pca$sdev^2)

pca


DF_mechanism_pca <- data.frame(title = colnames(count_matrix_standard))
DF_mechanism_pca <- left_join(DF_mechanism_pca, DF_mechanism, by="title")
DF_mechanism_pca[is.na(DF_mechanism_pca)] <- "Unknown"
DF_mechanism_pca <- select(DF_mechanism_pca,title,Mechnism)


library(factoextra)
# 绘制 PCA 双标图（biplot）
p_pca <- fviz_pca_biplot(pca, axes = c(1, 2),  # 选择主成分 PC1 和 PC2
                      geom.ind = c("point"),  # 样本点
                      geom.var = c("arrow", "text"),  # 变量箭头 + 标签
                      pointshape = 20,  # 设定点的形状
                      pointsize = 4,  # 设定点的大小
                      # 变量标签避免重叠，颜色设为黑色
                      label = "var", repel = TRUE, col.var = "black",  
                      labelsize = 0.5,  # 变量标签的大小
                      addEllipses = TRUE,  # 添加 95% 置信椭圆
                      ellipse.level = 0.95,  # 设置置信区间水平
                      col.ind = DF_mechanism_pca$Mechnism  # 按 `Subtype` 分类颜色
) +
  scale_color_manual(values =
                       colorRampPalette(brewer.pal(12,"Paired"))(6)) +  
  # 设置 X 轴标签
  labs(x = paste0("(PC1: ",round(var_explained[1] * 100,2),"%)"),  
       y = paste0("(PC2: ",round(var_explained[2] * 100,2),"%)")) + 
  # 设置 X 轴标题样式
  theme(panel.background=element_rect(fill='white',colour='black'),  
        axis.title.x = element_text(
          colour = "black",size = 12,
          margin = margin(t = 12)), 
        # 设置 Y 轴标题样式
        axis.title.y = element_text(
          colour = "black", size = 12, margin = margin(r = 12)), 
        axis.text = element_text(color = "black"),  # 设置轴刻度文字颜色
        plot.title = element_blank(),  # 不显示标题
        legend.title = element_blank(),  # 不显示图例标题
        legend.key = element_blank(),  # 使图例背景透明
        legend.text = element_text(color = "black", size = 9),  
        # 设置图例文本之间的水平间距
        legend.spacing.x = unit(0.06, 'cm'),  
        legend.key.width = unit(0.01, 'cm'),  # 设置图例的水平大小
        legend.key.height = unit(0.01, 'cm'),  # 设置图例的垂直大小
        legend.background = element_blank(),  # 使图例背景透明
        legend.position = c(1, 0),  # 设置图例在图的右下角
        legend.justification = c(1, 0))  # 使图例对齐到右下角
p_pca


#根据pcoa结果绘制不同cluster的DS组成特点####

setwd("D:/pku/20gaoyuan/10Defensefinder")

pcoadata2 <- read.csv("pcoadata.csv")
pcoadata2 <- pcoadata2[pcoadata2$genome %in% MAG_metadata$genome,]

pcoadata2$type <- ifelse(pcoadata2$PC1 < 0.125, "cluster12", "cluster34")
pcoadata2_cluster12 <- pcoadata2[pcoadata2$type == "cluster12",]

pcoadata2_cluster2_part1 <- pcoadata2_cluster12[pcoadata2_cluster12$PC1>=0 & pcoadata2_cluster12$PC2 <= 0.0625,]
pcoadata2_cluster2_part2 <- pcoadata2_cluster12[pcoadata2_cluster12$PC1>=-0.1 &pcoadata2_cluster12$PC1<0 & pcoadata2_cluster12$PC2 <= 0,]
pcoadata2_cluster2_part3 <- pcoadata2_cluster12[pcoadata2_cluster12$PC1>=-0.25 & pcoadata2_cluster12$PC1< -0.1 & pcoadata2_cluster12$PC2 <= -0.125,]
pcoadata2_cluster2 <- rbind(pcoadata2_cluster2_part1, pcoadata2_cluster2_part2, pcoadata2_cluster2_part3)
pcoadata2_cluster2$type <- "Cluster2"

pcoadata2_cluster1 <- pcoadata2_cluster12[!pcoadata2_cluster12$genome 
                                          %in% pcoadata2_cluster2$genome,]
pcoadata2_cluster1$type <- "Cluster1"

pcoadata2_cluster34 <- pcoadata2[pcoadata2$type == "cluster34",]
pcoadata2_cluster3 <- pcoadata2_cluster34[pcoadata2_cluster34$PC2>0.125,]
pcoadata2_cluster4 <- pcoadata2_cluster34[pcoadata2_cluster34$PC2<0.125,]
pcoadata2_cluster3$type <- "Cluster3"
pcoadata2_cluster4$type <- "Cluster4"

pcoadata2_cluster <- rbind(pcoadata2_cluster1,pcoadata2_cluster2,pcoadata2_cluster3,pcoadata2_cluster4)

dss <- count_matrix
dss$genome <- rownames(count_matrix)
dss <- dss[dss$genome %in% pcoadata2_cluster$genome,]
dss <- merge(dss, pcoadata2_cluster, by = "genome")
dss <- select(dss, genome, type, everything())
dss <- select(dss, -PC1, -PC2, -Region, -X)
#不同DS的均值/出现率

#均值
library(dplyr)
dss_average_cluster <- dss %>%
  select(-genome) %>%               # 移除 genome 列
  group_by(type) %>%                # 按 type 分组
  summarise(across(everything(),    # 对所有列计算均值
            ~mean(., na.rm = TRUE))) %>% # 使用匿名函数格式
  data.frame()
rownames(dss_average_cluster) <- dss_average_cluster$type
dss_average_cluster <- dss_average_cluster[,-1]

#绘制热图
library(pheatmap)

colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p_heatmap <- pheatmap(t(dss_average_cluster),
                #color = colorRampPalette(colors)(1000),
                #color = colorRampPalette(c('white','lightpink','red'))(1000),
                color = colors,
                cluster_cols = T,
                cluster_rows = F,
                #annotation_col = super_phyla[,-1],
                fontsize_row = 10
)
p_heatmap
setwd("D:/pku/20gaoyuan/10Defensefinder")

#出现率
dss_frequency_cluster <- dss %>%
  select(-genome) %>%               # 移除 genome 列
  group_by(type) %>%                # 按 type 分组
  summarise(across(everything(),    # 对所有列计算
            ~ mean(. != 0) * 100)) %>%   # 非零值比例转换为百分比
  data.frame(check.names = F)
rownames(dss_frequency_cluster) <- dss_frequency_cluster$type
dss_frequency_cluster <- dss_frequency_cluster[,-1]

#绘制热图
library(pheatmap)

colors <- c('white', colorRampPalette(c('#F1EDE4', '#3C4D61'))(9999))
p_heatmap <- pheatmap(t(dss_frequency_cluster),
                      #color = colorRampPalette(colors)(1000),
                      #color = colorRampPalette(c('white','lightpink','red'))(1000),
                      color = colors,
                      cluster_cols = T,
                      cluster_rows = F,
                      #annotation_col = super_phyla[,-1],
                      fontsize_row = 10
)
p_heatmap
setwd("D:/pku/20gaoyuan/10Defensefinder")

dss_frequency_cluster_filter <- dss_frequency_cluster[,colSums(dss_frequency_cluster)>3]

p_heatmap_65 <- pheatmap(t(dss_frequency_cluster_filter),
                      #color = colorRampPalette(colors)(1000),
                      #color = colorRampPalette(c('white','lightpink','red'))(1000),
                      color = colors,
                      cluster_cols = T,
                      cluster_rows = F,
                      #annotation_col = super_phyla[,-1],
                      fontsize_row = 10
)
p_heatmap_65
setwd("D:/pku/20gaoyuan/10Defensefinder")
ggsave("5_heatmap_DSfrequency.pdf",p_heatmap_65, width = 3, height = 8.5, units = "in")


#删除RM和SoFIC再次进行PCoA####
##PCoA######
count_matrix_standard <- count_matrix
count_matrix_standard <- count_matrix_standard[,colSums(count_matrix_standard)>=15]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]

count_matrix_standard <- count_matrix_standard[,!colnames(count_matrix_standard) %in% c("RM","SoFIC")]
count_matrix_standard <- count_matrix_standard[rowSums(count_matrix_standard) >0,]


library(vegan)
library(ape)
library(tidyverse)
library(magrittr)
library(RColorBrewer)

# 计算 Bray-Curtis 距离并进行 PCoA 分析
pcoa <- vegdist(count_matrix_standard, method = "bray") %>% 
  # 进行主坐标分析（PCoA），不进行负特征值校正
  pcoa(correction = "none")

summary(pcoa)

# 提取 PCoA 结果并合并分组信息
pcoadata <- data.frame(pcoa$vectors[,1], pcoa$vectors[,2]) %>%  
  # 获取 PC1 和 PC2 坐标
  set_colnames(c("PC1", "PC2")) %>%  # 重命名列名
  rownames_to_column(var = "genome") %>%  # 将行名转换为列 sample
  # 关联 Subtype 信息
  left_join(., MAG_metadata %>% 
              dplyr::select(genome, Region),by="genome")

pcoadata[is.na(pcoadata)] <- "Unknown"

#write.csv(pcoadata,"pcoadata_noRMSoFIC.csv")
pcoadata2 <- read.csv("pcoadata_noRMSoFIC.csv")

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

pcoadata2 <- read.csv("pcoadata_noRMSoFIC.csv")
pcoadata2 <- left_join(pcoadata2,mag_lineage_split,suffix = c("genome","genome"))
# 统计 phylum 列中每个门的计数
p_counts <- pcoadata2 %>%
  count(phylum)
# 找出计数小于 1000 的门
low_count_p <- p_counts %>%
  filter(n < 50) %>%
  pull(phylum)
# 将计数小于 1000 的门及其后续分类等级设置为 "others"
pcoadata2 <- pcoadata2 %>%
  mutate(across(phylum:species, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                         ifelse(phylum %in% low_count_p, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))
phylum <- table(pcoadata2$phylum) %>% as.data.frame()

#添加每个MAG的size大小
pcoadata3 <- pcoadata2 %>% left_join(., MAG_metadata, by="genome")
colnames(MAG_quality)[1] <- "genome"
pcoadata3 <- pcoadata3 %>% left_join(., MAG_quality, by="genome")
pcoadata3$MAG_size <- pcoadata3$Genome_Size/1000000
pcoadata3$magsize <- floor(pcoadata3$MAG_size)
pcoadata3$magsize <- factor(pcoadata3$magsize)

###1、根据ecosystem绘图####

# 颜色设置
color_palette <- colorRampPalette(brewer.pal(12, "Paired"))(6)

library(ggplot2)
#绘制pcoa散点图
plot_MAG_ecosystem <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = ecosystem, fill = ecosystem), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","s"= "#C4D8E9","Benthon"="#c38e70",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_ecosystem

###2、region着色####


#绘制pcoa散点图
plot_MAG_region <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Region.x, fill = Region.x), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_region

###3、medium着色####
#绘制pcoa散点图
plot_MAG_medium <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = medium, fill = medium), size = 0.5) +  
  #  scale_color_manual(values = color_palette) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_medium

###4、物种着色####
phylum <- table(pcoadata2$phylum) %>% as.data.frame()
color_palette_phylum <- colorRampPalette(brewer.pal(12, "Paired"))(10)


#绘制pcoa散点图
plot_MAG_phylum <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = phylum, fill = phylum), size = 0.5) +  
  scale_color_manual(values = color_palette_phylum) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background = element_blank(),
        #   element_rect(fill = 'NA',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_phylum

###5、用MAG size着色####


color_palette_magsize <- colorRampPalette(brewer.pal(12, "Paired"))(12)
#绘制pcoa散点图
plot_MAG_DS_magsize <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = magsize, fill = magsize), size = 0.5) +  
  scale_color_manual(values = color_palette_magsize) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  labs(x = paste0("(PC1:",
                  round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
       y = paste0("(PC2:",
                  round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_magsize


###6、用season着色####


color_palette_season <- colorRampPalette(brewer.pal(12, "Paired"))(4)
#绘制pcoa散点图
plot_MAG_DS_season <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Season, fill = Season), size = 0.5) +  
  scale_color_manual(values = color_palette_season) +  # 使用手动颜色填充
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_season

###7、用domain着色####

pcoadata3$domain <- factor(pcoadata3$domain, levels = c("d__Archaea","d__Bacteria"))
#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_domain <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = domain, fill = domain), size = 0.5) +  
  #  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
                                "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
                                "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
  geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_domain

###8、用Elevation着色####

#color_palette_domain <- colorRampPalette(brewer.pal(12, "Paired"))(2)
#绘制pcoa散点图
plot_MAG_DS_Elevation <- ggplot(pcoadata3, aes(PC1, PC2)) +
  # 颜色映射到 Subtype
  geom_point(aes(colour = Elevation, fill = Elevation), size = 0.5) +  
  #  scale_color_manual(values = color_palette_domain) +  # 使用手动颜色填充
  # scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF","Plankton"="#c7e9c0",
  #                               "LCJ"="#FFE2CE","Benthon"="#c38e70","Plateau"= "#C4D8E9",
  #                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b",
  #                               "d__Archaea" = "#c1121f", "d__Bacteria" = "#ffe1a8"
  # )) +  # 手动设置颜色
  # 添加 PC1 解释度
  # labs(x = paste0("(PC1:",
  #                 round(pcoa$values$Relative_eig[1]*100,2),"%)"), 
  #      y = paste0("(PC2:",
  #                 round(pcoa$values$Relative_eig[2]*100,2),"%)")) +
geom_vline(aes(xintercept = 0),linetype= "dotted") +  # 添加纵向虚线
  geom_hline(aes(yintercept = 0),linetype= "dotted") +  # 添加横向虚线
  # stat_ellipse(aes(label = Region, color=Region),
  #              geom = "textpath",  
  #              hjust = 0.3, level = 0.95, linewidth=0.5, 
  #              linetype=2, show.legend = F) +
  theme_bw() +  # 使用黑白主题
  # 设置背景颜色
  theme(panel.background =
          element_rect(fill = 'white',colour = 'black'), 
        # X 轴标题格式
        axis.title.x = 
          element_text(colour = "black",size = 12, 
                       margin = margin(t = 5), face = "bold"), 
        axis.title.y =
          element_text(colour = "black",size = 12,
                       margin = margin(r = 5), face = "bold"),
        axis.text = element_text(color = "black", face = "bold"), 
        plot.title = element_blank(),  # 移除标题
        legend.title = element_blank(),  # 移除图例标题
        legend.key = element_blank(),  # 图例键为空
        legend.text = element_text(
          color = "black",size = 9, face = "bold"),  # 图例文本格式
        legend.position = c(0.001, 0.999),
        legend.background = element_blank(),
        legend.justification = c(0.0001, 1))  # 调整图例位置
plot_MAG_DS_Elevation


library(Cairo)
library(cowplot)
cairo_pdf(filename = "5_DSpcoa_noRMSoFIC.pdf", width = 8,height = 8)
g <- plot_grid(plot_MAG_ecosystem,plot_MAG_region,plot_MAG_medium,
               plot_MAG_phylum,plot_MAG_DS_magsize,
               labels = letters[5], ncol = 2
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd("D:/pku/20gaoyuan/10Defensefinder")
dev.off()










