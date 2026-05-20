
setwd("D:/pku/20gaoyuan/01MAG_info")

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



# 绘制桑基图展示物种组成###################

library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)

setwd("D:/pku/20gaoyuan/01MAG_info")

df <- all_MAG_tax

# 统计 k 列中每个k的计数
k_counts <- df %>%
  count(Domain)

# 统计 phylum 列中每个门的计数
p_counts <- df %>%
  count(Phylum)
# 找出计数小于 10 的门
low_count_p <- p_counts %>%
  filter(n < 10) %>%
  pull(Phylum)
# 将计数小于 1000 的门及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(Phylum:Genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                              ifelse(Phylum %in% low_count_p, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))

# 统计 class 列中每个门的计数
c_counts <- df %>%
  count(Class)
# 找出计数小于 30 的class
low_count_c <- c_counts %>%
  filter(n < 30) %>%
  pull(Class)
# 将计数小于 1000 的c及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(Class:Genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                              ifelse(Class %in% low_count_c, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))

# 统计 order 列中每个门的计数
o_counts <- df %>%
  count(Order)
# 找出计数小于 75 的class
low_count_o <- o_counts %>%
  filter(n < 75) %>%
  pull(Order)
# 将计数小于 75 的order及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(Order:Genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                      ifelse(Order %in% low_count_o, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))


# 统计 family 列中每个门的计数
f_counts <- df %>%
  count(Family)
# 找出计数小于 55 的family
low_count_f <- f_counts %>%
  filter(n < 55) %>%
  pull(Family)
# 将计数小于 55 的family及其后续分类等级设置为 "others"
df <- df %>%
  mutate(across(Family:Genus, ~ ifelse(. == "unknown", .,  # 保留 "unknown"
                                      ifelse(Family %in% low_count_f, "Others", .) # 将低计数门及后续分类设为 "Others"
  )))



library(tidyverse)
library(ggsankeyfier) 
library(MetBrewer)
df1 <- df %>% select(2,3) %>% group_by(Domain,Phylum) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Domain", "Phylum"),
                      values_from = "n")
df2 <- df %>% select(3,4) %>% group_by(Phylum,Class) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Phylum", "Class"),
                      values_from = "n")
df3 <- df %>% select(4,5) %>% group_by(Class,Order) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Class", "Order"),
                      values_from = "n")

df4 <- df %>% select(5,6) %>% group_by(Order,Family) %>% count() %>%
  pivot_stages_longer(.,stages_from = c("Order", "Family"),
                      values_from = "n")

eukaryotes_sankey_4level <- 
  ggplot(data=df1,aes(x = stage,y =n,group = node,
                      edge_id = edge_id,connector = connector))+
  # 绘制第 1，2 层级
  geom_sankeyedge(aes(fill = node),
                  position = position_sankey(order ="ascending",v_space="auto",
                                             width = 0.01))+
  geom_sankeynode(aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="from"),
            aes(label = node),stat = "sankeynode",
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3,fontface="plain",color="black")+
  # 绘制第 2，3 层级
  geom_sankeyedge(data=df2,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df2,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df1 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  geom_text(data=df2 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending",nudge_x=0.05),
            hjust=0,size=3,color="black",fontface="plain")+
  # 绘制第 3，4 层级
  geom_sankeyedge(data=df3,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df3,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df3 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  # 绘制第 4，5 层级
  geom_sankeyedge(data=df4,aes(fill = node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.01))+
  geom_sankeynode(data=df4,aes(fill=node,color=node),
                  position = position_sankey(order = "ascending",v_space ="auto",
                                             width = 0.05))+
  geom_text(data=df4 %>% filter(connector=="to"),
            aes(label = node),stat = "sankeynode",angle=0,
            position = position_sankey(v_space ="auto",order="ascending", nudge_x=0.05),
            hjust=0,size=3,vjust=0.5,color="black",fontface="plain")+
  
  coord_cartesian(clip="off")+
  scale_x_discrete(position = "bottom")+
  # scale_fill_manual(values = met.brewer("Nizami")) +
  # scale_color_manual(values = met.brewer("Nizami")) +
  theme_void()+
  theme(plot.margin = margin(0,0,0,0,unit = "in"),
        axis.text.x=element_text(color="black",face="plain",size=10,
                                 margin = margin(b=0.1,unit = "cm")),
        legend.position="none")

eukaryotes_sankey_4level

setwd("D:/pku/20gaoyuan/01MAG_info")
ggsave("3_MAG_tax_sankey.pdf",eukaryotes_sankey_4level,width = 5.2,height = 5.2,units = 'in')
