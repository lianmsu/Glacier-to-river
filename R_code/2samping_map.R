
library(readxl)
library(ggspatial)
library(ggplot2)
library(cowplot)
library(colorspace)
library(sf)
library(dplyr)


### 设置工作路径
setwd("D:/pku/20gaoyuan/01MAG_info")

#metadata
MAG5455_metadata <- read.csv("5455MAG_metadata.csv")

#所有样品来自的国家的统计
sample_info_ecosystem <- MAG5455_metadata %>% count(ecosystem, name = "Count")

custom_colors <- c(
  "#9ecae1",  # 淡蓝色
  "#fdae6b",  # 淡橙色
  "#a1d99b",  # 淡绿色
  "#fc9272",  # 淡红色
  "#bcbddc",  # 淡紫色
  "#ffeda0",  # 淡黄色
  "#c6dbef",  # 浅蓝色
  "#fdd0a2",  # 浅橙色
  "#c7e9c0",  # 浅绿色
  "#fcbba1",  # 浅红色
  "#dadaeb",  # 浅紫色
  "#ffa69e"   # 浅黄色
)
# 绘制水平柱状图
sample_ecosystem <- 
  ggplot(sample_info_ecosystem, aes(x = ecosystem, y = Count, fill = ecosystem)) +
  geom_bar(stat = "identity",width = 0.8) +
  geom_text(aes(label = Count), hjust = -0.2, size = 3.5, color = "black") +  # 在柱子右侧显示数量
  coord_flip() +  # 将柱状图水平显示
  scale_y_continuous(limits = c(0, 3000), expand = c(0, 0)) +  # 设置 y 轴范围和扩展
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  labs(title = "",
       x = "",
       y = "Number of genomes") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # 标题居中，加粗，字体大小为16
    panel.grid.major = element_blank(),  # 去除主要网格线
    panel.grid.minor = element_blank(),  # 去除次要网格线
    axis.text.x = element_blank(),  # 去除纵轴刻度标签
    axis.ticks.x = element_blank(),  # 去除纵轴刻度线
    #axis.ticks.y = element_line(color = "black"),  # 保留 x 轴刻度线
    legend.position = "none"  # 不显示图例
  )
sample_ecosystem
setwd("D:/pku/20gaoyuan/01MAG_info")
ggsave("genome_ecosystem.pdf",sample_ecosystem,width = 4.5,height = 2.5,units = 'in')





#边界#############
# 加载必要的包
library(ggplot2)
library(sf)
library(ggspatial)
library(viridis)
### 设置工作路径
setwd("D:/pku/20gaoyuan/01MAG_info")

#metadata
MAG5455_metadata <- read.csv("5455MAG_metadata.csv")


### 样点数据包括：经度lon和纬度lat信息
scatter_df_tro <- st_as_sf(MAG5455_metadata, coords = c("Longitude", "Latitude"), crs = 4326) ## 设定参考坐标系
scatter_df_tro <- scatter_df_tro[,c("ecosystem","Region","geometry")] %>% distinct()

# 读取中国边界数据
china <- sf::read_sf('https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json')

# 假设dem_df已经准备好，包含x, y和elevation列


## 绘制中国大陆地图并添加采样点信息
china_map_ecosystem <- ggplot() +
  # 添加高程点
  # geom_sf(data = dem_sf, 
  #         aes(color = elevation),
  #         shape = 20,
  #         size = 0.5,
  #         alpha = 0.6) +
  # 设置高程颜色渐变
  # scale_color_gradientn(
  #   colors = c("#9BDFDF", "#FFE2CE", "#FB8C62"),
  #   name = "Elevation (m)"
  # ) +
  geom_sf(data = china, fill = 'NA', size = 0.2, color = 'black', lwd = 0.2) + # 中国地图轮廓
  # geom_sf(data = my_dili, aes(fill = QDABB1), size = 0.2, color = 'black', lwd = 0.5) + # 中国地图轮廓
  # scale_fill_manual(values = c("I"="#FCB2AF","Ⅱ"="#9BDFDF",
  #                              "Ⅲ"="#FFE2CE","Ⅳ"= "#C4D8E9",
  #                              "Ⅴ"= "#BEBCDF","Ⅵ"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +
  # 调整样点信息的外观
  geom_sf(data = scatter_df_tro, 
          aes(color = scatter_df_tro$ecosystem, fill = scatter_df_tro$ecosystem),
          shape = 21,           # 设置形状为带边框的圆形
          size = 1.5,           # 增大点的大小
          stroke = 0.5,         # 设置点边界的宽度
        #  color = 'blue',      # 边界颜色设置为黑色
        #  fill = '#00b4d8',        # 填充颜色为蓝色
          alpha = 0.7) +        # 设置透明度
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +  # 手动设置颜色
  # geom_sf(data = scatter_df_tro, 
  #         #aes(color=`Geo-environmental zone`),
  #         color = 'blue', 
  #         shape = 20,
  #         size = 0.8) + # 添加样点信息
 # geom_sf(data = china[c(21,32,35),], fill = 'NA', size = 0.2,lwd = 0.2, color = 'black') + 
  coord_sf(ylim = c(-2387082, 1654989), crs = '+proj=laea +lat_0=40 +lon_0=104') + #截取中国地图，不在大图上显示九段线区域
  # 添加比例尺
  annotation_scale(location = 'bl') +
  # 添加指北针
  annotation_north_arrow(location = 'tl', which_north = 'false', style = north_arrow_nautical()) +
  # 添加审图号并设定位置
  #labs(caption = '审图号: GS (2019) 1719号')+
  theme(plot.caption = element_text(hjust = 0, vjust = 0.0, size = 6.0)) + ## 设置主题的位置和大小，face = 'bold'为字体加粗
  # 设置坐标轴字体
  theme_linedraw()+
  theme(text = element_text(family = 'sans',size = 10, face = "bold"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.06), ## 设置轴标签大小
        legend.key = element_rect(fill = "NA"),
        #        legend.position = 'bottom',
        legend.position = "none")
china_map_ecosystem

## 绘制南海九段线小图
nine_map <- ggplot() +
  geom_sf(data = china, fill = 'NA', size = 0.2, color = 'black', lwd = 0.2) + # 中国地图轮廓
  # scale_fill_manual(values = c("I"="#FCB2AF","Ⅱ"="#9BDFDF",
  #                              "Ⅲ"="#FFE2CE","Ⅳ"= "#C4D8E9",
  #                              "Ⅴ"= "#BEBCDF","Ⅵ"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +
  # 调整样点信息的外观
  
  #geom_sf(data = china, fill = 'NA', size = 0.2, color = 'black') + 
  #geom_sf(data = nine_line, color = 'black', size = 0.5)+
  coord_sf(ylim = c(-4028017, -1577844), xlim = c(117131.4, 2115095), crs = "+proj=laea +lat_0=40 +lon_0=104")+
  theme_linedraw()+
  theme(
    # aspect.ratio = 1.25, #调节长宽比
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    #panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    #panel.border = element_rect(fill = NA, color = "grey10", linetype = 1, linewidth = 1.0),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    legend.position = "none"
  )
nine_map

## 使用cowplot包将中国大陆地图和南海九段线地图拼在一起
gg_inset_map_ecosystem <- ggdraw() +
  draw_plot(china_map_ecosystem,x = 0, y = 0, width = 1, height = 1) +
  draw_plot(nine_map, x = 0.885, y = 0.038, width = 0.1, height = 0.20)
gg_inset_map_ecosystem #Times字体格式无法保存

ggsave("2sample_map_ecosystem.pdf",gg_inset_map_ecosystem,width = 4,height = 3.3,units = 'in')



#Region#####
MAG_metadata <- read.csv("D:/pku/20gaoyuan/01MAG_info/5455MAG_metadata.csv")
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'Yes',"Plateau",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'No',"Plain",MAG_metadata$Plateau)
MAG_metadata$Plateau <- ifelse(MAG_metadata$Plateau == 'glacier',"Glacier",MAG_metadata$Plateau)
colnames(MAG_metadata)[17] <- 'Region'
### 样点数据包括：经度lon和纬度lat信息
scatter_df_tro_Region <- st_as_sf(MAG_metadata, coords = c("Longitude", "Latitude"), crs = 4326) ## 设定参考坐标系
scatter_df_tro_Region <- scatter_df_tro_Region[,c("ecosystem","Region","geometry")] %>% distinct()


## 绘制中国大陆地图并添加采样点信息
china_map_Region <- ggplot() +
  # 添加高程点
  # geom_sf(data = dem_sf, 
  #         aes(color = elevation),
  #         shape = 20,
  #         size = 0.5,
  #         alpha = 0.6) +
  # 设置高程颜色渐变
  # scale_color_gradientn(
  #   colors = c("#9BDFDF", "#FFE2CE", "#FB8C62"),
  #   name = "Elevation (m)"
  # ) +
geom_sf(data = china, fill = 'NA', size = 0.2, color = 'black', lwd = 0.2) + # 中国地图轮廓
  # geom_sf(data = my_dili, aes(fill = QDABB1), size = 0.2, color = 'black', lwd = 0.5) + # 中国地图轮廓
  # scale_fill_manual(values = c("I"="#FCB2AF","Ⅱ"="#9BDFDF",
  #                              "Ⅲ"="#FFE2CE","Ⅳ"= "#C4D8E9",
  #                              "Ⅴ"= "#BEBCDF","Ⅵ"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +
  # 调整样点信息的外观
  geom_sf(data = scatter_df_tro_Region, 
          aes(color = scatter_df_tro_Region$Region, fill = scatter_df_tro_Region$Region),
          shape = 21,           # 设置形状为带边框的圆形
          size = 1.5,           # 增大点的大小
          stroke = 0.5,         # 设置点边界的宽度
          #  color = 'blue',      # 边界颜色设置为黑色
          #  fill = '#00b4d8',        # 填充颜色为蓝色
          alpha = 0.7) +        # 设置透明度
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  # geom_sf(data = scatter_df_tro, 
  #         #aes(color=`Geo-environmental zone`),
  #         color = 'blue', 
  #         shape = 20,
  #         size = 0.8) + # 添加样点信息
  # geom_sf(data = china[c(21,32,35),], fill = 'NA', size = 0.2,lwd = 0.2, color = 'black') + 
  coord_sf(ylim = c(-2387082, 1654989), crs = '+proj=laea +lat_0=40 +lon_0=104') + #截取中国地图，不在大图上显示九段线区域
  # 添加比例尺
  annotation_scale(location = 'bl') +
  # 添加指北针
  annotation_north_arrow(location = 'tl', which_north = 'false', style = north_arrow_nautical()) +
  # 添加审图号并设定位置
  #labs(caption = '审图号: GS (2019) 1719号')+
  theme(plot.caption = element_text(hjust = 0, vjust = 0.0, size = 6.0)) + ## 设置主题的位置和大小，face = 'bold'为字体加粗
  # 设置坐标轴字体
  theme_linedraw()+
  theme(text = element_text(family = 'sans',size = 10, face = "bold"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.06), ## 设置轴标签大小
        legend.key = element_rect(fill = "NA"),
        #        legend.position = 'bottom',
        legend.position = "none")
china_map_Region

## 绘制南海九段线小图
nine_map <- ggplot() +
  geom_sf(data = china, fill = 'NA', size = 0.2, color = 'black', lwd = 0.2) + # 中国地图轮廓
  # scale_fill_manual(values = c("I"="#FCB2AF","Ⅱ"="#9BDFDF",
  #                              "Ⅲ"="#FFE2CE","Ⅳ"= "#C4D8E9",
  #                              "Ⅴ"= "#BEBCDF","Ⅵ"="#FB8C62","Ⅶ"="#FFF3CA"
  # )) +
  # 调整样点信息的外观
  
  #geom_sf(data = china, fill = 'NA', size = 0.2, color = 'black') + 
  #geom_sf(data = nine_line, color = 'black', size = 0.5)+
  coord_sf(ylim = c(-4028017, -1577844), xlim = c(117131.4, 2115095), crs = "+proj=laea +lat_0=40 +lon_0=104")+
  theme_linedraw()+
  theme(
    # aspect.ratio = 1.25, #调节长宽比
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    #panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    #panel.border = element_rect(fill = NA, color = "grey10", linetype = 1, linewidth = 1.0),
    plot.margin = unit(c(0, 0, 0, 0), "mm"),
    legend.position = "none"
  )
nine_map

## 使用cowplot包将中国大陆地图和南海九段线地图拼在一起
gg_inset_map_Region <- ggdraw() +
  draw_plot(china_map_Region,x = 0, y = 0, width = 1, height = 1) +
  draw_plot(nine_map, x = 0.885, y = 0.038, width = 0.1, height = 0.20)
gg_inset_map_Region #Times字体格式无法保存

ggsave("2sample_map_Region.pdf",gg_inset_map_Region,width = 4,height = 3.3,units = 'in')


#采样点的高程差#####################
elevation_longitude <- select(MAG_metadata,Longitude,Elevation,ecosystem,Region,Type)

# 去除所有列的完全重复行
elevation_unique <- elevation_longitude %>%
  distinct()

#env_factor <- read.csv("factor_all.csv")

# 加载ggplot2包
library(ggplot2)

plot_sampling_elevation_ecosystem <- 
# 绘制散点图
ggplot(elevation_unique, aes(x = Longitude, y = Elevation, color = ecosystem)) +
  geom_point(size = 1.5, alpha = 0.7) +  # 调整点的大小和透明度
  labs(
    x = "Longitude (°E)", 
    y = "Elevation (m)", 
    color = "Ecosystem"
  ) +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none")+  # 图例放在底部
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_sampling_elevation_ecosystem

#ggsave("2sample_map_ecosystemlegend.pdf",plot_sampling_elevation_ecosystem,width = 4,height = 3.3,units = 'in')



plot_sampling_elevation_Region <- 
  # 绘制散点图
  ggplot(elevation_unique, aes(x = Longitude, y = Elevation, color = Region)) +
  geom_point(size = 1.5, alpha = 0.7) +  # 调整点的大小和透明度
  labs(
    x = "Longitude (°E)", 
    y = "Elevation (m)", 
    color = "Ecosystem"
  ) +
  scale_color_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                                "LCJ"="#FFE2CE","Plateau"= "#C4D8E9",
                                "CJ"= "#BEBCDF","HH"="#FB8C62","Plain"="#d4e09b"
  )) +  # 手动设置颜色
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none")+  # 图例放在底部
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_sampling_elevation_Region

#ggsave("2sample_map_Regionlegend.pdf",plot_sampling_elevation_Region,width = 4,height = 3.3,units = 'in')

plot_sampling_elevation_medium <- 
  # 绘制散点图
  ggplot(elevation_unique, aes(x = Longitude, y = Elevation, color = Type)) +
  geom_point(size = 1.5, alpha = 0.7) +  # 调整点的大小和透明度
  labs(
    x = "Longitude (°E)", 
    y = "Elevation (m)", 
    color = "Ecosystem"
  ) +
  scale_color_manual(values =  c("Snow"="#E64B35FF","Ice"="#4DBBD5FF",
                                 "Cryoconite"="#00A087FF","water"= "#3C5488FF",
                                 "particle"= "#F39B7FFF","sediment"="#8491B4FF")) +  # 手动设置颜色
  theme_minimal() +  # 使用简洁主题
  theme(legend.position = "none")+  # 图例放在底部
  theme(
    panel.background = element_blank(),  # 设置背景为透明
    panel.grid.major = element_blank(),  # 去掉主要网格线
    panel.grid.minor = element_blank(),  # 去掉次要网格线
    #panel.border = element_blank(),  # 去掉边界线
    panel.border = element_rect(color = "grey20", fill = NA, size = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
    #axis.line = element_line(color = "grey20", size = 0.5),  # 设置坐标轴线条为深灰色，更加细致
    axis.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
    axis.text = element_text(size = 10, color = "grey20"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
    axis.text.x = element_text(angle = 0,vjust = 1,hjust = 0.5), #横轴旋转
    axis.ticks = element_line(color = "grey20", size = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
    legend.position = "none",  # 设置图例位置在右侧
    legend.title = element_text(face = "plain", size = 12, color = "grey20"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
    legend.text = element_text(size = 10, color = "grey20"),  # 设置图例文本字体大小为10，颜色为深灰色
    legend.background = element_blank(),  # 设置图例背景为透明
    legend.box.background = element_rect(color = "grey80", size = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
    legend.key = element_blank()  # 去掉图例键的背景
  )
plot_sampling_elevation_medium

#ggsave("2sample_map_mediumlegend.pdf",plot_sampling_elevation_medium,width = 4,height = 3.3,units = 'in')



library(Cairo)
library(cowplot)
cairo_pdf(filename = "2_sampling_elevation_ecosystem_Region_medium.pdf" , width = 8.55,height = 2.55)
g <- plot_grid(plot_sampling_elevation_ecosystem,plot_sampling_elevation_Region,plot_sampling_elevation_medium,
               labels = letters[3:4], ncol = 3
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd("D:/pku/20gaoyuan/01MAG_info")
dev.off()

#MAGs质量信息#########

MAG_quality <- read.delim2("quality_report_90_5_10kb.tsv")
MAG_quality <- MAG_quality[MAG_quality$Name %in% MAG_metadata$genome,]

MAG_quality$Completeness <- MAG_quality$Completeness %>% as.numeric()
MAG_quality$Contamination <- MAG_quality$Contamination %>% as.numeric()
MAG_quality$Coding_Density <- MAG_quality$Coding_Density %>% as.numeric()
MAG_quality$Contig_N50 <- MAG_quality$Contig_N50 %>% as.numeric()
MAG_quality$Genome_Size <- MAG_quality$Genome_Size %>% as.numeric()
MAG_quality$Total_Coding_Sequences <- MAG_quality$Total_Coding_Sequences %>% as.numeric()



mytheme <-  theme(
  panel.background = element_blank(),  # 设置背景为透明
  panel.grid.major = element_blank(),  # 去掉主要网格线
  panel.grid.minor = element_blank(),  # 去掉次要网格线
  panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),  # 添加边界线，颜色为深灰色，填充为透明，大小为0.5
  axis.title = element_text( lineheight = 12, color = "black"),  # 设置坐标轴标题字体为粗体，大小为12，颜色为深灰色
  axis.text = element_text(lineheight = 10, color = "black"),  # 设置坐标轴标签字体大小为10，颜色为深灰色
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),  # 横轴标签旋转
  axis.ticks = element_line(color = "black", linewidth = 0.5),  # 设置坐标轴刻度线颜色为深灰色，更加细致
  legend.position = 'none',
  #legend.position = c(0.2,0.75),  # 设置图例位置在右侧
  legend.title = element_text( lineheight = 12, color = "black"),  # 设置图例标题字体为粗体，大小为12，颜色为深灰色
  legend.text = element_text(lineheight = 10, color = "black"),  # 设置图例文本字体大小为10，颜色为深灰色
  legend.background = element_blank(),  # 设置图例背景为透明
  legend.box.background = element_rect(color = "grey80", linewidth = 0.5),  # 设置图例框背景为浅灰色，边框更加细致
  legend.key = element_blank()  # 去掉图例键的背景
)

completeness_density_plot <-
ggplot(MAG_quality, aes(x=Completeness)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue")+
  labs(
#    title = "Density Plot of MAG Completeness",
    x = "Completeness (%)",
    y = "Density"
  ) +
  mytheme

completeness_density_plot

Contamination_density_plot <-
  ggplot(MAG_quality, aes(x=Contamination)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue")+
  labs(
 #   title = "Density Plot of MAG Contamination",
    x = "Contamination (%)",
    y = "Density"
  ) +
  mytheme

Contamination_density_plot

Coding_Density_density_plot <-
  ggplot(MAG_quality, aes(x=Coding_Density)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue")+
  labs(
#    title = "Density Plot of MAG Coding Density",
    x = "Coding Density",
    y = "Density"
  ) +
  mytheme

Coding_Density_density_plot

Genome_size_density_plot <-
  ggplot(MAG_quality, aes(x=Genome_Size)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue")+
  labs(
#    title = "Density Plot of MAG Size",
    x = "MAG Size",
    y = "Density"
  ) +
  mytheme

Genome_size_density_plot


Total_Coding_Sequences_density_plot <-
  ggplot(MAG_quality, aes(x=Total_Coding_Sequences)) +
  geom_density(fill = "skyblue", alpha = 0.5, color = "blue")+
  labs(
 #   title = "Density Plot of MAG Total Coding Sequences",
    x = "Total Coding Sequences",
    y = "Density"
  ) +
  mytheme

Total_Coding_Sequences_density_plot


library(Cairo)
library(cowplot)
cairo_pdf(filename = "2_MAG_quality.pdf" , width = 2.78,height = 5.18)
g <- plot_grid(completeness_density_plot, Contamination_density_plot,
               Genome_size_density_plot, Total_Coding_Sequences_density_plot,
               Coding_Density_density_plot,
               labels = letters[5], ncol = 1
               #rel_heights = c(1, 1, 1.5),
               # label_x = 0.09,            # 标签x坐标（1为最右侧）
               #  label_y = 0.9             # 标签y坐标（1为最上方）
)
g
setwd("D:/pku/20gaoyuan/01MAG_info")
dev.off()






#高程###########
library(terra)
dem <- rast("D:/pku/20gaoyuan/01MAG_info/DEM数据-中国/DEM数据-中国/dem")
dem <- rast("D:/pku/20gaoyuan/01MAG_info/basins_dem/basins_dem")
plot(dem)

library(ggplot2)
library(viridis)

# 转换为数据框
dem_df <- as.data.frame(dem, xy = TRUE)
colnames(dem_df) <- c("x", "y", "elevation")

# 绘制高程图
ggplot(dem_df, aes(x, y, fill = elevation)) +
  geom_raster() +
  scale_fill_viridis(option = "terrain", name = "Elevation (m)") +
  coord_equal() +
  theme_minimal()

# 2. dem_df转为sf对象，并设置为WGS84
dem_sf <- dem_df %>%
  st_as_sf(coords = c("x", "y"), crs = 4326)  # 4326是WGS84


################################

set.seed(123)
dem_sf_sample <- dem_sf %>% sample_n(20000)  # 采样2000个点

# 修改地图代码
china_map <- ggplot() +
  # 添加高程点
  geom_sf(data = dem_sf_sample, 
          aes(color = elevation),
          shape = 20,
          size = 0.5,
          alpha = 0.6) +
  # 添加中国边界
  geom_sf(data = china, 
          fill = 'NA', 
          size = 0.2, 
          color = 'black', 
          lwd = 0.2) +
  # 添加样点信息
  geom_sf(data = scatter_df_tro, 
          aes(fill = ecosystem),
          shape = 21,
          size = 1.5,
          stroke = 0.5,
          color = 'NA',
          alpha = 0.7) +
  # 设置高程颜色渐变
  scale_color_gradientn(
    colors = c("#9BDFDF", "#FFE2CE", "#FB8C62"),
    name = "Elevation (m)"
  ) +
  # 设置生态系统颜色
  scale_fill_manual(values = c("YJ"="#FCB2AF","Glacier"="#9BDFDF",
                               "LCJ"="#FFE2CE","s"= "#C4D8E9",
                               "CJ"= "#BEBCDF","HH"="#FB8C62","Ⅶ"="#FFF3CA"
  )) +
  # 其他设置保持不变
  coord_sf(ylim = c(-2387082, 1654989), crs = '+proj=laea +lat_0=40 +lon_0=104') +
  annotation_scale(location = 'bl') +
  annotation_north_arrow(location = 'tl', which_north = 'false', style = north_arrow_nautical()) +
  theme_linedraw() +
  theme(text = element_text(family = 'sans', size = 10, face = "bold"),
        panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey80", linewidth = 0.06),
        legend.key = element_rect(fill = "NA"),
        legend.position = "right")  # 将图例位置改为右侧，以便显示高程图例

china_map


################

