# =============================================================================
# Glacier-to-river -- revision nc / Relative abundance of defense systems
# =============================================================================
# Weight the defense-system content by MAG relative abundance (TPM) to
# obtain community-level defensome profiles.
#
# Manuscript      : Revision
# Main outputs    : 05_defense_system_relative_abundance.RData
# Original script : 93_NC_comments/05_DS_relative_abundance/1_DS_relative_abundance.R
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

setwd(file.path(GTR_DATA, "93_NC_comments/05_DS_relative_abundance"))
library(dplyr)
library(stringr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(rstatix)
library(ggpubr)

# 样本信息及环境因子 ####
EF_info <- read.delim2("EF_info.tsv",sep = '\t')
EF_info$Region <- factor(EF_info$Region,
                         levels = c("Glacier","Plateau","Plain"))

# MAGs的相对丰度 ####
MAGs_4903_relative_abundance <- 
  read.delim2("all_abundance_tpm_bin.tsv")

# 将所有非Genome列转换为数值型
MAGs_4903_relative_abundance <- MAGs_4903_relative_abundance %>%
  mutate(across(-Genome, as.numeric))
# 1. 将第一列 'Genome' 设为行名，并移除该列
MAGs_4903_relative_abundance_2 <- MAGs_4903_relative_abundance
rownames(MAGs_4903_relative_abundance_2) <- MAGs_4903_relative_abundance_2$Genome
MAGs_4903_relative_abundance_2 <- MAGs_4903_relative_abundance_2[, -1]  # 删除 Genome 列
# 2. 转置（行名变列名，列名变行名）
t_MAGs_4903_relative_abundance_2 <- 
  as.data.frame(t(MAGs_4903_relative_abundance_2))
# 3. 将行名（原样本名）保存为新列 'Sample'
t_MAGs_4903_relative_abundance_2$Sample <- 
  rownames(t_MAGs_4903_relative_abundance_2)
# 4. 重新排列列，将 Sample 置于第一列
t_MAGs_4903_relative_abundance_2 <- 
  t_MAGs_4903_relative_abundance_2[, 
    c("Sample", setdiff(names(t_MAGs_4903_relative_abundance_2), "Sample"))]
# 5. 清除行名（以免冗余）
rownames(t_MAGs_4903_relative_abundance_2) <- NULL
# 处理样本名
t_MAGs_4903_relative_abundance_2$Sample <- 
  gsub(".sorted.TPM", "", t_MAGs_4903_relative_abundance_2$Sample)  # 可选：进一步简化
# sample_names <- gsub("_", " ", sample_names)      # 可选：用空格替代下划线
t_MAGs_4903_relative_abundance_2$Sample <- 
  sub("^[^.]*\\.", "", t_MAGs_4903_relative_abundance_2$Sample)

# DS 数据 ####
DF_gene_4903MAG <- read.delim2("merged_4903MAGs_defense_finder_genes.tsv")
DF_system_4903MAG <- read.delim2("merged_4903MAGs_defense_finder_systems.tsv")

DF_gene_4903MAG_DG <- 
  DF_gene_4903MAG[DF_gene_4903MAG$activity=="Defense",]
DF_system_4903MAG_DS <-
  DF_system_4903MAG[DF_system_4903MAG$activity=="Defense",]


# ==========================================================================
# DS 相对丰度整理 ####
# 思路：
#   用MAGs的TPM相对丰度值推导每个DS subtype在样本中的相对丰度。
#   逻辑链：
#     1. 每个MAG的TPM → 作为该MAG中每个DG（防御基因）的相对丰度
#     2. 每个DS（防御系统）的相对丰度 = 其所有DG相对丰度的平均值
#         → 由于同一MAG中一个DS的所有DG共享同一个TPM，平均值 = 该MAG的TPM
#     3. 每个样本中，某DS subtype的总相对丰度 =
#           Σ (每个MAG的TPM × 该MAG中该DS subtype的出现次数)
#     4. 汇总示例：
#        - MAG "CJ_14-10BD_WA-bin.108" 在WA样本中TPM=1047.867，
#          含1个DarTG系统 → 贡献 DarTG = 1047.867
#        - 若另一MAG含2个SoFic系统且TPM=500，则贡献 SoFic = 1000
#    最终结果：DS_relative_abundance_by_subtype.csv（311 DS subtypes × 682 样本）
#              DS_relative_abundance_by_subtype_transposed.csv（转置版本）
# ==========================================================================

# ---- 第1步：从sys_id中提取MAG名称 ----
# sys_id的格式为：CJ_14-10BD_WA-bin.108-144_DarTG_8530
# 其中bin.X-YYY表示MAG bin.X下的contig YYY，MAG名称为bin.X部分
DS_subtype_counts <- DF_system_4903MAG_DS |>
  mutate(
    MAG = str_extract(sys_id, "^.*?bin\\.\\d+")
  ) |>
  # 统计每个MAG中每种DS subtype的出现次数（一个MAG可能多次出现同一subtype）
  count(MAG, subtype, name = "n_DS")

# ---- 第2步：将DS计数与MAG相对丰度（宽格式）关联并计算 ----
# 将DS计数与MAG TPM数据按MAG名称左连接
# 然后对每个DS subtype，逐样本计算：
#   sum(MAG_TPM × n_DS)
#   即：该样本中所有MAG对该DS subtype的贡献之和

# 获取所有样本列名（排除"Genome"列）
sample_cols <- names(MAGs_4903_relative_abundance)[-1]

# 联合计算：对每个subtype，在所有样本中汇总（避免长格式转换，提高效率）
DS_abund_by_subtype <- DS_subtype_counts |>
  left_join(MAGs_4903_relative_abundance, by = c("MAG" = "Genome")) |>
  group_by(subtype) |>
  summarise(across(all_of(sample_cols),
                   ~ sum(.x * n_DS),
                   .names = "{.col}")) |>
  ungroup()

# ---- 第3步：保存结果 ----

# 格式1：DS subtypes × 样本（宽格式）
write.csv(DS_abund_by_subtype,
          file = "DS_relative_abundance_by_subtype.csv",
          row.names = FALSE)

# 格式2：样本 × DS subtypes（转置，便于与样本元数据联用）
DS_abund_t <- as.data.frame(t(DS_abund_by_subtype[, -1]))
names(DS_abund_t) <- DS_abund_by_subtype$subtype
DS_abund_t <- cbind(Sample = sample_cols, DS_abund_t)
write.csv(DS_abund_t,
          file = "DS_relative_abundance_by_subtype_transposed.csv",
          row.names = FALSE)

# ==========================================================================
# DS 相对丰度后续分析 ####
# 基于上一步得到的 DS_abund_by_subtype 进行 α 多样性、β 多样性和 Mantel 分析
# 参考：01_diversity_comparison/1_MAG_Phage_macrodiversity.R
#       04_environmental_factors/1_environmental_factor_effects.R
# ==========================================================================

# 避免 utf8.dll 加载失败导致的 tibble 打印问题
options(pillar.print_footnotes = FALSE)

# ---- 数据准备：将 DS 丰度转为样本×subtype 矩阵，并与 EF_info 关联 ----
DS_abund_matrix <- as.data.frame(t(DS_abund_by_subtype[, -1]))
rownames(DS_abund_matrix) <- names(DS_abund_by_subtype)[-1]
names(DS_abund_matrix) <- DS_abund_by_subtype$subtype

# 清理样本名：去掉前缀和 .sorted.TPM 后缀
# DS列名格式: "CJ_14.10BD_SAB.sorted.TPM" → 提取 "10BD_SAB"
sample_clean <- sub("\\.sorted\\.TPM$", "", rownames(DS_abund_matrix))
sample_clean <- sub("^[^.]+\\.", "", sample_clean)
rownames(DS_abund_matrix) <- sample_clean

# 按EF_info筛选共有样本
common_samples <- intersect(sample_clean, EF_info$SampleID)
cat("共有样本数:", length(common_samples), "\n")

DS_abund_matrix <- DS_abund_matrix[common_samples, ]

# 关联 EF_info 分组信息
DS_abund_meta <- data.frame(
  Sample = common_samples,
  Region = EF_info$Region[match(common_samples, EF_info$SampleID)]
)

# ==========================================================================
# α多样性分析 ####
# ==========================================================================
library(vegan)
library(tibble)

# 自定义 α 多样性计算函数（输入需为相对丰度，每行和为1）
calculate_alpha_diversity <- function(abund_matrix) {
  data.frame(
    Sample = rownames(abund_matrix),
    Observed_Richness = apply(abund_matrix, 1, function(x) sum(x > 0)),
    Shannon = diversity(abund_matrix, index = "shannon"),
    Simpson = diversity(abund_matrix, index = "simpson"),
    InvSimpson = diversity(abund_matrix, index = "invsimpson"),
    Pielou = apply(abund_matrix, 1, function(x) {
      p <- x[x > 0]
      S <- length(p)
      if (S <= 1) return(1)
      H <- -sum(p * log(p))
      H / log(S)
    })
  )
}

# 转换为相对丰度（每行和为1），并移除全零样本
DS_abund_rel <- DS_abund_matrix
zero_samples <- rowSums(DS_abund_rel) == 0
if (sum(zero_samples) > 0) {
  cat("移除", sum(zero_samples), "个全零样本:", 
      paste(rownames(DS_abund_rel)[zero_samples], collapse = ", "), "\n")
  DS_abund_rel <- DS_abund_rel[!zero_samples, ]
}
DS_abund_rel <- DS_abund_rel / rowSums(DS_abund_rel)

# 重新匹配分组信息（移除全零样本后）
DS_abund_meta <- DS_abund_meta[!zero_samples, ]

# 计算 α 多样性
DS_alpha <- calculate_alpha_diversity(DS_abund_rel)

# 合并分组信息
DS_alpha <- DS_alpha |>
  left_join(DS_abund_meta, by = "Sample")

# 设置 Region 因子水平
DS_alpha$Region <- factor(DS_alpha$Region,
                          levels = c("Glacier", "Plateau", "Plain"))

# ---- 统计检验（Wilcoxon，按 Region 两两比较）----
library(rstatix)

diversity_indices <- c("Pielou", "Shannon", "Simpson", "Observed_Richness")

DS_alpha_test <- DS_alpha |>
  pivot_longer(cols = all_of(diversity_indices),
               names_to = "index",
               values_to = "value") |>
  group_by(index) |>
  wilcox_test(value ~ Region) |>
  adjust_pvalue(p.col = "p", method = "bonferroni") |>
  add_significance(p.col = "p.adj") |>
  add_xy_position(x = "Region", dodge = 0.5) |>
  ungroup()

print(as.data.frame(DS_alpha_test))

# ---- 箱线图 ----
# 统一主题
my_theme <- theme(
  panel.background = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_rect(color = "grey20", fill = NA, linewidth = 0.5),
  axis.title = element_text(size = 12, color = "black"),
  axis.text = element_text(size = 10, color = "black"),
  axis.text.x = element_text(angle = 0, vjust = 0.8, hjust = 0.5),
  axis.ticks = element_line(color = "grey20", linewidth = 0.5),
  legend.position = "none"
)

# 颜色方案
region_colors <- c("Glacier" = "#9BDFDF",
                   "Plateau" = "#C4D8E9",
                   "Plain" = "#d4e09b")

# 各指数箱线图
p_DS_pielou <- ggplot(DS_alpha, aes(x = Region, y = Pielou)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(y = "Pielou evenness") +
  my_theme

p_DS_shannon <- ggplot(DS_alpha, aes(x = Region, y = Shannon)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(y = "Shannon diversity") +
  my_theme

p_DS_simpson <- ggplot(DS_alpha, aes(x = Region, y = Simpson)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(y = "Simpson diversity") +
  my_theme

p_DS_richness <- ggplot(DS_alpha, aes(x = Region, y = Observed_Richness)) +
  geom_boxplot(aes(fill = Region), alpha = 0.7, outlier.shape = NA) +
  geom_jitter(aes(color = Region), size = 1, width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(y = "Observed richness") +
  my_theme

# 组合图
library(cowplot)
g_DS_alpha <- plot_grid(p_DS_pielou, p_DS_shannon,
                        p_DS_simpson, p_DS_richness,
                        labels = letters[1:4], ncol = 2,
                        align = "v", axis = "lr")

ggsave(plot = g_DS_alpha,
       "2_DS_alpha_diversity.pdf",
       width = 6, height = 6, units = "in")

# ==========================================================================
# β多样性分析 ####
# ==========================================================================

# ---- Bray-Curtis 距离 ----
# 移除全零的样本行
nonempty <- rowSums(DS_abund_matrix) > 0
DS_bray <- vegdist(DS_abund_matrix[nonempty, ], method = "bray")
DS_bray_mat <- as.matrix(DS_bray)

# 匹配分组信息（过滤后）
beta_meta <- DS_abund_meta[nonempty, ]

# ---- PCoA ----
DS_pcoa <- cmdscale(DS_bray_mat, k = 3, eig = TRUE)
DS_pcoa_var <- DS_pcoa$eig / sum(DS_pcoa$eig[DS_pcoa$eig > 0]) * 100

DS_pcoa_df <- data.frame(
  Sample = rownames(DS_bray_mat),
  PC1 = DS_pcoa$points[, 1],
  PC2 = DS_pcoa$points[, 2],
  PC3 = DS_pcoa$points[, 3]
) |>
  left_join(DS_abund_meta, by = "Sample")

DS_pcoa_df$Region <- factor(DS_pcoa_df$Region,
                            levels = c("Glacier", "Plateau", "Plain"))

p_DS_pcoa <- ggplot(DS_pcoa_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Region, fill = Region), size = 2, alpha = 0.9) +
  stat_ellipse(aes(group = Region, color = Region),
               level = 0.68, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(x = sprintf("PCoA 1 (%.1f%%)", DS_pcoa_var[1]),
       y = sprintf("PCoA 2 (%.1f%%)", DS_pcoa_var[2])) +
  my_theme

# ---- NMDS ----
set.seed(943)
DS_nmds <- metaMDS(DS_bray_mat, k = 2, trymax = 100, distance = "bray")
cat("NMDS stress:", DS_nmds$stress, "\n")

DS_nmds_scores <- scores(DS_nmds, display = "sites")
DS_nmds_df <- data.frame(
  Sample = rownames(DS_nmds_scores),
  NMDS1 = DS_nmds_scores[, 1],
  NMDS2 = DS_nmds_scores[, 2]
) |>
  left_join(DS_abund_meta, by = "Sample")

DS_nmds_df$Region <- factor(DS_nmds_df$Region,
                            levels = c("Glacier", "Plateau", "Plain"))

p_DS_nmds <- ggplot(DS_nmds_df, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = Region, fill = Region), size = 2, alpha = 0.9) +
  stat_ellipse(aes(group = Region, color = Region),
               level = 0.68, linetype = "dashed", linewidth = 0.5) +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(subtitle = sprintf("stress = %.3f", DS_nmds$stress),
       x = "NMDS 1", y = "NMDS 2") +
  my_theme

# ---- ANOSIM ----
DS_anosim_region <- anosim(DS_bray,
                           grouping = beta_meta$Region,
                           permutations = 999)
print(as.data.frame(DS_anosim_region))

# 组合 β 多样性图
g_DS_beta <- plot_grid(p_DS_pcoa, p_DS_nmds,
                       labels = letters[1:2], ncol = 2,
                       align = "v", axis = "lr")

ggsave(plot = g_DS_beta,
       "2_DS_beta_diversity.pdf",
       width = 8, height = 4, units = "in")

# ==========================================================================
# Mantel 分析：DS 群落 vs 环境因子 ####
# 参考 1_environmental_factor_effects.R 中的 linkET 方法
# ==========================================================================
library(linkET)
library(FD)

# ---- 准备环境因子数据 ----
env_use <- EF_info |>
  select(SampleID, pH, NO3, DOC_TOC, Temperature, DO, Elevation,
         Heat, LE, NH4, srad, TP, River, Season, Type, Region) |>
  column_to_rownames("SampleID")

# ---- 准备 DS 丰度数据（与 env_use 取共有样本，移除全零行）----
DS_for_mantel <- DS_abund_matrix
common_mantel <- intersect(rownames(DS_for_mantel), rownames(env_use))

env_use <- env_use[common_mantel, ]
DS_for_mantel <- DS_for_mantel[common_mantel, ]

# 移除全零的 DS 样本行
nonzero_mantel <- rowSums(DS_for_mantel) > 0
if (sum(!nonzero_mantel) > 0) {
  cat("Mantel 分析中移除", sum(!nonzero_mantel), "个全零样本\n")
  DS_for_mantel <- DS_for_mantel[nonzero_mantel, ]
  env_use <- env_use[nonzero_mantel, ]
}

cat("Mantel 分析样本数:", length(common_mantel), "\n")

# ---- 执行 Mantel 检验 ----
DS_mantel_tbl <- mantel_test(
  spec = DS_for_mantel,
  env = env_use,
  spec_select = list(DS = 1:ncol(DS_for_mantel)),
  env_select = as.list(setNames(seq_len(ncol(env_use)),
                                colnames(env_use))),
  spec_dist = "bray",
  env_dist = "gower",
  na_omit = TRUE,
  permutations = 999
)

# ---- 整理 Mantel 结果 ----
DS_mantel_plot <- DS_mantel_tbl |>
  rename(from = env, to = spec, rd = r, pd = p) |>
  mutate(
    colour2 = case_when(
      abs(rd) < 0.3 | pd > 0.05 ~ "weak",
      rd > 0 ~ "positive",
      rd < 0 ~ "negative"
    ),
    linetype = ifelse(pd < 0.05, "solid", "dashed")
  )

# ---- 环境因子相关矩阵 ----
# 将数值型字符列转为 numeric（EF_info 中部分变量因读取格式为字符型）
num_vars <- c("pH", "NO3", "DOC_TOC", "Temperature", "DO", "Elevation",
              "Heat", "LE", "NH4", "srad", "TP")
env_cor <- env_use |>
  mutate(across(all_of(num_vars), ~ suppressWarnings(as.numeric(.)))) |>
  mutate(across(c(River, Season, Type, Region),
                ~ as.numeric(factor(.))))
env_correlation <- correlate(env_cor, method = "pearson",
                             use = "pairwise.complete.obs")

# ---- 绘图 ----
plot_DS_mantel <-
  qcorrplot(env_correlation,
            grid_col = "#00468BFF", "white", "#42B540FF",
            grid_size = 0.2,
            type = "upper", diag = FALSE) +
  geom_square() +
  scale_fill_gradientn(colours = c("#00468BFF", "white", "#42B540FF"),
                       limits = c(-1, 1)) +
  geom_mark(size = 4, only_mark = TRUE,
            sig_level = c(0.05, 0.01, 0.001),
            sig_thres = 0.05, colour = "white") +
  geom_couple(data = DS_mantel_plot,
              aes(colour = colour2, size = abs(rd), linetype = linetype),
              curvature = nice_curvature()) +
  scale_colour_manual(values = c(positive = "#ff9f1c",
                                 negative = "#377EB8",
                                 weak = "grey70")) +
  scale_size(range = c(0.5, 1)) +
  scale_linetype_identity()

plot_DS_mantel

ggsave(plot = plot_DS_mantel,
       "1_DS_mantel_environment.pdf",
       width = 6.2, height = 6.2, units = "in")

# ==========================================================================
# 保存工作空间
# ==========================================================================
save.image("./05_defense_system_relative_abundance.RData")
load("./05_defense_system_relative_abundance.RData")
# 清理
rm(sample_cols, sample_clean, common_samples, nonempty)

write.table(DS_mantel_plot, "DS_mantel_plot.tsv", sep = "\t",quote = F, row.names = F)
