# Glacier-to-river

> **Microbial evolutionary arms races across glacier-to-river gradients** 论文的代码仓库
> （Li P., Li Y., Ni J.）

本仓库收录了从公共宏基因组数据到论文全部图件的所有代码：上游生物信息学流程（MAG 目录、
移动遗传元件、防御/抗防御系统注释）与下游 R 统计分析和绘图。

## 目录结构

```
Glacier-to-river/
├── pipelines/          # 上游流程：SLURM 作业脚本 + conda 环境（shell / python）
│   ├── 01_genome_annotation/         MAG 重命名、Prokka/Prodigal、GTDB-Tk
│   ├── 02_mobile_genetic_elements/   geNomad、ISEScan、ICEfinder2、IntegronFinder、minced
│   ├── 03_defense_systems/           DefenseFinder 1/2、DefensePredictor
│   ├── 04_viral_contigs/             CheckV、BACPHLIP、BLAST、FastTree
│   ├── 05_adaptive_traits/           RGI/CARD、DIAMOND、Macrel、antiSMASH、OGT
│   └── 06_table_postprocessing/      防御岛相关中间表（python + awk）
├── analyses/           # 下游分析：R 统计与绘图
│   ├── 01_MAG_catalog/               图 1、图 S1、图 S7
│   ├── 02_defensome_landscape/       图 2、图 3、图 S2–S4
│   ├── 03_mobile_genetic_elements/   图 4、图 S5、S6、S10
│   ├── 04_defense_islands/           图 5、图 S11–S13
│   ├── 05_antidefensome/             图 6、图 S14–S16
│   ├── 06_adaptive_traits/           ARG / MRG / VFG / POP 降解等适应性特征
│   ├── 07_cross_study_comparison/    河流 vs 含水层（与本组另一项研究比较）
│   └── 08_revision_nc/               审稿过程中新增的分析
└── docs/               # 流程、目录、文件对照表、数据与软件说明
```

两层目录均按执行顺序编号；每个脚本开头都有一段头部注释，写明用途、对应图件、主要输出
文件以及它在原始工作目录中的位置。

## 如何复现

**1. 上游流程**：`pipelines/` 下是当时实际投递到集群的 SLURM 作业脚本（保留 `#SBATCH`
头），按 [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) 的顺序提交。脚本中
与运行环境相关的绝对路径已改为变量：设置 `LUSTRE_HOME`、`HOME_ROOT`、`DB_ROOT`（见每个
脚本 `#SBATCH` 头下方的 `# Site-specific settings` 段），改好队列名即可，conda 环境见
[docs/05_software_and_environment.md](docs/05_software_and_environment.md)。

**2. R 分析**：需要 R ≥ 4.2 及 docs/05 中列出的 R 包。脚本通过 `GTR_DATA` 定位输入数据，
默认取当前工作目录，无需改代码即可切换：

```r
Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
source("analyses/02_defensome_landscape/02_region_comparison.R")
```

每个模块内的脚本按编号顺序运行；部分脚本通过 `.RData` 向下游传递中间结果
（例如 `analyses/02_defensome_landscape/05_defense_system_composition_pcoa.RData`）。

## 数据

MAG、病毒目录以及防御/抗防御注释结果随论文以 Supplementary Data 形式提供；每张分析表
的来源见 [docs/04_data_availability.md](docs/04_data_availability.md)。原始宏基因组测序
数据在 NCBI SRA 的登录号见论文。

## 关于旧版仓库

旧版本的仓库是把 `*.sh`、`*.py`、`R_code/*.R` 平铺在一起的。本次整理为上面的分层结构：
没有丢弃任何分析（脚本做了重命名与路径可移植化，并补充了旧版缺失的 10 个分析，含 5 个
审稿意见分析），只删除了重复副本和空白草稿（6 个 shell 脚本各有两份，另有 4 个陈旧或空
文件）。新旧文件的对应关系见 [docs/03_file_mapping.md](docs/03_file_mapping.md)。

仓库中已不含任何绝对路径、个人主目录、账号名或本地 Windows 目录：所有与运行环境相关的
位置都由 `LUSTRE_HOME`、`HOME_ROOT`、`DB_ROOT`、`GTR_DATA`、`GTR_GW`、`GTR_DI_DIR`
环境变量控制（说明见
[docs/03](docs/03_file_mapping.md#site-specific-paths-and-how-they-are-hidden)）。

## 引用

Li, P., Li, Y. & Ni, J. Microbial evolutionary arms races across glacier-to-river
gradients.（论文信息将在正式发表后更新）

## 联系

通讯作者：倪晋仁（北京大学环境科学与工程学院）。具体联系方式见论文；代码相关问题请在
本仓库提交 issue。
