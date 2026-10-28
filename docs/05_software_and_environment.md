# Software and environment

## Cluster / job scheduler

All `pipelines/` scripts were submitted with **Slurm** (`#SBATCH` headers). The production
runs used two queues, `cpu` and `fat`, on a shared cluster; job names, memory requests and
thread counts are kept in the scripts as a record of the actual runs.

The site-specific absolute paths of those runs have been replaced by variables, so each
script reads them from its environment (or from the `# Site-specific settings` block just
below the `#SBATCH` header):

| Variable | Meaning | Default in the script |
| --- | --- | --- |
| `LUSTRE_HOME` | home directory on the large-storage filesystem (project directories) | `/path/to/lustre/home` |
| `HOME_ROOT` | home directory on the login-node filesystem | `/path/to/home` |
| `DB_ROOT` | reference databases (geNomad, CheckV, CARD, ...) | `${LUSTRE_HOME}/database` |
| `GTR_DI_DIR` | directory with the per-island DefenseFinder TSV files | `/path/to/project/10DefenseFinder2.0.0_5467MAGs/DIs` |

To re-run a script elsewhere, export the variables and adapt:

* `#SBATCH --partition=...` (queue name),
* `source activate <env>` — replace with `conda activate <env>` on modern conda,
* `# #SBATCH --reservation=...` — commented out in the repository (the reservation name
  of the production account is site-specific); uncomment and set your own if needed.

```bash
export LUSTRE_HOME=/your/lustre/home
export HOME_ROOT=/your/home
export DB_ROOT=/your/databases
sbatch pipelines/01_genome_annotation/05_gtdbtk_classify.sh
```

Python is `python3` with `pandas` only; the table post-processing scripts
(`pipelines/06_table_postprocessing/`) are plain scripts, no scheduler required.

## Conda environments (as used)

The environment names are the ones referenced by the scripts; the exact package versions
should be taken as "whatever was current at the time of the runs". Pin them from the
manuscript Methods if bit-wise reproducibility matters.

| Environment | Used by | Software |
| --- | --- | --- |
| `prokka` | `01_genome_annotation/02_prokka_annotate_MAGs.sh`, `04_prodigal_merged_MAGs.sh` | Prokka, Prodigal 2.6.3 (`module load prodigal/2.6.3-GCCcore-11.3.0`) |
| `gtdbtk-2.4.0-test` | `01_genome_annotation/05_gtdbtk_classify.sh` | GTDB-Tk 2.4.0 |
| `genomad` | `02_mobile_genetic_elements/01_genomad.sh` | geNomad |
| `isescan` | `02_mobile_genetic_elements/02_isescan.sh` | ISEScan |
| `icefinder` | `02_mobile_genetic_elements/03_icefinder_batch_1.sh`, `04_icefinder_batch_2.sh` | ICEfinder2 |
| `integron_finder` | `02_mobile_genetic_elements/05_integron_finder.sh` | IntegronFinder |
| *(none)* | `02_mobile_genetic_elements/06_minced_crispr.sh` | minced |
| `defensefinder` | `03_defense_systems/01_defensefinder_v1.sh` | DefenseFinder |
| `defensefinder2.0.0` | `03_defense_systems/02_defensefinder_v2.sh` | DefenseFinder 2.0.0 (`--preserve-raw`) |
| `defense_prefictor` | `03_defense_systems/03_defensepredictor_batch_1.sh`, `04_defensepredictor_batch_2.sh` | DefensePredictor (`module load gcc/7.2.0`) |
| `checkv` | `04_viral_contigs/01_checkv.sh` | CheckV (db v1.5 via `CHECKVDB`) |
| `BACPHLIP` | `04_viral_contigs/02_bacphlip.sh` | BACPHLIP |
| `blast2.16` | `04_viral_contigs/04_blastn_spacers_vs_phage_plasmid.sh` | BLAST+ 2.16 (`makeblastdb`, `blastn`) |
| `fasttree` | `04_viral_contigs/05_fasttree_acr_phylogeny.sh` | muscle, trimAl (`-gt 0.5`), FastTree |
| `RGI_CARD` | `05_adaptive_traits/01_rgi_card.sh` | RGI with the CARD database |
| `eggnog` | `05_adaptive_traits/02_diamond_functional_screening.sh` | DIAMOND (`-k 1 -e 1e-5 --id 60 --subject-cover 90`) |
| `macrel` | `05_adaptive_traits/04_macrel.sh` | Macrel |
| `antiSMASH` | `05_adaptive_traits/05_antismash_batch_1.sh` | antiSMASH |
| `OGT_prediction` | `05_adaptive_traits/06_ogt_prediction.sh` | OGT_prediction |
| *(seqkit)* | `06_table_postprocessing/05_DI_nondG_filter.sh` | seqkit, awk, comm |

Database paths (geNomad db, CheckV db, CARD, MRG, mibPOP, VFG, PaCRISPR) appear inline in
the scripts; they must be downloaded separately.

## R

The analyses in `analyses/` are plain R scripts (no package project, no `renv`). They were
developed and run under R 4.x; all 33 scripts parse cleanly under R 4.6.0
(`Rscript -e 'parse(file = f)'`). Packages used, by frequency:

* data handling: `dplyr`, `tidyr`, `stringr`, `tidyverse`, `magrittr`, `data.table`, `readxl`, `reshape2`, `reshape`, `tibble`
* plotting: `ggplot2`, `cowplot`, `Cairo`, `ggpubr`, `ggsignif`, `gghalves`, `ggsci`, `ggrepel`, `ggbeeswarm`, `gggenes`, `ggsankeyfier`, `RColorBrewer`, `viridis`, `MetBrewer`, `scales`, `treemap`, `treemapify`, `pheatmap`, `ComplexHeatmap`, `circlize`, `plotly`
* statistics / ecology: `vegan`, `rstatix`, `FD`, `linkET`
* phylogeny / trees: `ape`, `itol.toolkit`
* maps: `sf`, `ggspatial`, `terra`

`Cairo` needs the system `libcairo`; `sf`/`terra` need GDAL/GEOS/PROJ; `itol.toolkit` is
installed from CRAN.

A minimal install:

```r
install.packages(c("dplyr", "tidyr", "stringr", "tidyverse", "ggplot2", "cowplot",
                   "Cairo", "ggpubr", "rstatix", "vegan", "ape", "readxl",
                   "pheatmap", "RColorBrewer", "viridis"))
# optional, only for specific modules
install.packages(c("itol.toolkit", "gggenes", "gghalves", "ggsankeyfier",
                   "treemapify", "ComplexHeatmap", "circlize", "linkET", "FD"))
```

## Reproducibility notes

* Scripts are self-contained: no shared helper file, no `source()` between analyses.
* Intermediate state is passed through `.RData` files written into the data directory
  (`save.image("./<script name>.RData")`), so re-running a module overwrites them. See the
  dependency table at the end of [01_pipeline_overview.md](01_pipeline_overview.md).
* Figures are written as PDF with explicit `width`/`height` in inches; the sizes match the
  submission requirements of the manuscript.
