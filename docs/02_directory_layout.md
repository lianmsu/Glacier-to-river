# Directory layout

The repository is split into two layers plus documentation. Both layers use the same
conventions: a two-digit directory prefix gives the execution / narrative order, and a
two-digit file prefix gives the order inside a directory.

```
Glacier-to-river/
├── README.md                  English overview
├── README_zh.md               中文说明
├── docs/                      workflow, layout, mapping, data, software
├── pipelines/                 layer 1: bioinformatics, SLURM + conda
└── analyses/                  layer 2: statistics and figures, R
```

## `pipelines/` — upstream, ordered by data type

| Directory | Scope |
| --- | --- |
| `01_genome_annotation/` | contig renaming, Prokka/Prodigal gene prediction, GTDB-Tk taxonomy |
| `02_mobile_genetic_elements/` | geNomad (phage, plasmid), ISEScan (IS), ICEfinder2 (ICE/IME), IntegronFinder, minced (CRISPR) |
| `03_defense_systems/` | DefenseFinder (two passes) and DefensePredictor (two MAG batches) |
| `04_viral_contigs/` | CheckV, BACPHLIP, BLAST databases and spacer mapping, FastTree Acr phylogeny |
| `05_adaptive_traits/` | RGI/CARD, DIAMOND database screening, Macrel, antiSMASH, OGT prediction |
| `06_table_postprocessing/` | defense-island tables: python scripts and one awk/seqkit shell step |

The scripts are the SLURM job scripts as they were submitted, with their `#SBATCH`
headers and absolute paths preserved, so that the exact computational environment is
documented. Adapt paths and queue names before re-submitting.

## `analyses/` — downstream, ordered by manuscript narrative

| Directory | Manuscript | Content |
| --- | --- | --- |
| `01_MAG_catalog/` | Fig. 1; Figs. S1, S7 | MAG metadata, sampling map, taxonomy, bacterial/archaeal iTOL trees, supplementary MAG table |
| `02_defensome_landscape/` | Fig. 2, Fig. 3; Figs. S2-S4 | defensome overview, region/medium/season comparison, density vs genome size, DS composition, phylum patterns |
| `03_mobile_genetic_elements/` | Fig. 4; Figs. S5, S6, S10 | MGE counts and cargo, anti-defense genes on MGEs, representative ICE/IS/prophage elements |
| `04_defense_islands/` | Fig. 5; Figs. S11-S13 | DI characterization, anti-defense inside DIs, non-defense-gene annotation, DefenseFinder vs DefensePredictor |
| `05_antidefensome/` | Fig. 6; Figs. S14-S16 | phage antidefensome, representative plasmid/phage arms races, Acr phylogeny |
| `06_adaptive_traits/` | - | ARG / MRG / VFG / POP-degradation / antimicrobial-peptide traits per MAG |
| `07_cross_study_comparison/` | - | river/glacier defensome vs the companion aquifer study |
| `08_revision_nc/` | peer review | macro-diversity, VHR and prophage abundance, NCBI RefSeq defense islands, environmental factors, abundance-weighted defensome |

### File naming inside `analyses/`

`NN_<what it does>.R` — the number is the run order within the module, the rest says what
the script produces, not which software it calls. Examples:

* `02_defensome_landscape/02_region_comparison.R` — glacier vs plateau vs plain;
* `04_defense_islands/04_defensefinder_vs_defensepredictor.R` — tool comparison;
* `08_revision_nc/05_defense_system_relative_abundance.R` — abundance-weighted defensome.

### Script header

Every R script starts with the same block, which is also the fastest way to navigate the
repository:

```r
# Glacier-to-river -- <module> / <title>
# <one-sentence purpose>
#
# Manuscript      : Fig. 3a, Fig. 3e, Fig. S2
# Main outputs    : 2_genome_Region.pdf, ...
# Original script : 20gaoyuan/10Defensefinder/2defensefinder_res_Glacier_Plateau_Plain.R
GTR_DATA <- normalizePath(Sys.getenv("GTR_DATA", unset = "."), mustWork = FALSE)
```

## `docs/`

| File | Content |
| --- | --- |
| `01_pipeline_overview.md` | end-to-end workflow with per-stage inputs/outputs |
| `02_directory_layout.md` | this file |
| `03_file_mapping.md` | old path -> new path for every script |
| `04_data_availability.md` | per-script working directory and input tables |
| `05_software_and_environment.md` | conda environments, R packages, cluster notes |
