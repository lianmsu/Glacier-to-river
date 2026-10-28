# Glacier-to-river

Code accompanying **"Microbial evolutionary arms races across glacier-to-river gradients"**
(Li P., Li Y., Ni J.).

The repository contains every script used to go from public metagenomic data to the
figures of the manuscript: the upstream bioinformatic pipelines that build the MAG
catalogue, the mobile genetic element and defense/anti-defense annotation, and the R
analyses that produce the statistics and the figures.

> 中文版说明见 **[README_zh.md](README_zh.md)**。

---

## What is in here

```
Glacier-to-river/
├── pipelines/          # upstream bioinformatics (SLURM + conda, shell/python)
│   ├── 01_genome_annotation/         MAG renaming, Prokka/Prodigal, GTDB-Tk
│   ├── 02_mobile_genetic_elements/   geNomad, ISEScan, ICEfinder2, IntegronFinder, minced
│   ├── 03_defense_systems/           DefenseFinder 1/2, DefensePredictor
│   ├── 04_viral_contigs/             CheckV, BACPHLIP, BLAST, FastTree
│   ├── 05_adaptive_traits/           RGI/CARD, DIAMOND, Macrel, antiSMASH, OGT
│   └── 06_table_postprocessing/      defense-island tables (python + awk)
├── analyses/           # downstream statistics and figures (R)
│   ├── 01_MAG_catalog/               Fig. 1, Fig. S1, Fig. S7
│   ├── 02_defensome_landscape/       Fig. 2, Fig. 3, Figs. S2-S4
│   ├── 03_mobile_genetic_elements/   Fig. 4, Figs. S5, S6, S10
│   ├── 04_defense_islands/           Fig. 5, Figs. S11-S13
│   ├── 05_antidefensome/             Fig. 6, Figs. S14-S16
│   ├── 06_adaptive_traits/           ARG / MRG / VFG / POP-degradation traits
│   ├── 07_cross_study_comparison/    rivers vs aquifers (companion study)
│   └── 08_revision_nc/               analyses added during peer review
└── docs/               # workflow, layout, file mapping, data and software notes
```

Both layers are numbered in execution order; each file also carries a header with its
purpose, the manuscript figure it produces, its main output files and its original path
in the working copy.

## Documentation

| Document | Content |
| --- | --- |
| [docs/01_pipeline_overview.md](docs/01_pipeline_overview.md) | End-to-end workflow, stage by stage, with inputs and outputs |
| [docs/02_directory_layout.md](docs/02_directory_layout.md) | What lives where, and the naming conventions |
| [docs/03_file_mapping.md](docs/03_file_mapping.md) | Every script: original path -> new path, previous repo path, purpose |
| [docs/04_data_availability.md](docs/04_data_availability.md) | Input tables each analysis expects and where to obtain them |
| [docs/05_software_and_environment.md](docs/05_software_and_environment.md) | Cluster setup, conda environments, R packages |

## Reproducing the analyses

**1. Upstream pipelines.** The scripts under `pipelines/` are the SLURM job scripts used
for the production runs, wrapped per software (`#SBATCH` headers included). They are
reproduced as executed, except that the site-specific absolute paths were replaced by
variables: set `LUSTRE_HOME`, `HOME_ROOT` and `DB_ROOT` (see the
`# Site-specific settings` block of each script) and adapt the queue names, then submit in
the order given in [docs/01](docs/01_pipeline_overview.md), using the conda environments
listed in [docs/05](docs/05_software_and_environment.md).

**2. R analyses.** The scripts under `analyses/` require R >= 4.2 and the packages listed
in [docs/05](docs/05_software_and_environment.md). Input tables are read relative to
`GTR_DATA`, which defaults to the current working directory and can be overridden without
editing any script:

```r
Sys.setenv(GTR_DATA = "/path/to/plateau_defensome")
source("analyses/02_defensome_landscape/02_region_comparison.R")
```

Scripts are numbered within each module and are meant to be run in that order: some of
them hand over intermediate state through `.RData` files (e.g.
`analyses/02_defensome_landscape/05_defense_system_composition_pcoa.RData`), which are
written and re-read by later scripts.

## Data availability

The MAGs, the viral catalogue and the defense/anti-defense annotation tables are provided
as Supplementary Data with the manuscript. Sources and file names of every required input
table are listed in [docs/04_data_availability.md](docs/04_data_availability.md). The raw
metagenomic reads are available from the NCBI SRA under the accession given in the
manuscript.

## Note on the previous repository layout

An earlier version of this repository was a flat list of scripts (`*.sh`, `*.py`,
`R_code/*.R`). It has been reorganized into the layered structure above. No analysis was
dropped: the scripts were renamed and made path-portable, 10 analyses that were missing
from the old version were added (including the five peer-review analyses), and only stale
duplicates and empty drafts were removed (six shell scripts existed twice, plus four
obsolete or empty files). [docs/03_file_mapping.md](docs/03_file_mapping.md) maps every
old file to its new location.

No absolute path, user home directory, account name or Windows analysis folder remains in
the repository: all site-specific locations are handled by the `LUSTRE_HOME`,
`HOME_ROOT`, `DB_ROOT`, `GTR_DATA`, `GTR_GW` and `GTR_DI_DIR` environment variables
(see [docs/03](docs/03_file_mapping.md#site-specific-paths-and-how-they-are-hidden)).

## Citation

Li, P., Li, Y. & Ni, J. Microbial evolutionary arms races across glacier-to-river
gradients. (manuscript; citation details will be updated upon publication)

## Contact

Corresponding author: Prof. Jinren Ni, College of Environmental Sciences and Engineering,
Peking University. Correspondence details are given in the manuscript; for questions about
the code, please open an issue in this repository.
