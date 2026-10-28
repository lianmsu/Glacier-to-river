# Workflow overview

The analysis falls into two layers. The **upstream layer** (`pipelines/`) turns MAG
assemblies into annotation tables with SLURM job scripts, one per software. The
**downstream layer** (`analyses/`) turns those tables into statistics and figures with R.
The two layers are connected only through files, so any step can be re-run on its own.

```
                         MAG assemblies (5,467 bins; 5,455 high-quality)
                                        |
   [1] genome annotation ---------------+-------------------------------+
       contig renaming, Prokka/Prodigal, GTDB-Tk                          |
                                        |                                 |
   [2] mobile genetic elements ---------+   [3] defense systems ----------+
       geNomad (phage/plasmid)                 DefenseFinder 2.0          |
       ISEScan (IS)                            DefensePredictor           |
       ICEfinder2 (ICE/IME)                                               |
       IntegronFinder, minced (CRISPR)                                    |
                                        |                                 |
   [4] defense islands -----------------+   [5] viral contigs -----------+
       DI tables from DefenseFinder            CheckV, BACPHLIP           |
       DI / DG / non-DG ORF lists              BLAST spacer -> phage      |
       (python + awk)                          FastTree (Acr phylogeny)   |
                                        |                                 |
   [6] adaptive traits -----------------+                                 |
       RGI/CARD, DIAMOND, Macrel, antiSMASH, OGT                         |
                                        |                                 |
                                        v                                 v
   [7] R analyses: analyses/01_MAG_catalog ... analyses/08_revision_nc
       statistics, figures, supplementary tables
```

## Stage 1 — Genome annotation

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 1.1 | `01_rename_MAG_contigs.sh` | awk | `<MAG-id>-<n>` contig names (basis of every downstream gene ID) |
| 1.2 | `02_prokka_annotate_MAGs.sh` | Prokka (metagenome mode) | per-MAG GFF/FFN/FAA |
| 1.3 | `03_prodigal_MAGs.sh` | Prodigal | ORFs of the concatenated 5,455 MAG contigs |
| 1.4 | `04_prodigal_merged_MAGs.sh` | Prodigal | `merged_5467MAG_ORF.gff/.faa/.fna` |
| 1.5 | `05_gtdbtk_classify.sh` | GTDB-Tk 2.4.0 | `gtdbtk.bac120.summary.tsv`, `gtdbtk.ar53.summary.tsv` |

## Stage 2 — Mobile genetic elements

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 2.1 | `01_genomad.sh` | geNomad | `merged_5467MAG_virus_summary.tsv`, `merged_5467MAG_plasmid_summary.tsv` |
| 2.2 | `02_isescan.sh` | ISEScan | `5455MAG_filtered_completeIS.tsv` |
| 2.3 | `03_icefinder_batch_1.sh`, `04_icefinder_batch_2.sh` | ICEfinder2 | `processed_ice_sequences_ID.txt` (ICE/IME coordinates) |
| 2.4 | `05_integron_finder.sh` | IntegronFinder | integron calls (anti-Thoeris hits in integrons) |
| 2.5 | `06_minced_crispr.sh` | minced | CRISPR arrays and spacers |

## Stage 3 — Defense and anti-defense systems

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 3.1 | `01_defensefinder_v1.sh`, `02_defensefinder_v2.sh` | DefenseFinder (2.0.0) | `merged_5467MAG_defense_finder_systems.tsv`, `..._genes.tsv` |
| 3.2 | `03_defensepredictor_batch_1.sh`, `04_defensepredictor_batch_2.sh` | DefensePredictor | `defensepredictor_genecount.csv` |

`merged_5467MAG_defense_finder_systems.tsv` is the single most reused table of the study:
its `activity == "Defense"` / `"Antidefense"` rows feed every defensome, MGE and
defense-island analysis.

## Stage 4 — Defense islands

The defense islands are derived from the DefenseFinder output by a small python/awk chain
(`pipelines/06_table_postprocessing/`):

```
DIs/*.tsv                    (per-island DefenseFinder output)
   | 01_collect_DI_hits.py     keep Defense rows, standardise contig IDs
   v
5455MAG_defense_islands_summary.csv
   | 02_DI_gene_ids.py         expand records -> ORF IDs   (01all_DIorf.txt)
   | 03_DI_summary.py          per-island summary          (5455MAG_defense_islands.txt)
   | 04_DI_nondG_orf_ids.py    ORF IDs per island
   | 05_DI_nondG_filter.sh     subtract DG IDs, seqkit extraction, GFF filtering
   v
5455MAG_defense_islands_summary_no_anti.csv / _with_anti.csv
04non_DG_orf.faa, filtered_DI_genes.gff
```

## Stage 5 — Viral contigs

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 5.1 | `01_checkv.sh` | CheckV | completeness/contamination, quality tiers |
| 5.2 | `02_bacphlip.sh` | BACPHLIP | virulent vs temperate lifestyle |
| 5.3 | `03_makeblastdb.sh` | BLAST+ | phage and MAG-plasmid databases |
| 5.4 | `04_blastn_spacers_vs_phage_plasmid.sh` | BLAST+ | spacer/phage and spacer/plasmid hits (<=1 mismatch, >=95% coverage) |
| 5.5 | `05_fasttree_acr_phylogeny.sh` | muscle, trimAl, FastTree | `383acr_aln.trim.afa.treefile` |

The spacer hits of 5.4 establish the historical phage/plasmid-host links used in Fig. 6.

## Stage 6 — Adaptive traits

| # | Script | Software | Produces |
| --- | --- | --- | --- |
| 6.1 | `01_rgi_card.sh` | RGI/CARD | `rgiargs_2.txt` |
| 6.2 | `02_diamond_functional_screening.sh` | DIAMOND | ARG / MRG / mibPOP / VFG hits |
| 6.3 | `03_diamond_result_pretreatment.sh` | awk | per-MAG gene-count tables (`diamond_*_out3_type.csv`) |
| 6.4 | `04_macrel.sh` | Macrel | antimicrobial peptides |
| 6.5 | `05_antismash_batch_1.sh` | antiSMASH | biosynthetic gene clusters |
| 6.6 | `06_ogt_prediction.sh` | OGT_prediction | predicted growth temperature |

## Stage 7 — R analyses

Run each module in numbered order. Modules follow the manuscript narrative:

| Module | Manuscript | Key input tables |
| --- | --- | --- |
| `01_MAG_catalog` | Fig. 1; Figs. S1, S7 | GTDB-Tk summaries, MAG metadata, CheckM quality report |
| `02_defensome_landscape` | Fig. 2, Fig. 3; Figs. S2-S4 | DefenseFinder systems table, MAG metadata |
| `03_mobile_genetic_elements` | Fig. 4; Figs. S5, S6, S10 | geNomad, ISEScan, ICEfinder2 tables |
| `04_defense_islands` | Fig. 5; Figs. S11-S13 | DI summaries, DefenseFinder genes, DefensePredictor counts |
| `05_antidefensome` | Fig. 6; Figs. S14-S16 | viral catalogue, DefenseFinder (phage), BLAST spacer hits, Acr tree |
| `06_adaptive_traits` | - | RGI/CARD, DIAMOND tables |
| `07_cross_study_comparison` | - | companion groundwater/aquifer tables |
| `08_revision_nc` | peer-review revision | abundance tables, environmental factors, NCBI RefSeq DIs |

### Inter-script dependencies (`.RData` hand-over)

Some modules pass intermediate objects through `.RData` files. Run them in order; the
files are written next to the input data, in the directory the script `setwd()`s into:

| Produced by | Read by |
| --- | --- |
| `02_defensome_landscape/04_density_and_genome_size.R` | `08_revision_nc/02_VHR_and_prophage_abundance.R` |
| `02_defensome_landscape/05_defense_system_composition_pcoa.R` | `08_revision_nc/02_VHR_and_prophage_abundance.R` |
| `03_mobile_genetic_elements/01_MGE_cargo_and_defense_density.R` | `08_revision_nc/02_VHR_and_prophage_abundance.R` |

`08_revision_nc/02_VHR_and_prophage_abundance.R` assembles the final Fig. 2 from those
three intermediate states, so it must be run after them.
