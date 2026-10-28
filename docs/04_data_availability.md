# Data availability and input files

The repository contains **code only**. Every table the scripts read is
listed below so that the analyses can be re-run once the data are in place.

Two things determine where a script looks for its inputs:

* `GTR_DATA` — the root of your data directory (default: the original analysis root;
  override with `Sys.setenv(GTR_DATA = "...")`). Each script `setwd()`s into the
  sub-directory listed in the *Working directory* column, so that sub-directory must
  exist and hold the tables named in the *Reads* column.
* `GTR_GW` — only used by the cross-study comparison and by
  `analyses/04_defense_islands/03_nondG_functional_annotation.R`; it points to the
  companion aquifer/groundwater data directory.

The directory names below are the ones used during the study, i.e. the data directory
mirrors the working copy:

```
<GTR_DATA>/
├── 20gaoyuan/
│   ├── 01MAG_info/                 MAG metadata, GTDB-Tk, CheckM quality, sample factors
│   ├── 10Defensefinder/            DefenseFinder output, defense-system definitions
│   ├── 11_5455mag_tree/            bacterial MAG tree + iTOL annotation inputs
│   ├── 11_42mag_ar_tree/           archaeal MAG tree + iTOL annotation inputs
│   ├── 12_5455MAG_MGEs/            geNomad / ISEScan / ICEfinder / integron tables
│   ├── 13_special_MGE/             representative MGE GFFs
│   ├── 14_DIs/                     defense-island summaries and GFFs
│   ├── 15_antidefense/             representative plasmid/phage GFFs, Acr sequences
│   ├── 16DefensePredictor/, 17DI_nonDG/  DefensePredictor output, eggNOG annotations
│   ├── 20ARGMRGVFG/                RGI/CARD, DIAMOND, Macrel, antiSMASH tables
│   ├── 30allphage_antidefensome/   viral catalogue + phage DefenseFinder output
│   ├── 40plateau_groundwater_comparision/  aquifer comparison tables
│   └── Supplementary_tables/       supplementary-table inputs
└── 93_NC_comments/
    ├── 01_diversity_comparison/    abundance tables (丰度/Bins, 丰度/Virus), sample factors
    ├── 02_VHR_prophage_variation/  VHR table (Figure 2 source data)
    ├── 03_DIs_NCBI_RefSeq/         NCBI RefSeq genome list + DI summary
    ├── 04_environmental_factors/   environmental factor table
    └── 05_DS_relative_abundance/   MAG/phage abundance tables
```

## Where the data come from

| Data | Source |
| --- | --- |
| Metagenomic reads | NCBI SRA, accession given in the manuscript |
| MAG assemblies (5,467 bins; 5,455 high-quality) | this study; provided as Supplementary Data |
| GTDB-Tk, CheckM, Prodigal/Prokka output | regenerate with `pipelines/01_genome_annotation/` |
| geNomad / ISEScan / ICEfinder2 / IntegronFinder output | regenerate with `pipelines/02_mobile_genetic_elements/` |
| DefenseFinder / DefensePredictor output | regenerate with `pipelines/03_defense_systems/` |
| Viral catalogue, CheckV, BACPHLIP, spacer BLAST, Acr tree | regenerate with `pipelines/04_viral_contigs/` |
| ARG / MRG / VFG / Macrel / antiSMASH / OGT tables | regenerate with `pipelines/05_adaptive_traits/` |
| Defense-island summaries | regenerate with `pipelines/06_table_postprocessing/` |
| Abundance tables (`all_abundance_*.tsv`) | read mapping of the metagenomic reads to MAGs/phages |
| DefenseFinder system definitions (`df_List_Systems_lpw.csv`) | distributed with DefenseFinder 2.0.0 |
| `Acr_protein.xlsx` (339 validated Acr proteins) | PaCRISPR database |

## Per-script inputs


### `analyses/01_MAG_catalog/01_sample_and_MAG_metadata.R`

**Sample and MAG metadata assembly** — Split the MAG identifiers into region / site / bin and merge them with the sample metadata table.

* Working directory: `20gaoyuan/01MAG_info`
* Reads:
    * `mags_list.txt`
    * `factor_all.csv`
    * `info_corrected.txt`
* Writes: `in-memory metadata table`

### `analyses/01_MAG_catalog/02_sampling_map.R`

**Sampling map** — Maps of the sampling sites coloured by ecosystem and by region, sample elevation vs longitude, and the MAG quality panel.

* Working directory: `20gaoyuan/01MAG_info`
* Reads:
    * `MAG5455_metadata_new.tsv`
    * `5455MAG_metadata.csv`
    * `20gaoyuan/01MAG_info/5455MAG_metadata.csv`
    * `factor_all.csv`
    * `quality_report_90_5_10kb.tsv`
* Writes: `2sample_map_ecosystem.pdf`, `2sample_map_Region.pdf`, `2sample_map_medium*.pdf`, `genome_ecosystem.pdf`

### `analyses/01_MAG_catalog/03_MAG_taxonomy.R`

**MAG taxonomy** — Parse the GTDB-Tk bac120/ar53 classifications into domain-to-species columns and plot the taxonomic composition (sankey).

* Working directory: `20gaoyuan/01MAG_info`
* Reads:
    * `gtdbtk.bac120.summary.tsv`
    * `gtdbtk.ar53.summary.tsv`
* Writes: `3_MAG_tax_sankey.pdf`

### `analyses/01_MAG_catalog/04_MAG_tree_annotation_bacteria.R`

**Bacterial MAG tree annotation** — Build iTOL annotation files for the bacterial MAG phylogeny: region/ecosystem/medium rings plus defense-gene, defense-system, DI and MGE density heatmaps.

* Working directory: `20gaoyuan/11_5455mag_tree`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `20gaoyuan/01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
    * `ds_density_tax.csv`
    * `5455MAG_defense_islands_summary.csv`
    * `20gaoyuan/12_5455MAG_MGEs/genome_MGE_number.tsv`
* Writes: `iTOL annotation files (hub)`

### `analyses/01_MAG_catalog/05_MAG_tree_annotation_archaea.R`

**Archaeal MAG tree annotation** — Same iTOL annotation workflow as step 04, for the archaeal MAG tree.

* Working directory: `20gaoyuan/11_42mag_ar_tree`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `20gaoyuan/01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
    * `ds_density_tax.csv`
    * `5455MAG_defense_islands_summary.csv`
    * `20gaoyuan/12_5455MAG_MGEs/genome_MGE_number.tsv`
* Writes: `iTOL annotation files (hub)`

### `analyses/01_MAG_catalog/06_supplementary_MAG_table.R`

**Supplementary MAG table** — Merge MAG metadata with the CheckM quality report and the GTDB taxonomy to build the per-MAG supplementary table.

* Working directory: `20gaoyuan/Supplementary_tables`
* Reads:
    * `MAG5455_metadata.tsv`
    * `quality_report_90_5_10kb.tsv`
    * `gtdbtk.bac120.summary.tsv`
    * `gtdbtk.ar53.summary.tsv`
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `5455MAG_defense_islands_summary_no_anti.csv`
* Writes: `supplementary MAG table`

### `analyses/02_defensome_landscape/01_defense_system_overview.R`

**Defensome overview** — Build the defense-system count matrix from the DefenseFinder output and derive prevalence, per-MAG DS number/diversity and treemap abundance panels.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
    * `merged_5467MAG_defense_finder_genes.tsv`
* Writes: `defense-system matrices and overview figures`

### `analyses/02_defensome_landscape/02_region_comparison.R`

**Region comparison (glacier / plateau / plain)** — Compare DS number, density and diversity between the glacier, plateau and plain regions, with Wilcoxon/Mann-Whitney tests.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
* Writes: `2_genome_Region.pdf`, `2_two_barplot.pdf`, `2defensefinder_res_Glacier_Plateau_Plain.RData`

### `analyses/02_defensome_landscape/03_medium_comparison.R`

**Medium / lifestyle comparison** — Compare DS content across snow, ice, cryoconite, water, particle and sediment, and between planktonic and benthic lifestyles.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `20gaoyuan/01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
* Writes: `3genome_medium.pdf`, `2_genome_Region.pdf`

### `analyses/02_defensome_landscape/04_density_and_genome_size.R`

**DS density vs genome size and elevation** — Relate DS number and DS density to genome size and sample elevation; builds the density comparison panels.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
* Writes: `4defensefinder_res_density_compare.RData`

### `analyses/02_defensome_landscape/05_defense_system_composition_pcoa.R`

**Defense-system composition (PCoA / beta diversity)** — Ordination of MAGs by defense-system composition, alpha/beta diversity of the defensome and the DS-frequency heatmap.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
    * `pcoadata.csv`
    * `gtdbtk.bac120.summary.tsv`
    * `gtdbtk.ar53.summary.tsv`
    * `pcoadata_noRMSoFIC.csv`
* Writes: `5_heatmap_DSfrequency.pdf`, `5_MAG_DS_pcoa.RData`

### `analyses/02_defensome_landscape/06_taxonomic_composition.R`

**Taxonomic composition of MAGs** — Aggregate the MAG taxonomy into phylum-level composition across regions and ecosystems.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `gtdbtk.bac120.summary.tsv`
    * `gtdbtk.ar53.summary.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
* Writes: `7_tax_legend.pdf`

### `analyses/02_defensome_landscape/07_phylum_specialization.R`

**Defensome patterns across phyla** — Per-phylum DS number, density and diversity, plus defense-system frequency for Cyanobacteriota, Myxococcota and Gemmatimonadota.

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
    * `gtdbtk.bac120.summary.tsv`
    * `gtdbtk.ar53.summary.tsv`
* Writes: `6_phyla_dss_frequency.pdf`, `6_Cyanobacteriota_order_dss_frequency.pdf`, `6_Myxococcota_class_DS_frequency.pdf`, `6_Gemmatimonadota_phylum_DS_frequency.pdf`

### `analyses/02_defensome_landscape/08_season_comparison.R`

**Season comparison** — Compare DS number and density between sampling seasons (companion analysis to the region comparison).

* Working directory: `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
* Writes: `8genome_season.pdf`, `8defensefinder_res_season_compare.RData`

### `analyses/03_mobile_genetic_elements/01_MGE_cargo_and_defense_density.R`

**MGEs and the defensome they carry** — Count prophages, plasmids, ICE/IME and IS per MAG, compare DS density on chromosomes vs MGEs, and compute the O/E localisation ratios.

* Working directory: `20gaoyuan/12_5455MAG_MGEs`, `20gaoyuan/10Defensefinder`
* Reads:
    * `merged_5467MAG_virus_summary.tsv`
    * `merged_5467MAG_plasmid_summary.tsv`
    * `processed_ice_sequences_ID.txt`
    * `5455MAG_filtered_completeIS.tsv`
    * `all5455MAG_contig.integrons`
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `DG_gff.gff`
    * `DG_gff_DF.gff`
    * `MAG5455_metadata.tsv`
    * `quality_report_90_5_10kb.tsv`
    * `genome_MGE_length.tsv`
    * `.RData`
    * `genome_MGE_number.tsv`
* Writes: `1_p_MGE_number.pdf`, `1_mge_ds_density_comparison*.pdf`, `1_MGEnum_*.pdf`, `1_DSnum_glacier_plateau_plain_comprison.pdf`, `1_MGE_combination.RData`

### `analyses/03_mobile_genetic_elements/02_MGE_cargo_antidefense.R`

**Anti-defense genes carried by MGEs** — Count anti-defense genes per MGE type, their density and their O/E localisation ratios.

* Working directory: `20gaoyuan/12_5455MAG_MGEs`
* Reads:
    * `merged_5467MAG_virus_summary.tsv`
    * `merged_5467MAG_plasmid_summary.tsv`
    * `processed_ice_sequences_ID.txt`
    * `5455MAG_filtered_completeIS.tsv`
    * `all5455MAG_contig.integrons`
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `DG_gff.gff`
    * `DG_gff_DF.gff`
    * `MAG5455_metadata.tsv`
    * `quality_report_90_5_10kb.tsv`
    * `genome_MGE_length.tsv`
* Writes: `2_mge_antids_density_comparison.pdf`, `2_MGE_antidefense.RData`

### `analyses/03_mobile_genetic_elements/03_representative_MGEs.R`

**Representative MGEs** — Draw the representative ICE/IS elements adjacent to defense systems (type I, II and IIG R-M) and the CRISPR-Cas9 region of the special MGE.

* Working directory: `20gaoyuan/13_special_MGE`, `20gaoyuan/13_special_MGE/1ICE`, `20gaoyuan/13_special_MGE/2IS`, `20gaoyuan/13_special_MGE/3prophage`, `20gaoyuan/13_special_MGE/4CRISPR-Cas9`
* Reads:
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `DG_gff.gff`
    * `MAG5455_metadata.tsv`
    * `CJ_14-30WL_WA-bin.46-17.gff`
    * `CJ_14-30WL_WA-bin.46_contig_17_ICE2_gene.gff`
    * `Glacier-1611DDS-bin.12-3.gff`
    * `Glacier-1611DDS-bin.12_contig_3_ICE2_gene.gff`
    * `DG_gff_DF.gff`
    * `IS_DS_filtered_mix.tsv`
    * `selected_IS_DS_contig.gff`
    * `selected_prophage.tsv`
    * `selected_prophage_DS_contig.gff`
    * `final_prophage_marker.tsv`
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `sequence_lengths.tsv`
    * `MM_igenas9b.emapper.annotations_2.tsv`
* Writes: `p_final_contig_2ICE.pdf`, `1_p_cas9gene_number.pdf`, `1_p_protein_annotation.pdf`

### `analyses/04_defense_islands/01_defense_island_characterization.R`

**Defense-island characterization** — Number of DIs and DG density across the major lineages, host contribution by region/ecosystem/medium, DS composition of DIs and the DI O/E heatmap.

* Working directory: `20gaoyuan/14_DIs`
* Reads:
    * `5455MAG_defense_islands_summary_no_anti.csv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `../01MAG_info/MAG5455_metadata_new.tsv`
    * `../01MAG_info/quality_report_90_5_10kb.tsv`
    * `../10Defensefinder/gtdbtk.bac120.summary.tsv`
    * `../10Defensefinder/gtdbtk.ar53.summary.tsv`
    * `../10Defensefinder/merged_5467MAG_defense_finder_systems.tsv`
    * `5455MAG_defense_islands_summary_with_anti.csv`
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `../01MAG_info/MAG5455_metadata.tsv`
    * `selected_DI_contig.gff`
* Writes: `1_summary_stats_DI_n_normalized.pdf`, `1_DS_volcano_fisher.pdf`, `1_DIs.RData`

### `analyses/04_defense_islands/02_antidefense_in_defense_islands.R`

**Anti-defense genes inside defense islands** — Map anti-defense genes onto the defense islands and describe the representative DIs that carry both defense and anti-defense systems.

* Working directory: `20gaoyuan/15_antidefense/4antidefense_colocalization`
* Reads:
    * `3contig.id_anti.tsv`
    * `virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv`
    * `3contig.id.gff`
* Writes: `anti-defense co-localisation tables`

### `analyses/04_defense_islands/03_nondG_functional_annotation.R`

**Non-defense genes in defense islands** — Functional annotation (eggNOG/COG/KEGG/PFAM) of the non-defense genes inside defense islands and the DefensePredictor share of them.

* Working directory: `20gaoyuan/17DI_nonDG` and `GTR_GW`
* Reads:
    * `defensepredictor_genecount.csv`
    * `DF_DG_position.tsv`
    * `../01MAG_info/5455MAG_metadata.csv`
    * `5466mag_defensepredictor_processed_filtered.csv`
    * `DI_nonDG_emapper.annotations.tsv`
    * `filtered_DI_genes.gff`
    * `5455MAG_defense_islands_summary_no_anti.csv`
* Writes: `1_DP_DF_compare.pdf`, `1_DI_nonDG_eggnog.pdf`, `1_nonDG_DP_proportion.pdf`, `1DF_DP.RData`

### `analyses/04_defense_islands/04_defensefinder_vs_defensepredictor.R`

**DefenseFinder vs DefensePredictor** — Compare the defense genes reported by DefenseFinder and DefensePredictor across the 5,455 MAGs (counts, density, overlap).

* Working directory: -
* Reads:
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `defensepredictor_genecount.csv`
    * `DF_DG_position.tsv`
    * `5466mag_defensepredictor_processed_filtered.csv`
* Writes: `FigS6B_C_density.pdf`, `1_DP_DF_compare_3.pdf`

### `analyses/05_antidefensome/01_phage_antidefensome.R`

**Phage antidefensome** — Describe the viral catalog (completeness/length, taxonomy) and the defense and anti-defense systems encoded by 59,997 phage genomes.

* Working directory: `20gaoyuan/30allphage_antidefensome`
* Reads:
    * `virus_select_all_medium-10kb_high_complete_defense_finder_systems.tsv`
    * `virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv`
    * `virus_select_all_medium-10kb_high_complete_virus_summary_genomad.tsv`
    * `quality_summary_checkv.tsv`
    * `info_corrected.txt`
* Writes: `1_p_phage_quality.pdf`, `1_phage_tax_sankey.pdf`, `1_phage_tax_legend.pdf`, `1_allphage_antidefense.RData`

### `analyses/05_antidefensome/02_special_plasmid.R`

**Representative plasmid arms race** — Gene-level view of a plasmid encoding anti-RM genes (ArdC/ArdK/KlcA) targeted by the host type II-C CRISPR-Cas system.

* Working directory: `20gaoyuan/15_antidefense/1special_plasmid`
* Reads:
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `selected_plasmid_contig.gff`
    * `selected_plasmidhost_contig.gff`
* Writes: `special_plasmid_bind.tsv`, `1_special_plasmid.RData`

### `analyses/05_antidefensome/03_special_phage.R`

**Representative phage arms race** — Gene-level view of a phage encoding AcrVA3 against the host type I-F CRISPR-Cas system (linked by a CRISPR spacer match).

* Working directory: `20gaoyuan/15_antidefense/2special_phage`
* Reads:
    * `merged_5467MAG_defense_finder_genes.tsv`
    * `virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv`
    * `selected_phage_contig.gff`
    * `selected_phagehost_contig.gff`
* Writes: `special_phage_bind.tsv`, `1_special_phage.RData`

### `analyses/05_antidefensome/04_Acr_phylogeny.R`

**Anti-CRISPR phylogeny annotation** — Annotate the maximum-likelihood Acr tree (44 predicted + 339 validated Acr proteins) with iTOL colour strips.

* Working directory: `20gaoyuan/15_antidefense/3Acr_tree`
* Reads:
    * `./Acr_protein.xlsx`
    * `44_viral_acr.tsv`
    * `spacer_targeted_prophagewithAnti_CRISPR.txt`
    * `./Acr_pAcr_Sequence_Files/Acr_protein.xlsx`
    * `prophage_164Anti_CRISPR.txt`
* Writes: `iTOL annotation files (hub)`

### `analyses/06_adaptive_traits/01_ARG_MRG_VFG_traits.R`

**Adaptive traits of MAGs** — Combine ARG, MRG, VFG and POP-degradation gene counts with MAG metadata to relate mobile-element load to adaptive traits.

* Working directory: `20gaoyuan/20ARGMRGVFG`
* Reads:
    * `rgiargs_2.txt`
    * `diamond_MRG_out3_type.csv`
    * `diamond_VFG_out3_type.csv`
    * `diamond_POPDeg_out3_type.csv`
    * `MAG5455_metadata.tsv`
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `macrel.out.prediction`
    * `antismash_res_type.tsv`
* Writes: `per-MAG adaptive-trait tables`

### `analyses/07_cross_study_comparison/01_defense_system_comparison.R`

**Rivers vs aquifers: defense systems** — Compare the river/glacier defensome with the aquifer defensome of a companion study.

* Working directory: `20gaoyuan/40plateau_groundwater_comparision`
* Reads:
    * `merged_5467MAG_defense_finder_systems.tsv`
    * `20gaoyuan/01MAG_info/5455MAG_metadata.csv`
    * `quality_report_90_5_10kb.tsv`
    * `1626mag_ds_counts_matrix.csv`
    * `1626MAGs_90_5_N50100kb_score70.txt`
* Writes: `1_river_aquifer_p.pdf`

### `analyses/07_cross_study_comparison/02_antidefense_comparison.R`

**Rivers vs aquifers: phage anti-defense** — Compare the anti-defense gene content of river phages with groundwater phages.

* Working directory: `20gaoyuan/40plateau_groundwater_comparision`
* Reads:
    * `virus_select_all_medium-10kb_high_complete_defense_finder_systems.tsv`
    * `virus_select_all_medium-10kb_high_complete_defense_finder_genes.tsv`
    * `virus_select_all_medium-10kb_high_complete_virus_summary_genomad.tsv`
    * `quality_summary_checkv.tsv`
    * `GWVC.vOTUs_defense_finder_systems.tsv`
    * `GWVC.vOTUs.information.txt`
* Writes: `2_river_phage_aquifer_comparision.pdf`

### `analyses/08_revision_nc/01_MAG_and_phage_macrodiversity.R`

**MAG and phage macro-diversity** — Alpha/beta diversity of MAGs and phages along the glacier-to-river gradient (reviewer revision).

* Working directory: `93_NC_comments/01_diversity_comparison`
* Reads:
    * `./丰度/Bins/all_abundance_rpkm_bins.tsv`
    * `./info_corrected.txt`
    * `./丰度/Virus/all_abundance_rpkm.tsv`
    * `./丰度/Bins/all_abundance_tpm_bin.tsv`
* Writes: `1_MAG_Phage_macrodiversity.RData`, `1_MAG_Phage_macrodiversity_new.RData`, `NMDS/PCoA tables`

### `analyses/08_revision_nc/02_VHR_and_prophage_abundance.R`

**VHR and prophage variation** — Virus-to-host ratio (VHR) and prophage abundance across regions/mediums; regenerates the main Figure 2.

* Working directory: `93_NC_comments/02_VHR_prophage_variation`
* Reads:
    * `./VHR_from_NC_fig2c.xlsx`
    * `../01_diversity_comparison/info_corrected.txt`
    * `all_abundance_tpm.tsv`
    * `lifecycle_virus.csv`
    * `20gaoyuan/10Defensefinder/04_density_and_genome_size.RData`
    * `20gaoyuan/10Defensefinder/05_defense_system_composition_pcoa.RData`
    * `20gaoyuan/12_5455MAG_MGEs/01_MGE_cargo_and_defense_density.RData`
* Writes: `Figure2.pdf`, `1_VHR_prophage_abundance.RData`, `stat_tests_final.tsv`

### `analyses/08_revision_nc/03_defense_islands_in_NCBI_RefSeq.R`

**Defense islands in NCBI RefSeq genomes** — Repeat the defense-island analysis on 22,942 NCBI RefSeq genomes to test whether the observed DI patterns generalize.

* Working directory: `93_NC_comments/03_DIs_NCBI_RefSeq`
* Reads:
    * `22942genome_defense_islands_summary.csv`
    * `NCBI_refseq_22143_bacterial_genomes_20260305_withtax.tsv`
    * `NCBI_refseq_799_archaeal_genomes_20260305_withtax.tsv`
* Writes: `1_summary_stats_DI_n_normalized.pdf`, `1_DIs_NCBI_refseq_genome.RData`

### `analyses/08_revision_nc/04_environmental_factor_effects.R`

**Environmental drivers of the defensome** — Relate DS number/density and anti-defense content to the measured environmental factors.

* Working directory: `93_NC_comments/04_environmental_factors`
* Reads:
    * `info_corrected.txt`
    * `factor_all.csv`
    * `sourceDataFig1b_fromNC.tsv`
    * `../01_diversity_comparison/丰度/Bins/all_abundance_rpkm_bins.tsv`
    * `../01_diversity_comparison/丰度/Virus/all_abundance_rpkm.tsv`
    * `8_MAG_DS_type_num.tsv`
    * `MAG5455_metadata_new.tsv`
* Writes: `1_environmental_factor_effects.RData`

### `analyses/08_revision_nc/05_defense_system_relative_abundance.R`

**Relative abundance of defense systems** — Weight the defense-system content by MAG relative abundance (TPM) to obtain community-level defensome profiles.

* Working directory: `93_NC_comments/05_DS_relative_abundance`
* Reads:
    * `EF_info.tsv`
    * `all_abundance_tpm_bin.tsv`
    * `merged_4903MAGs_defense_finder_genes.tsv`
    * `merged_4903MAGs_defense_finder_systems.tsv`
* Writes: `1_DS_relative_abundance.RData`

