# File mapping

Every script of this repository, its original location in the working copy, and its
location in the previous (flat) version of this GitHub repository.

`-` means the script was not part of the previous repository version.

The scripts were **moved, renamed and made path-portable**; the statistical and
bioinformatic logic is unchanged. The only code edits in the R scripts are:

1. the hard-coded `setwd("D:/...")` and absolute Windows paths were replaced by
   `file.path(GTR_DATA, "<relative path>")`, where `GTR_DATA` defaults to the
   original analysis root and can be overridden with an environment variable;
2. the `.RData` artefacts written/read between scripts were renamed together
   with the scripts (`1_MGE_combination.RData` -> `01_MGE_cargo_and_defense_density.RData`, ...);
3. a short header (purpose, manuscript figure, outputs, original path) was added.


## Upstream pipelines

| New path | Original working-copy path | Previous repo path | Purpose |
| --- | --- | --- | --- |
| `pipelines/01_genome_annotation/01_rename_MAG_contigs.sh` | `20gaoyuan/submission/NC/code/MAGrename.sh` | `MAGrename.sh` | Rename every contig of each MAG assembly to `<MAG-id>-<n>`; this naming is the basis of all downstream gene IDs. |
| `pipelines/01_genome_annotation/02_prokka_annotate_MAGs.sh` | `20gaoyuan/submission/code_availability/01prokka_lpw.sh` | `01prokka_lpw.sh` | Prokka (metagenome mode) annotation of the 5,467 MAG assemblies, one output directory per MAG. |
| `pipelines/01_genome_annotation/03_prodigal_MAGs.sh` | `20gaoyuan/submission/NC/code/01mags_prodigal.sh` | `01mags_prodigal.sh` | Prodigal gene prediction on the concatenated contigs of the 5,455 high-quality MAGs (protein set reused by the defense-island analyses). |
| `pipelines/01_genome_annotation/04_prodigal_merged_MAGs.sh` | `20gaoyuan/submission/code_availability/01prodigal.sh` | `01prodigal.sh` | Prodigal gene prediction on the merged `merged_5467MAG.fa` contig set. |
| `pipelines/01_genome_annotation/05_gtdbtk_classify.sh` | `20gaoyuan/submission/NC/code/01GTDB_classify_1.sh` | `01GTDB_classify_1.sh` | GTDB-Tk `classify_wf` taxonomic assignment of the 5,455 MAGs (bac120 + ar53 summaries). |
| `pipelines/02_mobile_genetic_elements/01_genomad.sh` | `20gaoyuan/submission/code_availability/01genomad_lpw.sh` | `01genomad_lpw.sh` | geNomad end-to-end run on the merged MAG contigs: virus/provirus and plasmid detection. |
| `pipelines/02_mobile_genetic_elements/02_isescan.sh` | `20gaoyuan/submission/NC/code/ISEScan_5455MAGs.sh` | `ISEScan_5455MAGs.sh` | ISEScan detection of insertion sequences (IS) in the 5,455 MAG assemblies. |
| `pipelines/02_mobile_genetic_elements/03_icefinder_batch_1.sh` | `20gaoyuan/submission/NC/code/icefinder_1000.sh` | `icefinder_1000.sh` | ICEfinder2 detection of integrative conjugative elements (ICE/IME), MAG batch 1. |
| `pipelines/02_mobile_genetic_elements/04_icefinder_batch_2.sh` | `20gaoyuan/submission/code_availability/icefinder_11000.sh` | `icefinder_11000.sh` | ICEfinder2 detection of integrative conjugative elements (ICE/IME), MAG batch 2. |
| `pipelines/02_mobile_genetic_elements/05_integron_finder.sh` | `20gaoyuan/submission/NC/code/01integron_finder.sh` | `01integron_finder.sh` | IntegronFinder detection of integrons in the 5,455 MAG contigs. |
| `pipelines/02_mobile_genetic_elements/06_minced_crispr.sh` | `20gaoyuan/submission/code_availability/1minced.sh` | `1minced.sh` | minced detection of CRISPR arrays and spacer sequences (spacers are later mapped to phages/plasmids). |
| `pipelines/03_defense_systems/01_defensefinder_v1.sh` | `20gaoyuan/submission/code_availability/01DefenseFinder.sh` | `01DefenseFinder.sh` | DefenseFinder run on all MAG proteins (first pass, `res1`). |
| `pipelines/03_defense_systems/02_defensefinder_v2.sh` | `20gaoyuan/submission/NC/code/01DefenseFinder2.0.0.sh` | `01DefenseFinder2.0.0.sh` | DefenseFinder 2.0.0 run on the merged MAG proteins. Seeds the defense-system, defense-island and anti-defense analyses. |
| `pipelines/03_defense_systems/03_defensepredictor_batch_1.sh` | `20gaoyuan/submission/code_availability/01DefensePredictor_lpw.sh` | `01DefensePredictor_lpw.sh` | DefensePredictor (deep-learning defense gene finder), MAG batch 1 (1-1000). |
| `pipelines/03_defense_systems/04_defensepredictor_batch_2.sh` | `20gaoyuan/submission/NC/code/01DefensePredictor_lpw_1001_2000.sh` | `01DefensePredictor_lpw_1001_2000.sh` | DefensePredictor, MAG batch 2 (1001-2000). |
| `pipelines/04_viral_contigs/01_checkv.sh` | `20gaoyuan/submission/code_availability/01checkv_lpw.sh` | `01checkv_lpw.sh` | CheckV quality assessment of the viral contigs (completeness, contamination, provirus removal). |
| `pipelines/04_viral_contigs/02_bacphlip.sh` | `20gaoyuan/submission/code_availability/01bacphlip_lpw.sh` | `01bacphlip_lpw.sh` | BACPHLIP prediction of virulent vs temperate (lysogenic) lifestyle for each viral contig. |
| `pipelines/04_viral_contigs/03_makeblastdb.sh` | `20gaoyuan/submission/NC/code/2makeblastdb.sh` | `2makeblastdb.sh` | Build BLAST databases for the whole phage catalog and for the MAG-encoded plasmids. |
| `pipelines/04_viral_contigs/04_blastn_spacers_vs_phage_plasmid.sh` | `20gaoyuan/submission/code_availability/3blastn_allphage_plasmid.sh` | `3blastn_allphage_plasmid.sh` | Map CRISPR spacers against the phage and plasmid databases, keeping <=1 mismatch and >=95% coverage hits (historical phage/plasmid-host links). |
| `pipelines/04_viral_contigs/05_fasttree_acr_phylogeny.sh` | `20gaoyuan/submission/code_availability/2fasttree.sh` | `2fasttree.sh` | muscle alignment, trimAl trimming (-gt 0.5) and FastTree maximum-likelihood tree of anti-CRISPR (Acr) proteins. |
| `pipelines/05_adaptive_traits/01_rgi_card.sh` | `20gaoyuan/submission/code_availability/1RGI_CARD.sh` | `1RGI_CARD.sh` | RGI/CARD screening of the MAG proteins for antibiotic resistance genes (ARG). |
| `pipelines/05_adaptive_traits/02_diamond_functional_screening.sh` | `20gaoyuan/submission/code_availability/2mag_orf_diamond_li.sh` | `2mag_orf_diamond_li.sh` | DIAMOND blastp against the CARD (ARG), MRG, mibPOP and VFG databases (>=60% identity, >=90% subject coverage). |
| `pipelines/05_adaptive_traits/03_diamond_result_pretreatment.sh` | `20gaoyuan/submission/code_availability/3diamond_res_pretreatment.sh` | `3diamond_res_pretreatment.sh` | Annotate each DIAMOND hit with its MAG of origin and collapse hits into per-MAG gene-count tables. |
| `pipelines/05_adaptive_traits/04_macrel.sh` | `20gaoyuan/submission/code_availability/02macrel.sh` | `02macrel.sh` | Macrel prediction of antimicrobial peptides from the MAG protein set. |
| `pipelines/05_adaptive_traits/05_antismash_batch_1.sh` | `20gaoyuan/submission/code_availability/02antiSMASH_1-1000.sh` | `02antiSMASH_1-1000.sh` | antiSMASH detection of biosynthetic gene clusters, MAG batch 1. |
| `pipelines/05_adaptive_traits/06_ogt_prediction.sh` | `20gaoyuan/submission/NC/code/1OTG_prediction2.sh` | `1OTG_prediction2.sh` | OGT_prediction pipeline for optimal growth temperature of the MAGs. |
| `pipelines/06_table_postprocessing/01_collect_DI_hits.py` | `20gaoyuan/submission/NC/code/01DIs_results_3_useful.py` | `01DIs_results_3_useful.py` | Collect the DefenseFinder per-island TSV files, keep `activity == Defense` rows and standardise the contig IDs. |
| `pipelines/06_table_postprocessing/02_DI_gene_ids.py` | `20gaoyuan/submission/NC/code/01all_DI_orf.py` | `01all_DI_orf.py` | Expand each defense-island record into the individual ORF/gene IDs it contains. |
| `pipelines/06_table_postprocessing/03_DI_summary.py` | `20gaoyuan/submission/NC/code/02DIs_summary4_useful.py` | `02DIs_summary4_useful.py` | Parse `5455MAG_defense_islands.txt` into a tidy per-island summary (gene counts, defense types, boundaries). |
| `pipelines/06_table_postprocessing/04_DI_nondG_orf_ids.py` | `20gaoyuan/submission/NC/code/03DI_nonDG.sh (lines 1-33)` | `03DI_nonDG.sh` | Expand every defense-island record into the ORF IDs of its contig (Python block extracted from the original mixed scratch file). |
| `pipelines/06_table_postprocessing/05_DI_nondG_filter.sh` | `20gaoyuan/submission/NC/code/03DI_nonDG.sh (lines 35-117)` | `-` | Subtract the defense-gene IDs from the DI ORF list, extract the remaining ORFs with seqkit and filter the matching GFF records (shell block extracted from the same scratch file). |

## R analyses

| New path | Original working-copy path | Previous repo path | Manuscript |
| --- | --- | --- | --- |
| `analyses/01_MAG_catalog/01_sample_and_MAG_metadata.R` | `20gaoyuan/01MAG_info/1MAG_info.R` | `R_code/1MAG_info.R` | Fig. S1 |
| `analyses/01_MAG_catalog/02_sampling_map.R` | `20gaoyuan/01MAG_info/2samping_map.R` | `R_code/2samping_map.R` | Fig. S1a-c |
| `analyses/01_MAG_catalog/03_MAG_taxonomy.R` | `20gaoyuan/01MAG_info/3_MAG_tax.R` | `R_code/3_MAG_tax.R` | Fig. S1e, Fig. S7 |
| `analyses/01_MAG_catalog/04_MAG_tree_annotation_bacteria.R` | `20gaoyuan/11_5455mag_tree/1_MAG_tree_annotation.R` | `R_code/1_MAG_tree_annotation.R` | Fig. 1a |
| `analyses/01_MAG_catalog/05_MAG_tree_annotation_archaea.R` | `20gaoyuan/11_42mag_ar_tree/1_arMAG_tree_annotation.R` | `R_code/1_arMAG_tree_annotation.R` | Fig. 1b |
| `analyses/01_MAG_catalog/06_supplementary_MAG_table.R` | `20gaoyuan/Supplementary_tables/1MAG_info.R` | `-` | Supplementary Data |
| `analyses/02_defensome_landscape/01_defense_system_overview.R` | `20gaoyuan/10Defensefinder/1defensefinder_res.R` | `R_code/1defensefinder_res.R` | Fig. 2, Fig. S2, Fig. S3 |
| `analyses/02_defensome_landscape/02_region_comparison.R` | `20gaoyuan/10Defensefinder/2defensefinder_res_Glacier_Plateau_Plain.R` | `R_code/2defensefinder_res_Glacier_Plateau_Plain.R` | Fig. 3a, Fig. 3e, Fig. S2 |
| `analyses/02_defensome_landscape/03_medium_comparison.R` | `20gaoyuan/10Defensefinder/3defensefinder_res_plankton_benthon.R` | `R_code/3defensefinder_res_plankton_benthon.R` | Fig. 3b, Fig. 3d, Fig. S4 |
| `analyses/02_defensome_landscape/04_density_and_genome_size.R` | `20gaoyuan/10Defensefinder/4defensefinder_res_density_compare.R` | `R_code/4defensefinder_res_density_compare.R` | Fig. 2h-i, Fig. S2e, Fig. S4h-i, Fig. S5c-d |
| `analyses/02_defensome_landscape/05_defense_system_composition_pcoa.R` | `20gaoyuan/10Defensefinder/5_MAG_DS_pcoa.R` | `R_code/5_MAG_DS_pcoa.R` | Fig. 2e-f |
| `analyses/02_defensome_landscape/06_taxonomic_composition.R` | `20gaoyuan/10Defensefinder/7_tax_composition.R` | `R_code/7_tax_composition.R` | Fig. S7 |
| `analyses/02_defensome_landscape/07_phylum_specialization.R` | `20gaoyuan/10Defensefinder/6_specialtax_DS.R` | `R_code/6_specialtax_DS.R` | Fig. S8, Fig. S9 |
| `analyses/02_defensome_landscape/08_season_comparison.R` | `20gaoyuan/10Defensefinder/8defensefinder_res_season_compare.R` | `-` | - |
| `analyses/03_mobile_genetic_elements/01_MGE_cargo_and_defense_density.R` | `20gaoyuan/12_5455MAG_MGEs/1_MGE_combination.R` | `R_code/1_MGE_combination.R` | Fig. 4a-c, Fig. S5, Fig. S6 |
| `analyses/03_mobile_genetic_elements/02_MGE_cargo_antidefense.R` | `20gaoyuan/12_5455MAG_MGEs/2_MGE_antidefense.R` | `R_code/2_MGE_antidefense.R` | Fig. S15 |
| `analyses/03_mobile_genetic_elements/03_representative_MGEs.R` | `20gaoyuan/13_special_MGE/1_special_MGE.R` | `R_code/1_special_MGE.R` | Fig. 4d, Fig. S10 |
| `analyses/04_defense_islands/01_defense_island_characterization.R` | `20gaoyuan/14_DIs/1_DIs.R` | `R_code/1_DIs.R` | Fig. 5a-d, Fig. S11 |
| `analyses/04_defense_islands/02_antidefense_in_defense_islands.R` | `20gaoyuan/15_antidefense/4antidefense_colocalization/1_anti_colocalization.R` | `-` | Fig. 5e |
| `analyses/04_defense_islands/03_nondG_functional_annotation.R` | `20gaoyuan/17DI_nonDG/1DF_DP.R` | `-` | Fig. S12 |
| `analyses/04_defense_islands/04_defensefinder_vs_defensepredictor.R` | `20gaoyuan/16DefensePredictor/1DF_DP.R` | `R_code/1DF_DP.R` | Fig. S13 |
| `analyses/05_antidefensome/01_phage_antidefensome.R` | `20gaoyuan/30allphage_antidefensome/1_allphage_antidefense.R` | `R_code/1_allphage_antidefense.R` | Fig. S14, Fig. S15 |
| `analyses/05_antidefensome/02_special_plasmid.R` | `20gaoyuan/15_antidefense/1special_plasmid/1_special_plasmid.R` | `R_code/1_special_plasmid.R` | Fig. 6a |
| `analyses/05_antidefensome/03_special_phage.R` | `20gaoyuan/15_antidefense/2special_phage/1_special_phage.R` | `R_code/1_special_phage.R` | Fig. 6b |
| `analyses/05_antidefensome/04_Acr_phylogeny.R` | `20gaoyuan/15_antidefense/3Acr_tree/Acr_tree_annotation.R` | `R_code/Acr_tree_annotation.R` | Fig. S16 |
| `analyses/06_adaptive_traits/01_ARG_MRG_VFG_traits.R` | `20gaoyuan/20ARGMRGVFG/1_Adaptive_traits.R` | `R_code/1_Adaptive_traits.R` | - |
| `analyses/07_cross_study_comparison/01_defense_system_comparison.R` | `20gaoyuan/40plateau_groundwater_comparision/1_defense_comparision.R` | `R_code/1_defense_comparision.R` | - |
| `analyses/07_cross_study_comparison/02_antidefense_comparison.R` | `20gaoyuan/40plateau_groundwater_comparision/2_antidefense_comparision.R` | `R_code/2_antidefense_comparision.R` | - |
| `analyses/08_revision_nc/01_MAG_and_phage_macrodiversity.R` | `93_NC_comments/01_diversity_comparison/1_MAG_Phage_macrodiversity.R` | `-` | Revision |
| `analyses/08_revision_nc/02_VHR_and_prophage_abundance.R` | `93_NC_comments/02_VHR_prophage_variation/1_VHR_prophage_abundance.R` | `-` | Revision / Fig. 2 |
| `analyses/08_revision_nc/03_defense_islands_in_NCBI_RefSeq.R` | `93_NC_comments/03_DIs_NCBI_RefSeq/1_DIs_NCBI_refseq_genome.R` | `-` | Revision |
| `analyses/08_revision_nc/04_environmental_factor_effects.R` | `93_NC_comments/04_environmental_factors/1_environmental_factor_effects.R` | `-` | Revision |
| `analyses/08_revision_nc/05_defense_system_relative_abundance.R` | `93_NC_comments/05_DS_relative_abundance/1_DS_relative_abundance.R` | `-` | Revision |

## Inputs that were deliberately left out

| File | Why |
| --- | --- |
| `20gaoyuan/12_5455MAG_MGEs/2_MGE_DS.R` | empty stub (4 lines, `read.delim2("")`), never run |
| `20gaoyuan/16DefensePredictor/1DF_DP(1).R` | older copy of `1DF_DP.R`, superseded by the refactored version |
| `20gaoyuan/16DefensePredictor/03DI_nonDG(1).sh` | byte-identical duplicate of `03DI_nonDG.sh` |
| `92_codes/3_MAG_tax.R` | byte-identical duplicate of `20gaoyuan/01MAG_info/3_MAG_tax.R` |
| `20gaoyuan/submission/code_availability/*` vs `20gaoyuan/submission/NC/code/*` | six scripts existed twice (identical apart from line endings); only one copy is kept |

## Paths that are still site-specific

The cluster job scripts under `pipelines/` are reproduced verbatim: their
`/lustre/home/lipengwei2024phd/...` and `/home/lipengwei2024phd/...` paths are the
directories used for the production runs and have to be adapted to your own
cluster before submission. See [`05_software_and_environment.md`](05_software_and_environment.md).

