
mkdir -p /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/30minced_5455MAG/blastdb_allphage

makeblastdb -in /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau_phage/00virus_contigs/cdhit_before/virus_select_all_medium-10kb_high_complete.fa -dbtype nucl -out /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/30minced_5455MAG/blastdb_allphage/allphage.db -hash_index -max_file_sz '4GB'

#sz服务器无法makeblastdb建库，在bb建好copy过来的。

mkdir -p /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/30minced_5455MAG/blastdb_MAGplasmid

makeblastdb -in /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/20genomad_5467MAGs/res/merged_5467MAG_summary/merged_5467MAG_plasmid.fna -dbtype nucl -out /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/30minced_5455MAG/blastdb_MAGplasmid/MAGplasmid.db -hash_index -max_file_sz '4GB'

