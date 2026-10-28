#!/bin/bash
#SBATCH --job-name=blastn
#SBATCH --partition=cpu 
#SBATCH --nodes=1 
#SBATCH --get-user-env
#SBATCH --ntasks-per-node=24
#SBATCH --output=job-plasmid.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

source activate blast2.16

mkdir -p ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/2spacer2plasmid

blastn -db ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/blastdb_MAGplasmid/MAGplasmid.db -query ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/1minced_res/5455MAG_contig_withcas_minced_spacers.fa -out ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/2spacer2plasmid/spacers_plasmid.blastn.out -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' -num_threads 24 -evalue 1e-5

cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/2spacer2plasmid
cat spacers_plasmid.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > spacers_plasmid.blastn.1mis95cov.out
cat spacers_plasmid.blastn.1mis95cov.out | cut -f2 | sort | uniq > targeted.plasmid.id

mkdir -p ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/2spacer2allphage

blastn -db ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/blastdb_allphage/allphage.db -query ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/1minced_res/5455MAG_contig_withcas_minced_spacers.fa -out ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/2spacer2allphage/spacers_allphage.blastn.out -outfmt '6 qseqid sseqid pident evalue mismatch bitscore length qlen qcovs qcovhsp slen gapopen qstart qend sstart send' -num_threads 24 -evalue 1e-5

cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/30minced_5455MAG/2spacer2allphage
cat spacers_allphage.blastn.out | awk '{if($5<=1 && $9>=95) print $0}' > spacers_allphage.blastn.1mis95cov.out
cat spacers_allphage.blastn.1mis95cov.out | cut -f2 | sort | uniq > targeted.allphage.id


