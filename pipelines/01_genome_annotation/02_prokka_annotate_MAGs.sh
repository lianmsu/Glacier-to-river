#!/bin/bash
#SBATCH --job-name=p5467
#SBATCH --partition=cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01prokka_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

cp ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data/mags_list.txt ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/11prokka_5467MAGs

source activate prokka
cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/11prokka_5467MAGs
file=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/mags_list.txt
MAG_path=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data/10_5467MAGs_90_5_10kb

res_path=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res
mkdir -p ${res_path}

for i in $(cat $file)
do
prokka ${MAG_path}/"${i}.fa" --outdir ${res_path}/${i} --prefix ${i} --locustag ${i} --metagenome --kingdom Bacteria --cpus 28
done

