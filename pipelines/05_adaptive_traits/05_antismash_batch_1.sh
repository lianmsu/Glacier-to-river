#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J bgc1-1000
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=job.01antiSMASH_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/43antiSMASH_2
mkdir -p ./res
source activate antiSMASH

id_file=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/13DefensePredictor/mags_list.txt
start=1
end=1000
# 切片
ids=$(awk "NR >= $start && NR <= $end" ${id_file})

for i in ${ids}
do

mkdir -p ./res/${i}

antismash ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data/10_5467MAGs_90_5_10kb/${i}.fa \
	--taxon bacteria \
	--cpus 40 \
	--databases ${LUSTRE_HOME}/software/miniforge3/envs/antiSMASH/lib/python3.11/site-packages/antismash/databases \
	--output-dir ./res/${i} \
	--allow-long-headers \
	--genefinding-tool prodigal-m

done


