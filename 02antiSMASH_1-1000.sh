#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J bgc1-1000
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --output=job.01antiSMASH_lpw.%j.%N.out

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/43antiSMASH_2
mkdir -p ./res
source activate antiSMASH

id_file=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/13DefensePredictor/mags_list.txt
start=1
end=1000
# 切片
ids=$(awk "NR >= $start && NR <= $end" ${id_file})

for i in ${ids}
do

mkdir -p ./res/${i}

antismash /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/02metagenome_data/10_5467MAGs_90_5_10kb/${i}.fa \
	--taxon bacteria \
	--cpus 40 \
	--databases /lustre/home/lipengwei2024phd/software/miniforge3/envs/antiSMASH/lib/python3.11/site-packages/antismash/databases \
	--output-dir ./res/${i} \
	--allow-long-headers \
	--genefinding-tool prodigal-m

done


