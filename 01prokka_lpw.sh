#!/bin/bash
#SBATCH --job-name=p5467
#SBATCH --partition=cpu 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01prokka_lpw.%j.%N.out

cp /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/02metagenome_data/mags_list.txt /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs

source activate prokka
cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs
file=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/mags_list.txt
MAG_path=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/02metagenome_data/10_5467MAGs_90_5_10kb

res_path=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res
mkdir -p ${res_path}

for i in $(cat $file)
do
prokka ${MAG_path}/"${i}.fa" --outdir ${res_path}/${i} --prefix ${i} --locustag ${i} --metagenome --kingdom Bacteria --cpus 28
done

