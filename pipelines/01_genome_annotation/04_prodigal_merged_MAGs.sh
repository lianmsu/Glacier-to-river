#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J prodigal
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01prodigal_lpw.%j.%N.out

source activate prokka
# 开始计时
start_time=$(date +%s)

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/12prodigal_5467MAGs
mags=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/02metagenome_data/merged_5467MAG.fa
res_path=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/12prodigal_5467MAGs

echo "Running Prodigal"
prodigal -p meta -i "${mags}" \
		-f gff -o "${res_path}/merged_5467MAG_ORF.gff" \
		-m -d "${res_path}/merged_5467MAG_ORF.fna" \
		-a "${res_path}/merged_5467MAG_ORF.faa"




# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0