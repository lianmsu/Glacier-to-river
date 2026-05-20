#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J pd
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB
#SBATCH --mail-type=end

module load prodigal/2.6.3-GCCcore-11.3.0
# 开始计时
start_time=$(date +%s)

cd /home/lipengwei2024phd/05plateau/02_5455MAGprodigal
mags=/home/lipengwei2024phd/05plateau/02_5455MAGprodigal/all5455MAG_contig.fa
magorf=/home/lipengwei2024phd/05plateau/02_5455MAGprodigal

echo "Running Prodigal"
prodigal -p meta -i "${mags}" \
		-f gff -o "${magorf}/genes.gff" \
		-m -d "${magorf}/nucleotide_seq.fa" \
		-a "${magorf}/protein_seq.fa"




# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0