#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J gwvc2.0.0
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=500GB

# 开始计时
start_time=$(date +%s)


source activate defensefinder2.0.0
#########################################################################
# File Name: GO_DefenseFinder.sh
# Date: 17/4/2024
#########################################################################
cd /home/lipengwei2024phd/05plateau/40GWVC_DF2

#proteins=/home/lipengwei2024phd/03GW/03mag_orfs/all_27578_orf/all_protein_seq.fa
genome=/home/lipengwei2024phd/03GW_eu/60GWVC_vOTUs/GWVC.vOTUs.fna
magdefense=/home/lipengwei2024phd/05plateau/40GWVC_DF2/res2
echo "当前目录为：$(pwd)"
echo "Running DefenseFinder"

defense-finder run "$genome" -a --preserve-raw -o "$magdefense"

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
