#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J bac
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mail-type=end

# 记录开始时间（秒级时间戳）
start_time=$(date +%s)

#进入checkv并导入数据库目录
source activate BACPHLIP

virus_contig=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau_phage/00virus_contigs/cdhit_before/virus_select_all_medium-10kb_high_complete.fa
workdir=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau_phage/03BACPHLIP

cd ${workdir}
bacphlip -i ${virus_contig} --multi_fasta


end_time=$(date +%s)
cost_time=$(( end_time - start_time ))
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"


