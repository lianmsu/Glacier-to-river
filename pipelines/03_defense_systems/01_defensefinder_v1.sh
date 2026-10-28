#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J df2
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20


# 开始计时
start_time=$(date +%s)

source activate defensefinder
#########################################################################

#########################################################################

MAGs_path=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/02metagenome_data/merged_5467MAG.fa
defense_workdir=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/10DefenseFinder2.0.0_5467MAGs
output_dir="${defense_workdir}/res1"

mkdir -p ${output_dir}
cd ${defense_workdir}

defense-finder run -a ${MAGs_path} --preserve-raw -o "$output_dir"


# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
