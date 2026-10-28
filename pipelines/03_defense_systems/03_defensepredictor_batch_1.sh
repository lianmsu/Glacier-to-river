#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J DP1-1000
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.DefensePredictor_lpw.%j.%N.out

# 开始计时
start_time=$(date +%s)

module load gcc/7.2.0
source activate defense_prefictor
#########################################################################

#########################################################################
id_file=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/13DefensePredictor/mags_list.txt
out_path=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/13DefensePredictor/out

# 
start=1
end=1000

# 切片
ids=$(awk "NR >= $start && NR <= $end" ${id_file})

for i in ${ids}
do
gff=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res/${i}/${i}.gff
ffn=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res/${i}/${i}.ffn
faa=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/11prokka_5467MAGs/res/${i}/${i}.faa

defense_predictor --prokka_gff ${gff} --prokka_ffn ${ffn} --prokka_faa ${faa} --device cpu --output ${out_path}/${i}.csv

done

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
