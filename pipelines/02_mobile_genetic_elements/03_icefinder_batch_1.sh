#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J ICE_1000
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=400GB
#SBATCH --reservation=lipengwei2024phd_res2

# 开始计时
start_time=$(date +%s)
source activate icefinder
#########################################################################
# File Name: GO_DefenseFinder.sh
# Date: 17/4/2024
#########################################################################

#cp -r /home/lipengwei2024phd/05plateau/21icefinder_5455MAGs/ICEfinder2_moban /home/lipengwei2024phd/05plateau/21icefinder_5455MAGs/ICEfinder2_1000
ICE_workdir=/home/lipengwei2024phd/05plateau/21icefinder_5455MAGs/ICEfinder2_1000
cd ${ICE_workdir}
id_file=/home/lipengwei2024phd/05plateau/21icefinder_5455MAGs/5455MAG.id
MAG_dir=/home/lipengwei2024phd/05plateau/015455MAG_gtdb/5455MAG

start=1
end=1000
# 切片
ids=$(awk "NR >= $start && NR <= $end" ${id_file})

for i in ${ids}
do

python ICEfinder2.py -i ${MAG_dir}/"${i}.fa" -t Metagenome

done

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
