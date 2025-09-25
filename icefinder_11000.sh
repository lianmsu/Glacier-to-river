#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J ICE_11000
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=250GB
#SBATCH --reservation=lipengwei2024phd_res2

# 开始计时
start_time=$(date +%s)
source activate icefinder
#########################################################################
# File Name: GO_DefenseFinder.sh
# Date: 17/4/2024
#########################################################################

cp -r /home/lipengwei2024phd/03GW_prokaryotes/try/ICEfinder2_linux /home/lipengwei2024phd/03GW_prokaryotes/04ICEfinder/ICEfinder2_11000
ICE_workdir=/home/lipengwei2024phd/03GW_prokaryotes/04ICEfinder/ICEfinder2_11000
cd ${ICE_workdir}
id_file=/home/lipengwei2024phd/03GW_prokaryotes/04ICEfinder/all_19847MAGsid.txt
MAG_dir=/home/lipengwei2024phd/03GW/01ANI99MAGs

start=10001
end=11000
# 切片
ids=$(awk "NR >= $start && NR <= $end" ${id_file})

for i in ${ids}
do

python ICEfinder2.py -i ${MAG_dir}/${i} -t Metagenome

done

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
