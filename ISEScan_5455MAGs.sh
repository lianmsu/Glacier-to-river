#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J is5455
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB
#SBATCH --output=job.lpw-isescan5455.%j.%N.out
#SBATCH --reservation=lipengwei2024phd_res2

source activate isescan
cd /home/lipengwei2024phd/05plateau/23ISEscan_5455MAGs
id_file=/home/lipengwei2024phd/05plateau/23ISEscan_5455MAGs/5455MAG_path.txt
mkdir -p ./results
#start=15001
#end=20000
# 切片
#ids=$(awk "NR >= $start && NR <= $end" ${id_file})

#for i in ${ids}
#do

#isescan.py --seqfile /home/lipengwei2024phd/03GW/01ANI99MAGs/${i} --output ./results/${i} --nthread 24

#done
cd /home/lipengwei2024phd/05plateau/23ISEscan_5455MAGs
cat 5455MAG_path.txt | xargs -n 1 -P 8 -I{} isescan.py --seqfile {} --output ./results --nthread 6

