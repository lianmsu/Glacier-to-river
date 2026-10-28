#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J is5455
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB
#SBATCH --output=job.lpw-isescan5455.%j.%N.out
# #SBATCH --reservation=<your-reservation>   # uncomment if your cluster requires one
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

source activate isescan
cd ${HOME_ROOT}/05plateau/23ISEscan_5455MAGs
id_file=${HOME_ROOT}/05plateau/23ISEscan_5455MAGs/5455MAG_path.txt
mkdir -p ./results
#start=15001
#end=20000
# 切片
#ids=$(awk "NR >= $start && NR <= $end" ${id_file})

#for i in ${ids}
#do

#isescan.py --seqfile ${HOME_ROOT}/03GW/01ANI99MAGs/${i} --output ./results/${i} --nthread 24

#done
cd ${HOME_ROOT}/05plateau/23ISEscan_5455MAGs
cat 5455MAG_path.txt | xargs -n 1 -P 8 -I{} isescan.py --seqfile {} --output ./results --nthread 6

