#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J cv
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=end


#进入checkv并导入数据库目录
source activate checkv
export CHECKVDB=/lustre/home/lipengwei2024phd/database/checkv_db/checkv-db-v1.5

virus_contig=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau_phage/00virus_contigs/cdhit_before/virus_select_all_medium-10kb_high_complete.fa
workdir=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau_phage/02CheckV

#end-to-end：污染度、完整度分析和去除原噬菌体污染
checkv end_to_end ${virus_contig} ${workdir} -t $SLURM_NTASKS

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"


