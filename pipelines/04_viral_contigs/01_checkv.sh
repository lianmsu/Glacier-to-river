#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J cv
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=end
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
DB_ROOT="${DB_ROOT:-${LUSTRE_HOME}/database}"       # reference databases
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------


#进入checkv并导入数据库目录
source activate checkv
export CHECKVDB=${DB_ROOT}/checkv_db/checkv-db-v1.5

virus_contig=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau_phage/00virus_contigs/cdhit_before/virus_select_all_medium-10kb_high_complete.fa
workdir=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau_phage/02CheckV

#end-to-end：污染度、完整度分析和去除原噬菌体污染
checkv end_to_end ${virus_contig} ${workdir} -t $SLURM_NTASKS

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "execution time is $(($cost_time/60))min $(($cost_time%60))s"


