#!/bin/sh
#SBATCH --partition=cpu
#SBATCH -J genomadMAG
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH -o job.genomad_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
DB_ROOT="${DB_ROOT:-${LUSTRE_HOME}/database}"       # reference databases
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

# 开始计时
start_time=$(date +%s)


source activate genomad
genomad_db=${DB_ROOT}/genomad_db
#########################################################################
cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/20genomad_5467MAGs
MAG_contig=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data/merged_5467MAG.fa
workdir=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/20genomad_5467MAGs

genomad end-to-end --cleanup ${MAG_contig} ${workdir}/res ${genomad_db}

# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "脚本运行时间：$runtime 秒"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0



