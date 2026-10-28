#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J prodigal
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01prodigal_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

source activate prokka
# 开始计时
start_time=$(date +%s)

cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/12prodigal_5467MAGs
mags=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data/merged_5467MAG.fa
res_path=${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/12prodigal_5467MAGs

echo "Running Prodigal"
prodigal -p meta -i "${mags}" \
		-f gff -o "${res_path}/merged_5467MAG_ORF.gff" \
		-m -d "${res_path}/merged_5467MAG_ORF.fna" \
		-a "${res_path}/merged_5467MAG_ORF.faa"




# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0