#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J integron_5455
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB
# #SBATCH --reservation=<your-reservation>   # uncomment if your cluster requires one
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

# 开始计时
start_time=$(date +%s)



#########################################################################
cd ${HOME_ROOT}/05plateau/22integronfinder_5455MAGs
source activate integron_finder
#########################################################################


mags="${HOME_ROOT}/05plateau/22integronfinder_5455MAGs/all5455MAG_contig.fa"

integron_finder --local-max --func-annot ${mags} --cpu $SLURM_NTASKS --outdir ${HOME_ROOT}/05plateau/22integronfinder_5455MAGs/res --gbk --pdf


# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
