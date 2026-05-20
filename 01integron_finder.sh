#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J integron_5455
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB
#SBATCH --reservation=lipengwei2024phd_res2

# 开始计时
start_time=$(date +%s)



#########################################################################
cd /home/lipengwei2024phd/05plateau/22integronfinder_5455MAGs
source activate integron_finder
#########################################################################


mags="/home/lipengwei2024phd/05plateau/22integronfinder_5455MAGs/all5455MAG_contig.fa"

integron_finder --local-max --func-annot ${mags} --cpu $SLURM_NTASKS --outdir /home/lipengwei2024phd/05plateau/22integronfinder_5455MAGs/res --gbk --pdf


# 结束计时
end_time=$(date +%s)
# 计算运行时间
runtime=$((end_time - start_time))
echo "running time：$runtime s"
#输出完成消息
tput setaf 2; echo "Done!"; tput sgr0
