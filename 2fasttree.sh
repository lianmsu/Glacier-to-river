#!/bin/sh
#SBATCH -o job.tree503.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J tree503
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=24
#SBATCH --mem=250GB

source activate fasttree
cd /home/lipengwei2024phd/03GW/30Acr_tree
#1、先用muscle把protein序列对齐形成多序列比对（MSA）
muscle -align /home/lipengwei2024phd/03GW/30Acr_tree/all_Acr.faa -output 503acr_aln.afa

#2、然后用trimal修剪
trimal -in 503acr_aln.afa -gt 0.5 > 503acr_aln.trim.afa
#-gt 0.5代表删去有50%以上gap的位点。-gt 0.5这个是对应这篇文章方法中的参数 https://www.nature.com/articles/s41467-023-38400-0

#3、再用fasttree对修建后的MSA进行建树
FastTree 503acr_aln.trim.afa > 503aln.trim.tree