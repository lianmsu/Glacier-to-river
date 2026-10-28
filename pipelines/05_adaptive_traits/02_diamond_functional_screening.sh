#!/bin/sh
#SBATCH -o jobdia1060.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J diamond
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
#SBATCH --mail-type=end

source activate eggnog

orf="/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/5467MAG.prt_defensefinder.faa"
out=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG

#各种类型功能基因的数据库
ARGdb=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/database/ARG/CARD.dmnd
MRGdb=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/database/MRG/MRG.dmnd
POPDegdb=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/database/mibPOP/mibPOP.dmnd
VFGdb=/lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/database/VFG/VFDB_core.dmnd

query=${orf}

diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${ARGdb} -q ${query} -o ${out}/1ARG/diamond_ARG_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${MRGdb} -q ${query} -o ${out}/3MRG/diamond_MRG_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${POPDegdb} -q ${query} -o ${out}/4mibPOP/diamond_POPDeg_out --id 60 --subject-cover 90
diamond blastp -k 1 -e 0.00001 -p $SLURM_NTASKS -d ${VFGdb} -q ${query} -o ${out}/2VFG/diamond_VFG_out --id 60 --subject-cover 90
#--id 80 比较好
#The host system is detected to have 404 GB of RAM. It is recommended to increase the block size for better performance using these parameters : -b8 -c1
#--outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
#

#--id 80, 
echo "successfully"

