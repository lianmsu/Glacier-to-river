#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J rgicard
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.02RGI_CARD_lpw.%j.%N.out

cd /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG
source activate RGI_CARD

#rgi load --card_json /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/database/ARG/card-data/card.json --local

rgi main --input_sequence /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/5467MAG.prt_defensefinder.faa \
	--output_file /lustre/home/lipengwei2024phd/04Qinghai_Tibet_Plateau/42ARG/0ARG_RGI/rgiargs \
	--local \
	--clean \
	-t protein \
	--num_threads 28
