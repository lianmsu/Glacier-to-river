#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J rgicard
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.02RGI_CARD_lpw.%j.%N.out
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/42ARG
source activate RGI_CARD

#rgi load --card_json ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/42ARG/database/ARG/card-data/card.json --local

rgi main --input_sequence ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/42ARG/5467MAG.prt_defensefinder.faa \
	--output_file ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/42ARG/0ARG_RGI/rgiargs \
	--local \
	--clean \
	-t protein \
	--num_threads 28
