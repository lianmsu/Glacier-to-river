#!/bin/bash
#SBATCH -o job.gtdbtk.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J allbinrename
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=14
# ---------------------------------------------------------------------------
# Site-specific settings -- set these (or export them) before submitting the job.
# The absolute paths of the production runs have been replaced by these
# variables, so the script contains no site-specific path any more.
# ---------------------------------------------------------------------------
LUSTRE_HOME="${LUSTRE_HOME:-/path/to/lustre/home}"   # large-storage home (project + databases)
HOME_ROOT="${HOME_ROOT:-/path/to/home}"             # login-node home
# ---------------------------------------------------------------------------

cd ${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data
#给fa文件中的contig重命名###############################################################
input_dir="${LUSTRE_HOME}/04Qinghai_Tibet_Plateau/02metagenome_data/00all_bins"

for input_file in "$input_dir"/*.fa; do
    basename=$(basename "$input_file" .fa)
    awk -v basename="$basename" '/^>/ {print ">" basename "-" ++i; next} {print}' "$input_file" > tmpfile && mv tmpfile "$input_file"
    echo "已处理: $input_file"
done
###########################################################################################


