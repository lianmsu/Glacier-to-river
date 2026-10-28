#!/bin/bash
#SBATCH -o job.gtdbtk.%j.%N.out 
#SBATCH --partition=cpu
#SBATCH -J gtdb6
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB

source activate gtdbtk-2.4.0-test

#GMOC
binning_path=/home/lipengwei2024phd/05plateau/015455MAG_gtdb/5455MAG
out_path=/home/lipengwei2024phd/05plateau/015455MAG_gtdb/res6

start_time=$(date +%s)
mkdir -p ${out_path}
cd ${out_path}
#mkdir -p res
gtdbtk classify_wf --cpus $SLURM_NTASKS --pplacer_cpus $SLURM_NTASKS -x fa --genome_dir ${binning_path}/ --out_dir ${out_path} --skip_ani_screen
#--mash_db /home/lipengwei2024phd/05plateau/015455MAG_gtdb/mash

end_time=$(date +%s)
cost_time=$[ $end_time-$start_time ]
echo "running time is $(($cost_time/60)) min" >> run_time_gtdb_${SLURM_JOBID}.log



