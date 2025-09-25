#!/bin/sh
#SBATCH -o job.%j.%N.out
#SBATCH --partition=fat
#SBATCH -J amp
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --output=job.01macrel_lpw.%j.%N.out

source activate macrel
macrel peptides \
    --fasta extracted_proteins.faa \
    --output ./out_peptides \
	--threads 28

