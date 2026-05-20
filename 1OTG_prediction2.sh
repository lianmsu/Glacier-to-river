#!/usr/bin/bash
#SBATCH -o job-OGT.%j.%N.out
#SBATCH --partition=cpu
#SBATCH -J OGT
#SBATCH --get-user-env
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --mem=450GB

## this pipeline is used to predict CDS and annotate functions of them.

####################### Parameters ########################
# input, output and bin path
input=/home/lipengwei2024phd/05plateau/015455MAG_gtdb/5455MAG/* #input fasta file
out=/home/lipengwei2024phd/05plateau/30OGT_pre/res1 #output dir

mkdir -p ${out}
cd ${out}
#OGT_prediction
external_tools=/home/lipengwei2024phd/software/OGT_prediction/prediction/external_tools.txt #Create a file called external_tools.txt listing bedtools, tRNAscan-SE, barrnap, and prodigal tab separated from the absolute path for each executable.
ogt=/home/lipengwei2024phd/software/OGT_prediction/prediction/prediction_pipeline.py #prediction_pipeline.py
model=/home/lipengwei2024phd/software/OGT_prediction/data/prediction_demo/regression_models/ # OGT models path

source activate OGT_prediction

#
header="species\tsuperkingdom\tphylum\tclass\torder\tfamily"
echo -e "$header" >> tax_info.txt

for i in $input
do

#prefix="B6_bin.5" #prefix of genome
prefix=$(basename ${i} .fa)

taxonomy=$(awk -F"\t" -v prefix="$prefix" '{if($1==prefix) print $2}' /home/lipengwei2024phd/05plateau/30OGT_pre/gtdb_tax.tsv)

######################### OGT_Main ############################
real_out=$(realpath "$out")

#cd $real_out && 
if [ ! -f external_tools.txt ];then ln -s $external_tools;fi && \
if [ ! -d genomes/$prefix ];then mkdir -p genomes/$prefix;fi && \
gzip -c ${i} > genomes/$prefix/$prefix.fa.gz && \
echo -e "$prefix.fa.gz\t$prefix" >> genomes_info.txt && \
header="species\tsuperkingdom\tphylum\tclass\torder\tfamily" && \
tax=$(echo $taxonomy | sed 's/__;/__None;/g' | sed "s/d__\(.*\);p__\(.*\);c__\(.*\);o__\(.*\);f__\(.*\);g__.*/\1\t\2\t\3\t\4\t\5/") && \
echo -e "$prefix\t$tax" >> tax_info.txt 

done

python3 $ogt $model genomes_info.txt tax_info.txt && \
touch $real_out/OGT.all.done

#linux中，sed ""   与sed '' 效果不同。

