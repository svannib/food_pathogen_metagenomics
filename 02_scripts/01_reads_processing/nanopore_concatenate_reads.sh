#!/bin/bash
#SBATCH --job-name=nanopore_concatenate_reads       #This is the name of your job
#SBATCH --cpus-per-task=4                 #This is the number of cores reserved
#SBATCH --mem-per-cpu=8G              #This is the memory reserved per core.
#Total memory reserved: 32GB

#SBATCH --time=00:30:00        #This is the time that your task will run
#SBATCH --qos=30min           #You will run in this queue
#SBATCH --array=1-7

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=slurm_nanopore_concatenate     #This is the joined STDOUT and STDERR file

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

#load your required modules below
ml purge
ml BBMap/38.91-foss-2018b

folder=/scicore/home/egliadr/GROUP/projects/food_pathogen_metagenomics_VB
batch=library_prep_test_VB

sample_id=(BS_Z1 BS_Z2 BS_Z3 CC_Z1 CC_Z2 CC_Z3 HMW)
barcodes=(barcode02 barcode03 barcode04 barcode05 barcode06 barcode07 barcode08)

mkdir $folder/01_data/"$batch"/concatenated_reads/

for (( i=0; i<${#sample_id[@]}; i++ ))
do echo "${sample_id[$i]}"
echo "${barcodes[$i]}"

mkdir $folder/01_data/"$batch"/concatenated_reads/"${sample_id[$i]}"
zcat $folder/01_data/"$batch"/fastq_pass/${barcodes[$i]}/*.gz \
>> $folder/01_data/"$batch"/concatenated_reads/"${sample_id[$i]}"/"${sample_id[$i]}".fastq

reformat.sh -Xmx4g in=$folder/01_data/"$batch"/concatenated_reads/"${sample_id[$i]}"/"${sample_id[$i]}".fastq \
out=$folder/01_data/"$batch"/concatenated_reads/"${sample_id[$i]}"/"${sample_id[$i]}".fasta qin=33

done
