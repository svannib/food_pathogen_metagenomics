#!/usr/bin/env bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=30:00
#SBATCH --output=slurm_nanopore_concatenate
#SBATCH --job-name=nanopore_concatenate_reads

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

#load your required modules below
module load singularityce/3.10.2
singularity_path=/home/vbenve/share/bioinfo/singularity

# Declare working directory and sample batch as variables
folder=/net/cephfs/shares/amr.imm.uzh/data/projects/projects_VB/food_pathogen_metagenomics
batch=01_dna_extraction/batch_7
run=FC2

# Declare sample name and corresponding multiplexing barcode
sample_id=(7REFPE 7REFEC 7REFKM 7REFSE 7REFSL 7REFLM 7EMPTY)
barcodes=(barcode01 barcode02 barcode03 barcode04 barcode05 barcode06 barcode07)

#generate a new folder where the concatenated reads will be stored
mkdir -p $folder/01_data/"$batch"/"$run"/concatenated_reads/

# Loop trough the samples
for (( i=0; i<${#sample_id[@]}; i++ ))
do echo "Reads with "${barcodes[$i]}" concatenated in sample: "${sample_id[$i]}""

# Generate a folder for each sample and concatenate fastq reads there
mkdir -p $folder/01_data/"$batch"/"$run"/concatenated_reads/"${sample_id[$i]}"
zcat $folder/01_data/"$batch"/"$run"/fastq_pass/${barcodes[$i]}/*.gz \
>> $folder/01_data/"$batch"/"$run"/concatenated_reads/"${sample_id[$i]}"/"${sample_id[$i]}".fastq

# Convert the .fastq file to a .fasta
sed -n '1~4s/^@/>/p;2~4p' $folder/01_data/"$batch"/"$run"/concatenated_reads/"${sample_id[$i]}"/"${sample_id[$i]}".fastq \
> $folder/01_data/"$batch"/"$run"/concatenated_reads/"${sample_id[$i]}"/"${sample_id[$i]}".fasta

done


