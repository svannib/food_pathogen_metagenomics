#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=slurm_nanopore_kraken2
#SBATCH --job-name=nanopore_kraken2
#SBATCH --array=1-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#Import necessary modules
module load singularityce/3.10.2
singularity_path=/home/vbenve/share/bioinfo/singularity
folder=/shares/amr.imm.uzh/data/projects/projects_VB/food_pathogen_metagenomics
batch=01_dna_extraction/batch_7
run=FC2

# Define array of samples based on the folders with the trimmed reads
cd "$folder"/01_data/"$batch"/"$run"/trimmed_reads
folders=$(find . -maxdepth 1 -type d)
array=(mock *)
export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}

#Define useful variables
Kraken2=/shares/amr.imm.uzh/bioinfo/databases/Kraken_standard_db_20230605
raw="$folder"/01_data/"$batch"/"$run"/trimmed_reads/"$sample_id"/"$sample_id"_prokaryotic.fastq

# create directories for storing results
mkdir -p "$folder"/03_results/"$batch"/"$run"/kraken/$sample_id
cd "$folder"/03_results/"$batch"/"$run"/kraken/$sample_id

echo "classifying reads from "$sample_id" with kraken2"

#run kraken2 with the kraken2 database
srun $singularity_path/kraken2_2.1.2--pl5321h9f5acd7_3.sif kraken2 \
--db $Kraken2 --threads $SLURM_CPUS_PER_TASK \
--report "$sample_id"_kraken.k2report \
--report-minimizer-data --minimum-hit-groups 3 \
$raw > "$sample_id"_kraken.kraken2
