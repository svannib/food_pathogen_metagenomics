#!/usr/bin/env bash
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_nanopore_trim_align
#SBATCH --job-name=nanopore_trim_align
#SBATCH --array=1-7
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#load your required modules below
module load singularityce/3.10.2
singularity_path=/home/vbenve/share/bioinfo/singularity

# Indicate path to working directory, sample batch and raw reads
folder=/shares/amr.imm.uzh/data/projects/projects_VB/food_pathogen_metagenomics
batch=01_dna_extraction/batch_7
run=FC2

# Define array of samples based on the folders with the concatenated reads
cd "$folder"/01_data/"$batch"/"$run"/concatenated_reads
folders=$(find . -maxdepth 1 -type d)
array=(mock *)
export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}

# Define more useful variables
rawdata="$folder"/01_data/"$batch"/"$run"/concatenated_reads/"$sample_id"/"$sample_id".fastq #raw data fastq format
fasta="$folder"/01_data/"$batch"/"$run"/concatenated_reads/"$sample_id"/"$sample_id".fasta #raw data fasta format

#generate and access a new folder where the reads will be processed and the intermediate files will be stored
mkdir -p "$folder"/01_data/"$batch"/"$run"/trimmed_reads/$sample_id
cd "$folder"/01_data/"$batch"/"$run"/trimmed_reads/$sample_id

echo "Pipeline started for $sample_id" `date` > $sample_id.log

#Adapter trimming----------------------------------------------------------------------------------------------------
echo "Adapter trimming started" `date` >> $sample_id.log

# Adapter trimming, porechop takes 10'000 as samples to check adapter sets. Internal adapters will also be removed with --discard_middle.
srun $singularity_path/porechop_0.2.4--py310h30d9df9_3.sif porechop -t "$SLURM_CPUS_PER_TASK" -i $rawdata -v 1 -o "$sample_id"_trim.fastq --discard_middle

# Create a new directory and store results from nanoplot read QC there
mkdir -p nanoplot_trimmed_reads
srun $singularity_path/nanoplot_latest.sif NanoPlot --fastq "$sample_id"_trim.fastq -o nanoplot_trimmed_reads/

echo "Adapter trimming end" `date` >> $sample_id.log

# Host Removal----------------------------------------------------------------------------------------------------
echo "Host Removal started" `date` >> "$sample_id".log

# Map reads to reference host genome with minimap2
srun $singularity_path/minimap2_2.24--h7132678_1.sif minimap2 -ax map-ont \
"$folder"/01_data/00_reference_genomes/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna \
"$sample_id"_trim.fastq > $sample_id.sam

# After conversion to a binary alignment file, reads are sorted and tagged
srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools view -bS "$sample_id".sam > "$sample_id".bam
srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools sort "$sample_id".bam -T "$sample_id" -o "$sample_id".sorted.bam

# Reads are then indexed and statistics regarding alignment are generated
srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools index "$sample_id".sorted.bam
srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools flagstat "$sample_id".sorted.bam > "$sample_id"_hostaln.txt

#Retrieve unaligned sequences
srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools view -h -f 4 "$sample_id".sorted.bam -o "$sample_id"_unmapped.bam

# samtools is used to convert .BAM unaligned reads back to a fastq
srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools fastq -n "$sample_id"_unmapped.bam > "$sample_id"_prokaryotic.fastq

# Run nanoplot again on the reads not aligned to the reference genome
mkdir -p nanoplot_trimmed_prokaryotic_reads
srun $singularity_path/nanoplot_latest.sif NanoPlot --fastq "$sample_id"_prokaryotic.fastq -o nanoplot_trimmed_prokaryotic_reads

echo "Host Removal completed" `date`>> "$sample_id".log
echo "Pipeline finished for $sample_id" `date` > $sample_id.log
