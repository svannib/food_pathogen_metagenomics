#!/bin/bash
#SBATCH --job-name=nanopore_read_processing        #This is the name of your job
#SBATCH --cpus-per-task=10                  #This is the number of cores reserved
#SBATCH --mem-per-cpu=8G              #This is the memory reserved per core.
#Total memory reserved: 80GB

#SBATCH --time=24:00:00        #This is the time that your task will run
#SBATCH --qos=1day           #You will run in this queue
#SBATCH --array=1-6

# Paths to STDOUT or STDERR files should be absolute or relative to current working directory
#SBATCH --output=slurm_nanopore_trim_align_assemble_polish     #This is the joined STDOUT and STDERR file

#This job runs from the current working directory

#Remember:
#The variable $TMPDIR points to the local hard disks in the computing nodes.
#The variable $HOME points to your home directory.
#The variable $SLURM_JOBID stores the ID number of your job.

#load your required modules below
ml purge

# The source commands allows the cluster to access and activate personal conda environments
source /scicore/home/egliadr/benven0001/miniconda3/etc/profile.d/conda.sh

# the reads_processing environment contains porechop, nanoplot and minimap2
conda activate reads_processing

# Define array of samples
array=(mock A0CC A0F A1CC A1F REF0CC REF0F)
export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}

# Indicate path to working directory, sample batch and raw reads
folder=/scicore/home/egliadr/GROUP/projects/food_pathogen_metagenomics
batch=dna_extraction_test
run=20220519_run1
rawdata="$folder"/01_data/"$batch"/"$run"/concatenated_reads/"$sample_id"/"$sample_id".fastq #raw data fastq format
fasta="$folder"/01_data/"$batch"/"$run"/concatenated_reads/"$sample_id"/"$sample_id".fasta

#generate and access a new folder where the reads will be processed and the intermediate files will be stored
mkdir -p "$folder"/01_data/"$batch"/"$run"/trimmed_reads/$sample_id
cd "$folder"/01_data/"$batch"/"$run"/trimmed_reads/$sample_id

echo "Pipeline started for $sample_id" `date` > $sample_id.log

#Adapter trimming----------------------------------------------------------------------------------------------------
echo "Adapter trimming started" `date` >> $sample_id.log

# Adapter trimming, porechop takes 10'000 as samples to check adapter sets. Internal adapters will also be removed with --discard_middle.
#porechop -t "$SLURM_CPUS_PER_TASK" -i $rawdata -o "$sample_id"_trim.fastq --discard_middle

# Create a new directory and store results from nanoplot read QC there
mkdir -p nanoplot_trimmed_reads
NanoPlot --fastq "$sample_id"_trim.fastq -o nanoplot_trimmed_reads/

echo "Adapter trimming end" `date` >> $sample_id.log

# Host Removal----------------------------------------------------------------------------------------------------
echo "Host Removal started" `date` >> "$sample_id".log

# Map reads to reference host genome with minimap2
minimap2 -ax map-ont \
"$folder"/01_data/reference_genomes/GCF_016699485.2_bGalGal1.mat.broiler.GRCg7b_genomic.fna.gz \
"$sample_id"_trim.fastq > "$sample_id".sam

conda activate samtools
# After conversion to a binary alignment file, reads are sorted and tagged
samtools view -bS "$sample_id".sam > "$sample_id".bam
samtools sort "$sample_id".bam -T "$sample_id" -o "$sample_id".sorted.bam
# Reads are then indexed and statistics regarding alignment are generated
samtools index "$sample_id".sorted.bam
samtools flagstat "$sample_id".sorted.bam > "$sample_id"_hostaln.txt

#Retrieve unaligned sequences
samtools view -h -f 4 "$sample_id".sorted.bam -o "$sample_id"_unmapped.bam
samtools sort -n "$sample_id"_unmapped.bam -T "$sample_id"_unmapped -o "$sample_id"_sorted_unmapped.bam
rm "$sample_id"_unmapped.bam "$sample_id".sam "$sample_id".bam

# picard is used to convert .BAM alignment files back to a fastq
picard SamToFastq -I "$sample_id"_sorted_unmapped.bam -F "$sample_id"_UC_picard.fastq

# Run nanoplot again on the reads not aligned to the reference genome
conda activate reads_processing
mkdir -p nanoplot_trimmed_prokaryotic_reads
NanoPlot --fastq "$sample_id"_UC_picard.fastq -o nanoplot_trimmed_prokaryotic_reads

echo "Host Removal completed" `date`>> "$sample_id".log

#De novo assembly using flye----------------------------------------------------------------------------------------------------
echo "Denovo assembly started" `date` >> "$sample_id".log

# Assemble genomes from metagenomic reads
mkdir -p "$folder"/03_results/"$batch"/"$run"/assemblies/"$sample_id"/
cd "$folder"/03_results/"$batch"/"$run"/assemblies/"$sample_id"/

conda activate flye
flye --nano-raw "$folder"/01_data/"$batch"/"$run"/trimmed_reads/"$sample_id"/"$sample_id"_UC_picard.fastq  \
--out-dir "$sample_id".flye.assembly --threads "$SLURM_CPUS_PER_TASK" --meta

echo "#Assembly completed" `date` >> "$sample_id".log

#Assembly polishing with racon and medaka----------------------------------------------------------------------------------------------------
echo "Assembly polishing started" `date` >> "$sample_id".log

# Convert assemblies to .sam file
conda activate reads_processing
# map reads to assembly
minimap2 -ax map-ont ./"$sample_id".flye.assembly/assembly.fasta $rawdata > "$sample_id"_minimap.sam

#polish incorrect raw contigs
conda activate racon
racon -t "$SLURM_CPUS_PER_TASK" $rawdata "$sample_id"_minimap.sam ./"$sample_id".flye.assembly/assembly.fasta > "$sample_id"_racon_conse.fa

# Use medaka to create a consensus sequence through pileup of individual sequencing reads against a draft assembly
conda activate medaka
medaka_consensus -i $fasta -d "$sample_id"_racon_conse.fa -o "$sample_id"_medaka_consen  -t "$SLURM_CPUS_PER_TASK" -m r941_min_hac_g507

echo "Assembly polishing completed" `date` >> "$sample_id".log

#QUAST assembly quality control----------------------------------------------------------------------------------------------------
echo "QUAST QC started" `date` >> "$sample_id".log
conda activate quast

cd "$folder"/03_results/"$batch"/"$run"/assemblies/"$sample_id"
mkdir -p quast_reports

python /scicore/home/egliadr/benven0001/miniconda3/envs/quast/bin/metaquast.py -o ./quast_reports \
-r "$folder"/01_data/reference_genomes/ZymoBIOMICS.STD.refseq.v2/Genomes/ "$sample_id".flye.assembly/assembly.fasta

echo "QUAST QC finished" `date` >> "$sample_id".log

echo "Pipeline finished for $sample_id" `date` > $sample_id.log
