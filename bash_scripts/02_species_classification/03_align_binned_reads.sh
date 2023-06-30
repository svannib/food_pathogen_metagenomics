#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --output=slurm_nanopore_align_binned_reads
#SBATCH --job-name=nanopore_align_binned_reads
#SBATCH --array=1-1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

#Import necessary modules
module load singularityce/3.10.2
singularity_path=/home/vbenve/share/bioinfo/singularity

#Define useful variables
wd=/shares/amr.imm.uzh/data/projects/projects_VB/food_pathogen_metagenomics
batch=01_dna_extraction/batch_7
run=FC2

# Define array of samples
cd "$wd"/03_results/"$batch"/"$run"/kraken/
folders=$(find . -maxdepth 1 -type d)
array=(mock *)
export sample_id=${array["$SLURM_ARRAY_TASK_ID"]}

#Define useful variables
raw="$wd"/01_data/"$batch"/"$run"/trimmed_reads/"$sample_id"/"$sample_id"_trim.fastq
reference_genomes=$wd/01_data/00_reference_genomes/zymo_bacterial_genomes

#create a dictionary of species and taxonomy IDs
declare -A pathogens=(["Pseudomonas_aeruginosa"]=287 ["Escherichia_coli"]=562 ["Salmonella_enterica"]=28901 \
["Lactobacillus_fermentum"]=1613 ["Enterococcus_faecalis"]=1351 ["Staphylococcus_aureus"]=1280 \
["Listeria_monocytogenes"]=1639 ["Bacillus_subtilis"]=1423 ["Acinetobacter_baumanni"]=470 ["Klebsiella_pneumoniae"]=573 \
["Staphylococcus_lugdunensis"]=28035 ["Klebsiella_michiganensis"]=1134687 ["Klebsiella_oxytoca"]=571)

# move to the directory with the kraken results
cd $sample_id/binned_reads
echo "aligning reads from "$sample_id" with samtools"

#run a loop over the pathogen dictionary
for pathogen in "${!pathogens[@]}";
do
  # Align the binned reads to the reference genome
  echo $"$pathogen - ${pathogens[$pathogen]}";
  srun $singularity_path/minimap2_2.24--h7132678_1.sif minimap2 -ax map-ont \
  "$reference_genomes"/"$pathogen"_zymo.f* \
  "$sample_id"_"$pathogen".fa > "$sample_id"_"$pathogen".sam

  # convert the .sam alignment to a .bam alignment
  srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img \
  samtools view -bS "$sample_id"_"$pathogen".sam > "$sample_id"_"$pathogen".bam

  # Sort the reads
  srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img \
  samtools sort -o "$sample_id"_"$pathogen"_sorted.bam "$sample_id"_"$pathogen".bam
done

#delete empty files
find . -size 0 -delete

#Write the list of .bam file to a text file
ls *_sorted.bam > alignment_list.txt

#loop over the sorted bam files and write their coverage to a tab separated file
for sorted_bam_file in *_sorted.bam; do
    echo -ne "BAM File: $sorted_bam_file\t" >> alignment_coverage_depth.tsv
    srun $singularity_path/quay.io-biocontainers-samtools-1.13--h8c37831_0.img samtools coverage "$sorted_bam_file" >> alignment_coverage_depth.tsv
done