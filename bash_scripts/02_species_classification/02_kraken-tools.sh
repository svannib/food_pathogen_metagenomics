#!/usr/bin/env bash
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=slurm_nanopore_krakentools
#SBATCH --job-name=nanopore_krakentools
#SBATCH --array=1-7
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

# Define path to raw reads
raw="$wd"/01_data/"$batch"/"$run"/trimmed_reads/"$sample_id"/"$sample_id"_trim.fastq

#create a dictionary of species and taxonomy IDs
declare -A pathogens=(["Pseudomonas_aeruginosa"]=287 ["Escherichia_coli"]=562 ["Salmonella_enterica"]=28901 \
["Lactobacillus_fermentum"]=1613 ["Enterococcus_faecalis"]=1351 ["Staphylococcus_aureus"]=1280 \
["Listeria_monocytogenes"]=1639 ["Bacillus_subtilis"]=1423 ["Acinetobacter_baumanni"]=470 ["Klebsiella_pneumoniae"]=573 \
["Staphylococcus_lugdunensis"]=28035 ["Klebsiella_michiganensis"]=1134687 ["Klebsiella_oxytoca"]=571)

# move to the directory with the kraken results
cd $sample_id
echo "extracting reads from "$sample_id" with kraken-tools"

#run a kraken-tools loop to extract reads matching the species in the dictionary
for pathogen in "${!pathogens[@]}";
do
  echo $"$pathogen - ${pathogens[$pathogen]}";
  srun $singularity_path/krakentools_1.2--pyh5e36f6f_0.sif extract_kraken_reads.py \
  -k $PWD/"$sample_id"_minikraken.kraken2 \
  -s $raw -t ${pathogens[$pathogen]} --include-children \
  -r $PWD/"$sample_id"_minikraken.k2report -o "$sample_id"_"$pathogen".fa
done

#delete empty files
find . -size 0 -delete

# create and fill a subdirectory
mkdir binned_reads
mv *.fa binned_reads
