#!/bin/bash
#SBATCH --job-name=GTDBTK          #This is the name of your job
#SBATCH --cpus-per-task=1        #This is the number of cores reserved
#SBATCH --mem-per-cpu=300G       #This is the memory reserved per core.
#Total memory reserved: 300GB
#SBATCH --time=24:00:00    #This is the time that your task will run
#SBATCH --qos=1day     #You will run in this queue

ml purge
source /scicore/home/egliadr/benven0001/miniconda3/etc/profile.d/conda.sh
conda activate gtdbtk-2.1.0
gtdbtk test --out_dir /scicore/home/egliadr/GROUP/projects/food_pathogen_metagenomics/03_results/gtdb_test/ \
--cpus 1 --debug
