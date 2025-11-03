#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J run_ProphageMerger
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=xiaofen

eval "$(conda shell.bash hook)"
conda activate snakemake

python ProphageMerger.py -g test_data/genome_sequences.fa -o test_output
