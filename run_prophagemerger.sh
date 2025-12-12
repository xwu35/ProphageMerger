#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH -J run_ProphageMerger
#SBATCH --output=slurm-%x.%j.out
#SBATCH --error=slurm-%x.%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<your_email_address>

eval "$(conda shell.bash hook)"
conda activate snakemake

ProphageMerger.py --test_run
