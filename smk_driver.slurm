#!/bin/bash
#SBATCH --job-name=snakemake_driver
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --output=logs/snakemake_driver_%j.out
#SBATCH --error=logs/snakemake_driver_%j.err

module load GCC/12.3.0 OpenMPI/4.1.5 snakemake/8.4.2 Anaconda3/2024.02-1

snakemake -s workflow/Snakefile --profile profiles/slurm
