#!/usr/bin/env bash
#
#SBATCH -J nextflow_mega # Job name
#SBATCH --ntasks-per-node=4 # one core
#SBATCH -N 1 # on one node
#SBATCH -t 1-00:00 # days
#SBATCH --mem 20G
#SBATCH -o /standard/vol185/cphg_Manichaikul/users/csm6hg/err/nextflow_dna.out # Standard output
#SBATCH -e /standard/vol185/cphg_Manichaikul/users/csm6hg/err/nextflow_dna.err # Standard error
#SBATCH -p standard
#SBATCH --account manichaikul

# Modules to load
module load nextflow

# Working directory
cd /standard/vol185/cphg_Manichaikul/users/csm6hg/nextflow_dna

# Run nextflow
nextflow run main_TensorQTL.nf -profile slurm -resume
