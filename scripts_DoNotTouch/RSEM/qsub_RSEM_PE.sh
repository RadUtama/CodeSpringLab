#!/usr/bin/env bash
#SBATCH --job-name=rsem
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=8
#SBATCH --export=NONE

source ../scripts_DoNotTouch/RSEM/RSEM_PE.sh $1 $2 $3 $4 $5 $6

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 8 -cwd -o ../../csl_results/${7}/log/output_rsem.txt -e ../../csl_results/${7}/log/error_rsem.txt -V ../scripts_DoNotTouch/RSEM/RSEM_PE.sh $1 $2 $3 $4 $5 $6
