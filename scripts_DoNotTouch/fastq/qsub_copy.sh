#!/usr/bin/env bash
#SBATCH --job-name=copyFastq
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE

source ../scripts_DoNotTouch/fastq/copy.sh $1 $2 

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 1 -cwd -o ../../csl_results/${3}/log/output_copyFastq.txt -e ../../csl_results/${3}/log/error_copyFastq.txt -V ../scripts_DoNotTouch/fastq/copy.sh $1 $2
