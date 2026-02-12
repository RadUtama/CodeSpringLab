#!/usr/bin/env bash
#SBATCH --job-name=fastQC
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE

source ../scripts_DoNotTouch/FastQC/fastqc.sh $1 $2

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${3}/log/output_fastQC.txt -e ../../csl_results/${3}/log/error_fastQC.txt -V ../scripts_DoNotTouch/FastQC/fastqc.sh $1 $2
