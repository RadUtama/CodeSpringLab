#!/usr/bin/env bash
#SBATCH --job-name=listDir
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

source ../scripts_DoNotTouch/fastq/listdir.sh $1

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 1 -cwd -o ../../csl_results/${2}/log/output_listFastq.txt -e ../../csl_results/${2}log/error_listFastq.txt -V ../scripts_DoNotTouch/fastq/listdir.sh $1
