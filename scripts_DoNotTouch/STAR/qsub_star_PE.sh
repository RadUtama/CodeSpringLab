#!/usr/bin/env bash
#SBATCH --job-name=star
#SBATCH --mem-per-cpu=50G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE

source ../scripts_DoNotTouch/STAR/star_PE.sh $1 $2 $3 $4

# Ori script with grid qsub
#qsub -l mem_free=50G -pe threads 4 -cwd -o ../../csl_results/${5}/log/output_star.txt -e ../../csl_results/${5}/log/error_star.txt -V ../scripts_DoNotTouch/STAR/star_PE.sh $1 $2 $3 $4
