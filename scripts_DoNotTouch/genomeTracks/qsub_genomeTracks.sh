#!/usr/bin/env bash
#SBATCH --job-name=genomeTracks
#SBATCH --mem-per-cpu=5G
#SBATCH --cpus-per-task=4

source ../scripts_DoNotTouch/genomeTracks/genomeTracks.sh $1 $2 $3 $4 $5 $6

# Ori script with grid qsub
#qsub -l mem_free=5G -pe threads 4 -cwd -o ../../csl_results/${7}/log/output_genomeTracks.txt -e ../../csl_results/${7}/log/error_genomeTracks.txt -V ../scripts_DoNotTouch/genomeTracks/genomeTracks.sh $1 $2 $3 $4 $5 $6
