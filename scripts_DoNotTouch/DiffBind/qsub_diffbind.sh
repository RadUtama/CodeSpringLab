#!/usr/bin/env bash
#SBATCH --job-name=diffbind
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE

source ../scripts_DoNotTouch/DiffBind/diffbind.sh $1 $2 $3 $4 $5 $6 $7 $8

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${9}/log/output_diffbind.txt -e ../../csl_results/${9}/log/error_diffbind.txt -V ../scripts_DoNotTouch/DiffBind/diffbind.sh $1 $2 $3 $4 $5 $6 $7 $8
