#!/usr/bin/env bash
#SBATCH --job-name=deseq2
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=4
#SBATCH --export=NONE

source ../scripts_DoNotTouch/DESeq2/deseq2.sh $1 $2 $3 $4 $5 $6 $7

# Ori script with grid qsub
#qsub -l mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${8}/log/output_deseq2.txt -e ../../csl_results/${8}/log/error_deseq2.txt -V ../scripts_DoNotTouch/DESeq2/deseq2.sh $1 $2 $3 $4 $5 $6 $7
