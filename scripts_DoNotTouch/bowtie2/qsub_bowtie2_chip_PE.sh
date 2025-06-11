#!/usr/bin/env bash

qsub -l mem_free=50G -pe threads 8 -cwd -o ../../csl_results/${7}/log/output_bowtie2.txt -e ${1}Log.final.out -V ../scripts_DoNotTouch/bowtie2/bowtie2_chip_PE.sh $1 $2 $3 $4 $5 $6 $7

