#!/usr/bin/env bash

qsub -l mem_free=2G -pe threads 4 -cwd -o ../../csl_results/${8}/log/output_cutadapt.txt -e ../../csl_results/${8}/log/error_cutadapt.txt -V ../scripts_DoNotTouch/cutadapt_PE/cutadapt_PE.sh $1 $2 $3 $4 $5 $6 $7

