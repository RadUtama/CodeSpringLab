#!/usr/bin/env bash

qsub -l mem_free=5G -pe threads 8 -cwd -o ../../csl_results/${7}/log/output_macs2.txt -e ../../csl_results/${7}/log/error_macs2.txt -V ../scripts_DoNotTouch/MACS2/macs2_chip_SE.sh $1 $2 $3 $4 $5 $6 $8

