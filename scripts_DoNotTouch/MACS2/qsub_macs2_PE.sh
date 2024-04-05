#!/usr/bin/env bash

qsub -l m_mem_free=5G -pe threads 8 -cwd -o ../../csl_results/${8}/log/output_macs2.txt -e ../../csl_results/${8}/log/error_macs2.txt -V ../scripts_DoNotTouch/MACS2/macs2_PE.sh $1 $2 $3 $4 $5 $6 $7

