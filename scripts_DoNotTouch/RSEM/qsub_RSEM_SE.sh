#!/usr/bin/env bash

qsub -l m_mem_free=1G -pe threads 8 -cwd -o ../../csl_results/${7}/log/output_rsem.txt -e ../../csl_results/${7}/log/error_rsem.txt -V ../scripts_DoNotTouch/RSEM/RSEM_SE.sh $1 $2 $3 $4 $5 $6
