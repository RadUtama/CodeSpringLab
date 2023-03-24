#!/usr/bin/env bash

qsub -l m_mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${7}/log/output_deseq2.txt -e ../../csl_results/${7}/log/error_deseq2.txt -V ../scripts_DoNotTouch/DESeq2/deseq2.sh $1 $2 $3 $4 $5 $6

