#!/usr/bin/env bash

qsub -l m_mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${3}/log/output_fastQC.txt -e ../../csl_results/${3}/log/error_fastQC.txt -V ../scripts_DoNotTouch/FastQC/fastqc.sh $1 $2

