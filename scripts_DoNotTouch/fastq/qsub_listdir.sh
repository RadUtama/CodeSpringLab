#!/usr/bin/env bash

qsub -l m_mem_free=1G -pe threads 1 -cwd -o ../../csl_results/${2}/log/output_listFastq.txt -e ../../csl_results/${2}log/error_listFastq.txt -V ../scripts_DoNotTouch/fastq/listdir.sh $1

