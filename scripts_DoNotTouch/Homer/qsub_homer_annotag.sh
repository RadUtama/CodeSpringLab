#!/usr/bin/env bash


qsub  -l m_mem_free=50G -pe threads 4 -cwd -o ../../csl_results/${5}/log/output_homer_annotag.txt -e ../../csl_results/${5}/log/error_homer_annotag.txt homer_annotag.sh $1 $2 $3 $4


