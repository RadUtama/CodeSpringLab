#!/usr/bin/env bash


qsub  -l mem_free=50G -pe threads 4 -cwd -o ../../csl_results/${7}/log/output_homer_diffpeak.txt -e ../../csl_results/${7}/log/error_homer_diffpeak.txt homer_diffpeak.sh $1 $2 $3 $4 $5 $6


