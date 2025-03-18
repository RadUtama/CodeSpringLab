#!/usr/bin/env bash

qsub -l mem_free=5G -pe threads 4 -cwd -o ../../csl_results/${7}/log/output_genomeTracks.txt -e ../../csl_results/${7}/log/error_genomeTracks.txt -V ../scripts_DoNotTouch/genomeTracks/genomeTracks.sh $1 $2 $3 $4 $5 $6

