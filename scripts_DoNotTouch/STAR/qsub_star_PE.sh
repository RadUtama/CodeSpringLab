#!/usr/bin/env bash

qsub -l mem_free=50G -pe threads 4 -cwd -o ../../csl_results/${5}/log/output_star.txt -e ../../csl_results/${5}/log/error_star.txt -V ../scripts_DoNotTouch/STAR/star_PE.sh $1 $2 $3 $4

