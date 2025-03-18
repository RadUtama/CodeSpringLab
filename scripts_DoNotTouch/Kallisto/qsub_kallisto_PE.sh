#!/usr/bin/env bash

qsub -l mem_free=50G -pe threads 4 -cwd -o ../../csl_results/${5}/log/output_kallisto.txt -e ../../csl_results/${5}/log/error_kallisto.txt -V ../scripts_DoNotTouch/Kallisto/kallisto_PE.sh $1 $2 $3 $4

