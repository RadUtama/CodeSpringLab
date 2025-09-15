#!/usr/bin/env bash

qsub -l mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${9}/log/output_diffbind.txt -e ../../csl_results/${8}/log/error_diffbind.txt -V ../scripts_DoNotTouch/DiffBind/diffbind.sh $1 $2 $3 $4 $5 $6 $7 $8

