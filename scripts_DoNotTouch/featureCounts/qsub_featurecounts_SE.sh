#!/usr/bin/env bash

qsub -l mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${6}/log/output_featurecounts.txt -e ../../csl_results/${6}/log/error_featurecounts.txt -V ../scripts_DoNotTouch/featureCounts/featurecounts_SE.sh $1 $2 $3 $4 $5
