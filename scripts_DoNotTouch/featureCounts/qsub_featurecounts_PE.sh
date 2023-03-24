#!/usr/bin/env bash

qsub -l m_mem_free=1G -pe threads 4 -cwd -o ../../csl_results/${5}/log/output_featurecounts.txt -e ../../csl_results/${5}/log/error_featurecounts.txt -V ../scripts_DoNotTouch/featureCounts/featurecounts_PE.sh $1 $2 $3 $4

