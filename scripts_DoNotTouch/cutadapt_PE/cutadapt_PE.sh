#!/bin/sh

module load EBModules
module load Anaconda3/2022.05
conda init bash
conda activate cutadapt

cutadapt -j 4 -m $1 -a $2 -A $2 -o $3 -p $4 $5 $6

conda deactivate
