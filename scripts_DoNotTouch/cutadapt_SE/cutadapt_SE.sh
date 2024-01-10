#!/bin/sh

module load EBModules
module load cutadapt/4.4-GCCcore-12.2.0

#module load Anaconda3/2022.05
#conda init bash
#conda activate cutadapt

cutadapt -j 4 -m $1 -a $2 -o $3 $5

#conda deactivate
