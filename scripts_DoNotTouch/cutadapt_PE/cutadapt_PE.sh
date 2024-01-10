# !/bin/sh

module load EBModules
module load cutadapt/4.4-GCCcore-12.2.0

#module load Anaconda3/2022.05
#conda init bash
#conda activate cutadapt


cutadapt -j 4 -m $1 -a $2 -A $3 -o $4 -p $5 $6 $7

#conda deactivate
