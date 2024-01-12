# !/bin/bash

echo "Homer"

module load EBModules
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a

########### Homer differential peaks ###############

getDifferentialPeaksReplicates.pl -t ${3} -i ${2} -DESeq2 -genome ${4} -f 0.0001 -q 1 > ${1}/DiffPeak_${6}_vs_${5}(ref).txt

############ Homer annotate differential peaks ############

annotatePeaks.pl ${1}/DiffPeak_${6}_vs_${5}(ref).txt ${4} -raw > ${1}/DiffPeak_${6}_vs_${5}(ref)_annotated.txt


