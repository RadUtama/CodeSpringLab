# !/bin/bash

echo "Homer"

module load EBModules
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a

########### Homer create tag directory #######

makeTagDirectory ${2} ${1}Aligned.sortedByName.out.bed

############ Homer annotate peaks ############

annotatePeaks.pl ${3}_summits.bed ${4} -raw > ${3}_summits_annotated.txt


