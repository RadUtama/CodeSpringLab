

module load EBModules
module load MACS2/2.2.9.1-foss-2022b
## Alternative local installation:
## /grid/bsr/home/utama/.local/bin/macs2

macs2 callpeak --nomodel \
	-t ${2} -q ${7} \
	-f BED -g ${3} \
	--shift -100 --extsize 200 --bdg --call-summits \
	-n ${1} --keep-dup  all \
	--outdir ${5}

#module load Anaconda3/2023.03-1
#conda activate deeptools

module load EBModules
module load deepTools/3.5.2-foss-2022a

### TSS was used for workshop instead of center
    
computeMatrix reference-point -p 4 \
    --referencePoint TSS \
    -b 1000 -a 1000 \
    -R ${6} \
    -S ${9}/${1}/${1}Aligned.sortedByCoord_removeDup.out.bw \
    --skipZeros \
    -o ${5}/${1}_TSS.gz \
    --outFileSortedRegions ${5}/${1}_genes_TSS.bed

#computeMatrix reference-point -p 4 \
#    --referencePoint center \
#    -b 1000 -a 1000 \
#    -R ${6} \
#    -S ${5}/${1}Aligned.sortedByCoord_removeDup.out.bw \
#    --skipZeros \
#    -o ${5}/${1}_peakCenter.gz \
#    --outFileSortedRegions ${5}/${1}_genes_peakCenter.bed
     
     
plotHeatmap -m ${5}/${1}_TSS.gz \
    -out ${5}/${1}_heatmap_TSS.png \

#plotHeatmap -m ${5}/${1}_peakCenter.gz \
#    -out ${5}/${1}_heatmap_peakCenter.png \

#conda deactivate

###### Homer Annotate #########

module load EBModules
module load Anaconda3/2021.05
module load R/4.1.2-foss-2021a

export PATH=$PATH:/grid/bsr/data/data/utama/tools/homer/bin/

annotatePeaks.pl ${5}/${1}_peaks.xls ${8} > ${5}/${1}_peaks_annotated.txt
