

module load EBModules
module load MACS2/2.2.9.1-foss-2022b
## Alternative local installation:
## /grid/bsr/home/utama/.local/bin/macs2

macs2 callpeak --nomodel \
	-t ${2} \
	-f BED -g ${3} \
	--shift -100 --extsize 200 --bdg --call-summits \
	-n ${1} --keep-dup  all \
	--outdir ${5}

/grid/bsr/home/utama/bin/x86_64/bedGraphToBigWig \
	${5}/${1}_treat_pileup.bdg \
	${4} \
	${5}/${1}_treat_pileup.bw

#module load EBModules
#module load Anaconda3/2023.03-1

#echo 'modules'

#conda init bash
##exec $SHELL

#echo 'shell'

#conda activate deeptools
    
#echo 'conda'

#computeMatrix reference-point -p 4 \
#    --referencePoint TSS \
#    -b 1000 -a 1000 \
#    -R ${6} \
#    -S ${5}/${1}_treat_pileup.bw \
#    --skipZeros \
#    -o ${5}/${1}_TSS.gz \
#    --outFileSortedRegions ${5}/${1}_genes_TSS.bed

#echo 'TSS'

#computeMatrix reference-point -p 4 \
#    --referencePoint center \
#    -b 1000 -a 1000 \
#    -R ${6} \
#    -S ${5}/${1}_treat_pileup.bw \
#    --skipZeros \
#    -o ${5}/${1}_peakCenter.gz \
#    --outFileSortedRegions ${5}/${1}_genes_peakCenter.bed
     
#echo 'center'     
     
#plotHeatmap -m ${5}/${1}_TSS.gz \
#    -out ${5}/${1}_heatmap_TSS.png \

#echo 'heatmap tss'

#plotHeatmap -m ${5}/${1}_peakCenter.gz \
#    -out ${5}/${1}_heatmap_peakCenter.png \

#echo 'heatmap center'

#conda deactivate
