
module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

bowtie2  --very-sensitive --dovetail --threads 8 \
	--no-unal --no-mixed --phred33 \
	-x $2 \
	-U $3 > ${1}_temp1.sam

samtools view -h -q 30 -o ${1}_temp2.sam ${1}_temp1.sam
rm ${1}_temp1.sam

samtools sort -n -o ${1}_temp3.sam ${1}_temp2.sam
rm ${1}_temp2.sam

awk '( ($3 ~ /^chr/ && $3 != "chrM" && $3 != "chrUn") || (/^@/) )' ${1}_temp3.sam > ${1}Aligned.sortedByName.out.sam
rm ${1}_temp3.sam

samtools view -h -bS -o ${1}Aligned.sortedByName.out.bam ${1}Aligned.sortedByName.out.sam
rm ${1}Aligned.sortedByName.out.sam

samtools sort -o ${1}Aligned.sortedByCoord.out.bam ${1}Aligned.sortedByName.out.bam
rm ${1}Aligned.sortedByName.out.bam

samtools index -b ${1}Aligned.sortedByCoord.out.bam ${1}Aligned.sortedByCoord.out.bam.bai

#java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#REMOVE_DUPLICATES=true \
#I=${1}Aligned.sortedByName.out.bam \
#O=${1}Aligned.sortedByName_removeDup.out.bam \
#M=${1}_markedDup_metrics.txt

java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true \
I=${1}Aligned.sortedByCoord.out.bam \
O=${1}Aligned.sortedByCoord_removeDup.out.bam \
M=${1}_markedDup_metrics.txt

samtools index -b ${1}Aligned.sortedByCoord_removeDup.out.bam ${1}Aligned.sortedByCoord_removeDup.out.bam.bai

#bedtools bamtobed -i ${1}Aligned.sortedByName.out.bam > ${1}Aligned.sortedByName.out.bed
#bedtools bamtobed -i ${1}Aligned.sortedByName_removeDup.out.bam > ${1}Aligned.sortedByName_removeDup.out.bed
bedtools bamtobed -i ${1}Aligned.sortedByCoord.out.bam > ${1}Aligned.sortedByCoord.out.bed
bedtools bamtobed -i ${1}Aligned.sortedByCoord_removeDup.out.bam > ${1}Aligned.sortedByCoord_removeDup.out.bed

cat ../../csl_results/${7}/log/error_bowtie2.txt ${1}Log.final.out > ../../csl_results/${7}/log/error_bowtie2.txt

#####################################

module load EBModules
module load deepTools/3.5.2-foss-2022a

#--effectiveGenomeSize 2913022398 \ # Mice:2150570000; GRCh38:2913022398
### For male mice chrX should be ignored

#bamCoverage -b ${1}Aligned.sortedByCoord.out.bam \
#--normalizeUsing RPGC \
#--effectiveGenomeSize ${5} \
#--binSize 10 \
#--extendReads 200 \
#--ignoreForNormalization chrX chrM \
#--outFileFormat bedgraph \
#--outFileName ${1}Aligned.sortedByCoord.out.bdg

bamCoverage -b ${1}Aligned.sortedByCoord_removeDup.out.bam \
--normalizeUsing RPGC \
--effectiveGenomeSize ${5} \
--binSize 10 \
--extendReads 200 \
--ignoreForNormalization chrX chrM \
--outFileFormat bedgraph \
--outFileName ${1}Aligned.sortedByCoord_removeDup.out.bdg

awk '( $1 ~ /^chr/ && $1 != "chrM" && $1 != "chrUn" )' ${1}Aligned.sortedByCoord_removeDup.out.bdg > ${1}Aligned.sortedByCoord_removeDup_filtered.out.bdg

rm ${1}Aligned.sortedByCoord_removeDup.out.bdg

# The third line is chromsize in macs2 step
/grid/bsr/home/utama/bin/x86_64/bedGraphToBigWig \
${1}Aligned.sortedByCoord_removeDup_filtered.out.bdg \
${6} \
${1}Aligned.sortedByCoord_removeDup.out.bw

rm ${1}Aligned.sortedByCoord_removeDup_filtered.out.bdg
