
module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0
module load SAMtools/1.14-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0
module load picard/2.21.6-Java-11

bowtie2  --very-sensitive -k 10 -X 1000 --dovetail --threads 8 \
	-x $2 \
	-1 $3 -2 $4 > ${1}_temp1.sam

samtools view -h -q 30 -f 0x2 -o ${1}_temp2.sam ${1}_temp1.sam
rm ${1}_temp1.sam

samtools sort -n -o ${1}_temp3.sam ${1}_temp2.sam
rm ${1}_temp2.sam

awk '( ($3 ~ /^chr/ && $3 != "chrM" && $3 != "chrUn") || (/^@/) )' ${1}_temp3.sam > ${1}Aligned.sortedByName.out.sam
rm ${1}_temp3.sam

samtools view -h -bS -o ${1}Aligned.sortedByName.out.bam ${1}Aligned.sortedByName.out.sam
rm ${1}Aligned.sortedByName.out.sam

#java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#I=${1}Aligned.sortedByName.out.bam \
#O=${1}Aligned.sortedByName_removeDup.out.bam \
#M=${1}_markedDup_metrics.txt
# Too many duplicates removed
#java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
#REMOVE_DUPLICATES=true READ_NAME_REGEX=null \
#I=${1}Aligned.sortedByName.out.bam \
#O=${1}Aligned.sortedByName_removeDup.out.bam \
#M=${1}_markedDup_metrics.txt

java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
I=${1}Aligned.sortedByName.out.bam \
O=${1}_insert_size_metrics.txt \
H=${1}_insert_size_histogram.pdf \
M=0.5

#java -jar $EBROOTPICARD/picard.jar CollectInsertSizeMetrics \
#I=${1}Aligned.sortedByName_removeDup.out.bam \
#O=${1}_insert_size_metrics_removeDup.txt \
#H=${1}_insert_size_histogram_removeDup.pdf \
#M=0.5

pdftoppm -jpeg ${1}_insert_size_histogram.pdf ${1}_insert_size_histogram
mv ${1}_insert_size_histogram-1.jpg ${1}_insert_size_histogram.jpg
rm ${1}_insert_size_histogram.pdf

bedtools bamtobed -i ${1}Aligned.sortedByName.out.bam > ${1}Aligned.sortedByName.out.bed
#bedtools bamtobed -i ${1}Aligned.sortedByName_removeDup.out.bam > #${1}Aligned.sortedByName_removeDup.out.bed

cat ../../csl_results/${5}/log/error_bowtie2.txt ${1}Log.final.out > ../../csl_results/${5}/log/error_bowtie2.txt
