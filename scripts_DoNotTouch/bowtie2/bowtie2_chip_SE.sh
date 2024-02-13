
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


java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
REMOVE_DUPLICATES=true \
I=${1}Aligned.sortedByName.out.bam \
O=${1}Aligned.sortedByName_removeDup.out.bam \
M=${1}_markedDup_metrics.txt

bedtools bamtobed -i ${1}Aligned.sortedByName.out.bam > ${1}Aligned.sortedByName.out.bed
bedtools bamtobed -i ${1}Aligned.sortedByName_removeDup.out.bam > ${1}Aligned.sortedByName_removeDup.out.bed

cat ../../csl_results/${5}/log/error_bowtie2.txt ${1}Log.final.out > ../../csl_results/${5}/log/error_bowtie2.txt
