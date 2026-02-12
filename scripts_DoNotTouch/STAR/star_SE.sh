
module load EBModules
module load STAR/2.7.10a-GCC-10.3.0

ulimit -n 10000

STAR --runThreadN 4 \
--quantMode TranscriptomeSAM \
       --outFileNamePrefix $1 \
       --genomeLoad NoSharedMemory \
       --genomeDir $2 \
       --outSAMtype BAM SortedByCoordinate \
       --outFilterMismatchNmax "2" \
       --outFilterMultimapNmax "2" \
       --outSAMunmapped None \
       --outSAMstrandField None \
       --readFilesCommand zcat \
       --readFilesIn $3

module load EBModules
module load SAMtools/1.14-GCC-10.3.0

samtools index -b ${1}Aligned.sortedByCoord.out.bam ${1}Aligned.sortedByCoord.out.bam.bai
