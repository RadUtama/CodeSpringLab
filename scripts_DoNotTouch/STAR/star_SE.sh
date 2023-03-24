
module load EBModules
module load STAR/2.7.10a-GCC-10.3.0

STAR --runThreadN 4 \
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

