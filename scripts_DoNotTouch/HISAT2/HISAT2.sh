
hisat2 -p ${GALAXY_SLOTS:-1} \
	-x '/galaxy-central/tool-data/genomes/HISAT2/HISAT2_indexes/UCSC/hg38/hg38/genome' \
	-1 'input_f.fastq' \
	-2 'input_r.fastq' | \
samtools sort - -@ ${GALAXY_SLOTS:-1} \
	-l 6 \
	-o '/galaxy-central/database/files/000/139/dataset_139605.dat'