
module load EBModules-LegacyBNB
module load FastQC/0.11.8-Java-1.8

fastqc -t 4 -o $2 $1

module load EBModules
module load Bowtie2/2.4.4-GCC-10.3.0

/grid/bsr/data/data/utama/tools/FastQ-Screen-0.15.2/fastq_screen -threads 4 -outdir $2 $1

