#!/bin/sh

# This is a wrapper script for cutadapt; Use only from Galaxy

# $1 = min length
# $2 = adapter
# $3 = input file
# $4 = output file
# $5 = file type
# $6 = log file
# $7 = discard options

if [ -z "$6" ]; then
    echo "This script should be run from Galaxy" >&2
    exit 1
fi

MINLEN="$1"
ADAPT1="$2"
INFILE1="$3"
OUTFILE1="$4"
FILEEXT="$5"
LOGFILE="$6"
DISCARD="$7"

F1=`basename ${INFILE1}`
if [ "${FILEEXT}" == "fasta" ]; then
    F1="${F1}.fa"
    ln -s ${INFILE1} ${F1}
else
    F1="${F1}.fastq"
    ln -s ${INFILE1} ${F1}
fi    

cutadapt -m ${MINLEN} -a ${ADAPT1} $DISCARD -o ${OUTFILE1} ${F1} > ${LOGFILE} 2>&1

if [ $? -ne 0 ]; then
    echo "Error in cutadapt on read 1" >&2
    cat ${LOGFILE} >&2
    exit 1
fi

echo "Clipping completed without error"



bash /galaxy-central/tool-data/inhouse_tools/cshl_QC/cutadapt_se.sh '15' 'AGATCGGAAGAGC' '/galaxy-central/database/files/000/120/dataset_120744.dat' '/galaxy-central/database/files/000/142/dataset_142209.dat' 'fastqsanger' '/galaxy-central/database/files/000/142/dataset_142210.dat'