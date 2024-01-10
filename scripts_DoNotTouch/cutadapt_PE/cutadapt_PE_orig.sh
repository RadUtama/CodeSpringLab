#!/bin/sh

# This is a wrapper script for cutadapt; Use only from Galaxy

# $1 = min length
# $2 = first adapter
# $3 = second adapter
# $4 = input file 1
# $5 = input file 2
# $6 = output file 1
# $7 = output file 2
# $8 = file type
# $9 = log file
# $10 = discard options

if [ -z "$8" ]; then
    echo "This script should be run from Galaxy" >&2
    exit 1
fi

MINLEN="$1"
ADAPT1="$2"
ADAPT2="$3"
INFILE1="$4"
INFILE2="$5"
OUTFILE1="$6"
OUTFILE2="$7"
FILEEXT="$8"
LOGFILE="$9"
DISCARD="${10}"

F1=`basename ${INFILE1}`
F2=`basename ${INFILE2}`
if [ "${FILEEXT}" == "fasta" ]; then
    F1="${F1}.fa"
    ln -s ${INFILE1} ${F1}
    F2="${F2}.fa"
    ln -s ${INFILE2} ${F2}
    if [ $? -ne 0 ];then
        echo "Same file was used for paired-end clipping." >&2
        echo "Please ensure that each file corresponds to a single end of the paired-end read" >&2
        exit 1
    fi

else
    F1="${F1}.fastq"
    ln -s ${INFILE1} ${F1}
    F2="${F2}.fastq"
    ln -s ${INFILE2} ${F2}
    if [ $? -ne 0 ];then
        echo "Same file was used for paired-end clipping." >&2
        echo "Please ensure that each file corresponds to a single end of the paired-end read" >&2
        exit 1
    fi
fi    

if [ "${FILEEXT}" == "fasta" ]; then
    cutadapt -m ${MINLEN} -a ${ADAPT1} $DISCARD -o tmp.1.fa -p tmp.2.fa ${F1} ${F2} > ${LOGFILE} 2>&1
else
    cutadapt -m ${MINLEN} -a ${ADAPT1} $DISCARD -o tmp.1.fastq -p tmp.2.fastq ${F1} ${F2} > ${LOGFILE} 2>&1
fi

if [ $? -ne 0 ]; then
    echo "Error in cutadapt on read 1" >&2
    cat ${LOGFILE} >&2
    exit 1
fi

if [ "${FILEEXT}" == "fasta" ]; then
    cutadapt -m ${MINLEN} -a ${ADAPT2} $DISCARD -o ${OUTFILE2} -p ${OUTFILE1} tmp.2.fa tmp.1.fa >> ${LOGFILE} 2>&1
else
    cutadapt -m ${MINLEN} -a ${ADAPT2} $DISCARD -o ${OUTFILE2} -p ${OUTFILE1} tmp.2.fastq tmp.1.fastq >> ${LOGFILE} 2>&1
fi

if [ $? -ne 0 ]; then
    echo "Error in cutadapt on read 2" >&2
    cat ${LOGFILE} >&2
    exit 1
fi

echo "Clipping completed without error"




bash /galaxy-central/tool-data/inhouse_tools/cshl_QC/cutadapt_pe.sh '15' 'AGATCGGAAGAGC' 'AGATCGGAAGAGC' '/galaxy-central/database/files/000/142/dataset_142209.dat' '/galaxy-central/database/files/000/120/dataset_120745.dat' '/galaxy-central/database/files/000/142/dataset_142211.dat' '/galaxy-central/database/files/000/142/dataset_142212.dat' 'fastqsanger' '/galaxy-central/database/files/000/142/dataset_142213.dat '