# !/bin/bash
# Load modules: EBmodules, Anaconda3/2021.05,  R/4.1.2-foss-2021a 

echo "Homer"

module load EBModules
module load Anaconda3/2021.05

################## ALL PATHS #################

#SAMPLE_NAME=OE-462-PS
#BAM=./bam/${SAMPLE_NAME}.sortedByCoord.bam
#SAM=./bam/${SAMPLE_NAME}.sortedByCoord.sam
#OUT_TAG=./homer/tag/${SAMPLE_NAME}_tag

#TARGET1=./homer/tag/OE-394-PS_tag/
#TARGET2=./homer/tag/OE-399-PS_tag/
#TARGET3=./homer/tag/OE-462-PS_tag/

#CONTROL1=./homer/tag/OE-339-PS_tag/
#CONTROL2=./homer/tag/OE-349-PS_tag/
#CONTROL3=./homer/tag/OE-353-PS_tag/


#OUT=OlfacEpi_vs_SkinFibro_diffpeaks.txt
#OUT=N2A.miR_vs_N2A_diffpeaks.txt
#OUT=Omni_vs_OlfacEpi_diffpeaks.txt

########### Homer create tag directory #######

#exec $SHELL

#makeTagDirectory

makeTagDirectory ${2} ${1}

#makeTagDirectory ${OUT_TAG} ${SAM}

#############################################

#echo "success"

