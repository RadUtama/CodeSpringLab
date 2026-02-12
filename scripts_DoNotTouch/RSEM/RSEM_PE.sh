module load EBModules
module load RSeQC

#/grid/bsr/home/utama/.local/bin/infer_experiment.py -r $5 -i $1 > $4_strand.txt
infer_experiment.py -r $5 -i $1 > $4_strand.txt

fw=$(awk 'NR==5 {print $NF}' $4_strand.txt)
rv=$(awk 'NR==6 {print $NF}' $4_strand.txt)

diff=$(echo "scale=3 ;($fw - $rv)*100/($fw + $rv)" | bc)

if (( $(echo "${diff#-} < 50" | bc -l) ));then
        strand=none
        strand_idx=0
elif (( $(echo "$diff > 0" | bc -l) ));then
        strand=forward
        strand_idx=1
else
        strand=reverse
        strand_idx=2
fi

############ RSEM ################

module load EBModules
module load STAR/2.7.10a-GCC-10.3.0
module load RSEM/1.3.3-foss-2019b
module load SAMtools/1.16.1-GCC-11.3.0

ulimit -n 10000

rsem-calculate-expression --paired-end -p 8 --strandedness ${strand} \
	--bam $6 \
	$2 \
	$4

rm -rf ${4}.transcript.bam
rm -rf ${4}.stat

#################################

