
/grid/bsr/home/utama/.local/bin/infer_experiment.py -r $5 -i $1 > $4_strand.txt

fw=$(awk 'NR==5 {print $NF}' $4_strand.txt)
rv=$(awk 'NR==6 {print $NF}' $4_strand.txt)

diff=$(echo "scale=3 ;($fw - $rv)*100/($fw + $rv)" | bc)

if (( $(echo "${diff#-} < 50" | bc -l) ));then
        strand=unstranded
        strand_idx=0
elif (( $(echo "$diff > 0" | bc -l) ));then
        strand=forward
        strand_idx=1
else
        strand=reverse
        strand_idx=2
fi

module load EBModules
module load Subread

suffix=counts.txt

featureCounts -a $2 \
        -T 2 \
        -t "exon" \
        -g ${strand_idx} \
        -s 0 \
        -Q 12 \
        --minOverlap 1 \
        -o $4_${suffix} \
        $1

