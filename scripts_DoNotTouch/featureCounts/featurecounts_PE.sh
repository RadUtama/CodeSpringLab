
module load EBModules
module load Subread

suffix=counts.txt

featureCounts -a $2 \
        -T 2 \
        -p --countReadPairs \
        -t "exon" \
        -g $3 \
        -s 0 \
        -Q 12 \
        --minOverlap 1 \
        -C \
        -o $4_${suffix} \
        $1

