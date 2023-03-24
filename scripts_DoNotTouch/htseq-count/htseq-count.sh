
htseq-count --mode=union \
	--stranded=yes \
	--minaqual=10 \
	--type='exon' \
	--idattr='gene_id' \
	--order=pos \
	--format=bam \
	'/galaxy-central/database/files/000/120/dataset_120764.dat' "/galaxy-central/database/files/000/120/dataset_120775.dat" | \
awk '{if ($1 ~ "no_feature|ambiguous|too_low_aQual|not_aligned|alignment_not_unique") print $0 | "cat 1>&2"; else print $0}'