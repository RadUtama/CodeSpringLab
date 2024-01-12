

# GFFUtils

#gtf2bed FILE.gtf > FILE.bed

# pyGenomeTracks

# Good for tracks:
anno=/grid/bsr/data/data/utama/genome/hg38_p13_gencode/gencode.v42.chr_patch_hapl_scaff.annotation.gtf

label=treated2

treated2=~/csl_results/atac/data/macs2/treated_2/treated_2_treat_pileup.bw

treated2_bed=~/csl_results/atac/data/macs2/treated_2/treated_2_summits.bed

init=genomeTracks.ini
trackplot=genomeTracks.png

## On Conda

module load EBModules
module load Anaconda3/2023.03-1

exec $SHELL

conda activate gffutils
conda activate pygenometracks

make_tracks_file \
	--trackFiles ${control1} ${control2} ${control3} ${treated1} ${treated2} ${treated3} ${control1_bed} ${control2_bed} ${control3_bed} ${treated1_bed} ${treated2_bed} ${treated3_bed} ${anno} \
	-o ${init}
#make_tracks_file \
#	--trackFiles ${control1} ${control2} ${control3} ${treated1} ${treated2} ${treated3} ${control1_bed} ${control2_bed} ${control3_bed} ${treated1_bed} ${treated2_bed} ${treated3_bed} ${anno} \
#	-o ${init}

sed -i 's/labels\s=\sfalse/labels=true/g' ${init}
#sed -i 's/labels\s=\sfalse/labels=true/g' ${init}

pyGenomeTracks --tracks ${init} --region ${3} -o ${trackplot}
#pyGenomeTracks --tracks ${init} --region chr1:2700000-2800000 -o ${trackplot}

######## local installed but not piped well somehow

#/grid/bsr/home/utama/.local/bin/make_tracks_file \
#	--trackFiles ${control1} ${control2} ${control3} ${treated1} ${treated2} ${treated3} ${anno} \
#	-o ${init}

#/grid/bsr/home/utama/.local/bin/pyGenomeTracks --tracks ${init} --region chr1:1000000-4000000 -o ${trackplot}

conda deactivate

