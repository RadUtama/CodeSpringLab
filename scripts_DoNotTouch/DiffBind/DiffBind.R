
######################################################

.libPaths("/grid/bsr/data/data/oeldemer/4.3")
.libPaths()
getwd()
library("DiffBind")
library(dplyr)

args = commandArgs(trailingOnly=TRUE)

outpath <- args[1]
designpath <- args[2]
peakpath <- args[3]
refcond <- args[4]
compared <- args[5]
genome <- args[6]
bampath <- args[7]

design <- read.delim(paste(designpath,"/design_matrix.txt",sep=""),header=T,sep="\t",row.names = 1,check.names=FALSE)
design <- design %>% filter_all(any_vars(. %in% c(refcond,compared)))
design <- design[1:(ncol(design)-1)]

SampleID <- rownames(design)

diffcol <- sapply(design, function(x) any(x==refcond))
design <- as.data.frame(design[ , diffcol])
colnames(design) <- "Condition"

design$SampleID <- SampleID
design <- design[, c("SampleID","Condition")]

design_refcond <- design[design$Condition==refcond,] 
design_compared <- design[design$Condition==compared,]

design_refcond$Replicate <- 1:nrow(design_refcond)
design_compared$Replicate <- 1:nrow(design_compared)

design <- rbind(design_refcond,design_compared)

### Add rep3 just for example dataset ###
if (grepl("example_dataset",peakpath,fixed=TRUE)){
	design <- rbind(design,list("CD4_day1_rep3","day1",3),list("CD4_day3_rep3","day3",3))
}
#########################################

bam <- rep("Aligned.sortedByCoord_removeDup.out.bam",nrow(design))
narrowpeak <- rep("_peaks.narrowPeak",nrow(design))

inbamlist <- rep(bampath,nrow(design))
inpeaklist <- rep(peakpath,nrow(design))
                  
design$bamReads <- paste(inbamlist,"/",design$SampleID,"/",design$SampleID,bam,sep="")
design$Peaks <- paste(inpeaklist,"/",design$SampleID,"/",design$SampleID,narrowpeak,sep="")
design$PeakCaller <- rep("narrowpeak",,nrow(design))

names(design)
design

dbobject<-dba(sampleSheet=design)
dbobject

png(file=paste(outpath,"/diffbind_CorrHeatmap_byPeaks.png",sep=""), width = 1000, height = 800, res = 150)
plot(dbobject)
dev.off()

peakdata<-dba.show(dbobject)$Intervals
if (genome=="human"){
	dbobject<-dba.blacklist(dbobject,blacklist=DBA_BLACKLIST_HG38,greylist=FALSE)
} else if (genome=="mouse"){
	#dbobject<-dba.blacklist(dbobject,blacklist=paste(designpath,"/mm39-blacklist.bed",sep=""),greylist=FALSE)
}
dbobject

dbobject<-dba.count(dbobject)
dbobject

info<-dba.show(dbobject)
info

libsizes<-cbind(LibReads=info$Reads,FRiP=info$FRiP,PeakReads=round(info$Reads*info$FRiP))
rownames(libsizes)<-info$ID
libsizes

png(paste(outpath,"/diffbind_CorrHeatmap_byCounts.png",sep=""), width = 1000, height = 800, res = 150)
plot(dbobject)
dev.off()

dbobject<-dba.normalize(dbobject)
norm<-dba.normalize(dbobject, bRetrieve=TRUE)
norm

normlibs<-cbind(FullLibSize=norm$lib.sizes,NormFacs=norm$norm.factors,NormLibSize=round(norm$lib.sizes/norm$norm.factors))
rownames(normlibs)<-info$ID
normlibs

dbobject<-dba.contrast(dbobject,reorderMeta=list(Condition=refcond))
dbobject

dbobject<-dba.analyze(dbobject)
dba.show(dbobject,bContrasts=TRUE)

png(paste(outpath,"/diffbind_CorrHeatmap_byDiffPeaks.png",sep=""), width = 1000, height = 800, res = 150)
try(plot(dbobject,contrast=1))
dev.off()

dbobject.DB<-dba.report(dbobject,contrast=1)
outdifpeaks<-as.data.frame(dbobject.DB)
write.table(outdifpeaks,file=paste(outpath,"/DifferentialPeaks_",compared,"_vs_",refcond,"_ref.txt",sep=""),sep="\t",quote=F,row.names=F)

sum(dbobject.DB$Fold>0)
sum(dbobject.DB$Fold<0)

png(paste(outpath,"/diffbind_pca_byNormCounts.png",sep=""), width = 1000, height = 800, res = 150)
dba.plotPCA(dbobject,DBA_CONDITION,label=DBA_CONDITION)
dev.off()

png(paste(outpath,"/diffbind_pca_byDiffPeaks.png",sep=""), width = 1000, height = 800, res = 150)
try(dba.plotPCA(dbobject,contrast=1,label=DBA_CONDITION))
dev.off()

png(paste(outpath,"/diffbind_volcano_byDiffPeaks.png",sep=""), width = 1000, height = 800, res = 150)
dba.plotVolcano(dbobject)
dev.off()

png(paste(outpath,"/diffbind_boxplot_byDiffPeaks.png",sep=""), width = 1000, height = 800, res = 150)
try(dba.plotBox(dbobject))
dev.off()

pvals<-try(dba.plotBox(dbobject))
try(write.table(pvals,file=paste(outpath,"/pvals_boxplot_",compared,"_vs_",refcond,"_ref.txt",sep=""),sep="\t",quote=F,row.names=T))

hmap<-colorRampPalette(c("red","black","green"))(n=13)

png(paste(outpath,"/diffbind_SampleHeatmap_byDiffPeaks.png",sep=""), width = 1000, height = 800, res = 150)
try(dba.plotHeatmap(dbobject,contrast=1,correlations=FALSE,scale="row",colScheme=hmap))
dev.off()

profiles<-try(dba.plotProfile(dbobject))

png(paste(outpath,"/diffbind_ProfileHeatmap_bySampleContrast.png",sep=""),width=1000,height=800,res=150)
try(dba.plotProfile(profiles))
dev.off()

profiles<-try(dba.plotProfile(dbobject,merge=c(DBA_REPLICATE)))

png(paste(outpath,"/diffbind_ProfileHeatmap_byMergedContrast.png",sep=""),width=1000,height=800,res=150)
try(dba.plotProfile(profiles))
dev.off()

