library("DESeq2")

library(dplyr)
library(ggplot2)
#library(EnhancedVolcano)

args = commandArgs(trailingOnly=TRUE)

outpath <- args[3]
refcond <- args[4]
compared <- args[5]

#count <- read.delim(args[1],sep="",row.names = 1,check.names=FALSE)
count <- read.delim(args[1],sep="",row.names = NULL,check.names=FALSE)
modsum <- as.formula(paste(".~",colnames(count)[1],sep=""))
count <- aggregate(modsum,count,sum)
row.names(count) <- count[,1]
count[,1] <- NULL

design <- read.delim(args[2],header=T,sep="\t",row.names = 1,check.names=FALSE)
design <- design %>% filter_all(any_vars(. %in% c(refcond,compared)))
design <- design[1:(ncol(design)-1)]

overlap<-intersect(row.names(design),colnames(count))

count <- count[,overlap]
design_header<-colnames(design)
design <- as.data.frame(design[overlap,])
row.names(design)<-overlap
colnames(design)<-design_header

for (i in 1:length(colnames(design))){
    design[,i] <-factor(design[,i])
}

for (i in 1:length(colnames(design))){
  if (i == 1){
    model <- paste("~",colnames(design)[i])
  } else {
    model <- paste(model,"+",colnames(design)[i])
  }
}

formula <- as.formula(model)

dds <- DESeqDataSetFromMatrix(countData = count,colData = design,design = formula)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

### PCA Plot #####
for (i in 1:length(colnames(design))){

png(paste(outpath,'/pca_',colnames(design)[i],'.png',sep=""),units="in", width=5, height=5, res=300,type="cairo")
vsdata <- vst(dds, blind=FALSE)
pca<-plotPCA(vsdata, intgroup=colnames(design)[i])
pca<-pca+geom_text(label=rownames(design),color="black",nudge_x=1,nudge_y=0.7,size=2)
print(pca)
dev.off()
  
}
##################

for (i in 1:length(colnames(design))){
    if (refcond %in% design[,colnames(design)[i]]){
        dds@colData@listData[[i]]<-relevel(dds@colData@listData[[i]], ref = refcond)
        pheno_col<-colnames(design)[i]
    }
}

dds <- DESeq(dds)

normalized_counts<-as.data.frame(counts(dds, normalized=TRUE))
normalized_counts$DESCRIPTION <- 'na'
normalized_counts<-normalized_counts[,c(ncol(normalized_counts),1:ncol(normalized_counts)-1)]

#write.table(normalized_counts, file=paste(outpath,"/normalized_counts.txt",sep=""), sep="\t", quote=F, col.names=NA)
write.table(normalized_counts, file=paste(outpath,"/normalized_counts_",compared,"_vs_",refcond,"(ref).txt",sep=""), sep="\t", quote=F, col.names=NA)

resultsNames(dds)

res <- results(dds,name=paste(pheno_col,compared,"vs",refcond,sep="_"))
resLFC <- lfcShrink(dds, coef=paste(pheno_col,compared,"vs",refcond,sep="_"), type="apeglm")

### Volcano Plot (enhanced) ###
##png(paste(outpath,'/volcano_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=7, height=5, res=300)
##EnhancedVolcano(resLFC,
##                lab = rownames(resLFC),
##                x = 'log2FoldChange',
##                y = 'pvalue',
##                pCutoff = 1e-4,
##                FCcutoff = 1,
##                pointSize = 2.0,
##                labSize = 3.0)
##dev.off()
###########

######## MA Plot ###########
png(paste(outpath,'/MAplot_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=5, height=5, res=300,type="cairo")
plotMA(resLFC, ylim=c(-2,2))
dev.off()
############################

resOrdered <- resLFC[order(resLFC$pvalue),]
summary(resOrdered)

resOrdered <- as.data.frame(resOrdered)
resOrdered <- na.omit(resOrdered)

res_sig <- resOrdered[resOrdered$padj<0.05,]

write.table(res_sig,file=paste(outpath,"/DEG_fdr0.05_",compared,"_vs_",refcond,"(ref).txt",sep=""),sep="\t",quote=F) 
write.table(resOrdered,file=paste(outpath,"/DEG_",compared,"_vs_",refcond,"(ref).txt",sep=""),sep="\t",quote=F)
