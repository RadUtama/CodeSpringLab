library("DESeq2")

library(dplyr)
library(ggplot2)
#library(EnhancedVolcano)
library(ggrepel)

args = commandArgs(trailingOnly=TRUE)

outpath <- args[3]
refcond <- args[4]
compared <- args[5]
redundant <- args[6]

#count <- read.delim(args[1],sep="",row.names = 1,check.names=FALSE)
count <- read.delim(args[1],sep="",row.names = NULL,check.names=FALSE)
modsum <- as.formula(paste(".~",colnames(count)[1],sep=""))
count <- aggregate(modsum,count,sum)
row.names(count) <- count[,1]
count[,1] <- NULL

design <- read.delim(args[2],header=T,sep="\t",row.names = 1,check.names=FALSE)
if (redundant %in% colnames(design)){
    design <- design[,-match(redundant, colnames(design))]
}
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

png(paste(outpath,'/pca_',colnames(design)[i],'_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=5, height=5, res=300,type="cairo")
vsdata <- vst(dds, blind=FALSE)
pca<-plotPCA(vsdata, intgroup=colnames(design)[i])+theme(aspect.ratio = 1)
pca<-pca+geom_text_repel(label=rownames(design),color="black",size=2)+labs(title = "geom_text_repel()")
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
png(paste(outpath,'/volcano_',compared,"_vs_",refcond,'(ref).png',sep=""),units="in", width=7, height=5, res=300, ,type="cairo")

##EnhancedVolcano(resLFC,
##                lab = rownames(resLFC),
##                x = 'log2FoldChange',
##                y = 'pvalue',
##                pCutoff = 1e-4,
##                FCcutoff = 1,
##                pointSize = 2.0,
##                labSize = 3.0)

resLFC_data<-data.frame(resLFC)
resLFC_data$diffexpressed <- "NO"
resLFC_data$diffexpressed[resLFC_data$log2FoldChange > 0 & resLFC_data$padj < 0.05] <- "UP"
resLFC_data$diffexpressed[resLFC_data$log2FoldChange < 0 & resLFC_data$padj < 0.05] <- "DOWN"
resLFC_data$delabel <- NA
resLFC_data <- resLFC_data[order(resLFC_data$padj),]
resLFC_data$delabel[1:50]<-rownames(resLFC_data)[1:50]

ggplot(data = resLFC_data, aes(x = log2FoldChange, y = -log10(pvalue),col = diffexpressed,label=delabel))+geom_vline(xintercept = c(-1, 1), col = "red", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.0001), col = "red", linetype = 'dashed') + geom_point(size=0.5)+
  scale_color_manual(values = c("darkblue", "grey", "darkred"),
                     labels = c("Downregulated", "Not significant", "Upregulated"))+
  geom_text_repel(max.overlaps = Inf,color="black",size=1)

dev.off()
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
