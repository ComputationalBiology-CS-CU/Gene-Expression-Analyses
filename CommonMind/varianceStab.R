#### Script used to perform the variance stabilizing (affine log) transform on the gene expression data

library('DESeq2')
library("vsn")

d <- read.table("../GTEx_processed/GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized",header=TRUE,skip=2,stringsAsFactors=FALSE)
n <- d$Name
desc <- d$Description
d$Name <- NULL
d$Description <- NULL 
colData <- data.frame(V1=sample(x=c(1,2),replace=TRUE,size=ncol(d)),row.names=colnames(d))
dds <- as.matrix(d)


#dds <- DESeqDataSetFromMatrix(countData=dds, colData=colData, design=formula(~1))

t <- as.data.frame(log2(dds+1))
t <- data.frame(n,desc,t)
rownames(t)[1] <- "Name"
rownames(t)[2] <- "Description"
write.table(t,"./GTEx_counts_transformed.tsv",quote=FALSE,sep="\t",row.names=FALSE)

par(mfrow=c(1,2))
meanSdPlot(dds, ylim = c(0,2.5),main="Original expression")
meanSdPlot(log2(dds+1), ylim = c(0,2.5), main="log2(counts + 1)")

#t <- as.data.frame(log2(dds+1))
#rownames(t) <- rownames(d)
#colnames(t) <- colnames(d)




sprintf("Plotted log2 chart")
sprintf("Plotted original chart")

