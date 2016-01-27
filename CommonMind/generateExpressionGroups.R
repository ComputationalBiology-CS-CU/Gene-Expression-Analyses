#Script to create expression files for each affection status for context specific eqtl analyses

e <- read.table("../CommonMind/RNAseq/DorsolateralPrefrontalCortex/QuantitatedExpression/CMC_MSSM-Penn-Pitt_DLPFC_mRNA_IlluminaHiSeq2500_geneExpressionRaw.tsv",header=TRUE,sep="\t",row.names=1)

d <- read.table("../CommonMind/Clinical/disease_phenotypes.tsv",header=TRUE)

bp <- d$DLPFC_RNA_Sequencing_Sample_ID[d$Dx == "BP"]
scz <- d$DLPFC_RNA_Sequencing_Sample_ID[d$Dx == "SCZ"]
aff <- d$DLPFC_RNA_Sequencing_Sample_ID[d$Dx == "AFF"]

bp_e <- e[,bp]
scz_e <- e[,intersect(scz,colnames(e))]
aff_e <- e[,intersect(aff,colnames(e))]

write.table(bp_e,"./expr_BP.tsv",quote=FALSE,sep="\t")
write.table(scz_e,"./expr_SCZ.tsv",quote=FALSE,sep="\t")
write.table(aff_e,"./expr_AFF.tsv",quote=FALSE,sep="\t")




