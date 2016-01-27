#Converts phenotype files into covariate files for eQTL anaylyses. The covariates are separated by affection status

e <- read.table("../CommonMind/Clinical/CMC_MSSM-Penn-Pitt_Clinical.csv",header=TRUE,sep=",")
e <- e[,c("DLPFC_RNA_Sequencing_Sample_ID","Gender","Ethnicity","pH","Age_of_Death","Dx")]
e$Gender <- as.numeric(e$Gender)
e$Ethnicity <- as.numeric(e$Ethnicity)
control <- e[e$Dx == "Control",]
scz <- e[e$Dx == "SCZ",]
bp <- e[e$Dx == "BP",]
aff <- e[e$Dx == "AFF",]

e <- e[,c("DLPFC_RNA_Sequencing_Sample_ID","Gender","Ethnicity","pH","Age_of_Death")]
control <- control[,c("DLPFC_RNA_Sequencing_Sample_ID","Gender","Ethnicity","pH","Age_of_Death")]
scz <- scz[,c("DLPFC_RNA_Sequencing_Sample_ID","Gender","Ethnicity","pH","Age_of_Death")]
bp <- bp[,c("DLPFC_RNA_Sequencing_Sample_ID","Gender","Ethnicity","pH","Age_of_Death")]
aff <- aff[,c("DLPFC_RNA_Sequencing_Sample_ID","Gender","Ethnicity","pH","Age_of_Death")]

write.table(as.data.frame(t(e)),"covariates_ALL.tsv",sep="\t",quote=FALSE,col.names=FALSE)
write.table(as.data.frame(t(scz)),"covariates_SCZ.tsv",sep="\t",quote=FALSE,col.names=FALSE)
write.table(as.data.frame(t(bp)),"covariates_BP.tsv",sep="\t",quote=FALSE,col.names=FALSE)
write.table(as.data.frame(t(aff)),"covariates_AFF.tsv",sep="\t",quote=FALSE,col.names=FALSE)



