




#Read in the training dataset labels for each tissue
train_samples <- read.table("../tmp/train_samples.txt",sep="\t",header=FALSE,row.names=1)
row.names(train_samples) <- lapply(row.names(train_samples),function(x){gsub("-",".",x,fixed=TRUE)})
colnames(train_samples) <- "Tissue"
test_samples <- read.table("../tmp/test_samples.txt",sep="\t",header=FALSE,row.names=1)
row.names(test_samples) <- lapply(row.names(test_samples),function(x){gsub("-",".",x,fixed=TRUE)})
colnames(test_samples) <- "Tissue"

#Read in expr dataset
expr <- read.table("../GTEx_processed/GTEx_Data_2014-01-17_RNA-seq_RNA-SeQCv1.1.8_gene_rpkm.gct_processed_2_gene_normalized",header=TRUE,skip=2,stringsAsFactors=FALSE)
#Separate expr

train_expr <- expr[,colnames(expr) %in% c('Name',row.names(train_samples))]
test_expr <- expr[,colnames(expr) %in% c('Name',row.names(test_samples))]

train_expr = setNames(data.frame(t(train_expr[,-1])), train_expr[,1])
test_expr = setNames(data.frame(t(test_expr[,-1])), test_expr[,1])

train_expr <- data.frame(train_expr,train_samples[match(row.names(train_expr),row.names(train_samples)),"Tissue"])
test_expr <- data.frame(test_expr,test_samples[match(row.names(test_expr),row.names(test_samples)),"Tissue"])

colnames(train_expr)[length(colnames(train_expr))] <- "Tissue"
colnames(test_expr)[length(colnames(test_expr))] <- "Tissue"

Read in tss sites for each gene
tss <- read.table("../GTEx_processed/gencode.v18.genes.patched_contigs.gtf_gene_tss",header=FALSE)

Create data-frame for each gene and run model fitting
chrList <- vector()
corrList <- vector()
posList <- vector()

for(gene in expr[[1]]){
	gene.tss <- tss$V3[tss$V1 == gene]	
	gene.chr <- tss$V2[tss$V1 == gene]
	if(gene.chr == "MT" || gene.chr == "X" || gene.chr == "Y"){	
		next
	}
	chr.name <- paste("chr",gene.chr,sep="")	
	file_name = paste("../GTEx_processed/genotype_185_dosage_matrix_qc/",chr.name,"/genotypes_all.txt",sep="")
	geno <- read.table(file_name,sep="\t",header=TRUE)
	geno <- geno[geno[,"POS"] < (gene.tss + 1e6) & geno[,"POS"] > (gene.tss - 1e6),]
	curr_table <- train_expr[,c(gene,"Tissue")]
	ind <- sapply(row.names(curr_table),function(x){substr(x,1,9)})
	geno_train <- geno[,ind]
	curr_table <- data.frame(curr_table,t(geno_train))
	colnames(curr_table)[1] <- "Gene"	
	model <- lm(Gene ~ . ,,data=curr_table)

	}

	curr_test <- test_expr[c(gene,"Tissue")]	
	ind2  <- sapply(row.names(test_expr),function(x){substr(x,1,9)})
	geno_test <- geno[,ind2]		
	colnames(curr_test)[1] <- "Gene"
	p <- predict(model,newdata=data.frame(curr_test,t(geno_test)))
	corrList <- c(corrList,c(cor(p,test_expr[,gene])[1]))
	
	chrList <- c(chrList,gene.chr)
	posList <- c(posList,gene.tss)
	cat(sprintf("%d\t%d\t%f\n",gene.chr,gene.tss,cor(p,test_expr[,gene])))
	

}

final <- data.frame(chrList,posList,corrList)
colnames(final) <- c("CHR","BP","P")
final <- split(final,final$CHR)

j = 0
for(f in final){
	j <- j+1
	plot(f$BP,f$P,xlab="Position",ylab="Correlation",main=paste("Chromosome",j,sep=" "))
}


