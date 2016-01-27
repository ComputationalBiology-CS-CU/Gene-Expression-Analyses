#### Script adapted from tutorial on matrixEQTL website. FDR filter added

# R script to conduct cis eQTL analysis with MatrixEQTL package
# Necessary to format the input files according to specifications
# given in the website: http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
# Argument parameters for script call are as follows:
# 1) genotype
# 2) snp location
# 3) gene expression
# 4) gene location
# 5) Output file name
# 6) covariates (optional)
#
#
# Below are a list of set parameters that can be changed accordingly

args <- commandArgs(TRUE)

testModel = "linear" # Can be "ANOVA", or linear_cross as well
cis_pValueThres = 1 # Set to 0 to ignore SNP and gene locations
trans_pValueThresh = 0 # Set to 0 for cis eQTL only
cisDist = 1e6 # Distance to include cis SNPs
fdrThreshold = 0.01

# Load package
library("MatrixEQTL")


# Set Model
useModel = modelLINEAR

# Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = " "
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 0
snps$fileSkipColumns = 5
snps$fileSliceSize = 2000
snps$LoadFile(args[1])

# Load the gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"
gene$fileOmitCharacters = "NA"
gene$fileSkipRows = 1
gene$fileSkipColumns = 1
gene$fileSliceSize = 2000
gene$LoadFile(args[3])

# Load the covariates data
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"
cvrt$fileOmitCharacters = "NA"
cvrt$fileSkipRows = 1
cvrt$fileSkipColumns = 1
cvrt$fileSliceSize = 2000
cvrt$LoadFile(args[6])


snpspos = read.table(args[2], header=TRUE, stringsAsFactors = FALSE)
genepos = read.table(args[4], header=TRUE, stringsAsFactors = FALSE)

me = Matrix_eQTL_main(
	snps = snps,
	gene = gene,
	cvrt = cvrt,
	output_file_name = NULL,
	pvOutputThreshold = trans_pValueThresh,
	useModel = useModel,
	errorCovariance = numeric(),
	verbose = TRUE,
	output_file_name.cis = NULL,
	pvOutputThreshold.cis = cis_pValueThres,
	snpspos = snpspos,
	genepos = genepos,
	cisDist = cisDist,
	pvalue.hist = FALSE,
	min.pv.by.genesnp = FALSE,
	noFDRsaveMemory = FALSE);

#Plot
plot(me)

# Results
cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n')
cat('Detected local eQTLs:', '\n') 
final <- me$cis$eqtls[me$cis$eqtls[5] < fdrThreshold ,]
write.table(final,file=args[5],quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
final <- me$trans$eqtls[me$trans$eqtls[5] < fdrThreshold ,]
write.table(final,file=args[7],quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)
