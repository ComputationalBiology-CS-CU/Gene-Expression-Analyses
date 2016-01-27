library(peer)

library(qtl)

args = commandArgs(TRUE)

if(length(args) < 1) {
	cat("Usage: ./peer.R <genotype.file> <phenotype.file> <covariates.file>") 
}

cross <- read.cross(format="csvs",genfile=args[1],phefile=args[2],genotypes=c(0,1,2),sep="\t")

covariates <- read.table(args[3],skip=1,sep="\t",row.names=TRUE)

print("Read Tables")

model = PEER()

print("Created Model")

PEER_setNk(model,100)

print("Set Max Factors")

PEER_setPhenoMean(model, as.matrix(cross$pheno))

print("Set Pheno Mean")

PEER_setCovariates(model, as.matrix(covariates))

print("Set Covariates")

PEER_setNMax_iterations(model,250)

print("Set Max Iter")

PEER_update(model)

print("Finished Converging")

plot(precision)
