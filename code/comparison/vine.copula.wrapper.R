# Extract arguments
require(methods)
args = commandArgs(trailingOnly=TRUE)
if(length(args) < 1) { argsFile = "data/args.csv" } else { argsFile = args[1] }
if(length(args) < 2) { dataFile = "data/U.csv" } else { dataFile = args[2] }
if(length(args) < 3) { rvmFile = "data/rvm-params.csv" } else { rvmFile = args[3] }
if(length(args) < 4) { sampleFile = "data/samplesUt.csv" } else { sampleFile = args[4] }

# Load params
argsTable <- read.table(argsFile,
                 header = FALSE,
                 sep = ",")
nSamples = as.integer(argsTable[1,1])
vineType = as.integer(argsTable[1,2])
indTest = as.logical(argsTable[1,3])
pLevel = as.double(argsTable[1,4])
nWorkers = as.logical(argsTable[1,5])
if(NCOL(argsTable) > 5) {
    familySet = as.vector(argsTable[1,-(1:5)])
} else {
    familySet = NA
}
show(familySet)

# Load data
Xtable <- read.table(dataFile,
                 header = FALSE,
                 sep = ",")
Xt <- as.matrix(t(Xtable))

# Load library and source files for XMRF
require(VineCopula)

# Run program
RVM <- RVineStructureSelect(Xt, type=vineType, indeptest=indTest, level=pLevel, familyset=familySet, cores=nWorkers)
UtSample <- RVineSim( nSamples, RVM )

# Save results
write.table(rbind(RVM$Matrix, RVM$family, RVM$par, RVM$par2), rvmFile, row.names=FALSE, col.names=FALSE, sep=",")
write.table(UtSample, sampleFile, row.names=FALSE, col.names=FALSE, sep=",")
