## input
args = commandArgs(trailingOnly=TRUE)
print(args)
dname=args[1]
ctsFile = args[2]
q = as.integer(args[3])
out1 = args[4]
out2 = args[5]

start <- Sys.time()
suppressPackageStartupMessages({
    library(OUTRIDER)
    library(dplyr)
})
## load data
# ctsFile <- system.file('extdata', 'KremerNBaderSmall.tsv', package='OUTRIDER')
# ctsFile <-  '/media/eys/xwj/RNAseq/public_normal/df_cts_HC1157_corrupt_fc2_ngene100_nrep3.txt'
ctsTable <- read.table(ctsFile, check.names = FALSE)
ctsTable <- ctsTable[1:500, (ncol(ctsTable)-600+1):ncol(ctsTable)]

ods <- OutriderDataSet(countData=ctsTable)
ods <- filterExpression(ods, minCounts=TRUE, filterGenes=TRUE,)
ods <- estimateSizeFactors(ods)
ods <- findEncodingDim(ods, params = q, BPPARAM = SerialParam(),iterations=8)
# ods <- findEncodingDim(ods, params = q, implementation='pca',BPPARAM =MulticoreParam(1))

write.table(ods@metadata$encDimTable, file=out1, sep='\t',quote = FALSE)
saveRDS(ods, file=out2)

end <- Sys.time()
print(end-start)