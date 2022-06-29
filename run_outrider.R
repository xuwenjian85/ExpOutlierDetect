# --------------------------------------------
# simulation dataset
# randomly select corrupt target(gene, sample) pairs at frequency P_BENCH, 
# artificially inject with log(count+1) scaled raw matrix by shift +Z/-Zscore change
# input: a mask_matrix with -1/+1 injected
# output: counts injected according to mask
## input
args = commandArgs(trailingOnly=TRUE)
print(args)
ods_simu <- args[1]
q0 <- as.numeric(args[2])
ods_simu_outrider <- args[3]

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(dplyr)
})

setDTthreads(1)
ods <- readRDS(ods_simu)
ods <- OUTRIDER(ods, q = q0, BPPARAM = SerialParam(),iterations=8)

saveRDS(ods, ods_simu_outrider)

