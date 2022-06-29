# --------------------------------------------
# simulation dataset
# randomly select corrupt target(gene, sample) pairs at frequency P_BENCH, 
# artificially inject with log(count+1) scaled raw matrix by shift +Z/-Zscore change
# output: a mask_matrix with -1/+1 injected
print(paste('freq is', p))
print(paste('Zscore is', zScore))
print('input ods is')
print(ods)
# --------------------------------------------
# load samples with no aberrant genes
rawtable <- counts(ods)
print(counts(ods)[1:5,1:5])

#inject outliers
size <- dim(rawtable)[1] * (dim(rawtable)[2])

index_mat <- matrix(nrow=nrow(rawtable), 
        data=sample(c(0,1,-1), size, prob=c(1 - p, p/2, p/2), replace=TRUE))
print(dim(index_mat))

# create outrider object
mask_ods <- OutriderDataSet(countData=counts(ods), colData=colData(ods))
sizeFactors(mask_ods) <- sizeFactors(ods) 
# thetaCorrection(mask_ods) <- NULL
metadata(mask_ods)$optimalEncDim <- metadata(ods)$optimalEncDim
metadata(mask_ods)$encDimTable <- metadata(ods)$encDimTable
metadata(mask_ods)$dim <- dim(mask_ods)

counts(mask_ods) <- rawtable
print(mask_ods)

assay(mask_ods, "inj_mask", withDimnames=FALSE) <- index_mat

saveRDS(mask_ods, out_true_outliers)
