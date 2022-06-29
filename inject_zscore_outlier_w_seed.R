# --------------------------------------------
# simulation dataset
# randomly select corrupt target(gene, sample) pairs at frequency P_BENCH, 
# artificially inject with log(count+1) scaled raw matrix by shift +Z/-Zscore change
# input: a mask_matrix with -1/+1 injected
# output: counts injected according to mask
## input
args = commandArgs(trailingOnly=TRUE)
print(args)

in_true_outliers <- args[1]
zScore <- as.numeric(args[2])
out_inj_counts <- args[3]
inj <- 'both' ## high & low 

suppressPackageStartupMessages({
    library(OUTRIDER)
    library(dplyr)
})

setDTthreads(1)
ods <- readRDS(in_true_outliers)
k <- counts(ods)
index_mat <- assay(ods, "inj_mask")
size <- nrow(k) * ncol(k)
print(counts(ods)[1:3,1:5])

# --------------------------------------------
# check injection mode

if(inj=='low'){
    index_mat <- -abs(index_mat)
}
if(inj=='high'){
    index_mat <- abs(index_mat)
}

# --------------------------------------------
# inject outliers
 
# estimate from data if not simulated
sf <- DESeq2::estimateSizeFactorsForMatrix(k)
normtable <- t(t(k)/sf)
datasd <- matrix(rowSds(log2(normtable + 1)), nrow=nrow(ods), ncol=ncol(ods), byrow=FALSE)
lmu <- matrix(rowMeans(log2(normtable+1)), nrow=nrow(ods), ncol=ncol(ods), byrow=FALSE)

# use simulated values if present
if(all(c("true_mean", "true_sd") %in% assayNames(ods))){
    sf <- colData(ods)[['true_sizeFactor']]
    datasd <- assay(ods, "true_sd")
    lmu <- log2(assay(ods, "true_mean"))
}

# inject outliers
max_out <- 1E2 * min(max(k), .Machine$integer.max/1E3)
n_rejected <- 0
list_index <- which(index_mat != 0, arr.ind = TRUE)
for(i in seq_len(nrow(list_index))){
    row <- list_index[i,'row']
    col <- list_index[i,'col']
    fc <- zScore * datasd[row,col]
    clcount <- index_mat[row,col] * fc + lmu[row,col]
    
    #multiply size factor again
    art_out <- round(sf[col]*2^clcount)
    if(art_out < max_out){
        k[row,col] <- art_out
    }else{
        #remove super large outliers
        index_mat[row,col] <- 0 
        n_rejected <- n_rejected + 1
    }
}

print(paste('n_rejected is', n_rejected))
mode(k) <- "integer"

# create outrider object
ods_out <- copy(ods)
assay(ods_out, 'trueCounts') <- counts(ods)
assay(ods_out, 'inj_mask' )  <- index_mat # super large outliers may become
print(table(index_mat)) 
counts(ods_out) <- k

#save tables; save rawtable which is now the injected one.
saveRDS(ods_out, out_inj_counts)

