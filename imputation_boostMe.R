library(data.table)
library(rtracklayer)
library(bsseq)
library(boostme)
library(readxl)


## CpG methylation imputation of methylation values at low-coverage CpGs using boostme 
## 
## Input:  per-sample CpG-resolution BED files (*hg19_CpG_CovChromFiltered.sorted.bed)
## Output: boostme_imputed_data_merged.RData
##         — merged data frame of imputed methylation values across all samples,
##           used as input for epigenetic clock analysis (agePrediction_methylclocks.R)


################################################################################
## Impute CpG methylation with low coverage
################################################################################

# load data
input_files <- list.files(pattern = "hg19_CpG", full.names = TRUE)
names(input_files) <- sapply(input_files, function(x) strsplit(basename(x), "_hg19")[[1]][1])

input_bsseq <- lapply(names(input_files), function(x)
{
    print(x)
    data <- data.frame(fread(input_files[x], header = FALSE), stringsAsFactors = FALSE)
    colnames(data) <- c("chr", "start", "end", "cov", "meth")
    data$start <- data$start + 1
    data$end <- data$end - 1

    data_gr <- makeGRangesFromDataFrame(data[,1:3])
    meth <- as.matrix(data[,"meth",drop=FALSE])
    cov <- as.matrix(data[,"cov",drop=FALSE])
    colnames(meth) <- NULL
    colnames(cov) <- NULL

    bs <- BSseq(M = meth, Cov = cov, gr = data_gr, sampleNames = c(x))
    return(bs)
})

set.seed(42)

# impute methylation values with boostme 
imputed_data <- lapply(input_bsseq, function(x) return(boostme(x, minCov = 5, sampleAvg = FALSE, threads = 16)))

# convert imputed data to data frame and merge
imputed_data <- lapply(imputed_data, function(x)
{
    data <- data.frame(x)
    data$region <- rownames(data)
    return(data)
})

imputed_data_merged <- Reduce(function(x,y) merge(x = x, y = y, all.x = TRUE, all.y = TRUE), imputed_data)

coordinate_conversion <- data.frame(row.names = as.character(imputed_data_merged$region))
coordinate_conversion$chr <- as.character(sapply(rownames(coordinate_conversion), function(x) strsplit(x, ":")[[1]][1]))
coordinate_conversion$start <- sapply(rownames(coordinate_conversion), function(x) as.numeric(strsplit(x, ":")[[1]][2]) - 1)
coordinate_conversion$end <- sapply(rownames(coordinate_conversion), function(x) as.numeric(strsplit(x, ":")[[1]][2]) + 1)
coordinate_conversion$cpg <- paste(coordinate_conversion$chr, coordinate_conversion$start, coordinate_conversion$end, sep = "-")

rownames(imputed_data_merged) <- coordinate_conversion[as.character(imputed_data_merged$region),"cpg"]
imputed_data_merged$region <- NULL

save(imputed_data_merged, file = "results/tables/methylation_clock/boostme_imputed_data_merged.RData")
