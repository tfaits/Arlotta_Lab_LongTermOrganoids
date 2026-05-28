
library(tidyverse)
library(rtracklayer)
library(data.table)
library(xlsx)

## Epigenetic clock analysis
## Applies three DNA methylation-based age prediction clocks to WGBS data imputed to EPIC array probe positions.
##
## Clocks applied:
##   1. Multiple reference clocks via the methylclock R package
##      (Horvath)
##   2. Cortical Clock
##   3. Fetal Brain Clock
##
## Input:  boostme_imputed_data_merged.RData (from imputation_boostMe.R)
##         EPIC.hg19.manifest.tsv (probe coordinate annotation)




## READ DATA

# SAMPLE INFO
sample_info <- read.xlsx("data/sample_info/brain_organoids_WGBS.xlsx", sheetIndex = 1, check.names = F)


# EPIC PROBE ID ANNOTATION
EPIC_anno <- read.table("references/hg19/annotation/EPIC.hg19.manifest.tsv", sep = "\t", header = TRUE)
# remove NA rows
EPIC_anno <- EPIC_anno[!is.na(EPIC_anno$CpG_chrm), ]
# probe ID positions
EPIC_anno$epic_positions <- paste0(EPIC_anno$CpG_chrm, "-", EPIC_anno$CpG_beg, "-", EPIC_anno$CpG_end)


######################################################
### METHYLATION CLOCKS USING METHYLCLOCK R PACKAGE ###
######################################################

library(methylclock)


## CREATE DATA FRAME WITH EPIC IDS

# read imputed data
load(file = "results/tables/methylation_clock/boostme_imputed_data_merged.RData")
imputed_data_merged[which(imputed_data_merged > 1, arr.ind = T)] <- 1

# prepare dataframe
colnames(imputed_data_merged) <- gsub("\\.", "-", colnames(imputed_data_merged))
colnames(imputed_data_merged) <- sample_info$`Sample ID`[match(colnames(imputed_data_merged), sample_info$`WGBS library ID`)]


# subset imputed data to epic positions
probes_scored <- imputed_data_merged[rownames(imputed_data_merged) %in% EPIC_anno$epic_positions,]

# translate genomic position to probe ID
probes_scored$ProbeID = EPIC_anno[match(rownames(probes_scored), EPIC_anno$epic_positions),]$Probe_ID

saveRDS(probes_scored, 'results/tables/merged_methRates_wProbeIDs.rds')


## APPLY METHYLATION CLOCKS USING METHYLCLOCK R PACKAGE

predicted_ages_tbl <- NULL

# iterate over all samples and apply clock to each sample
for(sample in sample_info$`Sample ID`){
    print(paste("Start sample", sample))

    probes_scored_sample <- probes_scored[, c('ProbeID', sample, sample)] # select the column of the sample twice due to the input format of the methylclock function
    probes_scored_sample <- probes_scored_sample[rowSums(is.na(probes_scored_sample)) == 0,]

    prediced_ages <- DNAmAge(probes_scored_sample)


    print(paste0("Actual age: ", sample_info[sample_info$`Sample ID` == sample,]$`Time in culture`))
    print("Predicted ages:")
    print(prediced_ages)

    predicted_ages_tbl <- dplyr::bind_rows(predicted_ages_tbl, prediced_ages[1,])
}

saveRDS(predicted_ages_tbl, 'results/tables/refClocks/predictedAges_methylclock.rds')


######################################################
### CORTICAL CLOCK USING PROVIDED R FUNCTION #########
######################################################

source('resources/corticalClock/CorticalClock.r')


# apply cortical clock

# CorticalClock<-function(betas, ## betas = betas matrix (rownames=cpgs, colnames=IDs)
#                         pheno, ##  pheno file = file which contains IDs that match betas IDs and  contains actual Age col
#                         dir, ## directory where the ref file and coeffecients are saved
#                         IDcol, ## ID column which matches Betas IDs for your samples
#                         Agecol){  ## Age column

## prepare inputs

# betas
betas_matrix <- as.matrix(probes_scored[, !names(probes_scored) %in% "ProbeID"])
rownames(betas_matrix) <- probes_scored$ProbeID

# prepare pheno object for cortical clock

# translate time in culture to numeric values in months
sample_info$time_in_culture_numeric <- 0
sample_info[sample_info$`Time in culture` == "15 days",  ]$time_in_culture_numeric <- 0.5
sample_info[sample_info$`Time in culture` == "1 month", ]$time_in_culture_numeric <- 1
sample_info[sample_info$`Time in culture` == "3 month",  ]$time_in_culture_numeric <- 3
sample_info[sample_info$`Time in culture` == "6 month",  ]$time_in_culture_numeric <- 6
sample_info[sample_info$`Time in culture` == "9 month",  ]$time_in_culture_numeric <- 9
sample_info[sample_info$`Time in culture` == "1 year",   ]$time_in_culture_numeric <- 12
sample_info[sample_info$`Time in culture` == "1.5 year", ]$time_in_culture_numeric <- 18
sample_info[sample_info$`Time in culture` == "2 year",   ]$time_in_culture_numeric <- 24
sample_info[sample_info$`Time in culture` == "3 year",   ]$time_in_culture_numeric <- 36
sample_info[sample_info$`Time in culture` == "4 year",   ]$time_in_culture_numeric <- 48
sample_info[sample_info$`Time in culture` == "5 year",   ]$time_in_culture_numeric <- 60
sample_info[sample_info$`Time in culture` == "6 year",   ]$time_in_culture_numeric <- 72

# pheno object for cortical clock
pheno <- sample_info[,c("Sample ID", "time_in_culture_numeric")]
pheno$`time_in_culture_numeric` <- pheno$`time_in_culture_numeric`/12


# call cortical clock predictions
CorticalClock(betas = betas_matrix, pheno = pheno, dir = 'resources/corticalClock/', IDcol = "Sample ID", Agecol = "time_in_culture_numeric")

# results are saved as "CorticalPred.csv"

######################################################
### FETAL CLOCK USING PROVIDED R FUNCTION ############
######################################################


# The fetal brain clock function expects a matrix of beta values with CpGs as rows and samples as columns, and rownames as CpG IDs.
# This function transforms the data frame we have (with ProbeID as an additional column) into the required format.
make_betas <- function(ps) {
  m <- as.matrix(ps[, !names(ps) %in% "ProbeID"])
  rownames(m) <- ps$ProbeID
  m
}

betas_fetal <- make_betas(probes_scored)

## LOAD & APPLY FETAL BRAIN CLOCK

source('resources/fetalClock/FetalClockFunction.R')

# The FetalClock function returns predictions in days post-conception (DPC) by default. We convert to gestational weeks (GW) by dividing by 7.
reshape_pred <- function(pred_mat) {
  df <- as.data.frame(pred_mat)             # already samples × 1
  colnames(df) <- "pred_dpc"
  df$ID      <- rownames(df)
  df$pred_GW <- df$pred_dpc / 7
  df
}

# Apply the clock
fetal_pred_df <- reshape_pred(FetalClock(betas = betas_fetal, ageinyears = FALSE,
                                        dir   = 'resources/fetalClock/'))

# Print summary of predictions
cat("Fetal clock predictions (GW) — all samples:\n")
print(summary(fetal_pred_df$pred_GW))
