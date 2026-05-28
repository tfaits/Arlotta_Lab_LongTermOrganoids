
#!/usr/bin/bash


## bigWigAverageOverBed (UCSC tools) was used to compute mean methylation from CpG-resolution bigWig files over genomic feature BED files

## this is an examplary script to extract average methylation levels across PMDs/HMDs and DMVs 
## the same approach was followed to calculate average methylation levels across other features: bins, repeats, cdDMRs, super-enhancers


#############################
## paths
BW_DIR=data/bigwigs/
OUT_DIR=results/tables/HMD_PMD_DMV

#############################
UCSC_DIR=src/UCSCtools


## list of samples to process
SAMPLES=(
    V1
    V2
    V3
    V4
    V5
    V6
    V7
    V8
    V9
)


for smpl in "${SAMPLES[@]}"
do
    echo "Starting sample ${smpl}"

    ${UCSC_DIR}/./bigWigAverageOverBed ${BW_DIR}/${smpl}_hg19_CpG_CovChromFiltered.bw \
    data/annotation/DMVs_hg19.bed \
    ${OUT_DIR}/${smpl}_DMV_hg19.tab

    ${UCSC_DIR}/./bigWigAverageOverBed ${BW_DIR}/${smpl}_hg19_CpG_CovChromFiltered.bw \
    data/annotation/PMDs_hg19.bed \
    ${OUT_DIR}/${smpl}_PMD_hg19.tab

    ${UCSC_DIR}/./bigWigAverageOverBed ${BW_DIR}/${smpl}_hg19_CpG_CovChromFiltered.bw \
    data/annotation/HMDs_hg19.bed \
    ${OUT_DIR}/${smpl}_HMD_hg19.tab

done



