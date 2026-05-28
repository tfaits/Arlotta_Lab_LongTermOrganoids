#!/usr/bin/bash

### function to call DMRs between two organoid samples
call_DMRs() {
  local sample_1=$1
  local sample_2=$2

  # Create input table for metilene using bedtools unionbedg
  bedtools unionbedg -i data/beds/${sample_1}_hg19_CpG_CovChromFiltered.sorted.bed \
                         data/beds/${sample_2}_hg19_CpG_CovChromFiltered.sorted.bed \
                      -header -filler NA -names ${sample_1} ${sample_2} \
                      | cut -f1,3- | sed 's/end/pos/' \
                      > results/tables/metilene/unionbedg_methRates_${sample_1}_${sample_2}.tab

  # Call metilene
  metilene_v0.2-8/metilene_linux64 --mincpgs 10 --threads 4 --groupA "${sample_1}" --groupB "${sample_2}" \
                                                                          results/tables/metilene/unionbedg_methRates_${sample_1}_${sample_2}.tab 2>&1 > \
                                                                          results/tables/metilene/metilene_${sample_1}_${sample_2}.out | grep -v -i segmenting > results/tables/metilene/metilene_${sample_1}_${sample_2}.log

  # Filter significant DMRs and keep important columns
  less results/tables/metilene/metilene_${sample_1}_${sample_2}.out | perl -ane 'if($F[3]ne"NA" && $F[3]<0.05){print $_}' | cut -f1-6,9,10 | bedtools sort > results/tables/metilene/metilene_${sample_1}_${sample_2}_DMRs.bed

  # Split into hypo- and hyperDMRs
  less results/tables/metilene/metilene_${sample_1}_${sample_2}_DMRs.bed | perl -ane 'if($F[4]>0){print $_}' > results/tables/metilene/metilene_${sample_1}_${sample_2}_hypoDMRs.bed
  less results/tables/metilene/metilene_${sample_1}_${sample_2}_DMRs.bed | perl -ane 'if($F[4]<0){print $_}' > results/tables/metilene/metilene_${sample_1}_${sample_2}_hyperDMRs.bed
}

# call DMRS between timepoint 3 months vs. 5 years
call_DMRs V1 V9
