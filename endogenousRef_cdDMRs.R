# Load necessary libraries
library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(readxl)
library(tidyr)
require(RColorBrewer)


## Generates:
## - Ext Data Fig 1i
## - Ext Data Fig 1j

## Requires: bigWigAverageOverBed (UCSC tools) run upstream to produce
##   per-sample .tab files (see meanMethylation_HMDs_PMDs_DMVs.sh for a similar example) 

########################
### LOAD SAMPLE INFO ###
########################

sample_info <- read_excel("data/sample_info/brain_organoids_WGBS.xlsx")

# numeric time in culture (months)
Time_in_months <- c(
  "3 month" = 3,  "6 month" = 6,  "9 month" = 9,
  "1 year"  = 12, "1.5 year" = 18, "2 year"  = 24,
  "3 year"  = 36, "4 year"   = 48, "5 year"  = 60
)

sample_ids   <- sample_info$`Sample ID`
time_months  <- Time_in_months[sample_info$`Time in culture`]
names(time_months) <- sample_ids

########################
### LOAD cdDMRs ########
########################


# Load endogenous DMRs and their clusters from paper
endg_DMRs_df <- read_excel("resources/endogenousRef/cdDMRs.xlsx", sheet = 1)
endg_DMRs <- GRanges(seqnames = endg_DMRs_df$Chromosome,
                        ranges = IRanges(start = endg_DMRs_df$Start, end = endg_DMRs_df$End),
                        value = endg_DMRs_df$Value,
                        cluster = endg_DMRs_df$Cluster)

DMRs_adj <- import("resources/endogenousRef/endogenous_DMRs_adj.bed") # adjusted bed file with numbered regions to match bigWigAverageOverBed output 

# match regions in DMRs_adj with endg_DMRs_df
DMRs_adj$cluster <- endg_DMRs[match(DMRs_adj, endg_DMRs)]$cluster
DMRs_adj$cdDMR_value <- endg_DMRs[match(DMRs_adj, endg_DMRs)]$value

#############################################
### LOAD bigWigAverageOverBed OUTPUT ########
#############################################

# Load results from UCSC bigWigAverageOverBed 
load_ucsc_mean_bigwig <- function(sample_info, region_set_name, region_gr, mean_tab_file_folder){

  # create data frame to gather means from all samples 
  regions_scored <- data.frame(matrix(NA, nrow = length(region_gr), ncol = nrow(sample_info)))
  colnames(regions_scored) <- sample_info$`Sample ID` 
  rownames(regions_scored) <- region_gr$name

  # iterate over all samples and read results from UCSC tools
  for(sample in 1:nrow(sample_info)){
      sample_regions_scored <- read.table(paste0(mean_tab_file_folder, sample_info[sample,]$`Sample ID` , "_", region_set_name, "_hg19.tab"), 
                                  col.names = c("name", "size", "covered", "sum", "mean0", "mean")) 


      regions_scored[sample_regions_scored$name, sample_info[sample,]$`Sample ID`] <- sample_regions_scored$mean
  }
  return(regions_scored)
}

# read DMRs scored
DMRs_scored <- load_ucsc_mean_bigwig(sample_info = sample_info, region_set_name = "endogenousDMRs", region_gr = DMRs_adj, mean_tab_file_folder = "results/tables/endogenousDMRs/")

DMRs_scored$Cluster <- DMRs_adj$cluster[match(rownames(DMRs_scored), DMRs_adj$name)]


#####################################################################
##### COMPUTE PER-DMR LINEAR REGRESSION SLOPE OVER TIME #############
##################################################################### 

dmr_all_slopes <- data.frame()

DMRs_mat <- DMRs_scored[, c(sample_ids, "Cluster"), drop = FALSE]


for (cluster in sort(unique(DMRs_mat$Cluster))) {
  mat <- DMRs_mat[DMRs_mat$Cluster == cluster, sample_ids, drop = FALSE]

  slopes <- apply(mat, 1, function(y) {
    t <- time_months[colnames(mat)]
    ok <- !is.na(y) & !is.na(t)
    if (sum(ok) < 3) return(NA_real_)
    coef(lm(y[ok] ~ t[ok]))[2]  # slope per month
  })

  dmr_all_slopes <- rbind(dmr_all_slopes,
                          data.frame(Cluster = cluster,
                                     DMR = names(slopes),
                                     slope_per_month = slopes,
                                     stringsAsFactors = FALSE))
}

dmr_all_slopes <- dmr_all_slopes %>%
  dplyr::mutate(Direction = ifelse(slope_per_month > 0, "Hyper", "Hypo"))

# add coordinates 
dmr_all_slopes$region <- as.character(DMRs_adj[match(dmr_all_slopes$DMR, DMRs_adj$name),])


################################################################################
## SUMMARISE FRACTION OF DMRS EXCEEDING SLOPE THRESHOLD PER CLUSTER ############
################################################################################


# threshold: |slope| > 0.1 / 60 corresponds to ≥ 0.1 methylation change over
# 60 months (5 years), which is the span of the V1–V9 time series
X <- 0.1 / 60   # per-month slope threshold


cluster_colors <- brewer.pal(6, "Dark2")
names(cluster_colors) <- sort(unique(DMRs_adj$cluster))

cluster_frac <- dmr_all_slopes %>%
  dplyr::group_by(Cluster) %>%
  dplyr::summarize(
    Num_DMRs = sum(!is.na(slope_per_month)),
    Num_big  = sum(abs(slope_per_month) > X, na.rm = TRUE),
    Num_big_down  = sum((slope_per_month) < -X, na.rm = TRUE),
    Num_big_up  = sum((slope_per_month) > X, na.rm = TRUE),
    Fraction = Num_big / Num_DMRs,
    Fraction_down = Num_big_down/ Num_DMRs,
    Fraction_up = Num_big_up/ Num_DMRs,
    .groups = "drop"
  ) %>%
  dplyr::mutate(Cluster = factor(Cluster, levels = names(cluster_colors)))
  

################################################################################
## BARPLOT: FRACTION HYPER / HYPO PER CLUSTER ##################################
################################################################################


df_div <- cluster_frac %>%
  mutate(Up = Fraction_up,
         Down = -Fraction_down) %>%
  select(Cluster, Up, Down) %>%
  pivot_longer(cols = c(Up, Down), names_to = "Direction", values_to = "Fraction") 

ggplot(df_div, aes(x = Cluster, y = Fraction, fill = Cluster)) +
  geom_col(width = 0.8) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  theme_minimal() +
  theme(panel.background = element_blank(),
        axis.line = element_line(colour = "black", linewidth = rel(0.5)),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(x = "Cluster", y = "Fraction (Up / Down)", fill = NULL) +
  scale_fill_manual(values = cluster_colors) +
  coord_cartesian(ylim = c(-1, 1))

ggsave("fraction_DMRs_with_slopes_up_down_per_cluster_V1_to_V9.pdf", width = 6, height = 4)

################################################################################
## BOXPLOT OF SLOPES, SPLIT BY DIRECTION (HYPER/HYPO) AND CLUSTER ##############
################################################################################

ggplot(dmr_all_slopes %>% dplyr::filter(!is.na(slope_per_month)),
       aes(x = Cluster, y = slope_per_month, color = Cluster)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3) +
  geom_jitter(width = 0.15, alpha = 0.25, size = 0.5) +
  scale_color_manual(values = cluster_colors) +
  facet_wrap(~ Direction, nrow = 1, scales = "free_x") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "none") +
  labs(x = NULL, y = "Slope per month (methylation ~ time)")
  
ggsave("slope_per_month_by_cluster_hyper_hypo.pdf", width = 8, height = 3.5)

