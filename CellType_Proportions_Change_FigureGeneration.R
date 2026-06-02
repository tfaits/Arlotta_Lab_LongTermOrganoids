require(lme4)

# convert significance values into asterisks
makeAsterisk <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

media_cols <- c(
  'CDM4' = '#00B8A9',
  'Bphys' = '#F6416C',
  'APM' = '#F6416C'
)

umap.colors <- c(
  'aRG' = "#8dd3c7",
  'oRG' = "#f9e79d",
  'IP'= "#bebada",
  'PN' = "#fcbba1",
  'IN progenitors' = "#F0C4D0",
  'CFuPN' = "#80b1d3",
  'CPN' = "#fdb462",
  'Immature IN' = "#fa9fb5",
  'tRG' = "#9BD19B",
  'Glial precursors' = "#43A85E",
  'Astrocytes' = "#005A40",
  'OPC' = "#c3c483",
  'ChPL' = "#e1ce9a"
)

# use NBEM method to evaluate cell type proportional changes ----
# load meta data of age-specific integrated object from CDM4 & APM (BrainPhys/Bphys-based) scRNA-seq data (i.e. 6mo)
meta.6mo <- readRDS("6mo_harmonized_Metadata.RDS")

# generate contingency table of: Organoid x BroadType
myTab.all <- unclass(table(meta.6mo$Organoid, meta.6mo$BroadType))
myTab.all <- data.frame(myTab.all)
colnames(myTab.all) <- str_replace_all(colnames(myTab.all), "\\.", " ")
CellType.all <- colnames(myTab.all)

# Use the total cell count per organoid as an offset; add columns of treatment and genotype from sample ids
myTab.all$LibSize <- apply(myTab.all, 1, sum)
myTab.all$Treatment <- factor(sapply(strsplit(rownames(myTab.all), "_"), "[", 1), levels = c("CDM4", "Bphys"))
myTab.all$Genotype <- factor(sapply(strsplit(rownames(myTab.all), "_"), "[", 4))

# subset based on the shared genotypes between CDM4 and APM scRNA-seq data
genotype.tb <- table(myTab.all$Treatment, myTab.all$Genotype)
shared_genotypes <- colnames(genotype.tb)[apply(genotype.tb, 2, function(x) all(x > 0))]
myTab <- myTab.all[myTab.all$Genotype %in% shared_genotypes,]
myTab <- droplevels(myTab)

# create output directory
date <- str_replace_all(Sys.Date(), "-", "")
outDir <- paste0("visualization_df/CellTypeProp/", date, "/")
if(!dir.exists(outDir)){
  dir.create(outDir, recursive = TRUE)
}

Bphys.stats <- data.frame(CellType = character(),
                          pvalue = numeric())

print("Bphys vs. CDM4:")
anova.res <- list()
RandMod1 <- list()
TreatMod1 <- list()
for (i in 1:length(CellType.all)) {
  print(CellType.all[i])
  myTab$TypeOfInterest <- myTab[[CellType.all[i]]]
  
  # if the cell type is found in both conditions, proceed
  if (all(table(myTab$TypeOfInterest > 0, myTab$Treatment)['TRUE', ] > 0) == TRUE) {
    # if any sample has zero cell count in TypeOfInterest column, add small constant 1 to every TypeOfInterest. Each LibSize will increment by 1 accordingly
    if (any(myTab$TypeOfInterest == 0)) {
      myTab$TypeOfInterest <- myTab$TypeOfInterest + 1
      myTab$LibSize <- myTab$LibSize + 1
    }
    # First attempt at fitting the RandomModel and TreatmentModel
    tryCatch({
      RandomModel <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + (1|Genotype), data=myTab)
      TreatmentModel <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Genotype), data=myTab)
    }, error = function(e) {
      # If an error occurs, define control parameters and try fitting again
      cat("Encountered error. Refitting models with custom control parameters...\n")
      control_params <- glmerControl(optCtrl = list(maxfun = 10000), check.conv.grad = .makeCC("warning", tol = 3e-3))
      
      RandomModel <- glmer.nb(TypeOfInterest ~ offset(log(LibSize)) + (1|Genotype), data = myTab, control = control_params)
      TreatmentModel <- glmer.nb(TypeOfInterest ~ Treatment + offset(log(LibSize)) + (1|Genotype), data = myTab, control = control_params)
    })
    
    anova.res[[CellType.all[i]]] <- anova(RandomModel, TreatmentModel)
    Bphys.stats[i, ] <- c(CellType.all[i], scientific(anova.res[[CellType.all[i]]][2, 8]))
    RandMod1[[CellType.all[i]]] <- RandomModel
    TreatMod1[[CellType.all[i]]] <- TreatmentModel
  } else {
    print(paste0(CellType.all[[i]], " is only found in one condition, skip for testing.\n"))
    Bphys.stats[i, ] <- c(CellType.all[i], NA)
    next
  }
}
Bphys.stats$p.adjust <- scientific(p.adjust(Bphys.stats$pvalue, method="BH"))
Bphys.stats$asterisk <- makeAsterisk(Bphys.stats$p.adjust)
Bphys.stats

save(RandMod1, TreatMod1, anova.res, Bphys.stats, file = paste0(outDir, "6mo_NBME_Bphys_vs_CDM4.RData"))

# figure generation ----
df.all <- table(droplevels(meta.6mo$Organoid[meta.6mo$Genotype %in% shared_genotypes]),
                droplevels(meta.6mo$BroadType[meta.6mo$Genotype %in% shared_genotypes]))
df.all <- (df.all / apply(df.all, 1, sum)) * 100
df.all <- reshape2::melt(df.all)
colnames(df.all) <- c("Organoid", "Cell Type", "Percent of Cells")
df.all$Treatment <- sapply(strsplit(as.character(df.all$Organoid), "_"), "[", 1)
df.all$Treatment <- factor(df.all$Treatment, levels = c("CDM4", "Bphys"))
df.all
glimpse(df.all)

# assign broader categories to cell types for visual arrangements
df.all$category <- ""
df.all$category[df.all$`Cell Type` %in% c("aRG", "oRG", "IP", "PN", "IN progenitors")] <- "Category 1"
df.all$category[df.all$`Cell Type` %in% c("CFuPN", "CPN", "Immature IN")] <- "Category 2"
df.all$category[df.all$`Cell Type` %in% c("tRG", "Glial precursors", "Astrocytes", "OPC", "ChPL")] <- "Category 3"

df.all$`Cell Type` <- factor(df.all$`Cell Type`, levels = c("aRG", "oRG", "IP", "PN", "IN progenitors",
                                                            "CFuPN", "CPN", "Immature IN",
                                                            "tRG", "Glial precursors", "Astrocytes", "OPC", "ChPL"))

mean_df.all <-  df.all %>% 
  group_by(category, `Cell Type`, Treatment) %>%
  dplyr::summarise(mean= mean(`Percent of Cells`),
                   se = sd(`Percent of Cells`)/sqrt(n()))
mean_df.all
save(RandMod1, TreatMod1, Bphys.stats, df.all, mean_df.all, file = paste0(outDir, "6mo_NBME_Bphys_vs_CDM4.RData"))

# load statistics
stats1 <- Bphys.stats

# iterate and generate sub-barplots for each category
barplots <- list()
for (case in unique(mean_df.all$category)) {
  df <- df.all[df.all$category == case,]
  df$`Cell Type` <- droplevels(df$`Cell Type`)
  
  mean_df <- mean_df.all[mean_df.all$category == case,]
  mean_df$`Cell Type` <- droplevels(mean_df$`Cell Type`)
  
  strip <- strip_themed(background_x = elem_list_rect(fill = umap.colors[levels(mean_df$`Cell Type`)]),
                        text_x = elem_list_text(colour = c("white"), face = "bold", size = 10))
  
  composition.barplot <- ggplot(mean_df) +
    geom_bar(aes(x = `Treatment`, y = mean, fill = `Treatment`), stat="identity", alpha=0.8) +
    scale_fill_manual(values = media_cols) +
    geom_errorbar( aes(x = `Treatment`, ymin=mean-se, ymax=mean+se), width=0.2, colour="#6C6C6C", linewidth=.5) +
    geom_point(aes(x = `Treatment`, y = `Percent of Cells`, color = `Treatment`), size = 1.5, alpha=0.8, data = df, position = position_jitter(seed = 42,width=0.15, height=0)) +
    scale_color_manual(values = media_cols) +
    facet_grid2(. ~ `Cell Type`, strip = strip, scales="free_x") + 
    theme_classic() + 
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks=element_blank(),
          axis.ticks.x.bottom = element_blank(),
          axis.title.y = element_text(size = 12),
          panel.grid.major = element_blank(), 
          strip.background = element_rect(colour=NA, fill=NA),
          panel.spacing.x = unit(0, "null"),
          panel.border = element_rect(fill = NA, color = "white"),
          legend.title = element_text(size=13), #change legend title font size
          legend.text = element_text(size=12)) +
    labs(y = "Percent of Cells") +
    guides(color = "none") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(df$`Percent of Cells`) + 10)) 
  
  # for adding significance notations on grouped barplots at corresponding positions;
  # one way is to replace corresponding statistical values from a Wilcox_test tibble (used for anchor positions of stats) with the statistics
  stats.pos <- df %>%
    group_by(`Cell Type`) %>%
    wilcox_test(`Percent of Cells` ~ `Treatment`, ref.group = "CDM4") %>%
    adjust_pvalue(method = "BH") %>%
    add_significance() %>% 
    add_xy_position(x = "Treatment")
  
  for (col in c("statistic", "p", "p.adj", "p.adj.signif")) {
    stats.pos[[col]] <- NULL
  }
  
  stats.pos <- left_join(stats.pos, 
                                    stats1, by = c('Cell Type' = 'CellType'))
  stats.pos$`Cell Type` <- factor(stats.pos$`Cell Type`,
                                             levels = levels(mean_df$`Cell Type`))
  stats.pos$y.position <- stats.pos$y.position - 1
  
  composition.barplot <- composition.barplot + stat_pvalue_manual(stats.pos[stats.pos$asterisk != "ns",],
                                                                                        hide.ns = TRUE,
                                                                                        label = "{asterisk}",
                                                                                        label.size = 8,
                                                                                        bracket.nudge.y = -4,
                                                                                        remove.bracket = T)
  legend <- get_legend(composition.barplot)
  composition.barplot <- composition.barplot + NoLegend()
  
  barplots[[case]] <- composition.barplot
}

barplots.grid <- plot_grid(plotlist =  barplots, ncol = 3, rel_widths = c(.5, .3, .5))
combined.barplot <- plot_grid(barplots.grid, legend, ncol = 2, rel_widths = c(1, .1))
combined.barplot

# save the plot in svg format
ggsave(paste0("6mo_CellType_compositions_comparison_barplots_Bphys_vs_CDM4_with_signif_", str_replace_all(Sys.Date(), "-", ""), ".svg"), combined.barplot, path = outDir, width = 16.5, height = 4.4)





