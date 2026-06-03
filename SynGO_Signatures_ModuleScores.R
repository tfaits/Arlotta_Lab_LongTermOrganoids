require(Seurat)
require(qs)
require(lme4)
require(readxl)
require(dplyr)
require(plyr)
require(stringr)
require(reshape2)
require(rstatix)
require(openxlsx)
require(ggplot2)
require(ggpubr)
require(ggnewscale)

# convert significance values into asterisks
makeAsterisk <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

# color hex codes
media.border.cols <- c('CDM4' = "#5fb5c9",
                       'APM' = "#c593c2")
media.box.cols <- c('CDM4' = "#1a7c93",
                    'APM' = "#912883")
inner.genotype.cols <- c("H1" = "#e8d7a5", "PGP1" = "#ee9d23", "Mito210" = "#c86b29", "11a" = "#931007")
outer.genotype.cols <- c("H1" = "#ac9457", "PGP1" = "#ad6e29", "Mito210" = "#89421e", "11a" = "#410000")


# load signature gene sets retrieved from SynGO knowledge base ----
SynGO.ontologies <- read_excel("SynGO_bulk_download_release_20231201/syngo_ontologies.xlsx")

synapse.modules.all <- strsplit(SynGO.ontologies[SynGO.ontologies$name == 'synapse (GO:0045202)', ]$hgnc_symbol, ", ")[[1]]
presynapse.modules.all <- strsplit(SynGO.ontologies[SynGO.ontologies$name == 'presynapse (GO:0098793)', ]$hgnc_symbol, ", ")[[1]]
postsynapse.modules.all <- strsplit(SynGO.ontologies[SynGO.ontologies$name == 'postsynapse (GO:0098794)', ]$hgnc_symbol, ", ")[[1]]

# load age-specific integrated object from CDM4 & APM (BrainPhys/Bphys-based) scRNA-seq data (i.e. 6mo) ----
harmonized.6mo <- qread("v5_6mo_harmonized.qs")
levels(harmonized.6mo$Treatment) # should print out CDM4, Bphys in order

# run AddModuleScore ----
harmonized.6mo <- AddModuleScore(
  object = harmonized.6mo,
  features = list(synapse.modules.all),
  name = "Synapse.all"
)
harmonized.6mo <- AddModuleScore(
  object = harmonized.6mo,
  features = list(presynapse.modules.all),
  name = "Presynapse.all"
)
harmonized.6mo <- AddModuleScore(
  object = harmonized.6mo,
  features = list(postsynapse.modules.all),
  name = "Postsynapse.all"
)

harmonized.6mo[['Synapse.module']] <- harmonized.6mo$Synapse.all1
harmonized.6mo[['Presynapse.module']] <- harmonized.6mo$Presynapse.all1
harmonized.6mo[['Postsynapse.module']] <- harmonized.6mo$Postsynapse.all1
harmonized.6mo$Synapse.all1<- NULL
harmonized.6mo$Presynapse.all1 <- NULL
harmonized.6mo$Postsynapse.all1 <- NULL
qsave(harmonized.6mo, "v5_6mo_harmonized.qs")

# subset object based on the shared genotypes between 2 treatment conditions
genotype.tb <- table(harmonized.6mo$Treatment, harmonized.6mo$Genotype)
shared_genotypes <- colnames(genotype.tb)[apply(genotype.tb, 2, function(x) all(x > 0))]
sub1 <- harmonized.6mo@meta.data[harmonized.6mo$Genotype %in% shared_genotypes, ]
sub1 <- droplevels(sub1)

# arrange cell types into broader populations
lineages <- list('Neurons' = c("CFuPN", "CPN", "Migrating CPN", "PN", "Immature IN"),
                 'Neurogenic Progenitors' = c("aRG", "IN progenitors", "IP", "oRG", "tRG"),
                 'Astrocytes' = c('Astrocytes'),
                 'OPC' = c("OPC")
)

# generate graph data & run likelihood ratio test comparing linear mixed-effects models ----
date <- str_replace_all(Sys.Date(), "-", "")
outDir <- paste0("visualization_df/ModuleScores/", date, "/")
if(!dir.exists(outDir)){
  dir.create(outDir, recursive = TRUE)
}

age <- "6mo"
dfm.all <- list()
mean.dfm.all <- list()
lmer.stats <- list()
rm.list <- list()
tm.list <- list()
stats.betweenTreatments <- list()
for (lineage in names(lineages)) {
  print(paste0(lineage, " starting..."))
  lmer.stats[[lineage]] <- c()
  
  df.all <- sub1
  df.all <- df.all[df.all$BroadType %in% lineages[[lineage]], ]
  df.all <- droplevels(df.all)
  
  fileConn <- file(paste0(outDir, age, "_", lineage, "_lmer_summary.txt"), open = "wt")
  if (length(levels(df.all$Treatment)) > 1) {
    for (moduleOfInterest in c('Synapse.module', 'Presynapse.module', 'Postsynapse.module')) {
      print(paste0(moduleOfInterest, " starting..."))
      df <- df.all[c("Treatment", "Organoid", "Genotype", "BroadType", moduleOfInterest)]
      
      for (column in colnames(df)) {
        if (column != moduleOfInterest) {
          df[[column]] <- as.factor(df[[column]])
        }
      }
      dfm = reshape2::melt(df)
      dfm$Lineage <- lineage
      dfm$Module <- moduleOfInterest
      
      # average the module scores by each cell type per organoid
      mean.dfm <- dfm %>%
        group_by(Lineage, Module, Treatment, Organoid, Genotype) %>%
        dplyr::summarise(mean.module.score = mean(`value`)) %>%
        ungroup()
      
      dfm.all[[lineage]] <- rbind(dfm.all[[lineage]], dfm)
      mean.dfm.all[[lineage]] <- rbind(mean.dfm.all[[lineage]], mean.dfm)
      
      # run Mann-Whitney U test (only to provide a template dataframe to anchor positions of stats on the plot)
      stats.pos <- mean.dfm %>%
        wilcox_test(mean.module.score ~ Treatment) %>%
        add_xy_position(x = "Treatment")
      stats.pos$y.position <- stats.pos$y.position + 0.1
      stats.pos$Module <- moduleOfInterest
      stats.betweenTreatments[[lineage]] <- rbind(stats.betweenTreatments[[lineage]], stats.pos)
      
      # run lmer tests
      # use bobyqa optimizer to deal potential convergence warnings or when random effects variance is close to zero
      dfm <- dfm %>% 
        group_by(Organoid) %>% 
        mutate(obs_count = n()) %>% ungroup()
      dfm$weight_factor <- 1 / sqrt(dfm$obs_count)
      
      result <- tryCatch({
        random.model.weighted <- lmer(value ~ (1|Organoid) + (1|Genotype), data = dfm, weights = weight_factor, control = lmerControl(optimizer = "bobyqa"))
        treat.model.weighted <- lmer(value ~ Treatment + (1|Organoid) + (1|Genotype), data = dfm, weights = weight_factor, control = lmerControl(optimizer = "bobyqa"))
        rm.list[[paste0(lineage, " - ", moduleOfInterest, ": RandomModel")]] <- random.model.weighted
        tm.list[[paste0(lineage, " - ", moduleOfInterest, ": TreatmentModel")]] <- treat.model.weighted
        
        coef <- fixef(treat.model.weighted)['TreatmentBphys']
        OR <- exp(coef)
        CI <- confint(treat.model.weighted, parm="TreatmentBphys", method="Wald")
        CI_OR <- exp(CI)
        anov <- anova(random.model.weighted, treat.model.weighted)  
        pv <- anov$"Pr(>Chisq)"[2]
        
        row <- c(age, lineage, moduleOfInterest, "with weight_factor",
                 "value ~ (1|Organoid) + (1|Genotype)", "value ~ Treatment + (1|Organoid) + (1|Genotype)",
                 as.numeric(c(coef,OR,CI,CI_OR,pv)))
        lmer.stats[[lineage]] <- rbind(lmer.stats[[lineage]], row)
        
        writeLines(paste("\n\n### ", age, " ", lineage, " ", moduleOfInterest, " Formula: value ~ (1|Organoid) + (1|Genotype) (with weight_factor) ###"), fileConn)
        sink(fileConn, append = TRUE)
        cat("\n### Random Model Summary: ###\n")
        print(summary(random.model.weighted))
        
        cat("\n### Treatment Model Summary: ###\n")
        print(summary(treat.model.weighted))
        
        cat("\n### ANOVA Results ###\n")
        print(anov)
        
        writeLines("\n--------------------------------------------------------\n", fileConn)
      }, warning = function(w) {
        sink(fileConn, append = TRUE)
        cat(paste0(lineage, " encountered warnings:\n"))
        cat("Warning caught:", w$message, "\n")
        writeLines("\n--------------------------------------------------------\n", fileConn)
        return(NULL)
      }, error = function(e) {
        sink(fileConn, append = TRUE)
        cat(paste0(lineage, " encountered errors:\n"))
        cat("Error caught:", e$message, "\n")
        writeLines("\n--------------------------------------------------------\n", fileConn)
        return(NULL)
      })
    }
    
    # adjusting p values after multiple independent tests
    lmer.stats[[lineage]] <- data.frame(lmer.stats[[lineage]])
    colnames(lmer.stats[[lineage]]) <- c("Age", "Lineage", "Module", "weight_factor?", "rm.formula", "tm.formula", "coef","OR","CI_coef_low","CI_coef_high","CI_OR_low","CI_OR_high","pval")
    rownames(lmer.stats[[lineage]]) <- seq(1, nrow(lmer.stats[[lineage]]))
    
    lmer.stats[[lineage]] <- lmer.stats[[lineage]] %>%
      group_by(rm.formula) %>%
      mutate(p.adj = p.adjust(pval, method = "BH")) %>%
      ungroup()
    
    lmer.stats[[lineage]]$p.adj.signif <- makeAsterisk(lmer.stats[[lineage]]$p.adj)
    
    dfm.all[[lineage]]$Module <- factor(dfm.all[[lineage]]$Module, levels = c("Synapse.module", "Presynapse.module", "Postsynapse.module"))
    mean.dfm.all[[lineage]]$Module <- factor(mean.dfm.all[[lineage]]$Module, levels = c("Synapse.module", "Presynapse.module", "Postsynapse.module"))
    lmer.stats[[lineage]]$Module <- factor(lmer.stats[[lineage]]$Module, levels = c("Synapse.module", "Presynapse.module", "Postsynapse.module"))
    
    stats.betweenTreatments[[lineage]]$p.adj <- p.adjust(stats.betweenTreatments[[lineage]]$p, method = "BH")
    stats.betweenTreatments[[lineage]]$p.adj.signif <- makeAsterisk(stats.betweenTreatments[[lineage]]$p.adj)
    stats.betweenTreatments[[lineage]]$Module <- factor(stats.betweenTreatments[[lineage]]$Module, levels = c("Synapse.module", "Presynapse.module", "Postsynapse.module"))
  } else {
    stop(paste0("cannot run on ", lineage, " : observations only found in one level of Treatment."))
  }
  sink()
  close(fileConn)
}
save(dfm.all, mean.dfm.all, lmer.stats, rm.list, tm.list, stats.betweenTreatments, file = paste0(outDir, age, "_lmer.RData"))
write.xlsx(lmer.stats, paste0(outDir, age, "_lmer_results.xlsx"))

# draw violin plots by 2 layers (for easier loading in Adobe Illustrator)
vlnplots.6mo.background  <- list()
vlnplots.6mo.violin  <- list()
signif.df <- list() # stores rstatix stats' output df to anchor stats' positions; with statistical values replaced by lmer results

for (lineage in names(dfm.all)) {
  dfm <- dfm.all[[lineage]]
  mean.dfm <- mean.dfm.all[[lineage]]
  dfm$Treatment <- factor(mapvalues(as.character(dfm$Treatment), from = "Bphys", to = "APM"), levels = c("CDM4", "APM"))
  mean.dfm$Treatment <- factor(mapvalues(as.character(mean.dfm$Treatment), from = "Bphys", to = "APM"), levels = c("CDM4", "APM"))
  stats1 <- lmer.stats[[lineage]]
  
  # format positions where stats will be added
  stats1 <- stats1[, c('Module', 'p.adj', 'p.adj.signif')]
  stats1 <- join(stats.betweenTreatments[[lineage]][c(".y.", "group1", "group2", "n1", "n2", "y.position", "xmin", "xmax", "Module")], stats1, by = 'Module')
  stats1$group2 <- mapvalues(stats1$group2, from = "Bphys", to = "APM")
  signif.df[[lineage]] <- stats1
  
  p.background <- ggplot(dfm, aes(y=value, x = Treatment)) +
    geom_point(aes(x = `Treatment`, y = `value`), size = 0.1, data = dfm, position = position_jitter(seed = 42 , width = 0.2), color = '#DCDCDC') + 
    theme_bw() +
    facet_grid(. ~ Module) +
    ggtitle(lineage) +
    theme(plot.title = element_text(hjust = .5, size = 15),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(colour=NA, fill=NA),
          strip.text = element_text(size = 13),
          panel.border = element_rect(fill = NA, color = "white"),
          axis.text = element_text(size = 13),
          axis.title.y = element_text(size = 15)
    ) + 
    labs(y = "Module Score") +
    ylim(0, NA) +
    stat_pvalue_manual(stats1[stats1$p.adj.signif != "ns",],
                       hide.ns = T,
                       label = "{p.adj.signif}",
                       label.size = 8) +
    guides(fill = "none")
  
  p.violin <- ggplot(dfm, aes(y=value, x = Treatment)) +
    geom_violin(aes(color = Treatment), fill = "#f9f9f9", linewidth = 1, width = 0.6) + 
    scale_color_manual(values = scales::alpha(media.border.cols, 0.8)) +
    new_scale_color() +
    guides(
      color = guide_legend(order = 1)
    ) +
    geom_boxplot(aes(x = Treatment, y = mean.module.score, color = Treatment), outlier.shape = NA, data = mean.dfm, fill = NA, width = 0.25, linewidth = 1) +
    scale_color_manual(values = media.box.cols) +
    new_scale_color() +
    geom_point(aes(x = Treatment, y = mean.module.score, color = Genotype, fill = Genotype), 
               size = 2.5, shape = 21, stroke = .5,
               data = mean.dfm, position = position_jitter(seed = 42,width=0.1, height=0)) +
    scale_fill_manual(values = inner.genotype.cols) +
    scale_color_manual(values = outer.genotype.cols) +
    theme_bw() +
    facet_grid(. ~ Module) +
    ggtitle(lineage) +
    theme(plot.title = element_text(hjust = .5, size = 15),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x.bottom = element_blank(),
          panel.grid.major.x = element_blank(),
          strip.background = element_rect(colour=NA, fill=NA),
          strip.text = element_text(size = 13),
          panel.border = element_rect(fill = NA, color = "white"),
          axis.text = element_text(size = 13),
          axis.title.y = element_text(size = 15)
    ) + 
    labs(y = "Module Score") +
    ylim(0, NA) +
    stat_pvalue_manual(stats1[stats1$p.adj.signif != "ns",],
                       hide.ns = T,
                       label = "{p.adj.signif}",
                       label.size = 8)
  
  p.violin.legend <- get_legend(p.violin)
  p.violin <- p.violin + NoLegend()
  
  vlnplots.6mo.background[[lineage]] <- p.background
  vlnplots.6mo.violin[[lineage]] <- p.violin
  
  # export 2 layers into distinct svgs
  fileName <- paste0(age, "_SynGO_module_in_", lineage, "_", str_replace_all(Sys.Date(), "-", ""), ".svg")
  ggsave(paste0("Background_", fileName), p.background, path = outDir, width = 8, height = 4.3, device = "svg")
  ggsave(paste0("Violin_", fileName), p.violin, path = outDir, width = 8, height = 4.3, device = "svg")
}





