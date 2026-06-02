library(openxlsx)
library(ggplot2)
library(rstatix)
library(stats)
`%nin%` = Negate(`%in%`)
`%!in%` = Negate(`%in%`)

outDir <- "MEA_res/"

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = TRUE)
}

# load MEA data  ----
prism.PGP1 <- "PGP1_PRISM_data.xlsx"
prism.11a <- "11a_PRISM_data.xlsx"
prism.H1 <- "H1_PRISM_data.xlsx"

df.PGP1 <- read.xlsx(prism.PGP1)
df.11a <- read.xlsx(prism.11a)
df.H1 <- read.xlsx(prism.H1)

# filter out outliers
df.PGP1 <- df.PGP1[is.na(df.PGP1$Outlier) == TRUE, ] # only 1
df.11a <- df.11a[is.na(df.11a$Outlier) == TRUE, ] # only 1

# measures used in main figure panels
# note: H1 doesn't have measures for SpikeNumber and NeuronNumber
main.measures <- c("SpikeRate", "SpikeNumber", "NeuronNumber", "NetBurstRate")
# measures used in supplementary figure panels
suppl.measures <- c("ActiveNeurons", "BurstAmplitude", "BurstDuration", "BurstRate", "NetBurstDuration", "FanoFactor")

# Main ----
df.main <- rbind(df.H1[df.H1$Measurement %in% main.measures,], 
                 df.PGP1[df.PGP1$Measurement %in% main.measures,], 
                 df.11a[df.11a$Measurement %in% main.measures,])

df.main$Treatment <- factor(df.main$Treatment, levels = c("CDM4", "APM"))
df.main$Genotype <- factor(df.main$Genotype, levels = c("H1", "PGP1", "11a"))
df.main$Age <- factor(df.main$Age, levels = c(6, 9, 12))
df.main$Measurement <- factor(df.main$Measurement, levels = main.measures)

# function to get data summary for barplots
data_summary <- function(data, varname, groupnames){
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE),
      Qup = quantile(x[[col]], probs=0.75, na.rm=TRUE)[[1]],
      Qdown = quantile(x[[col]], probs=0.25, na.rm=TRUE)[[1]],
      median = median(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  return(data_sum)
}

# convert significance values into asterisks
makeAsterisk <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i]
}

## barplots ----
df.main.summary <- data_summary(df.main, varname = "Value", groupnames = c("Age","Treatment","Measurement"))

media.border.cols <- c('CDM4' = "#5fb5c9",
                       'APM' = "#c593c2")
media.box.cols <- c('CDM4' = "#1a7c93",
                    'APM' = "#912883")
genotype_cols <- c("H1" = "#E9D8A6", "PGP1" = "#EE9B00", "11a" = "darkred")

measurements <- main.measures
main.prism.plots <- list()
for (measure in measurements) {
  df.main.measure <- df.main %>%
    filter(Measurement == measure)
  df.main.measure <- droplevels(df.main.measure)
  
  p <- ggplot(droplevels(df.main.summary[df.main.summary$Measurement==measure, ]),
              aes(x=Age, y=mean)) +
    geom_bar(stat="identity", width=0.6, aes(fill = Treatment, color = Treatment)) +
    scale_fill_manual(name = "Treatment", values = scales::alpha(media.border.cols, 0.5)) +
    scale_color_manual(name = "Treatment", values=scales::alpha(media.box.cols, 0.8)) +
    facet_grid(~Treatment, space = "free") +
    new_scale_color() +
    theme(panel.background = element_blank(),
          axis.line = element_line(color="black"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13)) +
    xlab("Age (months)") +
    ylab(measure) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2) +
    geom_point(data=df.main.measure, aes(x=Age, y=Value, color=Genotype),
               position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.4, seed = 42),
               size = 2.5, inherit.aes = FALSE) +
    # buffer some blank space above the maximum value but also maintain y-axis breaks
    coord_cartesian(ylim = c(0, max(df.main.measure$Value) * 1.05)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
    scale_color_manual(values=genotype_cols)
  
  main.prism.plots[[measure]] <- p
  rm(df.main.measure)
}

for (measure in measurements) {
  ggsave(paste0("Main_", measure, "_", Sys.Date(),".svg"), main.prism.plots[[measure]], path = outDir, width = 7, height = 4.4)
}

## ANOVA (Type III) stats ----
contrasts(df.main$Treatment) <- contr.sum
contrasts(df.main$Age) <- contr.sum

df.main.subs <- list()
fileConn <- file(paste0(outDir, "Main_ANOVA_Type3.txt"), open = "wt")
for (measure in measurements) {
  df.main.subs[[measure]] <- df.main %>%
    filter(Measurement == measure)
  
  # check if there is an age group that has missing measurements in any treatment conditions
  tbl <- table(df.main.subs[[measure]]$Treatment, df.main.subs[[measure]]$Age)
  age.to.be.dropped <- colnames(tbl)[apply(tbl == 0, 2, any)]
  df.main.subs[[measure]] <- df.main.subs[[measure]][df.main.subs[[measure]]$Age %nin% age.to.be.dropped, ]
  df.main.subs[[measure]] <- droplevels(df.main.subs[[measure]])
  
  print(paste0("Computing ANOVA for: ", measure, "..."))
  
  # if there is only one level of Genotype, remove Genotype from fixed effects
  if (length(levels(df.main.subs[[measure]]$Genotype)) == 1) {
    res <- Anova(lm(Value ~ Treatment*Age,
                    contrasts=list(Treatment='contr.sum', Age ='contr.sum'),
                    data = df.main.subs[[measure]]),
                 type='III')
  } else {
    res <- Anova(lm(Value ~ Genotype + Treatment*Age,
                    contrasts=list(Treatment='contr.sum', Age ='contr.sum'),
                    data = df.main.subs[[measure]]),
                 type='III')
  }
  
  writeLines(paste("\n\n### ", measure, ": ", "###"), fileConn)
  sink(fileConn, append = TRUE)
  print(res)
  
  sink()
  writeLines("\n--------------------------------------------------------\n", fileConn)
}
close(fileConn)

## post-hoc Tukey stats ----
# perform pairwise comparisons across all Treatment x Age combinations, and then pick same-age comparisons.
Tukey.all <- list()
Tukey.df.main.all <- c()
fileConn2 <- file(paste0(outDir, "Main_Pairwise_Tukey.txt"), open = "wt")
for (measure in measurements) {
  print(paste0("Computing post-hoc Tukey tests for: ", measure, "..."))
  # if there is only one level of Genotype, remove Genotype from fixed effects
  if (length(levels(df.main.subs[[measure]]$Genotype)) == 1) {
    aov.model <- aov(Value ~ Treatment*Age, data = df.main.subs[[measure]])
  } else {
    aov.model <- aov(Value ~ Genotype + Treatment*Age, data = df.main.subs[[measure]])
  }
  
  Tukey <- TukeyHSD(aov.model, conf.level=.95)
  Tukey.all[[measure]] <- Tukey
  
  writeLines(paste("\n\n### ", measure, ": ", "###"), fileConn2)
  sink(fileConn2, append = TRUE)
  print(Tukey)
  
  sink()
  writeLines("\n--------------------------------------------------------\n", fileConn2)
  
  Tukey.df.main <- data.frame(
    "Measurement" = measure, 
    "Comparison" = names(Tukey$`Treatment:Age`[rownames(Tukey$`Treatment:Age`) %in% c("APM:12-CDM4:12","APM:9-CDM4:9","APM:6-CDM4:6"),4]),
    "p.adj" = as.numeric(Tukey$`Treatment:Age`[rownames(Tukey$`Treatment:Age`) %in% c("APM:12-CDM4:12","APM:9-CDM4:9","APM:6-CDM4:6"),4]))
  Tukey.df.main$p.adj.scientific.notation <- scientific(Tukey.df.main$p.adj, digits = 2)
  Tukey.df.main$p.adj.signif <- makeAsterisk(Tukey.df.main$p.adj)
  Tukey.df.main.all <- rbind(Tukey.df.main.all, Tukey.df.main)
}
close(fileConn2)
saveRDS(Tukey.df.main.all, paste0(outDir, "Main_Pairwise_Tukey_filtered.RDS"))
write.xlsx(Tukey.df.main.all, paste0(outDir, "Main_Pairwise_Tukey_filtered.xlsx"), rownames = F)

# Supplementary ----
df.suppl <- rbind(df.H1[df.H1$Measurement %in% suppl.measures,], 
                  df.PGP1[df.PGP1$Measurement %in% suppl.measures,], 
                  df.11a[df.11a$Measurement %in% suppl.measures,])

df.suppl$Treatment <- factor(df.suppl$Treatment, levels = c("CDM4", "APM"))
df.suppl$Genotype <- factor(df.suppl$Genotype, levels = c("H1", "PGP1", "11a"))
df.suppl$Age <- factor(df.suppl$Age, levels = c(6, 9, 12))
df.suppl$Measurement <- factor(df.suppl$Measurement, levels = suppl.measures)

## barplots ----
df.suppl.summary <- data_summary(df.suppl, varname = "Value", groupnames = c("Age","Treatment","Measurement"))

measurements <- suppl.measures
suppl.prism.plots <- list()
for (measure in measurements) {
  df.suppl.measure <- df.suppl %>%
    filter(Measurement == measure)
  df.suppl.measure <- droplevels(df.suppl.measure)
  
  p <- ggplot(df.suppl.summary[df.suppl.summary$Measurement==measure,], aes(x=Age, y=mean)) +
    geom_bar(stat="identity", width=0.6, aes(fill = Treatment, color = Treatment)) +
    scale_fill_manual(name = "Treatment", values = scales::alpha(media.border.cols, 0.5)) +
    scale_color_manual(name = "Treatment", values = scales::alpha(media.box.cols, 0.8)) +
    facet_grid(~Treatment, space = "free") +
    new_scale_color() +
    theme(panel.background = element_blank(),
          axis.line = element_line(color="black"),
          strip.background = element_blank(),
          strip.text = element_blank(),
          axis.title = element_text(size = 15),
          axis.text = element_text(size = 13)) +
    xlab("Age (months)") +
    ylab(measure) +
    geom_errorbar(aes(ymin=mean, ymax=mean+sd), width=.2) +
    geom_point(data=df.suppl.measure, aes(x=Age, y=Value, color=Genotype),
               position=position_jitterdodge(jitter.width = 0.25, dodge.width = 0.4, seed = 42),
               size = 2.5, inherit.aes = FALSE) +
    # buffer some blank space above the maximum value but also maintain y-axis breaks
    coord_cartesian(ylim = c(0, max(df.suppl.measure$Value) * 1.05)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0))) + 
    scale_color_manual(values=genotype_cols)
  
  suppl.prism.plots[[measure]] <- p
  rm(df.suppl.measure)
}

for (measure in measurements) {
  ggsave(paste0("Suppl_", measure, "_", Sys.Date(),".svg"), suppl.prism.plots[[measure]], path = outDir, width = 7, height = 4.4)
}

## ANOVA (Type III) stats ----
contrasts(df.suppl$Treatment) <- contr.sum
contrasts(df.suppl$Age) <- contr.sum

df.suppl.subs <- list()
fileConn <- file(paste0(outDir, "Suppl_ANOVA_Type3.txt"), open = "wt")
for (measure in measurements) {
  df.suppl.subs[[measure]] <- df.suppl %>%
    filter(Measurement == measure)
  
  # check if there is an age group that has missing measurements in any treatment conditions
  tbl <- table(df.suppl.subs[[measure]]$Treatment, df.suppl.subs[[measure]]$Age)
  age.to.be.dropped <- colnames(tbl)[apply(tbl == 0, 2, any)]
  df.suppl.subs[[measure]] <- df.suppl.subs[[measure]][df.suppl.subs[[measure]]$Age %nin% age.to.be.dropped, ]
  df.suppl.subs[[measure]] <- droplevels(df.suppl.subs[[measure]])
  
  print(paste0("Computing ANOVA for: ", measure, "..."))
  
  # if there is only one level of Genotype, remove Genotype from fixed effects
  if (length(levels(df.suppl.subs[[measure]]$Genotype)) == 1) {
    res <- Anova(lm(Value ~ Treatment*Age,
                    contrasts=list(Treatment='contr.sum', Age ='contr.sum'),
                    data = df.suppl.subs[[measure]]),
                 type='III')
  } else {
    res <- Anova(lm(Value ~ Genotype + Treatment*Age,
                    contrasts=list(Treatment='contr.sum', Age ='contr.sum'),
                    data = df.suppl.subs[[measure]]),
                 type='III')
  }
  
  writeLines(paste("\n\n### ", measure, ": ", "###"), fileConn)
  sink(fileConn, append = TRUE)
  print(res)
  
  sink()
  writeLines("\n--------------------------------------------------------\n", fileConn)
}
close(fileConn)

## post-hoc Tukey stats  ----
# perform pairwise comparisons across all Treatment x Age combinations, and then pick same-age comparisons.
Tukey.all <- list()
Tukey.df.suppl.all <- c()
fileConn2 <- file(paste0(outDir, "Suppl_Pairwise_Tukey.txt"), open = "wt")
for (measure in measurements) {
  print(paste0("Computing post-hoc Tukey tests for: ", measure, "..."))
  # if there is only one level of Genotype, remove Genotype from fixed effects
  if (length(levels(df.suppl.subs[[measure]]$Genotype)) == 1) {
    aov.model <- aov(Value ~ Treatment*Age, data = df.suppl.subs[[measure]])
  } else {
    aov.model <- aov(Value ~ Genotype + Treatment*Age, data = df.suppl.subs[[measure]])
  }
  
  Tukey <- TukeyHSD(aov.model, conf.level=.95)
  Tukey.all[[measure]] <- Tukey
  
  writeLines(paste("\n\n### ", measure, ": ", "###"), fileConn2)
  sink(fileConn2, append = TRUE)
  print(Tukey)
  
  sink()
  writeLines("\n--------------------------------------------------------\n", fileConn2)
  
  Tukey.df.suppl <- data.frame(
    "Measurement" = measure, 
    "Comparison" = names(Tukey$`Treatment:Age`[rownames(Tukey$`Treatment:Age`) %in% c("APM:12-CDM4:12","APM:9-CDM4:9","APM:6-CDM4:6"),4]),
    "p.adj" = as.numeric(Tukey$`Treatment:Age`[rownames(Tukey$`Treatment:Age`) %in% c("APM:12-CDM4:12","APM:9-CDM4:9","APM:6-CDM4:6"),4]))
  Tukey.df.suppl$p.adj.scientific.notation <- scientific(Tukey.df.suppl$p.adj, digits = 2)
  Tukey.df.suppl$p.adj.signif <- makeAsterisk(Tukey.df.suppl$p.adj)
  Tukey.df.suppl.all <- rbind(Tukey.df.suppl.all, Tukey.df.suppl)
}
close(fileConn2)
saveRDS(Tukey.df.suppl.all, paste0(outDir, "Suppl_Pairwise_Tukey_filtered.RDS"))
write.xlsx(Tukey.df.suppl.all, paste0(outDir, "Suppl_Pairwise_Tukey_filtered.xlsx"), rownames = F)








