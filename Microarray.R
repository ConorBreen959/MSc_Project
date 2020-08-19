#! /usr/bin/env R

## Load required packages

library(oligo)                              # Bioconductor package for preprocessing of Affymetrix oligonucleotide arrays
library(pd.mogene.1.1.st.v1)                # Microarray chip-specific design annotation package
library(mogene11sttranscriptcluster.db)     # Affymetrix transcript annotation data
library(limma)                              # Limma package for differential expression analysis
library(gplots)
library(tidyverse)
library(janitor)
library(kableExtra)
library(viridis)

## Preprocessing

# Read in platform design file data and set up celfiles object
pd <- read.AnnotatedDataFrame("data/pData.txt")
celfiles <- paste("data/", rownames(pd), sep="")

# Read raw data from CEL files based on platform design information
rawData <- read.celfiles(celfiles, phenoData=pd)

# Normalise raw data
normData <- rma(rawData, target = "core")

# Sanity check on normalised data
pData(normData)
fData(normData)

# Retrieve transcript annotation data
featureData(normData) <- getNetAffx(normData, "transcript")

# Filter normalised data by "main" category data
geneData <- normData[fData(normData)[, "category"] %in% "main", ]

## Linear Modelling

# Set up design object with sample data 
design <- cbind(C=as.numeric(pd$studyGroup=="C"),V=as.numeric(pd$studyGroup=="V"))

# Check design object
design

# Calculate relative quality weights for each array; based on how well expression values in each array follow the linear model
aw <- arrayWeights(geneData, design)

# Array-weighted linear model
fit <- lmFit(geneData, design, weight=aw)

# Contrasts matrix for sample comparisons and moderated F-statistic
cm <- makeContrasts(V-C, levels=design)
fit2 <- contrasts.fit(fit, cm)
fit2 <- eBayes(fit2)

## Calculate differentially expressed genes

# FDR-adjusted top table of DEGs
top <- topTable(fit2, adjust.method = 'fdr', coef=1, n=nrow(fit2))

# Filter DEGs by desired significance and fold change thresholds
sig <- top[top$adj.P.Val <= 0.05 & abs(top$logFC) >= 0.137, ];

# Retrieve expression values for significant genes
e <- exprs(geneData)
sig_exprs <- e[rownames(e) %in% as.character(sig$transcriptclusterid), ]

## Visualisations

# Write plots to results folder
pdf("results/Results.pdf")

# Setup
lab <- paste(pd$sampleID, pd$studyGroup)
colours <- as.numeric(factor(pd$studyGroup))+1

# Raw data box plot of expression
boxplot(rawData, target="core", main="Raw Data", ylab="log2(exprs)", names=lab, las=2)

# Normalised data box plot of expression
boxplot(geneData, main="Norm Data", ylab="log2(exprs)", names=lab, las=2)

# Expression density histogram
hist(geneData, main="Normalised Expression Density")

# Heatmap of expression correlation between sample groups
heatmap.2(cor(exprs(geneData))^2, trace="none", scale="none", margins=c(9,9), labRow=lab, labCol=lab)

# MDS plot of distance between samples
plotMDS(geneData, labels=lab, col=colours)

# Array weightings bar plot
barplot(aw)

# Volcano plot of fold change vs p-value for all DEGs, and significant DEGs marked in purple and yellow
ggplot(top, aes(x =  logFC, y = -log10(adj.P.Val))) +
  geom_point(color = "#20A387FF", size = 0.05) +
  geom_point(data = sig, aes(x = logFC, y = -log10(adj.P.Val), colour = logFC > 0), size = 0.05) +
  scale_x_continuous(breaks = seq(-1.5, 0.5, 0.5), minor_breaks = seq(-1.75,1,0.25)) +
  scale_y_continuous(breaks = seq(0, 6, 1)) +
  scale_colour_manual(name = 'PC1 > 0', values = setNames(c("#FDE725FF", "#481567ff"),c(T, F))) +
  xlab("Fold Change") + 
  ylab("P Value (log10)") +
  theme_minimal(base_size = 20) +
  theme(legend.position = "none")

# Expression heatmap
heatmap.2(sig_exprs, trace="none", scale="row", col="redgreen", cexRow=0.2, cexCol=0.7)
dev.off()

## Export Data

# Tidyverse trickery to set up dataframes for readability and export 
top <- top %>% separate(geneassignment, c('a','b','c','d','e'), sep = "//") %>% select(transcriptclusterid, seqname, start, stop, a, b, c, logFC, adj.P.Val) %>% rename(`Ensembl ID` = a, `MGI ID` = b, `Gene Name` = c)
sig <- sig %>% separate(geneassignment, c('a','b','c','d','e'), sep = "//") %>% select(transcriptclusterid, seqname, start, stop, a, b, c, logFC, adj.P.Val) %>% rename(`Ensembl ID` = a, `MGI ID` = b, `Gene Name` = c)

# Write excel files to results folder
write.table(top, "results/DEG_list.csv", sep = ",", row.names = F)
write.csv(sig, "results/sig_gene_list.csv", row.names = F)
write.csv(sig_exprs, "results/Sig_exprs.csv", row.names=F)
