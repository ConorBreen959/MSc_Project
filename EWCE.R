#! /usr/bin/env R

library(devtools)
library(biomaRt)
library(ggplot2)
library(cowplot)
library(limma)
library(readxl)
library(tidyverse)

# Install EWCE package using devtools
install_github("nathanskene/ewce")
library(EWCE)

## Preamble

# Load single-cell raw expression and annotation datafiles
exp <- read.csv("data/df_agg.csv",header=TRUE, row.names=1)
exp_mat <- data.matrix(exp, rownames.force = TRUE)
annot <- read.csv("data/annotation_agg.txt",header=TRUE, sep = "\t")

# Load and process desired gene lists; commands unique to this analysis pipeline and dataset
DEG_total_list <- read.csv("data/DEG_list.csv")
methyl_genes <- read.csv("data/Methylation_genes.csv")
DEG_list <- DEG_total_list[DEG_total_list$adj.P.Val <= 0.05 & abs(DEG_total_list$logFC) >= 0.137, ];
DEG_downreg <- DEG_list %>% filter(logFC < 0) %>% pull(Gene)
DEG_upreg <- DEG_list %>% filter(logFC > 0) %>% pull(Gene)
DMG_hypo <- methyl_genes %>% filter(Methylation < 0) %>% pull(Gene)
DMG_hyper <- methyl_genes %>% filter(Methylation > 0) %>% pull(Gene)

# Download gene symbol synonym file from MGI if not present
if(!file.exists("data/MRK_List2.rpt")){
  download.file("http://www.informatics.jax.org/downloads/reports/MRK_List2.rpt", destfile="data/MRK_List2.rpt")
}

# Fix bad MGI gene symbols in the expression matrix 
expData = fix.bad.mgi.symbols(exp_mat,mrk_file_path="data/MRK_List2.rpt")

#Set up comparison levels
l1=annot$ClusterName
l2=annot$TaxRank2
l3=annot$TaxRank4
annotLevels = list(l1=l1,l2=l2,l3=l3)

# Generate and load celltype data based on expression and annotation files, to be used in the enrichment test
fNames_Zeisel = generate.celltype.data(exp=expData,annotLevels,"Zeisel")
load(fNames_Zeisel[1])

# Get background geneset using bioMart
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
all_hgnc <- getBM(attributes = ("hgnc_symbol"), mart = human)
human.bg <- all_hgnc$hgnc_symbol


# Function set up to automate performing EWCE test. Genelist, ctd file, bootstrapping reps, and annotation level are required arguments.

perform.EWCE <- function(genelist, ctd, reps, annotLevel){
  # Convert genelist MGI symbols to human HGNC gene symbols using bioMart
  genes = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genelist, mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  human.hits <- genes$HGNC.symbol
  
  # Enrichment test
  results = bootstrap.enrichment.test(sct_data=ctd,hits=human.hits,bg=human.bg, 
                                          geneSizeControl=TRUE,genelistSpecies="human",sctSpecies="mouse", reps=reps, annotLevel=annotLevel)
  
  # Format results and write to file in results folder
  results <- results$results
  results <- merge(results, annot, by.x = 'CellType', by.y = 'ClusterName')
  results <- results[order(results$p),]
  write.csv(results, file = paste("results/EWCE_", deparse(substitute(genelist)),".csv",sep=""), row.names = F)
}

# Set up vector with desired genelists to test
genelists <- c(DEG_downreg, DEG_upreg, DMG_hyper, DMG_hypo)

# Perform EWCE analysis on gene lists

for(i in genelists) {
  perform.EWCE(i, ctd, 100, 1)
}

