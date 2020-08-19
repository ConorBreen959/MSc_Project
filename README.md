#### MSc. Research Project

##### "An analysis of the overlap between DNA methylation and gene expression patterns during maternal immune activation"



This directory contains the scripts used to process data for my masters research project which involved a comparison of gene expression and DNA metyhylation patterns in a mouse model of maternal immune activation. The final manuscript can be found [here](https://github.com/ConorBreen959/MSc_Project/blob/master/final_manuscript.pdf).



These scripts were run mainly through R to perform the three stages of analysis

1. [`Microarray.R`](https://github.com/ConorBreen959/MSc_Project/blob/master/Microarray.R) - To perform differential gene expression analysis on the microarray data used.
2. [`functional_enrichment.py`](https://github.com/ConorBreen959/MSc_Project/blob/master/functional_enrichment.py) - To find enrichment of gene ontology terms for each list of genes.
3. [`EWCE.R`](https://github.com/ConorBreen959/MSc_Project/blob/master/EWCE.R) - To perform expression-weighted cell type enrichment to identify neuronal cell types for which each gene list is more highly expressed.

-----

`Supplementary_Tables.ods` contains the differentially methylated and expressed gene lists used, their LogFC and adjusted P values, as well as results from each stage of the analysis.
