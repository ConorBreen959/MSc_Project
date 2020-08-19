# -*- coding: utf-8 -*-
from gprofiler import GProfiler
gp = GProfiler(return_dataframe=True)
import numpy as np
import scipy as sp
import pandas as pd
import os
os.chdir("/home/conor/Documents/Git_Repositories/MSc_Project")

DEG_list = pd.read_csv('data/DEG_list.csv')
methyl_genes = pd.read_csv('data/Methylation_genes.csv')

upreg = DEG_list[(DEG_list['adj.P.Val'] <= 0.1) & (DEG_list['logFC'] >= 0.1)]
upreg = upreg['Gene'].astype(str).tolist()
upreg[:] = map(str.strip,upreg)

downreg = DEG_list[(DEG_list['adj.P.Val'] <= 0.1) & (DEG_list['logFC'] <= -0.1)]
downreg = downreg['Gene'].astype(str).tolist()
downreg[:] = map(str.strip,downreg)

hyper = methyl_genes[methyl_genes['Methylation'] > 0]
hyper = hyper['Gene'].tolist()

hypo = methyl_genes[methyl_genes['Methylation'] < 0]
hypo = hypo['Gene'].tolist()

genelists = {'downreg': downreg, 'upreg': upreg, 'hyper': hyper, 'hypo': hypo}

for i in genelists:
    print("Calculating", i, "enrichment...")
    enrichment = gp.profile(genelists[i], organism='mmusculus', significance_threshold_method='fdr', measure_underrepresentation=False, sources=['GO:BP'])
    enrichment.sort_values('p_value').to_csv(str(i) + "_enrichment.csv")
    locals()[str(i) + "_enrichment"] = enrichment
    print("Done")
    
intersection_1 = pd.merge(downreg_enrichment, hyper_enrichment, how = 'inner', on = ['native'])
intersection_2 = pd.merge(upreg_enrichment, hypo_enrichment, how = 'inner', on = ['native'])
