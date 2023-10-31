"""
NI_A34780C_Ext_Fig_5_nasal_TRB.py

generates:
1. 2022_06_23_all_AIM_phenotypes_in_nasal_pubilc.tsv
2. 2022_06_23_all_AIM_phenotypes_in_nasal.tsv
"""

path_to_nasal = '/Users/kmayerbl/active/david_koelle/nasal/nasal_v2'
from hybrid.parse import parse_adaptive_v2
import os
import numpy as np
fs = [f for f in os.listdir(path_to_nasal) if f.endswith("tsv")]
print(fs)

from hybrid.parse import parse_adaptive_v2
import pandas as pd

d = dict()
for f in fs:
 print(f)
df_adaptive = pd.read_csv(os.path.join(path_to_nasal, f), sep = "\t")
df_parsed   = parse_adaptive_v2(df_adaptive)
df_parsed['filename'] = f
df_parsed = df_parsed.sort_values('frequency', ascending = False)
df_parsed['ptid'] = f.split("_")[0]
df_parsed['rank'] = [x for x in range(df_parsed.shape[0])]
d[f] = df_parsed

dfall = pd.concat(d.values())

# Now load S-reactive
sclones = pd.read_csv('plot_data/2022_06_23_all_AIM_phenotypes.tsv', sep = "\t")
sclones =  sclones[sclones['10x_enriched'] == True]
sclones['ptid'] = sclones['source_file'].apply(lambda x: x.split("_")[0])
sclones['v_b_gene'] = sclones['v_b_gene_x'].apply(lambda x : f'{x}*01')
sclones['j_b_gene'] = sclones['j_b_gene_x'].apply(lambda x : f'{x}*01')
sclones['cdr3_b_aa'] = sclones['cdr3_b_aa_x']

sclones['score'] = np.log(sclones['total_counts_cd8']) - np.log(sclones['total_counts_cd4']) 
sclones['total_cd_counts'] = sclones['total_counts_cd8'] + sclones['total_counts_cd4']
def cd4cd8_score2(score, total_cd_counts):
 if total_cd_counts < 10:
 r = "LOW_CD"
else:
 if score > 1:
 r = "CD8"
elif score < -1:
 r = "CD4"
else:
 r = "CD4CD8"
return r
sclones['cell_type'] = sclones.apply(lambda r: cd4cd8_score2(score = r['score'], total_cd_counts = r['total_cd_counts']), axis = 1)



dfall_S = dfall.merge(sclones[['10x_cdf_p_fdr_10x','v_b_gene', 'j_b_gene', 'cdr3_b_aa','ptid',
                              'total_counts_cd8','total_counts_cd4','cdr3_a_nt', 'v_a_gene', 'j_a_gene', 'cdr3_a_aa', 'cdr3_b_nt', 'cell_type']], 
                     how = "left", on = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa','ptid'])
dfall_S.to_csv('plot_data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv', sep = "\t")



dfall_S_public = dfall.merge(sclones[['10x_cdf_p_fdr_10x','v_b_gene', 'j_b_gene', 'cdr3_b_aa','ptid',
                                     'total_counts_cd8','total_counts_cd4','cdr3_a_nt', 'v_a_gene', 'j_a_gene', 'cdr3_a_aa', 'cdr3_b_nt', 'cell_type','ptid']], 
                            how = "left", on = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'])
dfall_S_public.to_csv('plot_data/2022_06_23_all_AIM_phenotypes_in_nasal_pubilc.tsv', sep = "\t")
