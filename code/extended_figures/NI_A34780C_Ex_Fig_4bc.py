# NA_A34780C_Ext_Fig_4bc.py
# Makes P525 and P581 AIM vs. bulk comparison files, without extra columns.
# This file is important for providing reasonable sized source data.

import os
import pandas as pd
d525 = pd.read_csv('/fh/fast/corey_l/esford3_kmayerbl_collab/tidy covid folder/15525.longitudinal_templates_only.tsv', sep ="\t")
d525cs = {
  'cdr3_b_nucseq':'cdr3_b_nucseq',
  'cdr3_b_aa':'cdr3_b_aa',
  'v_b_gene':'v_b_gene',
  'j_b_gene':'j_b_gene',
  '15525_3_TCRB_E01_pfreq':'E01' ,
  '15525_5_TCRB_E03_pfreq':'E03'}
df5 = d525[d525cs.keys()].rename(columns = d525cs)
df5['ptid'] = 525

d581= pd.read_csv('/fh/fast/corey_l/esford3_kmayerbl_collab/tidy covid folder/15581.longitudinal_templates_only.tsv', sep ="\t")
d581cs = {
  'cdr3_b_nucseq':'cdr3_b_nucseq',
  'cdr3_b_aa':'cdr3_b_aa',
  'v_b_gene':'v_b_gene',
  'j_b_gene':'j_b_gene',
  '15581_7_TCRB_E01_pfreq':'E01' ,
  '15581_9_TCRB_E03_pfreq':'E03'}
df1 = d581[d581cs.keys()].rename(columns = d581cs)
df1['ptid'] = 581

df = pd.concat([df5,df1])

df_sums = df[['E01','E03']].sum(axis = 1) > 0
df2 = df[df_sums]

# key file containing AIM information
key_cols = ['10x_pfreq','ptid','cdr3_b_aa','v_b_gene','j_b_gene','cdr3_b_nucseq']
key = pd.read_csv('/fh/fast/corey_l/esford3_kmayerbl_collab/software/NIA34780B/data/aggregate_phenotypes_unique_101223.csv')
key['v_b_gene'] = key['v_b_gene_x'].apply(lambda x: f"{x}*01")
key['j_b_gene'] = key['j_b_gene_x'].apply(lambda x: f"{x}*01")
key['cdr3_b_aa'] =key['cdr3_b_aa_x'] 
dfm = df2.merge(key[key_cols], how = "left", 
         left_on = ['ptid','cdr3_b_aa','v_b_gene','j_b_gene','cdr3_b_nucseq'],
         right_on = ['ptid','cdr3_b_aa','v_b_gene','j_b_gene','cdr3_b_nucseq'])


dfm.to_csv('/fh/fast/corey_l/esford3_kmayerbl_collab/software/NIA34780B/data/E01_E03_10x_match_525_581_cdr3vj.csv',index = False)
dfm[['ptid','E01','E03','10x_pfreq']].to_csv('/fh/fast/corey_l/esford3_kmayerbl_collab/software/NIA34780B/data/E01_E03_10x_match_525_581.csv',index = False)

dfaim = dfm[['ptid','E01','E03','10x_pfreq']][dfm['10x_pfreq'].notna()]
dfaim['count'] = 1
dfnm = dfm[['ptid','E01','E03','cdr3_b_aa']][dfm['10x_pfreq'].isna()].\
  groupby(['ptid','E01','E03']).count().\
  reset_index(drop = False).\
  rename(columns = {'cdr3_b_aa':'count'})

dfunique_plot = pd.concat([dfaim, dfnm]).sort_values(['ptid','count'], ascending = [True,False]).reset_index(drop = True)
dfunique_plot.to_csv('/fh/fast/corey_l/esford3_kmayerbl_collab/software/NIA34780B/data/E01_E03_10x_match_525_581_plot_unique.csv', index = False)







