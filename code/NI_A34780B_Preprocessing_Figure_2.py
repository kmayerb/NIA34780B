# Preprocessing_Figure_2_data.py

# Preproccessing Script -- Figure 2 inputs 
# This file creates inputs for Figure_2.R

# This script depends on code in the /hyrbid/ package 

# Key Outputs: 
# '2022_06_23_all_long_data.tsv', '2022_08_31_all_long_data.tsv'

# Xnote: These are highly similar outputs but, the 8-31-22 version 
# includes more information on each clones alpha chain receptor sequences.
# Xnote: This code was first as a ipython notebook written as 2022_06_02_CD4_CD8_longitudinal_plots.ipynb
# Xnote: This code was first as a a ipython noteboo written as 2022_08_31_CD4_CD8_longitudinal_plots.ipynb

##############################
##############################
# 2022_06_02/23 -- Version 
import os 
import pandas as pd
import numpy as np
import json
import os 
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe

r = 'phenotypes/'
frames = dict()
for f in os.listdir(r):
    if f.endswith('phenotype_unique.tsv'):
        print(f)
        df = pd.read_csv(os.path.join(r,f), sep = "\t")
        df['source_file'] = f
        df['ptid'] = f.split("_")[0]
        df['visit'] = f.split("_")[1]
        df['visit_code']= f.split("_")[3]
        frames[f] = df
all_E03_phenotypes = pd.concat([v for k,v in frames.items() if k.find("E03")!=-1])
all_E03_phenotypes['score'] = np.log(all_E03_phenotypes['total_counts_cd8']) - np.log(all_E03_phenotypes['total_counts_cd4']) 
all_E03_phenotypes['total_cd_counts'] = all_E03_phenotypes['total_counts_cd8'] + all_E03_phenotypes['total_counts_cd4']

# New cell type analysis, should include double positive
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
all_E03_phenotypes['cell_type'] = all_E03_phenotypes.apply(lambda r: cd4cd8_score2(score = r['score'], total_cd_counts = r['total_cd_counts']), axis = 1)
plots = []
summary_data = []
long_data_frames= []
count = 0
get_fold_changes = False
for ptid in all_E03_phenotypes['ptid'].value_counts().index:
    print(ptid)
    count = count+1
    phenotype_df = all_E03_phenotypes.query(f'ptid == "{ptid}"')
    os.system(f"head json/{ptid}_samples.json")

    json_file = f'json/{ptid}_samples.json'

    with open(json_file, 'r') as json_obj:
        samples = json.load(json_obj)

    print(f"CHECKING {ptid} SAMPLES EXIST")
    for sample_name, fp in samples.items():
        assert os.path.isfile(fp), f"{sample_name} FILE {fp} DOES NOT EXIST"
        print(f"SUCCESS  {sample_name} FILE {fp} EXISTS") 

   
    dft = make_wide_longitudinal_dataframe(
        sample_dictionary = samples, 
        dest = ".",
        ptid = None,
        write = False)
    
    if get_fold_changes:
        dftextra = compute_fold_changes_on_longitudinal_dataframe(
            dft = dft.copy(), 
            dest = ".",
            ptid = None,
            write = False)
    else:
        dftextra = dft.copy()
    # Combine Phenoytpe Information
    dftextra['v_b']  = dftextra['v_b_gene'].apply(lambda x: x.split("*")[0])
    dftextra['j_b']  = dftextra['j_b_gene'].apply(lambda x: x.split("*")[0])
    df_celltype = dftextra.merge(phenotype_df, how = "left", 
        left_on  = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b', 'j_b'],
        right_on = ['cdr3_b_nucseq', 'cdr3_b_aa_x', 'v_b_gene_x', 'j_b_gene_x'])

    # This is the critical line that needs to be modified to retain
    # JUNE 20th Modified df_celltype = df[df['cell_type'].notna()]
    # remove 10X cells which may not be enriched in spike specific pool
    df_celltype = df_celltype[(df_celltype['10x_cdf_p_fdr_10x'] < 0.05) & (df_celltype['10x_cdf_p_fdr_10x'].notna())]
    # confirm that we are using 10x enriched clones only
    print(df_celltype['10x_enriched'].all())
    assert df_celltype['10x_enriched'].all()
    dfp = df_celltype.copy()
    dfp['plot_clone_id'] = [i for i in range(dfp.shape[0])]
    # Get the value vars based on sample names
    value_vars =  [f"{x}_pfreq" for x in samples.keys()]
       
    print(value_vars)
    # Xnote: in this version plot columns only included v_b, 
    plotcols = ['plot_clone_id','cell_type','cd8vcd4_score','cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family','j_family',] + value_vars 
    print(plotcols)
    id_vars = [v for v in plotcols if v not in value_vars]
    df_long = pd.melt(dfp[plotcols], id_vars = id_vars, value_vars=value_vars)
    print(df_long.head())
    labels = {x:x.split("_")[3] for x in value_vars}
    print(labels)
    df_long['x_pos'] = df_long['variable'].apply(lambda x: labels.get(x))
    df_long['ptid'] = ptid
    long_data_frames.append(df_long)
    df_long_sum = df_long.groupby(['variable','cell_type'])['value'].mean().reset_index()
    df_long_max = df_long.groupby(['variable','cell_type'])['value'].max().reset_index()
    df_long_min = df_long.groupby(['variable','cell_type'])['value'].min().reset_index()
    df_long_sum['x_pos'] = df_long_sum['variable'].apply(lambda x: labels.get(x))
    df_long_max['x_pos'] = df_long_max['variable'].apply(lambda x: labels.get(x))
    df_long_min['x_pos'] = df_long_min['variable'].apply(lambda x: labels.get(x))
    df_long_sum['summary'] = "mean"
    df_long_min['summary'] = "min"
    df_long_max['summary'] = "max"
    
    summary_df = pd.concat([df_long_sum, df_long_max, df_long_min])
    summary_data.append(summary_df)
#     # We previously used plotnine, but now export and plot usign ggplot2
#     from plotnine import *
#     p = (ggplot(df_long, aes(x='x_pos', y='value')) 
#          + scale_y_log10()
#          + geom_point(size = .1) 
#          + geom_line(aes(group = 'plot_clone_id', color= 'cd8vcd4_score'), size = .1) 
#          + theme_bw() 
#          + facet_wrap('~cell_type', ncol = 1) 
#          + ggtitle(f"{ptid}: Trajectories of TCRB matching the E03 S-specific AIM single cells") 
#          + scale_color_gradient2(low = "green",mid = "white", high = "blue")
#          + annotation_logticks(sides = 'l')) 
#     p.save(f'figures/{ptid}_trajectory_cd4_cd8_clonal_frequency.pdf', height=5, width=8)
#     plots.append(p)
all_long_data_used = pd.concat(long_data_frames)
all_long_data_used.to_csv('2022_06_23_all_long_data.tsv', sep = "\t", index = False)
all_summary_data = pd.concat(summary_data)
all_summary_data.to_csv('2022_06_23_all_summary_data.tsv', sep = "\t", index = False)
############################################
# df_contract, df_expand               #####
# Focus on 12 PTIDs with complete data #####
############################################
# Xnote: was done originally in an ipython notebook : 2022_06_21_Expansion_Contraction.ipynb
import os 
import pandas as pd
r = 'phenotypes/'
frames = dict()
print("Loading Phenotypes_uniuqe.tsv")
for f in os.listdir(r):
    if f.endswith('phenotype_unique.tsv'):
        df = pd.read_csv(os.path.join(r,f), sep = "\t")
        df['source_file'] = f
        df['ptid'] = f.split("_")[0]
        df['visit'] = f.split("_")[1]
        df['visit_code']= f.split("_")[3]
        frames[f] = df
all_E03_phenotypes = pd.concat([v for k,v in frames.items() if k.find("E03")!=-1])
all_E03_phenotypes['10x_enriched'] == True 
import numpy as np
all_E03_phenotypes['score'] = np.log(all_E03_phenotypes['total_counts_cd8']) - np.log(all_E03_phenotypes['total_counts_cd4']) 
all_E03_phenotypes['total_cd_counts'] = all_E03_phenotypes['total_counts_cd8'] + all_E03_phenotypes['total_counts_cd4']
# New cell type analysis, should include double positive
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
all_E03_phenotypes['cell_type'] = all_E03_phenotypes.apply(lambda r: cd4cd8_score2(score = r['score'], total_cd_counts = r['total_cd_counts']), axis = 1)
# Looking at teh 12 PTID
import os 
import pandas as pd
import json
import os 
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe
results_expand = dict()
results_contract = dict()
for ptid in ['15525', '15527', '15531', '15577', '15581', '15669', '15673', '15684', '15754', '15761', '15836','15782']:
    
    get_fold_changes = True

    json_file = f'json/{ptid}_samples.json'

    with open(json_file, 'r') as json_obj:
        samples = json.load(json_obj)

    print(f"CHECKING {ptid} SAMPLES EXIST")
    for sample_name, fp in samples.items():
        assert os.path.isfile(fp), f"{sample_name} FILE {fp} DOES NOT EXIST"
        print(f"SUCCESS  {sample_name} FILE {fp} EXISTS") 


    dft = make_wide_longitudinal_dataframe(
        sample_dictionary = samples, 
        dest = ".",
        ptid = None,
        write = False)

    if get_fold_changes:
        dftextra = compute_fold_changes_on_longitudinal_dataframe(
            dft = dft.copy(), 
            dest = ".",
            ptid = None,
            write = False)
    else:
        dftextra = dft.copy()


    phenotype_df = all_E03_phenotypes.query(f'ptid == "{ptid}"')
    dftextra['v_b']  = dftextra['v_b_gene'].apply(lambda x: x.split("*")[0])
    dftextra['j_b']  = dftextra['j_b_gene'].apply(lambda x: x.split("*")[0])
    dftextra['longitudinal_id'] = [i for i in range(dftextra.shape[0])]
    df_celltype = dftextra.merge(phenotype_df, how = "left", 
        left_on  = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b', 'j_b'],
        right_on = ['cdr3_b_nucseq', 'cdr3_b_aa_x', 'v_b_gene_x', 'j_b_gene_x'])
    df_celltype['ptid'] = ptid

    e03_column = [x for x in dftextra.columns if x.find("E03") != -1][0]
    e01_column = [x for x in dftextra.columns if x.find("E01") != -1][0]
    # index on expansion
    ixe = df_celltype[f"{e03_column}_vac_class"] == "sig_expand" 
    # index on fold change > 2
    ixe2 = df_celltype[f"{e03_column}_vac_fc"] > 4

    ixc = df_celltype[f"{e03_column}_vac_class"] == "sig_contract" 
    ixc2 = df_celltype[f"{e03_column}_vac_fc"] < 1/4
    
    results_expand[ptid] = df_celltype[ixe2]
    results_contract[ptid] = df_celltype[ixc2]
    # Originally we previewed data with plotnine in python
    # Production graphics done in ggplot2 instead
    # do_plot = False    
    # if do_plot: 
    #     # index on fold chnage > 5
    #     ixe5 = df_celltype[f"{e03_column}_vac_fc"] > 5
    # 
    #     d2 = df_celltype[ixe&ixe2]
    #     d5 = df_celltype[ixe&ixe5]
    #     d2[f"{e01_column}_pfreq"] = d2[f"{e01_column}_pfreq"].apply(lambda x: x if x > 0 else 0.000001)
    #     d2[f"{e03_column}_pfreq"] = d2[f"{e03_column}_pfreq"].apply(lambda x: x if x > 0 else 0.000001)
    #     d5[f"{e01_column}_pfreq"] = d5[f"{e01_column}_pfreq"].apply(lambda x: x if x > 0 else 0.000001)
    #     d5[f"{e03_column}_pfreq"] = d5[f"{e03_column}_pfreq"].apply(lambda x: x if x > 0 else 0.000001)
    #     fig,plt = (ggplot(d5, aes(x= f"{e01_column}_pfreq", 
    #                    y = f"{e03_column}_pfreq",
    #                    color = 'cell_type' ))
    #          + geom_point() 
    #          + geom_point(data = d2, size = .1)
    #          + geom_abline() 
    #          + scale_y_log10() 
    #          + scale_x_log10()
    #          + coord_cartesian(xlim = (-6, -2), ylim = (-6, -2)) 
    #         ).draw(show=False, return_ggplot=True)
    #     fig.savefig(f"mock/{ptid}.pdf")

def count_expanded_clones(df, ptid , fc_threshold = 4):
    
    generic_cols = {'ptid': 'ptid',
     'cdr3_b_nucseq': 'cdr3_b_nucseq',
     'cdr3_b_aa': 'cdr3_b_aa',
     'v_b_gene': 'v_b_gene',
     'j_b_gene': 'j_b_gene',
     'v_family': 'v_family',
     'j_family': 'j_family',
     'v_b': 'v_b',
     'j_b': 'j_b',
     'longitudinal_id': 'longitudinal_id',
     'cdr3_a_nt': 'cdr3_a_nt',
     'v_a_gene': 'v_a_gene',
     'j_a_gene': 'j_a_gene',
     'cdr3_a_aa': 'cdr3_a_aa',
     'cdr3_b_nt': 'cdr3_b_nt',
     'v_b_gene_x': 'v_b_gene_x',
     'j_b_gene_x': 'j_b_gene_x',
     'cdr3_b_aa_x': 'cdr3_b_aa_x',
     'total_counts_cd8': 'total_counts_cd8',
     'total_counts_cd4': 'total_counts_cd4',
     'cell_type': 'cell_type',
     '10x_pfreq': '10x_pfreq',
     '10x_enriched': '10x_enriched',
     'clones_file_adpt': 'clones_file_adpt'}
    
    e03_column = [x for x in df.columns if x.find("E03") != -1][0]
    e01_column = [x for x in df.columns if x.find("E01") != -1][0]
    
    e00_columns = [x for x in df.columns if x.find("E00") != -1]
    e005_columns = [x for x in df.columns if x.find("E00.5") != -1]
    e01_columns = [x for x in df.columns if x.find("E01") != -1]
    e02_columns = [x for x in df.columns if x.find("E03") != -1]
    e03_columns = [x for x in df.columns if x.find("E03") != -1]
    e04_columns = [x for x in df.columns if x.find("E04") != -1]
    e05_columns = [x for x in df.columns if x.find("E05") != -1]
    
    e00_columns_rep  =  {x : x.split("_TCRB_")[1] for x in e00_columns}
    e005_columns_rep =  {x : x.split("_TCRB_")[1] for x in e005_columns}
    e01_columns_rep = {x : x.split("_TCRB_")[1] for x in e01_columns}
    e02_columns_rep = {x : x.split("_TCRB_")[1] for x in e01_columns}
    e03_columns_rep = {x : x.split("_TCRB_")[1] for x in e03_columns}
    e04_columns_rep = {x : x.split("_TCRB_")[1] for x in e04_columns}
    e05_columns_rep = {x : x.split("_TCRB_")[1] for x in e05_columns}
    
    rep_cols = {**generic_cols, **e00_columns_rep,  **e005_columns_rep,  **e01_columns_rep, **e02_columns_rep, **e03_columns_rep,**e04_columns_rep,**e05_columns_rep}
    df_generic = df[rep_cols.keys()].rename(columns = rep_cols)

    return(df_generic)

    
df_expand   = pd.concat([count_expanded_clones(df = v, ptid = k) for k,v in results_expand.items()])
df_contract = pd.concat([count_expanded_clones(df = v, ptid = k) for k,v in results_contract.items()])
df_expand.to_csv('2022_06_23_df_expand.tsv', sep = "\t", index = False)
df_contract.to_csv('2022_06_23_df_contract.tsv', sep = "\t", index = False)

##############################
##############################
##############################
##############################
# 2022_08_31 -- Slightly modified Version for Ext Figure 7
# We load all 17 phenotype tables, 'phenotype_unique' suffix
# refers to phenotypes that could be matched directly to bulk clone for
# confirmation and trajectory analysis

import os 
import pandas as pd
import numpy as np
import os 
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe
# Previously the phenotypes of AIM cells were generated. 
# Here we compile all TCRab across participants
# <all_E03_phenotypes>
r = 'phenotypes/'
frames = dict()
for f in os.listdir(r):
    if f.endswith('phenotype_unique.tsv'):
        print(f)
        df = pd.read_csv(os.path.join(r,f), sep = "\t")
        df['source_file'] = f
        df['ptid'] = f.split("_")[0]
        df['visit'] = f.split("_")[1]
        df['visit_code']= f.split("_")[3]
        frames[f] = df
        
all_E03_phenotypes = pd.concat([v for k,v in frames.items() if k.find("E03")!=-1])

# XNote: Cell type assignment can optionally included CD4/8 double positives
all_E03_phenotypes['score'] = np.log(all_E03_phenotypes['total_counts_cd8']) - np.log(all_E03_phenotypes['total_counts_cd4']) 
all_E03_phenotypes['total_cd_counts'] = all_E03_phenotypes['total_counts_cd8'] + all_E03_phenotypes['total_counts_cd4']
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
# Apply scoring
all_E03_phenotypes['cell_type'] = all_E03_phenotypes.apply(lambda r: cd4cd8_score2(score = r['score'], total_cd_counts = r['total_cd_counts']), axis = 1)
assert all_E03_phenotypes['ptid'].value_counts().shape[0] == 17

# Originally this script made plots in plotnine, but we have instead used this 
# data to make plots with ggplot2
plots = []
summary_data = []
long_data_frames= []
count = 0
get_fold_changes = False
for ptid in all_E03_phenotypes['ptid'].value_counts().index:
    print(ptid)
    count = count+1
    phenotype_df = all_E03_phenotypes.query(f'ptid == "{ptid}"')
    os.system(f"head json/{ptid}_samples.json")

    json_file = f'json/{ptid}_samples.json'

    with open(json_file, 'r') as json_obj:
        samples = json.load(json_obj)

    print(f"CHECKING {ptid} SAMPLES EXIST")
    for sample_name, fp in samples.items():
        assert os.path.isfile(fp), f"{sample_name} FILE {fp} DOES NOT EXIST"
        print(f"SUCCESS  {sample_name} FILE {fp} EXISTS") 

   
    dft = make_wide_longitudinal_dataframe(
        sample_dictionary = samples, 
        dest = ".",
        ptid = None,
        write = False)
    
    if get_fold_changes:
        dftextra = compute_fold_changes_on_longitudinal_dataframe(
            dft = dft.copy(), 
            dest = ".",
            ptid = None,
            write = False)
    else:
        dftextra = dft.copy()
    # Combine Phenoytpe Information
    dftextra['v_b']  = dftextra['v_b_gene'].apply(lambda x: x.split("*")[0])
    dftextra['j_b']  = dftextra['j_b_gene'].apply(lambda x: x.split("*")[0])
    df_celltype = dftextra.merge(phenotype_df, how = "left", 
        left_on  = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b', 'j_b'],
        right_on = ['cdr3_b_nucseq', 'cdr3_b_aa_x', 'v_b_gene_x', 'j_b_gene_x'])

    # remove 10X cells which may not be enriched in spike specific pool, applying most stringent
    df_celltype = df_celltype[(df_celltype['10x_cdf_p_fdr_10x'] < 0.05) & (df_celltype['10x_cdf_p_fdr_10x'].notna())]
    # confirm that we are using 10x enriched clones only
    print(df_celltype['10x_enriched'].all())
    assert df_celltype['10x_enriched'].all()
    dfp = df_celltype.copy()
    dfp['plot_clone_id'] = [i for i in range(dfp.shape[0])]
    # Get the value vars based on sample names
    value_vars =  [f"{x}_pfreq" for x in samples.keys()]
       
    print(value_vars)
    # Xnote: in this version plot columns only included v_b in this version we also include relevant alpha columns
    plotcols = ['plot_clone_id','cell_type','cd8vcd4_score','cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family','j_family','cdr3_a_aa', 'v_a_gene', 'j_a_gene'] + value_vars 
    print(plotcols)
    id_vars = [v for v in plotcols if v not in value_vars]
    df_long = pd.melt(dfp[plotcols], id_vars = id_vars, value_vars=value_vars)
    print(df_long.head())
    labels = {x:x.split("_")[3] for x in value_vars}
    print(labels)
    df_long['x_pos'] = df_long['variable'].apply(lambda x: labels.get(x))
    df_long['ptid'] = ptid
    long_data_frames.append(df_long)
    df_long_sum = df_long.groupby(['variable','cell_type'])['value'].mean().reset_index()
    df_long_max = df_long.groupby(['variable','cell_type'])['value'].max().reset_index()
    df_long_min = df_long.groupby(['variable','cell_type'])['value'].min().reset_index()
    df_long_sum['x_pos'] = df_long_sum['variable'].apply(lambda x: labels.get(x))
    df_long_max['x_pos'] = df_long_max['variable'].apply(lambda x: labels.get(x))
    df_long_min['x_pos'] = df_long_min['variable'].apply(lambda x: labels.get(x))
    df_long_sum['summary'] = "mean"
    df_long_min['summary'] = "min"
    df_long_max['summary'] = "max"
    
    summary_df = pd.concat([df_long_sum, df_long_max, df_long_min])
    summary_data.append(summary_df)
#   We skip plotting in Python as we instead us3 ggplot2 in R:
#   Code is left commented out a reference 
#     from plotnine import *
#     p = (ggplot(df_long, aes(x='x_pos', y='value')) 
#          + scale_y_log10()
#          + geom_point(size = .1) 
#          + geom_line(aes(group = 'plot_clone_id', color= 'cd8vcd4_score'), size = .1) 
#          + theme_bw() 
#          + facet_wrap('~cell_type', ncol = 1) 
#          + ggtitle(f"{ptid}: Trajectories of TCRB matching the E03 S-specific AIM single cells") 
#          + scale_color_gradient2(low = "green",mid = "white", high = "blue")
#          + annotation_logticks(sides = 'l')) 
#    
#     p.save(f'figures/{ptid}_trajectory_cd4_cd8_clonal_frequency.pdf', height=5, width=8)
#     plots.append(p)
# compile and output as long data for analysis
all_long_data_used = pd.concat(long_data_frames)
all_long_data_used['j_a_gene'] = all_long_data_used['j_a_gene'].apply(lambda x : f"{x}*01")
all_long_data_used['v_a_gene'] = all_long_data_used['v_a_gene'].apply(lambda x : f"{x}*01")
all_long_data_used.to_csv('2022_08_31_all_long_data.tsv', sep = "\t", index = False)
all_summary_data = pd.concat(summary_data)
all_summary_data.to_csv('2022_08_31_all_summary_data.tsv', sep = "\t", index = False)
