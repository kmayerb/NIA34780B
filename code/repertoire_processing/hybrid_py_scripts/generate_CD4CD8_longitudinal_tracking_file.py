"""
CD4CD8_longitudional_tracking.py

AIM+ (S-reactive clones were previously extracted from Cell Ranger outputs)

Here we produce a file that tracks each AIM+ clones longitudinal trajectories in earlier and later timepoints
based on matching TRB in participant matched deep bulk repertoires. 
"""
import os 
import pandas as pd
import numpy as np
import json
import os 
"""
github.com/kmayerbl/hybrid_t_cell_response
"""
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe
"""
We combine 17 "phenotype_unique" tables precomputed in an earlier phase
The 'phenotype_unique.tsv' suffix refers to phenotypes that could 
be matched directly to bulk clone for confirmation and trajectory analysis. 
Only the most abundant TCRab matching a possible TRB sequence is tracked, to 
avoid over-estimating the repertoire frequency of AIM+ cells.
"""
r = 'phenotypes/'
frames = dict()
for f in os.listdir(r):
    if f.endswith('phenotype_unique.tsv'):
        df = pd.read_csv(os.path.join(r,f), sep = "\t")
        df['source_file'] = f
        df['ptid'] = f.split("_")[0]
        df['visit'] = f.split("_")[1]
        df['visit_code']= f.split("_")[3]
        frames[f] = df
all_E03_phenotypes = pd.concat([v for k,v in frames.items() if k.find("E03")!=-1])
# New cell type analysis, should include double positive
"""
Use most current cell typing method to exclude cells with < 10 CD4/CD8 associated UMIs
"""
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
all_E03_phenotypes['cell_type'] = all_E03_phenotypes.apply(lambda r: cd4cd8_score2(score = r['score'], total_cd_counts = r['total_cd_counts']), axis = 1)
assert all_E03_phenotypes['ptid'].value_counts().shape[0] == 17
"""
* based on participant identification number we look up a json file with paths to all of the bulk repertoire files across multiple timepoints
* Here we align the AIM+ seuqences with a table of clonotype abundance at all available timepoints. 
* We limit the output to stringently enriched AIM+ cells
* We assemble a long form dataframe for use in plotting.
"""
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

    # remove 10X cells which may not be enriched in spike specific pool, applying most stringent criteria
    df_celltype = df_celltype[(df_celltype['10x_cdf_p_fdr_10x'] < 0.05) & (df_celltype['10x_cdf_p_fdr_10x'].notna())]
    # confirm that we are using 10x enriched clones only
    print(df_celltype['10x_enriched'].all())
    assert df_celltype['10x_enriched'].all()
    dfp = df_celltype.copy()
    dfp['plot_clone_id'] = [i for i in range(dfp.shape[0])]
    # Get the value vars based on sample names
    value_vars =  [f"{x}_pfreq" for x in samples.keys()]
       
    print(value_vars)
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
  
"""We combine intermediate output into a single long dataframe"""
all_long_data_used = pd.concat(long_data_frames)
all_long_data_used['j_a_gene'] = all_long_data_used['j_a_gene'].apply(lambda x : f"{x}*01")
all_long_data_used['v_a_gene'] = all_long_data_used['v_a_gene'].apply(lambda x : f"{x}*01")
all_long_data_used.to_csv('2022_09_13_all_long_data.tsv', sep = "\t", index = False)
all_summary_data = pd.concat(summary_data)
all_summary_data.to_csv('2022_09_13_all_summary_data.tsv', sep = "\t", index = False)

"""Visualization is performed in ggplot using script (plot/CD4_CD8_longitudinal_17_ptids.R)"""