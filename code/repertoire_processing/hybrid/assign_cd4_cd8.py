"""
assign_cd4_cd8.py
Barcode indexes a single droplet which have cell surface marker information
Clonotypes will map to many barcodes, so it may be useful to 
take distributional statistic of cell marker information.
"""
import os
import json
import scanpy as sc
import pandas as pd
import numpy as np
from tenextra.parse import select_likely_receptors
from hybrid.parse import parse_adaptive_v2
import statsmodels.stats.multitest as mt
from scipy.stats import binom
from statsmodels.stats.multitest import multipletests

def get_feature_names(matrix_h5):
    adata = sc.read_10x_h5(matrix_h5, gex_only= False)
    var_select = [name for name in adata.var_names if not name.startswith('Hash')]
    return(var_select)

def simple_call(x):
    if x >= 1:
        return "CD8"
    elif x <= -1:
        return "CD4"
    else:
        return None
       
def call_CD4_CD8_per_barcode(matrix_h5, var_select = ['CD4_TotalSeqC', 'CD8a_TotalSeqC']):
    """
    Here we make CD4 and CD8 calls based on antibody markers
    """
    adata = sc.read_10x_h5(matrix_h5, gex_only= False)
    adata.var_names_make_unique()
    adata = adata[:, var_select] 
    adata.var['cd8'] = adata.var_names.str.startswith('CD8a_TotalSeqC')
    adata.var['cd4'] = adata.var_names.str.startswith('CD4_TotalSeqC')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['cd8'], percent_top=None, log1p=False, inplace=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['cd4'], percent_top=None, log1p=False, inplace=True)
    obs_df = adata.obs.reset_index(drop = False).rename(columns = {'index':'barcode'})
    obs_df['cd8vcd4_score'] = obs_df.apply(lambda r: np.log(r['pct_counts_cd8']+1)-np.log(r['pct_counts_cd4']+1) , axis = 1)   
    obs_df['cell_type'] = obs_df['cd8vcd4_score'].apply(lambda x: simple_call(x))
    return obs_df

def fdr_valid_with_nas(pvals):
    """
    Parameters
    ----------
    pvals : 
        list or array of pvalues
    f1 : pd.Series 
        is the frequency in sample 1, used for checking detection > 0
    f2 : pd.Series 
        is the frequency in sample 2, , used for checking detection > 0

    Returns
    -------
    pd.Series with FDR values for the clones that are subject to valid comparison 

    Notes
    -----
    Get an fdr that is valid for two sampel comparison. 
    Avoid penalty of comparing all possible clones even those
    not seen in that sample
    """
    ix = pvals.notna()
    fdr_valid = mt.multipletests(pvals = pd.Series(pvals)[ix], method = "fdr_bh" )[1]
    df = pd.DataFrame({'pvals':pvals})
    df['fdr'] = np.NaN
    df.loc[ix,'fdr'] = fdr_valid
    return(df['fdr'])



def get_information_from_10X(
    clones_file_10X,
    matrix_file_10X,
    clones_file_adpt):
    
    # Here we deal with potential duplicates
    print("")
    clean_clones, all_clones, ct_chains, ct_pairs = \
        select_likely_receptors(
            f = clones_file_10X ,
            threshold_chains = 10)

    # we pull cd4 and cd8 barcodes
    cd4_cd8_barcodes = call_CD4_CD8_per_barcode(
        matrix_h5  = matrix_file_10X, 
        var_select = ['CD4_TotalSeqC', 'CD8a_TotalSeqC'])

    # Merge clean clones with total seq
    clean_clones_total_seq = clean_clones.merge(cd4_cd8_barcodes, how = "left", on = "barcode")
    icols = ['cdr3_a_nt', 'v_a_gene','j_a_gene','cdr3_a_aa', 'cdr3_b_nt','v_b_gene','j_b_gene','cdr3_b_aa']
    sumcols = ['total_counts_cd8','pct_counts_cd8' ,'total_counts_cd4','pct_counts_cd4','cd8vcd4_score']
    # store median value across all droplets form same clone
    a = clean_clones_total_seq.groupby(icols)[sumcols].median()
    # store barcodes counts in b
    b = clean_clones_total_seq.groupby(icols)['barcode'].size()
    clean_clones_total_seq_median = pd.concat([a,b], axis = 1).sort_values('barcode', ascending = False).reset_index(drop = False)
    clean_clones_total_seq_median['cell_type'] = clean_clones_total_seq_median['cd8vcd4_score'].apply(lambda x: simple_call(x))

    """Define a key to match TRBV-CDR3(AA)-TRBJ"""
    cc = clean_clones_total_seq_median
    cc['key'] = cc['v_b_gene'] + "+" +  cc['cdr3_b_aa'] + "+" + cc['j_b_gene']
    cc['N_10x'] = cc['barcode'].sum()
    cc['10x_count'] = cc['barcode']
    cc['10x_clone_id'] = [i for i in range(cc.shape[0])]
    """Load"""

    
    bulk_clones = parse_adaptive_v2(pd.read_csv(clones_file_adpt, sep = "\t"))
    bulk_clones = bulk_clones[(bulk_clones['v_b_gene'].notna()) & (bulk_clones['j_b_gene'].notna())].reset_index(drop = True)

    bulk_clones['v_b_gene'] = bulk_clones['v_b_gene'].apply(lambda x: x.split("*")[0])
    bulk_clones['j_b_gene'] = bulk_clones['j_b_gene'].apply(lambda x: x.split("*")[0])

    bulk_clones = bulk_clones.groupby(['cdr3_b_nucseq','v_b_gene','j_b_gene','cdr3_b_aa']).sum().reset_index()
    bulk_clones['adpt_clone_id'] = [i for i in range(bulk_clones.shape[0])]
    bulk_clones['key'] = bulk_clones['v_b_gene'] + "+" +  bulk_clones['cdr3_b_aa'] + "+" + bulk_clones['j_b_gene']
    cc_bulk = cc.merge(bulk_clones, how = "left", on = 'key')
    cc_bulk['nucleotide_level_match'] = \
        cc_bulk.apply(lambda r: r['cdr3_b_nucseq'].find(r['cdr3_b_nt']) != -1 if isinstance(r['cdr3_b_nt'], str) and isinstance(r['cdr3_b_nucseq'], str)  else False, axis = 1)
    cc_bulk = cc_bulk.query('nucleotide_level_match')

    """Subset rows that have a match in the bulk data"""
    """ Compute the productive frequency so we can compare to 10x"""
    cc_bulk['10x_pfreq'] = cc_bulk['10x_count'] / cc_bulk['barcode'].sum()
    cc_bulk[['10x_count', 'N_10x', 'productive_frequency','10x_pfreq']]
    cc_bulk['10x_enriched'] = cc_bulk.apply(lambda r: r['10x_pfreq']>r['productive_frequency'], axis =1) 

    cc_bulk['10x_binomial_cdf']    = cc_bulk.apply(lambda r : 1-binom.cdf( k =r['10x_count']-1, n = r['N_10x'], p  = r['productive_frequency']), axis  =1 )
    cc_bulk['10x_binomial_pmf']    = cc_bulk.apply(lambda r : binom.pmf( k =r['10x_count'], n = r['N_10x'], p  = r['productive_frequency']), axis  =1 )
    cc_bulk['10x_cdf_p_fdr_10x']   = fdr_valid_with_nas( pvals = cc_bulk['10x_binomial_cdf'])
    cc_bulk['10x_pmf_p_fdr_10x']   = fdr_valid_with_nas( pvals = cc_bulk['10x_binomial_pmf'])

    # take the most abundant 10x clone that matches 
    cc_bulk_unique = cc_bulk.sort_values('10x_count', ascending = False).groupby('adpt_clone_id').head(1).reset_index(drop = True)
    cc_bulk_unique['clones_file_10X'] =    os.path.basename(clones_file_10X)
    cc_bulk_unique['matrix_file_10X'] =    os.path.basename(matrix_file_10X)
    cc_bulk_unique['clones_file_adpt'] =    os.path.basename(clones_file_adpt)
    return cc_bulk_unique, cc_bulk







