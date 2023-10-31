"""
assemble_longitudinal_sample.py

"""
import pandas as pd
import os

import fishersapi
import statsmodels.stats.multitest as mt
import numpy as np
from hybrid.parse import parse_adaptive_v2


def fdr_valid(f1, f2, pvals):
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
    ix = (f1>0)|(f2>0)
    fdr_valid = mt.multipletests(pvals = pd.Series(pvals)[ix], method = "fdr_bh" )[1]
    df = pd.DataFrame({'f1':f1, 'f2':f2})
    df['fdr'] = np.NaN
    df.loc[ix,'fdr'] = fdr_valid
    return(df['fdr'])

def make_wide_longitudinal_dataframe(sample_dictionary, dest = ".", ptid = None, write = False, fdr_thr = 0.05):
    """
    Parameters
    ----------
    sample_dictionary : dict
        dictionary keyed on sample names, values are full filepaths to
        adaptive format
    dest : str
        place to write out if write is True 
    ptid : str or None
        not required
    fdr_thr : float 
        the level to use to determine if a clone has significantly expanded or contracted 
        default is 0.05 to accord with B>A analyses.

    Returns
    -------
    dft : pd.DataFrame
        Wide format longitudinal data
    """
    tcrdist_frames = dict()
    samples_in_order = list()
    for sample, fp in sample_dictionary.items():
        samples_in_order.append(sample)
        print(f"PARSING {sample} AT {fp}")
        assert os.path.isfile(fp)
        df = pd.read_csv(fp, sep = "\t")
        print(f"RAW DATAFRAME({os.path.basename(fp)}) -- {df.shape[0]} ROWS")
        try:
            df = parse_adaptive_v2(df)
        except KeyError:
            print("Adaptive V2 Parse Failed, Trying Earlier Version Parsing")
            df = parse_adaptive(df)
        df['sample_name'] = sample
        print(df.shape)
        print(f"PARSED DATAFRAME({os.path.basename(fp)}) -- {df.shape[0]} ROWS")
        tcrdist_frames[sample] = df

    dfall = pd.concat(tcrdist_frames.values())  
    dfall_s = dfall.groupby(['cdr3_b_nucseq', 'cdr3_b_aa',  'v_b_gene', 'j_b_gene', 'v_family','j_family', 'sample_name','vMaxResolved','jMaxResolved']).sum()
    #dfall_s = dfall.groupby(['cdr3_b_nucseq', 'cdr3_b_aa',  'v_b_gene', 'j_b_gene', 'v_family','j_family', 'sample_name']).sum()
    dfall_s = dfall_s.reset_index(drop = False)
    """<dft> This is a DataFrame with template counts, each timepoint as a unique sample"""
    dft = dfall_s.pivot(index=['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family','j_family', 'vMaxResolved','jMaxResolved'],
        columns='sample_name',values='templates')

    """<cs> numeric columns refering to number of templates in each sample, stored before we reset index"""
    # ensure we preserve the order from the json
    dft = dft[samples_in_order]
    cs = dft.columns.to_list()
    dft = dft.reset_index(drop = False)

    """replace na with zero in numeric template clolumns"""
    for col in cs:
        dft[col] = dft[col].fillna(0)

    """Compute pfreq (productive frequency) based on template divided by templates sum"""
    for col in cs:
        dft[f"{col}_pfreq"] = dft[col] / dft[col].sum()

    """This is an intermediate output, prior to addding additional information"""
    if write:
        outfile1 = os.path.join(dest, f"{ptid}.longitudinal_templates_only.tsv")
        print(f"OUTPUT {outfile1}")
        dft.to_csv(outfile1, sep = "\t", index = False)
    
    return dft

def compute_fold_changes_on_longitudinal_dataframe(dft, dest = ".", ptid = None, write = False, fdr_thr = 0.05):
    """
    Parameters
    ----------
    dft : pd.DataFrame
        Wide format longitudinal data from make_wide_longitudinal_dataframe().
        Rows are unique TRBV-CDR3b(nt) and columns are templates counts per
        samplename.
    Returns

    dft : pd.DataFrame

    """
    # We need to get the names of the samples columns, so we excluded index an _pfreq columns    
    index_cols = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family',
       'j_family', 'vMaxResolved', 'jMaxResolved']
    cs = [col for col in dft.columns if col not in index_cols]
    cs = [col for col in cs if col.find("_pfreq") == -1]

    """Compute templates sum"""
    for col in cs:
        dft[f"{col}_sum"] = dft[col].sum()

    """Temporal Fold Change is FC from one period to the next, what is the obesrved fold change from one period to the next"""
    psuedo = 1E-6
    for i in range(len(cs)-1):
        col_tupple = (cs[i+1], cs[i])
        print(col_tupple)
        dft[f'{cs[i+1]}_temporal_fc'] =  (dft[f'{cs[i+1]}_pfreq'] + psuedo) / (dft[f'{cs[i]}_pfreq'] + psuedo)  

    """Vaccine associated fold change is relative change to a pre-vaccine timepoint"""
    """add a small psuedo count to avoid division by zero"""    
    psuedo = 1E-6
    e01_column = [c for c in cs if c.find("E01") != -1]
    for i in range(len(cs)):
        col_tupple = ( e01_column[0], cs[i])
        print(col_tupple)
        dft[f'{cs[i]}_vac_fc'] =  (dft[f'{cs[i]}_pfreq'] + psuedo) / (dft[f'{e01_column[0]}_pfreq'] + psuedo)

    """To test for significance, perfor Fisher's Exact Tests of change relative to E01 (pre-vac) timepoint"""
    e01_column = [c for c in cs if c.find("E01") != -1]
    for i in range(len(cs)):
        print(f"COMPUTING EXPANSION FISHER'S EXACT TEST for {cs[i]} and {e01_column[0]} ")
        a = dft[f'{cs[i]}']
        b = dft[f'{cs[i]}_sum'] - dft[f'{cs[i]}']
        c = dft[f'{e01_column[0]}']
        d = dft[f'{e01_column[0]}_sum'] - dft[f'{e01_column[0]}']
        print("COMPUTING FISHER'S STATS WITH FISHERSAPI")
        odds, p = fishersapi.fishers_vec(a,b,c,d)
        dft[f'{cs[i]}_vac_OR']   = odds
        dft[f'{cs[i]}_vac_pval'] = p
        dft[f'{cs[i]}_vac_fdr']  = fdr_valid(f1 = dft[f"{cs[i]}"], f2 =dft[f'{e01_column[0]}'] , pvals = p)
        #fdrmt.multipletests(pvals = p, method = "fdr_bh" )[1]

    """Further test significantce of expansions or contractions from one period to another"""
    for i in range(len(cs)-1):
        print(f"COMPUTING TEMPORAL ONE PERIOD EXPANSION FISHER'S EXACT TEST for {cs[i+1]} and {cs[i]} ")
        a = dft[f'{cs[i+1]}']
        b = dft[f'{cs[i+1]}_sum'] - dft[f'{cs[i+1]}']
        c = dft[f'{cs[i]}']
        d = dft[f'{cs[i]}_sum'] - dft[f'{cs[i]}']

        print("COMPUTING FISHER'S STATS WITH FISHERSAPI")
        odds, p = fishersapi.fishers_vec(a,b,c,d)
        dft[f'{cs[i+1]}_temporal_OR']   = odds
        dft[f'{cs[i+1]}_temporal_pval'] = p
        dft[f'{cs[i+1]}_temporal_fdr']  = fdr_valid(f1 = dft[f"{cs[i]}"], f2 =dft[f'{cs[i+1]}'] , pvals = p)
        #fdrmt.multipletests(pvals = p, method = "fdr_bh" )[1]

    """ Finally let's dichotomize detection, providing a columns of whether the clone was even detected at each timepoint"""
    for i in range(len(cs)):
        dft[f'{cs[i]}_detect'] = dft[f'{cs[i]}'] > 0 

    
    """Finally we label a classification for the type of change"""
    print("CLASSIFYING SIGNIFICANT EXPANDERS AND CONTRACTORS RELATIVE TO PRIOR TIMEPOINT)")
    for i in range(1,len(cs)):
        col = cs[i]
        sig_exp       = (dft[f'{col}_temporal_OR'] > 1) & (dft[f'{col}_temporal_fdr'] < fdr_thr)
        sig_contract  = (dft[f'{col}_temporal_OR'] < 1) & (dft[f'{col}_temporal_fdr'] < fdr_thr)
        no_sig_change = (dft[f'{col}_temporal_fdr'] > fdr_thr) | (dft[f'{col}_temporal_fdr'].isna())
        detect     = dft[f'{col}_detect']
        import pandas as pd
        dftemp = pd.DataFrame({'a':sig_exp, 'b':sig_contract, 'c':no_sig_change, 'd':detect})
        #df.apply(lambda r: classify_dict.get((r['a'],r['b'],r['c'],r['d'])), axis = 1)
        classify_dict = {
        (True, False, False, True) : 'sig_expand',
        (False, True, False, True) : 'sig_contract',
        (False, True, False, False) : 'sig_contract',
        (False, False, True, True) : 'no_sig_change',
        (False, False, True, False) : 'no_detect'}
        dft[f'{col}_temporal_class'] = pd.Series([classify_dict.get(y) for y in [tuple(x) for x in dftemp.to_dict('split')['data'] ]])
    print("CLASSIFYING SIGNIFICANT EXPANDERS AND CONTRACTORS RELATIVE TO PRE-VACCINE (E01)")
    for i in range(len(cs)):
        col = cs[i]
        sig_exp       = (dft[f'{col}_vac_OR'] > 1) & (dft[f'{col}_vac_fdr'] < fdr_thr)
        sig_contract  = (dft[f'{col}_vac_OR'] < 1) & (dft[f'{col}_vac_fdr'] < fdr_thr)
        no_sig_change = (dft[f'{col}_vac_fdr'] > fdr_thr) | (dft[f'{col}_vac_fdr'].isna())
        detect     = dft[f'{col}_detect']
        import pandas as pd
        dftemp = pd.DataFrame({'a':sig_exp, 'b':sig_contract, 'c':no_sig_change, 'd':detect})
        #df.apply(lambda r: classify_dict.get((r['a'],r['b'],r['c'],r['d'])), axis = 1)
        classify_dict = {
        (True, False, False, True) : 'sig_expand',
        (False, True, False, True) : 'sig_contract',
        (False, True, False, False) : 'sig_contract',
        (False, False, True, True) : 'no_sig_change',
        (False, False, True, False) : 'no_detect'}
        dft[f'{col}_vac_class'] = pd.Series([classify_dict.get(y) for y in [tuple(x) for x in dftemp.to_dict('split')['data'] ]])
    
    """ We've now created a very large matrix with information that can be vizualized, this will be a 1GB before compression"""
    if write:
        outfile = os.path.join(dest, f"{ptid}.longitudinal_matrix.tsv")
        print(f"OUTPUT {outfile}")
        dft.to_csv(outfile, sep = "\t", index = False)
    
    return dft
