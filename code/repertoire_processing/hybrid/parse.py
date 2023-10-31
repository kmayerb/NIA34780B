from tcrdist.swap_gene_name import adaptive_to_imgt
import pandas as pd

def parse_adaptive_v2(df, use_max_resolved = True):
    """
    parse_adpative_v2
    
    Parameters
    ---------
    df : pd.DataFrame
        v2 version of Immune Access Adaptive Data
    use_max_resolved : bool
        if True, use the vMaxResolved and jMaxResolved columns for gene calling 

    Returns 
    -------
    dfx pd.DataFrame
    """
    cols = {'nucleotide': 'cdr3_b_nucseq',
        'aminoAcid': 'cdr3_b_aa',
        'count (templates/reads)': 'templates',
        'frequencyCount (%)': 'frequency',
        'vGeneName': 'v_b_gene_adpt',
        'jGeneName': 'j_b_gene_adpt',
        'vMaxResolved' : 'vMaxResolved',
        'jMaxResolved' : 'jMaxResolved'}

    df = df[df['sequenceStatus']=='In'].reset_index(drop = True)
    df = df[cols].rename(columns = cols)

    if use_max_resolved: 
        df['v_b_gene'] = df['vMaxResolved'].apply(lambda x: adaptive_to_imgt['human'].get(x.split("*")[0]) if isinstance(x,str) else None)
        df['j_b_gene'] = df['jMaxResolved'].apply(lambda x: adaptive_to_imgt['human'].get(x.split("*")[0])  if isinstance(x,str) else None)
    else:
        df['v_b_gene'] = df['v_b_gene_adpt'].apply(lambda x: adaptive_to_imgt['human'].get(x) if isinstance(x,str) else None)
        df['j_b_gene'] = df['j_b_gene_adpt'].apply(lambda x: adaptive_to_imgt['human'].get(x) if isinstance(x,str) else None)

    dfx = df[df['cdr3_b_aa'].notna()].reset_index(drop = True)
    dfx['v_family'] = dfx['v_b_gene'].apply(lambda x : x.split("*")[0].split("-")[0]  if isinstance(x,str) else None)
    dfx['j_family'] = dfx['j_b_gene'].apply(lambda x : x.split("*")[0].split("-")[0]  if isinstance(x,str) else None)
    dfx['productive_frequency'] = dfx['templates'] / dfx['templates'].sum()
    return dfx