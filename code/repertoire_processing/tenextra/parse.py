import os
import pandas as pd 
from functools import partial 

"""
Feb 28, 2022 

Assumption. Each T Cell has one physiologically 
relevant A:B TCR receptor. 

We consider the frequency of each unique chain (nucleotide-level), and 
the frequency of alpha-beta pairings within a pool of single cells.
We attempt to use that information to distinguish 'true' pairings from 
possible cell-free DNA artifacts. 

In those cases, if a cell has multiple alphas/betas pairings,
and no obvious contaminating receptor is found, the procedure 
then is to select, per single cell barcode, the A:B chains 
with the highest umi counts. 

V1. This script was written without option for multiple batches 
TODO: add a multibatch option ['batch_id', 'pool_id']
"""

def _generate_clone_stats_from_contig_annotation_file(
    f,
    p = None,
    sample_index_cols = ['arm', 'ptid', 'visit', 'stim', 'sortid'], 
    cell_index_cols = ['barcode', 'raw_clonotype_id'],
    tcra_cols       = ['v_a_gene', 'j_a_gene', 'cdr3_a_aa'],
    tcrb_cols       = ['v_b_gene', 'j_b_gene', 'cdr3_b_aa'],
    pairs_cols      = ['cdr3_a_nt', 'cdr3_b_nt'],
    chains_cols     = ['cdr3_nt', 'v_gene'],
    tcr_c_cols      = ['v_a_gene', 'j_a_gene', 'cdr3_a_aa', 'v_b_gene', 'j_b_gene', 'cdr3_b_aa', 'cdr3_a_nt', 'cdr3_b_nt']):
    """ 
    Parameters
    ----------
    p : str 
        path to file (optional)
    f : str
        filename for contig annoations file
    cell_index_cols : list
        list of columns that index each single cell (i.e. barcode)
    tcra_cols : list
        list of columns that index tcr alpha chain
    tcrb_cols : list
        lost of columns that index tcr beta chain 
    pairs_cols : list 
        list of columns used in a groupby to how often a particular alpha/beta pair occurs
    chains_cols : list
        list of columns used in groupby to count how often a particular chain appears
    

    Returns
    -------
    tr : pd.DataFrame 
        Dataframe with TRA and TRB joined together on <cell_index_cols>. 
    ct_chains : pd.DataFrame
        Dataframe with number of time each single chain occured
    ct_pairs : pd.DataFrame 
        Dataframe with the number of time each A/B pairing occured
    
    Notes: 

    In effect, we are joining all possible A and B receptors. Subsequently
    we append important columns
    'cts_a' : how often we see a particular alpha/pairing
    'cts_b'

    """
    if p is not None:
        filename = os.path.join(p, f)
    else:
        filename = f
    tcr_raw = pd.read_csv(filename, sep = ",")

    tra = tcr_raw.loc[tcr_raw['chain'] == 'TRA']
    trb = tcr_raw.loc[tcr_raw['chain'] == 'TRB']
    tr = pd.merge(tra, trb, on=cell_index_cols, how='outer', suffixes=('_a', '_b'))
    ren = {'v_gene_b':'v_b_gene',
        'v_gene_a':'v_a_gene',
        'j_gene_b':'j_b_gene',
        'j_gene_a':'j_a_gene',
        'cdr3_b':'cdr3_b_aa',
        'cdr3_a':'cdr3_a_aa',
        'cdr3_nt_b':'cdr3_b_nt',
        'cdr3_nt_a':'cdr3_a_nt'}
    tr = tr.rename(ren, axis=1)
    multis = tr.groupby(cell_index_cols)[['cdr3_a_nt', 'cdr3_b_nt']].nunique()
    multis.columns = ['cts_a', 'cts_b']
    tr = pd.merge(tr, multis, left_on=multis.index.names, right_index=True)

    ct_pairs = tr.groupby(pairs_cols)['barcode'].nunique().sort_values(ascending=False)
    ct_chains = tcr_raw.groupby(chains_cols)['barcode'].nunique().sort_values(ascending=False)
    ct_chains.name = 'ct_chains'
    ct_pairs.name = 'ct_pairs'

    return (tr, ct_chains, ct_pairs)

def _clean_multi_inframe(
    gby, 
    chain, 
    ct_chains, 
    ct_pairs, 
    threshold_chains = 10,
    threshold_pairs = 5):
    """
    gby : slice of DataFrame 
        defined by a groupby operation, such as , e.g., ['barcode', 'raw_clonotype_id'],
    chain : str
        'a' or 'b
    ct_chains : pd.DataFrame (names=['cdr3_nt', 'v_gene')
        This DataFrame provide info how oftem each individual receptor occurs
    ct_pairs : pd.DataFrame ( names=['cdr3_a_nt', 'cdr3_b_nt']), 
        This DataFrame provide info how often each receptor occurs
    threshold_chains : int
        When more than one receptor is present, if the specific receptor 
        is seen fewere more than <threshold_chains> times in the ct_chains lookup OR 
        if its seen at least n - <threshold_pairs> times as a pair then its a real pair.
    threshold_pairs : int
        When more than one receptor is present, if the specific receptor 
        is seen more than max_n times in the ct_chains lookup OR 
        if its seen at least n - <threshold_pairs> times as a pair then its a real pair.
    """
    assert chain in ['a','b']
    # <gby> is a DataFrame slice with all possible TRB:TRA combinations
    # <tmp> is left merge to ct_chains, adding the number of times the primary chain, receptor appears in the full batch
    tmp = pd.merge(gby, ct_chains, left_on=[f'cdr3_{chain}_nt', f'v_{chain}_gene'], right_index=True, how='left')
    # <tmp> a sceond left join to pairs adds counts of the specific pair
    tmp = pd.merge(tmp, ct_pairs, left_on=['cdr3_a_nt', 'cdr3_b_nt'], right_index=True, how='left')
    
    if ((tmp['ct_chains'] < threshold_chains ) & ((tmp['ct_chains'] < threshold_chains ) | tmp['ct_chains'].isnull())).all():
        """If all chains are seen <10 times as a chain then just pick the one with more UMIs, reads"""
        return tmp.sort_values(by=[f'umis_{chain}', f'reads_{chain}'], ascending=False).iloc[0]
    else:
        """If one of them is seen more than n > 10 times as a chain, if its seen at least n - 5 times
        as a pair then its a real pair."""
        ind = (tmp['ct_chains'] < threshold_chains ) | (tmp['ct_pairs'] >= (tmp['ct_chains'] - threshold_pairs))
        
        if ind.sum() >= 1:
            return tmp.loc[ind].sort_values(by=[f'umis_{chain}', f'reads_{chain}'], ascending=False).iloc[0]
        else:
            print(tmp[['cdr3_b_aa', 'cdr3_a_aa', 'reads_a', 'umis_a', 'ct_chains', 'ct_pairs']])
            return None

def _clean_multi_ab_inframe(gby, 
    ct_chains,
    ct_pairs,
    threshold_chains = 10,
    threshold_pairs = 5):
    """
    gby : slice of DataFrame 
        defined by a groupby operation, such as , e.g., ['barcode', 'raw_clonotype_id'],
    ct_chains : pd.DataFrame (names=['cdr3_nt', 'v_gene')
        This DataFrame provide info how oftem each individual receptor occurs
    ct_pairs : pd.DataFrame ( names=['cdr3_a_nt', 'cdr3_b_nt']), 
        This DataFrame provide info how often each receptor occurs
    """
    tmp = pd.merge(gby, ct_chains, left_on=['cdr3_a_nt', 'v_a_gene'], right_index=True, how='left')
    tmp = tmp.rename({'ct_chains':'ct_chains_a'}, axis=1)
    tmp = pd.merge(tmp, ct_chains, left_on=['cdr3_b_nt', 'v_b_gene'], right_index=True, how='left')
    tmp = tmp.rename({'ct_chains':'ct_chains_b'}, axis=1)
    tmp = pd.merge(tmp, ct_pairs, left_on=['cdr3_a_nt', 'cdr3_b_nt'], right_index=True, how='left')
    """If either chain is very common and the pair is uncommon then kick it out"""
    ind = ((tmp['ct_chains_a'] >= threshold_chains) | (tmp['ct_chains_b'] >= threshold_chains)) & ((tmp['ct_pairs'] < (tmp['ct_chains_a'] - threshold_pairs)) | (tmp['ct_pairs'] < (tmp['ct_chains_b'] - 5)))
    tmp = tmp.loc[~ind]
    if (~ind).sum() >= 1:
        return tmp.loc[~ind].sort_values(by=['umis_a', 'umis_b', 'reads_a', 'reads_b'], ascending=False).iloc[0]
    else:
        return None

def _all_multi_ab_inframe(gby, 
    ct_chains,
    ct_pairs):
    """
    gby : slice of DataFrame 
        defined by a groupby operation, such as , e.g., ['barcode', 'raw_clonotype_id'],
    ct_chains : pd.DataFrame (names=['cdr3_nt', 'v_gene')
        This DataFrame provide info how oftem each individual receptor occurs
    ct_pairs : pd.DataFrame ( names=['cdr3_a_nt', 'cdr3_b_nt']), 
        This DataFrame provide info how often each receptor occurs
    """
    tmp = pd.merge(gby, ct_chains, left_on=['cdr3_a_nt', 'v_a_gene'], right_index=True, how='left')
    tmp = tmp.rename({'ct_chains':'ct_chains_a'}, axis=1)
    tmp = pd.merge(tmp, ct_chains, left_on=['cdr3_b_nt', 'v_b_gene'], right_index=True, how='left')
    tmp = tmp.rename({'ct_chains':'ct_chains_b'}, axis=1)
    tmp = pd.merge(tmp, ct_pairs, left_on=['cdr3_a_nt', 'cdr3_b_nt'], right_index=True, how='left')
    return tmp

def _all_multi_inframe(
    gby, 
    chain, 
    ct_chains, 
    ct_pairs):
    """
    gby : slice of DataFrame 
        defined by a groupby operation, such as , e.g., ['barcode', 'raw_clonotype_id'],
    chain : str
        'a' or 'b
    ct_chains : pd.DataFrame (names=['cdr3_nt', 'v_gene')
        This DataFrame provide info how oftem each individual receptor occurs
    ct_pairs : pd.DataFrame ( names=['cdr3_a_nt', 'cdr3_b_nt']), 
        This DataFrame provide info how often each receptor occurs
    """
    # <gby> is a DataFrame slice with all possible TRB:TRA combinations
    # <tmp> is left merge to ct_chains, adding the number of times the primary chain, receptor appears in the full batch
    tmp = pd.merge(gby, ct_chains, left_on=[f'cdr3_{chain}_nt', f'v_{chain}_gene'], right_index=True, how='left')
    # <tmp> a sceond left join to pairs adds counts of the specific pair
    tmp = pd.merge(tmp, ct_pairs, left_on=['cdr3_a_nt', 'cdr3_b_nt'], right_index=True, how='left')
    return tmp


def select_likely_receptors( f,
                             p = None, 
                             threshold_chains = 10, 
                             threshold_pairs = 5):

    tr, ct_chains, ct_pairs = _generate_clone_stats_from_contig_annotation_file(p = p , f =f )
    # 0 <multi_ind> True if their are either more than 1 TRA or TRB receptor
    multi_ind = (tr[['cts_a', 'cts_b']] > 1).any(axis=1)
    # <single> DataFrame of all clones without duplicate chains
    singles = tr.loc[~multi_ind]
    # DataFrame <multi_a> multiple alpha chains
    multi_a = tr.loc[(tr['cts_a'] > 1) & (tr['cts_b'] <= 1)]
    # DataFrame <multi_b> multiple beta chains
    multi_b = tr.loc[(tr['cts_b'] > 1) & (tr['cts_a'] <= 1)]
    # DataFrame with multiple beta and alpha chains
    multi_ab = tr.loc[(tr['cts_a'] > 1) & (tr['cts_b'] > 1)]

    cell_index_cols = ['barcode', 'raw_clonotype_id']
    
    clean_a = multi_a.groupby(cell_index_cols).apply(partial(_clean_multi_inframe, 
        chain='a', 
        ct_chains = ct_chains, 
        ct_pairs = ct_pairs,
        threshold_chains = threshold_chains ,
        threshold_pairs = threshold_pairs )).reset_index(drop=True)
    
    clean_b = multi_b.groupby(cell_index_cols).apply(partial(_clean_multi_inframe,
        chain='b', 
        ct_chains = ct_chains, 
        ct_pairs = ct_pairs,
        threshold_chains = threshold_chains,
        threshold_pairs = threshold_pairs )).reset_index(drop=True)

    clean_ab = multi_ab.groupby(cell_index_cols).apply(_clean_multi_ab_inframe, 
        ct_chains = ct_chains, 
        ct_pairs = ct_pairs,
        threshold_chains = threshold_chains,
        threshold_pairs = threshold_pairs).reset_index(drop=True)
    
    clean_clones     =  pd.concat([singles, clean_a, clean_b, clean_ab], sort = True)

    """For diagnostic purposes"""
    all_clones = tr.groupby(cell_index_cols).apply(_all_multi_ab_inframe, 
        ct_chains = ct_chains, 
        ct_pairs = ct_pairs).reset_index(drop=True)


    return clean_clones, all_clones, ct_chains, ct_pairs






