"""
feasible_hla_walk.py
Reviewed on June 13, 2022, kmayerbl
"""
from inspect import classify_class_attrs
import os 
import sys
import pandas as pd
import numpy as np
import networkx as nx
from tcrdist.repertoire import TCRrep

def tcr_ab_to_network(
    tr, 
    edge_threshold = 120,
    columns_to_track = ['person_x', 
        'cdr3_a_aa', 'cdr3_b_aa','v_a_gene',
        'v_b_gene','j_a_gene','j_b_gene',
        'cell_type']):
    """
    Takes a TCRrep instance with paired alpha-beta data and returns a network DataFrame

    Parameters
    ----------
    tr : TCRrep
        tcrdist3 TCRrep instance, must contain all columns
        in <columns_to_track>
    edge_threshold : int
        TCRdist threshold for (default : 120)
    column_to_track : list
        List of columns to include df_net work DataFrame

    Returns
    -------
    df_net : pd.DataFrame

    G : nx.Graph

    """
    from tcrdist.public import _neighbors_fixed_radius
    import networkx as nx
    # Compute Network
    # <tr.pw_alpha_beta> is paired chain TCRdist.
    tr.pw_alpha_beta = tr.pw_beta + tr.pw_alpha
    # <network> initialize a list to populate with edges between TCRs.
    network = list()
    for i,n in enumerate(_neighbors_fixed_radius(tr.pw_alpha_beta, edge_threshold)):
        for j in n:
            if i != j:
                row = list()
                row.append(i)
                row.append(j)
                row.append(len(n)-1)
                row.append((tr.pw_alpha_beta )[i,j])
                for col in columns_to_track:
                    row.append(tr.clone_df[col].iloc[i])
                    row.append(tr.clone_df[col].iloc[j])
                network.append(row)                      # 'K_neighbors' - number of neighbors
    import numpy as np
    cols = ['node_1', 'node_2','K_neighbors','dist'] + np.concatenate([(f"{x}_1", f"{x}_2") for x in columns_to_track]).tolist()
    df_net = pd.DataFrame(network, columns = cols)
    df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold
    G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],
                                           'target' : df_net['node_2'],
                                           'weight' : df_net['weight']}))
    return(df_net, G)

def get_hla_dictionaries(
    tsv_file = '/fh/fast/gilbert_p/fg_data/koelle_covid/hla_file.tsv',
    cols = ['HLA-Aa', 'HLA-Ab', 'HLA-Ba', 'HLA-Bb', 'HLA-Ca', 'HLA-Cb', 'DQB1b', 'DPA1a', 'DPA1b', 'DPB1a', 'DPB1b'], 
    sample_key = 'sample'):
    """

    Parameters
    ----------
    tsv_file : str
        Full path to hla file which must have all columns in <cols> argument, must include 
    cols : list
        List of HLA column names
    sample_key : str
        the column name for the keys in the dictionary

    Returns
    -------
    sample_hla_dict : dictionary 
        from key -> to set of HLA alleles (A*02:01:01)
    sample_hla_dict_4d : dictionary 
        from key -> to set of HLA alleles (A*02:01)
    sample_hla_dict_4d : dictionary 
        from key -> to set of HLA alleles (A*02)
    
    """
    hla_file = pd.read_csv(tsv_file ,sep ='\t') 
    
    sample_hla_dict = dict()
    sample_hla_dict_2d = dict()
    sample_hla_dict_4d = dict()
    all_hla = list()
    for i,row in hla_file.iterrows():
        all_hla.extend(row[cols].to_list())
        sample_hla_dict[   row[sample_key]] = set(row[cols].to_list())
        sample_hla_dict_2d[row[sample_key]] = set([x.split(":")[0] for x in list(set(row[cols].to_list()))])
        sample_hla_dict_4d[row[sample_key]] = set([":".join(x.split(":")[0:2]) for x in list(set(row[cols].to_list()))])
    return(sample_hla_dict, sample_hla_dict_4d, sample_hla_dict_2d)
    

def compute_hla_feasibility(tr, G, sample_hla_dict, sample_hla_dict_4d, sample_hla_dict_2d, run_class_i= True, run_class_ii= True):
    """
    Add all of the HLA feasibility information
    Parameters
    ----------
    run_class_i: bool
        True, test related to class i HLA
    run_class_ii:
        True, test related to class ii HLA
    
    """
    ccs = [c for c in sorted(nx.connected_components(G), key=len, reverse=True)]
    clone_to_cc = dict()
    for cc in ccs: 
        for i in cc:
            clone_to_cc[i] = cc 
    tr.clone_df['connected_to']  = pd.Series(tr.clone_df.index).apply(lambda x :  clone_to_cc.get(x, {x})) 
    tr.clone_df['K_public_cc'] = tr.clone_df['connected_to'].apply(lambda x : pd.Series([tr.clone_df['person_x'].iloc[i] for i in x]).nunique()) 
    tr.clone_df['has_public_cc'] = tr.clone_df['connected_to'].apply(lambda x : pd.Series([tr.clone_df['person_x'].iloc[i] for i in x]).nunique() > 1) 
    # Index of connected pubic connected components
    tr.clone_df['connected_to_ranked'] = [sorted([(c,tr.pw_alpha_beta[i,c]) for c in cc], key = lambda x: x[1]) for i,cc in enumerate(tr.clone_df['connected_to'])]

    class_i  = lambda alleles :  set([a for a in alleles if a.startswith('A') or a.startswith('B') or a.startswith('C')  ])
    class_ii = lambda alleles :  set([a for a in alleles if not (a.startswith('A') or a.startswith('B') or a.startswith('C') ) ])

    # Class I 


    if run_class_i:
        sample_hla_dict_i    = {sample:class_i(a) for sample,a in sample_hla_dict.items()}
        sample_hla_dict_2d_i = {sample:class_i(a) for sample,a in sample_hla_dict_2d.items()}
        sample_hla_dict_4d_i = {sample:class_i(a) for sample,a in sample_hla_dict_4d.items()}
    
        
        feasible = dict()
        feasible2 = dict()
        feasible4 = dict()
        for i, row in tr.clone_df.iterrows():
            if i % 500 == 0 :
                print(f"COMPUTING HLA CLASS I FEASIBILITY FOR ROW {i} OF {tr.clone_df.shape[0]}")
            try:
                feasible_hla_set = sample_hla_dict_i[row['person_x']].copy()
                feasible_hla_set_2d = sample_hla_dict_2d_i[row['person_x']].copy()
                feasible_hla_set_4d = sample_hla_dict_4d_i[row['person_x']].copy()
            except KeyError:
                feasible[i] = None
                feasible2[i] = None 
                feasible4[i] = None 
            for x in row['connected_to_ranked']:
                person = tr.clone_df.iloc[x[0],]['person_x']
                try:
                    # Don't allow an empty set
                    if not feasible_hla_set.intersection(sample_hla_dict_i[person]) == set():
                        feasible_hla_set = feasible_hla_set.intersection(sample_hla_dict_i[person])

                    if not feasible_hla_set_2d.intersection(sample_hla_dict_2d_i[person]) == set():
                        feasible_hla_set_2d = feasible_hla_set_2d.intersection(sample_hla_dict_2d_i[person])

                    if not feasible_hla_set_4d.intersection(sample_hla_dict_4d_i[person]) == set():
                        feasible_hla_set_4d = feasible_hla_set_4d.intersection(sample_hla_dict_4d_i[person])
                except KeyError:
                    continue
            feasible[i]  = feasible_hla_set 
            feasible2[i] = feasible_hla_set_2d 
            feasible4[i] = feasible_hla_set_4d 

            tr.clone_df['feasible_hla_i'] =    [sorted(list(feasible[x])) for x in tr.clone_df.index]
            tr.clone_df['feasible_hla_2d_i'] = [sorted(list(feasible2[x])) for x in tr.clone_df.index]
            tr.clone_df['feasible_hla_4d_i'] = [sorted(list(feasible4[x])) for x in tr.clone_df.index]
    
    if run_class_ii:
        sample_hla_dict_ii    = {sample:class_ii(a) for sample,a in sample_hla_dict.items()}
        sample_hla_dict_2d_ii = {sample:class_ii(a) for sample,a in sample_hla_dict_2d.items()}
        sample_hla_dict_4d_ii = {sample:class_ii(a) for sample,a in sample_hla_dict_4d.items()}
        # Repeat for class ii
        feasible = dict()
        feasible2 = dict()
        feasible4 = dict()
        for i, row in tr.clone_df.iterrows():
            if i % 500 == 0 :
                print(f"COMPUTING HLA CLASS II FEASIBILITY FOR ROW {i} OF {tr.clone_df.shape[0]}")
            try:
                feasible_hla_set = sample_hla_dict_ii[row['person_x']].copy()
                feasible_hla_set_2d = sample_hla_dict_2d_ii[row['person_x']].copy()
                feasible_hla_set_4d = sample_hla_dict_4d_ii[row['person_x']].copy()
            except KeyError:
                feasible[i] = None
                feasible2[i] = None 
                feasible4[i] = None 
            for x in row['connected_to_ranked']:
                person = tr.clone_df.iloc[x[0],]['person_x']
                try:
                    # Don't allow an empty set
                    if not feasible_hla_set.intersection(sample_hla_dict_ii[person]) == set():
                        feasible_hla_set = feasible_hla_set.intersection(sample_hla_dict_ii[person])

                    if not feasible_hla_set_2d.intersection(sample_hla_dict_2d_ii[person]) == set():
                        feasible_hla_set_2d = feasible_hla_set_2d.intersection(sample_hla_dict_2d_ii[person])

                    if not feasible_hla_set_4d.intersection(sample_hla_dict_4d_ii[person]) == set():
                        feasible_hla_set_4d = feasible_hla_set_4d.intersection(sample_hla_dict_4d_ii[person])
                except KeyError:
                    continue
            feasible[i]  = feasible_hla_set 
            feasible2[i] = feasible_hla_set_2d 
            feasible4[i] = feasible_hla_set_4d 

            tr.clone_df['feasible_hla_ii'] =    [sorted(list(feasible[x])) for x in tr.clone_df.index]
            tr.clone_df['feasible_hla_2d_ii'] = [sorted(list(feasible2[x])) for x in tr.clone_df.index]
            tr.clone_df['feasible_hla_4d_ii'] = [sorted(list(feasible4[x])) for x in tr.clone_df.index]

    return tr

"""
# EXAMPLE 
# -------

import sys
import os 
import pandas as pd
from tcrdist.repertoire import TCRrep

sys.path.insert(0, '/fh/fast/gilbert_p/fg_data/koelle_covid/project_code/hybrid_t_cell_response') 
from hybrid.feasible import tcr_ab_to_network
from hybrid.feasible import get_hla_dictionaries
from hybrid.feasible import compute_hla_feasibility

clones_file = pd.read_csv('/fh/fast/gilbert_p/fg_data/koelle_covid/project_code/hybrid_t_cell_response/all_E03_phenotypes.tsv', sep = "\t")
clones_file['person_x'] = clones_file['source_file'].str.split('_').apply(lambda x: int(x[0]))
df  = clones_file.rename(columns = {'cdr3_b_aa_x':'cdr3_b_aa','v_b_gene_x':'v_b_gene','j_b_gene_x':'j_b_gene','hla_key':'person_x'})
df['v_a_gene'] = df['v_a_gene'].apply(lambda x : f"{x}*01")
df['j_a_gene'] = df['j_a_gene'].apply(lambda x : f"{x}*01")
df['v_b_gene'] = df['v_b_gene'].apply(lambda x : f"{x}*01")
df['j_b_gene'] = df['j_b_gene'].apply(lambda x : f"{x}*01")
df['cell_type'] = df['cell_type'].fillna('UNK')
clone_cols = ['person_x', 'cdr3_a_aa', 'cdr3_b_aa','v_a_gene','v_b_gene','j_a_gene','j_b_gene','cell_type']
df[clone_cols]
# Load a TCRrep insance:
tr = TCRrep(cell_df = df[clone_cols], chains = ['alpha','beta'], organism = 'human')
df_net, G = tcr_ab_to_network(tr = tr, edge_threshold = 120)
sample_hla_dict, sample_hla_dict_4d, sample_hla_dict_2d = get_hla_dictionaries(tsv_file = '/fh/fast/gilbert_p/fg_data/koelle_covid/hla_file.tsv')
sample_hla_dict.keys()
# Confirm that all person_x are present in the hla_dictionary, if not the next step will fail
assert pd.Series(tr.clone_df.person_x.unique().tolist()).isin(sample_hla_dict.keys()).all()
tr = compute_hla_feasibility(tr, G, sample_hla_dict, sample_hla_dict_4d, sample_hla_dict_2d)
tr.clone_df['feasible_hla_4d_i'] == ["A*03:01"]
ix  = tr.clone_df['feasible_hla_4d_i'].apply(lambda x : x == ['A*03:01']) 
ix2 = tr.clone_df['cell_type'] == "CD8"
tr.clone_df[ix&ix2]
tr.clone_df.to_csv('/fh/fast/gilbert_p/fg_data/koelle_covid/project_code/hybrid_t_cell_response/phenotypes/all_E03_phenotypes_hla_allele.tsv', sep = "\t", index = False)
"""
