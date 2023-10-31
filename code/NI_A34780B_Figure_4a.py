# Figure_4.py
# Code for making sequence similarty network plots
import sys
import os 
import pandas as pd
from tcrdist.repertoire import TCRrep
from hybrid.feasible import tcr_ab_to_network
from hybrid.feasible import get_hla_dictionaries
from hybrid.feasible import compute_hla_feasibility

# Step 1: Load data and remove TCRs that were not enriched over bulk repertoire frequencies
data_path =  './data/'
# Load a compiled file containing all phenotyped cells
clones_file = pd.read_csv(os.path.join(data_path, '2022_06_23_all_AIM_phenotypes.tsv'), sep = "\t")
clones_file['person_x'] = clones_file['source_file'].str.split('_').apply(lambda x: int(x[0]))
# Remove TCRs that were not enriched over bulk repertoire frequencies¶
    # For network clustering
    # matched bulk. We use this comparatively less stingent method since
    # we are detecting groups of pubic AIM+ similar TCRab ssquences. More stringent
    # statistical criteria for AIM+ enrichment is not used for graph based analaysis.
clones_file = clones_file[clones_file['10x_enriched'] == True].reset_index(drop = True)
# Format Data.Frame for use in tcrdist3
df  = clones_file.rename(columns = {'cdr3_b_aa_x':'cdr3_b_aa','v_b_gene_x':'v_b_gene','j_b_gene_x':'j_b_gene','hla_key':'person_x','barcode':'count'})
df['v_a_gene'] = df['v_a_gene'].apply(lambda x : f"{x}*01")
df['j_a_gene'] = df['j_a_gene'].apply(lambda x : f"{x}*01")
df['v_b_gene'] = df['v_b_gene'].apply(lambda x : f"{x}*01")
df['j_b_gene'] = df['j_b_gene'].apply(lambda x : f"{x}*01")
df['cell_type'] = df['cell_type'].fillna('UNK')

#Step 2: Compute HLA feasibility
# We try to infer feasibly HLA by walking the TCR sequence similarity network to 
# identify HLA-alleles that common among participants sharing graph edges 
# formed between connected components of similar sequence. 
# Load a TCRrep insance:
clone_cols = ['person_x', 'cdr3_a_aa', 'cdr3_b_aa','v_a_gene','v_b_gene','j_a_gene','j_b_gene','cell_type', 'total_counts_cd8', 'total_counts_cd4', 'count']
tr = TCRrep(cell_df = df[clone_cols], chains = ['alpha','beta'], organism = 'human')
df_net, G = tcr_ab_to_network(tr = tr, edge_threshold = 100)
sample_hla_dict, sample_hla_dict_4d, sample_hla_dict_2d = get_hla_dictionaries(tsv_file = os.path.join(data_path,'hla_file.tsv'))
sample_hla_dict.keys()
# Confirm that all person_x are present in the hla_dictionary, if not the next step will fail
assert pd.Series(tr.clone_df.person_x.unique().tolist()).isin(sample_hla_dict.keys()).all()
# Run hla feasibility routine, to find the interesction of feasible HLA within a clones neighborhood, walking the connected components
tr = compute_hla_feasibility(tr, G, sample_hla_dict, sample_hla_dict_2d, sample_hla_dict_4d)

# Step 3: Add meta-data on expansion, cell_type, and possible HLA restriction
cols =['v_a_gene', 'cdr3_a_aa', 'j_a_gene', 'v_b_gene', 'cdr3_b_aa', 'j_b_gene']
tr.clone_df['key_ab'] = tr.clone_df.apply(lambda r : "+".join(r[cols].to_list()), axis = 1)
# We take a moment to consider the number of S-reactive clones that matched a 
# unambiguosly expanded clone
import numpy as np
df_expand = pd.read_csv(os.path.join(data_path, '2022_06_23_df_expand.tsv'),sep = "\t") 
df_expand['v_a_gene'] = df_expand['v_a_gene'].apply(lambda x : f"{x}*01")
df_expand['j_a_gene'] = df_expand['j_a_gene'].apply(lambda x : f"{x}*01")
                      
# The unique features that define a key
complete_cases = df_expand[cols].apply(lambda r: np.all(r.notna()), axis = 1)
df_expand = df_expand[complete_cases].reset_index(drop = True)
df_expand['key_ab'] = df_expand.apply(lambda r : "+".join(r[cols].to_list()), axis = 1)

df_contract = pd.read_csv(os.path.join(data_path, '2022_06_23_df_contract.tsv'), sep = "\t") 
df_contract ['v_a_gene'] = df_contract['v_a_gene'].apply(lambda x : f"{x}*01")
df_contract ['j_a_gene'] = df_contract['j_a_gene'].apply(lambda x : f"{x}*01")
# The unique features that define a key
complete_cases = df_contract[cols].apply(lambda r: np.all(r.notna()), axis = 1)
df_contract = df_contract[complete_cases].reset_index(drop = True)
df_contract['key_ab'] = df_contract.apply(lambda r : "+".join(r[cols].to_list()), axis = 1)
#df_expand.merge(tr.clone_df, how = "inner", on = ["cdr3_b_aa",'c#dr3_a_aa','v_a_gene','j_a_gene','v_b_gene','j_b_gene'])
tr.clone_df['expanded'] = tr.clone_df['key_ab'].isin(df_expand['key_ab'].to_list())
tr.clone_df['contracted'] = tr.clone_df['key_ab'].isin(df_contract['key_ab'].to_list())

# CD4/CD8 tags
# assign CD4 and CD8 cell type based on median CD4 vs CD4 umi counts
import numpy as np
tr.clone_df['score'] = np.log(tr.clone_df['total_counts_cd8']) - np.log(tr.clone_df['total_counts_cd4']) 
tr.clone_df['total_cd_counts'] = tr.clone_df['total_counts_cd8'] + tr.clone_df['total_counts_cd4']
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
tr.clone_df['cell_type_2'] = tr.clone_df.apply(lambda r: cd4cd8_score2(score = r['score'], total_cd_counts = r['total_cd_counts']), axis = 1)
# replace old cell_type with this definition
tr.clone_df['cell_type'] = tr.clone_df['cell_type_2']
# HLA-association tags
# <tags> is a list of all the inferred common HLA groups
tags = list()
for i,r in tr.clone_df.iterrows():
    if len(r['feasible_hla_2d_i']) in [1] and r['cell_type'] == "CD8" :
        tags.append("_".join(r['feasible_hla_2d_i']))
    elif len(r['feasible_hla_2d_ii']) in  [1,2] and r['cell_type'] == "CD4" :
        tags.append("_".join(r['feasible_hla_2d_ii']))
    else: 
        tags.append("U")            
tr.clone_df['restriction_tag'] = tags    


# Step 4 Sequence Similarity Network
from tcrdist.public import _neighbors_fixed_radius
import pandas as pd
import os
# <edge_threshold> is used to define maximum distance to for a network edge.
edge_threshold = 100
# <tr.pw_alpha_beta> is paired chain TCRdist.
tr.pw_alpha_beta = tr.pw_beta + tr.pw_alpha
# <network> initialize a list to populate with edges between TCRs.
network = list()
for i,n in enumerate(_neighbors_fixed_radius(tr.pw_alpha_beta, edge_threshold)):
  for j in n:
      if i != j:
          network.append((
              i,                                 # 'node_1' - row index
              j,                                 # 'node_2' - column index
              (tr.pw_alpha_beta )[i,j],          # 'dist'- gets the distance between TCR(i,j)
              tr.clone_df['v_b_gene'].iloc[i],   # 'v_b_gene_1' - v beta gene of clone i
              tr.clone_df['v_b_gene'].iloc[j],   # 'v_b_gene_2' - v beta gene of clone j
              tr.clone_df['cdr3_b_aa'].iloc[i],  # 'cdr3_b_aa_1' - cdr3 beta of clone i
              tr.clone_df['cdr3_b_aa'].iloc[j],  # 'cdr3_b_aa_2' - cdr3 beta of clone j
              tr.clone_df['v_a_gene'].iloc[i],   # 'v_b_gene_1' - v alpha gene of clone i
              tr.clone_df['v_a_gene'].iloc[j],   # 'v_b_gene_2' - v alpha gene of clone j
              tr.clone_df['cdr3_a_aa'].iloc[i],  # 'cdr3_a_aa_1' - cdr3 alpha of clone i  
              tr.clone_df['cdr3_a_aa'].iloc[j],  # 'cdr3_a_aa_2' - cdr3 alpha of clone j  
              tr.clone_df['count'].iloc[i],      #  count 
              tr.clone_df['count'].iloc[j],      #  count
              tr.clone_df['total_counts_cd8'].iloc[i],      #  count 
              tr.clone_df['total_counts_cd8'].iloc[j],      #  count
              tr.clone_df['total_counts_cd4'].iloc[i],      #  count 
              tr.clone_df['total_counts_cd4'].iloc[j],      #  count
              tr.clone_df['person_x'].iloc[i],      #  count 
              tr.clone_df['person_x'].iloc[j],      #  count
              tr.clone_df['restriction_tag'].iloc[i],
              tr.clone_df['restriction_tag'].iloc[j],
              tr.clone_df['cell_type'].iloc[i],
              tr.clone_df['cell_type'].iloc[j],
              tr.clone_df['expanded'].iloc[i],
              tr.clone_df['expanded'].iloc[j],
              tr.clone_df['contracted'].iloc[i],
              tr.clone_df['contracted'].iloc[j],
              len(n)-1))                         # 'K_neighbors' - number of neighbors

cols = ['node_1', 'node_2', 'dist', 
      'v_b_gene_1', 'v_b_gene_2', 
      'cdr3_b_aa_1','cdr3_b_aa_2', 
      'v_a_gene_1', 'v_a_gene_2', 
      'cdr3_a_aa_1','cdr3_a_aa_2',
      'count_1', 'count_2',
      'total_counts_cd8_1', 'total_counts_cd8_2',
      'total_counts_cd4_1', 'total_counts_cd4_2',
      'subject_1', 'subject_2',
      'restriction_tag_1', 'restriction_tag_2',
      'cell_type_1','cell_type_2',
       'expanded_1','expanded_2',
       'contracted_1','contracted_2',
      'K_neighbors']
# Store the <network> edge list as a DataFrame.        
df_net = pd.DataFrame(network, columns = cols)
df_net['weight'] = (edge_threshold - df_net['dist'])/edge_threshold


# Graphically display of sequence similarity network
# This graph shows public (orange) and private edges (gray).
# Nodes are colored blue for CD8 and green for CD4.
# Cluster numbers are shown in purple which are removed later using Adobe Illustrator

################
### CD4/CD8 ####
################
import networkx as nx
from community import community_louvain
from tcrdist.html_colors import get_html_colors
from matplotlib import pyplot as plt
# <G> Initialize a networkx Graph instance from the columns of df_net.
G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],
                                          'target' : df_net['node_2'],
                                          'weight' : df_net['weight'],
                                          'count1' : df_net['count_1'],
                                          'count2' : df_net['count_2']}))

# get cluster definitions
G2 = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],'target' : df_net['node_2'], 'weight' :df_net['weight']}))
partition= community_louvain.best_partition(G2)
# Change partition such that cluster Id is in descending order based on community size 
partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size, 
                                                range(len(partitions_by_cluster_size)))}
partition = {k:partition_reorder.get(v) for k,v in partition.items()}

print(G)
print(nx.number_connected_components(G), "connected components")
plt.figure(1, figsize=(15, 15))
# layout graphs with positions using graphviz neato
pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
# color nodes the same in each connected subgraph
C = (G.subgraph(c) for c in nx.connected_components(G))

# extract nodes, counts, subjsect, cell_types, restrictions, expansion_status
ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()
sbs = df_net['subject_1'].to_list() + df_net['subject_2'].to_list()
tps = df_net['cell_type_1'].to_list() + df_net['cell_type_2'].to_list()
tgs = df_net['restriction_tag_1'].to_list() + df_net['restriction_tag_2'].to_list()
exs = df_net['expanded_1'].to_list() + df_net['expanded_2'].to_list()
# Dictionaries mapping nodes to counts, subjsect, cell_types, restrictions, expansion_status
node_cts = {node: x for node, x in zip(ns,cts)}
node_sbs = {node: x for node, x in zip(ns,sbs)}
node_tps = {node: x for node, x in zip(ns,tps)}
node_tgs = {node: x for node, x in zip(ns,tgs)}
node_exs = {node: x for node, x in zip(ns,exs)}

# tags for cell_types
tags =pd.Series(tps).value_counts().index
colors   = get_html_colors(len(pd.Series(tags).value_counts().index))
tag_color = {'CD8':'blue', 'CD4':'green', 'LOW_CD':'gray', 'CD4CD8':'black'}

for g in C:
    # For labeling elements of cluster:
    out_dict =dict()
    count = 0
    v_so_far = dict()
    for k,v in {i : partition.get(i) for i in g.nodes}.items():
        if v <= 300 and v_so_far.get(v, True):
            out_dict[k] =v 
        else:
            out_dict[k] = ""
        v_so_far[v] = False
    
    public_edges = [node_sbs.get(u) != node_sbs.get(v) for u,v in g.edges]
    public_edge_color = [{True :'orange', False : 'gray'}.get(x) for x in public_edges]
    #c = [random.random()] * nx.number_of_nodes(g)  # random color...
    nx.draw(g, pos, alpha= .5, 
        node_size = [5*node_cts.get(x) for x in g.nodes], 
        #node_color = [{True : "red", False: "black"}.get(node_exs.get(x)) for x in g.nodes],
        edge_color = public_edge_color,
        #node_color = [tag_color.get(node_tgs.get(x)) for x in g.nodes],
        node_color = [tag_color.get(node_tps.get(x)) for x in g.nodes],
        vmin=0.0, 
        vmax=1.0, 
        labels = out_dict,
        with_labels=True,
        font_size = 20, font_color = "violet")

# If saving figure:
# plt.savefig("figures/network_blue_cd8_green_cd4_black_CD4CD8.pdf")

################
### Jurkats ####
################

import networkx as nx
from community import community_louvain
from tcrdist.html_colors import get_html_colors
from matplotlib import pyplot as plt
# <G> Initialize a networkx Graph instance from the columns of df_net.
G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],
                                          'target' : df_net['node_2'],
                                          'weight' : df_net['weight'],
                                          'count1' : df_net['count_1'],
                                          'count2' : df_net['count_2']}))

# get cluster definitions
G2 = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],'target' : df_net['node_2'], 'weight' :df_net['weight']}))
partition= community_louvain.best_partition(G2)
# Change partition such that cluster Id is in descending order based on community size 
partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size, 
                                                range(len(partitions_by_cluster_size)))}
partition = {k:partition_reorder.get(v) for k,v in partition.items()}

print(G)
print(nx.number_connected_components(G), "connected components")
plt.figure(1, figsize=(15, 15))
# layout graphs with positions using graphviz neato
pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
# color nodes the same in each connected subgraph
C = (G.subgraph(c) for c in nx.connected_components(G))

ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()
sbs = df_net['subject_1'].to_list() + df_net['subject_2'].to_list()
tps = df_net['cell_type_1'].to_list() + df_net['cell_type_2'].to_list()
tgs = df_net['restriction_tag_1'].to_list() + df_net['restriction_tag_2'].to_list()
exs = df_net['expanded_1'].to_list() + df_net['expanded_2'].to_list()

node_cts = {node: x for node, x in zip(ns,cts)}
node_sbs = {node: x for node, x in zip(ns,sbs)}
node_tps = {node: x for node, x in zip(ns,tps)}
node_tgs = {node: x for node, x in zip(ns,tgs)}
node_exs = {node: x for node, x in zip(ns,exs)}


node_color = dict()
for n in ns:
    # Validate jurkat clonse
    if n in [2382, 2444, 2748, 2836, 2859, 2959, 3018]:
        node_color[n] = "red"
    else:
        node_color[n] = "black"


for g in C:
    
    # For labeling elements of cluster:
    out_dict =dict()
    count = 0
    v_so_far = dict()
    for k,v in {i : partition.get(i) for i in g.nodes}.items():
        if v <= 300 and v_so_far.get(v, True):
            out_dict[k] =v 
        else:
            out_dict[k] = ""
        v_so_far[v] = False
    
    public_edges = [node_sbs.get(u) != node_sbs.get(v) for u,v in g.edges]
    public_edge_color = [{True :'orange', False : 'gray'}.get(x) for x in public_edges]
    #c = [random.random()] * nx.number_of_nodes(g)  # random color...
    nx.draw(g, pos, alpha= .5, 
        node_size = [5*node_cts.get(x) for x in g.nodes], 
        edge_color = public_edge_color,
        node_color = [node_color.get(x) for x in g.nodes],
        vmin=0.0, 
        vmax=1.0, 
        labels = out_dict,
        with_labels=True,
        font_size = 20, font_color = "violet")

#plt.savefig("figures/network_blue_cloned_red_other_black.pdf")


######################
### HLA Inference ####
######################
# Assign colors to each HLA tags
tags =pd.Series(tgs).value_counts().index
colors   = get_html_colors(len(pd.Series(tags).value_counts().index))
tag_color = {t:c for t,c in zip(tags, colors)}
tag_color['U']= "gray"

cat_colors = ["#e41a1c","#377eb8","#4daf4a","#984ea3",
"#ff7f00","#ffff33","#a65628","#f781bf",
"#66c2a5","#fc8d62","#8da0cb","#e78ac3",
"#a6d854","#ffd92f","#e5c494","#b3b3b3",
"#8dd3c7","#ffffb3","#bebada","#fb8072",
"#80b1d3","#fdb462","#b3de69","#fccde5",
"#d9d9d9","#bc80bd","#ccebc5","#ffed6f"]
i = 0
tag_color_2 = dict()
for k,v in tag_color.items():
    tag_color_2[k] = cat_colors[i]
    i += 1
tag_color_2['U'] = 'black'


label_dict = {g: 'b' for g in G.nodes}
print(G)
print(nx.number_connected_components(G), "connected components")
plt.figure(1, figsize=(15, 15))
# layout graphs with positions using graphviz neato
pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
# color nodes the same in each connected subgraph
C = (G.subgraph(c) for c in nx.connected_components(G))

ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()
sbs = df_net['subject_1'].to_list() + df_net['subject_2'].to_list()
tps = df_net['cell_type_1'].to_list() + df_net['cell_type_2'].to_list()
tgs = df_net['restriction_tag_1'].to_list() + df_net['restriction_tag_2'].to_list()

tags =pd.Series(tgs).value_counts().index
node_cts = {node: x for node, x in zip(ns,cts)}
node_sbs = {node: x for node, x in zip(ns,sbs)}
node_tps = {node: x for node, x in zip(ns,tps)}
node_tgs = {node: x for node, x in zip(ns,tgs)}
colors   = get_html_colors(len(pd.Series(tags).value_counts().index))
tag_color = {t:c for t,c in zip(tags, colors)}
tag_color['U']= "gray"
for g in C:
    public_edges = [node_sbs.get(u) != node_sbs.get(v) for u,v in g.edges]
    public_edge_color = [{True :'orange', False : 'black'}.get(x) for x in public_edges]
    
    #c = [random.random()] * nx.number_of_nodes(g)  # random color...
    nx.draw(g, pos, alpha= .2, 
    node_size = [10*node_cts.get(x) for x in g.nodes], 
    node_color = [tag_color_2.get(node_tgs.get(x)) for x in g.nodes],
    #edge_color = public_edge_color,
    labels = {k:node_tgs.get(k) for k in g.nodes}, 
    vmin=0.0, 
    vmax=1.0, 
    with_labels=True,
    font_size = 7, font_color = "violet")

#plt.savefig('figures/hla_walk_labeled_network_16_x_16.pdf')

# HLA legend
from plotnine import *
legend_df = pd.DataFrame({'label':tag_color_2.keys(), 'color' : tag_color_2.values()} )
legend_df['x'] = legend_df.index

(ggplot(legend_df, aes(x = 1, y = 'label', fill = 'label'))
+ geom_point(size = 5, alpha = .5) 
+ scale_fill_manual(values = tag_color_2) 
+ theme_classic()
).save('figures/hla_legend_colors.pdf')

(ggplot(legend_df, aes(x = 1, y = 'label', fill = 'label'))
+ geom_point(size = 5, alpha = .5) 
+ scale_fill_manual(values = tag_color_2) 
+ theme_classic()
)

######################
### Cluster Number####
######################
import os 
import pandas as pd
import networkx as nx
import community.community_louvain as community_louvain
from tcrdist.repertoire import TCRrep
import networkx as nx

G = nx.from_pandas_edgelist(pd.DataFrame({'source' : df_net['node_1'],'target' : df_net['node_2'], 'weight' :df_net['weight']}))
partition= community_louvain.best_partition(G)
# Change partition such that cluster Id is in descending order based on community size 
partitions_by_cluster_size = list(pd.Series(partition.values()).value_counts().index)
partition_reorder = {id:rank for id,rank in zip(partitions_by_cluster_size, 
                                                range(len(partitions_by_cluster_size)))}
partition = {k:partition_reorder.get(v) for k,v in partition.items()}

from tcrdist.html_colors import get_html_colors
clusters = [i for i in pd.Series(partition.values()).value_counts().index]
colors   = get_html_colors(len(clusters))
cluster_to_color = {cluster:color for cluster,color, in zip(clusters,colors)}
print(G)
print(nx.number_connected_components(G), "connected components")

plt.figure(1, figsize=(16, 16))
# layout graphs with positions using graphviz neato
pos = nx.nx_agraph.graphviz_layout(G, prog="neato")
# color nodes the same in each connected subgraph
C = (G.subgraph(c) for c in nx.connected_components(G))
ns = df_net['node_1'].to_list() + df_net['node_2'].to_list()
cts = df_net['count_1'].to_list() + df_net['count_2'].to_list()
node_cts = {node: x for node, x in zip(ns,cts)}
for g in C:
    out_dict =dict()
    count = 0
    v_so_far = dict()
    for k,v in {i : partition.get(i) for i in g.nodes}.items():
        if v <= 100 and v_so_far.get(v, True):
            out_dict[k] =v 
        else:
            out_dict[k] = ""
        v_so_far[v] = False
    #c = [random.random()] * nx.number_of_nodes(g)  # random color...
    nx.draw(g, pos, alpha= .5, node_size = [5*node_cts.get(x) for x in g.nodes], 
            node_color = [cluster_to_color.get(partition.get(i)) for i in g.nodes],
            edge_color = "gray",
            labels = out_dict,
            vmin=0.0, vmax=1.0, with_labels=True, font_size = 20, font_color = "violet")
#plt.savefig("figures/network_cluster_number.pdf")

# legend
convert_html = ['IndianRed	#CD5C5C	rgb(205, 92, 92)','LightCoral	#F08080	rgb(240, 128, 128)','Salmon	#FA8072	rgb(250, 128, 114)','DarkSalmon	#E9967A	rgb(233, 150, 122)','LightSalmon	#FFA07A	rgb(255, 160, 122)','Crimson	#DC143C	rgb(220, 20, 60)','Red	#FF0000	rgb(255, 0, 0)','FireBrick	#B22222	rgb(178, 34, 34)','DarkRed	#8B0000	rgb(139, 0, 0)','Pink	#FFC0CB	rgb(255, 192, 203)','LightPink	#FFB6C1	rgb(255, 182, 193)','HotPink	#FF69B4	rgb(255, 105, 180)','DeepPink	#FF1493	rgb(255, 20, 147)','MediumVioletRed	#C71585	rgb(199, 21, 133)','PaleVioletRed	#DB7093	rgb(219, 112, 147)','LightSalmon	#FFA07A	rgb(255, 160, 122)','Coral	#FF7F50	rgb(255, 127, 80)','Tomato	#FF6347	rgb(255, 99, 71)','OrangeRed	#FF4500	rgb(255, 69, 0)','DarkOrange	#FF8C00	rgb(255, 140, 0)','Orange	#FFA500	rgb(255, 165, 0)','Gold	#FFD700	rgb(255, 215, 0)','Yellow	#FFFF00	rgb(255, 255, 0)','LightYellow	#FFFFE0	rgb(255, 255, 224)','LemonChiffon	#FFFACD	rgb(255, 250, 205)','LightGoldenrodYellow	#FAFAD2	rgb(250, 250, 210)','PapayaWhip	#FFEFD5	rgb(255, 239, 213)','Moccasin	#FFE4B5	rgb(255, 228, 181)','PeachPuff	#FFDAB9	rgb(255, 218, 185)','PaleGoldenrod	#EEE8AA	rgb(238, 232, 170)','Khaki	#F0E68C	rgb(240, 230, 140)','DarkKhaki	#BDB76B	rgb(189, 183, 107)','Lavender	#E6E6FA	rgb(230, 230, 250)','Thistle	#D8BFD8	rgb(216, 191, 216)','Plum	#DDA0DD	rgb(221, 160, 221)','Violet	#EE82EE	rgb(238, 130, 238)','Orchid	#DA70D6	rgb(218, 112, 214)','Fuchsia	#FF00FF	rgb(255, 0, 255)','Magenta	#FF00FF	rgb(255, 0, 255)','MediumOrchid	#BA55D3	rgb(186, 85, 211)','MediumPurple	#9370DB	rgb(147, 112, 219)','RebeccaPurple	#663399	rgb(102, 51, 153)','BlueViolet	#8A2BE2	rgb(138, 43, 226)','DarkViolet	#9400D3	rgb(148, 0, 211)','DarkOrchid	#9932CC	rgb(153, 50, 204)','DarkMagenta	#8B008B	rgb(139, 0, 139)','Purple	#800080	rgb(128, 0, 128)','Indigo	#4B0082	rgb(75, 0, 130)','SlateBlue	#6A5ACD	rgb(106, 90, 205)','DarkSlateBlue	#483D8B	rgb(72, 61, 139)','MediumSlateBlue	#7B68EE	rgb(123, 104, 238)','GreenYellow	#ADFF2F	rgb(173, 255, 47)','Chartreuse	#7FFF00	rgb(127, 255, 0)','LawnGreen	#7CFC00	rgb(124, 252, 0)','Lime	#00FF00	rgb(0, 255, 0)','LimeGreen	#32CD32	rgb(50, 205, 50)','PaleGreen	#98FB98	rgb(152, 251, 152)','LightGreen	#90EE90	rgb(144, 238, 144)','MediumSpringGreen	#00FA9A	rgb(0, 250, 154)','SpringGreen	#00FF7F	rgb(0, 255, 127)','MediumSeaGreen	#3CB371	rgb(60, 179, 113)','SeaGreen	#2E8B57	rgb(46, 139, 87)','ForestGreen	#228B22	rgb(34, 139, 34)','Green	#008000	rgb(0, 128, 0)','DarkGreen	#006400	rgb(0, 100, 0)','YellowGreen	#9ACD32	rgb(154, 205, 50)','OliveDrab	#6B8E23	rgb(107, 142, 35)','Olive	#808000	rgb(128, 128, 0)','DarkOliveGreen	#556B2F	rgb(85, 107, 47)','MediumAquamarine	#66CDAA	rgb(102, 205, 170)','DarkSeaGreen	#8FBC8B	rgb(143, 188, 139)','LightSeaGreen	#20B2AA	rgb(32, 178, 170)','DarkCyan	#008B8B	rgb(0, 139, 139)','Teal	#008080	rgb(0, 128, 128)','Aqua	#00FFFF	rgb(0, 255, 255)','Cyan	#00FFFF	rgb(0, 255, 255)','LightCyan	#E0FFFF	rgb(224, 255, 255)','PaleTurquoise	#AFEEEE	rgb(175, 238, 238)','Aquamarine	#7FFFD4	rgb(127, 255, 212)','Turquoise	#40E0D0	rgb(64, 224, 208)','MediumTurquoise	#48D1CC	rgb(72, 209, 204)','DarkTurquoise	#00CED1	rgb(0, 206, 209)','CadetBlue	#5F9EA0	rgb(95, 158, 160)','SteelBlue	#4682B4	rgb(70, 130, 180)','LightSteelBlue	#B0C4DE	rgb(176, 196, 222)','PowderBlue	#B0E0E6	rgb(176, 224, 230)','LightBlue	#ADD8E6	rgb(173, 216, 230)','SkyBlue	#87CEEB	rgb(135, 206, 235)','LightSkyBlue	#87CEFA	rgb(135, 206, 250)','DeepSkyBlue	#00BFFF	rgb(0, 191, 255)','DodgerBlue	#1E90FF	rgb(30, 144, 255)','CornflowerBlue	#6495ED	rgb(100, 149, 237)','MediumSlateBlue	#7B68EE	rgb(123, 104, 238)','RoyalBlue	#4169E1	rgb(65, 105, 225)','Blue	#0000FF	rgb(0, 0, 255)','MediumBlue	#0000CD	rgb(0, 0, 205)','DarkBlue	#00008B	rgb(0, 0, 139)','Navy	#000080	rgb(0, 0, 128)','MidnightBlue	#191970	rgb(25, 25, 112)','Cornsilk	#FFF8DC	rgb(255, 248, 220)','BlanchedAlmond	#FFEBCD	rgb(255, 235, 205)','Bisque	#FFE4C4	rgb(255, 228, 196)','NavajoWhite	#FFDEAD	rgb(255, 222, 173)','Wheat	#F5DEB3	rgb(245, 222, 179)','BurlyWood	#DEB887	rgb(222, 184, 135)','Tan	#D2B48C	rgb(210, 180, 140)','RosyBrown	#BC8F8F	rgb(188, 143, 143)','SandyBrown	#F4A460	rgb(244, 164, 96)','Goldenrod	#DAA520	rgb(218, 165, 32)','DarkGoldenrod	#B8860B	rgb(184, 134, 11)','Peru	#CD853F	rgb(205, 133, 63)','Chocolate	#D2691E	rgb(210, 105, 30)','SaddleBrown	#8B4513	rgb(139, 69, 19)','Sienna	#A0522D	rgb(160, 82, 45)','Brown	#A52A2A	rgb(165, 42, 42)','Maroon	#800000	rgb(128, 0, 0)','White	#FFFFFF	rgb(255, 255, 255)','Snow	#FFFAFA	rgb(255, 250, 250)','HoneyDew	#F0FFF0	rgb(240, 255, 240)','MintCream	#F5FFFA	rgb(245, 255, 250)','Azure	#F0FFFF	rgb(240, 255, 255)','AliceBlue	#F0F8FF	rgb(240, 248, 255)','GhostWhite	#F8F8FF	rgb(248, 248, 255)','WhiteSmoke	#F5F5F5	rgb(245, 245, 245)','SeaShell	#FFF5EE	rgb(255, 245, 238)','Beige	#F5F5DC	rgb(245, 245, 220)','OldLace	#FDF5E6	rgb(253, 245, 230)','FloralWhite	#FFFAF0	rgb(255, 250, 240)','Ivory	#FFFFF0	rgb(255, 255, 240)','AntiqueWhite	#FAEBD7	rgb(250, 235, 215)','Linen	#FAF0E6	rgb(250, 240, 230)','LavenderBlush	#FFF0F5	rgb(255, 240, 245)','MistyRose	#FFE4E1	rgb(255, 228, 225)','Gainsboro	#DCDCDC	rgb(220, 220, 220)','LightGray	#D3D3D3	rgb(211, 211, 211)','Silver	#C0C0C0	rgb(192, 192, 192)','DarkGray	#A9A9A9	rgb(169, 169, 169)','Gray	#808080	rgb(128, 128, 128)','DimGray	#696969	rgb(105, 105, 105)','LightSlateGray	#778899	rgb(119, 136, 153)','SlateGray	#708090	rgb(112, 128, 144)','DarkSlateGray	#2F4F4F	rgb(47, 79, 79)','Black	#000000	rgb(0, 0, 0)']
html_name_to_hex = dict()
for i in convert_html:
    html_name, hexname, rbg = (i.split("\t"))
    html_name_to_hex[html_name.lower()]=hexname
cluster_to_color_hex = {str(k):html_name_to_hex.get(v) for k,v in cluster_to_color.items()}

from plotnine import *
legend_df = pd.DataFrame({'label':cluster_to_color_hex.keys(), 'color' : cluster_to_color_hex.values()} )
legend_df['x'] = legend_df.index
legend_df['label'] = legend_df['label'].apply(lambda x: str(x))
legend_df['label_x'] = pd.Series(legend_df['label'].to_list(), dtype="category").cat.reorder_categories(legend_df['label'])
(ggplot(legend_df, aes(x = 1, y = 'x', fill = 'label_x'))
+ geom_point(size = 5, alpha = .5) 
+ scale_fill_manual(values = cluster_to_color_hex) 
+ theme_classic()
)

tr.clone_df['cluster_id'] = [str(partition.get(i)) if partition.get(i) is not None else None for i in tr.clone_df.index]
tr.clone_df['cluster_id_cat'] =  pd.Series(tr.clone_df['cluster_id'].to_list(), dtype="category").cat.reorder_categories(tr.clone_df['cluster_id'].value_counts().index)
#tr.clone_df.to_csv(os.path.join(data_path, '2022_08_14_all_E03_phenotypes_hla_allele_cluster_ids.tsv'), sep = "\t", index = False)
