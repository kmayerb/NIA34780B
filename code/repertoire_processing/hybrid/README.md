# hybrid_t_cell_response

Code supporting manuscript: 
**Memory CD8 T cells from prior SARS-CoV-2 infection dominate cellular immunity after mRNA vaccination**

## Contents

#### Functions

The `/hybrid/` folder contains functions used throughout the project are contained in 

#### Python Scripts

The `/project_py_scripts/` folder contains .py scripts resusing functions in the `/hybrid/` folder. 

#### JSON files

the `/json/` folder contains JSON files storing paths to  particiapnt related longitudinal files. 

#### Phenotype files

The `/phenotypes/` folder contains tables of information aggregated for identical clonotypes 
found in AIM+ cells. These tabels are generated with the (`project_py_sripts/generate_phenotype_files.py`).

#### Plotting Data

The `/plot_data/` folder contains intermediate outputs of the project's python scripts which are 
then used for vizualizaiton and statistical analysis.

#### Plotting R scripts

The `/project_R_scripts/` contain .R scripts used primarily for visualization and statistical testing

### Network Plot

The `/notebooks/` folder contains project notebook used for running tcrdist3 and constructing 
a sequence similarity network and background adjusted motifs.

## Usage snippets

### Example Parsing Adaptive Files 

Parsed dataframe contains only In-Frame CDR3 with valid IMGT TRBV and TRBJ gene calls and
pseudogenes (e.g., TCRB22-01) are also not considered. The results of this 
parsing function are sensitive to using the latest version of the tcrdist3, 
which has the most up to date mapping of adaptive gene names to IMGT names.

```python
from hybrid.parse import parse_adaptive_v2
import pandas as pd
f = '/fh/fast/gilbert_p/fg_data/koelle_covid/covid_tcr_compression/15761_4_TCRB.tsv'
df_adaptive = pd.read_csv(f, sep = "\t")
df_parsed   = parse_adaptive_v2(df_adaptive) 
```

### Example Assemble Longitudinal Data

A key feature of the study design was longitudinal repertoire sequencing after repeated antigen exposures.
Thus it is useful to track the abundance of each clone through the entire timecourse. Sample names and file
paths should be specified in a .json file (see examples in the `/json/` folder).



```python
import json
import os 
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe

ptid = '15673'
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

# optional to compute expansion statistics
dftextra = compute_fold_changes_on_longitudinal_dataframe(
    dft = dft.copy(), 
    dest = ".",
    ptid = None,
    write = False)
```

### Get Phenotype Information

For some samples we have phenotype information after AIM sorting and 
single cell sequencing samples with DNA-barcoded antibodies. 
In some case the TRBV-CDR3B-TRJV maps to multiple unique TCRab clones.
Thus both a "unique" output - retaining the match  most abundant TCRab clone -
and a "all" output -- retaining all matches -- are produced.

### Derive bulk phenotypes from a JSON file for 15 Files (May 23, 2022)

```python
from hybrid.assign_cd4_cd8 import get_information_from_10X
import os
import json

with open('json/2022_05_18_sc_samples.json', "r") as jsonfile:
    instructions = json.load(jsonfile)

for sample_name in instructions.keys():
    filepaths = prejson[sample_name]
    ptid = sample_name[0:5]
    bulk_visit = sample_name[5:]
    print(ptid, visit, filepaths)
    # Get phenotypes 
    phenotype_unique, phenotype_all = get_information_from_10X(
        clones_file_10X  = filepaths["filtered_contig_annotations.csv"]['filepath'],
        matrix_file_10X  = filepaths['raw_feature_bc_matrix.h5']['filepath'],
        clones_file_adpt = filepaths['clones_file_adpt']['filepath'])
    print(phenotype_unique)
    print(sample_name)
    bulk_filename = filepaths['clones_file_adpt']['filename'].replace(".tsv","")
    print(bulk_filename)
    phenotype_unique.to_csv(f'phenotypes/{bulk_filename}_{bulk_visit}.phenotype_unique.tsv', sep = "\t", index = False)
    phenotype_all.to_csv(f'phenotypes/{bulk_filename}_{bulk_visit}.phenotype_all.tsv', sep = "\t", index = False)
```


### What are the columns in the phenotype outputs

* `cdr3_a_nt` -- 10X TCRa nucleotide seqeunce
* `v_a_gene` -- 10X TCRa TRAV gene
* `j_a_gene` -- 10X TCRa TRAJ gene
* `cdr3_a_aa` --  10X TCRa CDR3 amino acid sequence
* `cdr3_b_nt`  -- 10X TCRb CDR3 nucleotide sequence
* `v_b_gene_x` -- 10X TCRb TRBV gene
* `j_b_gene_x` -- 10X TCRb TRBJ gene
* `cdr3_b_aa_x` -- 10X TCRb CDR3 amino sequence
* `total_counts_cd8` -- UMI counts associated with anti-CD8 antibodies
* `pct_counts_cd8` -- % UMI counts associated with anti-CD8 antibodies
* `total_counts_cd4` -- UMI counts associated with anti-CD4 antibodies
* `pct_counts_cd4` -- % UMI counts associated with anti-CD4 antibodies
* `cd8vcd4_score` -- np.log(r['pct_counts_cd8']+1)-np.log(r['pct_counts_cd4']+1)
* `barcode` -- NUMBER of droplet barcodes assigned to this TCRab clonotype
* `cell_type` -- CD8 if SCORE >=1, CD4 if SCORE <=1 , where SCORE = np.log(r['pct_counts_cd8']+1)-np.log(r['pct_counts_cd4']+1) 
* `key` -- TRBV-CDR3B-TRJV key for joining TCRab to bulk TCRb clone
* `N_10x` -- Total number of cells in the 10X reaction after Cell Ranger Filtering
* `10x_count` -- Same as barcode
* `10x_clone_id` -- 10X unique clone id for tracking duplicate join one to many
* `cdr3_b_nucseq` -- Bulk TCRb CDR3 nucleotide sequence
* `v_b_gene_y` -- Bulk TCRb TRBV nucleotide sequence
* `j_b_gene_y` -- Bulk TCRb TRBJ nucleotide sequence
* `cdr3_b_aa_y` -- Bulk TCRb CDR3 amino acid sequence
* `templates` -- Number of unique moelecular templates associated with that bulk clone
* `frequency` -- Overall frequency
* `productive_frequency` -- productive frequency
* `adpt_clone_id` -- Bulk 10X unique clone id for tracking duplicate join one to many
* `nucleotide_level_match` -- True if 10X CDR3b a subset of CDRb nucleotide sequence in bulk
* `10x_pfreq` -- Esimated productive frequency in the 10X AIM stimulated cell subset
* `10x_enriched` -- True, if 10x_pfreq > productive_frequency
* `10x_binomial_cdf` -- Binomial 1-CDF of observing  10x_count, conditional on productive frequency in bulk
* `10x_binomial_pmf` -- Binomial Probability Mass Function of observing 10x_count, conditional on productive frequency in bulk
* `10x_cdf_p_fdr_10x` -- FDR adjustment of p-value Binomial 1-CDF
* `10x_pmf_p_fdr_10x` -- FDR adjustment of p-value (PMF)
* `clones_file_10X` --  Reference to input clones 10X file
* `matrix_file_10X` -- Reference to input feature matrix 10X file 
* `clones_file_adpt` -- Reference to input adaptive file






### Consider HLA feasibility of S-specific 

```python
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
tr = compute_hla_feasibility(tr, G, sample_hla_dict, sample_hla_dict_2d, sample_hla_dict_4d)
tr.clone_df['feasible_hla_4d_i'] == ["A*03:01"]
ix  = tr.clone_df['feasible_hla_4d_i'].apply(lambda x : x == ['A*03:01']) 
ix2 = tr.clone_df['cell_type'] == "CD8"
tr.clone_df[ix&ix2]
tr.clone_df.to_csv('/fh/fast/gilbert_p/fg_data/koelle_covid/project_code/hybrid_t_cell_response/phenotypes/all_E03_phenotypes_hla_allele.tsv', sep = "\t", index = False)
```