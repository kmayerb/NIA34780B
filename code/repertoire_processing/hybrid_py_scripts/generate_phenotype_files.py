"""
This script combines information from Cell Ranger Output files and matching bulk TCRb 
repertoires to (i) assign CD4 and CD8 phenotypes to each clone and (ii) measure the 
degree to which the AIM+ clone is enriched by AIM sorting relative to its frequency
in the unenriched bulk repertoire.

The primary output is a set of phenotype_unique files that are used in 
later steps of the analysis to track AIM+ S-reactive clones abundance 
through a time series. 

The complexity of the analysis task is managed by a json file that contains 
groupings of relevant files
"""
from hybrid.assign_cd4_cd8 import get_information_from_10X
import os
import json

with open('json/2022_05_18_sc_samples.json', "r") as jsonfile:
    instructions = json.load(jsonfile)

for sample_name in instructions.keys():
    filepaths = instructions[sample_name]
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

"""
Phenotype unique columns
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
"""