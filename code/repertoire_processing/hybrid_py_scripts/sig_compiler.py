

# extract sig_expand clonotypes and merge with dfts to get measures of 10x capture
# second goal - grey dot graphs of the mystery cells 
# scheme - merge e03 file with dft file, then filter to just sig_exp clones from e03 

import os
import pandas as pd
import json

#10x info location
r1 = r"/fh/fast/corey_l/user/esford3/phenotypes"
#longitudinal file location 
r2 = r"/fh/fast/corey_l/user/esford3/output/dfts"

"""
#bot to help make the sample json
ptid_list = (15515,15525,15527,15531,
			15548,15577,15581,15630,
			15655,15669,15673,15684,
			15744,15754,15761,15782,15836)
for ptid in ptid_list:
	print(str(ptid)+"_10x': os.path.join(r1,'", "e03_10x_"+str(ptid)+"'),",
			str(ptid)+"_long': os.path.join(r2, '",(str(ptid)+".matrix_plus.tsv')"))

## notes about file sources (by file append)
phenotype_unique: files available online 
longitudinal_templates: via adaptive bulk data (online via immune access)
example code to generate longitudinal_templates is in hybrid/assemble_logitudinal_samples

matrix_plus: adaptive bulk data run w/ assemble_longitudinal_samples, followed by assignprevdetected 

"""
##or - 
import os
import json
import pandas as pd
path = r"/fh/fast/corey_l/user/esford3/phenotypes"
files = os.listdir(path)
files = sorted(files)
print(files)

#still needed a few formatting adjustments 
# problem - 5/17 have no matrix file (no E01 - need novel matrix file computations)
samples = {'15515':{'15515_10x': os.path.join(r1,'15515_4_TCRB_E03.phenotype_unique.tsv'), 
				'15515_long': os.path.join(r2, '15515.longitudinal_templates_only.tsv')},
			'15548':{'15548_10x': os.path.join(r1,'15548_7_TCRB_E03.phenotype_unique.tsv'), 
				'15548_long': os.path.join(r2, '15548.longitudinal_templates_only.tsv')},
			'15630':{'15630_10x': os.path.join(r1,'15630_3_TCRB_E03.phenotype_unique.tsv'),
				'15630_long': os.path.join(r2, '15630.longitudinal_templates_only.tsv')},
			'15655':{'15655_10x': os.path.join(r1,'15655_3_TCRB_E03.phenotype_unique.tsv'),
				'15655_long': os.path.join(r2, '15655.longitudinal_templates_only.tsv')},
			'15744':{'15744_10x': os.path.join(r1,'15744_3_TCRB_E03.phenotype_unique.tsv'), 
				'15744_long': os.path.join(r2, '15744.longitudinal_templates_only.tsv')}}

samples = {'15525':{'15525_10x': os.path.join(r1,'15525_5_TCRB_E03.phenotype_unique.tsv'), 
				'15525_long': os.path.join(r2, '15525.matrix_plus.tsv')},
			'15527':{'15527_10x': os.path.join(r1,'15527_6_TCRB_E03.phenotype_unique.tsv'), 
				'15527_long': os.path.join(r2, '15527.matrix_plus.tsv')},
			'15531':{'15531_10x': os.path.join(r1,'15531_6_TCRB_E03.phenotype_unique.tsv'), 
				'15531_long': os.path.join(r2, '15531.matrix_plus.tsv')},
			'15577':{'15577_10x': os.path.join(r1,'15577_4_TCRB_E03.phenotype_unique.tsv'), 
				'15577_long': os.path.join(r2, '15577.matrix_plus.tsv')},
			'15581':{'15581_10x': os.path.join(r1,'15581_9_TCRB_E03.phenotype_unique.tsv'), 
				'15581_long': os.path.join(r2, '15581.matrix_plus.tsv')},
			'15669':{'15669_10x': os.path.join(r1,'15669_4_TCRB_E03.phenotype_unique.tsv'), 
				'15669_long': os.path.join(r2, '15669.matrix_plus.tsv')},
			'15673':{'15673_10x': os.path.join(r1,'15673_8_TCRB_E03.phenotype_unique.tsv'),
				'15673_long': os.path.join(r2, '15673.matrix_plus.tsv')},
			'15684':{'15684_10x': os.path.join(r1,'15684_8_TCRB_E03.phenotype_unique.tsv'), 
				'15684_long': os.path.join(r2, '15684.matrix_plus.tsv')},
			'15754':{'15754_10x': os.path.join(r1,'15754_4_TCRB_E03.phenotype_unique.tsv'), 
				'15754_long': os.path.join(r2, '15754.matrix_plus.tsv')},
			'15761':{'15761_10x': os.path.join(r1,'15761_4_TCRB_E03.phenotype_unique.tsv'), 
				'15761_long': os.path.join(r2, '15761.matrix_plus.tsv')},
			'15836':{'15836_10x': os.path.join(r1,'15836_4_TCRB_E03.phenotype_unique.tsv'), 
				'15836_long': os.path.join(r2, '15836.matrix_plus.tsv')}}

## have to add 15782 last bc missing E05 
			'15782':{'15782_10x': os.path.join(r1,'15782_4_TCRB_E03.phenotype_unique.tsv'), 
				'15782_long': os.path.join(r2, '15782.matrix_plus.tsv')},
"""	
json - written with all samples, before distinction between long_ and matrix_ 	
json_object = json.dumps(samples, indent = 4) 
print(json_object)
with open('/fh/fast/corey_l/user/esford3/project_code/hybrid_t_cell_response/json_ef/sig_expand_merge.json', 'w') as outfile:
	outfile.write(json_object)
"""
#count=1

merge_list = list()
matrix_file = list()
index_cols = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family',
	'j_family', 'vMaxResolved', 'jMaxResolved']

for ptid in samples.keys():
	print(ptid)  
	x_filepath    = samples[ptid][f'{ptid}_10x']            
	matrix_filepath = samples[ptid][f'{ptid}_long']
	print(f"10X: {x_filepath} ")
	print(f"long : {matrix_filepath}")

	x_file = pd.read_csv(os.path.join(r1, x_filepath), sep = "\t")
	print(x_file.head(2))

	df = pd.read_csv(os.path.join(r2, matrix_filepath), sep = "\t")
	print(df.head(2))
	cs = [col for col in df.columns if col not in index_cols]
	e02_column = [c for c in cs if c.endswith("E02_vac_class") == True]
	e03_column = [c for c in cs if c.endswith("E03_vac_class") == True]
	e05_column = [c for c in cs if c.endswith("E05_vac_class") == True]
	e02_column = e02_column[0]
	e03_column = e03_column[0]
	e05_column = e05_column[0]

	#turn inclusion into boolean query
	df['include'] =   ( 
				((df[e02_column] == "sig_expand") | (df[e02_column] == "sig_contract")) |	
				((df[e03_column] == "sig_expand") | (df[e03_column] == "sig_contract")) |
				((df[e05_column] == "sig_expand") | (df[e05_column] == "sig_contract"))  
				)
	ix1=(df[e02_column] == "sig_expand")
	ix2=(df[e02_column] == "sig_contract")
	ix3=(df[e03_column] == "sig_expand")
	ix4=(df[e03_column] == "sig_contract")
	ix5=(df[e05_column] == "sig_expand")
	ix6=(df[e05_column] == "sig_contract")
	import numpy as np
	assert np.all(df['include']==(ix1 | ix2 | ix3 | ix4 | ix5 | ix6))

	dft = df.query('include == True')

	merge = pd.merge(dft, x_file.set_index('cdr3_b_aa_x'), left_on = ['cdr3_b_aa'], right_index=True)
	merge.to_csv(os.path.join(r2, f'{ptid}_sig_expand_cont_w10x_061722.tsv'), sep = "\t")

## for 15782 (no E05)

samples = {'15782':{'15782_10x': os.path.join(r1,'15782_4_TCRB_E03.phenotype_unique.tsv'), 
				'15782_long': os.path.join(r2, '15782.matrix_plus.tsv')}}

for ptid in samples.keys():
	print(ptid)  
	x_filepath    = samples[ptid][f'{ptid}_10x']            
	matrix_filepath = samples[ptid][f'{ptid}_long']
	print(f"10X: {x_filepath} ")
	print(f"long : {matrix_filepath}")

	x_file = pd.read_csv(os.path.join(r1, x_filepath), sep = "\t")
	print(x_file.head(2))

	df = pd.read_csv(os.path.join(r2, matrix_filepath), sep = "\t")
	print(df.head(2))
	cs = [col for col in df.columns if col not in index_cols]
	e02_column = [c for c in cs if c.endswith("E02_vac_class") == True]
	e03_column = [c for c in cs if c.endswith("E03_vac_class") == True]
	e02_column = e02_column[0]
	e03_column = e03_column[0]

	#turn inclusion into boolean query
	df['include'] =   ( 
				((df[e02_column] == "sig_expand") | (df[e02_column] == "sig_contract")) |	
				((df[e03_column] == "sig_expand") | (df[e03_column] == "sig_contract")) 
				)
	ix1=(df[e02_column] == "sig_expand")
	ix2=(df[e02_column] == "sig_contract")
	ix3=(df[e03_column] == "sig_expand")
	ix4=(df[e03_column] == "sig_contract")
	import numpy as np
	assert np.all(df['include']==(ix1 | ix2 | ix3 | ix4))

	dft = df.query('include == True')

	merge = pd.merge(dft, x_file.set_index('cdr3_b_aa_x'), left_on = ['cdr3_b_aa'], right_index=True)
	merge.to_csv(os.path.join(r2, f'{ptid}_sig_expand_cont_w10x_061722.tsv'), sep = "\t")


## make single files with either E03 sig_expand or sig_contract

import json
import os
import sys 
import pandas as pd
path = r"/fh/fast/corey_l/user/esford3/output/dfts/"
dest = r"/fh/fast/corey_l/user/esford3/output/dfts/"

## run this three times - for each vaccine time point 

## this version is for the phenotype files 
def generic_column_format(df, ptid , fc_threshold = 5):
    
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
    e02_columns = [x for x in df.columns if x.find("E02") != -1]
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
    df_generic['ptid'] = ptid
    #print(df_generic)
    return(df_generic)

### this version is for the longitudinal_matrix files 
def generic_column_format(df, ptid , fc_threshold = 4):
    
    generic_cols = {'cdr3_b_nucseq': 'cdr3_b_nucseq',
     'cdr3_b_aa': 'cdr3_b_aa',
     'v_b_gene': 'v_b_gene',
     'j_b_gene': 'j_b_gene',
     'v_family': 'v_family',
     'j_family': 'j_family',
     'vMaxResolved': 'vMaxResolved',
     'jMaxResolved': 'jMaxResolved',
     'predetected' : 'predetected',
     'mem' : 'mem'}
    
    e03_column = [x for x in df.columns if x.find("E03") != -1][0]
    e01_column = [x for x in df.columns if x.find("E01") != -1][0]
    
    e00_columns = [x for x in df.columns if x.find("E00") != -1]
    e005_columns = [x for x in df.columns if x.find("E00.5") != -1]
    e01_columns = [x for x in df.columns if x.find("E01") != -1]
    e02_columns = [x for x in df.columns if x.find("E02") != -1]
    e03_columns = [x for x in df.columns if x.find("E03") != -1]
    e04_columns = [x for x in df.columns if x.find("E04") != -1]
    e05_columns = [x for x in df.columns if x.find("E05") != -1]
    
    e00_columns_rep  =  {x : x.split("_TCRB_")[1] for x in e00_columns}
    e005_columns_rep =  {x : x.split("_TCRB_")[1] for x in e005_columns}
    e01_columns_rep = {x : x.split("_TCRB_")[1] for x in e01_columns}
    e02_columns_rep = {x : x.split("_TCRB_")[1] for x in e02_columns}
    e03_columns_rep = {x : x.split("_TCRB_")[1] for x in e03_columns}
    e04_columns_rep = {x : x.split("_TCRB_")[1] for x in e04_columns}
    e05_columns_rep = {x : x.split("_TCRB_")[1] for x in e05_columns}
    
    rep_cols = {**generic_cols, **e00_columns_rep,  **e005_columns_rep,  **e01_columns_rep, **e02_columns_rep, **e03_columns_rep,**e04_columns_rep,**e05_columns_rep}
    df_generic = df[rep_cols.keys()].rename(columns = rep_cols)
    df_generic['ptid'] = ptid
    #print(df_generic)
    return(df_generic)


files=os.listdir(path)
files=sorted(files)
print(files)

merge_list = list()
fdr_thr = 0.2

##sig_exp at E03 - have to exclude 15530
for file in files:
	if not file.endswith("matrix_plus.tsv"):
	  continue
	df = pd.read_csv(os.path.join(path, file), sep="\t")
	cs = [col for col in df.columns]
	e03_column = [c for c in cs if c.endswith("E03_vac_fdr") == True]
	e03_column = e03_column[0]
	df['include'] = ((df[e03_column] < fdr_thr))
	ptid = file.split(".", 1)[0]
	print(ptid)
	dft = df.query('include == True')
	df_generic = generic_column_format(dft, ptid = ptid)
	df_generic['ptid'] = ptid
	columns = df_generic.columns
	print(columns)
	merge_list.append(df_generic)

##sig_exp at E02 - have to exclude 15742 and 15758 (and 15530)
for file in files:
	if not file.endswith("matrix_plus.tsv"):
	  continue
	df = pd.read_csv(os.path.join(path, file), sep="\t")
	cs = [col for col in df.columns]
	e02_column = [c for c in cs if c.endswith("E02_vac_fdr") == True]
	e02_column = e02_column[0]
	df['include'] = ((df[e02_column] < fdr_thr))
	ptid = file.split(".", 1)[0]
	print(ptid)
	dft = df.query('include == True')
	df_generic = generic_column_format(dft, ptid = ptid)
	df_generic['ptid'] = ptid
	columns = df_generic.columns
	print(columns)
	merge_list.append(df_generic)

##sig_exp at E05 - exclude 15530, 15514, 15518, 15668, 15782, 15869
for file in files:
	if not file.endswith("matrix_plus.tsv"):
	  continue
	df = pd.read_csv(os.path.join(path, file), sep="\t")
	cs = [col for col in df.columns]
	e05_column = [c for c in cs if c.endswith("E05_vac_fdr") == True]
	e05_column = e05_column[0]
	df['include'] = ((df[e05_column] < fdr_thr))
	ptid = file.split(".", 1)[0]
	print(ptid)
	dft = df.query('include == True')
	df_generic = generic_column_format(dft, ptid = ptid)
	df_generic['ptid'] = ptid
	columns = df_generic.columns
	print(columns)
	merge_list.append(df_generic)

#for x in merge_list:
#	if not pd.Series(x.columns).value_counts().max()==1:
#		print(x, "you dummy")

all_merge = pd.concat(merge_list, axis=0)

index_cols = ['ptid', 'cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family',
	'j_family', 'vMaxResolved', 'jMaxResolved', 'predetected', 'mem']
cs = [col for col in all_merge.columns if col not in index_cols]
cs = [index_cols + sorted(cs)]
def flatten(d):
 v = [[i] if not isinstance(i, list) else flatten(i) for i in d]
 return [i for b in v for i in b]

cs = flatten(cs)

all_merge = all_merge.reindex(columns = cs)
all_merge.to_csv(os.path.join(dest, "all_e03_sig.2.tsv"), sep = "\t", index = False)
all_merge.to_csv(os.path.join(dest, "all_e02_sig.2.tsv"), sep = "\t", index = False)
all_merge.to_csv(os.path.join(dest, "all_e05_sig.2.tsv"), sep = "\t", index = False)
