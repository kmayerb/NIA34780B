

"""
## link CMV sequences to matrix_plus files

"""

def match_cmv(df, ptid = None, dest = ".", write = False): 

	list_of_cmv_keys = cmv_keys.values.tolist()
	df['key'] = df['v_b_gene'] + "+" + df['cdr3_b_aa'] + "+" + df['j_b_gene']
	df['_match_cmv'] = df['key'].isin(list_of_cmv_keys)

	print(f"ptid = {ptid}")
	print(df['_match_cmv'].value_counts())
	return df
	
	results = dict()
	results_list = list()
	results[f'{ptid}_match_cmv'] = sum(df['_match_cmv'])
	results_list.append((ptid, results[f'{ptid}_match_cmv']))
	results_df = pd.DataFrame(results_list, columns = ['ptid','matches'])
	return(results_df)

### second try 
def match_cmv(df, ptid = None, dest = ".", write = False):
	list_of_cmv_keys = cmv_keys.values.tolist()
	df['key'] = df['v_b_gene'] + "+" + df['cdr3_b_aa'] + "+" + df['j_b_gene']
	df['_match_cmv'] = df.apply(lambda x: int(x['key'] in list_of_cmv_keys), axis=1)

	print(f"ptid = {ptid}")
	print(df['_match_cmv'].value_counts())
	return df
	
	results = dict()
	results_list = list()
	results[f'{ptid}_match_cmv'] = sum(df['_match_cmv'])
	results_list.append((ptid, results[f'{ptid}_match_cmv']))
	results_df = pd.DataFrame(results_list, columns = ['ptid','matches'])
	return(results_df)

   """
   parameters
   	df - any longitudinal matrix file, 
   		here will use matrix_plus file with sig/not sig, memory/naive  
   	ptid - available from file name (in howto runfile)
   	dest = location for write file 
   	cmv_keys = in data_file (combination of VDJDB db entries based on multimer id studies
   		and Emerson/DeWitt computationally-defined sequences (767 unique sequences))
   		cmv_keys includes v gene (TCRdist format), cdr3b, j gene (TCRdist format)

   #test case
   import pandas as pd
   import os
   path = r"/fh/fast/corey_l/user/esford3/output/dfts/"
   dest = r"/fh/fast/corey_l/user/esford3/output/dfts/"
   file = "15518.matrix_plus.tsv"
   cmv_keys = pd.read_csv(os.path.join(path, "cmv_keys.csv"))
   df = pd.read_csv(os.path.join(path, file), sep="\t")
   ptid = file.split(".", 1)[0]

   """