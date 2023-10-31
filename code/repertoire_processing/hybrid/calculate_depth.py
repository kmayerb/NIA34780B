
"""
calculate_depth_breadth.py

"""
import sys
import os
sys.path.insert(0,os.path.abspath("/fh/fast/corey_l/user/esford3/project_code/hybrid_t_cell_response"))
#sys.path.append("/fh/fast/corey_l/user/esford3/project_code/hybrid_t_cell_response/hybrid")

import pandas as pd
import fishersapi
import statsmodels.stats.multitest as mt
import numpy as np
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe

"""
load dft matrix files
"""

import json
import os
import hybrid.parse
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe


path = r"/fh/fast/corey_l/user/esford3/output/dfts/"
dest = r"/fh/fast/corey_l/user/esford3/output/dfts/"

files=os.listdir(path)
files=sorted(files)
print(files)

for file in files:
   if not file.endswith("matrix.tsv"):
      continue
   df = pd.read_csv(os.path.join(path, file), sep="\t")
   ptid = file.split("_", 1)[0]

break

   dft = calculate_breadth_depth()

"""testing
import os
import pandas as pd
path = r"/fh/fast/corey_l/user/esford3/output/dfts/"
file = "15518.longitudinal_matrix.tsv"
df = pd.read_csv(os.path.join(path, file), sep="\t")
ptid = file.split("_", 1)[0]
"""

def calculate_depth_breadth(input = df, ptid = ptid, fold_change = 2)
    """
    Parameters
    ----------
    input : matrix file
         output from dftextra output (participants with E01 sample)
    dest : str
        place to write expanded dataframe 
    ptid : str or None
        interpolated from filename 
    fold_change : float 
        additional filter level 

    Returns
    -------
    matrixdb : pd.DataFrame
        Wide format longitudinal data
    """

    # We need to get the names of the samples columns, so we exclude index and _pfreq columns    
index_cols = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family',
       'j_family', 'vMaxResolved', 'jMaxResolved']
cs = [col for col in df.columns if col not in index_cols]
ca = [col for col in cs if col.find('_temporal_class') != -1]
cb = [col for col in cs if col.find("_temporal_fc") != -1]
cc = [col for col in cs if col.find("_vac_class") != -1]
cd = [col for col in cs if col.find("_vac_fc") != -1]
ce = [col for col in cs if col.find('_pfreq') != -1]
cf = [col for col in cs if col.find('_sum') != -1]
cg = [col for col in cs if col.find('_detect') != -1]
fold_change = 2 

## temporal_class first 
## calculate depth and breadth by identifying the number of clonotypes that are both 
## significantly expanded (by class call) and increased by >2 fold (by _fc_)
## calculate for each comparison 
## define the # expanding (_sigsum) here to pull into breadth calculation
## define the # of clonotypes in both samples of comparison (_inboth)
## calculate breadth for each temporal comparison _tempbreadth 

for i in range(len(ca)):
   col_tupple = (ca[i], cb[i], ce[i+1], ce[i])
   print(col_tupple)
   df[f'{ca[i]}_sigsum'] = len((df[f'{ca[i]}'] == 'sig_expand') & (df[f'{cb[i]}'] > fold_change))
   df[f'{ca[i]}_ineither'] = len((df[f'{cg[i+1]}'] == 'True') | (df[f'{cg[i]}'] == 'True')) 
   df[f'{ca[i]}_tempbreadth'] = (df[f'{ca[i]}_sigsum']) / (df[f'{ca[i]}_ineither'])

for i in range(len(ca)):
   col_tupple = (ca[i], cb[i], ce[i+1], ce[i])         
   if [(df[f'{ca[i]}'] == 'sig_expand') & (df[f'{cb[i]}'] > fold_change)] != -1:  
         df[f'{ca[i]}_tempdepth'] = df[f'{ce[i+1]}'] + df[f'{ce[i]}']
   else: 
         df[f'{ca[i]}_tempdepth'] = 0
 
## vax_class - same structure, just need to refer to E01 for _pfreq for depth
## sigsum = # of clonotypes that are sig expanding and meet _fc criteria (2)

e01_column = [c for c in ce if c.find("E01") != -1]

for i in range(len(cc)):
   col_tupple = (cc[i], cd[i], ce[i])
   print(col_tupple)
   df[f'{cc[i]}_sigsum'] = len((df[f'{cc[i]}'] == "sig_expand") & (df[f'{cd[i]}'] > fold_change))
   df[f'{cc[i]}_ineither'] = len((df[f'{cg[i]}'] == 'True') | (df[f'{e01_column[0]}'] == 'True')) 
   df[f'{cc[i]}_vaxbreadth'] = (df[f'{cc[i]}_sigsum']) / (df[f'{cc[i]}_ineither'])
     
for i in range(len(cc)):
   col_tupple = (cc[i], cd[i], ce[i], e01_column[0]) 
   print(col_tupple)   
   if [(df[f'{cc[i]}'] == 'sig_expand') & (df[f'{cd[i]}'] > fold_change)] != -1:   
         df[f'{cc[i]}_vaxdepth'] = df[f'{ce[i]}'] + df[f'{e01_column[0]}'] 
   else: 
         df[f'{cc[i]}_vaxdepth'] = 0


## make spreadsheet that puts ptid with _tempdepth and _vaxdepth
cs = [col for col in df.columns if col not in index_cols]
cg = [col for col in cs if col.find('depth') != -1]

dftemp = pd.DataFrame('a':ptid, )

dfoutput_cols = ['' ]

dft[f"{col}_sum"] = dft[col].sum()


## same with breadth - will be on sum 



