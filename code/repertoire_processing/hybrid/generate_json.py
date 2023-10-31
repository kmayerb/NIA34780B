"""
This is just for generating json files used in the program elsewhere
"""
import os
import pandas as pd
m = pd.read_csv('json/sample_key_041022.csv')
m['sample_name'] = m['person'].apply(lambda x : str(x)) + "_" \
    + m['sample'].apply(lambda x : str(x)) + "_TCRB_" + \
    m['extension code'].apply(lambda x: str(x) if isinstance(x, str) else "X")

m['filename'] = m['person'].apply(lambda x : str(x)) + "_" \
    + m['sample'].apply(lambda x : str(x)) + "_TCRB.tsv" 

rearly = '/fh/fast/gilbert_p/fg_data/koelle_covid/covid_tcr_compression'
rboost = '/fh/fast/gilbert_p/fg_data/koelle_covid/booster_repertoires_v2/'
samples = dict()
for p,sn,fn,ecode in zip(m['person'], m['sample_name'], m['filename'], m['extension code']):
    if ecode in ["E04", "E05"]:
        path = rboost 
    else: 
        path = rearly
    
    fp = os.path.join(path, fn)
    if os.path.isfile(fp):
        samples.setdefault(p, {})[sn] = fp
    #else:
    #    samples.setdefault(p, {})[sn] = None
    print(p,sn,fn,fp)

for k,v in samples.items():
    with open(f"json/{k}_samples.json", "w") as fh:
        fh.write(json.dumps(v, indent = 4) )

import json
json_object = json.dumps(samples, indent = 4) 
print(json_object)
with open("json/2022_04_12_all_samples.json", "w") as fh:
    fh.write(json_object )



