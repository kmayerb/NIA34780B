"""
generate_significant_expander_contractor_table.py

This code computes the significantly expanded and contracted clones comparing the 
E01 and E03 timpoints.

Complexity is managed in set of json file for each participant that records 
bulk TCRb repertoire file paths.
"""

import json
import os 
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe

get_fold_changes = True
ptids = ['15665','15518','15771','15545','15753','15668','15635','15706',
         '15869','15784','15514','15582','15653','15742','15773','15531',
         '15839','15577','15559','15837','15754','15527','15669','15782',
         '15763','15758','15581','15836','15673','15684','15525','15761','15530']
for ptid in do_these:
    if os.path.isfile(f'ab_fig/{ptid}_e01_e03.pdf'):
        print(f"skiping {ptid}")
        continue 

    os.system(f"head json/{ptid}_samples.json")

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

    if get_fold_changes:
        dftextra = compute_fold_changes_on_longitudinal_dataframe(
            dft = dft.copy(), 
            dest = ".",
            ptid = None,
            write = False)
    
    """
    Save dftextra to external directory or subset to significant clones
    """

