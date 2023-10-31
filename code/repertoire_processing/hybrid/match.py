"""
Code Example showing how to match to an external reference
"""
if __name__ == "__main__":
    # TWO EXAMPLES:
    # EXAMPLE 1, USE THIS ON ONE PTID
    import os 
    import sys 
    user = "kmayerbl"
    if user == 'esford3':
        git_path = "/fh/fast/corey_l/user/esford3/project_code/hybrid_t_cell_response"
        sys.path.insert(0,os.path.abspath(git_path))
    elif user == 'kmayerbl':
        git_path = "/fh/fast/gilbert_p/fg_data/koelle_covid/project_code/hybrid_t_cell_response"
        sys.path.insert(0,os.path.abspath(git_path))
    else:
        raise ValueError("User must be kmayerbl or esford3")


    import json
    from hybrid.parse import parse_adaptive_v2
    from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
    from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe
    from hybrid.calculate_depth_efficiently   import calculate_depth_breadth

    # IN the following block I reconstruct the dftextra from the raw files
    ptid = '15673'
    json_file = os.path.join(git_path, f'json/{ptid}_samples.json')

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

    import pandas as pd 
    cmv = pd.read_csv(os.path.join(git_path, 'reference', 'cmv_reference_sequences.csv'))
    # Define matching key. 
    # NOTE IN THIS CASE I DON'T TRUST GENE ASIGNMENT ACROSS MULTIPLE METHODS, SO WE DEFAULT TO FAMILY 
    # LEVEL ASSIGNMENT FOLLOWING METHOD BY DEWITT AND BRADLEY.
    # Get get V25 from TRBV25-1*01
    

    cmv['v_fam'] = cmv['v_b_gene'].apply(lambda x: x.split("*")[0].split("-")[0].replace("TRBV", "V"))
    cmv['key'] = cmv['v_fam'] + "," + cmv['cdr3_b_aa']
    # Here we keep track of all of the HLA alleles associated with a particular TCR key.
    cmv_key_to_hla_2 = dict() 
    for i,r in cmv.iterrows():
        if r['hla_2'] is not None:
            cmv_key_to_hla_2.setdefault(r['key'],[]).append(r['hla_2']) 
    # convert list of all viable to set so we can test for intersection wiht the ptid's allele set
    cmv_key_to_hla_2 = {k: set(v) for k,v in cmv_key_to_hla_2.items()}

    cmv_key_to_hla_4 = dict() 
    for i,r in cmv.iterrows():
        if r['hla_4'] is not None:
            cmv_key_to_hla_4.setdefault(r['key'],[]).append(r['hla_4']) 
    cmv_key_to_hla_4 = {k: set(v) for k,v in cmv_key_to_hla_4.items()}

    
    dft['v_fam'] = dft['v_b_gene'].apply(lambda x: x.split("*")[0].split("-")[0].replace("TRBV", "V"))
    dft['key'] = dft['v_fam'] + "," + dft['cdr3_b_aa']
    dft['_cmv_match'] = dft['key'].isin(cmv['key'])

    # I have made these dictionaries previously, so we can resuse that code
    # This provides set of all alleles in each ptid
    from hybrid.feasible import get_hla_dictionaries
    sample_hla_dict, sample_hla_dict_4d, sample_hla_dict_2d = sample_to_hla = get_hla_dictionaries(
        tsv_file = os.path.join(git_path, 'reference' ,'hla_file.tsv'),
        cols = ['HLA-Aa', 'HLA-Ab', 'HLA-Ba', 'HLA-Bb', 'HLA-Ca', 'HLA-Cb', 'DQB1b', 'DPA1a', 'DPA1b', 'DPB1a', 'DPB1b'], 
        sample_key = 'sample')
    # _cmv_hla_2, what is the 2-digit hla associated with the key
    dft['_cmv_hla_2'] = dft['key'].apply(lambda x : cmv_key_to_hla_2.get(x))
    # _cmv_hla_4, what is the 4-digit hla associated with the key
    dft['_cmv_hla_4'] = dft['key'].apply(lambda x : cmv_key_to_hla_4.get(x))
    # Now check if it matches
    # _cmv_hla_2_matches one of the ptid's HLA alleles (2-digit) Tue or False
    # Explanation. Here we lookup the HLA associated with the TCR key. Then we ask if it is in set of HLA allels for that ptid
    # See if there is an intersection between the HLAs associated with the TCR key and set of HLA allels in that ptid
    dft['_cmv_hla_2_match'] = dft['key'].apply(lambda x : cmv_key_to_hla_2.get(x).intersection(sample_hla_dict_2d.get(int(ptid))) != set() if cmv_key_to_hla_2.get(x) is not None else False)
    # _cmv_hla_2_matches one of the ptid's HLA alleles (2-digit) Tue or False
    dft['_cmv_hla_4_match'] = dft['key'].apply(lambda x : cmv_key_to_hla_4.get(x).intersection(sample_hla_dict_4d.get(int(ptid))) != set() if cmv_key_to_hla_4.get(x) is not None else False)
    
    # eg., if you want to check
    # dft['ptid_hla_2'] = ";".join(sorted(list(sample_hla_dict_2d.get(int(ptid)))))
    # dft[dft['_cmv_hla_2_match'] == True]
    # dft[dft['_cmv_match'] == True]










