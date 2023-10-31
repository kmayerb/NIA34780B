
"""
calculate_depth_breadth.py
"""

import os
import pandas as pd
def calculate_depth_breadth(df, ptid, fold_change = 4, verbose = True):
   """
   Parameters
   ----------
   df : pd.DataFrame 
      The DataFrame that is produced from 
         hybrid.assemble_longitudinal_samples.make_wide_longitudinal_dataframe()
         hybrid.assemble_longitudinal_samples.compute_fold_changes_on_longitudinal_dataframe()

      These pd.DataFrame can be read from precomputed output files (e.g., 15869.longitudinal_matrix.tsv ). 
         This was done only for people with E01 reference timepoint
   ptid : str
      supply a string for the ptid, so that the output dataframe is tagged by ptid
   fold_change : float 
      The increase (fold_change: default 2) or decrease (1/fold_change: default 0.5), filters 
      statistically significant change in counts 
   verbose: boolean
      If tru, prints out helpful information for understanding how the function is working

   Returns
   -------
   matrixdb : pd.DataFrame
       Wide format longitudinal data
   """
   # We need to get the names of the samples columns, so we exclude index and _pfreq columns    
   index_cols = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family','j_family', 'vMaxResolved', 'jMaxResolved']
   cs = [col for col in df.columns if col not in index_cols]
   ca = [col for col in cs if col.find('_temporal_class') != -1]
   cb = [col for col in cs if col.find("_temporal_fc") != -1]
   # CA and CB are one entry shorter than the other vectors, so we add a dummy
   ca.insert(0,None)
   cb.insert(0,None)
   cc = [col for col in cs if col.find("_vac_class") != -1]
   cd = [col for col in cs if col.find("_vac_fc") != -1]
   ce = [col for col in cs if col.find('_pfreq') != -1]
   cf = [col for col in cs if col.find('_sum') != -1]
   cg = [col for col in cs if col.find('_detect') != -1]

   assert len(ca) == len(cb)
   assert len(ca) == len(cc)
   assert len(ca) == len(cd)
   assert len(ca) == len(ce)
   assert len(ca) == len(cf)
   assert len(ca) == len(cg)
   
   if verbose:
      print('<ca> _temporal_class COLUMNS:')
      print(ca)
      print('<cb> _temporal_fc COLUMNS:')
      print(cb)
      print('<cc> _vac_class COLUMNS:')
      print(cc)
      print('<cd> _vac_fc COLUMNS:')
      print(cd)
      print('<ce> _pfreq COLUMNS:')
      print(ce)
      print('<cf> _sum COLUMNS:')
      print(cf)
      print('<cg> _detect COLUMNS:')
      print(cg)

   results = dict()
   results_list = list()
   for i in range(1,len(ca)): 
      # ca - _temporal_fc
      # cb - _vac_class (classification of significant expansion)
      # cc _vac_class
      # cd _vac_fc
      # ce - _pfreq COLUMNS
      # cg - _detect 

      # NAME THE SAMPLE <s>
      s = ca[i].replace('_temporal_class','')
      print(s)

      # We start by computing the breadth and depth of expansion for the "temporal_class". 
      # (That expansion from one sample to the next sample) 
      # Quanities of Interest
      #   _temporal_sigexpand_count
      #   _temporal_sigcontract_count
      #   _temporal_sigexpand_breadth
      #   _temporal_sigexpand_breadth2 (denominator is total clones in i, and i-1 sample)
      #   _temporal_sigcontract_breadth
      #   _temporal_sigexpand_pfreq
      #   _temporal_sigcontract_pfreq
      # 

      # First we create a boolean index of True for all significant 2X fold change expanders
      ix_expand   = (df[f'{ca[i]}'] == 'sig_expand') & (df[f'{cb[i]}'] > fold_change) 
      # Next we create a boolaen index of True for all signficiant 0.5X fold change contractors
      ix_contract = (df[f'{ca[i]}'] == 'sig_contract') & (df[f'{cb[i]}'] < 1/fold_change)  
      # Count how many signficant contractors , Sum of boolean is count of True values.
      results[f'{s}_temporal_sigexpand_count']     = sum(ix_expand)
      results[f'{s}_temporal_sigcontract_count']   = sum(ix_contract)
      # Compute the sum of productive frequencey (ce[i])
      results[f'{s}_temporal_sigexpand_pfreq']     = sum(df[f'{ce[i]}'][ix_expand])
      results[f'{s}_temporal_sigcontract_pfreq']   = sum(df[f'{ce[i]}'][ix_contract])
      # Breadth, here is defined as number detected clones (numerator) divided over the total number of unique clones in that sample (denominator)
      # cg column is detection collumn, so sum of Trues in this column is a denominator for our breadth measure. 
      # TODO: Is there 
      total_clones_dected_in_sample_i = sum(df[f'{cg[i]}'])
      total_clones_dected_in_sample_i_minus_1 = sum(df[f'{cg[i-1]}'])
      results[f'{s}_temporal_sigexpand_breadth']     = sum(ix_expand) / total_clones_dected_in_sample_i 
      results[f'{s}_temporal_sigcontract_breadth']   = sum(ix_contract) / total_clones_dected_in_sample_i 
      results[f'{s}_temporal_sigexpand_breadth2']    = sum(ix_expand) /   (total_clones_dected_in_sample_i + total_clones_dected_in_sample_i_minus_1)
      results[f'{s}_temporal_sigcontract_breadth2']  = sum(ix_contract) / (total_clones_dected_in_sample_i + total_clones_dected_in_sample_i_minus_1) 
      
      if verbose: 
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigexpand_count'],   '_temporal_sigexpand_count'     ))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigexpand_pfreq'],   '_temporal_sigexpand_pfreq'     ))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigexpand_breadth'], '_temporal_sigexpand_breadth'   ))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigexpand_breadth2'], '_temporal_sigexpand_breadth2'   ))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigcontract_count'],   '_temporal_sigcontract_count'   ))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigcontract_pfreq'],   '_temporal_sigcontract_pfreq'   ))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigcontract_breadth'], '_temporal_sigcontract_breadth'))
         print((i, s, ca[i], cb[i], cg[i], results[f'{s}_temporal_sigcontract_breadth2'], '_temporal_sigcontract_breadth2'))

      results_list.append((i, s, results[f'{s}_temporal_sigexpand_count'],      '_temporal_sigexpand_count'    ))
      results_list.append((i, s, results[f'{s}_temporal_sigexpand_pfreq'],      '_temporal_sigexpand_pfreq'    ))
      results_list.append((i, s, results[f'{s}_temporal_sigexpand_breadth'],    '_temporal_sigexpand_breadth'  ))
      results_list.append((i, s, results[f'{s}_temporal_sigexpand_breadth2'],   '_temporal_sigexpand_breadth2' ))
      results_list.append((i, s, results[f'{s}_temporal_sigcontract_count'],    '_temporal_sigcontract_count'  ))
      results_list.append((i, s, results[f'{s}_temporal_sigcontract_pfreq'],    '_temporal_sigcontract_pfreq'  ))
      results_list.append((i, s, results[f'{s}_temporal_sigcontract_breadth'],  '_temporal_sigcontract_breadth'))
      results_list.append((i, s, results[f'{s}_temporal_sigcontract_breadth2'], '_temporal_sigcontract_breadth2'))

   # For later computation we need to know total number of clones in the E01 reference 
   # sample 
   e01_detect_column = [x for x in cg if x.find('E01') !=-1][0]
   total_clones_dected_in_sample_E01 = sum(df[e01_detect_column])

   for i in range(len(cc)):
      # NAME THE SAMPLE <s>
      s = cc[i].replace('_vac_class','')
      print(s)
      # We start by computing the breadth and depth of expansion/contraction to E01 (pre-vaccine)
      # (That expansion from one sample to the next sample) 
      # Quanities of Interest
      #   _vac_sigexpand_count
      #   _vac_sigcontract_count
      #   _vac_sigexpand_breadth
      #   _vac_sigexpand_breadth2  (denominator is total clones in i, and E01 sample)
      #   _vac_sigcontract_breadth
      #   _vac_sigexpand_pfreq
      #   _vac_sigcontract_pfreq

      ix_expand   = (df[f'{cc[i]}'] == 'sig_expand') & (df[f'{cd[i]}'] > fold_change) 
      # Next we create a boolaen index of True for all signficiant 0.5X fold change contractors
      ix_contract = (df[f'{cc[i]}'] == 'sig_contract') & (df[f'{cd[i]}'] < 1/fold_change)  
      # Count how many signficant contractors , Sum of boolean is count of True values.
      results[f'{s}_vac_sigexpand_count']     = sum(ix_expand)
      results[f'{s}_vac_sigcontract_count']   = sum(ix_contract)
      # Compute the sum of productive frequencey (ce[i])
      results[f'{s}_vac_sigexpand_pfreq']     = sum(df[f'{ce[i]}'][ix_expand])
      results[f'{s}_vac_sigcontract_pfreq']   = sum(df[f'{ce[i]}'][ix_contract])
      # Breadth, here is defined as number detected clones (numerator) divided over the total number of unique clones in that sample (denominator)
      # cg column is detection collumn, so sum of Trues in this column is a denominator for our breadth measure. 
      # TODO: Is there 
      total_clones_dected_in_sample_i = sum(df[f'{cg[i]}'])
      results[f'{s}_vac_sigexpand_breadth']     = sum(ix_expand) / total_clones_dected_in_sample_i 
      results[f'{s}_vac_sigcontract_breadth']     = sum(ix_contract) / total_clones_dected_in_sample_i 
      results[f'{s}_vac_sigexpand_breadth2']     = sum(ix_expand) /     (total_clones_dected_in_sample_i + total_clones_dected_in_sample_E01 )
      results[f'{s}_vac_sigcontract_breadth2']     = sum(ix_contract) / (total_clones_dected_in_sample_i + total_clones_dected_in_sample_E01 )
      if verbose: 
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigexpand_count'],     '_vac_sigexpand_count'     ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigexpand_pfreq'],     '_vac_sigexpand_pfreq'     ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigexpand_breadth'],   '_vac_sigexpand_breadth'   ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigexpand_breadth2'],   '_vac_sigexpand_breadth2' ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigcontract_count'],   '_vac_sigcontract_count'   ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigcontract_pfreq'],   '_vac_sigcontract_pfreq'   ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigcontract_breadth'], '_vac_sigcontract_breadth' ))
         print((i, s, cc[i], cd[i], cg[i], results[f'{s}_vac_sigcontract_breadth2'], '_vac_sigcontract_breadth2' ))

      # Package the results, 
         # index
         # s - sample_name
         # value - computed from above
      results_list.append((i, s, results[f'{s}_vac_sigexpand_count'],      '_vac_sigexpand_count'     ))
      results_list.append((i, s, results[f'{s}_vac_sigexpand_pfreq'],      '_vac_sigexpand_pfreq'     ))
      results_list.append((i, s, results[f'{s}_vac_sigexpand_breadth'],    '_vac_sigexpand_breadth'   ))
      results_list.append((i, s, results[f'{s}_vac_sigexpand_breadth2'],    '_vac_sigexpand_breadth2' ))
      results_list.append((i, s, results[f'{s}_vac_sigcontract_count'],    '_vac_sigcontract_count'   ))
      results_list.append((i, s, results[f'{s}_vac_sigcontract_pfreq'],    '_vac_sigcontract_pfreq'   ))
      results_list.append((i, s, results[f'{s}_vac_sigcontract_breadth'],  '_vac_sigcontract_breadth' ))
      results_list.append((i, s, results[f'{s}_vac_sigcontract_breadth2'], '_vac_sigcontract_breadth2'))

   results_df = pd.DataFrame(results_list, columns = ['i','sample_name', 'value', 'variable'])
   results_df['ptid'] = ptid
   return(results_df)



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

   # optional to compute expansion statistics
   dftextra = compute_fold_changes_on_longitudinal_dataframe(
      dft = dft.copy(), 
      dest = ".",
      ptid = None,
      write = False)

   ptid_breadth_and_depth_df = calculate_depth_breadth(df = dftextra, ptid = ptid)
   print(ptid_breadth_and_depth_df)
   # EXAMPLE 2, PERFORM BREADTH AND DEPTH ON MULTIPLE PRECOMPUTED FILES
   # TODO: EMILY PLEASE TEST THIS I wasn't able to test this yet because I don't have read access
   
   import os 
   import sys 
   import pandas as pd
   from hybrid.calculate_depth_efficiently   import calculate_depth_breadth
   user = "esford3"
   if user == 'esford3':
      git_path = "/fh/fast/corey_l/user/esford3/project_code/hybrid_t_cell_response"
      sys.path.insert(0,os.path.abspath(git_path))
   elif user == 'kmayerbl':
      git_path = "/fh/fast/gilbert_p/fg_data/koelle_covid/project_code/hybrid_t_cell_response"
      sys.path.insert(0,os.path.abspath(git_path))
   else:
      raise ValueError("User must be kmayerbl or esford3")

   path = "/fh/fast/corey_l/user/esford3/output/dfts/"
   dest = "/fh/fast/corey_l/user/esford3/output/dfts/"
   fs = [f for f in os.listdir(path) if f.endswith('.longitudinal_matrix.tsv')]
   all_results = list()
   for file in fs:
      dftextra = pd.read_csv(os.path.join(path,file), sep = "\t")
      # extract ptid from filename 15869.longitudinal_matrix.tsv
      ptid = file.split(".")[0]
      ptid_breadth_and_depth_df = calculate_depth_breadth(df = dftextra, ptid = ptid, fold_change = 4)
      all_results.append(ptid_breadth_and_depth_df)
   
   # Put all the results into a single DataFrame
   all_ptid_breadth_and_depth_df = pd.concat(all_results)
   all_ptid_breadth_and_depth_df.to_csv(
      os.path.join(dest, "20220615_breadth_and_depth_summary.tsv"), 
      sep = "\t", 
      index = False)


