
"""
## assign memory/naive and previously_detected based on E00 and E01 columns 

"""

def assign_memory(df, ptid = None, dest = ".", write = False, memory_min = 2): 
   """
   Parameters
   ----------
   df : output from make_wide_longitudinal_array 
      mat  output from compute_fold_change 
   ptid: str or None
   
   memory_min : int
      default is zero (anything seen more than N times in pre-vaccine time point is considered memory)
      emily wants >2   
      kosh  wants >0 (Kosh concedes to emily if we have predetected catagory)
   dest : str
      filepath where to output result if write is True
   write : bool 
      if True, write out file as f"{ptid}.matrix_plus.tsv"
   Returns
   -------

   Example 
   -------
   #test case
   import pandas as pd
   import os
   path = r"/fh/fast/corey_l/user/esford3/output/dfts/"
   dest = r"/fh/fast/corey_l/user/esford3/output/dfts/"
   file = "15518.longitudinal_matrix.tsv"
   df = pd.read_csv(os.path.join(path, file), sep="\t")
   ptid = file.split("_", 1)[0]
   """
   ##first - define the columns needed (counts) 
   index_cols = ['cdr3_b_nucseq', 'cdr3_b_aa', 'v_b_gene', 'j_b_gene', 'v_family',
          'j_family', 'vMaxResolved', 'jMaxResolved']
   cs = [col for col in df.columns if col not in index_cols]
   cs = [col for col in cs if col.find("_pfreq") == -1]
   ## to make it more adaptable to matrix or non-matrix files, search for timepoints individually  
   ## NICE, THIS IS SMART!
   print("CLASSIFYING MEMORY OR NAIVE RELATIVE TO PRE VACCINE TIME POINTS")
   print(f"ptid = {ptid}")
   e00_column = [c for c in cs if c.endswith("E00") == True]
   e01_column = [c for c in cs if c.endswith("E01") == True]
   e02_column = [c for c in cs if c.endswith("E02") == True]
   e03_column = [c for c in cs if c.endswith("E03") == True]
   e05_column = [c for c in cs if c.endswith("E05") == True]

   ## to compute for longitudinal files - in case we want to look at E04 vs E05 
   # TODO : Emily explain logic of line above ??? 
   # If there is no e00 column, then use the E00.5
   if e00_column == []: 
      print("relying on e00.5 column")
      e00_column = [c for c in cs if c.endswith("E00.5") == True]
   # If there is no E01 colummn use e00.5 as e01    
   if e01_column == []: 
      print("relying on e00.5 column")
      e01_column = [c for c in cs if c.endswith("E00.5") == True]
   
   ## here - if a clonotype is detected at 3 or more copies in E00 or E01, we'll consider it likely to be a memory clonotype,
   ## and likewise if it is detected multiple times (E00 and E01, or E00.5 and E01), then it seems likely to be memory 
   ## can adjust these parameters here -   
   ## for less fuzzy distinction between memory vs naive, use just whether it was previously detected or not    

   # TODO: I believe if you saw it at all in E00 or E01 it is likely memory. 
   # So I made this part of the function memory_min = 0  
   # I don't think we need loops 

   # There is only on e00_column Right?
   e00_column = e00_column[0]
   e01_column = e01_column[0]
   # is predetected if obseved in either e01 or 
   df['is_predetected'] = ( (df[e00_column] > 0)  | (df[e01_column] > 0) ) 
   # convert boolean to string
   df['predetected']    =  df['is_predetected'].apply(lambda x: "previously_detected" if x else "new")
   df['is_memory'] = (
                     ( (df[e00_column] > memory_min) | (df[e01_column] > memory_min) ) | 
                     ( (df[e00_column] > 0)          & (df[e01_column] > 0) ) 
                     )
   # convert boolean to string
   df['mem']       =  df['is_memory'].apply(lambda x: "memory" if x else "naive")
   return df
   
   if write:
      outfile = os.path.join(dest, f"{ptid}.matrix_plus.tsv")
      print(f"OUTPUT {outfile}")
      df.to_csv(outfile, sep = "\t", index = False)

