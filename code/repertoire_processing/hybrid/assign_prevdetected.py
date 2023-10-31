
"""
## assign memory/naive and previously_detected based on E00 and E01 columns 

"""

def assign_memory(df, ptid = None, dest = ".", write = False): 

   """
   parameters
   df - output from make_wide_longitudinal_array 
   mat - output from compute_fold_change 

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

   print("CLASSIFYING MEMORY OR NAIVE RELATIVE TO PRE VACCINE TIME POINTS")
   print("ptid = ", ptid)
   e00_column = [c for c in cs if c.endswith("E00") == True]
   e01_column = [c for c in cs if c.endswith("E01") == True]
   e02_column = [c for c in cs if c.endswith("E02") == True]
   e03_column = [c for c in cs if c.endswith("E03") == True]
   e05_column = [c for c in cs if c.endswith("E05") == True]

   ## to compute for longitudinal files - in case we want to look at E04 vs E05 

   if e00_column == []: 
      print("relying on e00.5 column")
      e00_column = [c for c in cs if c.endswith("E00.5") == True]
   else: 
      if e01_column == []: 
         print("relying on e00.5 column")
         e01_column = [c for c in cs if c.endswith("E00.5") == True]
      else: 

   ## here - if a clonotype is detected at 3 or more copies in E00 or E01, we'll consider it likely to be a memory clonotype,
   ## and likewise if it is detected multiple times (E00 and E01, or E00.5 and E01), then it seems likely to be memory 
   ## can adjust these parameters here -   
   ## for less fuzzy distinction between memory vs naive, use just whether it was previously detected or not     
         for i in range(len(e00_column)):
               memory_df = (((df[f'{e00_column[i]}'] > 2) == True) | ((df[f'{e01_column[i]}'] > 2) == True) 
                        | ((df[f'{e00_column[i]}'] > 0) & (df[f'{e01_column[i]}'] > 0) == True))
               if  memory_df[i] != False:
                  df['mem'] = "memory"
               else:
                  df['mem'] = "naive"
               print(df.head(2))
         for i in range(len(e00_column)):
               prev_det = (((df[f'{e00_column[i]}'] > 0) == True) | ((df[f'{e01_column[i]}'] > 0) == True))
               if prev_det[i] != False: 
                  df['prev_detect'] = "previously_detected"
               else: 
                  df['prev_detect'] = "new"
               print(df.head(2))
               return df

         if write:
                 outfile = os.path.join(dest, f"{ptid}.matrix_plus.tsv")
                 print(f"OUTPUT {outfile}")
                 df.to_csv(outfile, sep = "\t", index = False)

