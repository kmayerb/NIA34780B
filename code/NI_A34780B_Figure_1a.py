"""
# Figure 1a.py

# E01-to-E03 Expansion plot for P673

# Note the original code to generate E01-to-E02, E01-to-E03 
# clonal expansion plots:
# /repertoire-processing/hybrid_py_scripts/generate_bulk_expander_contractor_plots.py

# Because this involves multiple large adaptive files, it must be run with access 
# to hybrid python modules:  /repertoire-processing/hybrid/
"""

import json
import os 
from hybrid.parse import parse_adaptive_v2
from hybrid.assemble_longitudinal_samples import make_wide_longitudinal_dataframe
from hybrid.assemble_longitudinal_samples import compute_fold_changes_on_longitudinal_dataframe
from plotnine import *

plots = []
summary_data = []
long_data_frames= []
count = 0
get_fold_changes = True

do_these = ['15673'] # ['15665','15518','15771','15545','15753','15668','15635','15706','15869','15784','15514','15582','15653','15742','15773','15531','15839','15577','15559','15837','15754','15527','15669','15782','15763','15758','15581','15836','15673','15684','15525','15761','15530']
for ptid in do_these:#

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

    ggdf = dftextra#[cols] #dftextra[['15673_6_TCRB_E01_pfreq', '15673_8_TCRB_E03_pfreq', '15673_8_TCRB_E03_vac_class','15673_8_TCRB_E03_vac_fc']]
    cols = dftextra.columns
    try:
        
        e01_pfreq_col = [k for k in cols if k.find("E01_pfreq") != -1][0]
        e02_pfreq_col = [k for k in cols if k.find("E02_pfreq") != -1][0]
        e03_pfreq_col = [k for k in cols if k.find("E03_pfreq") != -1][0]
        e02_vac_class = [k for k in cols if k.find("E02_vac_class") != -1][0]
        e03_vac_class = [k for k in cols if k.find("E03_vac_class") != -1][0]
        e02_vac_fc    = [k for k in cols if k.find("E02_vac_fc") != -1][0]
        e03_vac_fc    = [k for k in cols if k.find("E03_vac_fc") != -1][0]

        e03_temporal_fc    = [k for k in cols if k.find("E03_temporal_fc") != -1][0]
        e03_temporal_fdr    = [k for k in cols if k.find("E03_temporal_fdr") != -1][0]
    except IndexError:
        continue
    ggdf['E01'] = ggdf[e01_pfreq_col].apply(lambda x: x if x > 0 else 1E-6)
    ggdf['E02'] = ggdf[e02_pfreq_col].apply(lambda x: x if x > 0 else 1E-6)
    ggdf['E03'] = ggdf[e03_pfreq_col].apply(lambda x: x if x > 0 else 1E-6)
    i1 = ggdf[e03_vac_class] == "sig_expand"
    i2 = ggdf[e03_vac_fc] > 4
    i3 = ggdf[e03_vac_class] == "sig_contract"
    i4 = ggdf[e03_vac_fc] < .25
    ggdf13_expand = ggdf[i1&i2]
    ggdf13_contract = ggdf[i3&i4]
    ggdf13_tally = ggdf.groupby(['E01','E03']).count().reset_index(drop = False)
    
    #1,3, Check, for P673 
    #ggdf_expand = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_expand"').query('`15673_8_TCRB_E03_vac_fc` > 4')
    #ggdf_contract = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_contract"').query('`15673_8_TCRB_E03_vac_fc` < .25')
    
    
    """
    # ## FOR CREATING MANUSCRIPT FIGURE DATA FILES, 
    # First, Write out this files
    # ggdf13_tally.to_csv(f"ab_fig_REVIEW/{ptid}_ggdf13_tally.csv", index = False)
    # ggdf13_expand.to_csv(f"ab_fig_REVIEW/{ptid}_contract_ggdf13.csv", index = False)
    # ggdf13_contract.to_csv(f"ab_fig_REVIEW/{ptid}_expand_ggdf13.csv", index = False)
    # Second, open them in R, so that they can be written out 'figure_data_files/fig_1a_data.csv'
    # bind_rows(
    #   readr::read_csv('figure_data_files/15673_ggdf13_tally.csv') %>% 
    #     select(E01, E03, count=cdr3_b_aa)%>% 
    #     mutate(type = "all") %>%
    #     mutate(pubid = '673'),
    #   readr::read_csv('figure_data_files/15673_contract_ggdf13.csv')  %>% 
    #     select(E01, E03) %>% 
    #     mutate(count = 1) %>%
    #     mutate(type = "contract") %>% 
    #     mutate(pubid = '673') ,
    #   readr::read_csv('figure_data_files/15673_expand_ggdf13.csv')  %>% 
    #     select(E01, E03) %>% 
    #     mutate(count = 1) %>%
    #     mutate(type = "expand") %>% 
    #     mutate(pubid = '673') 
    # ) %>% 
    #   mutate(E01 = ifelse(E01 < 1E-6, 1E-6, E01)) %>% 
    #   mutate(E03 = ifelse(E03 < 1E-6, 1E-6, E03)) %>% 
    #   write.csv('figure_data_files/fig_1a_data.csv', row.names = F)
    #   #ggplot(aes(E01, E03)) + geom_point(aes(col = type)) + 
    #   #scale_y_log10()+ 
    #   #scale_x_log10()
    """
    
    i1 = ggdf[e02_vac_class] == "sig_expand"
    i2 = ggdf[e02_vac_fc] > 4
    i3 = ggdf[e02_vac_class] == "sig_contract"
    i4 = ggdf[e02_vac_fc] < .25
    ggdf12_expand = ggdf[i1&i2]
    ggdf12_contract = ggdf[i3&i4]
    #ggdf_expand = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_expand"').query('`15673_8_TCRB_E03_vac_fc` > 4')
    #ggdf_contract = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_contract"').query('`15673_8_TCRB_E03_vac_fc` < .25')
    ggdf12_tally = ggdf.groupby(['E01','E02']).count().reset_index(drop = False)

    # No used:
    i1 = ggdf[e03_temporal_fdr] < 0.05
    i2 = ggdf[e03_temporal_fc] > 4
    i3 = ggdf[e03_temporal_fdr] < 0.05
    i4 = ggdf[e03_temporal_fc] < 0.25
    ggdf23_expand = ggdf[i1&i2]
    ggdf23_contract = ggdf[i3&i4]
    #ggdf_expand = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_expand"').query('`15673_8_TCRB_E03_vac_fc` > 4')
    #ggdf_contract = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_contract"').query('`15673_8_TCRB_E03_vac_fc` < .25')
    ggdf23_tally = ggdf.groupby(['E02','E03']).count().reset_index(drop = False)

    # Plotnine (optional)
    # (ggplot(ggdf13_tally, aes('E01', 'E03' ) )
    # + geom_point(size = .05, color = "gray") 
    # + geom_point(data = ggdf13_expand, size = 1, shape = '2', color = 'red')
    # + geom_point(data = ggdf13_contract, size = 1, shape = '1',color = 'orange')
    # + theme_classic() 
    # + scale_y_log10()
    # + scale_x_log10() 
    # + geom_abline(linetype = "dashed" , color = "gray")
    # + annotation_logticks(sides = "lb", size = .2) 
    # + xlab("E01")
    # + ylab("E03")
    # + ggtitle(ptid[2:])
    # + coord_cartesian(xlim = (-6,-2), ylim = (-6,-2))
    # ).save(f'ab_fig/{ptid}_e01_e03.pdf')

    print(ptid)
    print("COMPLETED PLOTS")
