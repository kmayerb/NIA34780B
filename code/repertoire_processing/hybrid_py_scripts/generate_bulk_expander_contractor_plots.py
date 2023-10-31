"""
generate_bulk_expander_contractor_plots.py

Code to produce plot showing expander and contractor clones comparing two time points

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

do_these = ['15665','15518','15771','15545','15753','15668','15635','15706','15869','15784','15514','15582','15653','15742','15773','15531','15839','15577','15559','15837','15754','15527','15669','15782','15763','15758','15581','15836','15673','15684','15525','15761','15530']
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
    #1,3 
    #ggdf_expand = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_expand"').query('`15673_8_TCRB_E03_vac_fc` > 4')
    #ggdf_contract = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_contract"').query('`15673_8_TCRB_E03_vac_fc` < .25')
    ggdf13_tally = ggdf.groupby(['E01','E03']).count().reset_index(drop = False)

    i1 = ggdf[e02_vac_class] == "sig_expand"
    i2 = ggdf[e02_vac_fc] > 4
    i3 = ggdf[e02_vac_class] == "sig_contract"
    i4 = ggdf[e02_vac_fc] < .25
    ggdf12_expand = ggdf[i1&i2]
    ggdf12_contract = ggdf[i3&i4]
    #ggdf_expand = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_expand"').query('`15673_8_TCRB_E03_vac_fc` > 4')
    #ggdf_contract = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_contract"').query('`15673_8_TCRB_E03_vac_fc` < .25')
    ggdf12_tally = ggdf.groupby(['E01','E02']).count().reset_index(drop = False)


    i1 = ggdf[e03_temporal_fdr] < 0.05
    i2 = ggdf[e03_temporal_fc] > 4
    i3 = ggdf[e03_temporal_fdr] < 0.05
    i4 = ggdf[e03_temporal_fc] < 0.25
    ggdf23_expand = ggdf[i1&i2]
    ggdf23_contract = ggdf[i3&i4]
    #ggdf_expand = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_expand"').query('`15673_8_TCRB_E03_vac_fc` > 4')
    #ggdf_contract = ggdf.query('`15673_8_TCRB_E03_vac_class` == "sig_contract"').query('`15673_8_TCRB_E03_vac_fc` < .25')
    ggdf23_tally = ggdf.groupby(['E02','E03']).count().reset_index(drop = False)

    (ggplot(ggdf13_tally, aes('E01', 'E03' ) )
    + geom_point(size = .05, color = "gray") 
    + geom_point(data = ggdf13_expand, size = 1, shape = '2', color = 'red')
    + geom_point(data = ggdf13_contract, size = 1, shape = '1',color = 'orange')
    + theme_classic() 
    + scale_y_log10()
    + scale_x_log10() 
    + geom_abline(linetype = "dashed" , color = "gray")
    + annotation_logticks(sides = "lb", size = .2) 
    + xlab("E01")
    + ylab("E03")
    + ggtitle(ptid[2:])
    + coord_cartesian(xlim = (-6,-2), ylim = (-6,-2))
    ).save(f'ab_fig/{ptid}_e01_e03.pdf')

    (ggplot(ggdf12_tally, aes('E01', 'E02' ) )
    + geom_point(size = .05, color = "gray") 
    + theme_classic() 
    + scale_y_log10()
    + scale_x_log10() 
    + geom_point(data = ggdf12_expand, size = 1, shape = '2', color = 'red')
    + geom_point(data = ggdf12_contract, size = 1, shape = '1',color = 'orange')
    + geom_abline(linetype = "dashed" , color = "gray")
    + annotation_logticks(sides = "lb", size = .2) 
    + xlab("E01")
    + ylab("E02")
    + ggtitle(ptid[2:])
    + coord_cartesian(xlim = (-6,-2), ylim = (-6,-2))
    ).save(f'ab_fig/{ptid}_e01_e02.pdf')

    (ggplot(ggdf23_tally, aes('E02', 'E03' ) )
    + geom_point(size = .05, color = "gray") 
    + theme_classic() 
    + scale_y_log10()
    + scale_x_log10() 
    + geom_point(data = ggdf23_expand, size = 1, shape = '2', color = 'red')
    + geom_point(data = ggdf23_contract, size = 1, shape = '1',color = 'orange')
    + geom_abline(linetype = "dashed" , color = "gray")
    + annotation_logticks(sides = "lb", size = .2) 
    + xlab("E02")
    + ylab("E03")
    + ggtitle(ptid[2:])
    + coord_cartesian(xlim = (-6,-2), ylim = (-6,-2))
    ).save(f'ab_fig/{ptid}_e02_e03.pdf')

    print(ptid)
    print("COMPLETED PLOTS")

"""
Repeat for corner cases where E02 sample is not available
"""
plots = []
summary_data = []
long_data_frames= []
count = 0
get_fold_changes = True

do_these = ['15742','15758']
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

    ggdf = dftextra#[cols] #dftextra[['15673_6_TCRB_E01_pfreq', '15673_8_TCRB_E03_pfreq', '15673_8_TCRB_E03_vac_class','15673_8_TCRB_E03_vac_fc']]
    cols = dftextra.columns
    try:
        
        e01_pfreq_col = [k for k in cols if k.find("E01_pfreq") != -1][0]
        #e02_pfreq_col = [k for k in cols if k.find("E02_pfreq") != -1][0]
        e03_pfreq_col = [k for k in cols if k.find("E03_pfreq") != -1][0]
        #e02_vac_class = [k for k in cols if k.find("E02_vac_class") != -1][0]
        e03_vac_class = [k for k in cols if k.find("E03_vac_class") != -1][0]
        #e02_vac_fc    = [k for k in cols if k.find("E02_vac_fc") != -1][0]
        e03_vac_fc    = [k for k in cols if k.find("E03_vac_fc") != -1][0]

        e03_temporal_fc    = [k for k in cols if k.find("E03_temporal_fc") != -1][0]
        e03_temporal_fdr    = [k for k in cols if k.find("E03_temporal_fdr") != -1][0]
    except IndexError:
        continue
    ggdf['E01'] = ggdf[e01_pfreq_col].apply(lambda x: x if x > 0 else 1E-6)
    #ggdf['E02'] = ggdf[e02_pfreq_col].apply(lambda x: x if x > 0 else 1E-6)
    ggdf['E03'] = ggdf[e03_pfreq_col].apply(lambda x: x if x > 0 else 1E-6)
    i1 = ggdf[e03_vac_class] == "sig_expand"
    i2 = ggdf[e03_vac_fc] > 4
    i3 = ggdf[e03_vac_class] == "sig_contract"
    i4 = ggdf[e03_vac_fc] < .25
    ggdf13_expand = ggdf[i1&i2]
    ggdf13_contract = ggdf[i3&i4]
   
    ggdf13_tally = ggdf.groupby(['E01','E03']).count().reset_index(drop = False)

    i1 = ggdf[e03_temporal_fdr] < 0.05
    i2 = ggdf[e03_temporal_fc] > 4
    i3 = ggdf[e03_temporal_fdr] < 0.05
    i4 = ggdf[e03_temporal_fc] < 0.25

    (ggplot(ggdf13_tally, aes('E01', 'E03' ) )
    + geom_point(size = .05, color = "gray") 
    + geom_point(data = ggdf13_expand, size = 1, shape = '2', color = 'red')
    + geom_point(data = ggdf13_contract, size = 1, shape = '1',color = 'orange')
    + theme_classic() 
    + scale_y_log10()
    + scale_x_log10() 
    + geom_abline(linetype = "dashed" , color = "gray")
    + annotation_logticks(sides = "lb", size = .2) 
    + xlab("E01")
    + ylab("E03")
    + ggtitle(ptid[2:])
    + coord_cartesian(xlim = (-6,-2), ylim = (-6,-2))
    ).save(f'ab_fig/{ptid}_e01_e03.pdf')

print(ptid)
print("COMPLETED PLOTS")