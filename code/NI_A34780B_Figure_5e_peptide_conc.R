# Figure 5e.R
# NI_A347808B_Figure_1.R 
# FIGURE 1, AND ASSOCIATED EXTENDED FIGURE 3
# XNote: First Analysis: 2022_07_08_novel_expansion.R

# R DEPENDENCIES
require(ggplot2)
require(dplyr)
require(scales)

# INPUTS 
# set repo location as appropriate.
user = 'kmb'
if(user == 'esf') {
  repo_loc = '/Volumes/corey_l/esford3_kmayerbl_collab/software'
  repo = 'NIA34780B'
} else if(user == 'kmb') {
  repo_loc = '/Users/kmayerbl/active/david_koelle' ## if you're kosh
  repo = 'NI_A34780B'
} else {
  stop("set repo loc and repo manually")
}

# Inputs
peptide_dose_rsp_filename = file.path(repo_loc, repo, 
                             'data/Dose Response expt 2 Table CD8 origin TCRs.xls')


fig_5e_filename    = file.path(repo_loc, repo, 'figures/Fig_5e_NI_A34780C.pdf')
# OUTPUT FIGURE FILES
f5e_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_5e_data.csv')
# 2023-06-06-titration_plot.R
# concentration titration of KCY CD8 T cell

# Dependencies
require(ggplot2)
require(dplyr)
library(scales)
# Script
d = readxl::read_xls(peptide_dose_rsp_filename )
# Use more informative names for peptides.
swap = c('10mer' = 'S378-387\nKCYGVSPTKL',
         '9mer Left' = 'S378-386\nKCYGVSPTK',
         '9mer Right' = 'S379-387\nCYGVSPTKL',
         'No peptide' ="no peptide" )

swap2 = c('10mer' = 'KCYGVSPTKL',
         '9mer Left' = 'KCYGVSPTK',
         '9mer Right' = 'CYGVSPTKL',
         'No peptide' ="no peptide" )
# Use full description for names of the cloned TCRs.
swap_tcr = c('TCR1'= 'TCR1\nTRAV8-2-CVVSEKNTDKLIF-TRAJ34\nTRBV14-CASRRFGDTEAFF-TRBJ1-1',
             'TCR2'= 'TCR2\nTRAV13-1-CAARGVDAGGTSYGKLTF-TRAJ62\nTRBV19-CASSSIKASSYNEQFF-TRBJ2-1',
             'TCR3'= 'TCR3\nTRAV26-1-CIVTDNNAGNMLTF-TRAJ39\nTRBV9-CASSAWGGNPQHF-TRBJ1-5',
             'TCR4'= 'TCR4\nTRAV9-2-CALSDKKLTGGGNKLTF-TRAJ10\nTRBV20-1-CSARDLGGDTQYF-TRBJ2-3',
             'TCR8_1'='TCR8_1\nTRAV8-6-CAVSPRSNDYKLSF-TRAJ20\nTRBV20-1-CSARSWGSETQYF-TRBJ2-5',
             'TCR8_2'='TCR8_2\nTRAV9-2-CAVSARSNDYKLSF-TRAJ20\nTRBV3-1-CASRPLGEETQYF-TRBJ2-5')
# Create figure
gg_conc = 
  d %>% 
  select(tcr = TCR, 
         peptide = Peptide, 
         conc =Concentration, 
         a03 = `LCL (A*03:01+ or A*03:01-)`,
         percent = `Percent of CTV(+) cells that are mNeonGreen (+)`) %>% 
  mutate(peptide_name = swap[peptide]) %>%
  mutate(peptide_name2 = swap2[peptide]) %>%
  mutate(tcr_name = swap_tcr[tcr]) %>%
  filter(a03 == "A3+ LCL")%>%
  filter(!is.na(conc)) %>%
  mutate(conc = ifelse(conc == 0, 1E-7, conc)) %>%
  mutate(sym = ifelse(conc < 1E-6, "no peptide neg. control", "peptide")) %>%
  filter(stringr::str_starts(string = tcr, pattern = "TCR")) %>% 
  ggplot(aes(x = conc, y= percent, group = peptide_name, col = peptide_name)) + 
  geom_line()+
  geom_point() + 
  scale_x_continuous(trans = 'log10',
                     breaks = trans_breaks('log10', function(x) 10^x),
                     labels = trans_format('log10', math_format(10^.x)),
                     limits =c(1E-5,1) )+
  facet_wrap(~tcr,ncol = 3) +
  scale_color_manual("",values = c("orange","darkgray","darkblue")) + 
  xlab("peptide concetration micrograms per ml") + 
  ylab("% mNeonGreen positive") + 
  theme_bw() + 
  theme(panel.grid = element_blank())+
  theme(strip.background = element_blank()) + 
  annotation_logticks(side = "b") + 
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 0))
gg_conc
pdf(fig_5e_filename, width = 6, height = 4.5)
gg_conc
dev.off()

gg_conc$data %>% 
  select(tcr = tcr, 
         peptide= peptide_name2 , 
         stim = sym, 
         value = percent) %>% 
  mutate(variable = "percent_parent_mneongreen") %>% 
  write.csv(f5e_fdata, row.names= F)

