# NI_A347808B_Figure_2.R 
# FIGURE 2, AND ASSOCIATED EXTENDED FIGURES
# R DEPENDENCIES
require(ggplot2)
require(dplyr)
require(scales)

# INPUTS 
# Input data used for plotting comes from a Python script, with a copy
# 'ipython/2022_06_02_CD4_CD8_longitudinal_plots.ipynb'
# set repo location as appropriate.
# INPUTS
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

filename_f1 = file.path(
  repo_loc, repo, 
  'data/2022_06_23_all_long_data.tsv')

filename_f2 = file.path(
  repo_loc, repo,
  'data/2022_06_23_all_summary_data.tsv')

# OUTPUTS
# Figure 2c, 2d, Extended Data Figure 3: 
# Trajectories of the vaccine expanded TRB clonotypes 
# that significantly increased from E01 to E03 
cd_colors=c("CD4"= "#009966","CD8"="#333399")
#CMYK, green and blue
global_font_size = 8
global_line_size = .2

fig_2c_cd4_filename = file.path(repo_loc,'NI_A34780B/figures/Fig_2c_cd4_NI_A34780C.pdf') #
fig_2c_cd8_filename = file.path(repo_loc,'NI_A34780B/figures/Fig_2c_cd8_NI_A34780C.pdf') #
fig_2d_cd4_filename = file.path(repo_loc,'NI_A34780B/figures/Fig_2d_cd4_NI_A34780C.pdf') #
fig_2d_cd8_filename = file.path(repo_loc,'NI_A34780B/figures/Fig_2d_cd8_NI_A34780C.pdf') #
ex_fig_6a_filename = file.path(repo_loc,'NI_A34780B/figures/Ex_Fig_6a_NI_A34780C.pdf')

f2c_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_2c_data.csv')
f2d_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_2d_data.csv')
ef6_fdata= file.path(repo_loc,repo,'figure_data_files/NI_A34780C_ex_fig_6_data.csv')
#f2g_fdata=file.path(repo_loc,repo, 'figure_data_files/fig_2g_data.csv')

# <d> for data
# read in the data, exclude those cells which had too few CD4 or CD8
# UMIs to be identified as either cell type 
d = readr::read_tsv(filename_f1) %>% 
  filter(cell_type != "LOW_CD") 

ds = readr::read_tsv(filename_f2) %>% 
  filter(cell_type != "LOW_CD") %>%
  filter(summary == "mean")  %>% 
  tidyr::separate(variable, sep = "_", into = c("ptid",'v','chain','visit','pfreq')) %>% 
  mutate(ptid = as.numeric(ptid))

# rank the ptids
rank_order_ptid_cd8_breadth_E03 = d %>% group_by(x_pos, ptid, cell_type) %>% 
  tally() %>%
  filter(x_pos == 'E03', cell_type == "CD8") %>% 
  arrange(desc(n)) %>% pull(ptid) 

# Count number of clones per PTID at each time (x_pos)
n_per_ptid_xpos = d %>% filter(value > 0) %>% group_by(x_pos, ptid, cell_type) %>% tally()  %>% 
  mutate(ptid = factor(ptid, levels = rank_order_ptid_cd8_breadth_E03 ))

# Compute number of clones at E03 that were not found 
# at any of the pre-vaccine visit
E03_denovo_per_ptid = d %>% select(plot_clone_id, cell_type, x_pos, value, ptid) %>% 
  group_by(plot_clone_id, ptid, x_pos,cell_type) %>% 
  slice(1) %>%
  tidyr::spread(key = x_pos, value = value) %>% 
  replace(is.na(.), 0) %>% 
  mutate(E03_denovo_count = ifelse(E03 > 0 & (E00 == 0 & `E00.5`== 0 & E01 == 0), T,F)) %>% 
  group_by(cell_type, ptid) %>% 
  summarise(n = sum(E03 > 0), n_E03_denovo = sum(E03_denovo_count)) %>% 
  mutate(percent_E03_denovo = n_E03_denovo/n) %>%
  mutate(x_pos= "E03") %>%
  mutate(ptid = factor(ptid, levels = rank_order_ptid_cd8_breadth_E03 )) %>% 
  mutate(percent_E03_denovo = round(percent_E03_denovo, 2)) %>%
  mutate(percent_E03_denovo = paste0(100*round(percent_E03_denovo, 2), "%"))

# Compute number of clones at E03 that were not found 
# at E02
E03_new_after_E02_ptid = d %>% select(plot_clone_id, cell_type, x_pos, value, ptid) %>% 
  #filter(cell_type == "CD8") %>%
  group_by(plot_clone_id, ptid, x_pos,cell_type) %>% 
  slice(1) %>%
  tidyr::spread(key = x_pos, value = value) %>% 
  replace(is.na(.), 0) %>% 
  mutate(E03_new_E02_count = ifelse(E03 > 0 & (E02 == 0), T,F)) %>% 
  group_by(cell_type, ptid) %>% 
  summarise(n = sum(E03 > 0) , n_E03_new_E02 = sum(E03_new_E02_count)) %>% 
  mutate(percent_E03_new_E02 = n_E03_new_E02/n) %>%
  mutate(x_pos= "E03") %>%
  mutate(ptid = factor(ptid, levels = rank_order_ptid_cd8_breadth_E03 )) %>% 
  mutate(percent_E03_new_E02 = round(percent_E03_new_E02, 2)) %>%
  mutate(percent_E03_new_E02= paste0(100*round(percent_E03_new_E02, 2), "%"))

# Compute number of clones at E02 that 
# were not found at any of the pre-vaccinte timepoints
E02_denovo_per_ptid = d %>% select(plot_clone_id, cell_type, x_pos, value, ptid) %>% 
  #filter(cell_type == "CD8") %>%
  group_by(plot_clone_id, ptid, x_pos,cell_type) %>% 
  slice(1) %>%
  tidyr::spread(key = x_pos, value = value) %>% 
  replace(is.na(.), 0) %>% 
  mutate(E02_denovo_count = ifelse(E02 > 0 & (E00 == 0 & `E00.5`== 0 & E01 == 0), T,F)) %>% 
  group_by(cell_type, ptid) %>% 
  summarise(n = sum(E02 > 0) , n_E02_denovo = sum(E02_denovo_count)) %>% 
  mutate(percent_E02_denovo = n_E02_denovo/n) %>%
  mutate(x_pos= "E02") %>%
  mutate(ptid = factor(ptid, levels = rank_order_ptid_cd8_breadth_E03 )) %>% 
  mutate(percent_E02_denovo = paste0(100*round(percent_E02_denovo, 2), "%"))

E02_denovo_per_ptid  %>% write.csv()
E03_denovo_per_ptid  %>% write.csv()
E03_new_after_E02_ptid %>% write.csv()

# create dictionary from long to short ptid
d$ptid_short = d$ptid %>% sub(pattern = "^15", replace = "P")
lab_over = unique(d$ptid)
new <- unique(d$ptid) %>% sub(pattern = "^15", replace = "P")
names(new) = unique(d$ptid)



# Only consider CD4 and CD8 T cells
gg_all_cd4 =d %>% filter(cell_type == "CD4") %>% 
  ggplot(., aes(x = x_pos, y = value, col = cell_type)) + 
  geom_point(size = .05, alpha = .4) + 
  geom_line(aes(group = plot_clone_id), size = .1, alpha = .2) + 
  #geom_line(aes(group = ptid), data = ds, size = .75, alpha = .9)+
  facet_grid(ptid~cell_type, labeller = labeller(ptid = new)) + 
  scale_y_log10(expand = c(.05,0), 
                breaks = c(10^-6,10^-5,10^-4,10^-3,10^-2, 10^1),
                labels = c('ND',expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)),
                limits = c(1E-6,1))+
  scale_color_manual(values = cd_colors) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = global_font_size)) + 
  theme(axis.text.x  = element_text(size = global_font_size, angle = 90)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size)) + 
  theme(strip.text.y = element_text(size = global_font_size, angle = 0))+
  #theme(axis.line.y = element_blank()) + 
  theme(strip.text = element_text(size = 6)) + 
  theme(strip.background = element_blank()) + 
  ylab("CD4+ TS TRB\n(frequency)") + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2) + 
  geom_text(data = n_per_ptid_xpos %>% filter(cell_type== "CD4") , 
            aes(x = x_pos, y = .5, label = n), size = 1.5, angle = 0, 
            nudge_x = .4, col = "black") + 
  geom_text(data = E03_denovo_per_ptid %>% filter(cell_type== "CD4"), 
            aes(x = x_pos, y = .05, label = percent_E03_denovo), 
            size = 1.5, angle = 0, nudge_x = .5, col ="black") + 
  geom_text(data = E02_denovo_per_ptid %>% filter(cell_type== "CD4"), 
            aes(x = x_pos, y = .05, label = percent_E02_denovo), 
            size = 1.5, angle = 0, nudge_x = .5, col ="black")+
  geom_text(data = E03_new_after_E02_ptid %>% filter(cell_type== "CD4"),
            aes(x = x_pos, y = .007, label = percent_E03_new_E02), 
            size = 1.5, angle = 0, nudge_x = .5, col ="#666666")+
  theme(axis.ticks.y = element_blank())  

gg_all_cd8 =d %>% filter(cell_type == "CD8") %>% 
  ggplot(., aes(x = x_pos, y = value, col = cell_type)) + 
  geom_point(size = .05, alpha = .4) + 
  geom_line(aes(group = plot_clone_id), size = .1, alpha = .2) + 
  facet_grid(ptid~cell_type, labeller = labeller(ptid = new)) + 
  scale_y_log10(expand = c(.05,0), 
                breaks = c(10^-6,
                           10^-5,
                           10^-4,
                           10^-3,
                           10^-2, 
                           10^1),
                labels = c('ND',
                           expression(10^-5),
                           expression(10^-4),
                           expression(10^-3),
                           expression(10^-2),
                           expression(10^-1)),
                limits = c(1E-6,1))+
  scale_color_manual(values = cd_colors) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = global_font_size)) + 
  theme(axis.text.x  = element_text(size = global_font_size, angle = 90)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size)) + 
  theme(strip.text.y = element_text(size = global_font_size, angle = 0))+
  #theme(axis.line.y = element_blank()) + 
  theme(strip.text = element_text(size = 6)) + 
  theme(strip.background = element_blank()) + 
  ylab("CD8+ TS TRB\n(frequency)") + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2) + 
  geom_text(data = n_per_ptid_xpos %>% filter(cell_type== "CD8") , 
            aes(x = x_pos, y = .5, label = n), size = 1.5, angle = 0, 
            nudge_x = .4, col = "black") + 
  geom_text(data = E03_denovo_per_ptid %>% filter(cell_type== "CD8"), 
            aes(x = x_pos, y = .05, label = percent_E03_denovo), 
            size = 1.5, angle = 0, nudge_x = .5, col ="black") + 
  geom_text(data = E02_denovo_per_ptid %>% filter(cell_type== "CD8"), 
            aes(x = x_pos, y = .05, label = percent_E02_denovo), 
            size = 1.5, angle = 0, nudge_x = .5, col ="black")+
  geom_text(data = E03_new_after_E02_ptid %>% filter(cell_type== "CD8"),
            aes(x = x_pos, y = .007, label = percent_E03_new_E02), 
            size = 1.5, angle = 0, nudge_x = .5, col ="#666666")+
  theme(axis.ticks.y = element_blank())  

# This is the plot that is used at Extended Data Fig 6, 
# 10^-6 changed to ND
pdf(ex_fig_6a_filename, height = 10, width = 5.5)
gridExtra::grid.arrange(gg_all_cd4, gg_all_cd8, ncol = 2) 
dev.off()
# Write data
bind_rows(gg_all_cd4$data %>% 
  select(pubid = ptid_short, 
         cell_type,
         visit = x_pos, 
         value, 
         cdr3_b_nucseq,
         cdr3_b_aa,
         v_b_gene,
         j_b_gene),
gg_all_cd8$data %>% 
  select(pubid = ptid_short, 
         cell_type,
         visit = x_pos, 
         value, 
         cdr3_b_nucseq,
         cdr3_b_aa,
         v_b_gene,
         j_b_gene)) %>%
  mutate(pubid = stringr::str_remove(pubid, "P")) %>%
  write.csv(ef6_fdata, row.names = F)
         

# Only consider CD4 and CD8 T cells
gg_fig_2c_cd4 =d %>% filter(cell_type == "CD4") %>% 
  filter(ptid %in% c(15673, 15836)) %>%
  ggplot(., aes(x = x_pos, y = value, col = cell_type)) + 
  geom_point(size = .05, alpha = .4) + 
  geom_line(aes(group = plot_clone_id), size = .1, alpha = .2) + 
  #geom_line(aes(group = ptid), data = ds, size = .75, alpha = .9)+
  facet_grid(~ptid, labeller = labeller(ptid = new)) + 
  scale_y_log10(expand = c(.05,0), 
                breaks = c(10^-6,10^-5,10^-4,10^-3,10^-2, 10^-1),
                labels = c('ND',expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)),
                limits = c(1E-6,1))+
  scale_color_manual(values = cd_colors) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = global_font_size)) + 
  theme(axis.text.x  = element_text(size = global_font_size, angle = 90)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size)) + 
  theme(strip.text.y = element_text(size = global_font_size, angle = 0))+
  #theme(axis.line.y = element_blank()) + 
  theme(strip.text = element_text(size = 6)) + 
  theme(strip.background = element_blank()) + 
  ylab("CD4+ TS TRB\n(frequency)") + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2) + 
  geom_text(data = n_per_ptid_xpos %>% filter(cell_type== "CD4") %>% filter(ptid %in% c(15673, 15836)), 
            aes(x = x_pos, y = .5, label = n), size = 2, angle = 0, 
            nudge_x = .4, col = "black") + 
  geom_text(data = E03_denovo_per_ptid %>% filter(cell_type== "CD4") %>% filter(ptid %in% c(15673, 15836)),  
            aes(x = x_pos, y = .05, label = percent_E03_denovo), 
            size = 2, angle = 0, nudge_x = .5, col ="black") + 
  geom_text(data = E02_denovo_per_ptid %>% filter(cell_type== "CD4") %>% filter(ptid %in% c(15673, 15836)),  
            aes(x = x_pos, y = .05, label = percent_E02_denovo), 
            size = 2, angle = 0, nudge_x = .5, col ="black")+
  geom_text(data = E03_new_after_E02_ptid %>% filter(cell_type== "CD4") %>% filter(ptid %in% c(15673, 15836)), 
            aes(x = x_pos, y = .007, label = percent_E03_new_E02), 
            size = 2, angle = 0, nudge_x = .5, col ="#666666")+
  theme(axis.ticks.y = element_blank())  

gg_fig_2c_cd8 =d %>% filter(cell_type == "CD8") %>% 
  filter(ptid %in% c(15673, 15836)) %>%
  ggplot(., aes(x = x_pos, y = value, col = cell_type)) + 
  geom_point(size = .05, alpha = .4) + 
  geom_line(aes(group = plot_clone_id), size = .1, alpha = .2) + 
  #geom_line(aes(group = ptid), data = ds, size = .75, alpha = .9)+
  facet_grid(~ptid, labeller = labeller(ptid = new)) + 
  scale_y_log10(expand = c(.05,0), 
                breaks = c(10^-6,10^-5,10^-4,10^-3,10^-2, 10^-1),
                labels = c('ND',expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1)),
                limits = c(1E-6,1))+
  scale_color_manual(values = cd_colors) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = global_font_size)) + 
  theme(axis.text.x  = element_text(size = global_font_size, angle = 90)) + 
  theme(axis.title.x = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size)) + 
  theme(strip.text.y = element_text(size = global_font_size, angle = 0))+
  #theme(axis.line.y = element_blank()) + 
  theme(strip.text = element_text(size = 6)) + 
  theme(strip.background = element_blank()) + 
  ylab("CD8+ TS TRB\n(frequency)") + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2) + 
  geom_text(data = n_per_ptid_xpos %>% filter(cell_type== "CD8") %>% filter(ptid %in% c(15673, 15836)), 
            aes(x = x_pos, y = .5, label = n), size = 2, angle = 0, 
            nudge_x = .4, col = "black") + 
  geom_text(data = E03_denovo_per_ptid %>% filter(cell_type== "CD8") %>% filter(ptid %in% c(15673, 15836)),  
            aes(x = x_pos, y = .05, label = percent_E03_denovo), 
            size =2, angle = 0, nudge_x = .5, col ="black") + 
  geom_text(data = E02_denovo_per_ptid %>% filter(cell_type== "CD8") %>% filter(ptid %in% c(15673, 15836)),  
            aes(x = x_pos, y = .05, label = percent_E02_denovo), 
            size = 2, angle = 0, nudge_x = .5, col ="black")+
  geom_text(data = E03_new_after_E02_ptid %>% filter(cell_type== "CD8") %>% filter(ptid %in% c(15673, 15836)), 
            aes(x = x_pos, y = .007, label = percent_E03_new_E02), 
            size = 2, angle = 0, nudge_x = .5, col ="#666666")+
  theme(axis.ticks.y = element_blank())  


# FIGURE 2C
pdf(fig_2c_cd4_filename, height = 1.5, width = 5)
gg_fig_2c_cd4
dev.off()

pdf(fig_2c_cd8_filename, height = 1.5, width = 5)
gg_fig_2c_cd8
dev.off()

bind_rows(gg_fig_2c_cd4$data %>% 
            select(pubid = ptid_short, 
                   cell_type,
                   visit = x_pos, 
                   value, 
                   cdr3_b_nucseq,
                   cdr3_b_aa,
                   v_b_gene,
                   j_b_gene),
          gg_fig_2c_cd8$data %>% 
            select(pubid = ptid_short, 
                   cell_type,
                   visit = x_pos, 
                   value, 
                   cdr3_b_nucseq,
                   cdr3_b_aa,
                   v_b_gene,
                   j_b_gene)) %>%
  mutate(pubid = stringr::str_remove(pubid, "P")) %>%
  write.csv(f2c_fdata, row.names = F)


# FIGURE 2D
gg_mean_cd4 = ds %>%
  filter(cell_type == "CD4") %>% 
  filter(x_pos != "E04") %>% 
  ggplot(. , aes(x = x_pos, y = value, col = cell_type)) + 
  geom_line(aes(group = ptid),size = .4, alpha = .5)+
  scale_y_log10(expand = c(.05,.15), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1E-6, 1E-3))+
  scale_color_manual(values =cd_colors) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = global_font_size)) + 
  theme(axis.text.x  = element_text(size = global_font_size, angle = 90)) + 
  theme(axis.title.y = element_text(size = global_font_size)) + 
  theme(axis.title.x = element_blank()) + 
  theme(strip.text.y = element_text(size = global_font_size, angle = 0))+
  theme(strip.text = element_text(size = global_font_size)) + 
  theme(strip.background = element_blank()) + 
  #theme(axis.text = element_line(size = global_line_size))+
  ylab(bquote('CD4+ TS TRB\n(mean frequency)')) + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2)+
  theme(axis.ticks.y = element_blank())


gg_mean_cd8 = ds %>%
  filter(cell_type == "CD8") %>% 
  filter(x_pos != "E04") %>% 
  ggplot(. , aes(x = x_pos, y = value, col = cell_type)) + 
  geom_line(aes(group = ptid),size = .4, alpha = .5)+
  scale_y_log10(expand = c(.05,.15), 
                breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1E-6, 1E-3))+
  scale_color_manual(values =cd_colors) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = global_font_size)) + 
  theme(axis.text.x  = element_text(size = global_font_size, angle = 90)) + 
  theme(axis.title.y = element_text(size = global_font_size)) + 
  theme(axis.title.x = element_blank()) + 
  theme(strip.text.y = element_text(size = global_font_size, angle = 0))+
  theme(strip.text = element_text(size = global_font_size)) + 
  theme(strip.background = element_blank()) + 
  #theme(axis.text = element_line(size = global_line_size))+
  ylab(bquote('CD8+ TS TRB\n(mean frequency)')) + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2)+
  theme(axis.ticks.y = element_blank())

# FIGURE 2D CD4+ Ts (mean frequency)
pdf(fig_2d_cd4_filename, height = 1.5, width = 2.25)
gg_mean_cd4
dev.off()
# FIGURE 2D CD8+ Ts (mean frequency)
pdf(fig_2d_cd8_filename , height = 1.5, width =2.25)
gg_mean_cd8
dev.off()

bind_rows(gg_mean_cd4$data %>% 
            select(pubid = ptid, 
                   visit, 
                   value,
                   variable = summary) %>% 
            mutate(cell_type = "CD4"),
          gg_mean_cd8$data %>% 
            select(pubid = ptid, 
                   visit, 
                   value,
                   variable = summary) %>% 
            mutate(cell_type = "CD8")) %>%
  mutate(pubid = stringr::str_replace(pubid, pattern ="^15",replacement="")) %>%
  write.csv(f2d_fdata, row.names = F)
