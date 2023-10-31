# NI_A347808C_Figure_1.R 
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

filename_e03_sig = file.path(repo_loc, repo, 
                             'data/all_e03_sig.2.tsv')
filename_total_clones_per_sample = file.path(repo_loc, repo,
                                             'data/number_of_unique_productive_templates_per_sample.tsv')

# OUTPUTS
# Figure 1b-h
# Extended Data Figure 3: 
  # Trajectories of the vaccine expanded TRB clonotypes 
  # that significantly increased from E01 to E03 
global_font_size = 8
global_line_size = .2

# XNote: 3 participants excluded because of documented breakthrough infection 
breakthrough_ptid = c(15545,15669,15664)
e02_missing_ptid  = c(15742,15758)

fig_1b_filename    = file.path(repo_loc, repo, 'figures/Fig_1b_NI_A34780C.pdf')
fig_1c_filename    = file.path(repo_loc, repo, 'figures/Fig_1c_NI_A34780C.pdf')
fig_1c_p2_filename = file.path(repo_loc, repo, 'figures/Fig_1c_part2_NI_A34780C.pdf')
fig_1d_filename    = file.path(repo_loc, repo, 'figures/Fig_1d_NI_A34780C.pdf')#e
fig_1e_filename    = file.path(repo_loc, repo, 'figures/Fig_1e_NI_A34780C.pdf')#f
fig_1f_filename    = file.path(repo_loc, repo, 'figures/Fig_1f_NI_A34780C.pdf')#g
fig_1g_filename    = file.path(repo_loc, repo, 'figures/Fig_1g_NI_A34780C.pdf')#h
fig_1h_filename    = file.path(repo_loc, repo, 'figures/Fig_1h_NI_A34780C.pdf')#h
ex_fig_3a_filename = file.path(repo_loc, repo, 'figures/Ex_Fig_3a_NI_A34780C.pdf')
ex_fig_3b_filename = file.path(repo_loc, repo, 'figures/Ex_Fig_3b_NI_A34780C.pdf')

# OUTPUT FIGURE FILES
f1b_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1b_data.csv')
f1c_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1c_data.csv')
f1d_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1d_data.csv')
f1e_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1e_data.csv')
f1f_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1f_data.csv')
f1g_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1g_data.csv')
f1h_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_1h_data.csv')
ef3a_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_ex_fig_3a_data.csv')
ef3b_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_ex_fig_3b_data.csv')


# CODE
# LOAD DATA INCLUDING SIGNIFICANTLY VACCINE-EXPANDED CLONES E03 vs.E01
d = readr::read_tsv(filename_e03_sig)
# XNote: gg_581_837
# c(581, 837) representative participants
gg_fig_1b = d %>% 
  mutate(ptid = stringr::str_remove(ptid, pattern = "15")) %>% 
  filter(ptid %in% c(581, 837)) %>%  #representative participants
  select(ptid,  E03_vac_fdr,  E03_vac_fc, 
         predetected, E00_pfreq, 
         E00.5_pfreq, E01_pfreq, 
         E02_pfreq, E03_pfreq, E05_pfreq,       cdr3_b_aa,  v_b_gene,  j_b_gene,  cdr3_b_nucseq) %>% 
  filter(., E03_vac_fdr <0.05 & E03_vac_fc >4) %>%
  mutate(uid = seq_along(.$ptid)) %>%
  tidyr::gather(variable, value, -ptid, -uid, -E03_vac_fdr, -E03_vac_fc, -predetected,-cdr3_b_aa,  -v_b_gene,  -j_b_gene,  -cdr3_b_nucseq) %>%
  mutate(xpos= gsub(variable, pattern = "_pfreq", replacement ="")) %>% 
  mutate(value = ifelse(value < 1E-6, 1E-6, value)) %>% 
  filter(xpos != "E00.5") %>%
  filter(!(xpos == "E00" & predetected == "new")) %>% 
  filter(!ptid %in% c('15742','15758')) %>% 
  ggplot(aes(x =xpos, y = value, group = uid, col = predetected)) +
  geom_line(size = .2) + 
  scale_y_log10(breaks = c(10^-6,10^-5,10^-4,10^-3,10^-2),
                labels = c('ND',expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2)),
                #breaks = trans_breaks("log10", function(x) 10^x), 
                #labels = trans_format("log10", math_format(10^.x)),
                limits = c(1E-6, 5E-2)) + 
  theme_classic() + 
  facet_wrap(~ptid, scale= "free_y",ncol = 3) + 
  scale_color_manual(values = c("orange", "#00000030")) + 
  theme(legend.position = "none") + 
  geom_point(size = .05) + 
  theme(strip.background = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  annotation_logticks(side = "l" , size = global_line_size) + 
  theme(panel.spacing = unit(.65, "lines"))+
  ylab("expanded TRB productive frequency") 

# FIGURE 1B
pdf(fig_1b_filename, width = 4.25, height = 2)
gg_fig_1b 
dev.off()

# Write figure data
gg_fig_1b$data %>% 
  select(pubid = ptid, 
        visit = xpos,
         variable,
         value,
         e03_vac_fdr = E03_vac_fdr,
         e03_vac_fc = E03_vac_fc,
         group = predetected,
         uid,
         cdr3_b_aa,  v_b_gene,  j_b_gene,  cdr3_b_nucseq) %>% 
  write.csv(f1b_fdata,row.names = F)

# # EXTENDED FIGURE 3A, Relating to Figure 1 
# <gg_all> 
# <gg_ef3a> Extended Figure 3A
gg_efig_3a = d %>% 
  select(ptid,  E03_vac_fdr,  E03_vac_fc, predetected, 
         E00_pfreq, E00.5_pfreq, E01_pfreq, E02_pfreq, E03_pfreq,E05_pfreq,
         cdr3_b_aa,  v_b_gene,  j_b_gene,  cdr3_b_nucseq) %>% 
  filter(!ptid %in% e02_missing_ptid)%>%
  mutate(ptid = as.character(ptid) )%>%
  mutate(ptid = stringr::str_replace(ptid, pattern = "15",replacement="P")) %>%
  filter(., E03_vac_fdr <0.05 & E03_vac_fc >4) %>%
  mutate(uid = seq_along(.$ptid)) %>%
  tidyr::gather(variable, value, -ptid, -uid, -E03_vac_fdr, -E03_vac_fc, -predetected,-cdr3_b_aa,  -v_b_gene,  -j_b_gene,  -cdr3_b_nucseq) %>%
  mutate(xpos= gsub(variable, pattern = "_pfreq", replacement ="")) %>% 
  mutate(value = ifelse(value < 1E-6, 1E-6, value)) %>% 
  filter(xpos != "E00.5") %>%
  filter(!(xpos == "E00" & predetected == "new")) %>% 
  filter(!ptid %in% c('742','758')) %>% # Sample without E02 timepoint
  arrange(desc(predetected))%>%
  ggplot(aes(x =xpos, y = value, group = uid, col = predetected)) +
  geom_line(alpha = .4, size = .1) + 
  scale_y_log10(breaks = c(10^-6,10^-5,10^-4,10^-3,10^-2),
                labels = c('ND',expression(10^-5),expression(10^-4),expression(10^-3),expression(10^-2)),#trans_format("log10", math_format(10^.x)),
                limits = c(1E-6, 1E-2)) + 
  theme_classic() + 
  facet_wrap(~ptid,scales = "free", ncol = 4) +
  annotation_logticks(side = "l" , size = .1) + 
  scale_color_manual(values = c("orange", "black")) + 
  theme(legend.position = "none") + 
  geom_point(size = .05) + 
  theme(strip.background = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  annotation_logticks(side = "l" , size = global_line_size) + 
  theme(panel.spacing = unit(0, "lines"))+
  ylab("")+xlab("") 

# EXTENDED FIGURE 3A
# XNote: Original File Title: pdf('figures/full_cohort_novel_expanders_with_32.pdf', width = 7, height = 10)
pdf(ex_fig_3a_filename , width = 7, height = 10)
gg_efig_3a
dev.off()

gg_efig_3a$data %>% 
  select(pubid = ptid, 
         visit = xpos,
         variable,
         value,
         e03_vac_fdr = E03_vac_fdr,
         e03_vac_fc = E03_vac_fc,
         group = predetected,
         uid,
         cdr3_b_aa,  v_b_gene,  j_b_gene,  cdr3_b_nucseq) %>% 
  mutate(pubid = stringr::str_remove(pubid, pattern = "P"))%>%
  write.csv(ef3a_fdata, row.names = F)


# <gg_ef3b> EXTENDED FIGURE 3B
gg_efig_3b = d %>% 
  select(ptid,  E03_vac_fdr,  
         E03_vac_fc, predetected, 
         E00_pfreq, E00.5_pfreq, 
         E01_pfreq, E02_pfreq,
         E03_pfreq,E05_pfreq) %>% 
  filter(!ptid %in% e02_missing_ptid)%>% #pull(ptid) %>% unique()
  mutate(ptid = as.character(ptid) )%>%
  mutate(ptid = stringr::str_remove(ptid, pattern = "15")) %>% 
  filter(E01_pfreq == 0) %>%
  filter(E03_pfreq > 0) %>%
  filter(!is.na(E05_pfreq)) %>% 
  mutate(go_away_e05 = ifelse(E05_pfreq == 0, "Dissapear", "Persist")) %>% 
  group_by(ptid, go_away_e05, predetected) %>% 
  tally() %>% 
  tidyr::spread(key = go_away_e05, value = n) %>% 
  mutate(Dissapear = ifelse(is.na(Dissapear),0, Dissapear))%>%
  mutate(Persist = ifelse(is.na(Persist), 0, Persist))%>%
  mutate(percent_dis = 100-(100*Dissapear/(Dissapear +Persist))) %>% 
  select(ptid, predetected, percent_dis) %>% 
  ggplot(aes(x = predetected, y = percent_dis)) + 
  geom_boxplot(width = .25, aes(fill =predetected), outlier.size = 0, alpha = .5) + 
  geom_point(size = .5, col = "black") + 
  theme_classic() + 
  ylab("% Retained") + 
  xlab("")+
  scale_fill_manual(values = c("orange","black")) + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(angle= 90))+ 
  ggpubr::stat_compare_means(method = "wilcox", size = 2)

# EXTENDED FIGURE 3B
# XNote: Originally "figures/percent_e05_retained.pdf"
pdf(ex_fig_3b_filename, width = 1.5, height = 2)
gg_efig_3b + theme(axis.text.x = element_blank())
dev.off()


gg_efig_3b$data %>% 
  select(pubid = ptid, 
         group = predetected, 
         percent_e01_to_e03xp_at_e05 = percent_dis) %>% 
  write.csv(ef3b_fdata, row.names = F)
gg_efig_3b$data %>% group_by(predetected) %>% tally()
# FIGURE 1C
# "TRB productive frequency (sum)"
gg_fig_1c = d %>% 
  select(ptid,  E03_vac_fdr,  E03_vac_fc, predetected, 
         E00_pfreq, E00.5_pfreq, E01_pfreq, E02_pfreq, E03_pfreq,E05_pfreq) %>% 
  filter(., E03_vac_fdr <0.05 & E03_vac_fc >4) %>%
  mutate(uid = seq_along(.$ptid)) %>%
  tidyr::gather(variable, value, -ptid, -uid, -E03_vac_fdr, -E03_vac_fc, -predetected) %>%
  filter(!((ptid %in% breakthrough_ptid) & (variable == 'E05_pfreq'))) %>% #### <-
  mutate(xpos= gsub(variable, pattern = "_pfreq", replacement ="")) %>% 
  filter(xpos != "E00.5") %>%
  filter(!(xpos == "E00" & predetected == "new")) %>% 
  filter(!ptid %in% c('15742','15758')) %>% 
  mutate(value= ifelse(value ==1E-6, 0, value)) %>%
  group_by(ptid, xpos, predetected) %>% 
  summarise(n = sum(value > 1E-6 ), sum= sum(value)) %>%
  mutate(sum= ifelse(sum <1E-6, 1E-6, sum)) %>%
  mutate(group = paste0(ptid,predetected)) %>% 
  mutate(predetected2 = ifelse(predetected == "new", "Post-vaccination", "Post-infection")) %>%
  ggplot(aes(x =xpos, y = sum, group= group, col = predetected)) +
  geom_line(size = .3) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  facet_wrap(~predetected2, scale = "free_y") +
  theme_classic() + 
  scale_color_manual(values = c("#FFA50070", "#00000050")) + 
  theme(legend.position = "none") + 
  geom_point( pch = 20, size = .5) + 
  theme(strip.background = element_blank()) + 
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  annotation_logticks(side = "l" , size = global_line_size) + 
  theme(panel.spacing = unit(.65, "lines"))+
  ylab("TRB productive frequency (sum)")

# make a copy as ggpf
ggpf = gg_fig_1c

ggpf$data %>% group_by(predetected,xpos) %>% tally()

# FIGURE 1C
# XNote: gg_fig_1c = ggpf ('figures/fig_emily_1_productive_frequency.pdf')
pdf(fig_1c_filename, width = 4.25, height = 2)
gg_fig_1c
dev.off()

gg_fig_1c$data %>% 
  mutate(ptid = stringr::str_replace(ptid, pattern = "^15",replacement= ""))%>%
  select(pubid = ptid, 
         visit = xpos, 
         group = predetected2, 
         sum_freq = sum ) %>% 
  write.csv(f1c_fdata, row.names = F)

gg_fig_1c$data %>% 
  mutate(ptid = stringr::str_replace(ptid, pattern = "^15",replacement= ""))%>%
  select(pubid = ptid, 
         visit = xpos, 
         group = predetected2, 
         sum_freq = sum) %>%
  group_by(group, visit) %>% 
  filter(visit %in% c("E02","E03")) %>%
  summarise(max = 100*max(sum_freq, na.rm= T), min=100*min(sum_freq, na.rm= T) )

gg_fig_1c$data %>% 
  mutate(ptid = stringr::str_replace(ptid, pattern = "^15",replacement= ""))%>%
  select(pubid = ptid, 
         visit = xpos, 
         group = predetected2, 
         sum_freq = sum) %>%
  group_by(visit) %>% 
  filter(visit %in% c("E02","E03")) %>%
  summarise(max = 100*max(sum_freq, na.rm= T), min=100*min(sum_freq, na.rm= T) )


# FIGURE 1C part 2, make corresponding boxplot for the right margin
gg_fig_1c_part2 = gg_fig_1c$data %>% 
  filter(xpos %in% c("E02","E03","E05")) %>%
  ggplot(aes(x = xpos, y = sum)) + 
  geom_boxplot(outlier.shape = 21, outlier.size =.3, outlier.color = NA, aes(fill = predetected), size = .2, alpha = .5, width = .5) + 
  #geom_point( size = .1,pch = 20, aes(col = predetected), position = position_jitterdodge(jitter.width = .1, jitter.height = 0))+ 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + 
  scale_color_manual(values = c("orange", "black")) + 
  scale_fill_manual(values = c("orange", "black")) + 
  theme(legend.position = "none") + 
  theme(legend.position = "none") + 
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  annotation_logticks(side = "l" , size = global_line_size) + 
  theme(panel.spacing = unit(.1, "lines"))+
  ylab("TRB productive frequency (sum)")


gg_fig_1c_part2$data %>% 
  mutate(ptid = stringr::str_replace(ptid, pattern = "^15",replacement= ""))%>%
  select(pubid = ptid, 
         visit = xpos, 
         group = predetected2, 
         sum_freq = sum ) %>% 
  write.csv(f1d_fdata, row.names = F)
  


# FIGURE 1D
pdf(fig_1d_filename, width = 1.25, height = 1.75)
gg_fig_1c_part2 
dev.off()

# for testing we pull the data in wide form
x = gg_fig_1c$data %>% 
  filter(xpos %in% c("E01","E02","E03","E05")) %>% 
  select(ptid, xpos, predetected, sum) %>% 
  tidyr::spread(key = xpos, value = sum)
# xn - orange - "post-vaccination"
xn = x %>% filter(predetected == "new")
# xp - orange - "post-infection"
xp = x %>% filter(predetected != "new")
# "post-infection" comparisons between timepoints
wilcox.test(xp$E01, xp$E02, paired = T)
# 1.863e-09 ****
wilcox.test(xp$E02, xp$E03, paired = T)
# p-value = 1.863e-09 ****
wilcox.test(xp$E02, xp$E05, paired = T)
#p-value = 0.004339 ***
wilcox.test(xp$E03, xp$E05, paired = T)
# p-value = 0.06504 ns
# wilcox.test(xp$E02, xp$E05, paired = T)

# "post-vaccination" comparisons between timepoints
wilcox.test(xn$E02, xn$E03, paired = T)
#V = 0, p-value = 1.863e-09 ****
wilcox.test(xn$E02, xn$E05, paired = T)
#V = 36, p-value = 0.001118 ***
wilcox.test(xn$E03, xn$E05, paired = T)
#V = 237, p-value = 0.00166 ***

# Compare "post-infection" with "post-vaccination"
xnp = xn %>% left_join(xp, by = "ptid") 
wilcox.test(xnp$E02.x, xnp$E02.y, paired = T)
#p-value = 1.061e-05 ****
wilcox.test(xnp$E03.x, xnp$E03.y, paired = T)
#p-value = 3.79e-06 ****
wilcox.test(xnp$E05.x, xnp$E05.y, paired = T)
#p-value = 1.311e-05 ****
wilcox.test(xnp$E03.x, xnp$E05.y, paired = T)
# p-value = 0.04844 *

# WE LOAD FILE <breadth>
# num_unique_productive_templates used as denominator to compute breadth.
# breadth = unique templates / total unique templaates
breadth= readr::read_tsv(filename_total_clones_per_sample)%>% 
  mutate(xpos = visit)
#breadth= readr::read_tsv("number_of_unique_productive_templates_per_sample.tsv") %>% 
#  mutate(xpos = visit)
# First depicted as line plot but, not in figure
gg_fig_1f_lines = d %>% select(ptid,  
                   E03_vac_fdr,  E03_vac_fc, predetected, 
                   E00_pfreq, E00.5_pfreq, E01_pfreq, 
                   E02_pfreq, E03_pfreq, E05_pfreq) %>% 
  filter(., E03_vac_fdr <0.05 & E03_vac_fc >4) %>%
  mutate(uid = seq_along(.$ptid)) %>%
  tidyr::gather(variable, value, -ptid, -uid, -E03_vac_fdr, -E03_vac_fc, -predetected) %>%
  filter(!((ptid %in% breakthrough_ptid) & (variable == 'E05_pfreq'))) %>%
  mutate(xpos= gsub(variable, pattern = "_pfreq", replacement ="")) %>% 
  mutate(value = ifelse(value < 1E-6, 1E-6, value)) %>% 
  filter(xpos != "E00.5") %>%
  filter(!(xpos == "E00" & predetected == "new")) %>% 
  filter(!ptid %in% c('15742','15758')) %>% # No E02 timepoint
  group_by(ptid, xpos, predetected) %>% 
  summarise(n = sum(value > 1E-6 ), mean = mean(value), sum= sum(value)) %>%
  mutate(mean = ifelse(mean < 1E-6, 1E-6, mean)) %>%
  mutate(group = paste0(ptid,predetected)) %>% 
  left_join(breadth ,  by = c("ptid","xpos")) %>%
  mutate(breadth = n /num_unique_productive_templates) %>% 
  ggplot(aes(x =xpos, y = breadth, group= group, col = predetected)) +
  geom_line(alpha = .5, size = .51) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x))) + 
  facet_wrap(~predetected) +
  theme_classic() + 
  annotation_logticks(side = "l" , size = .1) + 
  scale_color_manual(values = c("orange", "#00000030")) + 
  theme(legend.position = "none") + 
  geom_point( size = .5,pch = 20)+ 
  theme(legend.position = "none") + 
  geom_point( pch = 20, size = .5) + 
  theme(strip.background = element_blank()) +
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  annotation_logticks(side = "l" , size = global_line_size) + 
  ylab("Clonal Breadth") + 
  xlab("")

# we reference ggb below, so make a copy
ggb = gg_fig_1f_lines

# FIGURE 1F
gg_fig_1f = gg_fig_1f_lines$data %>% 
  filter(xpos %in% c("E02","E03","E05")) %>%
  ggplot(aes(x = xpos, y = breadth)) + 
  geom_boxplot(outlier.shape =21, 
               outlier.color = NA,
               outlier.size = .5, 
               aes(fill = predetected), size = .2, alpha = .5, width = .5) + 
  #geom_point( size = .5,pch = 20, aes(col = predetected), position = position_jitterdodge(jitter.width = .1, jitter.height = 0))+ 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", math_format(10^.x))) + 
  theme_classic() + 
 
  scale_color_manual(values = c("orange", "black")) + 
  scale_fill_manual(values = c("orange", "black")) + 
  theme(legend.position = "none") + 
  theme(legend.position = "none") + 
  theme(strip.background = element_blank()) +
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  annotation_logticks(side = "l" , size = .2) + 
  #theme(axis.tick = element_line(size = global_line_size)) +
  ylab("E03 expanded TRB\n(clonal breadth)") + 
  xlab("")

# FIGURE 1F
#XNote: pdf("figures/fig_emily_1_boxplot_breadth.pdf", width = 1.5, height = 2)
pdf(fig_1f_filename, width = 1.75, height = 2)
gg_fig_1f
dev.off()


gg_fig_1f$data %>% 
  mutate(ptid = stringr::str_replace(ptid, pattern = "^15",replacement= ""))%>%
  select(pubid = ptid, 
         visit = xpos, 
         group = predetected, 
         breadth ) %>% 
  write.csv(f1f_fdata, row.names = F)
  

# Pull data for tests tests
x = gg_fig_1f_lines$data %>% 
  filter(xpos %in% c("E02","E03","E05")) %>% 
  select(ptid, xpos, predetected, breadth) %>% 
  tidyr::spread(key = xpos, value = breadth)
xn = x %>% filter(predetected == "new")
xp = x %>% filter(predetected != "new")
wilcox.test(xp$E02, xp$E03, paired = T)
# p-value = 0.3387 , ns
wilcox.test(xp$E02, xp$E05, paired = T)
#  = 214, p-value = 0.01956
wilcox.test(xp$E03, xp$E05, paired = T)
#V = 208, p-value = 0.03267
wilcox.test(xn$E02, xn$E03, paired = T)
#V = 131, p-value = 0.03643
wilcox.test(xn$E02, xn$E05, paired = T)
#V = 219, p-value = 0.01229
wilcox.test(xn$E03, xn$E05, paired = T)
#V = 231, p-value = 0.003453

xnp = left_join(xn, xp, by = "ptid")
wilcox.test(xnp$E02.x, xnp$E02.y, paired = T)
#p-value = 2.689e-05 ****
wilcox.test(xnp$E03.x, xnp$E03.y, paired = T)
#V = 93, p-value = 0.003223, **
wilcox.test(xnp$E05.x, xnp$E05.y, paired = T)
#V = 41, p-value = 0.002136 **

# FIGURE 1G
#my_comparisons <- list( c("E02", "E03"), c("E02", "E05"), c("E03", "E05") )
gg_fig_1g_box_percent_memory = ggb$data %>% 
  filter(xpos %in% c("E02","E03","E05")) %>%
  select(ptid, xpos, n, predetected) %>%
  tidyr::spread(key =predetected, value = n) %>% 
  mutate(percent_memory = previously_detected / (new + previously_detected)) %>%
  ggplot(aes(x = xpos, y = 100* percent_memory )) +
  geom_line(aes(group = ptid), size = .1, col = "gray")+
  geom_boxplot(outlier.size = NA, outlier.shape = NA, size = .2, alpha = .5, width = .5) + 
  geom_point( color = "black", size = .5,pch = 20, position = position_jitter(width = .01, height = 0))+
  theme_classic() + 
  theme(legend.position = "none") + 
  theme(legend.position = "none") + 
  theme(strip.background = element_blank()) +
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text = element_text(size  = global_font_size)) + 
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  ylab("expanded TRB\nderived from memory (%)") + 
  scale_color_gradient(low = "orange", high = "black")

# FIGURE 1G
#pdf("figures/fig_emily_1_boxplot_percent_memory.pdf", width = 1.5, height = 2)
pdf(fig_1g_filename,width = 1.5, height = 2)
gg_fig_1g_box_percent_memory 
dev.off()

gg_fig_1g_box_percent_memory$data %>% 
  mutate(ptid = stringr::str_replace(ptid, pattern = "^15",replacement= ""))%>%
  select(pubid = ptid, 
         visit = xpos, 
         post_vac = new, 
         post_inf = previously_detected, 
         percent_post_inf = percent_memory ) %>% 
  write.csv(f1g_fdata, row.names = F)

# Looking at percent memory changes at each time point
x = gg_fig_1f_lines$data %>% 
  filter(xpos %in% c("E02","E03","E05")) %>%
  select(ptid, xpos, n, predetected) %>%
  tidyr::spread(key =predetected, value = n) %>% 
  mutate(percent_memory = previously_detected / (new + previously_detected)) %>% 
  select(ptid, xpos, percent_memory) %>%
  tidyr::spread(key = xpos, value = percent_memory)
# Paired tests, of memory
wilcox.test(x$E02, x$E03, paired = T)
#p-value = 6.206e-06
wilcox.test(x$E02, x$E05, paired = T)
# p-value = 0.05625

# <PTID ORDER> DETERMINE THE PARTICIPANT ORDER BASED ON BREADTH OF EXPANDED CLONES AT E03
ggb$data %>%   filter(xpos=="E03") %>%
  select(ptid, n, predetected) %>% 
  tidyr::spread(key = predetected, value = n) %>% 
  mutate(p = new/previously_detected) %>% 
  arrange(desc(p)) %>% pull(ptid) %>% as.character()-> ptid_order

# substitute P for 15
ptid_order = stringr::str_replace(ptid_order, "15","P")

# FIGURE 1E
gg_fig_1e = ggpf$data %>% 
  filter(xpos %in% c("E03"))%>%
  mutate(ptid = stringr::str_replace(ptid, "15","P")) %>%
  ggplot(aes(x = forcats::fct_relevel(factor(ptid),ptid_order), y = sum)) + 
  geom_bar(aes(fill = predetected),stat = 'identity', position = "dodge", width = .5) + 
  theme_classic() + 
  scale_color_manual(values = c("orange", "black")) + 
  scale_fill_manual(values = c("orange", "black")) + 
  ylab("Expanded TRB repertoire fraction (%)") + 
  theme(legend.position = "none") + 
  theme(strip.background = element_blank()) +
  theme(strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.y = element_text(size  = global_font_size)) +
  theme(axis.text.x = element_text(angle = 90, size =global_font_size))+
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  coord_cartesian(ylim = c(0, 0.25))

#FIGURE 1E
pdf(fig_1e_filename, width = 4, height = 2)
gg_fig_1e
dev.off()

gg_fig_1e$data %>% 
  select(pubid = ptid, 
         visit = xpos, 
         group = predetected2, 
         sum_freq = sum) %>% 
  write.csv(f1e_fdata, row.names = F)

#FIGURE 1H
# <a> are barplots based on unique clones
gg_fig_1h = ggb$data %>% 
  filter(xpos %in% c('E02',"E03"))%>%
  mutate(ptid = stringr::str_replace(ptid, "15","P")) %>%
  ggplot(aes(x = forcats::fct_relevel(factor(ptid), ptid_order), y = n)) + 
  geom_bar(aes(fill = predetected),stat = 'identity', position="fill") + 
  facet_wrap(~xpos,scale = 'free',ncol = 2 )+
  theme_classic() + 
  scale_color_manual(values = c("orange", "black")) + 
  scale_fill_manual(values = c("orange", "black")) + 
  theme(axis.text.x = element_text(angle = 90, size = 5))+
  ylab("Unique expanded\nTRB (%)") + 
  xlab("")+
  theme(legend.position = "none") + 
  theme(strip.background = element_blank()) +
  theme(axis.title = element_text(size = 6)) +
  theme(axis.text = element_text(size = 6)) + 
  scale_y_continuous(expand =c(0,0), breaks = c(0,.25,.5,.75,1), labels = c("","25","50","75","100")) + 
  geom_text(data= ggb$data %>% filter(predetected == "new") %>%  
            mutate(ptid = stringr::str_replace(ptid, "15","P")) %>%
            filter(xpos %in% c('E02',"E03")) , 
            aes(y = .9, x = factor(ptid), label = n), size = 2.5, angle = 90, col = "black")+
  geom_text(data= ggb$data %>% filter(predetected == "previously_detected")%>% 
              mutate(ptid = stringr::str_replace(ptid, "15","P")) %>%
              filter(xpos %in% c('E02',"E03")) , 
            aes(y = .1, x = factor(ptid), label = n), size = 2.5, angle = 90, col = "white") + 
  theme(strip.background = element_blank()) +
  theme(strip.text = element_blank()) +
  theme(axis.title.y = element_text(size = global_font_size))+
  theme(axis.title.x = element_blank())+
  theme(axis.text.y = element_text(size  = global_font_size)) +
  theme(axis.text.x = element_text(angle = 90, size =global_font_size))+
  theme(strip.text = element_text(size = global_font_size)) +
  theme(axis.line = element_line(size = global_line_size)) +
  theme(panel.spacing = unit(2, "lines"))

#FIGURE 1H
pdf(fig_1h_filename, width = 8, height = 1.25)
gg_fig_1h 
dev.off()

gg_fig_1h$data %>% 
  mutate(group =  ifelse(predetected == "new", "Post-vaccination", "Post-infection")) %>%
  select(pubid = ptid, 
         visit = xpos, 
         group,
         n = n) %>% 
  write.csv(f1h_fdata, row.names = F)


# END FIGURE 1 #