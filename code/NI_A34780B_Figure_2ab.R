# NI_A34780B_Figure_2a.R # 

# DEPENDENCIES
require(dplyr)
require(ggplot2)
require(scales)
require(stringr)

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

# Expanded was based on both FDR < 0.05 and log[2](FC) > 2

filename_f1 = file.path(
  repo_loc, repo,
  '/data/2022_06_23_df_expand.tsv')
filename_f2 = file.path(
  repo_loc, repo,
  '/data/2022_06_23_df_contract.tsv')
filename_f3 = file.path(
  repo_loc, repo,
  'data/2022_06_23_all_long_data.tsv')

# PARAMETERS
fc = 4                        # FC = 4, Log2 Fold Change > 2
my_ptids = c('15673','15836') # Example Participants
global_font_size = 8

# OUTPUTS
global_font_size = 8
figure_2a_filename = file.path(repo_loc, repo, 'figures/Fig_2a_NI_A34780C.pdf')
figure_2b_filename = file.path(repo_loc, repo, 'figures/Fig_2b_NI_A34780C.pdf')
ex_fig_4_filename = file.path(repo_loc, repo, 'figures/Ex_Fig_4_NI_A34780C.pdf')


f2a_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_2a_data.csv')
f2b_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_2b_data.csv')
ef4a_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_ex_fig_4a_data.csv')



# PLOT SCRIPT

# Expanded Clones
# Expanded/Contracted clones include bulk TRB clones, 
# only some match TRB in AIM+ TCRab Ts 
df = readr::read_tsv(filename_f1) %>% 
  mutate(key_b = paste0(v_b_gene, "+", cdr3_b_aa, "+", j_b_gene))

# Contracted Clones
dfc = readr::read_tsv(filename_f2) %>% 
  mutate(key_b = paste0(v_b_gene, "+", cdr3_b_aa, "+", j_b_gene))

# Read in all longitudinal clones, remove those those that expanded or contracted
sig_change_keys = c(df %>% filter(E03_vac_fc > fc) %>% pull(key_b), 
                    dfc %>% filter(E03_vac_fc > 1/fc) %>% pull(key_b))
d = readr::read_tsv(filename_f3) %>% 
  filter(cell_type != "LOW_CD") %>% 
  mutate(key_b = paste0(v_b_gene, "+", cdr3_b_aa, "+", j_b_gene))
d = d[! d$key_b %in% sig_change_keys ,]


# Make a boundary diagram, approximate only
x = c(seq(-7, -3, by = .1) ,seq(from = -4, to = -1 , by = 1))
x =10^x
ypos = 2E-5 + 4*x
yneg = .25*x - 5E-6
boundary = data.frame(x,ypos, yneg)
boundary$yneg

# Count number of expanded AIM+ clones by person
df_expansions = df %>% 
  filter(E03_vac_fc > fc) %>%
  group_by(ptid, longitudinal_id) %>% 
  slice(1) %>%
  mutate(E01_pfreq = ifelse(E01_pfreq < 1E-6, 1E-6,  E01_pfreq )) %>% 
  mutate(E03_pfreq = ifelse( E03_pfreq < 1E-6, 1E-6,  E03_pfreq )) %>% 
  mutate(ptid2 = stringr::str_replace(ptid, "15","P"))
# Count number of contracted AIM+ clones by person
df_contractions = dfc %>% 
  filter(E03_vac_fc < 1/fc) %>%
  group_by(ptid, longitudinal_id) %>% 
  slice(1) %>%
  mutate(E01_pfreq = ifelse(E01_pfreq < 1E-6, 1E-6,  E01_pfreq )) %>% 
  mutate(E03_pfreq = ifelse( E03_pfreq < 1E-6, 1E-6,  E03_pfreq )) %>% 
  mutate(ptid2 = stringr::str_replace(ptid, "15","P"))

# Tablulate number of expanded clones
number_of_expanded = df %>% 
  group_by(ptid, longitudinal_id) %>% 
  slice(1) %>% 
  filter(E03_vac_fc > fc) %>%
  group_by(ptid) %>% 
  tally() %>% 
  arrange(desc(n))
# check for all 12 ptids
stopifnot(length(number_of_expanded$ptid)==12)
d_non_change = d %>% filter(x_pos %in% c("E01","E03")) %>% 
  filter(value > 0) %>%
  group_by(ptid, plot_clone_id, cell_type, x_pos,  cdr3_b_nucseq, cdr3_b_aa,v_b_gene, j_b_gene ) %>% 
  slice(1) %>% 
  ungroup()%>%
  select(ptid ,cell_type, value, x_pos,  plot_clone_id,cdr3_b_nucseq, cdr3_b_aa,v_b_gene, j_b_gene) %>% 
  tidyr::spread(key = x_pos, value = value) %>% 
  mutate(E01_pfreq = ifelse(is.na(E01), 1E-6, E01 )) %>%
  mutate(E03_pfreq = ifelse(is.na(E03), 1E-6, E03)) %>% 
  filter(ptid %in% number_of_expanded$ptid) %>%
  mutate(check_fc = E03_pfreq/E01_pfreq) %>% 
  filter(check_fc < 100) %>%
  mutate(ptid2 = stringr::str_replace(ptid, "15","P"))



# Number of expanded clones per person
df_exp = df %>% 
  mutate(cell_type2 = ifelse(is.na(cell_type), "UNK", cell_type)) %>% 
  group_by(ptid,cell_type2) %>% tally() %>% 
  ungroup()%>%
  select(ptid, cell_type2, n) %>%
  tidyr::spread(key = cell_type2, value = n) %>%
  mutate(ptid2 = stringr::str_replace(ptid, "15","P"))

# Number of contraced clones per person
df_con = dfc %>% 
  mutate(cell_type2 = ifelse(is.na(cell_type), "UNK", cell_type)) %>% 
  group_by(ptid,cell_type2) %>% tally() %>% 
  ungroup()%>%
  select(ptid, cell_type2, n) %>%
  tidyr::spread(key = cell_type2, value = n) %>%
  mutate(ptid2 = stringr::str_replace(ptid, "15","P"))

d_n_non_expand = d_non_change %>%
  group_by(ptid, cell_type) %>% 
  tally() %>% 
  tidyr::spread(key = cell_type, value = n) %>%
  mutate(ptid2 = stringr::str_replace(ptid, "15","P"))

# Figure 2a 
## SELECT ONLY TWO PTIDS
ggx_2ptid_selected = ggplot(df_expansions %>% filter(ptid %in% my_ptids) , aes(x = E01_pfreq, y = E03_pfreq)) + 
  geom_point(data = df_expansions %>% filter(ptid %in% my_ptids),
             aes(col = cell_type), pch = 2, size = .6, alpha = .5)+
  geom_point(data = df_contractions %>% filter(ptid %in% my_ptids), aes(col = cell_type), size = .6, pch = 6, alpha = .5)+
  geom_point(data = d_non_change %>% filter(ptid %in% my_ptids), aes(x = E01_pfreq, y = E03_pfreq, col = cell_type), size =.6, alpha = .5, pch = 1) +
  scale_x_log10(expand = c(.05,.05), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(expand = c(.05,.15), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  theme_classic() + 
  theme(axis.text.x = element_text(size = 6)) + 
  theme(axis.text.y = element_text(size = 6)) +
  theme(axis.title.x = element_text(size = 6)) + 
  theme(axis.title.y = element_text(size = 6)) +
  annotation_logticks(sides = "lb", size = .2, 
                      short = unit(0.005, "npc"),
                      mid = unit(0.0075, "npc"),
                      long = unit(0.01, "npc"),) + 
  coord_cartesian(xlim = c(1E-6, .1), ylim = c(1E-6, .1)) + 
  geom_abline(linetype = "dashed", size = .1) + 
  scale_color_manual("",values = c("CD4"="#009966", "CD8"="#333399","CD4CD8"="black", "LOW_CD" = "purple"), na.value = "darkgray") + 
  facet_wrap(~forcats::fct_infreq(factor(ptid2)), scale = "free", ncol = 2) + 
  theme(legend.position = "none") + 
  geom_line(data = boundary, aes(x = x, y = ypos), col = "red", size = .1, linetype = 'dashed')+
  geom_line(data = boundary, aes(x = x, y = yneg), col = "red", size = .1, linetype = 'dashed')+
  ylab("") + xlab("") + 
  #geom_text(data = df_exp%>% filter(ptid %in% my_ptids), aes(x= 1E-6, y = 2E-1, label = CD8 ), col = "blue",size =2) +
  #geom_text(data = df_exp%>% filter(ptid %in% my_ptids), aes(x= 1E-6, y = 6E-2, label = CD4 ), col = "darkgreen",size =2)+
  #geom_text(data = df_exp%>% filter(ptid %in% my_ptids), aes(x= 1E-6, y = 2E-2, label = UNK ), col = "gray",size =2)+
  #geom_text(data = df_con%>% filter(ptid %in% my_ptids), aes(x= 1E-2, y = 2E-5, label = CD8 ), col = "blue",size =2) +
  #geom_text(data = df_con%>% filter(ptid %in% my_ptids), aes(x= 1E-2, y = 6E-6, label = CD4 ), col = "darkgreen",size =2)+
  #geom_text(data = df_con%>% filter(ptid %in% my_ptids), aes(x= 1E-2, y = 2E-6, label = UNK ), col = "gray",size =2)+
  #geom_text(data = d_n_non_expand%>% filter(ptid %in% my_ptids), aes(x= 5E-4, y = 1E-3, label = CD8 ), col = "blue",size =2) +
  #geom_text(data = d_n_non_expand%>% filter(ptid %in% my_ptids), aes(x= 5E-4, y = 2E-4, label = CD4 ), col = "darkgreen",size =2)
  theme(axis.text.y = element_text(size = global_font_size )) + 
  theme(axis.text.x = element_text(size = global_font_size , angle = 0)) + 
  theme(axis.title.x = element_text(size = global_font_size)) + 
  theme(axis.title.y = element_text(size = global_font_size)) +
  ylab("E03 TRB (frequency)") + xlab("E01 TRB (frequency)") + 
  theme(axis.line = element_line(size = .2)) + 
  theme(strip.text = element_text(size = global_font_size)) + 
  theme(strip.background = element_blank())  + 
  theme(axis.ticks = element_blank()) + 
  geom_point(data = d_non_change%>% filter(ptid %in% my_ptids), aes(x = E01_pfreq, y = E03_pfreq, col = cell_type), size =.2, alpha = .4, pch = 1) 

# Figure 2a 
pdf(figure_2a_filename , width = 4.5, height = 2.25)
ggx_2ptid_selected
dev.off()

# Data
bind_rows(
 df_expansions %>%
    filter(ptid %in% my_ptids) %>%
    mutate(ptid = stringr::str_replace(ptid, pattern = "^15", ""))%>%
    select(pubid = ptid,
           cdr3_b_aa, 
           v_b_gene,
           j_b_gene,
           cdr3_b_nucseq,
           longitudinal_id, 
           E01_pfreq, 
           E03_pfreq, 
           cell_type) %>%
    mutate(data_type = "expansion"),
  
  df_contractions %>% 
    filter(ptid %in% my_ptids) %>%
    mutate(ptid = stringr::str_replace(ptid, pattern = "^15", ""))%>%
    
    select(pubid = ptid,
           uid= longitudinal_id, 
           E01_pfreq, 
           E03_pfreq, 
           cell_type) %>%
    mutate(data_type = "contraction"),
  
  d_non_change %>% 
    filter(ptid %in% my_ptids) %>% 
    mutate(ptid = stringr::str_replace(ptid, pattern = "^15", ""))%>%

    select(pubid = ptid,
           cdr3_b_aa, 
           v_b_gene,
           j_b_gene,
           cdr3_b_nucseq,
           uid= plot_clone_id,
           E01_pfreq, 
           E03_pfreq, 
           cell_type) %>%
    mutate(data_type = "no_sig_change")
) %>% 
  write.csv(f2a_fdata, row.names = F)





df_expansions %>% group_by(ptid) %>% 
  tally() %>% 
  arrange(desc(n)) %>% pull(ptid) %>%
  stringr::str_replace(., "15", "P")->ptid_order

gg_up=df %>% group_by(ptid, longitudinal_id) %>% 
  slice(1) %>% 
  filter(cell_type %in% c("CD4","CD8") | is.na(cell_type)) %>%
  filter(E03_vac_fc > fc) %>%
  group_by(ptid,cell_type) %>% 
  tally()%>%
  mutate(ptid2 = stringr::str_replace(ptid, "15", "P")) %>% 
  filter(ptid2 %in% ptid_order) %>%
  ggplot(., aes(x = factor(ptid2, levels = ptid_order), y = n, fill= cell_type)) + geom_bar(stat = "identity", width = .5) + 
  scale_fill_manual("",values = c("CD4"="#009966", "CD8"="#333399","CD4CD8"="black", "LOW_CD" = "purple"), na.value = "gray") + 
  theme_classic() + 
  xlab("\n")+
  ylab(paste0("Exexpanded clones [log2(FC) < 2]") ) + 
  theme(legend.position = "None") + 
  theme(axis.text.x = element_text(angle = 90, size = 6)) + 
  theme(axis.title = element_text(size = global_font_size))+
  scale_y_continuous(expand = c(0, 0)) 

gg_non_expand = d_non_change %>%
  group_by(ptid, cell_type)%>% 
  tally() %>% 
  filter(cell_type %in% c("CD4","CD8") | is.na(cell_type)) %>%
  mutate(ptid2 = stringr::str_replace(ptid, "15", "P")) %>% 
  filter(ptid2 %in% ptid_order) %>%
  ggplot(., aes(x = factor(ptid2, levels = ptid_order), y = n, fill= cell_type)) + geom_bar(stat = "identity", width = .5) + 
  scale_fill_manual("",values = c("CD4"="#009966", "CD8"="#333399","CD4CD8"="black", "LOW_CD" = "purple"), na.value = "gray") + 
  theme_classic() + 
  xlab("\n")+
  ylab(paste0("Unexpanded clones [log2(FC) < 2]") ) + 
  theme(legend.position = "None") + 
  theme(axis.text.x = element_text(angle = 90, size = 6)) + 
  theme(axis.title = element_text(size = global_font_size)) +
  scale_y_continuous(expand = c(0, 0)) 

pdf(figure_2b_filename, width = 3, height =2.25)
gridExtra::grid.arrange(gg_up + coord_cartesian(ylim= c(0,400)), 
                        gg_non_expand + coord_cartesian(ylim= c(0,400)), ncol = 2)
dev.off()

# export data types
bind_rows(
  gg_up$data %>% 
    ungroup %>%
    select(pubid = ptid2, 
           cell_type, 
           n) %>% 
    mutate(data_type = "expansion_count"),
  gg_non_expand$data %>% 
    ungroup %>%
    select(pubid = ptid2, 
           cell_type, 
           n) %>% 
    mutate(data_type = "no_sig_change_count")) %>%
write.csv(f2b_fdata, row.names = F)
  

# Extended Figure 5
ggx4 = ggplot(df_expansions ,aes(x = E01_pfreq, y = E03_pfreq)) + 
  
  geom_point(aes(col = cell_type), pch = 2, size = .6, alpha = .5)+
  geom_point(data = df_contractions, aes(col = cell_type), size = .6, pch = 6, alpha = .5)+
  geom_point(data = d_non_change, aes(x = E01_pfreq, y = E03_pfreq, col = cell_type), size =.4, alpha = .5, pch = 1) +
  scale_x_log10(expand = c(.05,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_y_log10(expand = c(.05,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  theme_classic() + 
  annotation_logticks(sides = "lb", size = .1, 
                      short = unit(0.005, "npc"),
                      mid = unit(0.0075, "npc"),
                      long = unit(0.01, "npc"),) + 
  coord_cartesian(xlim = c(1E-6, .1), ylim = c(1E-6, .1)) + 
  geom_abline(linetype = "dashed", size = .1) + 
  scale_color_manual("",values = c("CD4"="#009966", "CD8"="#333399","CD4CD8"="black", "LOW_CD" = "purple"), na.value = "darkgray") + 
  facet_wrap(~forcats::fct_infreq(factor(ptid2)), scales = "free",ncol = 3) +
  theme(legend.position = "none") + 
  geom_line(data = boundary, aes(x = x, y = ypos), col = "red", size = .1, linetype = 'dashed')+
  geom_line(data = boundary, aes(x = x, y = yneg), col = "red", size = .1, linetype = 'dashed')+
  #geom_text(data = df_exp, aes(x= 5E-4, y = 6E-2, label = CD8 ), col = "blue",size =2) +
  #geom_text(data = df_exp, aes(x= 5E-4, y = 2E-2, label = CD4 ), col = "darkgreen",size =2)+
  #geom_text(data = df_exp, aes(x= 5E-4, y = 1E-2, label = UNK ), col = "gray",size =2)+
  #geom_text(data = df_con, aes(x= 5E-4, y = 2E-5, label = CD8 ), col = "blue",size =2) +
  #geom_text(data = df_con, aes(x= 5E-4, y = 6E-6, label = CD4 ), col = "darkgreen",size =2)+
  #geom_text(data = df_con, aes(x= 5E-4, y = 2E-6, label = UNK ), col = "gray",size =2)+
  #geom_text(data = d_n_non_expand, aes(x= 5E-4, y = 1E-3, label = CD8 ), col = "blue",size =2) +
  #geom_text(data = d_n_non_expand, aes(x= 5E-4, y = 2E-4, label = CD4 ), col = "darkgreen",size =2)->ggx4
  theme(axis.text.y = element_text(size = global_font_size )) + 
  theme(axis.text.x = element_text(size = global_font_size , angle = 0)) + 
  theme(axis.title.x = element_text(size = global_font_size)) + 
  theme(axis.title.y = element_text(size = global_font_size)) +
  ylab("E03 TRB (frequency)") + xlab("E01 TRB (frequency)") + 
  theme(strip.text = element_text(size =  global_font_size)) + 
  theme(strip.background = element_blank())  + 
  theme(axis.ticks = element_blank()) 


pdf(ex_fig_4_filename , width = 6, height = 8)
ggx4
dev.off()


# Data
bind_rows(
  df_expansions %>%
    mutate(ptid = stringr::str_replace(ptid, pattern = "^15", ""))%>%
    select(pubid = ptid,
           cdr3_b_aa, 
           v_b_gene,
           j_b_gene,
           cdr3_b_nucseq,
           longitudinal_id, 
           E01_pfreq, 
           E03_pfreq, 
           cell_type) %>%
    mutate(data_type = "expansion"),
  
  df_contractions %>% 
    mutate(ptid = stringr::str_replace(ptid, pattern = "^15", ""))%>%
    
    select(pubid = ptid,
           uid= longitudinal_id, 
           E01_pfreq, 
           E03_pfreq, 
           cell_type) %>%
    mutate(data_type = "contraction"),
  
  d_non_change %>% 
    mutate(ptid = stringr::str_replace(ptid, pattern = "^15", ""))%>%
    select(pubid = ptid,
           cdr3_b_aa, 
           v_b_gene,
           j_b_gene,
           cdr3_b_nucseq,
           uid= plot_clone_id,
           E01_pfreq, 
           E03_pfreq, 
           cell_type) %>%
    mutate(data_type = "no_sig_change")
) %>% 
  write.csv(ef4a_fdata, row.names = F)


