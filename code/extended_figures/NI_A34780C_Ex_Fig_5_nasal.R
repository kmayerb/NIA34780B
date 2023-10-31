# Extended Data Figure Showing CD69+CD137+, TRB matches in E05 
# Nasal repertoires.

# NI_A34780B_Ext_Fig_5_nasal_TRB.R

# OUTPUTS 
# Ext Data Fig 5A
# Ext Data Fig 5B

# ex_fig5a_data
# ex_fig5b_data

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

filename_nasal1= file.path(repo_loc, repo, 
                           'data/2022_06_23_all_AIM_phenotypes_in_nasal_pubilc.tsv')
filename_nasal2= file.path(repo_loc, repo, 
                           'data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv')
ex_fig5a_data =  file.path(repo_loc, repo, 
                          'figure_data_files/ex_fig5a_nasal_data.csv')
ex_fig5b_data =  file.path(repo_loc, repo, 
                           'figure_data_files/ex_fig5b_nasal_data.csv')
#d = readr::read_tsv('/Users/kmayerbl/active/david_koelle/hybrid_t_cell_response/plot_data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv') %>% 
my_ptids = c('15754', '15548', '15782', '15577', '15525', '15515', '15630',
             '15684', '15669', '15673','15761', '15527', '15581', '15655', 
             '15744', '15836', '15531')

d = readr::read_tsv(filename_nasal2) %>% 
  mutate(ptid = as.character(ptid))%>%
  filter(ptid %in% my_ptids)

d = d %>% mutate(ptid = paste0('P', ptid))
length(d$ptid %>% unique())

table(d$cell_type, d$ptid)
d$short_code = stringr::str_extract(d$filename, pattern = '15[0-9]{3}')%>% 
  stringr::str_replace(., pattern = "^15", replacement = "P")
dct = d %>% filter(!is.na(cell_type)) %>% filter(`10x_cdf_p_fdr_10x` < 0.05)
dtct = dct %>% filter(cell_type %in% c('CD4', 'CD8'))

length(15673 %in% unique(d$ptid_x))

gall = ggplot(d, aes(x = rank, y = frequency)) + 
  geom_point(size = .1, color = 'gray', pch = 19) + 
  facet_wrap(~short_code,  scales = "free", ncol = 3) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)),
                limits = c(1E-5, 1E-1))+
  theme_classic() + 
  annotation_logticks(side = "l",size = .2)+
  geom_point(data = dtct %>% filter(cell_type == "CD4"), 
             size = 1.5,  aes(color = cell_type, y= 1E-5), pch = "|") + 
  geom_point(data = dtct %>% filter(cell_type == "CD8"), 
             size = 1.5,  aes(color = cell_type, y= 3E-5), pch = "|") + 
  scale_color_manual(values = c("CD4"="#009966", "CD8"="#333399")) + 
  ylab("E05 Nasal TRB\n(frequency)") + xlab("E05 TRB (rank abundance)") + 
  theme_classic()+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+ 
  theme(strip.text = element_text(size = 7))+
  theme(axis.title = element_text(size = 7))+
  theme(axis.text = element_text(size = 7))+
  theme(axis.line = element_line(size = .1)) +
  theme(panel.border = element_blank())+
  theme(panel.spacing = grid::unit(1, "lines"))

pdf('figures/NIA34780C_Ex_Fig_5_nasal_rank3_all__2.pdf', width =6, height = 9)
gall
dev.off()

bind_rows(
  gall$data %>% select(pubid = short_code, rank, frequency) %>%
    mutate(type = "nasal"), 
  dtct %>% select(pubid = short_code, rank ,frequency, type = cell_type)
) %>% 
  write.csv(ex_fig5a_data, row.names = F)





my_ptids = c('15673')
#d = readr::read_tsv('/Users/kmayerbl/active/david_koelle/hybrid_t_cell_response/plot_data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv') %>% 
d = readr::read_tsv(filename_nasal2) %>% 
  filter(ptid%in% my_ptids)
d$short_code = stringr::str_extract(d$filename, pattern = '15[0-9]{3}')%>% stringr::str_remove(., pattern = "^15")

draw = readr::read_tsv('/Users/kmayerbl/active/david_koelle/hybrid_t_cell_response/plot_data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv') 

draw %>% filter(!is.na(cell_type)) %>% 
  pull(ptid) %>% unique() %>% length()


# Validated TCRs
cloned_tcrs = tibble(
  'label'    = c('TCR1','TCR2','TCR3','TCR4','TCR5','TCR8.1','TCR8.2','TCR14'),
  'v_b_gene' = c('TRBV14*01','TRBV19*01','TRBV9*01','TRBV20-1*01','TRBV5-1*01','TRBV20-1*01','TRBV3-1*01','TRBV20-1*01'),
  'j_b_gene' = c('TRBJ1-1*01','TRBJ2-1*01','TRBJ1-5*01','TRBJ2-3*01','TRBJ1-1*01','TRBJ2-5*01','TRBJ2-5*01','TRBJ2-3*01'),
  'cdr3_b_aa'= c('CASRRFGDTEAFF','CASSIKASSYNEQFF','CASSAWGGNQPQHF','CSARDLGGDTQYF','CASKDSLNTEAFF','CSARSWGSETQYF','CASRPLGEETQYF','CSARVIGSASLQYF'),
  'v_a_gene' = c('TRAV8-2','TRAV13-1','TRAV26-1','TRAV9-2','TRAV12-2','TRAV8-6','TRAV8-6','TRAV13-1'),
  'j_a_gene' = c('TRAJ34','TRAJ52','TRAJ39','TRAJ10','TRAJ34','TRAJ20','TRAJ20','TRAJ48'),
  'cdr3_a_aa'= c('CVVSEKNTDKLIF','CAARGVDAGGTSYGKLTF','CIVTDNNAGNMLTF','CALSDKKLTGGGNKLTF','CAVKENTDKLIF','CAVSPRSNDYKLSF','CAVSARSNDYKLSF','CAASSNFGNEKLTF'),
)

d_confirmed = d %>% 
  left_join(cloned_tcrs, by = c("v_b_gene","cdr3_b_aa", 'j_b_gene','cdr3_a_aa', 'v_a_gene','j_a_gene')) %>%
  filter(!is.na(label)) %>% 
  filter(ptid == '15673') %>% 
  arrange(desc(templates)) %>% 
  group_by(label) %>% 
  slice(1)

d$`10x_cdf_p_fdr_10x` 
dct = d %>% filter(!is.na(cell_type)) %>% filter(`10x_cdf_p_fdr_10x` < 0.05)
dim(dct)
dtct = dct %>% filter(cell_type %in% c('CD4', 'CD8'))
dct100 = dct %>% filter(rank <100)

require(scales)
g673 = ggplot(d, aes(x = rank, y = frequency)) + 
  geom_point(size = .1, color = 'black') + 
  facet_wrap(~short_code, scale= "free_x") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = trans_format("log10", math_format(10^.x)))+
  theme_classic() + 
  annotation_logticks(side = "l",size = .2)+
  geom_point(data = dct, size = 1.5, aes(color = cell_type), pch = 1) + 
  scale_color_manual(values = c("CD8"="blue" , "CD4"= "darkgreen")) + 
  ylab("") + xlab("") + 
  geom_text(data= d_confirmed, aes(label = label), col = "blue", size = 2, nudge_y = -.1, nudge_x = -.1)+
  theme(strip.background = element_blank()) +
  theme(legend.position = "none") + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))

pdf('figures/nasal_rank_673.pdf', width = 4, height = 4)
g673
dev.off()


bind_rows(
  g673$data %>% 
    select(pubid = short_code, rank, frequency) %>%
    mutate(type = "nasal"), 
  d %>% 
    select(pubid = short_code, rank ,frequency, type = cell_type)
) %>% 
  write.csv(ex_fig5b_data, row.names = F)

  
  
























###########################
# Additional Code Note used
# ex_fig5_data 
# d %>% filter(ptid_x == 15525)
# 
# dtct %>% filter(cell_type == "CD4") %>% filter(!is.na(`10x_cdf_p_fdr_10x`)) %>% dim()
# 
# #/Users/kmayerbl/active/david_koelle/hybrid_t_cell_response/plot_data/2022_06_23_all_AIM_phenotypes_in_nasal_pubilc.tsv
# #/Users/kmayerbl/active/david_koelle/hybrid_t_cell_response/plot_data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv'
# 
# # filename --  generate_nose_rank_match.py
# d = readr::read_tsv(filename_nasal1) %>% 
#   filter(ptid_x%in% c('15754', '15548', '15782', '15577', '15525', '15515', '15630',
#                       '15684', '15669', '15761', '15527', '15581', '15655', '15673',
#                       '15744', '15836', '15531'))
# 
# dct = d %>% filter(!is.na(cell_type)) %>% filter(`10x_cdf_p_fdr_10x` < 0.05)
# dim(dct)
# dtct = dct %>% filter(cell_type %in% c('CD4', 'CD8'))
# dct100 = dct %>% filter(rank <100)
# ggplot(d, aes(x = rank, y = frequency)) + 
#   geom_point(size = .1, color = 'gray') + facet_wrap(~filename, scale= "free_x") + 
#   scale_y_log10() + 
#   theme_classic() + 
#   annotation_logticks(side = "l")+
#   geom_point(data = dct, size = 1.5, aes(color = cell_type), pch = 1) + 
#   scale_color_manual(values = c("CD8"="blue" , "CD4"= "darkgreen"))  +
#   ggrepel::geom_text_repel(data = dct100, aes(label = cdr3_b_aa, color = cell_type),size = 1.5)
# 
# 
# d %>% group_by(cell_type, filename) %>%
#   summarise(sum_freq = sum(frequency)) %>%
#   tidyr::spread(key = cell_type, value = sum_freq) %>% 
#   arrange(desc(CD8))
# 
# 
# 
# 
# 
# # 
# # 
# # #d = readr::read_tsv('/Users/kmayerbl/active/david_koelle/hybrid_t_cell_response/plot_data/2022_06_23_all_AIM_phenotypes_in_nasal.tsv') %>% 
# # my_ptids = c('15754', '15548', '15782', '15577', '15525', '15515', '15630',
# #              '15684', '15669', '15761', '15527', '15581', '15655', 
# #              '15744', '15836', '15531')
# # d = readr::read_tsv(filename_nasal1) %>% 
# #   mutate(ptid = ptid_x)%>%
# #   filter(ptid_x%in% my_ptids)
# # table(d$cell_type, d$ptid)
# # d$short_code = stringr::str_extract(d$filename, pattern = '15[0-9]{3}')%>% stringr::str_remove(., pattern = "^15")
# # dct = d %>% filter(!is.na(cell_type)) %>% filter(`10x_cdf_p_fdr_10x` < 0.05)
# # dtct = dct %>% filter(cell_type %in% c('CD4', 'CD8'))
# # 
# # d = d %>% mutate(ptid = paste0('P', ptid))
# # 
# # gall = ggplot(d, aes(x = rank, y = frequency)) + 
# #   geom_point(size = .1, color = 'black', pch = 19) + 
# #   facet_wrap(~short_code, ncol = 4, scales = "free") + 
# #   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
# #                 labels = trans_format("log10", math_format(10^.x)),
# #                 limits = c(1E-5, 1E-1))+
# #   theme_classic() + 
# #   annotation_logticks(side = "l",size = .2)+
# #   geom_point(data = dtct, size = 1.5, aes(color = cell_type), pch = 1) + 
# #   scale_color_manual(values = c("CD8"="blue" , "CD4"= "darkgreen")) + 
# #   ylab("") + xlab("") + 
# #   theme(strip.background = element_blank()) +
# #   theme(legend.position = "none") + 
# #   theme(panel.border = element_rect(colour = "black", fill=NA, size=1))
# # 
# # pdf('figures/nasal_rank_all.pdf', width = 6, height = 6)
# # gall
# # dev.off()
# # 
# # 
# # 
# # 
# # rankme <- function(x){
# #   x = arrange(x, desc(frequency))
# #   x = mutate(x, rank2 = seq_along(x$frequency))
# #   return(x)
# # }
# # 
# # d2 = d %>% group_by(filename) %>% group_split() %>% 
# #   purrr::map(.,  ~rankme(.x)) %>% 
# #   do.call(rbind, .)
# # 
# # ggplot(d2, aes(x = rank2, y = frequency)) + geom_point(size = .1) + facet_wrap(~filename) + 
# #   scale_y_log10()
# 
