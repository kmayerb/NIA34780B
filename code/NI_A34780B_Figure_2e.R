# fig_2e_filename    = file.path(repo_loc,'NI_A34780B/figures/Fig_2e_NI_A34780B.pdf') 
# March 15, 2023
# /Volumes/fg_data/koelle_covid/2023_project_code/hybrid_t_cell_response/project_R_scripts
# see: 2023-03-08-median-expansion.R on server
# Reviewer 3 : The observation that mainly CD8+ T cells, but not CD4+ T cells expand upon vaccination after prior infection is interesting. Can the authors please plot the mean or median fold expansion of SARS-CoV-2 reactive CD4+ and CD8+ TCRs in Figure 2 as well? Since CD8+ T cells expand so strongly upon vaccination could this outcompete CD4+ responses and maybe even negatively affect protective immune responses following too many booster vaccinations?

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

filename_f1 = file.path(
  repo_loc, repo,
  '/data/SupTable7.longitudinal_matrix.csv')


# OUTPUTS
global_font_size = 8
figure_2e_part1_filename = file.path(repo_loc, repo, 'figures/Fig_2e_part1_NI_A34780C.pdf')
figure_2e_part2_filename = file.path(repo_loc, repo, 'figures/Fig_2e_part2_NI_A34780C.pdf')
f2e_figdata = file.path(repo_loc, repo, 'figure_data_files/NI_A34780C_fig_2e_data.csv')


# Xnote: ref_dir =  '/fh/fast/gilbert_p/fg_data/koelle_covid/2023_project_code/hybrid_t_cell_response/reference/'
# Xnote: f = 'SupTable7.longitudinal_matrix.csv'

df = readr::read_csv(filename_f1) %>% 
  filter(cell_type_2_y%in% c("CD4","CD8")) %>%
  filter(`10x_cdf_p_fdr_10x` < 0.05) %>% 
  mutate(person_x = paste0("P", person_x))

gg3 = ggplot(df , aes(x = factor(person_x), 
                      y = E03_vac_fc,
                      fill = cell_type_2_y,
                      alpha = .5)) + 
  geom_boxplot(outlier.size = .1, width = .5, size = .1) +
  theme_classic() +
  scale_y_continuous(limits = c(0.1, 2048),
                     trans='log2',
                     breaks=c(2,8,32,128,512,2048))+ #trans_breaks('log2', function(x) 2^x))+
  #labels=trans_format('log2', math_format(2^.x)))+
  scale_fill_manual(values = c("CD4"="#009966", "CD8"="#333399")) + 
  ylab("Fold change\n(E01 vs. E03)") + 
  xlab("") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size =global_font_size, angle = 90))+
  theme(axis.text.y = element_text(size =global_font_size))+
  theme(axis.title = element_text(size =global_font_size)) 
gg3

gg2 = ggplot(df , aes(x = factor(person_x), 
                      y = E02_vac_fc,
                      fill = cell_type_2_y,
                      alpha = .5)) + 
  geom_boxplot(outlier.size = .1, width = .5, size = .1) +
  theme_classic() +
  scale_y_continuous(
    limits = c(0.1, 2048),
    trans='log2',
    breaks=c(2,8,32,128,512,2048))+ #trans_breaks('log2', function(x) 2^x))+
  #labels=trans_format('log2', math_format(2^.x)))+
  scale_fill_manual(values = c("CD4"="#009966", "CD8"="#333399")) + 
  ylab("Fold change\n(E01 vs. E02)") +
  xlab("") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  theme(legend.position = "none") + 
  theme(axis.text.x = element_text(size = global_font_size, angle = 90))+
  theme(axis.text.y = element_text(size = global_font_size))+
  theme(axis.title = element_text(size = global_font_size)) 


pdf(figure_2e_part1_filename, width = 3, height = 1.5)
gg2 
dev.off()

pdf(figure_2e_part2_filename, width = 3, height = 1.5)
gg3
dev.off()

dput(names(df))
df %>% 
  mutate(person_x = stringr::str_replace(person_x, pattern = "P",replacement = ''))%>%
  select(pubid = person_x, 
         E02_vac_fc,
         E03_vac_fc,
         clone_id, 
         v_a_gene, j_a_gene, cdr3_a_aa, v_b_gene, j_b_gene, cdr3_b_aa, cdr3_a_nt, cdr3_b_nt) %>% 
  write.csv(f2e_figdata, row.names = F)
  

