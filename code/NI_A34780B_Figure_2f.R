# NI_A34780B_Figure_2f.R

# DEPENDENCIES
require(dplyr)
require(ggplot2)
require(scales)
require(tidyr)
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

filename_ics_summary = file.path(
  repo_loc, repo,
  'data/E01_E02_E03_ics_summary.tsv')

# OUTPUTS
global_font_size = 8
figure_2f_filename = file.path(repo_loc, repo, 'figures/Fig_2f_NI_A34780C.pdf')
fig_2f_figdata = file.path(repo_loc, repo, 'figure_data_files/NI_A34780C_fig_2f_data.csv')


# SCRIPT
data = readr::read_tsv(filename_ics_summary)
data = data %>% 
  select(subject = SubjShort, 
         visit = COLLECTION, 
         antigen   = ANTIGEN,
         `CD4 [2+]` = `CD4[2+]`,
         `CD8 [IFNg+]` = `CD8[Any_IFNg]`,
         `CD4 [IL-2+]` = `CD4[Any_IL2]`) %>%
  tidyr::gather(key, value, -subject, -visit, -antigen) %>% 
  filter(antigen == "peptide") 

# Figure 2f (CD4+ T cells)
gg4 = data %>% 
  filter(key == "CD4 [IL-2+]") %>%
  ggplot(aes(x = visit, y =value)) + 
  geom_point(col = "#009966", size = .35, alpha = .75) + 
  geom_line(aes(group = subject), col = "#009966", size = .25, alpha = .5) + 
  theme_classic() +
  scale_x_discrete(expand = c(0.05,0.05))+
  scale_y_log10(limits = c(0.001, 100), 
                breaks = c(0.001, 0.01, 0.1, 1, 10,100),
                labels = c(0.001, 0.01, 0.1, 1, 10,100)) + 
  ylab("% CD4+ T cells\nIL-2+")

# Figure 2f (CD8+ T cells)
gg8 = data %>% 
  filter(key == "CD8 [IFNg+]") %>%
  ggplot(aes(x = visit, y =value)) + 
  geom_point(col = "#333399", size = .35, alpha = .75) + 
  geom_line(aes(group = subject), col = "#333399", size = .25, alpha = .5) + 
  theme_classic() +
  scale_x_discrete(expand = c(0.05,0.05))+
  scale_y_log10(limits = c(0.001, 100), 
                breaks = c(0.001, 0.01, 0.1, 1, 10,100),
                labels = c(0.001, 0.01, 0.1, 1, 10,100)) + 
  #annotation_logticks(size= .3, side= "l") + 
  ylab("% CD8+ T cells\nINFg+")

gg8
extra_theme = theme(axis.text.y = element_text(size = global_font_size, angle = 0), 
                    axis.text.x = element_text(size = global_font_size, angle = 90), 
                    axis.title = element_text(size = global_font_size))

# Median values
gg4$data %>% group_by(visit) %>% 
  summarise(n = n(), median = median(value)) %>% 
  write.csv()
# "","visit","n","median"
# "1","E01",7,0.25157
# "2","E02",14,0.336555
# "3","E03",14,0.240735
gg8$data %>% group_by(visit) %>% 
  summarise(n = n(), median = median(value)) %>% 
  write.csv()
# "","visit","n","median"
# "1","E01",7,0.01227
# "2","E02",14,0.175
# "3","E03",14,0.45736
# Figure 2f 
pdf(figure_2f_filename, width = 2.25, height =1.5)
gridExtra::grid.arrange(
  gg4 + 
    ggpubr::stat_compare_means(comparisons = list(c("E02", "E03")), 
                               method = "wilcox",  label ="p.signif", 
                               paired = T, size = 4) + 
    extra_theme,
  gg8 + 
    ggpubr::stat_compare_means(comparisons = list(c("E02", "E03")), 
                               label ="p.signif", method = "wilcox", 
                               paired = T, size = 4) + 
    extra_theme,
  ncol = 2)
dev.off()

data %>% 
  filter(ANTIGEN == "peptide")%>%
  select(pubid = SUBJECT, 
         visit = COLLECTION,
         cd8_any_ifng = `CD8[Any_IFNg]`, 
         cd4_any_il2 = `CD4[Any_IL2]`) %>% 
  mutate(variable = "percent_parent_pos") %>% 
  write.csv(fig_2f_figdata, row.names = F)


# confirm paired wilcox tests CD8+/CD4+ E02 to E03
wide8 = gg8$data %>% select(subject, visit, value) %>% 
  tidyr::spread(key =visit, value = value) %>% 
  select(E01, E02, E03)
wide4 = gg4$data %>% select(subject, visit, value) %>% 
  tidyr::spread(key =visit, value = value)


gg8$data %>% group_by(visit) %>% summarise(median(value))
gg4$data %>% group_by(visit) %>% summarise(median(value))


print("tests")
wilcox.test(wide4$E01, wide4$E02, paired = T)
wilcox.test(wide4$E03, wide4$E02, paired = T)
wilcox.test(wide8$E01, wide8$E02, paired = T)
wilcox.test(wide8$E03, wide8$E02, paired = T)


