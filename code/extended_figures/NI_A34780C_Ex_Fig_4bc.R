# NI_A34780C_Ex_Fig_4bc.R
# Uses inputs from NI_A34780C_Ex_Fig_4bc.py

# XNote: Was Ext 5 now Ext 4

## clear env
rm(list=ls())
## PACKAGES
library(readr)
library(ggplot2)

# INPUTS 
# set repo location as appropriate.
user = 'kmb'
if(user == 'esf') {
  repo_loc = '/Volumes/fh/fast/corey_l/esford3_kmayerbl_collab/software' ## if  you're emily
  repo = 'NIA34780B'
} else if(user == 'kmb') {
  repo_loc = '/Users/kmayerbl/active/david_koelle' ## if you're kosh
  repo = 'NI_A34780B'
} else {
  stop("set repo loc and repo manually") ## delete these and set repo loc, repo for your data
}

# inputs
E01_E03_525_581_filename = file.path(repo_loc, repo, 'data/E01_E03_10x_match_525_581_plot_unique.csv')

## OUTPUTS
ex_fig_5bc_filename = file.path(repo_loc, repo, 'figures/Ex_Fig_4bc_NI_A34780C.pdf')
ex_fig_5bc_filename_legend = file.path(repo_loc, repo, 'figures/Ex_Fig_4bc_NI_A34780C_legend.pdf')
# Figure Data
ef5b_fdata = file.path(repo_loc, repo, 'figure_data_files',' NI_A34780C_ex_fig_4b.data.csv')
ef5c_fdata = file.path(repo_loc, repo, 'figure_data_files',' NI_A34780C_ex_fig_4c.data.csv')


# all E01,E03 compbinations
d=readr::read_csv(E01_E03_525_581_filename)%>%
  mutate(E03= ifelse(E03 < 1E-6, 1E-6, E03)) %>% 
  mutate(E01=ifelse(E01 < 1E-6, 1E-6, E01))

daim = d %>% filter(!is.na(`10x_pfreq`)) 

b <- ggplot(d %>% filter(ptid == 525), 
            aes(x = E01, y = E03)) +
  geom_point(pch = 19, size = .1, col = "gray") + 
  geom_point(pch = 15, data = daim %>% filter(ptid == 525), 
             aes(col = log10(`10x_pfreq`)), size = .75)+
  scale_y_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) + 
  scale_x_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) +
  theme_classic()+
  annotation_logticks(side = 'lb', size = .2)+
  theme(axis.text =element_text(size = 8))+
  theme(axis.title =element_text(size = 8)) +
  theme(strip.background = element_blank()) +
  scale_color_viridis("", option = "plasma") +
  facet_wrap(~ptid) +
  theme(legend.position = "top")

c <- ggplot(d %>% filter(ptid == 581), 
            aes(x = E01, y = E03)) +
  geom_point(pch = 19,size = .1, col = "gray") + 
  geom_point(pch = 15, data = daim %>% filter(ptid == 581), 
             aes(col = log10(`10x_pfreq`)), size = .75)+
  scale_y_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) + 
  scale_x_log10(breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) +
  theme_classic()+
  annotation_logticks(side = 'lb', size = .2)+
  theme(axis.text =element_text(size = 8))+
  theme(axis.title =element_text(size = 8)) +
  theme(strip.background = element_blank()) +
  scale_color_viridis("", option = "plasma") +
  facet_wrap(~ptid)+
  theme(legend.position = "top")


# create legend
pdf(ex_fig_5bc_filename_legend, width = 4, height = 2.5)
gridExtra::grid.arrange(b,c, ncol = 2)
dev.off()

pdf(ex_fig_5bc_filename, width = 3.5, height = 1.75)
gridExtra::grid.arrange(b+theme(legend.position = "none"),
                        c+theme(legend.position = "none") , ncol = 2)
dev.off()


b$data %>% 
  select(pubid = ptid, 
         e01 = E01, 
         e03 = E03, 
         `10x_pfreq`,
         count) %>%
         write.csv(ef5b_fdata, row.names = FALSE)

c$data %>% 
  select(pubid = ptid, 
         e01 = E01, 
         e03 = E03, 
         `10x_pfreq`,
         count) %>%
  write.csv(ef5c_fdata, row.names = FALSE)


