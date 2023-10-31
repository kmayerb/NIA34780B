### R code for Figure 6. 
# files: 

## packages
library(readr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)
library(data.table)
library(grDevices)
library(stringr)
library(cowplot)
library(corrplot)

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

filename_clin_corr = file.path(repo_loc, repo, 'data/all_clinical_correlates_101123.csv')
df <- read.csv(filename_clin_corr)

##OUTPUTS
fig_6a_filename    = file.path(repo_loc, repo, 'figures/Fig_6a_NI_A34780C.pdf')
fig_6b_filename    = file.path(repo_loc, repo, 'figures/Fig_6b_NI_A34780C.pdf')
fig_6c_filename    = file.path(repo_loc, repo, 'figures/Fig_6c_NI_A34780C.pdf')
fig_6d_filename    = file.path(repo_loc, repo, 'figures/Fig_6d_NI_A34780C.pdf')
fig_6e_filename    = file.path(repo_loc, repo, 'figures/Fig_6e_NI_A34780C.pdf')
fig_6f_filename    = file.path(repo_loc, repo, 'figures/Fig_6f_NI_A34780C.pdf')
fig_6g_filename    = file.path(repo_loc, repo, 'figures/Fig_6g_NI_A34780C.pdf')


f6a_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6a_data.csv')
f6b_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6b_data.csv')
f6c_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6c_data.csv')
f6d_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6d_data.csv')
f6e_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6e_data.csv')
f6f_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6f_data.csv')
f6g_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6g_data.csv')
f6h_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_6h_data.csv')


## Fig 4a: Diagnostic CD4+ spike breadth vs hospitalization status
# get minimum (non-zero) value for log transformation
k = df$breadth_spike_classii 
min(k[k>0 & !is.na(k)])
# 6.84e-06
my_comparisons <- list(c("E01", "E02"), c("E01", "E03"),  c("E01", "E05"),
                       c("E00", "E01"), c("E00", "E02"), c("E00", "E03"), c("E00", "E05"))
# brii breadth class ii S
brii <- df %>% 
  select(subject_id, extension, breadth_spike_classii, Hospitalized) %>%
  filter(!extension %in% c("E00.5","E04")) %>%
  mutate(breadth_spike_classii = ifelse(breadth_spike_classii< 6.84e-06, 6.84e-06, breadth_spike_classii)) %>%
  ggplot(aes(fill = factor(Hospitalized), x = extension, y = breadth_spike_classii)) + 
  geom_line(aes(group = subject_id ), size = .2, alpha = .5, col = "gray") +
  geom_boxplot(alpha = .2, size = .25, width = .7) +
  geom_point(aes(fill= factor(Hospitalized)), size = .5,  pch = 21, position = position_jitterdodge(jitter.height = 0, jitter.width = 0)) +
  scale_fill_manual(values = c("#0C7BDC","#FFC20A")) + 
  theme_classic() +
  scale_y_log10(expand = c(.1,.1))+
  annotation_logticks(side = "l", size = .2) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,  method = "wilcox", size = 2, 
                             label = "p.signif", tip.length = 0, step.increase = .06, hjust = -1.5)+
  ggpubr::stat_compare_means(method = "wilcox", size = 3, label = "p.signif", label.y.npc = 0, col = "black")+
  theme(axis.title.x = element_blank()) +  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(1E-5, 5E-3)) + 
  ylab(expression(~CD4^+ ~" "~T^S~"cells diagnostic breadth"))
brii

# FIGURE 4A
pdf(fig_6a_filename, width = 3, height =2.5)
print(brii)
dev.off()

brii$data %>% 
  select(pubid=subject_id, 
         visit = extension, 
         breadth_spike_classii,
         severe = Hospitalized) %>%
  write.csv(f6a_fdata,row.names = F)


## Fig 4b: Expansion breadth vs hospitalization status
# get minimum (non-zero) value for log transformation
k = df$breadth_exp
min(k[k>0 & !is.na(k)])
# 1.367273e-05
my_comparisons <- list(c("E02","E03"), c("E03","E05"))

br_exp = df %>% select(subject_id, extension, breadth_exp, Hospitalized) %>%
  filter(!extension %in% c("E00", "E01", "E00.5","E04")) %>%
  filter(!is.na(breadth_exp)) %>%
  mutate(breadth_exp = ifelse(breadth_exp <1.367273e-05, 1.367273e-05, breadth_exp)) %>%
  ggplot(aes(fill = factor(Hospitalized), x = extension, y = breadth_exp)) + 
  #geom_line(aes(group = subject_id ), size = .1, alpha = .5, col = "gray") +
  geom_boxplot(alpha = .2, size = .25, width = .7, outlier.shape = NA) +
  geom_point(aes(fill= factor(Hospitalized)),size = .5,  pch = 21, 
             position = position_jitterdodge(jitter.height = 0, jitter.width = 0)) +
  scale_fill_manual(values = c("#0C7BDC","#FFC20A"), 
                    labels = c("Not hospitalized", "Hospitalized"),
                    name = "Hospitalization \nstatus") + 
  theme_classic() +
  scale_y_log10(expand = c(.1,.1))+
  annotation_logticks(side = "l", size = .2) + 
  #ggpubr::stat_compare_means(comparisons = my_comparisons, paired = FALSE, na.rm = FALSE, 
  #                           method = "wilcox", size = 2, label = "p.signif", tip.length = 0)+
  ggpubr::stat_compare_means(method = "wilcox", size = 3, label = "p.signif", label.y.npc = 0, col = "black")+
  theme(axis.title.x = element_blank()) +  #theme(legend.position = "none") + 
  coord_cartesian(ylim = c(1E-5, 1E-2)) + 
  ylab("Expanded clones breadth")
br_exp

pdf(fig_6b_filename, width = 4, height = 3)
print(br_exp)
dev.off()

br_exp$data %>% 
  select(pubid=subject_id, 
         visit = extension, 
         breadth_exp,
         severe = Hospitalized) %>%
  write.csv(f6b_fdata,row.names = F)


## Figure 4C - Spearman correlation 
library(corrplot)
plot <- filter(df, subject_id!="845") %>% 
  select(
    extension,
    ptid = subject_id,
    B  = breadth_adapt,
    D    = depth_adapt,
    B_S         = breadth_spike,
    D_S          = depth_spike,
    B_NS         = breadth_nonspike,
    D_NS        = depth_nonspike,
    B_8        = breadth_classi,
    D_8      = depth_classi,
    B_4       = breadth_classii,
    D_4        = depth_classii, 
    B_S_8      = breadth_spike_classi,
    D_S_8    = depth_spike_classi,
    B_S_4    = breadth_spike_classii,
    D_S_4    = depth_spike_classii,
    B_NS_8   = breadth_nonspike_classi,
    D_NS_8   = depth_nonspike_classi,
    B_NS_4   = breadth_nonspike_classii,
    D_NS_4    = depth_nonspike_classii,
    `CMV exposure`   = cmv_status,
    `Hospitalization`= Hospitalized,
    `NT50`           = NIAID_neuts) 
plot$Hospitalization <- as.numeric(plot$Hospitalization)
plot$`CMV exposure` <- as.numeric(plot$`CMV exposure`)
str(plot)
plot <- filter(plot, extension == "E00") %>% select(-extension, -ptid)
plot = plot %>% 
  filter(complete.cases(.))
correlation_matrix <- cor(plot, method = "spearman") ##put in all time points
correlation_p <- cor.mtest(plot, conf.level = 0.95, method = "spearman")

corrplot::corrplot(correlation_matrix, p.mat = correlation_p$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", 
                   method = "square", type = "upper", tl.col = "black", tl.cex = 0.7, cl.cex = 0.7, 
                   pch.col = 'grey20', pch.cex = 0.7, order = "original", diag = F)

##save as pdf - device size, 'E00 corrplot.pdf' (3.91*4.71)
pdf(fig_6c_filename, width = 4, height = 4)
corrplot::corrplot(correlation_matrix, p.mat = correlation_p$p, sig.level = c(0.001, 0.01, 0.05), insig = "label_sig", 
                   method = "square", type = "upper", tl.col = "black", tl.cex = 0.7, cl.cex = 0.7, 
                   pch.col = 'grey20', pch.cex = 0.7, order = "original", diag = F)
dev.off()


# Write out data
df %>% 
  mutate(pubid = subject_id, 
         visit = extension, 
         severe = Hospitalized )%>%
  dplyr::select(c("pubid", "visit", "severe", "breadth_adapt", "depth_adapt", 
                  "breadth_spike", "depth_spike", "breadth_nonspike", "depth_nonspike", 
                  "breadth_classi", "depth_classi", "breadth_classii", "depth_classii", 
                  "breadth_spike_classi", "depth_spike_classi", "breadth_spike_classii", 
                  "depth_spike_classii", "breadth_nonspike_classi", "depth_nonspike_classi", 
                  "breadth_nonspike_classii", "depth_nonspike_classii", "cmv_status", 
                  "NIAID_neuts", "breadth_exp")) %>% 
  write.csv(f6c_fdata, row.names = F)


# Fig 4d: diagnostic CD4+ non-spike breadth 
#  Find min non-zero breadth
k = df$breadth_spike_classii
c2_S_min = k[k>0 & !is.na(k)] %>% min()
k = df$breadth_nonspike_classii
c2_NS_min = k[k>0 & !is.na(k)] %>% min()

# Figure 4d
my_comparisons <- list(c("E00", "E05"), c("E00", "E03"), c("E00", "E02"), c("E00", "E01"), 
                       c("E01","E05"), c("E01", "E03"), c("E01", "E02"))
brSNS = df %>% select(subject_id, extension, breadth_nonspike_classii, breadth_spike_classii, Hospitalized)%>%
  filter(!extension %in% c("E00.5","E04"))%>%
  mutate(breadth_spike_classii= ifelse(breadth_spike_classii<c2_S_min , c2_S_min, breadth_spike_classii)) %>%
  mutate(breadth_nonspike_classii = ifelse(breadth_nonspike_classii< c2_NS_min, c2_NS_min, breadth_nonspike_classii)) %>% 
  tidyr::gather(key , value, -Hospitalized, -extension, -subject_id) %>% 
  #mutate(key2 = paste0(Hospitalized, key)) %>% 
  mutate(key3 = c("breadth_spike_classii"="CD4+ S", 
                  "breadth_nonspike_classii"="CD4+ NS")[key]) %>%
  mutate(spike = c("breadth_spike_classii"="S", 
                   "breadth_nonspike_classii"="NS")[key]) %>%
  mutate(key3 = factor(key3, levels = c("CD4+ S", "CD4+ NS"))) %>% 
  #filter(spike == "NS") %>%
  ggplot(aes(x = extension, y = value, fill = key3)) +
  geom_boxplot(width = .8, outlier.size =.2 , alpha = .5) + scale_y_log10() + 
  geom_point(size = .5,  pch = 21, position = position_jitterdodge(jitter.height = 0, jitter.width = 0)) +
  scale_fill_manual("",values = c("red", "gray"))+ 
  theme_classic() + 
  annotation_logticks(side = "l", size = .3) + 
  theme() + 
  theme(axis.title.x = element_blank())+
  ylab(expression(~CD4^+ ~" "~T^S~" & "~T^NS~ "cells diagnostic breadth")) + 
  ggpubr::stat_compare_means(comparisons = my_comparisons,  method = "wilcox", size = 2, 
                             label = "p.signif", tip.length = 0, step.increase = .06, vjust = 1)+
  coord_cartesian(ylim= c(6E-6, 1E-3))
brSNS

pdf(fig_6d_filename, width = 4, height = 3)
print(brSNS)
dev.off()


brSNS$data %>% 
  select(pubid = subject_id, 
         visit = extension, 
         severe = Hospitalized, 
         value, 
         variable = key, 
         variable2 = key3) %>% 
  write.csv(f6d_fdata, row.names = F)


# STAT COMPARISIONS, TO SUPPORT FIGURE 4d 
## (spike to spike and non-spike to non-spike - 
my_comparisons <- list(c("E00", "E05"), c("E00", "E03"), c("E00", "E02"), c("E00", "E01"), 
                       c("E01","E05"), c("E01", "E03"), c("E01", "E02"))

my_comparisons <- list(  #c("E00", "E01"),c("E00", "E02"),c("E00", "E03"),c("E00", "E05"),
  c("E01", "E02"), c("E01", "E03"),  c("E01","E05"))#c("E02","E03"),c("E02","E05"))
brS = df %>% select(subject_id, extension, breadth_nonspike_classii, breadth_spike_classii, Hospitalized)%>%
  filter(!extension %in% c("E00.5","E04"))%>%
  mutate(breadth_spike_classii= ifelse(breadth_spike_classii<c2_S_min, c2_S_min, breadth_spike_classii)) %>%
  mutate(breadth_nonspike_classii = ifelse(breadth_nonspike_classii< c2_NS_min, c2_NS_min, breadth_nonspike_classii)) %>% 
  tidyr::gather(key , value, -Hospitalized, -extension, -subject_id) %>% 
  #mutate(key2 = paste0(Hospitalized, key)) %>% 
  mutate(key3 = c("breadth_spike_classii"="CD4+ S", 
                  "breadth_nonspike_classii"="CD4+ NS")[key]) %>%
  mutate(spike = c("breadth_spike_classii"="S", 
                   "breadth_nonspike_classii"="NS")[key]) %>%
  mutate(key3 = factor(key3, levels = c("CD4+ S", "CD4+ NS"))) %>% 
  filter(spike == "S") %>%
  ggplot(aes(x = extension, y = value, fill = key3)) +
  geom_boxplot(width = .5, outlier.size =.2 , alpha = .5) + scale_y_log10() + 
  geom_point(size = .5,  pch = 21, position = position_jitterdodge(jitter.height = 0, jitter.width = 0)) +
  scale_fill_manual("",values = c("red","gray"))+ 
  theme_classic() + 
  annotation_logticks(side = "l") + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank())+
  ylab("CD4+ breadth")+
  ggpubr::stat_compare_means(comparisons = my_comparisons,  method = "wilcox", size = 2, label = "p.signif",tip.length = 0, step.increase = .06, vjust = 1)+
  #ggpubr::stat_compare_means(method = "wilcox", size = 3, label = "p.signif", label.y.npc = 0, col = "red" )  + 
  facet_wrap(~spike)
brS

my_comparisons <- list(  c("E01", "E02"), c("E01", "E03"),  c("E01","E05"))
brNS = df %>% select(subject_id, extension, breadth_nonspike_classii, breadth_spike_classii, Hospitalized)%>%
  filter(!extension %in% c("E00.5","E04"))%>%
  mutate(breadth_spike_classii= ifelse(breadth_spike_classii<6.84e-06, 6.84e-06 ,breadth_spike_classii)) %>%
  mutate(breadth_nonspike_classii = ifelse(breadth_nonspike_classii< 9.45e-06, 9.45e-06, breadth_nonspike_classii)) %>% 
  tidyr::gather(key , value, -Hospitalized, -extension, -subject_id) %>% 
  #mutate(key2 = paste0(Hospitalized, key)) %>% 
  mutate(key3 = c("breadth_spike_classii"="CD4+ S", 
                  "breadth_nonspike_classii"="CD4+ NS")[key]) %>%
  mutate(spike = c("breadth_spike_classii"="S", 
                   "breadth_nonspike_classii"="NS")[key]) %>%
  mutate(key3 = factor(key3, levels = c("CD4+ S", "CD4+ NS"))) %>% 
  filter(spike == "NS") %>%
  ggplot(aes(x = extension, y = value, fill = key3)) +
  geom_boxplot(width = .5, outlier.size =.2 , alpha = .5) + scale_y_log10() + 
  geom_point(size = .5,  pch = 21, position = position_jitterdodge(jitter.height = 0, jitter.width = 0)) +
  scale_fill_manual("",values = c("gray"))+ 
  theme_classic() + 
  annotation_logticks(side = "l") + 
  theme(legend.position = "none") + 
  theme(axis.title.x = element_blank())+
  ylab("CD4+ breadth")+
  ggpubr::stat_compare_means(comparisons = my_comparisons,  method = "wilcox", size = 2, label = "p.signif",tip.length = 0, step.increase = .06, vjust = 1)+
  #ggpubr::stat_compare_means(method = "wilcox", size = 3, label = "p.signif", label.y.npc = 0, col = "red" )  + 
  facet_wrap(~spike)
brNS






## Fig. 4e
# <md> summarizes how many people had each level neutralization at each point 
# <per_ext> gives you total #s of each
per_ext = df %>% filter(!extension %in% c("E00.5","E04")) %>% 
  group_by(extension, Hospitalized) %>% tally()

md  = df %>% filter(!extension %in% c("E00.5","E04")) %>% 
  group_by(extension, Hospitalized, NIAID_neuts) %>% 
  summarise(n1=n()) %>% 
  left_join(per_ext) %>% 
  mutate(p = 100*n1/n)

df$subject_id <- as.numeric(df$subject_id)
ann_text <- data.frame(extension = "E05", NIAID_neut = 220, lab = "845", Hospitalized = 1)

plot <- ggplot(md, aes(x = extension, y = NIAID_neuts, size = p, col = factor(Hospitalized))) +
  geom_point()+
  scale_y_continuous(trans='log2', breaks = sort(unique(df$NIAID_neuts))) + 
  theme_classic() + 
  geom_line(data = df %>% filter(!extension %in% c("E00.5","E04")), 
            aes(x = extension, y = NIAID_neuts, group = subject_id, ), 
            size = .4, alpha = .2) + 
  scale_y_continuous(trans = "log2", breaks = c(20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480),
                     labels = c("20", "40", "80", "160", "320", "640", "1280", "2560", "5120", "10240", "20480"), 
                     limits = c(20, 20480),
                     name = "NT50") +
  scale_x_discrete(breaks = c("E00", "E01", "E02", "E03", "E05"),
                   labels = c("E00", "E01", "E02", "E03", "E05"),
                   name = "Visit") +
  scale_size(range = c(0,3)) +
  #scale_linewidth(range = c(0.1-0.3)) + 
  theme_classic() + 
  labs(title = "Neutralization titer over time") + 
  scale_color_manual("",values = c("#0C7BDC","#FFC20A"), 
                     labels = c("Not hospitalized", "Hospitalized")) + 
  facet_wrap(~Hospitalized, ncol = 1) + 
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  guides(size = "none") +
  geom_hline(yintercept = 640, linetype = "dashed") 
pplot <- plot + geom_text(data = ann_text, aes(x=extension, y=NIAID_neut, label = lab), size = 3, color = "gold")
pplot

pdf(fig_6e_filename, width = 4, height = 4)
print(pplot)
dev.off()

df %>% filter(!extension %in% c("E00.5","E04")) %>% 
  select(pubid = subject_id, 
         visit = extension,
         value= NIAID_neuts) %>% 
  mutate(variable =  'NIAID_neuts') %>%
  write.csv(f6e_fdata, row.names = F)



## Fig. 4f-g
## compare Nt50 and diagnostic CD4 breadth
#f: 
x <- df %>%
  filter(extension == "E00") %>%
  select(subject_id, NIAID_neuts)
y <- df %>%
  filter(extension == "E00") %>%
  select(subject_id, breadth_classii)
test <- drop_na(merge(x,y, by = "subject_id"))
test <- filter(test, subject_id!="845")
test$breadth_classii[test$breadth_classii == 0] <- 0.000016700

a <- ggplot(test, aes(x = breadth_classii, y = NIAID_neuts)) + 
  geom_point() + scale_y_continuous(trans = "log2", breaks = c(20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480),
                                    labels = c("20", "40", "80", "160", "320", "640", "1280", "2560", "5120", "10240", "20480"), 
                                    limits = c(20,20480),
                                    name = "NT50 (E00)") +
  scale_x_continuous(limits = c(0.0002, 0.0018), 
                     name = "Diagnostic CD4 breadth (E00)",
                     labels = c('0.0005', '0.0010', '0.0015'),
                     breaks = c(0.0005, 0.001, 0.0015)) + 
  geom_smooth(method = "loess", span = 1, se = F, color = "grey") + 
  theme_classic() 
test_a = test

x <- df %>%
  filter(extension == "E02") %>%
  select(subject_id, NIAID_neuts)
y <- df %>%
  filter(extension == "E00") %>%
  select(subject_id, breadth_classii)
test <- drop_na(merge(x,y, by = "subject_id"))
test <- filter(test, subject_id!="845")
test$breadth_classii[test$breadth_classii == 0] <- 0.000016700

test_b = test
b <- ggplot(test, aes(x = breadth_classii, y = NIAID_neuts)) + 
  geom_point() + scale_y_continuous(trans = "log2", breaks = c(20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480),
                                    labels = c("20", "40", "80", "160", "320", "640", "1280", "2560", "5120", "10240", "20480"), 
                                    limits = c(20,20480),
                                    name = "NT50 (E02)") +
  scale_x_continuous(limits = c(0.0002, 0.0018), 
                     name = "Diagnostic CD4 breadth (E00)",
                     labels = c('0.0005', '0.0010', '0.0015'),
                     breaks = c(0.0005, 0.001, 0.0015)) + 
  geom_smooth(method = "loess", span = 1, se = F, color = "grey") + 
  theme_classic() 

both_plots <- plot_grid(a, b, ncol = 2)
both_plots
pdf(fig_6f_filename, width = 5, height = 2)
print(both_plots)
dev.off()

bind_rows(
  test_a %>% mutate(visit_neut = 'E00', visit_breadth = 'E00', variable_breadth = 'breadth_classii'),
  test_b %>% mutate(visit_neut = 'E00', visit_breadth = 'E02', variable_breadth = 'breadth_classii'),
) %>% 
  write.csv(f6f_fdata, row.names = F)



## Fig 4g - Neutralization at E00 vs E02
x <- df %>%
  filter(extension == "E00") %>%
  select(subject_id, NIAID_neuts)
y <- df %>%
  filter(extension == "E02") %>%
  select(subject_id, NIAID_neuts)
test <- drop_na(merge(x,y, by = "subject_id"))
test <- filter(test, subject_id!="845")
test_c = test
plot <- ggplot(test, aes(x = NIAID_neuts.x, y = NIAID_neuts.y)) + 
  geom_point() + scale_y_continuous(trans = "log2", breaks = c(20, 40, 80, 160, 320, 640, 1280, 2560, 5120, 10240, 20480),
                                    labels = c("20", "40", "80", "160", "320", "640", "1280", "2560", "5120", "10240", "20480"), 
                                    limits = c(20,20480),
                                    name = "NT50 (E02)") +
  scale_x_continuous(trans = "log2", breaks = c(20, 40, 80, 160, 320, 640, 1280),
                     labels = c("20", "40", "80", "160", "320", "640", "1280"), 
                     limits = c(20,1280),
                     name = "NT50 (E00)") +
  geom_smooth(method = "loess", span = 1, se = F, color = "grey") + 
  theme_classic() 
plot

pdf(fig_6g_filename, width = 3, height = 2)
print(plot)
dev.off()



df %>%
  filter(subject_id!="845") %>%
  filter(extension %in% c("E00","E02")) %>%
  select(pubid = subject_id, visit = extension,value = NIAID_neuts ) %>% 
  mutate(variable = 'NIAID_neuts') %>% 
  select(pubid, visit, value, variable) %>% 
  write.csv(f6g_fdata, row.names = F)
