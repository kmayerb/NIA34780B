# NI_A34780B_Figure_3_Trajectory.R

# Xnote: Prev NI_A34780B_Extended_Figure_7.R # 
# Trajectory Analysis 
# August 31, 2022
# Code used for trajectory analysis.
require(ggplot2)
require(dplyr)
library(scales)
require(cluster)
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
set.seed(0623)
# Added TCRa information  
filename_f3 = file.path(
  repo_loc, repo,
  'data/2022_08_31_all_long_data.tsv')

# Added TCRa columns
#f1 = 'plot_data/2022_06_02_all_long_data.tsv'
#f1 = 'plot_data/2022_08_31_all_long_data.tsv'

# Outputs
global_font_size = 8
figure_3a_filename = file.path(repo_loc, repo, 'figures/Fig_3a_NI_A34780C.pdf')
figure_3b_filename = file.path(repo_loc, repo, 'figures/Fig_3b_NI_A34780C.pdf')
figure_3e_filename = file.path(repo_loc, repo, 'figures/Fig_3e_NI_A34780C.pdf')
f3a_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_3a_data.csv')
f3b_fdata=file.path(repo_loc,repo, 'figure_data_files/NI_A34780C_fig_3b_data.csv')


# Script
d = readr::read_tsv(filename_f3) %>% 
  filter(cell_type != "LOW_CD") %>% 
  filter(x_pos  %in% c("E00","E01","E02","E03","E05"))


names(d)
table(d$cell_type)
d %>% filter(cell_type!= "CD4CD8") %>% 
  ggplot(., aes(x = x_pos, y = value, col = cell_type))+ 
  geom_line(aes(group = plot_clone_id), size = .1, alpha = .5) + 
  #geom_line(aes(group = ptid), data = ds, size = .75, alpha = .9)+
  facet_grid(ptid~cell_type) + 
  scale_y_log10(expand = c(.05,0), labels = trans_format("log10", math_format(10^.x))) +
  #scale_color_viridis_d(option = "inferno", begin = 0, end = .7) +
  #scale_y_log10(expand = c(.05,0), breaks = trans_breaks("log10", function(x) 10^x), labels = trans_format("log10", math_format(10^.x)))+
  scale_color_manual(values = c("darkgreen", "blue")) + 
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .1)+
  theme(axis.text.y = element_text(size = 7)) + 
  theme(axis.text.x = element_text(size = 7, angle = 90)) + 
  #theme(axis.line.y = element_blank()) + 
  theme(strip.text = element_text(size = 8)) + 
  theme(strip.background = element_blank()) + 
  xlab("")+ylab("") + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2)  

# construction of matrix of CD4 and CD8 T cells.
d %>% filter(cell_type!= "CD4CD8") %>% 
  mutate(plot_clone_id_ptid = paste0(plot_clone_id, "_",ptid))%>%
  mutate(value = ifelse(value<=0, 1E-6, value)) %>%
  mutate(log10_value = log10(value)) %>% 
  filter(x_pos %in% c("E00", "E01","E02","E03","E05")) %>% 
  select(x_pos, log10_value, plot_clone_id_ptid, ptid, cell_type,
         cdr3_b_nucseq, cdr3_b_aa, v_b_gene, j_b_gene, 
         cdr3_a_aa, v_a_gene, j_a_gene) %>% 
  tidyr::spread(key = x_pos, value = log10_value) %>% 
  filter(complete.cases(.)) %>% 
  select(c("E00", "E01","E02","E03","E05")) -> m

d %>% filter(cell_type!= "CD4CD8") %>% 
  mutate(plot_clone_id_ptid = paste0(plot_clone_id, "_",ptid))%>%
  mutate(value = ifelse(value<=0, 1E-6, value)) %>%
  mutate(log10_value = log10(value)) %>% 
  filter(x_pos %in% c("E00", "E01","E02","E03","E05")) %>% 
  select(x_pos, log10_value, plot_clone_id_ptid, ptid, cell_type, 
         cdr3_b_nucseq, cdr3_b_aa, v_b_gene, j_b_gene,
         cdr3_a_aa, v_a_gene, j_a_gene) %>% 
  tidyr::spread(key = x_pos, value = log10_value) %>% 
  filter(complete.cases(.)) %>% 
  select(-c("E00", "E01","E02","E03","E05")) -> meta

# check the dimensions
stopifnot(length(unique(d$ptid)) == 17)
stopifnot(length(unique(meta$ptid)) == 11)

# Construct Matrices
# <m> E00,01,02,03,05
m = as.matrix(m)
# <mdiff> - log difference rowwise
mdiff = t(apply(m, 1, function(x) diff(x)))
# <cosine> - matrix diff
cos_mat      = lsa::cosine(t(m))
cos_mat_diff = lsa::cosine(t(mdiff))
# convert to distances
cos_dist      = as.dist(1- cos_mat)
cos_dist_diff = as.dist(1- cos_mat_diff)
# hierarchical clusters
cos_dist_hc = hclust(cos_dist, method = "ward.D2")
cos_dist_diff_hc = hclust(cos_dist_diff, method = "ward.D2")
# viz 
plot(cos_dist_hc, labels = FALSE)
plot(cos_dist_diff_hc, labels = FALSE)
# gap statistic
hclusCut <- function(x, k, cos_mat){
  print(k)
  cos_dist=as.dist(1- cos_mat)
  hc = hclust(cos_dist, method="ward.D2")
  return(list(cluster = cutree(hc, k=k) ))
}
# Boot-strap the gap statistic
gsk      <- clusGap(m, FUN = hclusCut, cos_mat = cos_mat, K.max = 15, B = 10)
gsk_diff <- clusGap(mdiff, FUN = hclusCut, cos_mat = cos_mat_diff, K.max = 15, B = 10)
plot(gsk)
plot(gsk_diff)

# ggplot gap statistic m
# --- #
gap = as.data.frame(gsk$Tab)
gap$x <- seq_along(gap$gap)
ggplot(gap, aes(x = x, y = gap)) + 
  geom_point(size = .7) + 
  theme_bw() + 
  geom_segment(aes(x = x, y = gap, yend = gap+2*SE.sim, xend= x))+
  geom_segment(aes(x = x, y = gap, yend = gap-2*SE.sim, xend= x))+
  scale_x_continuous(breaks = seq(1,15)) + 
  xlab("Cluster Number (K)") + 
  ylab(expression(Gap[k]))

# ggplot gap statistic mdiff
# --- #
gap_diff = as.data.frame(gsk_diff$Tab)
gap_diff$x <- seq_along(gap$gap)

pdf(figure_3e_filename, width= 3, height =3)
ggplot(gap_diff %>% filter(x < 11), aes(x = x, y = gap)) + 
  geom_point(size = .7) + 
  theme_bw() + 
  geom_segment(aes(x = x, y = gap, yend = gap+2*SE.sim, xend= x))+
  geom_segment(aes(x = x, y = gap, yend = gap-2*SE.sim, xend= x))+
  scale_x_continuous(breaks = seq(1,15)) + 
  xlab("Cluster Number (K)") + 
  ylab(expression(Gap[k])) + 
  theme_classic()
dev.off()
## CUSTOM PLOTTING FOR K=5 ##
plot(gsk)
meta$cluster = cutree(cos_dist_hc,k = 5)
meta$cluster_global5 = cutree(cos_dist_hc,k = 5)

plot(gsk_diff)
meta$cluster_diff6 = cutree(cos_dist_diff_hc,k  = 6)
meta$cluster_diff11 = cutree(cos_dist_diff_hc,k = 11)
# THIS CAN BE USED FOR FURTHER EXTERNAL ANALYSIS 

# MAKE SUPORTING FIGURE
df      = as.data.frame(m)
meta
df2     = cbind(meta,df)
# STORE DIFFERENCE BASED CLUSTERING
#df2diff = as.data.frame(mdiff)
#df2diff = cbind(meta,df2diff)
#df2diff %>% write.csv("umap_input_diff.csv")

names(df2)
dim(df2)
table(df2$cell_type, df2$cluster)

# MAKE LONG FORM FOR PLOTTING
df3 = df2 %>% 
  tidyr::gather(key = "x_pos", value= "log10_value",
                -plot_clone_id_ptid, 
                -ptid, 
                -cell_type, 
                -cluster, 
                -cluster_global5,
                -cluster_diff6, 
                -cluster_diff11, 
                -cdr3_b_nucseq, 
                -cdr3_b_aa, 
                -v_b_gene, 
                -j_b_gene,
                -cdr3_a_aa,
                -v_a_gene,
                -j_a_gene) 
  


#COMPUTE MEAN
df3 %>% filter(cell_type!= "CD4CD8") %>% 
  group_by(cell_type, cluster, x_pos) %>% 
  mutate(value = 10^log10_value) %>%
  summarize(value = mean(value)) -> df3_mean

cluster_order=c('1','4','5','3','2')
t1 = df3 %>% filter(cell_type!= "CD4CD8") %>% 
  mutate(value = 10^log10_value) %>%
  ggplot(., aes(x = x_pos, y = value, col =forcats::fct_relevel(factor(cluster),cluster_order  ))) + 
  geom_point(size = .05, alpha = .4) + 
  geom_line(aes(group = plot_clone_id_ptid), size = .1, alpha = .25) + 
  geom_line(data = df3_mean, aes(group = cluster), size = 1, alpha = 1, color = "white") + 
  geom_line(data = df3_mean, aes(group = cluster), size = .5, alpha = 1, linetype = 'dashed') + 
  facet_grid(cell_type~forcats::fct_relevel(factor(cluster), cluster_order ), scales= "free") + 
  scale_y_log10(expand = c(.05,0), labels = trans_format("log10", math_format(10^.x)),
                limits= c(1E-6, 1)) +
  scale_color_viridis_d(option = "inferno", begin = 0, end = .7) +
  theme_classic()+
  theme(legend.position = "none") + 
  annotation_logticks(side = "l", size = .2)+
  theme(axis.text.y = element_text(size = 8)) + 
  theme(axis.text.x = element_text(size = 8, angle = 90)) + 
  theme(axis.line = element_line(size = .2))+
  theme(strip.text = element_text(size = 8)) + 
  theme(strip.background = element_blank()) + 
  xlab("")+ylab("") + 
  geom_vline(aes(xintercept = "E02"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E03"), linetype = "dashed", size =.2)+
  geom_vline(aes(xintercept = "E05"), linetype = "dashed", size =.2)+
  theme(panel.spacing = unit(.5, "lines"))

new <- unique(df3$ptid) %>% sub(pattern = "^15", replace = "P")
names(new) = unique(df3$ptid)

t2 = df3 %>% filter(cell_type!= "CD4CD8") %>% 
  mutate(cell_type= c("CD4"="Ts CD4+","CD8"= "Ts CD8+")[cell_type]) %>%
  ggplot(., aes(x = cell_type, fill = forcats::fct_relevel(factor(cluster), cluster_order ))) + 
  geom_bar(size = .05, alpha = .9, stat = "count", position = "fill") + 
  facet_wrap(~ptid,#scales= "free", 
             labeller = labeller(ptid = new), nrow = 1)+
  scale_y_continuous(expand = c(0,0),labels = scales::percent_format()) +
  scale_fill_viridis_d(option = "inferno", begin = 0, end = .7) +
  theme_classic() +
  theme(strip.background =  element_blank(),
        strip.text.y = element_blank(), 
        strip.text.x = element_text(size = 8)) +
  theme(legend.position = "none") + 
  theme(axis.title = element_blank()) + 
  theme(axis.text.x = element_text(size = 8)) + 
  theme(axis.text.x = element_text(size = 8, angle = 90)) + 
  theme(panel.spacing = unit(0, "lines"))

df3 %>%  filter(cell_type!= "CD4CD8") %>%
  group_by(cell_type, cluster, x_pos) %>% 
  tally() %>% filter(x_pos== 'E03', cell_type == "CD4") %>%
  write.csv()
# 1112/(3+1112+308+28)
# 308/(3+1112+308+28)

pdf(figure_3a_filename, width = 7, height = 2)
t1
dev.off()
pdf(figure_3b_filename, width = 6.5, height = 2)
t2
dev.off()

# WRITE FIGURE DATA
cluster_order=c('1'="1",'4'="2",'5'="3",'3'="4",'2'="5")
t1$data$group = cluster_order[t1$data$cluster]
t1$data %>% 
  select(pubid = ptid,
         cell_type,
         group, 
         visit = x_pos,
         log10_value,
         cdr3_a_aa,
         v_a_gene,
         j_a_gene,
         cdr3_b_aa,
         v_b_gene,
         j_b_gene,
         cdr3_b_nucseq) %>% 
  mutate(pubid = stringr::str_replace(pubid, pattern = "^15", "")) %>% 
  write.csv(f3a_fdata, row.names = F)

t2$data$group = cluster_order[t2$data$cluster]
t2$data %>% group_by(group, cell_type, ptid, x_pos) %>% 
  tally() %>% 
  filter(x_pos =="E03") %>%
  select(pubid = ptid, 
         cell_type,
         n = n) %>% 
  mutate(pubid = stringr::str_replace(pubid, pattern = "^15", "")) %>% 
  write.csv(f3b_fdata, row.names = F)
  









## Output clonotype trajectory cluster
# df2 %>% write.csv("2022_08_31_trajectory_clusters_with_clonoytpes.csv")

