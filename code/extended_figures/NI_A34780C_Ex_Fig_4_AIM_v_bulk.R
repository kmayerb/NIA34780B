# NI_A34780C_Ex_Fig_4_AIM_v_bulk.R

# Code for extended data figure 4a,b,c,d,e, showing AIM+ enrichment 
# for high confidence antigen specificity determination
# Xnote: 20231022 - ESF code
# Xnote: Previously was Ext 5, now Ext 4
## clear env
rm(list=ls())
# Set repo location as appropriate.
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
## OUTPUTS
ex_fig_5a_filename  = file.path(repo_loc,  repo, 'figures/Ex_Fig_4a_NI_A34780C.pdf')
ex_fig_5bc_filename = file.path(repo_loc,  repo, 'figures/Ex_Fig_4bc_NI_A34780C.pdf')
ex_fig_5de_filename = file.path(repo_loc,  repo, 'figures/Ex_Fig_4de_NI_A34780C.pdf')

### R code for Extended Figure 
#setwd('/Volumes/corey_l/esford3_kmayerbl_collab/software')

## PACKAGES
library(readr)
library(tidyverse)
library(dplyr)
library(tidyr)
library(readxl)
library(data.table)
library(grDevices)
library(gplots)
library(ggplot2)
library(viridis)
library(cowplot)
library(ggbeeswarm)

# INPUTS 
## phenotype file is an aggregation of all 10x CD4/CD8 citeseq typing results
phenotype_file = file.path(repo_loc, repo, 'data/aggregate_phenotypes_unique_101223.csv')
## all sig_e03 is an aggregation of all bulk expanded clonotypes (by TRB-seq)
sig_e03_file = file.path(repo_loc, repo, 'data/all_e03_sig.2.tsv')

## this set uses files where adaptive nucleotide sequence was joined by substring matching
## to the 10x nucleotide sequences
repo_loc2 = '/Volumes/fh/fast/corey_l/esford3_kmayerbl_collab/software' ## if  you're emily
repo2 = 'NIA34780B'
l525_file = file.path(repo_loc2, repo2, 'data/e03_10x_15525.csv')
l581_file = file.path(repo_loc2, repo2, 'data/e03_10x_15581.csv')



## EXTENDED DATA FIG 5a
all_10x <- read.csv(phenotype_file) %>%
  mutate(cell_type = ifelse((cell_type_2 == "CD4CD8" | cell_type_2 == "LOW_CD"), "not defined", cell_type_2)) %>%
  rename('cdr3_b_aa' = 'cdr3_b_aa_x') %>%
  mutate(above1 = ifelse(X10x_pfreq> productive_frequency, "enriched", "not_enriched")) %>%
  mutate(above = ifelse(X10x_cdf_p_fdr_10x < 0.05, "stat_enriched", "stat_not_enriched")) %>%
  mutate(enrich_type = ifelse(above1 == "not_enriched", 
                            'not_enriched', 
                            ifelse(above == "stat_not_enriched", "stat_not_enriched", "enriched"))) %>%
  select(cdr3_b_aa, above, ptid, X10x_pfreq, productive_frequency, cdr3_b_nt, cdr3_b_nucseq, 
         cell_type, enrich_type)
#limit to persons with E01-E03 and 10x 
list <- c('525', '527', '531', '577', '581', '669', '673', '684', '754', '761', '782', '836')
all_10x <- filter(all_10x, ptid %in% list)
all_10x = all_10x %>% mutate(ptid2 = paste0("P", ptid))
#all_10x$enrich_type == all_10x$en
## check 10x vs productive frequency 
ggplot(all_10x, aes(y = X10x_pfreq, x = productive_frequency, fill = cell_type, col= enrich_type)) + 
  geom_point(pch = 21) + 
  scale_y_log10() + 
  scale_x_log10() + 
  scale_fill_manual(values = c("CD4"="green", "CD8"="darkblue", "not definted"="gray"))+
  scale_color_manual(values = c("stat_not_enriched"="orange", "not_enriched"="red", "enriched"="black")) 

## confirm that adaptive and 10x sequences have equivalent nucleotide sequences (based on available sequence length)
## phenotype file is matched by v+cdr3+j key due to nucleotide output being different lengths
colnames(all_10x)
match_vector <- rep(FALSE, nrow(all_10x))
# Loop through the rows of mismatches and check for substring matches
for (i in 1:nrow(all_10x)) {
  match_vector[i] <- grepl(all_10x$cdr3_b_nt[i], all_10x$cdr3_b_nucseq[i])
}
table(match_vector) ## 4304 TRUE merging is exact for each clonotype

## graph 
cell_color <- c("CD4" = "#009966", "CD8" = "#333399", "not defined" = "purple")
table(all_10x$cell_type)
a <- ggplot(data = all_10x, aes(x = productive_frequency, y = X10x_pfreq, group = cell_type)) +
  geom_abline(linetype = 'dashed') + 
  geom_point(data = subset(all_10x, enrich_type == "not_enriched"), aes(fill = cell_type), color = "red", alpha = 0.3, size = .5, shape = 20, show.legend = F) + 
  geom_point(data = subset(all_10x, enrich_type == "stat_not_enriched"), aes(fill = cell_type), color = "orange", alpha = 0.3, size = .5, shape = 20, show.legend = F) + 
  geom_point(data = subset(all_10x, enrich_type == "enriched"), aes(fill = cell_type), shape = 21, color = "black", 
             stroke = 0, color = NA, size = 1, show.legend = T, alpha = .5) + 
  theme_bw() + 
  scale_fill_manual(values = cell_color, name = "cell type") + 
                    #labels = c('CD4+TS', 'CD4+TS', 
                    #'No sig. enrichment', 'No enrichment')) +
  labs(x="E03 TRB clonotype (frequency)", y="E03 AIM-scTCRseq (frequency)") + 
  scale_y_log10(breaks = c(1E-1, 1E-2, 1E-3, 1E-4),
                limits = c(1E-4, 3E-1),
                labels = c(expression(~'10'^'-1'), expression(~'10'^'-2'), 
                           expression(~'10'^'-3'), expression(~'10'^'-4'))) + 
  scale_x_log10(breaks = c(1E-1, 1E-2, 1E-3, 1E-4, 1E-5, 1E-6),
                limits = c(1E-6, 1E-1),
                labels = c(expression(~'10'^'-1'), expression(~'10'^'-2'), 
                           expression(~'10'^'-3'), expression(~'10'^'-4'), expression(~'10'^'-5'),  "ND")) + 
  annotation_logticks(sides = "bl", size = .2) + 
  theme(title = element_text(size = 8), 
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 8),
        axis.title = element_text(size = 8),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.line = element_line(color = "black")) + 
  facet_wrap(~ptid2, scales = "free", ncol = 4) +
  guides(shape = guide_legend(override.aes = list(size = 5))) + 
  theme(legend.position = "top")


pdf(ex_fig_5a_filename, width = 6.75, height = 7.5)
a
dev.off()

#__________________________________________________________________________________________


## EXTENDED DATA FIG 5b,c 
## load in 10x results to re-match against TRB clonotypes by cdr3_b_nt sequence 
## with additional phenotype data
key <- read.csv(phenotype_file) %>%
  filter(`X10x_enriched` == TRUE) %>%
  select(c("cdr3_b_nt", "X10x_count", "X10x_pfreq", "ptid", "cell_type_2"))

e525 <- readr::read_csv(l525_file) 
k = e525$X15525_3_TCRB_E01 
min(k[k>0 & !is.na(k)]) ## 1.8e-06
k = e525$X15525_5_TCRB_E03 
min(k[k>0 & !is.na(k)]) ## 2.2e-06

df <- e525 %>%
  mutate_all(~ifelse(.=="", NA ,.)) %>%
  select(c("X", "X15525_3_TCRB_E01", "X15525_5_TCRB_E03", "cdr3_b_nt")) %>%
  left_join(subset(key, ptid == "525"), by = "cdr3_b_nt") %>% ## join to pre-matched 10x sequence
  filter(!duplicated(rleidv(., cols = c("X15525_3_TCRB_E01", "X15525_5_TCRB_E03", 'cdr3_b_nt')))) %>% ## make this more feasible to graph
  mutate(e01 = ifelse(X15525_3_TCRB_E01 == 0, 1e-6, X15525_3_TCRB_E01)) %>% ## set non-detect value below min
  mutate(e03 = ifelse(X15525_5_TCRB_E03 == 0, 1e-6, X15525_5_TCRB_E03)) %>%
  select(e01, e03, X10x_count, X10x_pfreq, cdr3_b_nt)

table(!duplicated(e525$cdr3_b_nt)) ##190 
table(!duplicated(df$cdr3_b_nt)) ##190 (make sure no overlappers are lost)

b <- ggplot(df, aes(x = e01, y = e03)) +
  geom_point(data = df %>% filter(is.na(X10x_count)), alpha = 0.2, size=2, shape=18, color='black', show.legend = F) +
  geom_point(data = df %>% filter(!is.na(X10x_count)), aes(color = log10(X10x_pfreq)), size=3, shape = 18, show.legend = T)  +
  labs(title = "P525") + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  scale_color_viridis(name = expression(AIM^"+"~" frequency"~log[10]), option = "plasma") +
  labs(x="E01", y="E03") + 
  scale_y_log10(breaks = c(10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-6'), expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) + 
  scale_x_log10(breaks = c(10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-6'), expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) + 
  annotation_logticks(sides = "bl") + 
  coord_cartesian(xlim = c(1E-6, 1E-1), ylim = c(1E-6, 1E-1))
b

e581 <- readr::read_csv(l581_file) 
k = e581$X15581_7_TCRB_E01 
min(k[k>0 & !is.na(k)]) ## 4.7e-06
k = e581$X15581_9_TCRB_E03 
min(k[k>0 & !is.na(k)]) ## 2.4e-06
remove(k)

df <- e581 %>%
  mutate_all(~ifelse(.=="", NA ,.)) %>%
  select(c("X", "X15581_7_TCRB_E01", "X15581_9_TCRB_E03", "cdr3_b_nt")) %>%
  left_join(subset(key, ptid == "581"), by = "cdr3_b_nt") %>% ## join to pre-matched 10x sequence
  filter(!duplicated(rleidv(., cols = c("X15581_7_TCRB_E01", "X15581_9_TCRB_E03", 'cdr3_b_nt')))) %>% ## make this more feasible to graph
  mutate(e01 = ifelse(X15581_7_TCRB_E01 == 0, 1e-6, X15581_7_TCRB_E01)) %>% ## set non-detect value below min
  mutate(e03 = ifelse(X15581_9_TCRB_E03 == 0, 1e-6, X15581_9_TCRB_E03)) %>%
  select(e01, e03, X10x_count, X10x_pfreq, cdr3_b_nt)

table(!duplicated(e581$cdr3_b_nt)) ##339 
table(!duplicated(df$cdr3_b_nt)) ##339 (make sure no overlappers are lost)

c <- ggplot(df, aes(x = e01, y = e03)) +
  geom_point(data = df %>% filter(is.na(X10x_count)), alpha = 0.2, size=2, shape=18, color='black', show.legend = F) +
  geom_point(data = df %>% filter(!is.na(X10x_count)), aes(color = log10(X10x_pfreq)), size=3, shape = 18, show.legend = T)  +
  labs(title = "P581") + theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5, size = 12)) +
  scale_color_viridis(name = expression(AIM^"+"~" frequency"~log[10]), option = "plasma") +
  labs(x="E01", y="E03") + 
  scale_y_log10(breaks = c(10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-6'), expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) + 
  scale_x_log10(breaks = c(10e-7, 10e-6, 10e-5, 10e-4, 10e-3, 10e-2),
                limits = c(1E-6, 1E-1),
                labels = c("ND", expression(~'10'^'-6'), expression(~'10'^'-5'), expression(~'10'^'-4'),  expression(~'10'^'-3'),  expression(~'10'^'-2'))) + 
  annotation_logticks(sides = "bl") + 
  coord_cartesian(xlim = c(1E-6, 1E-1), ylim = c(1E-6, 1E-1))
c

plots <- plot_grid(b, c, 
                   align = "v",
                   axis = "b",
                   nrow = 2,
                   ncol = 1)

pdf(ex_fig_5bc_filename, width = 5, height = 7)
print(plots)
dev.off()

#__________________________________________________________________________________________
## EXTENDED DATA FIG 5d, e
## question - of the significantly expanded clonotypes, does frequency account for why some were captured in AIM vs not
## start with the 10x-bulk matched phenotype file and add in the unmatched sig_expanded clonotypes as comparison - 
## match by key to avoid non-cdr3-region nucleotide mismatching 

key <- read.csv(phenotype_file) %>%
  filter(`X10x_enriched` == TRUE) %>%
  mutate(key = paste0(v_b_gene_y, "+", cdr3_b_aa_y, "+", j_b_gene_y)) %>% 
  select(c("key", "cdr3_b_nucseq", "X10x_count", "ptid"))
list <- unique(key$ptid) #N = 17

all <-read_tsv(sig_e03_file) %>% 
  filter(E03_vac_class == "sig_expand") %>%
  mutate(ptid = sub("^\\d{2}", "", .$ptid)) %>%
  filter(ptid %in% list) %>% 
  mutate(ptid = as.numeric(ptid)) %>%
  mutate(key = paste0(gsub("\\*01", "", v_b_gene), "+", cdr3_b_aa, "+", gsub("\\*01", "", j_b_gene))) %>%
  select(c("ptid", 'key', "cdr3_b_nucseq", "predetected", "E03_pfreq", "E03_vac_fc"))

count_color <- c("TRUE" = "#2D708EFF", "FALSE" = "#C2DF23FF")

all_merge <- left_join(all, key, by = c("cdr3_b_nucseq", "ptid")) 
## fewer - not all of the clonotypes in the phenotype file are sig_exp 

all_merge$detected <- !is.na(all_merge$X10x_count)
table(all_merge$detected)

d <- ggplot(data = all_merge, aes(E03_pfreq, group = detected, fill = detected)) + 
  geom_density(aes(fill = detected), alpha = 0.5) + 
  labs(x = 'Frequency of detected clonotypes', y="Fraction detected clonotypes") + 
  scale_y_continuous() + 
  scale_x_log10(breaks = c(1E-1, 1E-2, 1E-3, 1E-4, 1E-5),
                limits = c(1E-5, 1E-1),
                labels = c(expression(~'10'^'-1'), expression(~'10'^'-2'), expression(~'10'^'-3'), expression(~'10'^'-4'), expression(~'10'^'-5'))) + 
  scale_fill_manual(values = count_color, name = "Detected \n by AIM-scTCRseq") +   theme_classic()+
  theme(plot.title = element_text(size = 8, hjust = 0.5), legend.position = c(0.8,0.7))+
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 8))

e <- ggplot(data = all_merge, aes(E03_vac_fc, group = detected, fill = detected)) + 
  geom_density(aes(fill = detected), alpha = 0.5) + 
  labs(x = 'Fold change after primary vaccination', y="Fraction detected clonotypes") + 
  scale_y_continuous(limits = c(0,1)) + 
  scale_x_log10(breaks = c(1E-2, 1E-1, 1E1, 1E2, 1E3, 1E4),
                limits = c(0.01, 1E4),
                labels = c(expression(~'10'^'-2'), expression(~'10'^'-1'), expression(~'10'^'1'), expression(~'10'^'2'), expression(~'10'^'3'), expression(~'10'^'4'))) + 
  scale_fill_manual(values = count_color, guide = NULL) +   theme_classic()+
  theme(plot.title = element_text(size = 8, hjust = 0.5)) + 
  theme(axis.text = element_text(size = 8))+
  theme(axis.title = element_text(size = 8))  


plots <- plot_grid(d, e,
                   align = "v",
                   axis = "b",
                   nrow = 2,
                   ncol = 1)
plots

pdf(ex_fig_5de_filename, width = 2, height = 3)
gridExtra::grid.arrange(d,e)
dev.off()

