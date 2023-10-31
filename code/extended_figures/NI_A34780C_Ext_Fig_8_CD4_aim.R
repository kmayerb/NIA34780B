# NI_A34780C_Ext_Fig_8.R
# Alternative CD4+ AIM markers
# 20231022

##################
##################
##################
# XNOTE: fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/2023_01_28_visualize.R
require(ggplot2)
require(dplyr)
# SPECIFY USER TO GET PREFERRED PATH
user = 'kosh'
if (user == 'emily'){
  path1 = '/fh/fast/corey_l/esford3_kmayerbl_collab/2023_tcrb/'
  path2 = '/fh/fast/corey_l/esford3_kmayerbl_collab/2023_tcrb/2023_tcrb/2023_tcrb_REVIEW/'
  #outpath_fig = "/fh/fast/corey_l/esford3_kmayerbl_collab/2023_tcrb/figures/"
  outpath_fig = "/fh/fast/corey_l/esford3_kmayerbl_collab/2023_tcrb/figures_REVIEW/"
}else if (user == 'kosh') {
  path1 = '/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/'
  path2 = '/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/2023_tcrb_REVIEW/'
  #outpath_fig = "/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/figures/"
  outpath_fig = "/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/figures_REVIEW/"
}
# CREATE A LIST TO STORE RESULTS BY PARTICIPANT

store = list()
ggst = list()
store_df = list()
store_dfnr = list()
store_df4 = list()
store_df37 = list()
store_df34 = list()
store_df3734 = list()


for (p in c('673', '525', '684', '754', '761', '782', '836')){
  # GET FILENAME AND CORRECT PATH
  f = paste0(path1,p,'_dfp_wide.tsv')
  # LOAD Data., replace zero with 1/1M, and look at enrichment factor
  df = read.csv(f, sep ="\t", stringsAsFactors = F) %>% 
    mutate(e01_pbmc = ifelse(e01_pbmc == 0, 1E-6, e01_pbmc),
           e03_pbmc = ifelse(e03_pbmc == 0, 1E-6, e03_pbmc))
  head(df[,c('pubid','total_counts_cd4')], 2)
  df[!is.na(df$total_counts_cd4),]$total_counts_cd4
  
  # IS TCRB ENRICHED IN CD4?
  df$enrich_cd4                   = df$e03_cd4 > df$e03_pbmc
  # IS TCRB ENRICHED IN CD69/CD137?
  df$enrich_cd37                  = df$e03_cd4_cd69_cd137     > df$e03_pbmc
  # IS TCRB ENRICHED IN CD134/CD154
  df$enrich_cd34                  = df$e03_cd4_cd154_or_cd134 > df$e03_pbmc
  # IS TCRB ENRICHED IN BOTH CD4, and CD69/CD137?
  df$enrich_cd4_AND_cd37          = df$enrich_cd4 & df$enrich_cd37  
  # IS TCRB ENRICHED IN BOTH CD4, and CD134|CD154
  df$enrich_cd4_AND_cd34          = df$enrich_cd4 & df$enrich_cd34  
  # IS TCRB ENRICHED IN CD4, and, CD69/CD137, CD134|CD154
  df$enrich_cd4_AND_cd37_AND_cd34 = (df$enrich_cd4_AND_cd37 ) & (df$enrich_cd4_AND_cd34)
  
  store[[p]] = data.frame( 
    pubid = p,
    num_cd4 = sum(df$enrich_cd4),
    num_cd4_AND_cd69137 = sum(df$enrich_cd4_AND_cd37), 
    num_cd4_AND_cd134OR154 = sum(df$enrich_cd4_AND_cd34), 
    num_cd4_AND_cd69_AND_cd134OR154 = sum(df$enrich_cd4_AND_cd37_AND_cd34),
    num_cd69137_no_cd4_filter =  sum( df$enrich_cd37 ), 
    num_cd134OR154_no_cd4_filter = sum( df$enrich_cd34) ,
    num_cd69137_AND_cd134OR154_no_cd4_filter = sum( df$enrich_cd37 & df$enrich_cd34)
  )
  
  
  
  # PLOTTING
  df4 = df %>% filter(enrich_cd4)
  df37 = df %>% filter(enrich_cd4_AND_cd37 )
  df34 = df %>% filter(enrich_cd4_AND_cd34)
  df3734 = df %>% filter(enrich_cd4_AND_cd37_AND_cd34)
  
  df_nr = df %>% group_by(e01_pbmc, 
                          e03_pbmc) %>% 
    summarise(count = n())
  store_df[[p]] = df
  store_dfnr[[p]] = df_nr
  store_df4[[p]] = df4 
  store_df37[[p]] = df37
  store_df34[[p]] = df34
  store_df3734[[p]] = df3734
  
  
  # REDUCE REDUNDANT VECTORS
  df_reduce = df %>% group_by(e01_pbmc, 
                              e03_pbmc) %>% 
    slice(1)
  
  df4_reduce = df4 %>% group_by(e01_pbmc, 
                                e03_pbmc) %>% 
    slice(1)
  require(scales)
  
  gg1 = ggplot(df4_reduce, aes(x = e01_pbmc, y = e03_pbmc))+ 
    geom_point(data = df_reduce, size = .5, col = "gray")+
    geom_point(data = df4_reduce, size= .5) + 
    geom_point(data = df37, pch = 1, size= 1, col = "green") + 
    geom_abline(aes(slope =1, intercept = 0), col = "gray", linetype = "dashed")+
    theme_classic() + 
    scale_y_continuous(trans='log10',
                       limits = c(1e-6, .1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x)))+
    scale_x_continuous(trans='log10',
                       limits = c(1e-6, .1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x)))+
    annotation_logticks(side = 'lb', size = .3) + 
    xlab("E01")+
    ylab("E03") + 
    ggtitle(paste0("PUBID: " , p, "\n(CD4+ AND CD169+ AND CD137+)"))+
    theme(plot.title = element_text(size=6))
  
  
  gg2 = ggplot(df4, aes(x = e01_pbmc, y = e03_pbmc))+
    geom_point(data = df_reduce, size = .5, col = "gray")+
    geom_point(data = df4_reduce, size= .5) + 
    geom_point(data = df34, pch = 1, size= 1, col = "violet") + 
    geom_point(data = df3734, pch = 1, size= 1, col = "orange") +
    geom_abline(aes(slope =1, intercept = 0), col = "gray", linetype = "dashed")+
    theme_classic() + 
    scale_y_continuous(trans='log10',
                       limits = c(1e-6, .1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x)))+
    scale_x_continuous(trans='log10',
                       limits = c(1e-6, .1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x)))+
    annotation_logticks(side = 'lb', size = .3) + 
    xlab("E01")+
    ylab("E03") + 
    ggtitle(paste0("PUBID: ", p, " \n(CD4+ AND CD134+ AND/OR CD154+)")) + 
    theme(plot.title = element_text(size=6)) 
  
  
  pdf(paste0(outpath_fig, p, "CD4_CONDITIONAL_CD137_and_CD13454.pdf"), width = 8, height = 4)
  gridExtra::grid.arrange(gg1,gg2, ncol = 2)
  dev.off()
  ggst[[p]][['cd37']] = gg1
  ggst[[p]][['cd34']] = gg2
  
}

### MAKE DATA FILES ###
store_dfnr_all   = list()
store_df4_all    = list()
store_df37_all   = list()
store_df34_all   = list()
store_df3734_all = list()

for (p in c('673', '525', '684', '754', '761', '782', '836')){
  print(p)
  
  store_dfnr_all[[p]]   = store_df[[p]] %>% 
    group_by(e01_pbmc, e03_pbmc) %>% 
    summarise(count=n()) %>% 
    mutate(data_type = "PBMC") %>%
    mutate(pubid = p) %>% 
    ungroup()
  
  store_df4_all[[p]] = store_df4[[p]] %>% 
    group_by(e01_pbmc, e03_pbmc) %>% 
    summarise(count=n()) %>% 
    mutate(data_type = "CD4+") %>%
    mutate(pubid = p) %>% 
    ungroup()
  
  store_df37_all[[p]]   =  store_df37[[p]] %>% 
    group_by(e01_pbmc, e03_pbmc) %>% 
    summarise(count=n()) %>% 
    mutate(data_type = "CD4+CD69+CD37+") %>%
    mutate(pubid = p) %>% 
    ungroup()
  
  store_df34_all[[p]]   = store_df34[[p]] %>% 
    group_by(e01_pbmc, e03_pbmc) %>% 
    summarise(count=n()) %>% 
    mutate(data_type = "CD4+CD69+/CD37+CD134+/CD154+") %>%
    mutate(pubid = p) %>% 
    ungroup()
  
  store_df3734_all[[p]] = store_df3734[[p]] %>% 
    group_by(e01_pbmc, e03_pbmc) %>% 
    summarise(count=n()) %>% 
    mutate(data_type = "CD4+CD69+CD137+ANDCD134+/CD154+") %>%
    mutate(pubid = p) %>% 
    ungroup()
  
}


store_dfnr_all_fdata   = do.call(rbind,store_dfnr_all)
store_df4_all_fdata    = do.call(rbind,store_df4_all)
store_df37_all_fdata   = do.call(rbind,store_df37_all)
store_df34_all_fdata   = do.call(rbind,store_df34_all)
store_df3734_all_fdata = do.call(rbind,store_df3734_all)

figdata_cd4_e01_e03_markers = bind_rows(
  store_dfnr_all_fdata,
  store_df4_all_fdata,
  store_df37_all_fdata, 
  store_df34_all_fdata,
  store_df3734_all_fdata)

figdata_cd4_e01_e03_markers %>% 
  ggplot(aes(x = e01_pbmc ,y =  e03_pbmc  , col = data_type)) + 
  geom_point(size = .001) + 
  scale_y_log10()+ 
  scale_x_log10() + 
  facet_wrap(~pubid)

# WRITE OUT FIGURE FILE
fig_9_fdata = paste0(repo_loc, repo, "figure_data_file/ex_fig9a_data.csv")
figdata_cd4_e01_e03_markers %>% 
  write.csv( file.path(path2 = '/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/2023_tcrb_REVIEW/ex_fig9a_data.csv'), row.names = F)
# COPY OVER MANUALLY FORM SERVER


# Extra simple theme
etheme = theme(axis.text = element_text(size = 6)) + 
  theme(axis.title = element_text(size = 6))
pdf(paste0(outpath_fig, "ALL_PTID_CD4_CONDITIONAL_CD137_and_CD13454_large.pdf"), width = 8, height =3.25*7)
print(gridExtra::grid.arrange(
  ggst[['525']][['cd37']]+etheme, ggst[['525']][['cd34']]+etheme,
  ggst[['673']][['cd37']]+etheme, ggst[['673']][['cd34']]+etheme,
  ggst[['684']][['cd37']]+etheme, ggst[['684']][['cd34']]+etheme,
  ggst[['754']][['cd37']]+etheme, ggst[['754']][['cd34']]+etheme,
  ggst[['761']][['cd37']]+etheme, ggst[['761']][['cd34']]+etheme,
  ggst[['782']][['cd37']]+etheme, ggst[['782']][['cd34']]+etheme,
  ggst[['836']][['cd37']]+etheme, ggst[['836']][['cd34']]+etheme,
  ncol = 2))
dev.off()


pdf(paste0(outpath_fig, "ALL_PTID_CD4_CONDITIONAL_CD137_and_CD13454_small.pdf"), width = 4.25, height =2*7)
print(gridExtra::grid.arrange(
  ggst[['525']][['cd37']]+etheme, ggst[['525']][['cd34']]+etheme,
  ggst[['673']][['cd37']]+etheme, ggst[['673']][['cd34']]+etheme,
  ggst[['684']][['cd37']]+etheme, ggst[['684']][['cd34']]+etheme,
  ggst[['754']][['cd37']]+etheme, ggst[['754']][['cd34']]+etheme,
  ggst[['761']][['cd37']]+etheme, ggst[['761']][['cd34']]+etheme,
  ggst[['782']][['cd37']]+etheme, ggst[['782']][['cd34']]+etheme,
  ggst[['836']][['cd37']]+etheme, ggst[['836']][['cd34']]+etheme,
  ncol = 2))
dev.off()


pdf(paste0(outpath_fig, "ALL_PTID_CD4_CONDITIONAL_CD137_and_CD13454_small_wide.pdf"), width = 8.5, height =2*4)
print(gridExtra::grid.arrange(
  ggst[['525']][['cd37']]+etheme, ggst[['525']][['cd34']]+etheme,  
  ggst[['673']][['cd37']]+etheme, ggst[['673']][['cd34']]+etheme,
  ggst[['684']][['cd37']]+etheme, ggst[['684']][['cd34']]+etheme,
  ggst[['754']][['cd37']]+etheme, ggst[['754']][['cd34']]+etheme,
  ggst[['761']][['cd37']]+etheme, ggst[['761']][['cd34']]+etheme,
  ggst[['782']][['cd37']]+etheme, ggst[['782']][['cd34']]+etheme,
  ggst[['836']][['cd37']]+etheme, ggst[['836']][['cd34']]+etheme,
  ncol = 4))
dev.off()
# COMBINE A LIST OF DATAFRAMES INTO A SINGLE DATA.FRAME
store_all = do.call(rbind, store)
store_all %>% write.csv()
store_all %>% readr::write_csv(paste0(path2, "2022_02_01_marker_summary.csv"))

##################
##################
##################

# XNOTE: fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/2023_03_06_visualize.R
# March 9, 2023
# fh/fast/gilbert_p/koell_kovid/2023_tcrb/2023_03_06_visualize.R
# Extended Data Figures 12
# We want to create and extended data figure that shows
# "classic" and "augmented" CD4+ AIM+ longitudinal trajectories.


require(readr)
require(dplyr)
require(ggplot2)
outpath = '/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb'
store = list()
store_fc = list()
# THIS FIRST BLOCK OPENS AND PROCESSES EACH FILE
# STORING in store_fc, fold changes of each clones for each of 7 ptids
# STORING in store a ggplot graph for each of 7 ptids
for (pubid in c('673','525','684','754','761','782','836')){
  fp = file.path(outpath, paste0(pubid,'.longitudinal_freq_cd4_activated.tsv') )
  df = readr::read_tsv(print(fp))
  #break
  # IS TCRB ENRICHED IN CD4?
  df$enrich_cd4                   = df$e03_cd4 > df$e03_pbmc
  # IS TCRB ENRICHED IN CD69/CD137?
  df$enrich_cd37                  = df$e03_cd4_cd69_cd137     > df$e03_pbmc
  # IS TCRB ENRICHED IN CD134/CD154
  df$enrich_cd34                  = df$e03_cd4_cd154_or_cd134 > df$e03_pbmc
  # IS TCRB ENRICHED IN BOTH CD4, and CD69/CD137?
  df$enrich_cd4_AND_cd37          = df$enrich_cd4 & df$enrich_cd37  
  # IS TCRB ENRICHED IN BOTH CD4, and CD134|CD154
  df$enrich_cd4_AND_cd34          = df$enrich_cd4 & df$enrich_cd34  
  # IS TCRB ENRICHED IN CD4, and, CD69/CD137, CD134|CD154
  df$enrich_cd4_AND_cd37_AND_cd34 = (df$enrich_cd4_AND_cd37 ) & (df$enrich_cd4_AND_cd34)
  require(scales)
  
  
  ns = names(df)[grepl(names(df), pattern ="_pfreq")]
  rns = ns %>% gsub(pattern = "_pfreq", replacement =  "")
  names(rns) = ns 
  
  fc_df = bind_rows(
    df %>% 
      filter(enrich_cd4) %>% 
      filter(enrich_cd37) %>% 
      select(  c('cdr3_b_nucseq', ns)) %>% 
      mutate(fc_e03 = (E03_pfreq + 1e-6) / (E01_pfreq + 1e-6) ) %>%
      mutate(fc_e02 = (E02_pfreq + 1e-6) / (E01_pfreq + 1e-6) ) %>% 
      mutate(type = "CD69CD137 Double Pos.\n")%>% 
      mutate(pubid = pubid) %>%
      select(type,  pubid, fc_e03 , fc_e02 ),
    df %>% filter(enrich_cd4) %>% 
      filter(enrich_cd34) %>% 
      select(  c('cdr3_b_nucseq', ns)) %>% 
      mutate(fc_e03 = (E03_pfreq + 1e-6) / (E01_pfreq + 1e-6) ) %>%
      mutate(fc_e02 = (E02_pfreq + 1e-6) / (E01_pfreq + 1e-6) ) %>% 
      mutate(type = "CD69CD137 Single Pos.\n CD134 and/or CD154") %>%
      mutate(pubid = pubid) %>%
      select(type,  pubid, fc_e03 , fc_e02 ))
  # STORE FOLD_CHANGE INFORMATION  
  store_fc[[pubid]] = fc_df
  
  classic = df %>% filter(enrich_cd4) %>% 
    filter(enrich_cd37) %>% 
    select(  c('cdr3_b_nucseq', ns)) %>%
    filter(E02_pfreq > 4*E01_pfreq | E03_pfreq > 4*E01_pfreq) %>%
    tidyr::gather(key, value, -cdr3_b_nucseq) %>% 
    mutate(key = rns[key]) %>%
    mutate(type = "CD69CD137 Double Pos.\n") %>% 
    mutate(type2 = "type1")
  
  additional =  df %>% filter(enrich_cd4) %>% 
    filter(enrich_cd34) %>% 
    select(  c('cdr3_b_nucseq', ns)) %>%
    filter(E02_pfreq > 4*E01_pfreq | E03_pfreq > 4*E01_pfreq) %>%
    tidyr::gather(key, value, -cdr3_b_nucseq) %>% 
    mutate(key = rns[key]) %>%
    mutate(type = "CD69CD137 Single Pos.\n CD134 and/or CD154") %>% 
    mutate(type2 = "type2" )
  
  both =  df %>% filter(enrich_cd4) %>% 
    filter(enrich_cd4_AND_cd37_AND_cd34) %>% 
    select(  c('cdr3_b_nucseq', ns)) %>%
    filter(E02_pfreq > 4*E01_pfreq | E03_pfreq > 4*E01_pfreq) %>%
    tidyr::gather(key, value, -cdr3_b_nucseq) %>% 
    mutate(key = rns[key]) %>%
    mutate(type = "both") %>% 
    mutate(type2 = "type2")
  
  
  ca = bind_rows(classic, additional, both) %>%
    mutate(type = factor(type,levels = c("CD69CD137 Double Pos.\n","CD69CD137 Single Pos.\n CD134 and/or CD154", "both")))
  
  levels(ca$type)
  
  gg= ggplot(ca, aes(x = key, y =value, col = type)) + 
    geom_point(size = .5) + 
    geom_line(aes(group = cdr3_b_nucseq), alpha = .2) + 
    theme_classic() + 
    scale_color_manual("",values = c("CD69CD137 Double Pos.\n"= "green",
                                     "CD69CD137 Single Pos.\n CD134 and/or CD154" = "violet",
                                     "both" = "orange"))+
    scale_y_continuous(trans='log10',
                       limits = c(1e-6, .1),
                       breaks=trans_breaks('log10', function(x) 10^x),
                       labels=trans_format('log10', math_format(10^.x)))+
    annotation_logticks(side = "l", size = .2)+ 
    theme(legend.position = "none")+
    facet_wrap(~type2) + 
    ggtitle(pubid)+ 
    ylab("") + xlab("") + 
    theme(strip.text = element_blank(), strip.background = element_blank()) + 
    theme(title = element_text(size = 7))+
    theme(axis.text.y = element_text(size = 7))+
    theme(axis.text.x = element_text(size = 7, angle = 0, hjust = 0, vjust = .5)) 
  
  store[[pubid]] =gg
}

#combine all fold change into a single data.frame 
fc_all = do.call(rbind,store_fc)
fc_all


# plot fold change distribution as boxplots
fc2 = ggplot(fc_all, aes(x =factor(pubid), fill = type, y = fc_e02)) + 
  geom_boxplot(outlier.size = .3, width = .5) +
  theme_classic() +
  scale_y_continuous(trans='log2',
                     limits = c(1e-2, 1e2),
                     breaks=c(0.0625,0.125, 0.25,0.5,0, 2, 4,8,16, 32, 64),#trans_breaks('log2', function(x) 2^x),
                     labels=trans_format('log2', math_format(2^.x)))+
  scale_fill_manual("",values = c("CD69CD137 Double Pos.\n"= "green",
                                  "CD69CD137 Single Pos.\n CD134 and/or CD154" = "violet"
  ))+
  ylab("log2 fold change\n(pre-vac vs. post-vac 1)") +
  xlab("") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 7))+
  theme(axis.text.x = element_text(size = 7))+
  theme(axis.title = element_text(size = 7)) 
fc2
fc3 = ggplot(fc_all, aes(x =factor(pubid), fill = type, y = fc_e03)) + 
  geom_boxplot(outlier.size = .3, width = .5) +
  theme_classic() +
  scale_y_continuous(trans='log2',
                     limits = c(1e-2, 1e2),
                     breaks=c(0.0625,0.125, 0.25,0.5,0, 2, 4,8,16, 32, 64),#trans_breaks('log2', function(x) 2^x),
                     labels=trans_format('log2', math_format(2^.x)))+
  #scale_y_continuous(trans='log10',
  #                   limits = c(1e-2, 1e2),
  #                   breaks=trans_breaks('log10', function(x) 10^x),
  #                   labels=trans_format('log10', math_format(10^.x)))+
  scale_fill_manual("",values = c("CD69CD137 Double Pos.\n"= "green",
                                  "CD69CD137 Single Pos.\n CD134 and/or CD154" = "violet"
  ))+
  ylab("log2 fold change\n(pre-vac vs. post-vac 2)") + 
  xlab("") +
  geom_hline(aes(yintercept = 1), linetype = "dashed") + 
  theme(legend.position = "none") + 
  theme(axis.text.y = element_text(size = 7))+
  theme(axis.text.x = element_text(size = 7))+
  theme(axis.title = element_text(size = 7)) 

#pdf('/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/figures/CD4_activated_all_4x_exp.pdf',
#    width = 8, height = 10)
lm = matrix(c(1,2,
              3,4,
              5,6,
              7,NA, 
              8,9,
              8,9), byrow = TRUE, ncol = 2)
gridExtra::grid.arrange(
  store[['525']]+theme(plot.title = element_text(hjust = 0.5)) ,
  store[['761']]+theme(plot.title = element_text(hjust = 0.5)),
  store[['673']]+theme(plot.title = element_text(hjust = 0.5)),
  store[['782']]+theme(plot.title = element_text(hjust = 0.5)),
  store[['684']]+theme(plot.title = element_text(hjust = 0.5)),
  store[['836']]+theme(plot.title = element_text(hjust = 0.5)),
  store[['754']]+theme(plot.title = element_text(hjust = 0.5)),
  fc2 + geom_hline(aes(yintercept = 4), size = .3, linetype = "dashed", col = "red"),
  fc3 + geom_hline(aes(yintercept = 4), size = .3, linetype = "dashed", col = "red"),
  layout_matrix = lm)
#dev.off()

# 
# 
# pdf('/fh/fast/gilbert_p/fg_data/koelle_covid/2023_tcrb/figures/CD4_activated_all_clones.pdf',
#     width = 8, height = 10)
# gridExtra::grid.arrange(
#              store[['525']] ,
#              store[['761']],
#              store[['673']],
#              store[['782']],
#              store[['684']],
#              store[['836']],
#              store[['754']],
#              ncol = 2)
# dev.off()


