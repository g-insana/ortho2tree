#!/usr/bin/env Rscript --vanilla

################
##
## 24-Nov-2023
##
## clade_stats.R -f output_gc231004
##
## plot proposed difference lengths by taxon and canonical db type
##
## reads in output_gc231004 file with cummulative changes
##
## 03-Nov-2023
## uses standard strategy for --pub supression of command line doc
##

library('ggplot2')
library('getopt')
library('optparse')
library('stringr')
library('RColorBrewer')
library('reshape2')

library(patchwork)

options(width=500)

option_list = list(
	    make_option(c("-f","--file"),type='character',action='store',help='output_gc file REQUIRED'),
	    make_option(c("-p","--pdf"),type='character',action='store',default='fig7cp', help='optional pdf file name [default= \"%default\"]'),
	    make_option(c("-P","--pub"),action='store_true',default=FALSE, help='remove plot documentation')
	    )


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt$pdf)

if (is.null(opt$file)) {
    print_help(opt_parser)
    cat(" clade_stats.R -f clade_cnts.tab\n\n")
    quit(save='no',status=1)
}

rel_map = c('2022_05'="2023_01 (8)",'2023_01'="2023_02 (13)",'2023_02'="2023_03 (13)",'2023_03'="2023_04 (13*)",'2023_04'="2023_05 (13**)",'2023_05'="2024_01 (35)")

print(names(rel_map))

shape_map = setNames(c(1, rep(2,4), 0),names(rel_map))

## globals
def_shape=2 ## triangle

args<- commandArgs(trailingOnly=TRUE)

p.name<-basename(get_Rscript_filename())
plabel=paste(c(p.name,"\n",args,"\n",date()),collapse=' ', sep=' ')

plot_name = "fig7cp"
if (! is.null(opt$pdf)) {
   plot_name = opt$pdf
}

################
## names for output plots and text files

fig_file=paste0(plot_name,"_clade_size")
if (opt$pub) {
   fig_file = paste0(fig_file,"_pub")
}
fig_file = paste0(fig_file,".pdf")

fig_file

arg_files<-strsplit(opt$file,',',fixed=TRUE)[[1]]

clade_df = NULL
rel_list = c()

ix <- 1
for (tab_file in arg_files) {

    rel_name = str_extract(tab_file,'[0-9]+_[0-9]+')

    print(sprintf("%s : %s : %s",tab_file, rel_name, rel_map[rel_name]))

    rel_list = append(rel_list,rel_name)

    tmp_df <- read.table(tab_file,header=TRUE,sep='\t')
    ## colnames(tmp_df) = c('pthr_id','n_taxa','MANEstatus','c_type')
    tmp_df$rel_name = rel_name

    clade_df = rbind(clade_df, tmp_df)

    ix = ix + 1
}

clade_df$c_type <- factor(clade_df$c_type,levels=c('conf','prop'),ordered=TRUE)
clade_df$rel_name <- factor(clade_df$rel_name)

print("All clades")
aggregate(clade_df, n_taxa ~ rel_name + c_type, function(x) quantile(x,probs=seq(0,1,0.1)))

n_by_rel <- aggregate(clade_df, n_taxa ~ rel_name + c_type, NROW)

print("MANE_good clades")
aggregate(clade_df[clade_df$MANEstatus=='MANE_good',], n_taxa ~ rel_name + c_type, function(x) quantile(x,probs=seq(0,1,0.1)))

n_by_rel_M <- aggregate(clade_df[clade_df$MANEstatus=='MANE_good',], n_taxa ~ rel_name + c_type, NROW)

clade_N <- merge(n_by_rel, n_by_rel_M,by=c('rel_name','c_type'))

colnames(clade_N) = c('rel_name','c_type','N_all','N_mane')
clade_N

## colnames(clade_N) = c('rel_name','N_all','N_mane')

rel_label_conf <- sprintf("%s %s",rel_map,format(clade_N[clade_N$c_type=='conf',]$N_all,big.mark=','))
rel_label_prop <- sprintf("%s %s",rel_map,format(clade_N[clade_N$c_type=='prop',]$N_all,big.mark=','))

rel_label_conf_m <- sprintf("%s %s",rel_map,format(clade_N[clade_N$c_type=='conf',]$N_mane,big.mark=','))
rel_label_prop_m <- sprintf("%s %s",rel_map,format(clade_N[clade_N$c_type=='prop',]$N_mane,big.mark=','))

## set the colors : enough for the number of releases
a.colors <- brewer.pal(7,'Dark2') # 'Dark2', 'Set2', 'Paired'
a.colors[6] <- a.colors[7]
a.colors <- setNames(a.colors, names(rel_map))


s.color_conf <- scale_color_manual(values=a.colors, labels=rel_label_conf)
s.color_prop <- scale_color_manual(values=a.colors, labels=rel_label_prop)
s.color_conf_m <- scale_color_manual(values=a.colors, labels=rel_label_conf_m)
s.color_prop_m <- scale_color_manual(values=a.colors, labels=rel_label_prop_m)

s.shape_conf <- scale_shape_manual(values=shape_map, labels=rel_label_conf)
s.shape_prop <- scale_shape_manual(values=shape_map, labels=rel_label_prop)
s.shape_conf_m <- scale_shape_manual(values=shape_map, labels=rel_label_conf_m)
s.shape_prop_m <- scale_shape_manual(values=shape_map, labels=rel_label_prop_m)

y_lab8k=c(0,seq(2000,8000,2000))

y_lab4k=c(0,seq(1000,4000,1000))
y_lab4k_lr=c("0    ","1,000 ","2,000","3,000","4,000")

y_lab2k=c(0,seq(500,2000,500))
y_lab2k_lr=c("0    ","500 ","1,000","1,500","2,000")

s.y.8K <- scale_y_continuous(limits=c(0,8000), breaks=y_lab8k,labels=format(y_lab8k,big.mark=','))

s.y.4K <- scale_y_continuous(limits=c(0,4000), breaks=y_lab4k,labels=format(y_lab4k,big.mark=','))
s.y.4Kr <- scale_y_continuous(limits=c(0,4000), breaks=y_lab4k,labels=y_lab4k_lr,position='right')

s.y.2K <- scale_y_continuous(limits=c(0,2000),breaks=y_lab2k,labels=format(y_lab2k,big.mark=','))
s.y.2Kr <- scale_y_continuous(limits=c(0,2000),breaks=y_lab2k,labels=y_lab2k_lr,position='right')

theme_set(theme_linedraw(base_size=14))
theme.base <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0),
  plot.caption=element_text(size=6,hjust=0),
  axis.text.x = element_text(size=12),
  legend.text=element_text(size=9),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.text.align = 0,
  legend.title=element_blank(),
  legend.margin=margin(c(0.2,5,5))
)

theme.leg = theme.base + theme(legend.position=c(0.98,0.98),legend.justification=c('right','top'))
theme.leg_m = theme.base + theme(legend.position=c(0.98,0.98),legend.justification=c('right','top'))

c_breaks = seq(1:35)

p_clades_conf <- ggplot(clade_df[clade_df$c_type=='conf',],aes(x=n_taxa, color=rel_name, shape=rel_name)) +
  geom_freqpoly(binwidth=1, breaks=c_breaks)

p_clades_conf <- p_clades_conf +
  geom_point(stat='bin',aes(y=after_stat(count)),breaks=c_breaks) +
  theme.leg + s.color_conf + s.shape_conf +
  theme(axis.text.y=element_text(size=12)) +
  geom_vline(xintercept=8, linetype='longdash') + 
  geom_vline(xintercept=13, linetype='dashed') +
  geom_vline(xintercept=35, linetype='dotdash') +
  xlab('proteomes / clade') + ylab('number of clades') +
  s.y.8K + labs(subtitle='A. all confirmed clades')

p_clades_prop <-ggplot(clade_df[clade_df$c_type=='prop',],aes(x=n_taxa, color=rel_name, shape=rel_name)) +
  geom_freqpoly(binwidth=1, breaks=c_breaks)

p_clades_prop = p_clades_prop + 
  geom_point(stat='bin',aes(y=after_stat(count)),breaks=c_breaks) +
  theme.leg + s.color_prop + s.shape_prop +
  geom_vline(xintercept=8, linetype='longdash') + 
  geom_vline(xintercept=13, linetype='dashed') +
  geom_vline(xintercept=35, linetype='dotdash',color='black') +
  xlab('proteomes / clade') + ylab('number of clades') +
  s.y.4Kr + labs(subtitle='B. all proposed clades') +
  theme(axis.text.y=element_text(size=12,hjust=0))

clade_df_M = clade_df[clade_df$MANEstatus=='MANE_good',]

p_clades_conf_M <-ggplot(clade_df_M[clade_df_M$c_type=='conf',],aes(x=n_taxa, color=rel_name, shape=rel_name)) +
  geom_freqpoly(binwidth=1, breaks=c_breaks)

p_clades_conf_M = p_clades_conf_M +  geom_point(stat='bin',aes(y=after_stat(count)),breaks=c_breaks) +
  theme.leg_m + s.color_conf_m + s.shape_conf_m + s.y.4K +
  theme(axis.text.y=element_text(size=12)) +
  geom_vline(xintercept=8, linetype='longdash') + 
  geom_vline(xintercept=13, linetype='dashed') +
  geom_vline(xintercept=35, linetype='dotdash') +
  xlab('proteomes / clade') + ylab('number of clades') + labs(subtitle='C. confirmed MANE good clades')

p_clades_prop_M <-ggplot(clade_df_M[clade_df_M$c_type=='prop',],aes(x=n_taxa, color=rel_name, shape=rel_name)) +
  geom_freqpoly(binwidth=1, breaks=c_breaks)

p_clades_prop_M = p_clades_prop_M +  geom_point(stat='bin',aes(y=after_stat(count)),breaks=c_breaks) +
  theme.leg_m + s.color_prop_m + s.shape_prop_m + s.y.2Kr +
  theme(axis.text.y=element_text(size=12,hjust=0)) +
  geom_vline(xintercept=8, linetype='longdash') + 
  geom_vline(xintercept=13, linetype='dashed') +
  geom_vline(xintercept=35, linetype='dotdash') +
  xlab('proteomes / clade') + ylab('number of clades') + labs(subtitle='D. proposed MANE good clades')

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (! opt$pub) {
  p_all <- ((p_clades_conf | p_clades_prop) / (p_clades_conf_M | p_clades_prop_M)) / doc_panel
  p_all <- p_all  + plot_layout(heights=c(5, 5.0, 0.1))

  ggsave(file=fig_file, plot=p_all, width=8.0, height=7.5)
} else {
  p_all <- ((p_clades_conf | p_clades_prop) / (p_clades_conf_M | p_clades_prop_M)) + plot_layout(heights=c(5.0, 5.0))
  ggsave(file=fig_file, plot=p_all, width=8.0, height=7.0)
}

