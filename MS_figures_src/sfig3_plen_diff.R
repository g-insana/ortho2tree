#!/usr/bin/env Rscript --vanilla

################
##
## 24-Nov-2023
##
## plen_diff1r.R -f output_gc231004
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
library('ggallin')   ## get log axis scaler for negative numbers
library('reshape2')

## library(cowplot)
library(patchwork)
library('RColorBrewer')
library('Cairo')

option_list = list(
	    make_option(c("-f","--file"),type='character',action='store',help='output_gc file REQUIRED'),
	    make_option(c("-p","--pdf"),type='character',action='store',default='fig4', help='optional pdf file name [default= \"%default\"]'),
	    make_option(c("-P","--pub"),action='store_true',default=FALSE, help='remove plot documentation')
	    )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$file)) {
    print_help(opt_parser)
    cat(" plen_diff.R  -f output_gc231004 \n\n")
    quit(save='no',status=1)
}

args<- commandArgs(trailingOnly=TRUE)

p.name<-basename(get_Rscript_filename())
plabel=paste(c(p.name,args,"\n",date()),collapse=' ', sep=' ')

plot_name = "fig4"
if (! is.null(opt$pdf)) {
   plot_name = opt$pdf
}

################
## names for output plots and text files

fig_file=paste0(plot_name,"_ldelta")
if (opt$pub) {
   fig_file = paste0(fig_file,"_pub")
}
fig_file = paste0(fig_file,".pdf")

fig_file

################

read_fields <- c("pthr_id","taxon","canon_acc","canon_len","canon_md5","prop_acc","prop_len","prop_md5","rank_score","canon_cost","prop_cost","n_sp","n_tax","n_canon","wn_canon","scaled_prop_f","scaled_p_cost","p_len_diff","clade","clade_members","MANEstatus","release","dataset")

save_fields <- c("pthr_id","taxon","canon_acc","canon_len","prop_acc","prop_len","canon_cost","prop_cost","n_sp","n_tax","clade_members","MANEstatus","release","dataset")

target_taxa = toupper(c('human','gorgo','mouse','rat','bovin'))
taxa_order = setNames(1:5,target_taxa)

################
## read in the data

prop_df<-read.table(file=opt$file,header=FALSE,sep='\t',quote='',col.names=read_fields)
prop_df <- prop_df[prop_df$taxon %in% target_taxa,save_fields]

prop_df$l_delta <- prop_df$prop_len - prop_df$canon_len
prop_df$delta_pos <- ifelse(prop_df$l_delta > 0, TRUE, FALSE)

prop_df_tr = prop_df[grepl('^tr\\|',prop_df$canon_acc,perl=TRUE),]
prop_df_tr$delta_pos <- ifelse(prop_df_tr$l_delta > 0, TRUE, FALSE)

print(sprintf("length all: %d; tr/tr: %d",NROW(prop_df),NROW(prop_df_tr)))

del_breaks = c(-5000,-2000,-1000,-500,-200,-100,-50,-20,-10,-5,-2,-1,1,2,5,10,20,50,100,200,500,1000,2000,5000)
pdel_breaks = c(-1000,-100,-10,-1,1,10,100,1000)
##pdel_break_lab = expression(-10^3,-10^2,-10^1,10^1,10^2,10^3)
pdel_break_lab = c('-1000','-100','-10','-1','1','10','100','1000')

## set the colors
a.colors <- brewer.pal(7,'Dark2') # 'Dark2', 'Set2', 'Paired'

cl = length(target_taxa)+1
a.colors[cl] <- a.colors[7]
a.colors <- a.colors[1:cl]

a.colors <- setNames(a.colors,c('GORGO','NA1','MOUSE','RAT','BOVIN','HUMAN'))

theme_set(theme_linedraw(base_size=14))
theme.base <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0),
  plot.caption=element_text(size=6,hjust=0),
  axis.text.x = element_text(size=9),
  axis.text.y = element_text(size=12),
  legend.text=element_text(size=9,family='Menlo'),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.text.align = 0,
  legend.title=element_blank(),
  legend.margin=margin(c(0.2,5,5))
)

  theme.leg = theme.base + theme(legend.position=c(0.98,0.98), legend.justification=c('right','top'))

## pdf(file=fig1i_file,width=8, height=4)

## build count structure by taxa

prop_df$taxon = factor(prop_df$taxon, levels=target_taxa, ordered=TRUE)
prop_df_tr$taxon = factor(prop_df_tr$taxon, levels=target_taxa, ordered=TRUE)

ldelta_sum = aggregate(cbind(N=l_delta) ~ taxon, data=prop_df ,FUN=NROW)
ldelta_sum_pos = aggregate(cbind(Npos=l_delta) ~ taxon, data=prop_df[prop_df$delta_pos,],FUN=NROW)

ldelta_sum <- merge(ldelta_sum, ldelta_sum_pos, by='taxon')
ldelta_sum$Nneg = ldelta_sum$N - ldelta_sum$Npos

##taxa_labels = sprintf("%s  \n  N=%d",tolower(ldelta_sum$taxon),ldelta_sum$N)
taxa_labels = sprintf("%-5s \u0394>0:%4d\n      \u0394<0:%4d",tolower(ldelta_sum$taxon),ldelta_sum$Npos, ldelta_sum$Nneg)
##taxa_labels = sprintf("%-5s N+:%4d\n      N\u2013:%4d",tolower(ldelta_sum$taxon),ldelta_sum$Npos, ldelta_sum$Nneg)


taxa_labels = setNames(taxa_labels, ldelta_sum$taxon)

taxa_labels = taxa_labels[order(taxa_order[names(taxa_labels)])]

print("all sums")
ldelta_sum
taxa_labels

ldelta_sum_tr = aggregate(cbind(N=l_delta) ~ taxon, data=prop_df_tr ,FUN=NROW)
ldelta_sum_tr_pos = aggregate(cbind(Npos=l_delta) ~ taxon, data=prop_df_tr[prop_df_tr$delta_pos,],FUN=NROW)

print("ldelta_sum_tr")
ldelta_sum_tr
ldelta_sum_tr_pos

ldelta_sum_tr <- merge(ldelta_sum_tr, ldelta_sum_tr_pos, by='taxon', all.x=TRUE, all.y=TRUE)
ldelta_sum_tr$Npos = ifelse(is.na(ldelta_sum_tr$Npos),0,ldelta_sum_tr$Npos)

ldelta_sum_tr$Nneg = ldelta_sum_tr$N - ldelta_sum_tr$Npos

ldelta_sum_tr

taxa_labels_tr = sprintf("%-5s \u0394>0:%4d\n      \u0394<0:%4d",tolower(ldelta_sum_tr$taxon),ldelta_sum_tr$Npos, ldelta_sum_tr$Nneg)
## taxa_labels_tr = sprintf("%-5s N+:%4d\n      N\u2013:%4d",tolower(ldelta_sum_tr$taxon),ldelta_sum_tr$Npos, ldelta_sum_tr$Nneg)

taxa_labels_tr = setNames(taxa_labels_tr, ldelta_sum_tr$taxa)


prop_df_tr_gt0 <- prop_df_tr[prop_df_tr$l_delta > 0,]
ldelta_sum_tr0 <- do.call(data.frame,aggregate(prop_df_tr_gt0$l_delta, by=list(prop_df_tr_gt0$taxon),FUN=NROW))
colnames(ldelta_sum_tr0)=c("taxon","N")

ldelta_merge <- merge(ldelta_sum_tr, ldelta_sum_tr0, by='taxon')
ldelta_merge$diff = ldelta_merge$N.x - ldelta_merge$N.y
print("total(N.x) vs >0(N.y) counts")
ldelta_merge

colSums(ldelta_merge[,2:4])

## this must be done after aggregated counts are available
##
s.color <- scale_color_manual(values=a.colors, labels=taxa_labels)
s.color_tr <- scale_color_manual(values=a.colors, labels=taxa_labels_tr)

################
## do the plots

xaxis_title="(prop. - canon.) length"
yaxis_title="number of changes"

################
## all changes
p_delta <-ggplot(prop_df,aes(x=l_delta,color=taxon)) +
  geom_freqpoly(breaks=del_breaks)

h_delta <- layer_data(p_delta)

p_delta <- p_delta +
  geom_point(stat='bin',aes(y=after_stat(count)),breaks=del_breaks,shape=1) +
  theme.leg +   
  s.color +
  scale_x_continuous(name=xaxis_title,trans=pseudolog10_trans, breaks=pdel_breaks, labels=pdel_break_lab) +
##  scale_x_continuous(name=xaxis_title, breaks=pdel_breaks, labels=pdel_break_lab) +
  scale_y_continuous(name=yaxis_title) +
  labs(subtitle='A. all changes')

h_delta <- h_delta[,c('group','count','x','xmin')]
h_delta$taxon = target_taxa[h_delta$group]
h_delta = h_delta[,c('taxon','x','xmin','count')]
w_delta <- dcast(h_delta,x + xmin ~ taxon, value.var='count')

print("all prop counts")
w_delta

## print(h_delta)

################
## canon_tr only
p_delta_tr <-ggplot(prop_df_tr,aes(x=l_delta, color=taxon)) +
  geom_freqpoly(breaks=del_breaks)

h_delta_tr <- layer_data(p_delta_tr)

p_delta_tr <- p_delta_tr +
  geom_point(stat='bin',aes(y=after_stat(count)),breaks=del_breaks,shape=1) +
  theme.leg + 
  s.color_tr +
  scale_x_continuous(name=xaxis_title, trans=pseudolog10_trans, breaks=pdel_breaks, labels=pdel_break_lab) +
##  scale_x_continuous(name=xaxis_title, breaks=pdel_breaks, labels=pdel_break_lab) +
  scale_y_continuous(name=yaxis_title, position='right') +
  labs(subtitle='B. Trembl/Trembl only')

h_delta_tr <- h_delta_tr[,c('group','count','x','xmin')]
h_delta_tr$taxon = target_taxa[h_delta_tr$group]
h_delta_tr <- h_delta_tr[,c('taxon','x','xmin','count')]

w_delta_tr <- dcast(h_delta_tr,x + xmin ~ taxon, value.var='count')
print("tr:tr prop counts")
w_delta_tr



doc_panel = ggplot() + theme_void() + labs(caption=plabel)

p_2deltas <- p_delta + p_delta_tr

if (! opt$pub) {
  p_2deltas <- p_2deltas / doc_panel + plot_layout(heights=c(5, 0.1))
  ggsave(file=fig_file, plot=p_2deltas, width=7.5, height=5.5, device=cairo_pdf)
} else {
  ggsave(file=fig_file, plot=p_2deltas, width=7.5, height=5.0, device=cairo_pdf)
}

