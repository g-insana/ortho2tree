#!/usr/bin/env Rscript --vanilla

################
##
## 24-Nov-2023
##
## sfig1_cost_hist.R -f output_gc231004
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
library('reshape2')

## library(cowplot)
library(patchwork)
library('RColorBrewer')
library('Cairo')

option_list = list(
	    make_option(c("-f","--files"),type='character',action='store',help='output_changes,output_confirm file list REQUIRED'),
	    make_option(c("-p","--pdf"),type='character',action='store',default='supp_fig1', help='optional pdf file name [default= \"%default\"]'),
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

plot_name = "suppl_fig1"
if (! is.null(opt$pdf)) {
   plot_name = opt$pdf
}

################
## globals
min_cost_nz = 0.001

################
## names for output plots and text files

fig_file=paste0(plot_name,"_cost")

if (opt$pub) {
   fig_file = paste0(fig_file,"_pub")
}
fig_file = paste0(fig_file,".pdf")

fig_file

################

chng_fields <- c("pthr_id","taxon","canon_acc","canon_len","prop_acc","prop_len","rank_score","canon_cost","prop_cost","n_sp","n_tax","n_canon","wn_canon","scaled_prop_f","scaled_p_cost","p_len_diff","clade","clade_members","MANEstatus")

conf_fields <- c("pthr_id","canon_cost","n_sp","n_tax","clade","clade_members","MANEstatus")

save_fields_conf <- c("pthr_id","canon_cost","MANEstatus")
save_fields_chng <- c("pthr_id","prop_cost","MANEstatus")

################
## read in the data

cost_df = NULL
for (tab_file in strsplit(opt$files,',',fixed=TRUE)[[1]]) {

    print(paste("tab_file:",tab_file))
    tmp_df<-read.table(file=tab_file,header=FALSE,sep='\t',quote='')
    if (NCOL(tmp_df) == length(chng_fields)) {
         colnames(tmp_df) = chng_fields
         tmp_df <- tmp_df[,save_fields_chng]
         tmp_df$type = 'changed'
    } else {
         colnames(tmp_df) = conf_fields

         tmp_df <- tmp_df[,save_fields_conf]
         colnames(tmp_df) = save_fields_chng

         tmp_df$type='confirmed'
    }
    cost_df = rbind(cost_df, tmp_df)
}

cost_df$type = factor(cost_df$type, levels=c('confirmed','changed'), ordered=TRUE)
cost_df$prop_cost_nz = ifelse(cost_df$prop_cost < min_cost_nz,min_cost_nz,cost_df$prop_cost)

q10tile <-function(x) {
  quantile(x,probs=seq(0.0,1.0,0.05))
}

aggregate(prop_cost ~ type, data=cost_df, FUN=q10tile)

cost_rows <- aggregate(cbind(nrow=prop_cost) ~ type, data=cost_df, FUN=NROW)
cost_rows

theme_set(theme_linedraw(base_size=14))
theme.base <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0),
  plot.caption=element_text(size=6,hjust=0),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  legend.text=element_text(size=9),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.text.align = 0,
  legend.title=element_blank(),
  legend.margin=margin(c(0.2,5,5))
)

theme.leg = theme.base + theme(legend.position=c(0.98,0.02), legend.justification=c('right','bottom'))

source("set_color_a90.R")
## associate pair names with symbols

a.colors=setNames(a.colors[1:2],c('confirmed','changed'))
sb.color = scale_color_manual(values=a.colors,labels=sprintf("%s\n  N=%s",cost_rows$type,format(cost_rows$nrow,big.mark=',')))

sb_x.log10 = c(0.001,0.002,0.005,0.010,0.020)
label_x.log10 = c('0.0','0.002','0.005','0.010','0.020')
scale.x_log10 = scale_x_log10('clade cost',breaks=sb_x.log10,labels=label_x.log10, limits=c(0.001,0.025))
scale.x_lin = scale_x_continuous('clade cost',limits=c(0.000,0.025))

## scale.y = scale_y_continuous('fraction of clades',limits=c(0.0,1.0),breaks=seq(0,1.0,0.2),labels=sprintf("%.1f",seq(0,1.0,0.2))) +
s.y_breaks4 = seq(0.4,1.0,0.1)
s.y_breaks10 = seq(0.0,1.0,0.2)
scale.y4 = scale_y_continuous('fraction of clades',limits=c(0.4,1.0),breaks=s.y_breaks4,labels=sprintf("%.1f",s.y_breaks4))
scale.y10 = scale_y_continuous('fraction of clades',limits=c(0.0,1.0),breaks=s.y_breaks10,labels=sprintf("%.1f",s.y_breaks10))

################
## all changes
p_cost <-ggplot(cost_df,aes(x=prop_cost, color=type)) +
  scale.y10 + scale.x_lin + theme.leg + sb.color +
##  geom_freqpoly(position='identity')
  stat_ecdf(geom='step',pad=FALSE)
##   stat_bin(aes(y= cumsum(after_stat(count))),geom='step')
  

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (! opt$pub) {
  f_plot <- p_cost / doc_panel + plot_layout(heights=c(5, 0.1))
  ggsave(file=fig_file, plot=f_plot, width=4.5, height=4.5, device=cairo_pdf)
} else {
  ggsave(file=fig_file, plot=p_cost, width=4.5, height=4.0, device=cairo_pdf)
}

