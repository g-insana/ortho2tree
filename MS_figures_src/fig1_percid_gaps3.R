#!/usr/bin/env Rscript --vanilla

################
## reads in two files: pair.times with tag<>time<>class
##                     *.gap_summ that reports query/subj/alen/percid/gaps
##                     
## also needs a title for output file names
## 
## creates one plot panels (in one file): of gap count vs percent identity quantilea
##

## 5-Jan-2023 --
##
## modified to make more modular, putting cut/cumm procedures in functions

## 22-Dec-2023 --
## go from cummumlative plot to percentile count plot

## 03-Nov-2023
## uses standard strategy for --pub supression of command line doc
##

library('ggplot2')
library('ggtext')
library(tidyr)
library('getopt')
library('optparse')
library('plotly')
library('Cairo')

library(patchwork)
library('RColorBrewer')

stat_labels = c('g_med'="median gap length",'g_max'='maximum gap length','g_3q'='3rd-quartile gap length')
stat_names = names(stat_labels)
stat_names_str = paste(stat_names,collapse=', ')

## print(stat_names_str)

################
## globals
min_qlen = 100  ## minimum query length

min_stat_nz = 0.1
min_gap_cnt_fn = 0.05
g_gap_offset = 0.05
pct_fact = 100.0
long_gap_thresh = 50

option_list = list(
	    make_option(c("-t","--times"),action='store',help='times file REQUIRED'),
	    make_option(c("-f","--files"),type='character',action='store',help='comma separated *.summ, files REQUIRED'),
	    make_option(c("-Y","--type"),action='store_true',help='has type: canon_HUMAN_V_MOUSE',default=FALSE),
	    make_option(c("-p","--pdf"),type='character',action='store',default='fig1', help='optional pdf file name [default= \"%default\"]'),
	    make_option(c("-s","--stat"),type='character',action='store',default='g_max',help=paste(sprintf('statistic to plot [%s];',stat_names_str),"default=%default")),
	    make_option(c("-P","--pub"),action='store_true',default=FALSE,help='remove plot doc for publication'),
            make_option(c("-T","--txt"),action='store_true',default=FALSE,help='create text files of high gaps'),
            make_option(c("-L","--ptly"),action='store',type='character',default='',help='plotly file name'),
            make_option(c("--Cgapl"),action='store',type='numeric',default=5,help='panel C threshold'),
            make_option(c("--Bthresh"),action='store',type='numeric',default=0.9,help='panel B threshold')
	)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print(opt)

args<- commandArgs(trailingOnly=TRUE)

if (is.null(opt$files) || is.null(opt$times)) {

    ## print_help(opt_parser)
    cat(" percid_up23_gapl3t.R -t times.file -f ,... \n\n")

  my_args = readline(prompt="enter command line arguments: ")
  my_args_l <- strsplit(my_args,' ')
  opt = parse_args(opt_parser,args=my_args_l[[1]]);
  ##quit(save='no',status=1)
  args <- my_args_l[[1]]
}

p.name<-basename(get_Rscript_filename())
plabel=paste(c(p.name,"\n",args,"\n",date()),collapse=' ', sep=' ')

time_file <- opt$time
tab_file_str <- opt$files

comp_name = "fig1"
if (! is.null(opt$pdf)) {
   comp_name = opt$pdf
}

stat_labels = c('g_med'="median gap length",'g_max'='maximum gap length','g_3q'='3rd-quartile gap length')

target_samp_rate = 100

ygap_label = stat_labels[opt$stat]
ygap_label_gcnt = sprintf("percent of alignments with gaps \u2265 %d",opt$Cgapl)
cat(sprintf("%s : %s\n",opt$stat,ygap_label))

################
## names for output plots and text files

top_gaps_file=paste0(comp_name,"_top_gaps.txt")
top_gaps_file2=paste0(comp_name,"_summ_top_gaps.txt")

fig_file=paste0(comp_name,"_id_gaps")
if (opt$pub) {
   fig_file = paste0(fig_file, "_pub")
}
fig_file=paste0(fig_file,'.pdf')

################
## column names for _raw.summ files

save_fields <- c("tag","qseqid","pair","s_pair","qlen","sseqid","slen","percid","alen","n_gaps","q_start","q_end",
                 "n_gaps","g_min","g_1q","g_med","g_3q","g_max","g_max_fn","n_gaps_out","end_diff")

################
## get pair_times

pair_times<-read.table(time_file,header=FALSE,sep='\t',quote='', stringsAsFactors=TRUE)
colnames(pair_times) = c('pair','time','class','class_label')

pair_times$pair<-factor(pair_times$pair,levels=pair_times$pair,ordered=TRUE)
tax.order=pair_times$pair

print("tax.order:")
tax.order

tax.levels=c('rod','mam','vrt','yst','bct','arch')

## pair_times$class<-factor(pair_times$class,levels=tax.levels,ordered=TRUE)

class.order=levels(pair_times$class)
print("class.order")
print(class.order)

## str(pair_times)

################
## read in the data

up_hit1<-NULL

## print(paste0("files: ",tab_file_str))
arg_files<-strsplit(tab_file_str,',',fixed=TRUE)[[1]]

for (tab_file in arg_files) {

  tmp_df<-read.table(file=tab_file,header=TRUE,sep='\t',quote='')

  print(sprintf("file: %s: N: %d",tab_file,length(tmp_df$tag)))

  ## this version just expects human_v_mouse_10090,
  ## no type, or "ref/bad/???"
  ## split tag into type_query_v_subj_staxid

  ## need percid 0.0 - 1.0, not 0.0 - 100.0

  if (max(tmp_df$percid) > 2.0) {
    tmp_df$percid <- tmp_df$percid/100.0
  }   

  tag_list <- strsplit(tmp_df$tag,'_')
  tag_mat <- matrix(unlist(tag_list),byrow=TRUE,ncol=length(tag_list[[1]]))
      
  if (! opt$type) {
    tmp_df$q_name = tag_mat[,1]
    tmp_df$s_name = tag_mat[,3]
    tmp_df$s_taxid = tag_mat[,4]
    tmp_df$a_type = 'canon'
  } else if (tag_mat[1][1] == 'ecoli') {
    tmp_df$q_name = tag_mat[,1]
    tmp_df$s_name = tag_mat[,3]
    tmp_df$s_taxid = tag_mat[,4]
    tmp_df$a_type = 'canon'
  } else {
  ## tag is of the form: 'canon_HUMAN_v_BOVIN' or 'prop_HUMAN_v_BOVIN'
  ## extract fields 3,4
    tmp_df$a_type = tag_mat[,1]
    tmp_df$q_name = tolower(tag_mat[,2])
    tmp_df$s_name = tolower(tag_mat[,4])
  }

  tmp_df$pair = paste0(tmp_df$q_name,'_v_',tmp_df$s_name)


  ## tmp_df$s_pair = ifelse (tmp_df$q_name == 'human', sprintf("%s:%s",tmp_df$s_name,tmp_df$q_name), sprintf("%s:%s",tmp_df$q_name,tmp_df$s_name))

  tmp_df$s_pair = sprintf("%s:%s",tmp_df$q_name,tmp_df$s_name)

  tmp_df$end_diff <- tmp_df$qlen - tmp_df$q_end

  ## print(head(tmp_df))

  ## only longer alignments
  up_hit1 <- rbind(up_hit1,tmp_df[tmp_df$qlen>min_qlen,save_fields])
}

## only target taxa pairs
up_hit1 <- up_hit1[up_hit1$pair %in% tax.order,]
## make zero gaps min_stat_nz
up_hit1$stat_nz = ifelse(up_hit1[,opt$stat] < min_stat_nz, min_stat_nz, up_hit1[,opt$stat])

up_hit1$pair = factor(up_hit1$pair,levels=tax.order, ordered=T)

up_hit1_e90 <- up_hit1[up_hit1$percid > opt$Bthresh,]

################
## all done reading/parsing data
## up_hit1 has gap statistics for all hits, up_hit1_e90 for percid > opt$Bthresh

################
##
## build up_hit1_cnts data structure for >5 gap counts vs percid percentile

gap_cnt_vs_percid <- function(up_hit1, uniq_pairs, gap_stat_thresh, g_gap_offset) {

    ################
    ## for each pair, divide in intervals between 0.50 and 1.0, and count the number of gaps longer than c(5)
    ## 
    ## id quantiles are independent of pair, so simply mark them:

    q_breaks=seq(0.50,1.0,2*g_gap_offset)

    up_hit1 <- within(up_hit1,id_q <- cut(percid, q_breaks, include.lowest=TRUE, labels=FALSE))
    up_hit1$id_q_percid <- q_breaks[up_hit1$id_q]

    ## print(head(up_hit1))

    up_hit1_cnts <- NULL
    ## gap_stat_thresh <- c(5,20,100)
    ## gap_stat_thresh <- c(5,20)

    for (this_pair in uniq_pairs) {

      up_pair_sub <- up_hit1[up_hit1$pair == this_pair,]
      for (gap_stat in gap_stat_thresh) {

        tmp_alns <- aggregate(cbind(a_cnt=tag) ~ id_q_percid,data=up_pair_sub,NROW)
        tmp_gaps <- aggregate(cbind(g_cnt=tag) ~ id_q_percid,data=up_pair_sub[up_pair_sub$g_max >= gap_stat,],NROW)

        tmp_alns$pair = this_pair
        tmp_gaps$pair = this_pair

	zero_intervals_a <- setdiff(q_breaks, tmp_alns$id_q_percid)
	zero_intervals_a <- zero_intervals_a[! is.na(zero_intervals_a)]
    	zero_intervals_a <- zero_intervals_a[order(zero_intervals_a)]

    	## print(paste("length_a:",length(zero_intervals_a)))
    	## print(zero_intervals_a)

    	zero_intervals_g <- setdiff(q_breaks, tmp_gaps$id_q_percid)

	zero_intervals_g <- zero_intervals_g[! is.na(zero_intervals_g)]
    	zero_intervals_g <- zero_intervals_g[order(zero_intervals_g)]
    	tmp_zeros_a <- data.frame(id_q_percid=zero_intervals_a, a_cnt=0,pair=this_pair)
    	tmp_zeros_a <- tmp_zeros_a[! is.na(tmp_zeros_a$id_q_percid),]
    	## print("tmp_zeros_a")
    	## print(tmp_zeros_a)

    	tmp_a_zeros <- rbind(tmp_alns,tmp_zeros_a)
    	tmp_a_zeros <- tmp_a_zeros[order(tmp_a_zeros$id_q_percid),]
    	## print(tmp_a_zeros)

    	tmp_zeros_g <- data.frame(id_q_percid=zero_intervals_g, g_cnt=0,pair=this_pair)
    	tmp_zeros_g <- tmp_zeros_g[! is.na(tmp_zeros_g$id_q_percid),]

    	tmp_gaps$g_cut = gap_stat
    	tmp_zeros_g$g_cut = gap_stat

    	tmp_g_zeros <- rbind(tmp_gaps,tmp_zeros_g)
    	tmp_g_zeros <- tmp_g_zeros[order(tmp_g_zeros$id_q_percid),]

    	tmp_g_zeros$a_cnt <- tmp_a_zeros$a_cnt

    	tmp_g_zeros$g_fract <- ifelse(tmp_g_zeros$a_cnt > 0, 100 * tmp_g_zeros$g_cnt/tmp_g_zeros$a_cnt, 0)

    	## print(tmp_g_zeros)

    	up_hit1_cnts <- rbind(up_hit1_cnts,tmp_g_zeros)
     }
  }

  return(up_hit1_cnts)
}

uniq_pairs = unique(up_hit1$pair)
gap_stat_thresh <- c(opt$Cgapl)

print(aggregate(cbind(NROW=tag) ~ pair, data=up_hit1, FUN=NROW))

up_hit1_cnts <- gap_cnt_vs_percid(up_hit1, uniq_pairs, gap_stat_thresh, g_gap_offset)

## print(head(up_hit1_cnts))

up_hit1_cnts$g_fract_nz = ifelse(up_hit1_cnts$g_fract > 0, up_hit1_cnts$g_fract, min_gap_cnt_fn)
up_hit1_cnts$pair = factor(up_hit1_cnts$pair,levels=tax.order, ordered=T)
up_hit1_cnts$g_cut = factor(up_hit1_cnts$g_cut,levels=gap_stat_thresh, ordered=T)

################
## set up cummulative percid columns

up_percid_cumm = NULL
up_gap90_cumm = NULL
top_gaps_df = NULL

leg.labels_id = c()
leg.labels_gap=c()
leg.labels_gcnt=c()

plot_columns = c('pair','s_pair','class','qseqid','sseqid','percid','g_med','g_max','stat_nz')

for (pair in uniq_pairs) {

    this_pair = pair_times[pair_times$pair == pair,]
    
    ## do statistics for cummulative percid
    pair_subset = up_hit1[up_hit1$pair==pair,]
    pair_subset$class = this_pair$class
    
    pair_cumm = pair_subset[order(pair_subset$percid),plot_columns]
    pair_cumm$row_num = 1:nrow(pair_cumm)
    pair_cumm$fract_row = pair_cumm$row_num/nrow(pair_cumm)
    pair_cumm$c_fract_row = 1.0 - pair_cumm$fract_row
    
    pair_cumm$point_label = sprintf("%s\n%s\n%s\npercid: %.4f stat_nz:%.0f",
    pair_cumm$pair,pair_cumm$qseqid,pair_cumm$sseqid,pair_cumm$percid,pair_cumm$stat_nz)
    
    ## set up label for id with number of hits
    
    # s_pair = gsub('_v_',':',pair)
    s_pair <- pair_subset[1,]$s_pair

    ##  this_label = sprintf("%s\n %.0f My; N = %d",s_pair, this_pair$time,length(pair_cumm$row_num))
    N_form = format(length(pair_cumm$row_num),big.mark=',')
    this_label = sprintf("%s: %.0f My  \n N = %s",s_pair, this_pair$time,N_form)
    this_label = setNames(this_label,pair)
    ##  print(this_label)
    
    leg.labels_id <- append(leg.labels_id,this_label)
    
    ## print(tail(pair_cumm))
    
    samp_rate = trunc(length(pair_cumm$row_num)/target_samp_rate + 0.5)
    
    pair_cumm_samp = pair_cumm[pair_cumm$row_num %% samp_rate ==0,]
    
    print(sprintf("pc: %s samp: %d len: %d",s_pair,samp_rate,NROW(pair_cumm_samp)))
    
    ## statistics for cumm gaps > 90% id
    
    pair_subset_e90 = up_hit1_e90[up_hit1_e90$pair==pair,]
    pair_subset_e90$class = this_pair$class
    
    top_gaps = pair_subset_e90[order(-pair_subset_e90$stat_nz),c('pair','qseqid','qlen','sseqid','slen','percid','n_gaps','g_med','g_max','stat_nz')]

    ## print("head(top_gaps")
    ## print(head(top_gaps))

    if (nrow(top_gaps[top_gaps$stat_nz > long_gap_thresh,]) > 25) {
        top_gaps = top_gaps[top_gaps$stat_nz > long_gap_thresh,]
    } else {
        top_gaps = top_gaps[1:50,]
    }

    top_gaps25 = top_gaps[1:min(nrow(top_gaps),25),]

    top_gaps_df = rbind(top_gaps_df, top_gaps)
    ## print(top_gaps[,c('qseqid','sseqid','percid','g_max')],row.names=FALSE)
    ## write.table(top_gaps,gaps_fd,sep='\t',quote=FALSE,row.names=FALSE)
    
    pair_cumm_gap90 = pair_subset_e90[order(pair_subset_e90$stat_nz),plot_columns]
    
    pair_cumm_gap90$row_num = 1:nrow(pair_cumm_gap90)
    pair_cumm_gap90$fract_row = pair_cumm_gap90$row_num/nrow(pair_cumm_gap90)
    pair_cumm_gap90$c_fract_row = 1.0 - pair_cumm_gap90$fract_row
    
    pair_cumm_gap90$point_label = sprintf("%s\n%s\n%s\npercid: %.4f stat_nz:%.0f",
    pair_cumm_gap90$pair,pair_cumm_gap90$qseqid,pair_cumm_gap90$sseqid,pair_cumm_gap90$percid,pair_cumm_gap90$stat_nz)
    
    ##  this_label = sprintf("%s\n N = %d; %.0f%%",s_pair, length(pair_cumm_gap90$row_num),100.0*length(pair_cumm_gap90$row_num)/length(pair_cumm$row_num))
    N_form = format(length(pair_cumm_gap90$row_num),big.mark=',')
    ## older version
    ## this_label = sprintf("%s\n N = %s; %.0f%%",s_pair, N_form,100.0*length(pair_cumm_gap90$row_num)/length(pair_cumm$row_num))
    
    ## this_label = sprintf("%s (%.1f%%)  \n  %s &ge;20; %s &ge;5",s_pair, 
    ## 100.0*length(pair_cumm_gap90$row_num)/length(pair_cumm$row_num),
    ## format(length(pair_cumm_gap90[pair_cumm_gap90$stat_nz >= 20,]$row_num),big.mark=','),
    ## format(length(pair_cumm_gap90[pair_cumm_gap90$stat_nz >= 5,]$row_num),big.mark=','))
    
    this_label_ex = sprintf("%s (%.1f%%)  \n  %s>100 %s >=20; %s >=5; %s >= 1",s_pair, 
    100.0*length(pair_cumm_gap90$row_num)/length(pair_cumm$row_num),
    format(length(pair_cumm_gap90[pair_cumm_gap90$stat_nz > 100,]$row_num),big.mark=','),
    format(length(pair_cumm_gap90[pair_cumm_gap90$stat_nz >= 20,]$row_num),big.mark=','),
    format(length(pair_cumm_gap90[pair_cumm_gap90$stat_nz >= 5,]$row_num),big.mark=','),
    format(length(pair_cumm_gap90[pair_cumm_gap90$stat_nz >= 1,]$row_num),big.mark=','))
    
    gall_len = length(pair_cumm$row_num)
    g90_len = length(pair_cumm_gap90$row_num)
    this_label = sprintf("%s  \n N= %s (%.0f%%)",s_pair, 
       format(g90_len,big.mark=','),
       100.0*g90_len/gall_len)

    this_label = setNames(this_label,pair)

    print(this_label_ex)

    leg.labels_gap <- append(leg.labels_gap,this_label)
    
    ## print(tail(pair_cumm_gap90))
    
    ge5_len = length(pair_cumm_gap90[pair_cumm_gap90$g_max >= opt$Cgapl,]$row_num)
    if (ge5_len/g90_len < 0.1) {
        gcnt_fmt = "%s  \n N(>%.0f%%)=%s (%.1f%%)"
    } else {
        gcnt_fmt = "%s  \n N(>%.0f%%)=%s (%.0f%%)"
    }

    ## gcnt_fmt = "%s  \n N(>%.0f%%)=%s (%s%%)"
    ## ge5_fn_str = formatC(100*ge5_len/g90_len, digits=2, format='fg')

    gcnt_label = sprintf(gcnt_fmt,s_pair,opt$Bthresh*100,format(ge5_len,big.mark=','),100.0 * ge5_len/g90_len)
    ## gcnt_label = sprintf(gcnt_fmt,s_pair,opt$Bthresh*100,format(ge5_len,big.mark=','),ge5_fn_str)
    

    gcnt_label = setNames(gcnt_label, pair)
    leg.labels_gcnt <- append(leg.labels_gcnt, gcnt_label)

    samp_rate = trunc(length(pair_cumm_gap90$row_num)/target_samp_rate + 0.5)
    
    ## pair_cumm_gap90_samp = pair_cumm_gap90[(pair_cumm_gap90$row_num %% samp_rate==0),]
    pair_cumm_gap90_samp = pair_cumm_gap90[(pair_cumm_gap90$row_num %% samp_rate==0) | (pair_cumm_gap90$qseqid %in% top_gaps25$qseqid),]
    ## cat(sprintf("%s : %0.f\n",pair, length(pair_cumm_gap90_samp$row_num)))
    
    print(sprintf("pc90: %s samp: %d len: %0.f",s_pair,samp_rate,NROW(pair_cumm_gap90_samp)))
    
    up_percid_cumm = rbind(up_percid_cumm,pair_cumm_samp)
    up_gap90_cumm = rbind(up_gap90_cumm,pair_cumm_gap90_samp)
}

## leg.labels_id
## leg.labels_gap

up_percid_cumm$pair = factor(up_percid_cumm$pair,levels=tax.order, ordered=T)
up_percid_cumm$class = factor(up_percid_cumm$class,levels=class.order, ordered=T)

up_gap90_cumm$pair = factor(up_gap90_cumm$pair,levels=tax.order, ordered=T)
up_gap90_cumm$class = factor(up_gap90_cumm$class,levels=class.order, ordered=T)

## remove zeros for log plot

## print("head(up_percid_cumm)")
## head(up_percid_cumm)


################################################################
## set up shapes, colors, and scales

source("set_color_a90.R")
## associate pair names with symbols

pt_names = as.character(pair_times$pair)

a.colors = a.colors[1:length(pt_names)]
a.colors = setNames(a.colors,pt_names)
a.colors['yeast_v_sacar'] = "grey1"
a.colors['ecoli_v_salty'] = "grey1"
print("a.colors")
print(a.colors)

print("pair_times$class")
print(pair_times$class)
print("tax.sym/tax.sym[pair_times$class]")
print(tax.sym)
## print(tax.sym[as.character(pair_times$class)])

print("pt_names")
print(pt_names)

p.sym=setNames(tax.sym[as.character(pair_times$class)],pt_names)
print("p.sym")
print(p.sym)

g.syms = c(0,1,5)

pg.sym = setNames(g.syms, gap_stat_thresh)
## p.sym

a.colors = a.colors[1:length(pair_times$class)]
a.colors = setNames(a.colors,pt_names)

## a.colors

sb.color_id   <- scale_color_manual(name='id_legend', values=a.colors, labels=leg.labels_id)
sb.color_g <- scale_color_manual(name='gap_legend', values=a.colors, labels=leg.labels_gap)
sb.color_gcnt <- scale_color_manual(name='gcnt_legend', values=a.colors, labels=leg.labels_gcnt)

sb.shape_id <- scale_shape_manual(name='id_legend', values=p.sym, labels=leg.labels_id)
sb.shape_g  <- scale_shape_manual(name='gap_legend', values=p.sym, labels=leg.labels_gap)
sb.shape_gcnt <- scale_shape_manual(name='gcnt_legend', values=p.sym, labels=leg.labels_gcnt)

p.alpha = setNames(rep(1.0,length(pt_names)),pt_names)
p.alpha['ecoli_v_salty'] = 0.33
p.alpha['yeast_v_sacar'] = 0.33

sb.alpha_id <- scale_alpha_manual(name='id_legend', values=p.alpha, labels=leg.labels_id)
sb.alpha_g <- scale_alpha_manual(name='gap_legend', values=p.alpha, labels=leg.labels_gap)
## sb.alpha_gcnt <- scale_alpha_manual(name='gcnt_legend', values=p.alpha, labels=leg.labels_gcnt)

sym.size = 1
p.size = setNames(rep(sym.size, length(pt_names)),pt_names)
p.size['ecoli_v_salty'] = 1.33
p.size['yeast_v_sacar'] = 1.33
sb.size_id <- scale_size_manual(name='id_legend', values=p.size, labels=leg.labels_id)
sb.size_g <- scale_size_manual(name='gap_legend', values=p.size, labels=leg.labels_gap)

scalex.id <- scale_x_continuous(breaks=seq(0.0,1.0,0.2),limits=c(0.0,1.0))
scalex.gap <- scale_x_continuous(breaks=seq(0.0,1.0,0.2),limits=c(0.0,1.02))
scalex.gap_cnt <- scale_x_continuous(breaks=pct_fact*seq(0.5,1.0,0.1),limits=pct_fact*c(0.5,1.0))

## y-axis for identity
scaley.id <- scale_y_continuous("percent identity", breaks=pct_fact*seq(0.4, 1.0, 0.1), limits=pct_fact*c(0.4,1.0))

scaley.gap_r <- scale_y_continuous("number of gaps", position='right')

l.yvals <- c(0.1,1.0,5.0, 10.0,20.0,50.0,100.0,200.0,500.0,1000.0,2000.0)
l.yvals100 <- c(0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0)
ly.labels <- c("0","1","5","10","20","50","100","200","500","1000","2000")
ly.labels100 <- c("0","0.1","0.2","0.5","1","2","5","10","20","50","100")
## scaley.log_g <- scale_y_log10(ygap_label, breaks=l.yvals,minor_breaks=l.yvals,labels=ly.labels,limits=c(0.1,1000.0))
scaley.log_g <- scale_y_log10(ygap_label, breaks=l.yvals,minor_breaks=l.yvals,labels=ly.labels,limits=c(0.1,1000.0))
scaley.log_gr <- scale_y_log10(ygap_label, breaks=l.yvals,minor_breaks=l.yvals,labels=ly.labels,position='right',limits=c(0.1,2000.0))
scaley.log_gcnt <- scale_y_log10(ygap_label_gcnt, breaks=l.yvals100,minor_breaks=l.yvals100,limits=c(0.05,100.0))
scaley.log_gcntr <- scale_y_log10(ygap_label_gcnt, breaks=l.yvals100,minor_breaks=l.yvals100,labels=ly.labels100,position='right',limits=c(0.05,100.0))


################################################################
## themes

theme_set(theme_linedraw(base_size=14))
theme.base <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0),
  plot.caption=element_text(size=6,hjust=0),
  axis.text.x = element_text(size=12),
  axis.text.y = element_text(size=12),
  axis.title = element_markdown(size=12),
  ## legend.text=element_text(size=9),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.text.align = 0,
  legend.title=element_blank(),
  legend.text=element_markdown(size=9),
  legend.margin=margin(c(0.2,5,5))
)

theme.leg_id = theme.base + theme(legend.position=c(0.98,0.02),legend.justification=c('right','bottom'))
theme.leg_gap = theme.base + theme(legend.position=c(0.02,0.98))
theme.leg_cnt = theme.base + theme(legend.position=c(0.02,0.02),legend.justification=c('left','bottom'))

size.guide = guides(color = guide_legend(override.aes = list(size = 2)))

################################################################
## plots

p_id <-ggplot(up_percid_cumm,aes(x=fract_row, y=pct_fact * percid)) +
theme.leg_id + scaley.id + scalex.id + ## scale_x_reverse() + 
  geom_hline(yintercept=100*opt$Bthresh,linetype='longdash')+
  geom_point(aes(color=pair,shape=pair,alpha=pair),size=sym.size) + 
  sb.color_id + sb.shape_id + sb.alpha_id + size.guide +
  labs(x='fraction of alignments',
      subtitle='A.')

p_gap90r <-ggplot(up_gap90_cumm,aes(x=fract_row, y=stat_nz)) +
  theme.leg_gap + scalex.gap + scaley.log_gr + 
  geom_hline(yintercept=opt$Cgapl,linetype='dotdash')+
  annotate('rect',xmin=0.98,xmax=1.019999,ymin=50,ymax=2000,alpha=0.0,color='red') +
  geom_point(aes(color=pair,shape=pair,alpha=pair,size=pair)) + 
  sb.color_g + sb.shape_g + sb.alpha_g + sb.size_g + size.guide +
  labs(x=sprintf('f\'n of alignments > %.0f%% identical',100.0*opt$Bthresh),
      subtitle='B.')

print(up_hit1_cnts)

p_gapcnt_r <- ggplot(up_hit1_cnts,aes(x=pct_fact*(id_q_percid+g_gap_offset), y=g_fract_nz, color=pair, shape=pair)) +
  sb.color_gcnt + scalex.gap_cnt + sb.shape_gcnt + theme.leg_cnt + scaley.log_gcntr +
  labs(x='percent identity', 
      subtitle='C.') + 
  geom_vline(xintercept=100.0*opt$Bthresh,linetype='longdash')+
  geom_point() + ##  geom_step(position=position_nudge(x= -g_gap_offset))
  geom_line()

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

p_id_gap <- p_id + p_gap90r + p_gapcnt_r

if (! opt$pub) {
  p_id_gap <- p_id_gap / doc_panel + plot_layout(heights=c(10, 0.1))
  ggsave(file=fig_file, plot=p_id_gap, width=12.5, height=4.6, device=cairo_pdf)
} else {
  ggsave(file=fig_file, plot=p_id_gap, width=12.5, height=4.5, device=cairo_pdf)
}

if (opt$txt) {
   gaps_fd = file(top_gaps_file,open='w')

   write(paste0("# ",paste0(colnames(top_gaps_df),collapse='\t',sep='')),gaps_fd)
   for (pair in tax.order) {
       these_rows = top_gaps_df[top_gaps_df$pair==pair,]
       write(sprintf(">%s %.0f",pair,NROW(these_rows)),gaps_fd)
       write.table(these_rows,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,gaps_fd)
   }      
   close(gaps_fd)
}
