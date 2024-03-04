#!/usr/bin/env Rscript --vanilla

################
## reads in two files: pair.times with tag<>time<>class
##                     *.gap_summ that reports query/subj/alen/percid/gaps
##                     
## also needs a title for output file names
## 
## creates 3 plot files: identity qq-plot, gaps qq-plot (at 90% identity), combined id/gaps qq-plot
##

## 5-Jan-2023
## modified to plot gap > 5 cnt vs percid percentile

## 20-Nov-2023
## modified to use cummulative plots of sampled percent identity raw data
##
## 03-Nov-2023
## uses standard strategy for --pub supression of command line doc

## one of the challenges of this plot is to get the colors/symbols
## right, while also getting the labels right.  It is easy to get the colors right -- just match them to the pair (human_v_gorgo, etc)
## but if you do that, you do not get the labels in the legends correct.  To get those correct, you need to match the colors to the pair_type (human_v_gorgo_canon)
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

option_list = list(
	    make_option(c("-t","--times"),action='store',help='times file REQUIRED'),
	    make_option(c("-f","--files"),type='character',action='store',help='comma separated *.summ, files REQUIRED'),
	    make_option(c("-Y","--type"),action='store_true',help='has type: canon_HUMAN_V_MOUSE',default=FALSE),
	    make_option(c("-p","--pdf"),type='character',action='store',default='fig5', help='optional pdf file name [default= \"%default\"]'),
      make_option(c("-s","--stat"),type='character',action='store',default='g_max',help=paste(sprintf('statistic to plot [%s];',stat_names_str),"default=%default")),
      make_option(c("-P","--pub"),action='store_true',default=FALSE,help='remove plot doc for publication'),
      make_option(c("-T","--txt"),action='store_true',default=FALSE,help='create text files of high gaps'),
      make_option(c("-L","--ptly"),action='store',type='character',default='',help='plotly file name'),
      make_option(c("-S","--supp"),action='store_true',default=FALSE,help='make separate supplemental file')
	    )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

## globals
target_samp_rate = 200
min_qlen = 100  ## minimum query length
min_stat_nz = 0.1
min_gap_cnt_fn = 0.05

long_gap_thresh = 50

p.name<-basename(get_Rscript_filename())

if (is.null(opt$files) || is.null(opt$times)) {

    print_help(opt_parser)

    cat(sprintf(" %s -t times.file -f ,... \n\n",p.name))

    quit(save='no',status=1)
}

args<- commandArgs(trailingOnly=TRUE)

plabel=paste(c(p.name,"\n",args,"\n",date()),collapse=' ', sep=' ')
## plabel
time_file <- opt$time
tab_file_str <- opt$files

stat_labels = c('g_med'="median gap length",'g_max'='maximum gap length','g_3q'='3rd-quartile gap length')

ygap_label = stat_labels[opt$stat]
## cat(sprintf("%s : %s\n",opt$stat,ygap_label))

ygap_label_gcnt = "percent of alignments with gaps \u22655"

################
## names for output plots and text files

plot_name = 'fig5'
if (! is.null(opt$pdf)) {
  plot_name = opt$pdf
}

fig_file_b=paste0(plot_name,"_id_gapsF3")

top_gaps_file=paste0(plot_name,"_top_gaps.txt")

if (opt$pub) {
   fig_file = paste0(fig_file_b, "_pub.pdf")
} else {
   fig_file = paste0(fig_file_b,'.pdf')
}

if (opt$supp) {
  if (opt$pub) {
     fig_file = paste0(fig_file_b,'_main_pub.pdf')
     fig_file_supp = paste0(fig_file_b,'_supp_pub.pdf')
  } else {
     fig_file = paste0(fig_file_b,'_main.pdf')
     fig_file_supp = paste0(fig_file_b,'_supp.pdf')
  }  
}

fig_file
## fig_file_supp

################
## column names for _raw.summ files

save_fields <- c("tag","pthr_id","qseqid","pair","s_pair","a_type","pair_type","qlen","sseqid","slen","percid","alen","n_gaps","q_start","q_end",
                 "n_gaps","g_min","g_1q","g_med","g_3q","g_max","g_max_fn","n_gaps_out","end_diff")

################
## get pair_times

pair_times<-read.table(time_file,header=FALSE,sep='\t',quote='', stringsAsFactors=TRUE)
colnames(pair_times) = c('pair','time','class','class_label')

##pair_times = pair_times[order(pair_times$time),]
## print("pair_times")
## pair_times

pair_times$pair<-factor(pair_times$pair,levels=pair_times$pair,ordered=TRUE)
tax.order=pair_times$pair

## print("tax.order:")
## tax.order

tax.levels=c('rod','mam','vrt','yst','bct','arch')

pair_times$class<-factor(pair_times$class,levels=tax.levels,ordered=TRUE)

class.order=levels(pair_times$class)
## print("class.order")
## print(class.order)

## str(pair_times)

################
## read in the data

up_hit1<-NULL

## print(paste0("files: ",tab_file_str))
arg_files<-strsplit(tab_file_str,',',fixed=TRUE)[[1]]

## read in the data, generate $a_type when needed

for (tab_file in arg_files) {

  print(paste("file:",tab_file))

  tmp_df<-read.table(file=tab_file,header=TRUE,sep='\t',quote='')

  ## this version just expects human_v_mouse_10090,
  ## no type, or "ref/bad/???"
  ## split tag into type_query_v_subj_staxid

  ## need percid 0.0 - 1.0, not 0.0 - 100.0

  if (max(tmp_df$percid) > 2.0) {
    tmp_df$percid <- tmp_df$percid/100.0
  }   

  tag_list <- strsplit(tmp_df$tag,'_')
  tag_mat <- matrix(unlist(tag_list),byrow=TRUE,ncol=length(tag_list[[1]]))
      
  ## print("tag_mat")
  ## print(tag_mat)

  if (! opt$type) {
    tmp_df$q_name = tag_mat[,1]
    tmp_df$s_name = tag_mat[,3]
    tmp_df$s_taxid = tag_mat[,4]
    tmp_df$a_type = 'canon'
  } else if (tag_mat[1][1] == 'ecoli' | tag_mat[1][1] == 'yeast') {
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
  tmp_df$pair_type = sprintf("%s_%s",tmp_df$pair,tmp_df$a_type)

  tmp_df$s_pair = sprintf("%s:%s",tmp_df$q_name,tmp_df$s_name)

  tmp_df$end_diff <- tmp_df$qlen - tmp_df$q_end

  ## print(head(tmp_df))

  ## only longer alignments
  up_hit1 <- rbind(up_hit1,tmp_df[tmp_df$qlen>min_qlen,save_fields])
}

## make two duplicate sets of ecoli_v_salty data if it is there:
if ('ecoli_v_salty' %in% unique(up_hit1$pair)){
  tmp_df1 = up_hit1[up_hit1$pair=='ecoli_v_salty',]
  tmp_df2 = tmp_df1

  tmp_df1$a_type = 'prop'
  tmp_df1$pair_type = sprintf("%s_%s",tmp_df1$pair,'prop')
  up_hit1 <- rbind(up_hit1,tmp_df1[tmp_df1$qlen>min_qlen,save_fields])

  tmp_df2$a_type = 'conf'
  tmp_df2$pair_type = sprintf("%s_%s",tmp_df2$pair,'conf')
  up_hit1 <- rbind(up_hit1,tmp_df2[tmp_df2$qlen>min_qlen,save_fields])
}

if ('yeast_v_sacar' %in% unique(up_hit1$pair)){
  tmp_df1 = up_hit1[up_hit1$pair=='yeast_v_sacar',]
  tmp_df2 = tmp_df1

  tmp_df1$a_type = 'prop'
  tmp_df1$pair_type = sprintf("%s_%s",tmp_df1$pair,'prop')
  up_hit1 <- rbind(up_hit1,tmp_df1[tmp_df1$qlen>min_qlen,save_fields])

  tmp_df2$a_type = 'conf'
  tmp_df2$pair_type = sprintf("%s_%s",tmp_df2$pair,'conf')
  up_hit1 <- rbind(up_hit1,tmp_df2[tmp_df2$qlen>min_qlen,save_fields])
}

up_hit1$un_aln=up_hit1$q_start+up_hit1$end_diff
up_hit1$nn_gaps=up_hit1$n_gaps/(up_hit1$n_gaps+up_hit1$alen)

## subset based on tax_pair
up_hit1_e <- up_hit1[up_hit1$pair %in% tax.order,]
up_hit1_e$pair = factor(up_hit1_e$pair,levels=tax.order, ordered=T)

####
## make sure no zero gaps for opt$stat that will be plotted
##
up_hit1_e$stat_nz = ifelse(up_hit1_e[,opt$stat] > min_stat_nz, up_hit1_e[,opt$stat], min_stat_nz)
## print("head(up_hit1_e)")
## print(head(up_hit1_e))

up_hit1_e$pair = factor(up_hit1_e$pair,levels=tax.order, ordered=T)

up_hit1_e90 <- up_hit1_e[up_hit1_e$percid > 0.9,]
## up_hit1_e90 <- up_hit1_e

################
## do the gap vs percid percentile, but do it separately for canonical, proposed

################
##
## build up_hit1_cnts data structure for >5 gap counts vs percid percentile

## g_plot_offset=0.050
g_plot_offset=0.05

gap_cnt_vs_percid <- function(up_hit1, uniq_pairs, gap_stat_thresh, fract_bin) {

    ################
    ## for each pair, divide in intervals between 0.50 and 1.0, and count the number of gaps longer than c(5)
    ## 
    ## id quantiles are independent of pair, so simply mark them:

    q_breaks=seq(0.50,1.0,fract_bin)

    up_hit1 <- within(up_hit1,id_q <- cut(percid, q_breaks, include.lowest=FALSE, labels=FALSE))
    up_hit1$id_q_percid <- q_breaks[up_hit1$id_q]

    ## print("head(up_hit1)")
    ## print(head(up_hit1))

    up_hit1_cnts <- NULL

    for (this_pair in uniq_pairs) {

        up_pair_sub <- up_hit1[up_hit1$pair == this_pair,]

## 	print(paste("This pair:",this_pair))
##	print(head(up_pair_sub))

        for (gap_stat in gap_stat_thresh) {
    
           tmp_alns <- aggregate(cbind(a_cnt=tag) ~ id_q_percid,data=up_pair_sub,NROW)
       
           tmp_alns$pair = this_pair

           zero_intervals_a <- setdiff(q_breaks, tmp_alns$id_q_percid)
           zero_intervals_a <- zero_intervals_a[! is.na(zero_intervals_a)]
           zero_intervals_a <- zero_intervals_a[order(zero_intervals_a)]
       
           ## print(paste("alns length_a:",length(zero_intervals_a)))
           ## print(zero_intervals_a)
       
	   if (length(up_pair_sub[up_pair_sub$g_max >= gap_stat,]) > 0) {
	       tmp_gaps <- aggregate(cbind(g_cnt=tag) ~ id_q_percid,data=up_pair_sub[up_pair_sub$g_max >= gap_stat,],NROW)
               tmp_gaps$pair = this_pair

               zero_intervals_g <- setdiff(q_breaks, tmp_gaps$id_q_percid)
               zero_intervals_g <- zero_intervals_g[! is.na(zero_intervals_g)]
               zero_intervals_g <- zero_intervals_g[order(zero_intervals_g)]

	       tmp_gaps$g_cut = gap_stat

	   } else {
	       ## print(sprintf("%s no > %d",this_pair, gap_stat))

	       tmp_gaps <- NULL
	       zero_intervals_g <- q_breaks
	   }
	       
           tmp_zeros_a <- data.frame(id_q_percid=zero_intervals_a, a_cnt=0,pair=this_pair)
           tmp_zeros_a <- tmp_zeros_a[! is.na(tmp_zeros_a$id_q_percid),]
           ## print("tmp_zeros_a")
           ## print(tmp_zeros_a)
       
           tmp_a_zeros <- rbind(tmp_alns,tmp_zeros_a)
           tmp_a_zeros <- tmp_a_zeros[order(tmp_a_zeros$id_q_percid),]
           ## print(tmp_a_zeros)
       
           tmp_zeros_g <- data.frame(id_q_percid=zero_intervals_g, g_cnt=0,pair=this_pair)
           tmp_zeros_g <- tmp_zeros_g[! is.na(tmp_zeros_g$id_q_percid),]
       
           tmp_zeros_g$g_cut = gap_stat
       
           tmp_g_zeros <- rbind(tmp_gaps,tmp_zeros_g)
           tmp_g_zeros <- tmp_g_zeros[order(tmp_g_zeros$id_q_percid),]

           tmp_g_zeros$a_cnt <- tmp_a_zeros$a_cnt

           ## print("tmp_g_zeros")
           ## print(tmp_g_zeros)

           tmp_g_zeros$g_fract <- ifelse(tmp_g_zeros$a_cnt > 0, 100 * tmp_g_zeros$g_cnt/tmp_g_zeros$a_cnt, 0)
       
           ## print(tmp_g_zeros)
       
           up_hit1_cnts <- rbind(up_hit1_cnts,tmp_g_zeros)
        }
    }

    return(up_hit1_cnts)
}

uniq_pairs = unique(up_hit1_e$pair)
## print(uniq_pairs)
## print(is.na(uniq_pairs))
uniq_pairs = uniq_pairs[! is.na(uniq_pairs)]
gap_stat_thresh = c(5)

## print(head(up_hit1))
## print(head(up_hit1[up_hit1$a_type=='canon',]))
## print(head(up_hit1[up_hit1$a_type=='prop',]))

canon_hit1_cnts <- gap_cnt_vs_percid(up_hit1_e[up_hit1_e$a_type=='canon',], uniq_pairs, gap_stat_thresh, g_plot_offset*2)
canon_hit1_cnts$pair_type = paste0(canon_hit1_cnts$pair,'_canon')

## print("canon_hit1_cnts")
## print(head(canon_hit1_cnts))

prop_hit1_cnts <- gap_cnt_vs_percid(up_hit1_e[up_hit1_e$a_type=='prop',], uniq_pairs, gap_stat_thresh, g_plot_offset*2)
prop_hit1_cnts$pair_type = paste0(prop_hit1_cnts$pair,'_prop')

## print("prop_hit1_cnts")
## print(head(prop_hit1_cnts))

conf_hit1_cnts <- gap_cnt_vs_percid(up_hit1_e[up_hit1_e$a_type=='conf',], uniq_pairs, gap_stat_thresh, g_plot_offset*2)
conf_hit1_cnts$pair_type = paste0(conf_hit1_cnts$pair,'_conf')

## print("conf_hit1_cnts")
## print(head(conf_hit1_cnts))

canon_hit1_cnts <- canon_hit1_cnts[canon_hit1_cnts$g_cnt > 0,]

canon_hit1_cnts$g_fract_nz = ifelse(canon_hit1_cnts$g_fract > 0, canon_hit1_cnts$g_fract, min_gap_cnt_fn)
canon_hit1_cnts$pair = factor(canon_hit1_cnts$pair,levels=tax.order, ordered=T)
canon_hit1_cnts$g_cut = factor(canon_hit1_cnts$g_cut,levels=gap_stat_thresh, ordered=T)

prop_hit1_cnts <- prop_hit1_cnts[prop_hit1_cnts$g_cnt > 0,]

prop_hit1_cnts$g_fract_nz = ifelse(prop_hit1_cnts$g_fract > 0, prop_hit1_cnts$g_fract, min_gap_cnt_fn)
prop_hit1_cnts$pair = factor(prop_hit1_cnts$pair,levels=tax.order, ordered=T)
prop_hit1_cnts$g_cut = factor(prop_hit1_cnts$g_cut,levels=gap_stat_thresh, ordered=T)

conf_hit1_cnts <- conf_hit1_cnts[conf_hit1_cnts$g_cnt > 0,]

conf_hit1_cnts$g_fract_nz = ifelse(conf_hit1_cnts$g_fract > 0, conf_hit1_cnts$g_fract, min_gap_cnt_fn)
conf_hit1_cnts$pair = factor(conf_hit1_cnts$pair,levels=tax.order, ordered=T)
conf_hit1_cnts$g_cut = factor(conf_hit1_cnts$g_cut,levels=gap_stat_thresh, ordered=T)

################
## set up cummulative percid columns by pair, a_type

up_percid_cumm = NULL
up_gap90_cumm = NULL
top_gaps_df = NULL

leg.labels_id = c()
leg.labels_gap=c()
leg.labels_gcnt=c()
leg.labels_names = c()
pair.names = c()

## print(head(up_hit1_e[,c('pair','a_type')]))
## print(unique(up_hit1_e$a_type))

plot_columns = c('pair','class','a_type','pair','pair_type','pthr_id','qseqid','sseqid','percid','g_med','g_max','stat_nz')

for (pair in unique(up_hit1_e$pair)) {

  pair.names = append(pair.names, pair)
  this_pair = pair_times[pair_times$pair == pair,]

  for (a_type in c('canon','prop','conf')) {

    ## print(sprintf("for %s %s\n",pair, a_type))

    this_pair_atype = sprintf("%s_%s",pair,a_type)
    leg.labels_names = append(leg.labels_names,this_pair_atype)
    ## print(this_pair_atype)
    
    ## subset on pair, a_type
    pair_subset = up_hit1_e[up_hit1_e$pair==pair & up_hit1_e$a_type == a_type,]

    ## print(paste("pair_subset",pair, a_type))
    pair_subset$class = this_pair$class
    
##    print("head(pair_subset)")
##    print(head(pair_subset))

    pair_cumm = pair_subset[order(pair_subset$percid),plot_columns]
    pair_cumm$row_num = 1:nrow(pair_cumm)
    pair_cumm$fract_row = pair_cumm$row_num/nrow(pair_cumm)
    pair_cumm$c_fract_row = 1.0 - pair_cumm$fract_row
    
##    print("head(pair_cumm)")
##    print(head(pair_cumm))

    ## set up label for id with number of hits
    
    ## s_pair = gsub('_v_',':',pair)
    s_pair <- pair_subset[1,]$s_pair
##    print(paste("s_pair:",s_pair))

##    print(paste("id--pair", pair,'atype:',a_type," length_e90:",length(pair_cumm$row_num)))

    this_label = sprintf("%s  \n N = %s",s_pair, format(length(pair_cumm$row_num),big.mark=','))
    this_label = setNames(this_label,this_pair_atype)
##    print(this_label)
    
    leg.labels_id <- append(leg.labels_id,this_label)
    
    samp_rate = trunc(length(pair_cumm$row_num)/target_samp_rate + 0.5)
    
    if (a_type == 'canon' || a_type=='conf' || pair=='ecoli_v_salty' || pair=='yeast_v_sacar') {
        pair_cumm_samp = pair_cumm[pair_cumm$row_num %% samp_rate ==0,]

	## now we have a sampled set of canon, get the same sampled set of prop for _id and gaps
	canon_pthr_ids = pair_cumm_samp$pthr_id
    } else {
        pair_cumm_samp = pair_cumm[pair_cumm$pthr_id %in% canon_pthr_ids,]
    }
    pair_subset_e90 = up_hit1_e90[up_hit1_e90$pair==pair & up_hit1_e90$a_type == a_type,]
    pair_subset_e90$class = this_pair$class
    
    pair_cumm_gap90 = pair_subset_e90[order(pair_subset_e90$stat_nz),plot_columns]
    pair_cumm_gap90$row_num = 1:nrow(pair_cumm_gap90)
    pair_cumm_gap90$fract_row = pair_cumm_gap90$row_num/nrow(pair_cumm_gap90)
    pair_cumm_gap90$c_fract_row = 1.0 - pair_cumm_gap90$fract_row
    
    top_gaps = pair_subset_e90[order(-pair_subset_e90$stat_nz),c('a_type','pair','qseqid','qlen','sseqid','slen','percid','n_gaps','g_med','g_max','stat_nz')]

    if (nrow(top_gaps[top_gaps$stat_nz > long_gap_thresh,]) > 25) {
        top_gaps = top_gaps[top_gaps$stat_nz > long_gap_thresh,]
    } else {
        top_gaps = top_gaps[1:50,]
    }

    top_gaps25 = top_gaps[1:min(nrow(top_gaps),25),]

    top_gaps_df = rbind(top_gaps_df, top_gaps)

    g90_len = length(pair_cumm_gap90$row_num)

    this_label = sprintf("%s  \n N= %s (%.0f%%)",s_pair, 
       format(g90_len,big.mark=','),
       100.0*length(pair_cumm_gap90$row_num)/length(pair_cumm$row_num))

    this_label = setNames(this_label,this_pair_atype)
    leg.labels_gap <- append(leg.labels_gap,this_label)
    
    ge5_len = length(pair_cumm_gap90[pair_cumm_gap90$g_max >= 5,]$row_num)
    
    if (ge5_len/g90_len < 0.1) {
        g90_fmt = "%s  \n N(>90%%)=%s (%.1f%%)"
    } else {
        g90_fmt = "%s  \n N(>90%%)=%s (%.0f%%)"
    }

    gcnt_label = sprintf(g90_fmt,s_pair,ge5_len,100.0*ge5_len/g90_len)
    gcnt_label = setNames(gcnt_label, this_pair_atype)
    leg.labels_gcnt <- append(leg.labels_gcnt, gcnt_label)

    if (! opt$type) {
      samp_rate = 50
    } else {
      samp_rate = 2
    }
    
    if (length(pair_cumm_gap90$row_num)/samp_rate > target_samp_rate) {
       samp_rate = trunc(length(pair_cumm_gap90$row_num)/target_samp_rate + 0.5)
    }
    
    pair_cumm_gap90_samp = pair_cumm_gap90[(pair_cumm_gap90$row_num %% samp_rate==0) | (pair_cumm_gap90$qseqid %in% top_gaps25$qseqid),]

    if (pair == 'ecoli_v_salty' || pair == 'yeast_v_sacar' || a_type == 'conf') {
        pair_cumm_gap90_samp = pair_cumm_gap90[(pair_cumm_gap90$row_num %% samp_rate == 0)| (pair_cumm_gap90$qseqid %in% top_gaps25$qseqid) ,]
    } else {
        pair_cumm_gap90_samp = pair_cumm_gap90[pair_cumm_gap90$pthr_id %in% canon_pthr_ids,]
    }
    
    ## cat(sprintf("%s %s: %d\n",pair, a_type, length(pair_cumm_gap90_samp$row_num)))
    
    up_percid_cumm = rbind(up_percid_cumm,pair_cumm_samp)
    
    up_gap90_cumm = rbind(up_gap90_cumm,pair_cumm_gap90_samp)
  }
}

## print("unique pair_type")
## print(unique(up_gap90_cumm$pair_type))

pair.names = unique(pair.names)
## pair_times = pair_times[pair_times$name %in% pair.names,]
## print("pair.names")
## print(pair.names)

pair.names.N  = length(tax.order)


leg.labels_levels = tax.order
leg.labels_c = paste0(tax.order, '_canon')
leg.labels_p = paste0(tax.order, '_prop')
leg.labels_cnf = paste0(tax.order, '_conf')

leg.labels_levels3 = append(leg.labels_c, leg.labels_p)
leg.labels_levels3 = append(leg.labels_levels3, leg.labels_cnf)

## print("leg.labels_levels3")
## print(leg.labels_levels3)

leg.labels_names = factor(leg.labels_names,levels=leg.labels_levels3,ordered=TRUE)

## print("leg.labels_levels3/names/id/gap")
## leg.labels_levels3
## print("names:")
## leg.labels_names
## print("id/gap")
## leg.labels_id
## leg.labels_gap
## leg.labels_gcnt

up_percid_cumm$pair = factor(up_percid_cumm$pair,levels=tax.order, ordered=T)
up_percid_cumm$pair_type = factor(up_percid_cumm$pair_type,levels=leg.labels_levels3, ordered=T)
up_percid_cumm$class = factor(up_percid_cumm$class,levels=class.order, ordered=T)

canon_hit1_cnts$pair_type = factor(canon_hit1_cnts$pair_type,levels=leg.labels_levels3, ordered=T)
prop_hit1_cnts$pair_type = factor(prop_hit1_cnts$pair_type,levels=leg.labels_levels3, ordered=T)
conf_hit1_cnts$pair_type = factor(conf_hit1_cnts$pair_type,levels=leg.labels_levels3, ordered=T)

## remove zeros for log plot
up_gap90_cumm$stat_nz = up_gap90_cumm$stat_nz

if (nrow(up_gap90_cumm[up_gap90_cumm$stat_nz<0.1,]) > 0) {
   up_gap90_cumm[!is.na(up_gap90_cumm$stat_nz) & up_gap90_cumm$stat_nz<0.1,]$stat_nz <- 0.1
}

up_gap90_cumm$pair = factor(up_gap90_cumm$pair,levels=tax.order, ordered=T)
up_gap90_cumm$pair_type = factor(up_gap90_cumm$pair_type,levels=leg.labels_levels3, ordered=T)
up_gap90_cumm$class = factor(up_gap90_cumm$class,levels=class.order, ordered=T)

source("set_color_a90.R")

pt_names=as.character(pair_times$pair)

a.colors = a.colors[1:length(pt_names)]
a.colors = setNames(a.colors,pt_names)
a.colors['ecoli_v_salty'] = "grey1"
a.colors['yeast_v_sacar'] = "grey1"

a.colors_c = setNames(a.colors, paste0(pt_names,'_canon'))
a.colors_p = setNames(a.colors, paste0(pt_names,'_prop'))
a.colors_cnf = setNames(a.colors, paste0(pt_names,'_conf'))

a.colors3 = append(a.colors_c, a.colors_p)
a.colors3 = append(a.colors3, a.colors_cnf)

## associate pair names with symbols

pt_names_c = paste0(pt_names,"_canon")
pt_names_p = paste0(pt_names,"_prop")
pt_names_cnf = paste0(pt_names,"_conf")

## p.sym=setNames(tax.sym[pair_times$class],pt_names)

p.sym_c = setNames(tax.sym[as.character(pair_times$class)],pt_names_c)
p.sym_p = setNames(tax.sym[as.character(pair_times$class)],pt_names_p)
p.sym_cnf = setNames(tax.sym[as.character(pair_times$class)],pt_names_cnf)

p.sym3 = append(p.sym_c, p.sym_p)
p.sym3 = append(p.sym3, p.sym_cnf)

p.lines = setNames(rep('solid',length(pt_names)),pt_names)
p.lines['ecoli_v_salty'] = 'longdash'
p.lines['yeast_v_sacar'] = 'longdash'

p.lines_c = p.lines
p.lines_c = setNames(p.lines_c,pt_names_c)
p.lines_p = p.lines
p.lines_p = setNames(p.lines_p,pt_names_p)
p.lines_cnf = p.lines
p.lines_cnf = setNames(p.lines_cnf,pt_names_cnf)

p.lines3 = append(p.lines_c, p.lines_p)
p.lines3 = append(p.lines3, p.lines_cnf)

# p.sym_g['ecoli_v_salty_prop'] = 17

sb.color_id   <- scale_color_manual(name='id_legend', values=a.colors3, labels=leg.labels_id)
sb.color_g <- scale_color_manual(name='gap_legend', values=a.colors3, labels=leg.labels_gap)
sb.color_gcnt <- scale_color_manual(name='gcnt_legend', values=a.colors3, labels=leg.labels_gcnt)

sb.shape_id <- scale_shape_manual(name='id_legend', values=p.sym3, labels=leg.labels_id)
sb.shape_g  <- scale_shape_manual(name='gap_legend', values=p.sym3, labels=leg.labels_gap)
sb.shape_gcnt <- scale_shape_manual(name='gcnt_legend', values=p.sym3, labels=leg.labels_gcnt)

p.alpha = setNames(rep(1.0,length(pt_names)),pt_names)
p.alpha['ecoli_v_salty'] = 0.33
p.alpha['yeast_v_sacar'] = 0.33

p.alpha_id <- p.alpha
p.alpha_c <- setNames(p.alpha_id, pt_names_c)
p.alpha_p <- setNames(p.alpha, pt_names_p)
p.alpha_cnf <- setNames(p.alpha, pt_names_cnf)

p.alpha3 <- append(p.alpha_c, p.alpha_p)
p.alpha3 <- append(p.alpha3, p.alpha_cnf)

sb.alpha_id <- scale_alpha_manual(name='id_legend', values=p.alpha3, labels=leg.labels_id)
sb.alpha_g <- scale_alpha_manual(name='gap_legend', values=p.alpha3, labels=leg.labels_gap)
sb.alpha_gcnt <- scale_alpha_manual(name='gcnt_legend', values=p.alpha3, labels=leg.labels_gcnt)

sym.size = 1
p.size = setNames(rep(sym.size, length(pt_names)),pt_names)
p.size['ecoli_v_salty'] = 1.33
p.size['yeast_v_sacar'] = 1.33

p.size_c <- setNames(p.size, pt_names_c)
p.size_p <- setNames(p.size, pt_names_p)
p.size_cnf <- setNames(p.size, pt_names_cnf)

p.size3 = append(p.size_c, p.size_p)
p.size3 = append(p.size3, p.size_cnf)

sb.size_g <- scale_size_manual(name='gap_legend', values=p.size3, labels=leg.labels_gap)

sb.line_gcnt <- scale_linetype_manual(name='gcnt_legend',values=p.lines3, labels=leg.labels_gcnt)

scalex.id <- scale_x_continuous(breaks=seq(0.0,1.0,0.2),limits=c(0.0,1.0))
scalex.gap <- scale_x_continuous(breaks=seq(0.0,1.0,0.2),limits=c(0.0,1.02))
scalex.gap_cnt <- scale_x_continuous(breaks=100*seq(0.5,1.0,0.1),limits=100*c(0.5,1.0))

## y-axis for identity
scaley.id <- scale_y_continuous("percent identity", breaks=100.0*seq(0.4, 1.0, 0.1), limits=100*c(0.4,1.0))

scaley.gap <- scale_y_log10(ygap_label)

l.yvals <- c(0.1,1.0,5.0,10.0,20.0,50.0,100.0)
l.yvals100 <- c(0.05,0.1,0.2,0.5,1.0,2.0,5.0,10.0,20.0,50.0,100.0)
ly.labels <- c("0","1","5","10","20","50","100")
ly.labels100 <- c("0","0.1","0.2","0.5","1","2","5","10","20","50","100")

## scaley.log_g <- scale_y_log10(ygap_label, breaks=l.yvals,minor_breaks=l.yvals,labels=ly.labels,limits=c(0.1,1000.0))
scaley.log_g <- scale_y_log10(ygap_label, breaks=l.yvals,minor_breaks=l.yvals,labels=ly.labels,limits=c(0.1,1000.0))
scaley.log_gr <- scale_y_log10(ygap_label, breaks=l.yvals,minor_breaks=l.yvals,labels=ly.labels,position='right',limits=c(0.1,1000.0))

scaley.log_gcnt <- scale_y_log10(ygap_label_gcnt, breaks=l.yvals100,minor_breaks=l.yvals100,limits=c(0.05,100.0))
scaley.log_gcntr <- scale_y_log10(ygap_label_gcnt, breaks=l.yvals100,minor_breaks=l.yvals100,labels=ly.labels100,position='right',limits=c(0.05,100.0))

theme_set(theme_linedraw(base_size=14))
theme.base <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_markdown(face='plain', hjust=0),
  plot.caption=element_markdown(size=6,hjust=0),
  axis.text.x = element_markdown(size=12),
  axis.text.y = element_markdown(size=12),
  legend.text=element_markdown(size=9),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.text.align = 0,
  legend.title=element_blank(),
  legend.margin=margin(c(0.2,5,5)),
  strip.background=element_rect(fill='white'), strip.text=element_text(color='black'))

theme.leg_id = theme.base + theme(legend.position=c(0.98,0.02),legend.justification=c('right','bottom'))
theme.leg_gap = theme.base + theme(legend.position=c(0.02,0.98))
theme.leg_cnt = theme.base + theme(legend.position=c(0.02,0.02),legend.justification=c('left','bottom'))

size.guide = guides(color = guide_legend(override.aes = list(size = 2)))

## pdf(file=fig1i_file,width=8, height=4)

## p_id_c -- canon
p_id_c <-ggplot(up_percid_cumm[up_percid_cumm$a_type=='canon',],aes(x=fract_row, y=percid*100.0)) +
  theme.leg_id + scaley.id + scalex.id + ## scale_x_reverse() + 
  geom_hline(yintercept=90,linetype='longdash')+
  geom_point(aes(color=pair_type,shape=pair_type,alpha=pair_type),size=sym.size) + 
  sb.color_id + sb.shape_id + sb.alpha_id + size.guide +
  labs( x='fraction of alignments',
        subtitle='A. canonical isoforms (identity)'
	)

## p_id_p -- prop
p_id_p <-ggplot(up_percid_cumm[up_percid_cumm$a_type=='prop',],aes(x=fract_row, y=percid*100.0)) +
  theme.leg_id + scaley.id + scalex.id + ## scale_x_reverse() + 
  geom_hline(yintercept=90,linetype='longdash')+
  geom_point(aes(color=pair_type,shape=pair_type,alpha=pair_type),size=sym.size) + 
  sb.color_id + sb.shape_id + sb.alpha_id + size.guide +
##  xlab('fraction of alignments')
  labs( x='fraction of alignments',
        subtitle='D. proposed isoforms (identity)'
	)

## p_id_cnf -- conf

fig5_sub_a = 'G. confirmed isoforms (identity)'
fig5_sub_b = 'H. confirmed isoforms (gaps)'
fig5_sub_c = 'I. confirmed isoforms (\u2265 5 gap count)'

if (opt$supp) {
   fig5_sub_a = 'A. confirmed isoforms (identity)'
   fig5_sub_b = 'B. confirmed isoforms (gaps)'
   fig5_sub_c = 'C. confirmed isoforms (\u2265 5 gap count)'
}


p_id_cnf <-ggplot(up_percid_cumm[up_percid_cumm$a_type=='conf',],aes(x=fract_row, y=percid*100.0)) +
  theme.leg_id + scaley.id + scalex.id + ## scale_x_reverse() + 
  geom_hline(yintercept=90,linetype='longdash')+
  geom_point(aes(color=pair_type,shape=pair_type,alpha=pair_type),size=sym.size) + 
  sb.color_id + sb.shape_id + sb.alpha_id + size.guide +
##  xlab('fraction of alignments')
  labs( x='fraction of alignments',
        subtitle=fig5_sub_a
	)

## 90% gaps with axis label on right
## p_gap90r_c -- canon

p_gap_alpha=0.5

p_gap90r_c <-ggplot(up_gap90_cumm[up_gap90_cumm$a_type=='canon',],aes(x=fract_row, y=stat_nz)) +
  theme.leg_gap + scaley.log_gr + scalex.gap + ## scale_x_reverse() + 
  geom_point(aes(color=pair_type,shape=pair_type,alpha=pair_type,size=pair_type)) + 
  annotate('rect',xmin=0.98,xmax=1.019999,ymin=50,ymax=999,alpha=0.0,color='red') +
  geom_hline(yintercept=5,linetype='dotdash')+
  sb.color_g + sb.shape_g + sb.alpha_g + sb.size_g + size.guide +
  labs( x='f\'n of alignments > 90% identical',
        subtitle='B. canonical isoforms (gaps)'
	)

p_gap90r_p <-ggplot(up_gap90_cumm[up_gap90_cumm$a_type=='prop',],aes(x=fract_row, y=stat_nz)) +
  theme.leg_gap + scaley.log_gr + scalex.gap + ## scale_x_reverse() + 
  geom_point(aes(color=pair_type,shape=pair_type,alpha=pair_type,size=pair_type)) + 
  annotate('rect',xmin=0.98,xmax=1.019999,ymin=50,ymax=999,alpha=0.0,color='red') +
  geom_hline(yintercept=5,linetype='dotdash')+
  sb.color_g + sb.shape_g + sb.alpha_g + sb.size_g + size.guide +
  labs( x='f\'n of alignments > 90% identical',
        subtitle='E. proposed isoforms (gaps)'
	)

p_gap90r_cnf <-ggplot(up_gap90_cumm[up_gap90_cumm$a_type=='conf',],aes(x=fract_row, y=stat_nz)) +
  theme.leg_gap + scaley.log_gr + scalex.gap + ## scale_x_reverse() + 
  geom_point(aes(color=pair_type,shape=pair_type,alpha=pair_type,size=pair_type)) + 
  annotate('rect',xmin=0.98,xmax=1.019999,ymin=50,ymax=999,alpha=0.0,color='red') +
  geom_hline(yintercept=5,linetype='dotdash')+
  sb.color_g + sb.shape_g + sb.alpha_g + sb.size_g + size.guide +
  labs( x='f\'n of alignments > 90% identical',
        subtitle=fig5_sub_b
	)

p_gapcnt_r_c <- ggplot(canon_hit1_cnts,aes(x=100.0*(id_q_percid+g_plot_offset), y=g_fract_nz, color=pair_type, shape=pair_type, linetype=pair_type)) +
  sb.color_gcnt + scalex.gap_cnt + sb.shape_gcnt + theme.leg_cnt + scaley.log_gcntr + sb.line_gcnt +
  labs(x='percent identity', 
      subtitle='C. canonical (\u2265 5 gap count)') + 
  geom_vline(xintercept=90,linetype='longdash')+
  geom_point() + geom_line()

p_gapcnt_r_p <- ggplot(prop_hit1_cnts,aes(x=100.0*(id_q_percid+g_plot_offset), y=g_fract_nz, color=pair_type, shape=pair_type, linetype=pair_type)) +
  sb.color_gcnt + scalex.gap_cnt + sb.shape_gcnt + theme.leg_cnt + scaley.log_gcntr + sb.line_gcnt +
  geom_vline(xintercept=90,linetype='longdash')+
  geom_point() + geom_line() + 
  labs(x='fraction identical', subtitle='F. proposed (\u2265 5 gap count)')

p_gapcnt_r_cnf <- ggplot(conf_hit1_cnts,aes(x=100.0*(id_q_percid+g_plot_offset), y=g_fract_nz, color=pair_type, shape=pair_type, linetype=pair_type)) +
  sb.color_gcnt + scalex.gap_cnt + sb.shape_gcnt + theme.leg_cnt + scaley.log_gcntr + sb.line_gcnt +
  geom_vline(xintercept=90,linetype='longdash')+
  geom_point() + geom_line() +
  labs(x='percent identity', subtitle=fig5_sub_c)

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

## p_id_gap <- ((p_id_c | p_gap90r_c) / (p_id_p | p_gap90r_p)) + plot_annotation(tag_level='A',tag_suffix='.')

if (opt$supp) {
   p_id_gap <- (p_id_c | p_gap90r_c | p_gapcnt_r_c) / (p_id_p | p_gap90r_p |  p_gapcnt_r_p)
   p_id_gap_supp <- (p_id_cnf | p_gap90r_cnf | p_gapcnt_r_cnf)
} else {
   p_id_gap <- ((p_id_c | p_gap90r_c | p_gapcnt_r_c) / (p_id_p | p_gap90r_p |  p_gapcnt_r_p) / (p_id_cnf | p_gap90r_cnf | p_gapcnt_r_cnf))
}

if (opt$supp) {
  if (! opt$pub) {
    p_id_gap <- p_id_gap / doc_panel + plot_layout(heights=c(10, 10, 0.1))
    ggsave(file=fig_file, plot=p_id_gap, width=12.5, height=8.5, device=cairo_pdf)
    p_id_gap_supp <- p_id_gap_supp / doc_panel + plot_layout(heights=c(10, 0.1))
    ggsave(file=fig_file_supp, plot=p_id_gap_supp, width=12.5, height=4.6, device=cairo_pdf)
  } else {
    ggsave(file=fig_file, plot=p_id_gap, width=12.5, height=8.4, device=cairo_pdf)
    ggsave(file=fig_file_supp, plot=p_id_gap_supp, width=12.5, height=4.5, device=cairo_pdf)
  }
} else {
  if (! opt$pub) {
    p_id_gap <- p_id_gap / doc_panel + plot_layout(heights=c(10, 10, 10, 0.1))
    ggsave(file=fig_file, plot=p_id_gap, width=12.5, height=12.5, device=cairo_pdf)
  } else {
    ggsave(file=fig_file, plot=p_id_gap, width=12.5, height=12.0, device=cairo_pdf)
  }
}

if (opt$txt) {
   gaps_fd = file(top_gaps_file,open='w')

   write(paste0("# ",paste0(colnames(top_gaps_df),collapse='\t',sep='')),gaps_fd)

   for (a_type in c('canon','prop','conf')) {
       for (pair in tax.order) {
           these_rows = top_gaps_df[top_gaps_df$pair==pair & top_gaps_df$a_type==a_type,]
           write(sprintf(">%s %s %d",a_type,pair,NROW(these_rows)),gaps_fd)
           write.table(these_rows,sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE,gaps_fd)
       }      
   }      

   close(gaps_fd)
}
