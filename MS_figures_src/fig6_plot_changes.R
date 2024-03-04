#!/usr/bin/env Rscript --vanilla

################
##
## 08-Dec-2023
##
## plot_changes.R history_new_changes_byspecies.tsv
##
## plot history of ortho2tree proposed changes
##

library('ggplot2')
library('dplyr')
library('optparse')
library('reshape2')

library(patchwork)
library('RColorBrewer')

option_list = list(
	    make_option(c("--file","-f"),type='character',action='store',default='',help='history file'),
	    make_option(c("-p","--pdf"),type='character',action='store',default='fig6', help='optional pdf file name [default= \"%default\"]'),
	    make_option(c("-P","--pub"),action='store_true',default=FALSE, help='remove plot documentation')
	    )

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (opt$file=='') {
    print_help(opt_parser)
    cat(" plen_diff.R  output_gc231004 \n\n")
    quit(save='no',status=1)
}

d_args<- commandArgs(trailingOnly=TRUE)

p.name<-'plot_changes.R'
plabel=paste(c(p.name,d_args,"\n",date()),collapse=' ', sep=' ')

plot_name = "fig6"
if (! is.null(opt$pdf)) {
   plot_name = opt$pdf
}

ensembl_rel = c('2023_01'='107', '2023_02'='108','2023_03'='108','2023_04'='109','2023_05'='109','2024_01'='109')

## 107* new dog
## 108->109 new horse


################
## names for output plots and text files

fig_file=paste0(plot_name,"_cumm_prop")
if (opt$pub) {
   fig_file = paste0(fig_file,"_pub")
}
fig_file = paste0(fig_file,".pdf")

fig_file

################

read_fields <- c('release','cnt','taxon')

################
## read in the data

prop_df<-read.table(file=opt$file,header=FALSE,sep='\t',quote='',col.names=read_fields)

prop_df$taxon = tolower(prop_df$taxon)

## prop_df <- aggregate(cnt ~ taxon, prop_df, FUN=cumsum)

total_df <- prop_df
total_df <- total_df %>% group_by(release) %>% dplyr::mutate(totcnt=sum(cnt))
total_df <- total_df[total_df$taxon=='human',]

total_df$taxon = 'total'
total_df$cnt = total_df$totcnt
total_df$cumcnt <- cumsum(total_df$totcnt)

prop_df <- prop_df %>% group_by(taxon) %>% dplyr::mutate(cumcnt=cumsum(cnt))

target_taxa = unique(prop_df$taxon)

## print(target_taxa)

names.4 = tolower(c('GORGO','MOUSE','RAT','BOVIN'))
names.5 = tolower(c('GORGO','MOUSE','RAT','BOVIN','HUMAN'))
names.p3 = tolower(c('CANLF','MONDO','PANTR'))
names.8 = c(names.5, names.p3)

names.mam_x = tolower(c('FELCA','HORSE','MACMU','ORNAN','PIG'))

taxa_levels = c('human',names.4, names.p3, names.mam_x)
print("taxa_levels")
print(taxa_levels)

## print("prop_df")
## prop_df

## here, we need to match the colors (as much as possible) with Fig 1, 4, and 5.

## colors for the 5 taxa used in the other plots:
## 
## 
a.colors <- brewer.pal(8,'Dark2') # 'Dark2', 'Set2', 'Paired'
a.color_save2 = a.colors[2]
a.color_save6 <- a.colors[6]
a.colors[6] <- a.colors[7]

a.colors[2:5] <- a.colors[3:6]
a.colors[6] <- a.color_save2
a.colors[7] <- a.color_save6

a.colors8 <- setNames(a.colors,names.8)

## set the colors
b.colors5 <- brewer.pal(5,'Dark2') # 'Dark2', 'Set2', 'Paired'
b.colors5 <- setNames(b.colors5, names.mam_x)

a.colors13 <- append(a.colors8, b.colors5)

## print(a.colors13)

print("head(prop_df)")
head(prop_df)

## a.colors13

a.sym <- c(rep(1,8),rep(2,8))
a.sym <- a.sym[1:length(target_taxa)]
a.sym <- setNames(a.sym, taxa_levels)
a.sym['human'] = 0

s.color <- scale_color_manual(values=a.colors13, labels=taxa_levels)
s.sym <- scale_shape_manual(values=a.sym, labels=taxa_levels)

prop_df$taxon = factor(prop_df$taxon, levels=taxa_levels, ordered=TRUE)

# rel_labels = c('2023_01'='2023_01\n(8 org.)','2023_02'='2023_02\n(13 org.)','2023_03'='2023_03\n(13 org.)','2023_04'='2023_04 \n(13 org.*)','2023_05'='2023_05  \n(13 org.**)','2024_01'='2024_01\n(35 org.)')

rel_labels = c('2023_01'='2023_01\n8 org. ','2023_02'='2023_02\n13 org.','2023_03'='2023_03\n13 org.','2023_04'='2023_04\n13 org.*','2023_05'='2023_05\n13 org.**','2024_01'='2024_01\n35 org.')

rel_labels1= c('2023_01'='2023_01 (8)','2023_02'='2023_02 (13)','2023_03'='2023_03 (13)','2023_04'='2023_04 (13*)','2023_05'='2023_05 (13**)','2024_01'='2024_01 (35)')

scale.x = scale_x_discrete('release',labels=rel_labels)


theme_set(theme_linedraw(base_size=14))
theme.base <- theme(panel.background=element_rect(colour='black', linewidth=1.0),
  panel.grid.major=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  panel.grid.minor=element_line(colour='darkgrey', linewidth=0.4,linetype='dashed'),
  plot.title=element_text(face='plain', hjust=0),
  plot.caption=element_text(size=6,hjust=0),
  axis.text.x = element_text(size=9,angle=45,vjust=1,hjust=1),
  axis.text.y = element_text(size=12),
  legend.text=element_text(size=9),
  legend.key=element_blank(),
  legend.background=element_rect(fill='white', color='black',linetype='solid',linewidth=0.4),
  legend.justification=c(0,1),
  legend.text.align = 0,
  legend.title=element_blank(),
  legend.margin=margin(c(0.5,5,5,5)))

theme.leg = theme.base ## + theme(legend.position=c(0.10,0.95))

## text.x_axes = sprintf("%s\n%s",unique(prop_df$release),ensembl_rel[unique(prop_df$release)])
## print(text.x_axes)


################
## convert wide (ldelta_t, l_delta_t_tr) to long (ldelta_df, ldelta_df_tr)

xaxis_title="release number"
yaxis_title="changes"

yA.breaks = c(0,1000,2000,3000,4000)
yA.breaks_l = format(yA.breaks,big.mark=',')

yB.breaks = c(0,10000,20000,30000)
yB.breaks_l = format(yB.breaks,big.mark=',')

################
## all changes
p_delta <-ggplot(prop_df,aes(x=release, y=cumcnt, color=taxon, group=taxon)) +
  theme.leg +   geom_point(aes(shape=taxon)) +
  geom_line() +
  s.color + s.sym +
  scale_y_continuous(name=yaxis_title,breaks=yA.breaks, labels=yA.breaks_l) + 
  labs(subtitle='A. changes by taxon')

p_total <-ggplot(total_df,aes(x=release, y=cumcnt, group=taxon)) +
  theme.leg +   geom_point(shape=1) +
  geom_line() + 
  scale_y_continuous(name=yaxis_title,limits=c(0,30000),breaks=yB.breaks, labels=yB.breaks_l) +
  scale.x +
  labs(subtitle='B. total changes')

doc_panel = ggplot() + theme_void() + labs(caption=plabel)

if (! opt$pub) {
  f_plot <- p_delta / p_total / doc_panel + plot_layout(heights=c(6.5, 2, 0.1))
  ggsave(file=fig_file, plot=f_plot, width=5, height=7.5)
} else {
  f_plot <- p_delta / p_total + plot_layout(heights=c(6.5, 2))
  ggsave(file=fig_file, plot=f_plot, width=5, height=7.0)
}

