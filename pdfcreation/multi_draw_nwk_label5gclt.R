#!/usr/bin/env Rscript --vanilla

################
## multi_draw_nwk_label5gclt.r tstamp pdf_dir tree_data aln_data lab_data .lab_suff
## reads gene_centric file if sys.env(GC_FILE) is defined and file exists
##
## multi_draw_nwk_label5.r tstamp pdf_dir tree_data aln_data lab_data .lab_suff
##
## derived from:
## multi_draw_nwk_label4d.r list_file pdf_dir tree_data aln_data lab_data .lab_suff
##

## (1) given tstamp, reads output_changes{tstamp} and output_confirmed{tstamp}
##     to find pthr_id, accs changed and confirmed
## (2) for all changes, report ntax, canonical, proposed clade cost, and 
##     produce label (S1, S2, etc)
## (3) for confirmed, report ntax, clade cost, produce label (C1, C2, etc)
##
## (4) if only confirmed, write to _mltmC.pdf, otherwise _mltm.pdf
##

## multi_draw_nwk_label4d.r -- read a output file, build trees for each
## of the rows in the file
##
## ../draw_nwk_label4d.r query_subject_msa.nwk_l  (also looks for .aln_faX, produces _msa2g.pdf)
## ../draw_nwk_label4d.r PTHR012345.nwk 4d PTHR012345.lab PTHR012345.lab2 PTHR012345.faX
##
################

## install.packages()
library("ggplot2", quietly=TRUE)
## library("phangorn", warn.conflicts=FALSE, quietly=TRUE)
library("ape", quietly=TRUE)
## library("grid")

## BiocManager::install()
library("Biostrings", warn.conflicts=FALSE, quietly=TRUE)
library("treeio", warn.conflicts=FALSE, quietly=TRUE)
library("ggtree", warn.conflicts=FALSE, quietly=TRUE)
## library("ggplotify")

color_map = c('prop_iso/MANE_sel'='aquamarine3', 'prop_iso'='green3', 'prop_iso0'='darkorange','prop_iso2'='darkgreen','canon_good'='blue', 'canon_good/MANE_sel'='blue4','sugg_canon'='green3', 'canon'='black', 'canon_orig'='black', 'other_iso'='grey','MANE_sel'='blue','canon_bad'='red', 'curr_canon'='red', 'canon_bad/MANE_sel'= 'orange', 'canon/MANE_sel'='purple', 'sp'='black', 'tr'='black', 'iso'='grey','iso/MANE_sel'='orange')

msa.alpha=0.1

## read_tree_plot() creates a single plot

read_tree_plot <- function(pthr_id, tree, tax_labels, aln_dir, pdf_dir, pdf_suff, title_str) {

    faX.file <- paste0(aln_dir,'/',pthr_id,'.faX')

    ## print(c(pdf_suff, label_file, faX.file))

    ## print(tree)
    ## str(tree)
    ## print(tree$tip.label)
    ## length(tree$tip.label)

    ltz <- tree$edge.length < 0.0
    tree$edge.length[ltz] = 0.00001

    p_tree = ggtree(tree) %<+% tax_labels + geom_tiplab(aes(label=label_str,color=type),offset=0.004) + scale_color_manual(values=color_map) + theme(legend.position='none') ## + geom_label(aes(label=rank),fill='lightgreen',vjust=1.0) + ## geom_label(aes(label=C), hjust=-.2, vjust=1., alpha=msa.alpha) +
    geom_tippoint(aes(shape=rank),position=position_nudge(x=0.002),fill='lightgreen')

    if (file.exists(faX.file)) {
      fa<-treeio::read.fasta(faX.file)
      fa <- as.matrix(fa)
      ## print(sprintf("%s faX file read",faX.file))
      ## print(labels(fa))
      ## print(str(fa))
      return (msaplot(p_tree,faX.file,width=1.0,offset=0.2, height=0.50) + ggtitle(title_str))
    } else {
      print(sprintf("%s faX file not found",faX.file))
      return(p_tree + ggtitle(title_str))
    }
}

args <- commandArgs(trailingOnly=TRUE)

tstamp = args[1]
## multi_file = 'output.pthr100'

change_file=paste0('output_changes',tstamp)
confirm_file=paste0('output_confirm',tstamp)

pdf_dir = 'pdf_data'
if (! is.na(args[2])) {
   pdf_dir = args[2]
}

if (! file.exists(pdf_dir)){
    dir.create(pdf_dir)
}

tree_dir = 'tree_data'
if (! is.na(args[3])) {
   tree_dir = args[3]
}

aln_dir = 'aln_data'
if (! is.na(args[4])) {
   aln_dir = args[4]
}

label_dir = 'lab_data'
if (! is.na(args[5])) {
   label_dir = args[5]
}

label_suff = '.lab_ltm'
if (! is.na(args[6])) {
   label_suff = args[6]
}
label_suff0 = '.lab_lt'

gc_df <- NULL
gc_file <- Sys.getenv("GC_FILE")
if (gc_file != "" && file.exists(gc_file)) {

   gc_df <- read.csv(gc_file,header=TRUE)
   gc_df <- gc_df[,c('acc','groupid')]

   ## genecentric can have duplicate acc/groupid rows
   ## for later joins, they must be unique, so do it now

   gc_df <- unique(gc_df)
   ## head(gc_df)
}

change_col_names = c("pthr_id","taxon","canon_acc","canon_len","prop_acc","prop_len","rank_score","canon_cost","prop_cost","n_sp","n_tax","n_canon", "wn_canon", "scaled_prop_f","scaled_p_cost","p_len_diff","clade","prop_canon", "MANE")

confirm_col_names = c("pthr_id","canon_cost","n_sp","n_tax","clade","clade_members", "MANE")

## cat(paste(change_col_names),'\n')
## cat(length(change_col_names),'\n')
change_data <- read.table(change_file, header=FALSE, sep='\t', stringsAsFactors=FALSE,col.names=change_col_names)
## cat(length(change_data[1,]),'\n')
## head(change_data)

## cat(paste(confirm_col_names),'\n')
## cat(length(confirm_col_names),'\n')

if (file.exists(confirm_file)) {
   confirm_data <- read.table(confirm_file, header=FALSE, sep='\t', stringsAsFactors=FALSE,col.names=confirm_col_names)
} else {
   confirm_data = NULL
}
## cat(length(confirm_data[1,]),'\n')
## head(confirm_data)

if (NROW(confirm_data) > 0) {
  confirm_data$rank_score = 100.0
  confirm_data$prop_cost = confirm_data$canon_cost
}

u_change_pthr = unique(change_data$pthr_id)

cat("uniq change pthr: ",length(u_change_pthr),"\n")

u_conf_pthr = unique(confirm_data$pthr_id)
cat("uniq conf pthr: ",length(u_conf_pthr),"\n")

u_out_data = unique(c(u_change_pthr, u_conf_pthr))

cat("uniq change/conf union pthr: ",length(u_out_data),"\n")
cat("change not in conf (pthr): ",length(setdiff(u_change_pthr, u_conf_pthr)),"\n")
cat("conf not in change (pthr): ",length(setdiff(u_conf_pthr, u_change_pthr)),"\n")

## head(u_out_data)

for (i in 1:length(u_out_data)) {

    pthr_id = u_out_data[i]
    ## print(pthr_id)

    change_out_data = change_data[change_data$pthr_id==pthr_id,]

    ## build up the title string for changes
    title_str = ''

    taxa_clade_map = NULL

    change_flag = FALSE

    ## collect changed clades

    change_out_clades = unique(change_out_data$clade)

    if (length(change_out_clades) > 0) {
        change_flag = TRUE
        for (j in 1:length(change_out_clades)) {
    	    cl_uniq_data = change_out_data[change_out_data$clade==change_out_clades[j],]

	    ## cat("cl_uniq_data:\n")
	    ## print(head(cl_uniq_data))

	    ## extract canonical acc's from changes
	    canon_changed_df = do.call(rbind,strsplit(cl_uniq_data$canon_acc,'|',fixed=TRUE))
	    colnames(canon_changed_df) = c('db','acc')
	    
	    ## make canon_df to accumulate canonicals that were changed (will be added to this_taxa_clade_map)
	    canon_df = data.frame('taxa'=sprintf("%s:%s",cl_uniq_data$taxon,canon_changed_df[,2]),'clade'='0','label'=sprintf("S%d",j),'tx_name'=cl_uniq_data$taxon,'acc'=canon_changed_df[,2])

	    ##cat("canon_changed_df 1:\n")
	    ## print(canon_changed_df)

	    ## get info on proposed canonicals in clade
   	    this_out_data = cl_uniq_data[1,]
	    these_taxa = strsplit(this_out_data$prop_canon,',')[[1]]
	    this_taxa_clade_map = data.frame('taxa'=these_taxa, 'clade'=change_out_clades[j],'label'=sprintf('S%d',j))

	    this_map_df = do.call(rbind,strsplit(this_taxa_clade_map$taxa,':'))
	    colnames(this_map_df) = c('tx_name','acc')
	    this_taxa_clade_map = cbind(this_taxa_clade_map,this_map_df)

	    ## add bad canonicals (canon_df) to this_taxa_clade_map
	    this_taxa_clade_map = rbind(canon_df, this_taxa_clade_map)

	    ## add all labeled taxa for tree
	    taxa_clade_map <- rbind(taxa_clade_map,this_taxa_clade_map)

	    ## add this clade to title
	    this_title <- sprintf(" S%d: n_tax: %d rank_score: %.2f canon_cost: %.4f prop_cost: %.4f\n",
               	         j,this_out_data$n_tax,this_out_data$rank_score,this_out_data$canon_cost,this_out_data$prop_cost)
	    title_str = paste0(title_str,this_title)
	}
    }

    ## collect confirmed clades

    confirm_out_data = confirm_data[confirm_data$pthr_id==pthr_id,]

    confirm_out_clades = unique(confirm_out_data$clade)
    if (length(confirm_out_clades) > 0) {
        confirm_out_clades = unique(confirm_out_data$clade)
        for (j in 1:length(confirm_out_clades)) {
    	    cl_uniq_data = confirm_out_data[confirm_out_data$clade==confirm_out_clades[j],]
    	    this_out_data = cl_uniq_data[1,]

	    ## extract canonicals
	    these_taxa = strsplit(this_out_data$clade_members,',')[[1]]
	    this_taxa_clade_map = data.frame('taxa'=these_taxa, 'clade'=change_out_clades[j],'label'=sprintf('C%d',j))

	    ## extract accessions 
	    this_map_df = do.call(rbind,strsplit(this_taxa_clade_map$taxa,':'))
	    colnames(this_map_df) = c('tx_name','acc')
	    taxa_clade_map <- rbind(taxa_clade_map,cbind(this_taxa_clade_map,this_map_df))

	    ## add confirmed clade to title
	    this_title <- sprintf(" C%d: n_tax: %d canon_cost: %.4f\n", j,this_out_data$n_tax,this_out_data$canon_cost)
	    title_str = paste0(title_str,this_title)
        }
    }

    ## cat('taxa_clade_map')
    ## print(taxa_clade_map)

    ## get tree file for plotting
    tree_file <- paste0(tree_dir,'/',pthr_id,'.nwk')
    tree <- read.newick(tree_file)

    tree_tax_labels = do.call(rbind,strsplit(tree$tip.label,'|',fixed=TRUE))
    ## print(tree_tax_labels)
    tx_name_len = length(unique(tree_tax_labels[,3]))

    tree_label_str = sprintf("%s : %d taxa %d organisms (%s)\n",pthr_id,length(tree$tip.label),tx_name_len,tstamp)
    cat(tree_label_str)
    cat(title_str)

    ## if changes, then 'mltm', if ALL confirmed, 'mltmC'
    pdf_suff='mltm'
    if (! change_flag) {
       pdf_suff='mltmC'
    }

    ## get label_file of sequence lengths, canon/prop_canon status

    label_file <- paste0(label_dir,'/',pthr_id,label_suff)
    label_file0 <- paste0(label_dir,'/',pthr_id,label_suff0) 

    if (file.exists(label_file)) {
        tax_labels = read.table(label_file, header=TRUE,sep='\t',quote='',stringsAsFactors=TRUE)
        ## print(tax_labels)
    } else if (file.exists(label_file0)) {
        ## print(sprintf("No %s label - trying %s",label_file, label_file0))
        tax_labels = read.table(label_file0, header=TRUE,sep='\t',quote='',stringsAsFactors=TRUE)
	## pdf_suff = 'ltm'
    } else {
        print(sprintf("No %s label",label_file0))
        next
    }

    ## associate the suggested change index and the confirmed index with taxa
    ## one way would be to do a merge/join between the labels from the unique clades and the label_file

    ## extract acc from tax_labels
    labels <- strsplit(as.character(tax_labels$label),'|',fixed=TRUE)
    tax_label_df = do.call(rbind,labels)
    colnames(tax_label_df) = c('db','acc','tx_name')

    ## print(head(tax_label_df))

    ## add db, acc, tx_name back to tax_labels
    tax_labels <- cbind(tax_labels,tax_label_df)

    ## merge (left-join) tax_labels with taxa_clade_map (clades with suggestions) including all tax_labels (all.x=TRUE)
    tax_labels_cl <- merge(tax_labels, taxa_clade_map, by='acc',all.x=TRUE)

    ## get gc_groupids for acc's for tax_labels
    ## first check if it's in the .lab_lt file?

    if ('groupid' %in% colnames(tax_labels)) {
        taxa_clade_map <- merge(tax_labels[,c('groupid','acc')],taxa_clade_map,by='acc')
	## cat(sprintf(" %i  %s using 'groupid' in tax_labels\n",i,pthr_id),file=stderr())
	for (gc_id in unique(taxa_clade_map$groupid)) {
	    this_gc_set = taxa_clade_map[taxa_clade_map$groupid== gc_id,]
       	    this_label = unique(this_gc_set$label)
	    if (length(this_label) > 1) {
	       cat("*** Warning *** %s more than one gc_groupid\n",pthr_id)
	       cat("  ",this_gc_set$acc,"\n")
	       cat("  ",this_gc_set$label,"\n")
            }	       
	    else if (length(this_label) == 0) {
	       cat("*** Warning *** %s zero gc_groupid\n",pthr_id)
	       cat("  ",this_gc_set$acc,"\n")
	       cat("  ",this_gc_set$label,"\n")
            } else {
	       tax_labels_cl[tax_labels_cl$groupid == gc_id,]$label.y=this_label[1]
	    }	       
	}
    }
    ##
    ## if not, do we have a gene_centric file?
    ##
    else if (nrow(gc_df) > 0) {

        tax_labels_cl <- merge(gc_df[,c('groupid','acc')],tax_labels_cl,by='acc',all.y=TRUE)

	if (sum(is.na(tax_labels_cl$groupid))>0) {
	    tax_labels_cl[is.na(tax_labels_cl$groupid),]$groupid=0
	}

        taxa_clade_map <- merge(gc_df[,c('groupid','acc')],taxa_clade_map,by='acc')
	## print(taxa_clade_map)

	for (gc_id in unique(taxa_clade_map$groupid)) {
	    this_gc_set = taxa_clade_map[taxa_clade_map$groupid== gc_id,]
       	    this_label = unique(this_gc_set$label)
	    if (length(this_label) > 1) {
	       cat("*** Warning *** %s more than one gc_groupid\n",pthr_id)
	       cat("  ",this_gc_set$acc,"\n")
	       cat("  ",this_gc_set$label,"\n")
	    } else if (length(this_label) == 0) {
	       cat("*** Warning *** %s zero gc_groupid\n",pthr_id)
	       cat("  ",this_gc_set$acc,"\n")
	       cat("  ",this_gc_set$label,"\n")
            } else {
	       tax_labels_cl[tax_labels_cl$groupid == gc_id,]$label.y=this_label[1]
	       ## cat("gc_id: ",gc_id,'\n')
	       ## these_data = tax_labels_cl[tax_labels_cl$groupid == gc_id,]
	       ## print(these_data[,c('groupid','acc','label.y')])
	    }	       
	}
    }

    ## convert NA label.y to '0'
    if (sum(is.na(tax_labels_cl$label.y) > 0)) {
        tax_labels_cl[is.na(tax_labels_cl$label.y),]$label.y='0'
    }

    ## print(tax_labels[,c('acc','label')])
    ## print(tax_labels_cl[,c('acc','label.y')])

    ## cat('tax_labels not labels_cl\n')
    ## print(setdiff(tax_labels[,'acc'],tax_labels_cl[,'acc']))
    ## cat('tax_labels_cl not labels\n')
    ## print(setdiff(tax_labels_cl[,'acc'],tax_labels[,'acc']))

    ## cat("uniq len(tax_labels.acc):",length(unique(tax_labels$acc)),"\n")
    ## cat("uniq len(tax_label_cl.acc):",length(unique(tax_labels_cl$acc)),"\n")

    dup_label_cl = tax_labels_cl[duplicated(tax_labels_cl$acc),]
    if (nrow(dup_label_cl) > 0) {
        cat(sprintf("*** dup: %s\n",pthr_id))
	print(dup_label_cl)
	tax_labels_cl <- unique(tax_labels_cl)
    }

    ## build label string with taxon_string, length, clade_label (label.y)
    ## cat("len(tax_labels.label):",length(tax_labels$label),"\n")
    ## cat("len(tax_label_cl.y):",length(tax_labels_cl$label.y),"\n")

    if (length(tax_labels$label) != length(tax_labels_cl$label.y)) {
       print(paste("mismatched length: ",pthr_id))
       cat("len(tax_labels.label):",length(tax_labels$label),"\n")
       cat("len(tax_label_cl.y):",length(tax_labels_cl$label.y),"\n")
       print(tax_labels)
       print("\n")
       print(tax_labels_cl)
    }

    tax_labels$label_str = sprintf("%s [%d:%s]",tax_labels$label,tax_labels$len,tax_labels_cl$label.y)

    tree_fig=paste0(pdf_dir,'/',pthr_id,'_',pdf_suff,'.pdf')
    pdf(file=tree_fig, width=9, height = 0.3 * length(tree$tip.label))

    ## title string has list of clade scores, add PTHR_id and total taxa before
    title_str <- paste0(tree_label_str,title_str)

    print(read_tree_plot(pthr_id, tree, tax_labels, aln_dir, tree_dir, pdf_suff, title_str ))

    dev.off()

    cat(sprintf("%d %s printed\n",i,tree_fig))
}
