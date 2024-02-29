#!/usr/bin/env bash
#  22-June-2023
## modified produce one tree with multi_draw_nwk_tree_label5.R
## modified to capture options
##

## mtl_pdf_builer.sh qfomam tstamp pdf_data

BINDIR=$(pwd)

if [ "X$1" == "X" ]; then
    >&2 echo "You need to provide a set name, e.g. qfomam" 
    exit 22
fi
set=$1
if [ ! -d "${set}/" ]; then
  >&2 echo "No such set dir"
  exit 2
fi

if [ "X$2" == "X" ]; then
    >&2 echo "You need to provide a time_stamp, e.g. 230609" 
    exit 22
fi
tstamp=$2

if [ "X$3" == "X" ]; then
    pdf_data='pdf_data'
else
    pdf_data=$3
fi

if [ "X$4" == "X" ]; then
    tree_data='tree_data'
else
    tree_data=$4
fi

if [ "X$5" == "X" ]; then
    aln_data='aln_data'
else
    aln_data=$5
fi

if [ "X$6" == "X" ]; then
    lab_data='lab_data'
else
    lab_data=$6
fi

if [ "X$7" == "X" ]; then
    lab_suff='.lab_ltm'
else
    lab_suff=$7
fi

set -u
set -e
cd ${set}/
Rscript --vanilla $BINDIR/pdfcreation/multi_draw_nwk_label5gclt.R ${tstamp} ${pdf_data} ${tree_data} ${aln_data} ${lab_data} ${lab_suff} >${set}_multi_draw_gclt.log 2>${set}_multi_draw_gclt.err
