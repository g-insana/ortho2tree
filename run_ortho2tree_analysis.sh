#!/bin/bash
if [ "X$1" == "X" ]; then
  >&2 echo "You need to provide as first argument the name of the set, e.g. 5taxa"
  exit 22
fi
set=$1
set -e
set -u
outstamp=$(date +%y%m%d)
mkdir -p ${set}

python ortho2tree.py -set ${set} -outstamp ${outstamp} 2>${set}_log${outstamp}

# for pdf creation:
# pdfcreation/mtl_pdf_builder.sh ${set} ${outstamp}
# NOTE: requires R environment with appropriate packages, see the .R file under pdfcreation
