#!/usr/bin/env python
# coding: utf-8
"""
module providing functions
* build_tree_for_orthogroup
* get_three_for_orthogroup
"""
import os
import re
from multiprocessing import current_process
from Bio import Phylo
from Bio import AlignIO
from config_muscle import ALIGN_FORMAT
from o2t_utils import eprint, get_orthologs_df_from_pantherid, sequences2fasta
from o2t_alignment import get_alignment
from o2t_alignment2tree import alignment2tree
from o2t_tree2ndata import tree_to_ndata, get_leaf_acc

REFASTAHEADER_ADDISO = re.compile(r"^>(sp|tr)(|.*)\n")  # to add iso after sp/tr


def get_tree_for_orthogroup(
    orthoid,
    orthologs=None,
    ortho_df=None,
    config=None,
    get_sequences_fn=None,
    verbose=False,
):
    """
    used by scan_ndata and build_tree_for_orthogroup
    """
    nwk_tree = None
    tree_fname = os.path.join(config["tree_data_dir"], orthoid + ".nwk")

    # try to read previously created tree, if there is
    if config["cache_trees_flag"]:
        if os.path.isfile(tree_fname):
            # eprint("Reading cached tree for {}".format(orthoid)) #debug
            nwk_tree = Phylo.read(tree_fname, "newick")

    if nwk_tree is None:
        align = None
        aln_fname = os.path.join(config["aln_data_dir"], orthoid + ".aln")

        # try to read previously cached alignment, if there is
        if config["cache_alignments_flag"]:
            if os.path.isfile(aln_fname):
                # eprint("Reading cached alignment for {}".format(orthoid)) #debug
                align = AlignIO.read(aln_fname, ALIGN_FORMAT)  # muscle5.1

        # if we don't have alignment, let's gather sequences and align them now
        if align is None:
            # eprint("Creating new alignment for {}".format(orthoid)) #debug

            # get sequences
            if orthologs is None:  # can override list
                orthologs_df = get_orthologs_df_from_pantherid(
                    orthoid, ortho_df, remove_outliers=config["remove_outliers"]
                )
                orthologs = orthologs_df.index.to_list()
            else:
                orthologs_df = ortho_df.loc[orthologs]
            if not orthologs:
                eprint("No orthologs!")
                return [None] * 3
            # eprint("{} collecting sequences for {}".format(workerid, orthologs))
            sequences_dict = get_sequences_fn(orthologs, orthoid=orthoid)
            # eprint("{} got sequences for {}".format(workerid, set(sequences_dict.keys())))
            if not sequences_dict:
                eprint("No sequences!")
                return [None] * 3

            if verbose:
                outliers = ortho_df[
                    (ortho_df["pantherid"] == orthoid) & (ortho_df["outlier"])
                ].index.to_list()
                if outliers:
                    eprint(
                        "{}: INFO removed {} outliers: {}".format(
                            orthoid, len(outliers), outliers
                        )
                    )

            # add _iso to acc in fasta header if not canonical
            canonical_status = orthologs_df["is_canonical"].to_dict()  # dict acc -> T/F
            for acc in sequences_dict:
                if not canonical_status[acc]:  # then add iso
                    sequences_dict[acc] = REFASTAHEADER_ADDISO.sub(
                        r">\1_iso\2\n", sequences_dict[acc]
                    )

            if len(sequences_dict) < config["min_taxa_threshold"]:
                # heuristic: number of different organisms cannot be higher than number of sequences
                eprint(
                    "{} WARN: too few sequences remaining, no point aligning".format(
                        orthoid
                    )
                )
                return [None] * 3

            # check number of different taxa: if under threshold, skip
            organisms_dict = orthologs_df["org"].to_dict()  # dict acc -> oscode
            organisms = set(organisms_dict.values())
            if len(organisms) < config["min_taxa_threshold"]:
                eprint(
                    "{} WARN: not enough different taxa ({}) in sequences, no point continuing".format(
                        orthoid, len(organisms)
                    )
                )
                return [None] * 3

            # seq2fasta
            fasta_txt = sequences2fasta(sequences_dict, organisms_dict)
            if fasta_txt == "":
                eprint("{} ERROR: No fasta!".format(orthoid))
                return [None] * 3
            # eprint("prepared fasta: {}".format(fasta_txt)) #debug

            # call muscle
            align = get_alignment(fasta_txt)  # , textoutput=True)
            if align is None:
                eprint("{} ERROR: No alignment!".format(orthoid))
                return [None] * 3

            if config["cache_alignments_flag"]:  # if we use cache
                # write alignment file
                with open(aln_fname, "w", encoding="utf-8") as afd:
                    AlignIO.write(align, afd, ALIGN_FORMAT)

        # eprint(align) #debug short
        # AlignIO.write(align, sys.stderr, ALIGN_FORMAT) #debug full

        # build a gap-distance based evolutionary tree aka calc_pdist_g2.py
        nwk_tree = alignment2tree(align, dist_only=False, verbose=False)
        if nwk_tree is None:
            eprint("{} ERROR: No tree!".format(orthoid))
            return [None] * 3

        if config["cache_trees_flag"]:  # if we use a cache
            # write tree file
            with open(tree_fname, "w", encoding="utf-8") as tfd:
                Phylo.write(nwk_tree, tfd, "newick")

    return nwk_tree


def build_tree_for_orthogroup(
    orthoid,
    orthologs=None,
    ortho_df=None,
    config=None,
    verbose=False,
    debuginfo=False,
    draw_tree=False,
    get_sequences_fn=None,
):
    """
    build a phylo tree for a given orthoid (or for a given list of sequence accessions)
    """

    def _get_seqlen_from_acc(acc):
        """
        gets stored sequence length from dataframe or retrieve sequence (and then compute its length)
        return -1 if cannot get
        """
        seqlen = -1
        if acc in ortho_df.index:
            seqlen = ortho_df.loc[acc]["seqlen"]
        else:
            sequences = get_sequences_fn([acc])
            if sequences:
                seqlen = len(sequences[acc])
        return seqlen

    if current_process().name == "MainProcess":
        workerid = ""
    else:  # for debug of multithreading
        workerid = "[#{}]".format(current_process().name.split("-")[1])
    if verbose:
        eprint("{} building tree for {}".format(workerid, orthoid))

    nwk_tree = get_tree_for_orthogroup(
        orthoid,
        orthologs=orthologs,
        ortho_df=ortho_df,
        config=config,
        get_sequences_fn=get_sequences_fn,
        verbose=verbose,
    )

    # here we have nwk_tree
    # eprint(nwk_tree.format('newick'))
    if draw_tree and debuginfo:
        Phylo.draw(nwk_tree)  # debug, original tree

    # report zero-cost clades aka tree_cost2r.py
    l_tree, n_data, labels = tree_to_ndata(
        nwk_tree,
        ortho_df=ortho_df,
        target_id=orthoid,
        filename_prefix=orthoid,
        config=config,
        verbose=verbose,
        debuginfo=debuginfo,
        drop_taxa=True,
        create_ndata_file=True,
    )
    if debuginfo:
        eprint("labelled tree: {}".format(l_tree.format("newick")))  # debug
        eprint("n_data: {}".format(n_data))
        eprint("labels/ranks: {}".format(labels))

    if draw_tree:
        # l_tree_ori = copy.deepcopy(l_tree) #if n_data tree also needed
        # append sequence length information
        for cl in l_tree.find_clades(terminal=True):
            seqlen = _get_seqlen_from_acc(get_leaf_acc(cl.name))
            cl.name = cl.name + " [{}]".format(seqlen)
            # eprint(repr(cl), cl.name)
        Phylo.draw(l_tree)

    return l_tree, n_data, labels
