#!/usr/bin/env python
# coding: utf-8
"""
module providing functions related to output files
"""
from glob import glob
from o2t_utils import eprint, print_subfiles, delete_files


def output_headers(config=None):
    """
    prepare output headers with
    treebuilding parameters, cladescanning parameters,
    scoring function weights and fields names
    """
    wgts = config["suggestion_ranking_weights"]
    wgt_keys = list(wgts.keys())
    params_header = "#[trees] panther:{} up:{} aln_min_taxa:{} clade_min_taxa:{} tree_max_cost:{} tree_drop_fact:{}\n".format(
        config["panther_version"],
        config["up_release"],
        config["min_taxa_threshold"],
        config["taxa_threshold"],
        config["tree_max_cost"],
        config["tree_drop_fact"],
    )
    params_header += "#[suggestions] cost_diff:{} min_taxa:{} min_canon:{} max_cost:{} only_zero_cost:{} min_delta:{} min_delta_sp:{}\n".format(
        config["suggestion_score_difference"],
        config["suggestion_taxa_threshold"],
        config["suggestion_min_canon"],
        config["suggestion_max_clade_cost"],
        config["suggestion_only_zero_cost"],
        config["min_delta_threshold"],
        config["min_delta_sp_threshold"],
    )
    wgts_str = ["%s:%.1f" % (k, wgts[k]) for k in wgt_keys]
    params_header += "#[scoring] " + " ".join(wgts_str) + "\n"

    out_fields = [
        "pthr_id",
        "taxon",
        "canon_acc",
        "canon_len",
        "prop_acc",
        "prop_len",
        "rank_score",
        "canon_cost",
        "prop_cost",
    ]
    suggestions_fields = "\t".join(
        out_fields + wgt_keys + ["clade", "clade_members", "MANEstatus"]
    )
    confirmed_fields = "\t".join(
        [
            "pthr_id",
            "canon_cost",
            "n_sp",
            "n_tax",
            "clade",
            "clade_members",
            "MANEstatus",
        ]
    )
    gc_fields = (
        suggestions_fields.replace("canon_len\t", "canon_len\tcanon_md5\t").replace(
            "prop_len\t", "prop_len\tprop_md5\t"
        )
        + "\trelease\tdataset"
    )

    headers = {}
    headers["confirm"] = "#Canonical clades confirmed\n{}# {}\n".format(
        params_header, confirmed_fields
    )
    headers["gc"] = "#Cumulative Changes Proposed\n{}# {}\n".format(
        params_header, gc_fields
    )
    headers["skipped"] = "#Groups skipped\n{}# {}\n".format(
        params_header, "pthr_id\tclade\treason\tinfo"
    )
    headers["changes"] = "#Changes Proposed\n{}# {}\n".format(
        params_header, suggestions_fields
    )
    headers["conflict"] = "#oldcanon\treplacement\n"

    return headers


def combine_and_print_output(filename, key, header):
    """
    combine the temporary output files produced by each thread to create final output files
    """
    tempfiles = glob("{}..*".format(filename))
    if tempfiles:
        # combine the chunks
        eprint(
            "* Combining {} chunks into final output files for {}".format(
                len(tempfiles), key
            )
        )
        with open(filename, "w", encoding="utf-8") as fh:
            fh.write(header)
            print_subfiles(tempfiles, fh)
        delete_files(tempfiles)  # cleanup
    else:
        eprint("NOTICE: no data for {}".format(key))


def clean_up_tempfiles(filename_prefix):
    """
    delete temporary output files (produced by each thread)
    """
    delete_files(glob("{}..*".format(filename_prefix)))
