#!/usr/bin/env python
# coding: utf-8
"""
module providing the function scan_ndata_file
"""
import re
import os
import sys
import math
import time
from Bio import AlignIO
from Bio.Seq import Seq  # for aln2fastax
from .o2t_utils import eprint, get_orthologs_df_from_pantherid, elapsed_time
from .o2t_gc_integration import check_altgroup_suggestion, check_suggestion_against_prev
from .config_muscle import ALIGN_FORMAT
from .o2t_tree2ndata import selected_clade_cost
from .o2t_buildtree import get_tree_for_orthogroup


# helper functions for scan_ndata_file
def aln2fastax(orthoid, aln_data_dir):
    """
    to create faX files used in pdf creation
    """
    REALL2K = re.compile(r"[A-Za-z]")  # to convert all aa to K
    aln_filename = os.path.join(aln_data_dir, orthoid + ".aln")
    fax_filename = os.path.join(aln_data_dir, orthoid + ".faX")
    if not os.path.exists(aln_filename):
        eprint(
            "{}  ERROR: no alignment file found. Has it been deleted or not cached? Cannot create faX".format(
                orthoid
            )
        )
        return
    with open(aln_filename, "r", encoding="utf-8") as aln_fh:
        aln = AlignIO.read(aln_fh, ALIGN_FORMAT)
    for record in aln:
        this_seq = REALL2K.sub("K", str(record.seq))
        record.seq = Seq(this_seq)
    with open(fax_filename, "w", encoding="utf-8") as ffd:
        AlignIO.write(aln, ffd, format="fasta")


def parse_nd3_line(line, nd3_fields):
    data_fields = re.split(r" [:\!-][:=] ", line)
    n_info = dict(zip(nd3_fields, data_fields))

    (prop_cost, prop_ntax) = n_info["p_cost_n"].split("/")
    (xc_canon_cost, xc_canon_ntax) = n_info["xc_cost_n"].split("/")
    # (pr_canon_cost, pr_canon_ntax) = n_info['pr_cost_n'].split('/')
    # (gc_canon_cost, gc_canon_ntax) = n_info['gc_cost_n'].split('/')

    prop_acc_info = {}

    prop_iso = {}
    prop_sp_iso = []
    prop_tr_iso = []
    prop_ref_iso = []

    for acc_cost in n_info["p_acc_info"].split(","):
        split1 = acc_cost.split(":")
        split2 = split1[1].split(" ")
        split3 = split2[0].split("|")
        taxa = split1[0]
        db = split3[0]
        acc = split3[1]
        prop_acc_info[taxa] = {"db": db, "acc": acc, "cost": float(split2[1])}
        if db == "sp_iso":
            prop_iso[taxa] = {"p_db": db, "p_acc": acc}
            prop_sp_iso.append(taxa)
        elif db == "tr_iso":
            prop_iso[taxa] = {"p_db": db, "p_acc": acc}
            prop_tr_iso.append(taxa)
        elif db == "ref_iso":
            prop_iso[taxa] = {"p_db": db, "p_acc": acc}
            prop_ref_iso.append(taxa)

    canon_accs = {}
    if n_info["xc_acc_info"]:
        for info in n_info["xc_acc_info"].split(","):
            (c_taxa, c_acc) = info.split(":")
            canon_accs[c_taxa] = c_acc
            if c_taxa in prop_iso:
                prop_iso[c_taxa]["canon_acc"] = c_acc
    else:
        eprint("ERROR! empty xc_acc_info in {}!!!".format(n_info))
        return {"clade": n_info["clade"]}

    return {
        "clade": n_info["clade"],
        "prop_cost": float(prop_cost),
        "prop_ntax": int(prop_ntax),
        "canon_cost": float(xc_canon_cost),
        "canon_ntax": int(xc_canon_ntax),
        "prop_acc_info": prop_acc_info,
        "canon_accs": canon_accs,
        "prop_iso": prop_iso,
        "prop_sp_iso": prop_sp_iso,
        "prop_tr_iso": prop_tr_iso,
        "prop_ref_iso": prop_ref_iso,
    }


def calc_rank(wts, vals):
    """
    calc_rank generates a score using a dictionary of weights and
    values.  Each value is multiplied by its corresponding weight, and
    the sum is returned.
    """

    tot_score = 0.0

    for wt in wts.keys():
        # if verbose: eprint("DEBUG: {}: {} * {}".format(wt, wts[wt], vals[wt])) #DEBUG
        tot_score += wts[wt] * vals[wt]

    return tot_score


def scale_exp(value, lowest):
    if value <= lowest:
        return 1.0
    scaled = math.exp(-float(value))
    return min(scaled, 1.0)


def format_confirmed(orthologs_df, orthoid, prop_clade):
    """
    prepare output line for confirmed canonical solution
    """
    clade_taxa = prop_clade["prop_acc_info"].keys()
    # taxa and accessions from the clade to which suggestion belongs to
    clade_members = ",".join(
        [x + ":" + prop_clade["prop_acc_info"][x]["acc"] for x in clade_taxa]
    )
    confirmed_output = "{}\t{}\t{}\t{}\t{}\t{}".format(
        orthoid,  # orthoid
        # canon_cost (cost of canons)
        prop_clade["canon_cost"],
        # n_sp (number of sp entries)
        prop_clade["n_sp"],
        # n_tax (number of taxa)
        prop_clade["n_tax"],
        # clade (identifier from n_data file)
        prop_clade["clade"],
        clade_members,
    )  # members of the clade
    # MANEtagging of confirmed
    if "HUMAN" in clade_taxa:
        if orthologs_df.loc[prop_clade["prop_acc_info"]["HUMAN"]["acc"]]["has_mane"]:
            confirmed_output += "\tMANE_good"
        else:
            confirmed_output += "\tMANE_bad"
    else:
        if orthologs_df["has_mane"].sum() > 0:
            confirmed_output += (
                "\tMANE_bad"  # MANE tags something which is not in our clade
            )
        else:
            confirmed_output += "\tNAM"

    return confirmed_output


def format_changed(
    orthologs_df, orthoid, prop_clade, canon_acc, p_acc, suggestion_ranking_weights
):
    """
    prepare output line for proposed replacement
    """
    # currently consisting in 19 columns:

    # first 4 columns: orthoid, taxon, canon db+acc and length

    suggestion_output = "{}\t{}\t{}\t{}".format(
        orthoid,
        # taxon
        orthologs_df.loc[canon_acc].org,
        "{}|{}".format(
            # canon_acc
            orthologs_df.loc[canon_acc].entry_type,
            canon_acc,
        ),
        # canon_len
        orthologs_df.loc[canon_acc].seqlen,
    )
    # variant for gc_output, adding also md5sum information
    suggestion_output_md5 = suggestion_output + "\t{}".format(
        orthologs_df.loc[canon_acc].md5sum
    )  # canon_md5

    # other 2 columns: suggested replacement db+acc and length
    prop_acc_info = "\t{}\t{}".format(
        "{}|{}".format(orthologs_df.loc[p_acc].entry_type, p_acc),  # prop_acc
        # prop_len
        orthologs_df.loc[p_acc].seqlen,
    )
    suggestion_output += prop_acc_info

    # variant for gc_output
    suggestion_output_md5 += prop_acc_info + "\t{}".format(
        orthologs_df.loc[p_acc].md5sum
    )  # prop_md5

    # 3 more columns: rank_score, canon_cost, prop_cost
    score_and_costs_info = "\t{:.2f}\t{:.4f}\t{:.4f}\t".format(
        prop_clade["prop_iso_score"],  # rank_score
        # canon_cost
        prop_clade["canon_cost"],
        # prop_cost
        prop_clade["prop_cost"],
    )
    suggestion_output += score_and_costs_info
    suggestion_output_md5 += score_and_costs_info

    # 7 weight values used in score computation: n_sp n_tax n_canon wn_canon scaled_prop_f scaled_p_cost p_len_diff
    wt_strs = []
    for k in suggestion_ranking_weights.keys():
        if isinstance(prop_clade[k], float):
            wt_strs.append("%.4f" % (prop_clade[k]))
        else:
            wt_strs.append(str(prop_clade[k]))
    suggestion_output += "\t".join(wt_strs)
    suggestion_output_md5 += "\t".join(wt_strs)

    # 2 more columns: clade identifier (from n_data file) and members
    # taxa and accessions from the clade to which suggestion belongs to:
    clade_taxa = prop_clade["prop_acc_info"].keys()
    clade_members = ",".join(
        [x + ":" + prop_clade["prop_acc_info"][x]["acc"] for x in clade_taxa]
    )
    clade_info = "\t{}\t{}".format(
        prop_clade["clade"], clade_members  # clade
    )  # clade_members
    suggestion_output += clade_info
    suggestion_output_md5 += clade_info

    # 1 final column: MANEtagging of suggestions
    mane_marked = False
    if orthologs_df.loc[canon_acc].org == "HUMAN":
        if orthologs_df.loc[canon_acc]["has_mane"]:
            # we are suggesting to replace a MANE canonical
            suggestion_output += "\tMANE_bad"
            suggestion_output_md5 += "\tMANE_bad"
            mane_marked = True
        elif orthologs_df.loc[p_acc]["has_mane"]:
            # we are suggesting a replacement which has MANE tag
            suggestion_output += "\tMANE_good"
            suggestion_output_md5 += "\tMANE_good"
            mane_marked = True
    if not mane_marked:
        if (
            "HUMAN" in clade_taxa
            and orthologs_df.loc[prop_clade["prop_acc_info"]["HUMAN"]["acc"]][
                "has_mane"
            ]
        ):
            suggestion_output += "\tMANE_good"
            suggestion_output_md5 += "\tMANE_good"
        else:
            if orthologs_df["has_mane"].sum() > 0:
                # MANE tags something which is not in our clade
                suggestion_output += "\tMANE_bad"
                suggestion_output_md5 += "\tMANE_bad"
            else:
                suggestion_output += "\tNAM"
                suggestion_output_md5 += "\tNAM"

    return suggestion_output, suggestion_output_md5


def tset_suggestions(
    prop_clade, suggestions_by_acc, suggestions_by_groupid, orthologs_df, do_set=False
):
    """
    test and set suggestions to avoid duplicates
    this is in order not to re-use something that from a previous clade has been selected
    """
    have_acc = where_seen_clade = False
    prop_acc_info = prop_clade["prop_acc_info"]
    prop_iso = prop_clade["prop_iso"]
    # check each taxa
    for p_taxa in prop_acc_info.keys():
        prop_canon_acc = prop_clade["canon_accs"][p_taxa]
        this_groupid = orthologs_df.loc[prop_canon_acc].groupid
        acc = prop_acc_info[p_taxa]["acc"]

        if acc in suggestions_by_acc:
            where_seen_clade = "by acc: {}".format(acc)
            have_acc = True
        elif this_groupid in suggestions_by_groupid:
            where_seen_clade = "by groupid: {}={}".format(acc, this_groupid)
            have_acc = True

        # update suggestion IF do_set=True or canonical
        if do_set or (p_taxa not in prop_iso):
            suggestions_by_acc[prop_canon_acc] = p_taxa
            suggestions_by_groupid[this_groupid] = p_taxa

    return have_acc, where_seen_clade


def scan_ndata_file(
    orthoid,
    ortho_df=None,
    prevgc_df=None,
    config=None,
    get_sequences_fn=None,
    verbose=False,
):
    """
    identify best clades from n_data files to propose canonical replacement
    suggestions or confirm clade canonicals
    """
    suggestion_ranking_weights = config["suggestion_ranking_weights"]
    process_start_time = time.time()

    gc_output_text = ""  # all cumulative gc
    suggestions_output_text = ""  # let's print all suggestions from one file together
    confirmed_output_text = (
        ""  # in order to create output with confirmed canonical clades
    )
    skipped_output_text = ""  # in order to create output for skipped groups
    conflict_output_text = ""  # in order to create output for conflicts

    nd3_fields = (
        "clade",
        "p_cost_n",
        "p_acc_info",
        "xc_cost_n",
        "xc_acc_info",
        "pc_cost_n",
        "pc_acc_info",
        "gc_cost_n",
        "gc_acc_info",
    )

    if orthoid is None:
        eprint("ERROR! no orthoid given to scan_ndata_file()")
        sys.exit(22)

    ndata_fname = os.path.join(config["n_data_dir"], orthoid + ".n_data2")
    if not os.path.isfile(ndata_fname):
        # GCTODO: maybe we have previous gc cumulative for that orthoid? if so, first '' should change
        return [""] * 5

    # for sequence lengths
    orthologs_df = get_orthologs_df_from_pantherid(
        orthoid, ortho_df, remove_outliers=True
    )[
        [
            "org",
            "groupid",
            "seqlen",
            "md5sum",
            "entry_type",
            "is_canonical",
            "is_fragment",
            "has_mane",
        ]
    ]
    orthologs_dict = orthologs_df.to_dict(orient="index")

    n_clades = []
    with open(ndata_fname, "r", encoding="utf-8") as fd:
        for line in fd:
            line = line.strip("\n")
            if line.startswith("#"):
                continue
            if line.startswith(">"):
                # no need, we have orthoid already
                # pthr_id = line[1:]
                continue
            else:
                n_clades.append(parse_nd3_line(line, nd3_fields))

    better_clades = []
    prot_lengths = {}

    for ix, n_clade in enumerate(n_clades):
        clade_name = n_clade["clade"]
        if "canon_cost" not in n_clade:
            # SKIPPED due to errors in clade data
            eprint("CASE 0: {}".format(orthoid))
            skipped_output_text += "{}\t{}\tC0\t-\n".format(orthoid, clade_name)
            continue
        delta = n_clade["canon_cost"] - n_clade["prop_cost"]
        (prop_ntax, canon_ntax) = (n_clade["prop_ntax"], n_clade["canon_ntax"])

        info_str = "%s\t%0.4f\t%0.4f\t%d\t%d" % (
            orthoid,
            delta,
            n_clade["prop_cost"],
            prop_ntax,
            canon_ntax,
        )
        if verbose:
            eprint("* checking clade {}: {}".format(clade_name, info_str))
        if delta < -0.0001 and prop_ntax != canon_ntax:
            # SKIPPED due to trembl canonical fragment
            if verbose:
                eprint("CASE 1: {}".format(orthoid))
            skipped_output_text += (
                "{}\t{}\tC1negdelta+orphaniso\tΔ:{:.4f}iso:{}\n".format(
                    orthoid, clade_name, delta, "TODO"
                )
            )
            if ix == 0:  # first one in the list
                if verbose:
                    eprint(
                        "**neg delta %.4f/n_tax prop_ntax: %d canon_ntax: %d**\t%s\n"
                        % (
                            delta,
                            prop_ntax,
                            canon_ntax,
                            info_str,
                        )
                    )
            if verbose:
                eprint(
                    "{}: delta < 0: suggestion is worse than existing canon and missing canonical".format(
                        clade_name
                    )
                )
        elif delta < -0.0001:  # bad delta, but matching prop_ntax/canon_ntax
            if verbose:
                eprint("CASE 2: {}".format(orthoid))
            skipped_output_text += "{}\t{}\tC2negdelta>confirmed?\tΔ:{:.4f}\n".format(
                orthoid, clade_name, delta
            )
            # CONFIRMED?
            if ix == 0:
                if verbose:
                    eprint(
                        "**neg delta %.4f **\t%s\n"
                        % (
                            delta,
                            info_str,
                        )
                    )
            if verbose:
                eprint(
                    "{}: delta < 0: suggestion is worse than existing canon".format(
                        clade_name
                    )
                )
        elif prop_ntax != canon_ntax:
            if verbose:
                eprint("CASE 3: {}".format(orthoid))
            # should we print the orphan isoform as third field? TODO
            skipped_output_text += "{}\t{}\tC3orphaniso\tiso:{}\n".format(
                orthoid, clade_name, "-"
            )
            if verbose:
                eprint(
                    "**n_tax prop_ntax: %d canon_ntax: %d**\t%s\n"
                    % (
                        prop_ntax,
                        canon_ntax,
                        info_str,
                    )
                )
        else:
            # delta is >= 0.0
            # this test should only be done after we check for a zero-cost canonical set
            if delta < config["suggestion_score_difference"]:  # delta is too small
                if verbose:
                    eprint("CASE 4: {}".format(orthoid))
                skipped_output_text += "{}\t{}\tC4delta_too_small\tΔ:{:.4f}\n".format(
                    orthoid, clade_name, delta
                )
                # CONFIRMED IF cost is low otherwise SKIPPED
                if verbose:
                    eprint(
                        "delta ({}) < suggestion score ({})".format(
                            delta, config["suggestion_score_difference"]
                        )
                    )  # debug
                continue

            # is zero-cost proposal required??
            if config["suggestion_only_zero_cost"] and n_clade["prop_cost"] > 0.00:
                if verbose:
                    eprint("CASE 5: {}".format(orthoid))
                skipped_output_text += (
                    "{}\t{}\tC5sugg_only_zero\tprop_cost:{:.2f}\n".format(
                        orthoid, clade_name, n_clade["prop_cost"]
                    )
                )
                if verbose:
                    eprint("n_clade prop_cost > 0")  # debug
                continue

            # is proposed cost low enough?
            if config["suggestion_max_clade_cost"] < n_clade["prop_cost"]:
                if verbose:
                    eprint("CASE 6: {}".format(orthoid))
                skipped_output_text += (
                    "{}\t{}\tC6sugg_too_costly\tprop_cost:{:.2f}\n".format(
                        orthoid, clade_name, n_clade["prop_cost"]
                    )
                )
                if verbose:
                    eprint(
                        "config['suggestion_max_clade_cost'] ({}) < n_clade prop_cost ({})".format(
                            # debug
                            config["suggestion_max_clade_cost"],
                            n_clade["prop_cost"],
                        )
                    )
                continue

            # is number of taxa high enough?
            if config["suggestion_taxa_threshold"] > n_clade["prop_ntax"]:
                if verbose:
                    eprint("CASE 7: {}".format(orthoid))
                skipped_output_text += (
                    "{}\t{}\tC7not_enough_taxa\tprop_ntax:{}\n".format(
                        orthoid, clade_name, n_clade["prop_ntax"]
                    )
                )
                if verbose:
                    eprint(
                        "config['suggestion_taxa_threshold'] < n_clade prop_ntax"
                    )  # debug
                continue

            # possible good candidate, get lengths if we don't have them
            if verbose:
                eprint("possible good candidate clade: {}".format(clade_name))
            # first for substitutes

            prop_accs = []
            prop_iso_accs = []
            prop_iso_len = []

            for taxa in n_clade["prop_acc_info"].keys():
                this_prot = n_clade["prop_acc_info"][taxa]
                this_prot["taxon"] = taxa
                prop_accs.append(this_prot)

                if this_prot["acc"] not in prot_lengths:
                    prot_lengths[this_prot["acc"]] = orthologs_dict[this_prot["acc"]][
                        "seqlen"
                    ]

                if re.search(r"iso", this_prot["db"]):
                    prop_iso_accs.append(this_prot)
                    prop_iso_len.append(prot_lengths[this_prot["acc"]])

            # count canonicals in clade
            (sp_cnt, can_cnt, w_can_cnt) = (0, 0, 0)
            for p_acc in prop_accs:
                if not re.search("iso", p_acc["db"]):
                    if p_acc["taxon"] in config["suggestion_taxon_weights"]:
                        w_can_cnt += config["suggestion_taxon_weights"][p_acc["taxon"]]
                    else:
                        w_can_cnt += config["suggestion_taxon_weight_default"]
                if p_acc["db"] == "sp":
                    sp_cnt += 1
                    can_cnt += 1
                elif p_acc["db"] == "tr":
                    can_cnt += 1

            if can_cnt < config["suggestion_min_canon"]:
                if verbose:
                    eprint("CASE 8: {}".format(orthoid))
                skipped_output_text += (
                    "{}\t{}\tC8low_cano_count\tcan_count:{}\n".format(
                        orthoid, clade_name, can_cnt
                    )
                )
                if verbose:
                    eprint(
                        "{}: canonical count {} lower than min_canon {}".format(
                            clade_name, can_cnt, config["suggestion_min_canon"]
                        )
                    )  # debug
                continue

            # do the same for canonicals
            canon_lengths = []
            for taxa in n_clade["canon_accs"].keys():
                this_acc = n_clade["canon_accs"][taxa]
                if this_acc not in prot_lengths:
                    prot_lengths[this_acc] = orthologs_dict[this_acc]["seqlen"]
                canon_lengths.append(prot_lengths[this_acc])

            # now do some scoring:
            # (a) # of swissprot
            # (b) # of canonicals
            # canon median length
            # prop length diff

            canon_lengths.sort()
            # if verbose: eprint("DEBUG: canon_lengths:", canon_lengths) #DEBUG
            canon_len = len(canon_lengths)
            canon_half = int(canon_len / 2)
            if canon_len % 2 == 1:
                med_canon_length = float(canon_lengths[canon_half])
            else:
                med_canon_length = (
                    float(canon_lengths[canon_half - 1] + canon_lengths[canon_half])
                    / 2.0
                )

            prop_len_diff = 0.0
            for p_len in prop_iso_len:
                prop_len_diff += abs(float(p_len) - med_canon_length)

            if len(prop_iso_len) > 0:
                prop_len_diff /= float(len(prop_iso_len))
            else:
                prop_len_diff = 0.000

            # prop_len_diff needs to be rescaled so that no difference is good, and a large difference is bad

            prop_len_diff_score = scale_exp(float(prop_len_diff), 0.0)
            scaled_p_cost = scale_exp(float(n_clade["prop_cost"]), 0.0)
            scaled_prop_f = scale_exp(
                float(len(prop_iso_accs)) / float(n_clade["prop_ntax"]), 0.0
            )

            # now we need a scoring function

            rank_vals = {
                "n_sp": sp_cnt,
                "n_canon": can_cnt,
                "wn_canon": w_can_cnt,
                "p_len_diff": prop_len_diff_score,
                "n_iso": len(n_clade["prop_iso"]),
                "n_clade": n_clade,
                "n_tax": n_clade["prop_ntax"],
                "scaled_prop_f": scaled_prop_f,
                "scaled_p_cost": scaled_p_cost,
            }

            for k in config["suggestion_ranking_weights"].keys():
                n_clade[k] = rank_vals[k]

            n_clade["prop_iso_score"] = calc_rank(
                config["suggestion_ranking_weights"], rank_vals
            )

            better_clades.append(n_clade)

    if len(better_clades) and verbose:
        # eprint("\nbetter clades to analyse: {}".format([x['clade'] for x in sorted(better_clades, key=lambda x: -x['prop_iso_score'])]))
        eprint("\nbetter clades to analyse:")
        out_list = (
            "clade",
            "n_tax",
            "wn_canon",
            "n_canon",
            "scaled_p_cost",
            "scaled_prop_f",
            "p_len_diff",
            "prop_iso_score",
            "prop_cost",
        )
        for x in sorted(better_clades, key=lambda x: -x["prop_iso_score"]):
            out_vals = [x[y] for y in out_list]
            # eprint("{} {} {:.1f} {} {:.4f} {:.4f} {:.4f} {:.2f} {:.4f}".format(x[y] for y in out_list))
            eprint(
                "{} {} {:.1f} {} {:.4f} {:.4f} {:.4f} {:.2f} {:.4f}".format(*out_vals)
            )

    # dict to avoid printing repeated suggestions: only print the best one
    suggestions_by_acc = {}
    # needs to before checking confirmed clades so that a confirmed acc is not re-used
    # should be by genecentric_id, rather than acc, so that all the isoforms are also excluded
    suggestions_by_groupid = {}

    suggestions = (
        {}
    )  # dict to avoid printing repeated suggestions: only print the best one
    suggestion_clade2rank = {}

    # get the original in_tree for later analyses, outside a loop
    this_tree = get_tree_for_orthogroup(
        orthoid,
        ortho_df=ortho_df,
        config=config,
        get_sequences_fn=get_sequences_fn,
        verbose=verbose,
    )

    # to hold clades which have been reviewed and for which either CONFIRM or CHANGE has been proposed
    used_clades = set()
    for prop_clade in sorted(better_clades, key=lambda x: -x["prop_iso_score"]):
        pc_info = prop_clade["prop_acc_info"]
        used_clades_key = ":".join(sorted([pc_info[x]["acc"] for x in pc_info.keys()]))

        # if prop_clade['clade'] in used_clades: #skip this clade as it has already been used
        if (
            used_clades_key in used_clades
        ):  # skip this clade as it has already been used
            if verbose:
                eprint(
                    "* skipping clade {}: {}/{} [already used]".format(
                        prop_clade["clade"],
                        prop_clade["prop_cost"],
                        prop_clade["prop_ntax"],
                    )
                )
            continue

        if verbose:
            eprint(
                "* analysing clade {} :: {}/{}".format(
                    prop_clade["clade"],
                    prop_clade["prop_cost"],
                    prop_clade["prop_ntax"],
                )
            )

        if prop_clade["n_tax"] == prop_clade["n_canon"]:
            if verbose:
                eprint("  * all members of clade are canonical, skipping")
            have_acc, where_seen_clade = tset_suggestions(
                prop_clade,
                suggestions_by_acc,
                suggestions_by_groupid,
                orthologs_df,
                True,
            )
            if not have_acc:
                if verbose:
                    eprint("  * printing as confirmed {}".format(prop_clade["clade"]))
                used_clades.add(used_clades_key)  # adding to used_clades
                confirmed_output_text += (
                    format_confirmed(orthologs_df, orthoid, prop_clade) + "\n"
                )
            continue
        if prop_clade["prop_cost"] >= prop_clade["canon_cost"]:
            if verbose:
                eprint(
                    "  * canon_cost: %.4f <= prop_cost: %.4f skipping"
                    % (prop_clade["canon_cost"], prop_clade["prop_cost"])
                )

            # replace prop_clade with canon_accs
            this_clade = pc_info
            for tax in prop_clade["prop_iso"].keys():
                this_clade[tax]["acc"] = prop_clade["canon_accs"][tax]
                this_clade[tax]["db"] = re.sub("_iso", "", this_clade[tax]["db"])

            have_acc, where_seen_clade = tset_suggestions(
                prop_clade,
                suggestions_by_acc,
                suggestions_by_groupid,
                orthologs_df,
                True,
            )
            if not have_acc:
                if verbose:
                    eprint("  * printing as confirmed {}".format(prop_clade["clade"]))
                used_clades.add(used_clades_key)  # adding to used_clades
                confirmed_output_text += (
                    format_confirmed(orthologs_df, orthoid, prop_clade) + "\n"
                )
            continue

        # changed clades new/old accs are in n_clade['prop_iso']

        # skip if this acc is already in a chosen canonical solution
        have_acc, where_seen_clade = tset_suggestions(
            prop_clade, suggestions_by_acc, suggestions_by_groupid, orthologs_df, False
        )
        if have_acc:
            if verbose:
                eprint("  * acc already in chosen canonical solution")
            continue

        prop_iso = prop_clade["prop_iso"]
        # if verbose: eprint("clade {} prop_acc_info keys {}".format(prop_clade['clade'], prop_acc_info.keys())) #debug
        # if verbose: eprint("* analysing clade {}, prop_iso: {}".format(prop_clade['clade'], prop_iso))
        # if verbose: eprint("* taxa: {}".format(prop_acc_info.keys())) #debug
        # check for proposed changes

        prop_hybrid = False

        suggest_conflict = False
        new_suggestions = []

        for p_taxa in pc_info.keys():
            # if (p_taxa in prop_iso and p_taxa not in prop_change): no need to check p_taxa in prop_change

            if p_taxa in prop_iso:  # prop_iso only has non-canonical taxa
                ptaxa_cost = pc_info[p_taxa]["cost"]
                canon_acc = prop_clade["canon_accs"][p_taxa]
                # get database of canonical (can't use prop_iso, need canon)
                canon_db = orthologs_df.loc[canon_acc].entry_type
                canon_groupid = orthologs_df.loc[canon_acc].groupid

                if (
                    canon_acc in suggestions_by_acc
                    or canon_groupid in suggestions_by_groupid
                ):
                    suggest_conflict = True
                    eprint("{} CONFLICTtrue for {}".format(orthoid, canon_acc))
                    continue

                new_suggestions.append(
                    {"acc": canon_acc, "groupid": canon_groupid, "taxa": p_taxa}
                )

                # is the change large enough that SwissProt does not matter?
                canon_is_sp = "sp" == canon_db
                iso_is_sp = prop_iso[p_taxa]["p_db"][:2] == "sp"
                iso_is_tr = prop_iso[p_taxa]["p_db"][:2] == "tr"

                if ptaxa_cost >= config["min_delta_sp_threshold"]:
                    # because delta_sp_threshold is always > delta_threshold we check this first
                    continue
                elif not canon_is_sp and ptaxa_cost >= config["min_delta_threshold"]:
                    continue
                else:  # keep canon
                    # keeping canonical means:
                    prop_hybrid = True  # record that we need a new clade cost

                    # replacing isoform in prop_acc_info
                    prop_clade["prop_acc_info"][p_taxa] = {
                        "db": canon_db,
                        "acc": canon_acc,
                        "cost": 0.0,
                    }

                    # removing from prop_iso
                    del prop_clade["prop_iso"][p_taxa]

                    # removing from prop_sp_iso
                    if iso_is_sp:
                        if p_taxa in prop_clade["prop_sp_iso"]:
                            prop_clade["prop_sp_iso"].remove(p_taxa)
                        else:
                            eprint(
                                " cannot remove {} from prop_clade['prop_sp_iso']".format(
                                    p_taxa
                                )
                            )
                    elif iso_is_tr:
                        if p_taxa in prop_clade["prop_tr_iso"]:
                            prop_clade["prop_tr_iso"].remove(p_taxa)
                        else:
                            eprint(
                                " cannot remove {} from prop_clade['prop_tr_iso']".format(
                                    p_taxa
                                )
                            )
                    else:
                        if p_taxa in prop_clade["prop_ref_iso"]:
                            prop_clade["prop_ref_iso"].remove(p_taxa)
                        else:
                            eprint(
                                " cannot remove {} from prop_clade['prop_ref_iso']".format(
                                    p_taxa
                                )
                            )

                    # update n_canon
                    prop_clade["n_canon"] += 1
                    # update n_sp
                    prop_clade["n_sp"] += 1

        # check to see if we have a suggest_conflict -- if not, then add to suggestions
        if not suggest_conflict:
            for suggest in new_suggestions:
                suggestions_by_acc[suggest["acc"]] = suggest["taxa"]
                suggestions_by_groupid[suggest["groupid"]] = suggest["taxa"]
            if verbose:
                eprint(
                    "{} CONFLICTfalse;adding sugg_by_acc: {} sugg_by_groupid: {}".format(
                        orthoid, suggestions_by_acc, suggestions_by_groupid
                    )
                )
        else:
            continue

        # if we have reverted from a prop_iso to a canon, then we need to recalculate the tree cost
        # but only if there are still proposed isoforms
        # we could avoid calculating the common_clade and cost if they are all canonical, since it has been done
        ##
        if prop_hybrid and prop_clade["n_tax"] > prop_clade["n_canon"]:
            # extract proposed taxa:
            this_acc_dict = {}
            this_prop_acc_info = prop_clade["prop_acc_info"]

            # this_acc_dict is used to build the list of acc's used
            # to find the clade and cost -- but it is already
            # available in prop_clade['prop_acc_info'] (this_prop_acc_info)

            this_acc_dict = {
                this_prop_acc_info[x]["acc"]: True for x in this_prop_acc_info.keys()
            }

            leaf_list = []
            n_taxa = len(list(this_acc_dict.keys()))
            for leaf in this_tree.get_terminals():
                (db, acc, id_oscode) = leaf.name.split("|")
                if acc in this_acc_dict:
                    leaf_list.append(leaf)
                if len(leaf_list) == n_taxa:
                    break

            common_clade = this_tree.common_ancestor(leaf_list)
            this_cost = selected_clade_cost(this_tree, common_clade, leaf_list)

            ##
            prop_clade["prop_cost"] = this_cost
            prop_clade["clade"] = common_clade.name

        # this section used to be above the previous for loop, but it
        # has been moved because the prop_clade['n_canon'] can be
        # increased in that loop

        if prop_clade["n_tax"] == prop_clade["n_canon"]:
            # eprint("sugg_by_acc: {}\nsugg_by_groupid: {}".format(suggestions_by_acc, suggestions_by_groupid)) #debug
            have_acc, where_seen_clade = tset_suggestions(
                prop_clade,
                suggestions_by_acc,
                suggestions_by_groupid,
                orthologs_df,
                True,
            )
            if verbose:
                eprint(
                    "  * all members of clade are (now) canonical (n_canon={}), skipping?".format(
                        prop_clade["n_canon"]
                    )
                )
            # if (not have_acc) and verbose: eprint("  * printing as confirmed {}".format(prop_clade['clade']))
            used_clades.add(used_clades_key)  # adding to used_clades
            confirmed_output_text += (
                format_confirmed(orthologs_df, orthoid, prop_clade) + "\n"
            )
            # elif verbose:
            # eprint("  * skipping to print as confirmed {} because have_acc is True as acc was seen in {}; sugg_by_acc: {}; sugg_by_groupid: {}".format(prop_clade['clade'], where_seen_clade, suggestions_by_acc, suggestions_by_groupid))
            continue

        # put suggestions in output
        for p_taxa in prop_clade["prop_acc_info"].keys():
            this_acc_info = prop_clade["prop_acc_info"][p_taxa]

            # only make suggestions if non-canonical
            if p_taxa not in prop_iso:
                if verbose:
                    eprint("{}, has canon: {}".format(p_taxa, this_acc_info["acc"]))
                continue

            # this has been superceded by suggestions_by_acc and suggestions_by_groupid
            # check to see if this canon has already been used up
            p_acc = this_acc_info["acc"]
            canon_acc = prop_clade["canon_accs"][p_taxa]
            if canon_acc in suggestions:
                # we already have a suggestion to replace this canon, let's check if it was for the same iso as replacement:
                if p_acc in suggestions[canon_acc]:
                    # this suggestion was printed already, skipping
                    if verbose:
                        eprint(
                            "{}, skipping suggesting {}->{} as we already printed this".format(
                                p_taxa, canon_acc, p_acc
                            )
                        )
                    continue  # do not put in output
                suggestions[canon_acc].add(p_acc)
            else:
                suggestions[canon_acc] = {p_acc}

            if prop_clade["clade"] in suggestion_clade2rank:
                suggestion_rank_level = suggestion_clade2rank[prop_clade["clade"]]
            else:
                suggestion_rank_level = len(suggestion_clade2rank) + 1
                suggestion_clade2rank[prop_clade["clade"]] = suggestion_rank_level

            if verbose:
                eprint(
                    "{}: {}, printing suggestion {}->{} at ranklevel {}".format(
                        orthoid, p_taxa, canon_acc, p_acc, suggestion_rank_level
                    )
                )
            used_clades.add(used_clades_key)  # adding to used_clades
            formatted_change, formatted_change_md5 = format_changed(
                orthologs_df,
                orthoid,
                prop_clade,
                canon_acc,
                p_acc,
                suggestion_ranking_weights,
            )
            suggestions_output_text += formatted_change + "\n"

            if not prevgc_df.empty: # if we have a prevgc_df
                # first check if a suggestion involving a different orthoid existed (different family assignment)
                conflict_text = check_altgroup_suggestion(
                    prevgc_df,
                    p_taxa,
                    canon_acc,
                    p_acc,
                    orthoid,
                    config=config,
                    verbose=True,
                )
                if conflict_text: #this one for obsolete ones belonging to different orthoid
                    # to be able to remove conflicts even in multi-thread operation (the marking of conflict will not be shared across parallel workers)
                    conflict_output_text += conflict_text
                # then check the new_sugg against prev_sugg for possible conflicts or flipflips:
                print_new_sugg, conflict_text = check_suggestion_against_prev(
                    prevgc_df,
                    p_taxa,
                    canon_acc,
                    p_acc,
                    orthoid,
                    config=config,
                    verbose=True,
                )
                if print_new_sugg:
                    # if all ok (if check returns True) we'll print the new_sugg also to gc_output
                    gc_output_text += formatted_change_md5 + "\t{}\t{}\n".format(
                        config["up_release"], config["dataset_name"]
                    )
                if conflict_text: #this one for all other cases
                    # to be able to remove conflicts even in multi-thread operation (the marking of conflict will not be shared across parallel workers)
                    conflict_output_text += conflict_text
            else:  # simply print the new_sugg
                gc_output_text += formatted_change_md5 + "\t{}\t{}\n".format(
                    config["up_release"], config["dataset_name"]
                )

            # update labelling info for all accessions in the suggested clade
            # accessions from the clade to which suggestion belongs to
            clade_accs = [
                prop_clade["prop_acc_info"][x]["acc"]
                for x in prop_clade["prop_acc_info"].keys()
            ]
            for acc in clade_accs:
                if "accrank" not in orthologs_dict[acc]:
                    orthologs_dict[acc]["accrank"] = str(suggestion_rank_level)
                else:
                    if "acctype" not in orthologs_dict[acc]:
                        orthologs_dict[acc]["accrank"] = str(suggestion_rank_level)
                if acc == p_acc:  # the suggestion
                    orthologs_dict[p_acc]["acctype"] = "prop_iso"
                    orthologs_dict[p_acc]["accrank"] = str(suggestion_rank_level)
                else:
                    if "acctype" not in orthologs_dict[acc]:
                        orthologs_dict[acc]["acctype"] = (
                            "canon_good"
                            if orthologs_dict[acc]["is_canonical"]
                            else "iso"
                        )

            # specifically mark the canon being replaced
            if "acctype" not in orthologs_dict[canon_acc]:
                orthologs_dict[canon_acc]["acctype"] = "canon_bad"
            if "accrank" not in orthologs_dict[canon_acc]:
                orthologs_dict[canon_acc]["accrank"] = str(suggestion_rank_level)

    # if there are no suggestions but one or more confirmed canonical solution(s), prepare labelling information
    if (
        suggestions_output_text == ""
        and confirmed_output_text != ""
        and config["create_pdf_files4confirmed_flag"]
    ):
        canonical_rank_level = 1
        # go through all canonical solutions
        for canonical_solution in confirmed_output_text[:-1].split("\n"):
            clade_members = canonical_solution.split("\t")[5]
            for taxon_acc in clade_members.split(","):
                taxon, acc = taxon_acc.split(":")
                # eprint(canonical_rank_level, taxon, acc) #debug
                orthologs_dict[acc]["accrank"] = str(canonical_rank_level)
                orthologs_dict[acc]["acctype"] = "canon_good"
            canonical_rank_level += 1

    # if there are suggestions (or only confirmations), prepare data for pdf creation
    if suggestions_output_text != "" or (
        confirmed_output_text != "" and config["create_pdf_files4confirmed_flag"]
    ):
        if config["create_lablt_files_flag"]:  # true if create_pdf_files_flag
            label_output = ""
            lab2_fname = os.path.join(config["lab_data_dir"], orthoid + ".lab_lt")
            label_output += "label\ttype\tlen\trank\tgroupid\n"  # header
            for acc in orthologs_dict.keys():
                acclabel = "{}{}|{}|{}".format(
                    orthologs_dict[acc]["entry_type"],
                    "" if orthologs_dict[acc]["is_canonical"] else "_iso",
                    acc,
                    orthologs_dict[acc]["org"],
                )
                if "acctype" in orthologs_dict[acc]:
                    acctype = orthologs_dict[acc]["acctype"]
                else:
                    acctype = "canon" if orthologs_dict[acc]["is_canonical"] else "iso"
                acclen = orthologs_dict[acc]["seqlen"]
                if "accrank" in orthologs_dict[acc]:
                    accrank = orthologs_dict[acc]["accrank"]
                else:
                    accrank = "0"  # no rank level
                if orthologs_dict[acc]["has_mane"]:
                    manetag = "/MANE_sel"
                else:
                    manetag = ""
                label_output += "{}\t{}{}\t{}\t{}\t{}\n".format(
                    acclabel,
                    acctype,
                    manetag,
                    acclen,
                    accrank,
                    orthologs_dict[acc]["groupid"],
                )

            with open(lab2_fname, "w", encoding="utf-8") as fd:
                fd.write(label_output)
            if verbose:
                eprint("\n" + label_output)

    if config["create_faX_files_flag"]:  # true if create_pdf_files_flag
        if (
            suggestions_output_text != "" or confirmed_output_text != ""
        ):  # if we created lab_lt files
            # to create .faX file in aln_data dir for pdf creation
            aln2fastax(orthoid, config["aln_data_dir"])

    if suggestions_output_text == "" and confirmed_output_text == "":
        if skipped_output_text == "":
            eprint("uncaught SKIPPED: {}".format(orthoid))
    elif suggestions_output_text != "" or confirmed_output_text != "":
        skipped_output_text = ""  # do not print skipped if we have confirm or change

    process_name = "o"  # scanning n_data to produce output
    process_duration = time.time() - process_start_time
    semaphore_file = os.path.join(config["semaphores_dir"], orthoid + ".done")
    with open(semaphore_file, "a") as fh:
        fh.write("{}\t{}\n".format(process_name, process_duration))
    if verbose:
        eprint("scan time: {}".format(elapsed_time(process_start_time)))

    # if nothing, returns ('', '', '', '', '')
    return (
        gc_output_text,
        suggestions_output_text,
        confirmed_output_text,
        skipped_output_text,
        conflict_output_text,
    )
