#!/usr/bin/env python
# coding: utf-8
"""
module providing trees manipulation functions
and in particular the function tree_to_ndata()
to identify low cost clades
"""
import os
import re
import sys
import copy
from .o2t_utils import eprint, get_orthologs_df_from_pantherid


REIDSTRIP = re.compile(r"\|\w*$")  # to remove oscode from ids
TAXA_RANK = {
    "HUMAN": 0,
    "MOUSE": 1,
    "RAT": 2,
    "BOVIN": 4,
    "CANLF": 5,
    "PIG": 7,
    "FELCA": 7,
    "MONDO": 6,
    "MACMU": 4,
    "GORGO": 4,
    "PANTR": 4,
    "ARATH": 0,
    "WHEAT": 1,
    "TOBAC": 3,
    "MEDTR": 4,
    "GOSHI": 5,
    "MAIZE": 6,
    "HORVV": 7,
    "ORYSJ": 8,
}


# tree leaves manipulation functions
def get_leaves(tree):
    leaves = []
    for clade in tree.find_clades(terminal=True):
        leaves.append(clade.name)

    return leaves


def get_leaf_taxa(leaf, verbose=False):
    """
    takes a seqid, possibly sp|P09488|GSTM1_HUMAN
    and returns oscode
    """
    fields = leaf.split("|")
    if (
        len(fields) == 3
    ):  # sp|P09488|GSTM1_HUMAN; tr|Q12345|HUMAN; ref_iso|NP_012345|HUMAN
        taxa = fields[2]
        if taxa.find("_") != -1:
            taxa = taxa.split("_")[1]  # full oscode
    else:
        if verbose:
            eprint(" *** Not leaf!! %s\n" % leaf)
        taxa = None
    return taxa


def get_leaf_acc(leaf_name, verbose=False):
    fields = leaf_name.split("|")
    if (
        len(fields) == 3
    ):  # sp|P09488|GSTM1_HUMAN; tr|Q12345|HUMAN; ref_iso|NP_012345|HUMAN
        acc = fields[1]
    else:
        if verbose:
            eprint(" *** Not leaf!! %s\n" % leaf_name)
        acc = None
    return acc


def get_leaf_acc_db(leaf_name, verbose=False):
    fields = leaf_name.split("|")
    db = None
    if (
        len(fields) == 3
    ):  # sp|P09488|GSTM1_HUMAN; tr|Q12345|HUMAN; ref_iso|NP_012345|HUMAN
        db = fields[0]  # sp, tr, sp_iso, tr_iso, ref_iso
        acc = fields[1]
    else:
        if verbose:
            eprint(" *** Not leaf!! %s\n" % leaf_name)
        acc = None
    return acc, db


def get_leaf_acc_db_taxon(leaf_name, verbose=False):
    fields = leaf_name.split("|")
    db = acc = taxon = None
    if len(fields) == 2:  # NP_12345|HUMAN
        acc = fields[0]
        db = "ref"
        taxon = fields[1]
    elif len(fields) == 3:  # sp|P09488|GSTM1_HUMAN; tr|Q12345|HUMAN
        db = fields[0]  # sp, tr, sp_iso, tr_iso
        acc = fields[1]
        taxon = fields[2]
        if taxon.find("_") != -1:
            taxon = taxon.split("_")[1]  # full oscode
    else:
        if verbose:
            eprint(" *** Not leaf!! %s\n" % leaf_name)
    return acc, db, taxon


def get_leaf_taxa_set(leaves, taxa_set={}):
    """
    returns taxa_set, where keys are distinct taxa, and values are lists of leaves with that taxa
    """
    for leaf in leaves:
        if isinstance(leaf, str):
            taxa = get_leaf_taxa(leaf)
        else:
            taxa = get_leaf_taxa(leaf.name)

        if taxa:
            if taxa not in taxa_set:
                taxa_set[taxa] = [leaf]
            else:
                taxa_set[taxa].append(leaf)

    return taxa_set


def get_accs_from_subtree(subtree):
    """
    get the accessions in a subtree
    """
    accs = []
    for leaf in subtree["leaves"]:
        accs.append(get_leaf_acc(leaf.name))
    return accs


def selected_clade_cost(tree, start_clade, selected_leaves):
    """
    find a path from start_clade to each of the selected_leaves,
    building a structure that ensures that each node on the path is only
    counted once
    """
    t_cost = 0.0

    clade_seen = {}

    for leaf in selected_leaves:
        this_trace = tree.trace(start_clade, leaf)
        # sys.stderr.write(repr(this_trace),"\n"))
        for t_clade in this_trace[1:]:
            if t_clade.name not in clade_seen:
                t_cost += t_clade.branch_length
                clade_seen[t_clade.name] = True

    return t_cost


def tree_clade_cost(tree, sel_clade=None, canonicals=None, debuginfo=False, orthoid=""):
    """
    return the cost of a particular clade
    (not a real clade, but a set of leaves, one per taxa)
    or a specified set of accessions
    """
    if sel_clade is None and canonicals is None:
        eprint("you must pass me a sel_clade and a list of canonical accessions")
        sys.exit(1)
    accs = get_accs_from_subtree(sel_clade)  # target accs
    if not accs:  # empty
        eprint("{} ERROR: cannot retrieve accessions for given clade?!".format(orthoid))
        return "NA"

    canonical_accs = canonicals.values()
    if len(accs) != len(canonical_accs):
        eprint(
            "{} WARNING: number of canonical accs not equal to number of members of the set!".format(
                orthoid
            )
        )
        # in this case we must skip computing canon_diff for the taxon missing canonical

    target_leaves = {}
    canonical_leaves = {}

    # get leaves for target and for canonicals
    for leaf in tree.get_terminals():
        # eprint("leaf name is {}".format(leaf.name)) #debug
        (db, acc, id_oscode) = leaf.name.split("|")
        if acc in accs:
            target_leaves[get_leaf_taxa(leaf.name)] = leaf
        if acc in canonical_accs:
            canonical_leaves[get_leaf_taxa(leaf.name)] = leaf

    # paranoid sanity checks:
    if len(target_leaves) != len(accs):
        eprint("{} ERROR: could not retrieve leaves for the targets?!".format(orthoid))
        return "NA"
    if len(canonical_leaves) != len(canonical_accs):
        eprint(
            "{} ERROR: could not retrieve leaves for the canonicals?!".format(orthoid)
        )
        return "NA"

    target_taxa_list = target_leaves.keys()

    # procedure to compute costs:
    # have the targets, find the for that clade
    # find the LCA for the selected leaves
    # find the cost to the selected leaves from the LCA

    # 1) compute cost of the selected clade
    s_common_clade = tree.common_ancestor(target_leaves.values())
    target_cost = selected_clade_cost(tree, s_common_clade, target_leaves.values())
    if debuginfo:
        eprint("targets: %s cost:%.4f\n" % (",".join(accs), target_cost))

    # 2) compute cost of the canonicals
    c_common_clade = tree.common_ancestor(canonical_leaves.values())
    canonical_cost = selected_clade_cost(
        tree, c_common_clade, canonical_leaves.values()
    )
    if debuginfo:
        eprint(
            "canonicals: %s cost:%.4f\n" % (",".join(canonical_accs), canonical_cost)
        )

    # 3) compute the targets cost if replacing a single taxon with canonical
    # (resulting differential cost written into leaf.canon_diff)
    set_to_try = target_leaves.copy()  # take a copy of the target leaves

    for taxon in target_leaves.keys():
        target_leaf = target_leaves[taxon]
        if taxon not in canonical_leaves:
            # no canonical to replace (case of orphan isoform), set as 0 and continue
            target_leaf.canon_diff = 0.0
            continue
        canonical_leaf = canonical_leaves[taxon]
        target_acc = target_leaf.name
        canonical_acc = canonical_leaf.name
        if target_acc == canonical_acc:  # target is the canonical
            target_leaf.canon_diff = 0.0
        else:
            # replacing targets with canonicals one at a time:
            old_leaf = set_to_try[taxon]
            # replace one target with the corresponding canonical
            set_to_try[taxon] = canonical_leaf
            common_clade = tree.common_ancestor(set_to_try.values())
            set_cost = selected_clade_cost(tree, common_clade, set_to_try.values())

            set_to_try[taxon] = old_leaf

            diff_cost = target_leaf.canon_diff = set_cost - target_cost
            if debuginfo:
                eprint(
                    "set_to_try(%s): %s cost:%.4f\n"
                    % (
                        taxon,
                        ",".join([x.name for x in set_to_try.values()]),
                        diff_cost,
                    )
                )
            if diff_cost < 0.000001:
                if debuginfo:
                    eprint(
                        " *** replacing isoform: {} {} {}\n".format(
                            taxon, target_leaves[taxon], canonical_leaf
                        )
                    )
                target_leaves[taxon] = canonical_leaf
                target_leaves[taxon].canon_diff = 0.0

    # save the cost information, and appropriate taxon/db/acc/_oscodes in selected_clades for later display
    return {
        "leaves": target_leaves.values(),
        "cost": target_cost,
        "canonical_cost": canonical_cost,
        "canonical_dict": canonicals,
        "clade": s_common_clade,
        "taxa": target_taxa_list,
    }


def get_tree_canonicals(leaves, taxa_list_set, orthogroup):
    """
    Get the canonicals for all taxa in a tree
    In case more than one canonical per taxon, return them as set
    """
    canonicals = {}
    repeated_count = 0
    repeated_taxa = {}
    for leaf_name in leaves:
        acc, db, taxon = get_leaf_acc_db_taxon(leaf_name)
        if db in ["sp", "tr"]:
            if taxon in canonicals:
                if taxon in repeated_taxa:
                    repeated_taxa[taxon].add(acc)
                else:
                    repeated_taxa[taxon] = {canonicals[taxon], acc}
                # eprint("WARNING: we already have canonical for {} ({}), now replacing with {}".format(taxon, canonicals[taxon], acc))
                repeated_count += 1
            else:
                canonicals[taxon] = acc
    canonicals_taxa_set = set(canonicals.keys())
    if canonicals_taxa_set != taxa_list_set:
        eprint(
            "{} WARN: missing canonical(s) for {}".format(
                orthogroup, sorted(taxa_list_set.difference(canonicals_taxa_set))
            )
        )
    return canonicals, repeated_count, repeated_taxa


def get_gc_clade_canonicals(sel_clade, orthologs_df):
    """
    Return for each taxa the corresponding canonical from the same genecentric group as the accession in the clade
    """
    gc_clade_canonicals = {}
    # only the canonicals
    orthologs_canonicals = orthologs_df[orthologs_df["is_canonical"]]

    for leaf in sel_clade["leaves"]:
        acc, db, taxon = get_leaf_acc_db_taxon(leaf.name)
        if db in ["sp", "tr"]:  # already canonical, copy it as it is
            gc_clade_canonicals[taxon] = acc
        else:  # isoform, look for corresponding canonical in same gc
            can = None
            if (
                acc not in orthologs_df.index
            ):  # paranoid check! acc SHOULD be in the orthologs_df!
                eprint("ERROR! acc {} not in the df".format(acc))
                return gc_clade_canonicals
            acc_groupid = orthologs_df.loc[acc]["groupid"]
            # find the canonical with same groupid
            a = orthologs_canonicals[
                orthologs_canonicals["groupid"] == acc_groupid
            ].index
            if len(a):
                can = a.to_list()[0]
            if can is not None:
                gc_clade_canonicals[taxon] = can
            else:
                eprint("WARNING: orphan isoform '{}' without canonical".format(acc))

    return gc_clade_canonicals


def lookup_by_names(tree):
    names = {}
    for clade in tree.find_clades():
        if clade.name is None:
            clade.name = "Root"

        if clade.name in names:
            raise ValueError("Duplicate key: %s" % clade.name)
        else:
            names[clade.name] = clade

    return names


def lookup_by_taxa(names_d):
    """
    return dict of clades by taxa name
    """
    taxa_clades = {}

    for this_name in names_d.keys():
        this_leaf = names_d[this_name]
        this_taxon = get_leaf_taxa(this_leaf.name)
        if not this_taxon:
            continue
        if this_taxon in taxa_clades:
            taxa_clades[this_taxon].append(this_leaf)
        else:
            taxa_clades[this_taxon] = [this_leaf]


def bt_gap_dist(s1, s2):
    if len(s1) != len(s2):
        sys.stderr.write(
            "*** simple_dist lengths differ %s: %d :: %s: %d\n"
            % (s1.id, len(s1), s2.id, len(s2))
        )
        return 1.0


# helper functions of tree_to_ndata
def wrp_is_semiterminal(self_clade):
    """Check if all direct descendents are terminal."""
    if self_clade.root.is_terminal():
        return False
    for clade in self_clade.root.clades:
        if clade.is_terminal():
            return True
    return False


def wrp_collapse(self_clade, target=None, **kwargs):
    """Delete target from the tree, relinking its children to its parent.

    :returns: the parent clade.

    """
    path = self_clade.get_path(target, **kwargs)
    if not path:
        raise ValueError("couldn't collapse %s in this tree" % (target or kwargs))
    if len(path) == 1:
        parent = self_clade.root
    else:
        parent = path[-2]
    popped = parent.clades.pop(parent.clades.index(path[-1]))
    extra_length = popped.branch_length or 0
    extra_length /= len(popped)

    for child in popped:
        child.branch_length += extra_length
    parent.clades.extend(popped.clades)
    return parent


def wrp_collapse_all(self_clade, target=None, **kwargs):
    # 21-Sept-2023 wrp
    # after collapse, turn all very short terminal branches into zero-length branches
    # this solves a problem where zero length branches somehow have short lengths after collapse
    # for t_clade in self_clade.get_terminals():
    #     if (t_clade.branch_length < 0.00001):
    #         t_clade.branch_length = 0.0
    # Read the iterable into a list to protect against in-place changes
    matches = list(self_clade.find_clades(target, False, "level", **kwargs))
    if not matches:
        # No matching nodes to collapse
        return
    # Skip the root node -- it can't be collapsed
    if matches[0] == self_clade.root:
        matches.pop(0)
    for clade in matches:
        wrp_collapse(self_clade, clade)


def check_taxa(name):
    """figures out which species the taxa is from, based off a string of the name"""
    if name[2] == "_":
        letter = name.split("|")[-1][0].upper()
    else:
        letter = name.split("_")[-1][0].upper()
    return letter


def all_parents(tree):
    """find all parents and save in dict"""
    parents = {}
    for clade in tree.find_clades(order="level"):
        for child in clade:
            if not child.name:
                parents["Root"] = clade
            else:
                parents[child.name] = clade
    return parents


def get_db_rank(name):
    """give rank to clade name based on db"""
    if re.search(r"iso", name):
        rank = 2
    elif name[0:2] == "sp":
        rank = 0
    elif name[0:2] == "tr":
        rank = 1
    else:
        rank = 3

    return rank


def sort_good_taxa(old_taxa_set):
    """given a list of taxa, sort by branch_length and db_rank"""
    new_taxa_set = {}

    for clade in old_taxa_set.keys():
        new_taxa_set[clade] = {}
        for tx in old_taxa_set[clade].keys():
            new_taxa_set[clade][tx] = sorted(
                old_taxa_set[clade][tx],
                key=lambda x: (x.branch_length, get_db_rank(x.name)),
            )

    return new_taxa_set


def get_taxa_rank(name, taxa=None):
    if taxa is None:
        taxa = get_leaf_taxa(name)

    if taxa and taxa in TAXA_RANK:
        return TAXA_RANK[taxa]
    else:
        TAXA_RANK[taxa] = 10
        return 10


def db_rank(leaves, n_taxa):
    """produce a score reflecting the number of canonical, isoform, refseq sequences"""
    db_score = 0
    for leaf in leaves:
        if re.search(r"iso", leaf.name):
            db_score += n_taxa * get_taxa_rank(leaf.name)
        elif leaf.name[0:2] == "sp" or leaf.name[0:2] == "tr":
            db_score += get_taxa_rank(leaf.name)
        else:
            db_score += n_taxa * n_taxa * get_taxa_rank(leaf.name)

    return db_score


def select_best_leaves(clade_set, verbose=False):
    """
    given a clade with multiple children,
    select a set of children from each taxa preferring
    canonical/isoform/refseq
    """
    selected_leaves = []

    for taxa in clade_set.keys():
        if isinstance(clade_set[taxa], list):
            if len(clade_set[taxa]) > 1:
                rank_child_list = []
                for child in clade_set[taxa]:
                    rank_child_list.append(
                        {
                            "taxa": taxa,
                            "clade": child,
                            "rank": get_db_rank(child.name),
                        }
                    )
                rank_child_list.sort(
                    key=lambda x: (x["clade"].branch_length, x["rank"])
                )
                selected_leaves.append(rank_child_list[0]["clade"])
            else:
                selected_leaves.append(clade_set[taxa][0])
        else:
            selected_leaves.append(clade_set[taxa])

    return selected_leaves


def list_missing_taxa(all_taxa, found_taxa):
    """return list of missing taxa"""
    missing_taxa = []
    for taxon in all_taxa:
        if taxon not in found_taxa:
            missing_taxa.append(taxon)

    return missing_taxa


def print_ndata(select):
    """
    given a properly formatted selected clade/set, print it out in n_data2 format

    Example of output line and explanation of the fields:
    Inner8 :: 0.0000/3 :: CANLF:tr_iso|A0A5F4D1K6 0.2535,HUMAN:sp_iso|Q9UM11-3 0.1835,MOUSE:tr_iso|D3YTV2 0.1806 :: 0.0796/3 :: CANLF:A0A5F4C1U6,HUMAN:Q9UM11,MOUSE:Q9R1K5
    CLADENAME :: CLADECOST/CLADETAXANUM :: TAXON:db|ACC diff_cost,TAXON:db|ACC diff_cost[...] :: CANONCOST/CANONTAXANUM :: TAXON:CANONACC,TAXON:CANONACC[...]
    """
    (clade, leaves, cost) = (select["clade"], select["leaves"], select["cost"])
    if (cost < 0.000):
        cost = 0.000
    gc_canonicals_cost = select["canonical_cost"]
    gc_clade_canonicals = select["canonical_dict"]

    leaf_info = []
    for leaf in leaves:
        taxon = get_leaf_taxa(leaf.name)
        # leaf_info.append("{}:{} {:.4f}".format(taxon, REIDSTRIP.sub('|_',leaf.name), leaf.branch_length)) #CHANGED (stripped oscode)
        # leaf_info.append("{}:{} {:.4f}".format(taxon, REIDSTRIP.sub('',leaf.name), leaf.branch_length)) #OLD: printing branch_length
        leaf_info.append(
            "{}:{} {:.4f}".format(taxon, REIDSTRIP.sub("", leaf.name), leaf.canon_diff)
        )  # NEW: replaced branch_length with canon_diff

    return "{} :: {:.4f}/{} :: {} :: {:.4f}/{} :: {}\n".format(
        clade.name,
        cost,
        len(leaves),
        ",".join(sorted(leaf_info)),
        gc_canonicals_cost,
        len(gc_clade_canonicals),
        ",".join([k + ":" + v for k, v in sorted(gc_clade_canonicals.items())]),
    )


def remove_to_improve(
    in_tree,
    t_clade,
    level,
    max_rlevel,
    min_taxa,
    new_selected,
    old_tested,
    dropped,
    tree_drop_fact=1.5,
    verbose=False,
):
    """
    recursively check to see if removing a taxon improves score
    adds to new_selected list with t_clade/t_cost
    returns new_selected, but also updates old_tested
    """
    zero_eps = 0.0000001

    if level > max_rlevel:
        return new_selected

    level += 1

    if t_clade["cost"] <= zero_eps:
        return new_selected

    if len(t_clade["leaves"]) > min_taxa:
        best_cost = t_clade["cost"]
        best_miss = None
        best_out = None
        this_dropped = [x for x in dropped]
        this_dropped.append("")

        # check for improvement leaving out each taxa
        for t_leaf in t_clade["leaves"]:
            this_dropped[-1] = t_leaf.name
            t_selected = [x for x in t_clade["leaves"] if x.name != t_leaf.name]

            # build ordered set of acc's to avoid duplication
            t_select_names = ":".join(
                [t.name for t in sorted(t_selected, key=lambda x: x.name)]
            )
            if t_select_names in old_tested:
                continue
            else:
                old_tested[t_select_names] = True

            t_parent = in_tree.common_ancestor(t_selected)
            t_cost = selected_clade_cost(in_tree, t_parent, t_selected)

            nt_cost = t_cost / float(len(t_selected))
            if verbose:
                eprint(
                    "|- D [%d:%d]: %s : %.5f %.5f : %s"
                    % (
                        level,
                        len(t_selected),
                        ",".join(this_dropped),
                        t_cost,
                        nt_cost,
                        t_select_names,
                    )
                )

            if t_cost < best_cost:
                # removing a taxa has improved the cost
                best_cost = t_cost
                best_miss = [x for x in t_selected]
                best_out = t_leaf.name
                best_taxa_out = get_leaf_taxa(t_leaf.name)

        # only save lower-taxa solution if it is better
        # this check should keep zero-cost clades from being improved
        # -- they should be improved if canonical available
        # drop_thresh =
        # max(max_cost/5.0,(t_clade['cost']/float(len(t_clade['leaves']))))
        # * tree_drop_fact if (best_out and t_clade['cost'] -
        # best_cost > drop_thresh):
        drop_thresh = 1.0 / tree_drop_fact

        n_taxa = len(t_clade["taxa"])
        if (
            best_out
            and t_clade["cost"] > 0.001
            and best_cost * n_taxa < t_clade["cost"] * (n_taxa - 1)
        ):
            this_dropped[-1] = best_out

            if verbose:
                eprint(
                    "|- [%d:%d] remove %s : %.4f - %.4f > %.4f"
                    % (
                        level,
                        len(t_clade["taxa"]),
                        ",".join(this_dropped),
                        t_clade["cost"],
                        best_cost,
                        drop_thresh,
                    )
                )

            # this has not worked, because taxa are not what is in best_out, should have been best_taxa_out
            t_taxa = [x for x in t_clade["taxa"] if x != best_taxa_out]

            # check to see if the new set of taxa uses a different clade

            common_clade = in_tree.common_ancestor(best_miss)
            t_cost = selected_clade_cost(in_tree, common_clade, best_miss)
            db_score = db_rank(best_miss, n_taxa)

            t_clade1 = {
                "leaves": best_miss,
                "cost": t_cost,
                "clade": common_clade,
                "db_rank": db_score,
                "taxa": t_taxa,
            }

            new_selected.append(t_clade1)

            # we have an n-1 solution, how about n-2??
            new_selected = remove_to_improve(
                in_tree,
                t_clade1,
                level,
                max_rlevel,
                min_taxa,
                new_selected,
                old_tested,
                this_dropped,
                tree_drop_fact=tree_drop_fact,
                verbose=verbose,
            )

    return new_selected


def tree_to_ndata(
    in_tree,
    ortho_df=None,
    config={},
    filename_prefix=None,
    target_id=None,
    verbose=False,
    debuginfo=False,
    create_lab_file=False,
    create_ndata_file=False,
    drop_taxa=True,
):
    # max_cost: exclude solutions with costs higher (default=0.01)
    max_cost = config.get("tree_max_cost", 0.01)
    tree_drop_fact = config.get("tree_drop_fact", 1.5)

    if create_lab_file or create_ndata_file:
        if filename_prefix is None:
            raise Exception(
                "    => ERROR: you need to specify a filename_prefix if you want to create files"
            )
            return
        # output files
        lab2_fname = os.path.join(config["lab_data_dir"], filename_prefix + ".lab2")
        ndata_fname = os.path.join(config["n_data_dir"], filename_prefix + ".n_data2")

    # output variables
    n_data_in = ""
    labels = ""

    if target_id is None:
        eprint("No target id\n")

    selected_clades = []

    # we are given a tree:
    in_clade_names = lookup_by_names(in_tree)

    l_rank = 0
    for sel in selected_clades:
        l_rank += 1

    zero_tree = copy.deepcopy(in_tree)
    cost_tree = copy.deepcopy(in_tree)

    zero_eps = 0.0000001
    wrp_collapse_all(zero_tree, lambda c: c.branch_length <= zero_eps)

    wrp_collapse_all(cost_tree, lambda c: c.branch_length <= max_cost)

    leaves = get_leaves(in_tree)
    # important: init with new hash
    all_taxa_set = get_leaf_taxa_set(leaves, {})

    # for canonicals reference cost computations:
    orthologs_df = get_orthologs_df_from_pantherid(
        target_id, ortho_df, remove_outliers=True
    )[
        ["groupid", "is_canonical"]
    ]  # for get_gc_clade_canonicals()

    min_taxa = n_taxa = len(all_taxa_set.keys())
    min_taxa = min(min_taxa, config["taxa_threshold"])

    # first check to see if we have clades in zero_tree with >= min_taxa different taxa; if so, we are done

    # 13-Sept-2023 (wrp) -- check to make certain that collapsed clades have < zero_eps length to parent
    zc_taxa = scan_clade_taxa(zero_tree, zero_eps, min_taxa)

    # zc_taxa is a dict with keys:preterminal clades and values:dicts with keys:taxa and values: lists of clades
    # sort from largest number of taxa to smallest
    zc_clade_srt = sorted(
        zc_taxa.keys(), key=lambda tx: len(zc_taxa[tx].keys()), reverse=True
    )

    # check for clades with >= min_taxa
    good_clades = []

    # starting checking zero cost clades from zero cost trees
    for zc in zc_clade_srt:
        if len(zc_taxa[zc].keys()) >= min_taxa:
            good_clades.append({"taxa": zc_taxa[zc], "clade": zc, "phase": "zero"})

    if debuginfo:
        eprint("** %s in_tree cost %s %.2f\n" % (target_id, in_tree, max_cost))
        # eprint(" % s\n"%(str(cost_tree)))

    # 13-Sept-2023 (wrp)
    # like the similar search of zero_tree, we need to skip long branches in collapsed clades
    cc_taxa = scan_clade_taxa(cost_tree, max_cost, min_taxa)

    # sort by how many taxa
    cc_clade_srt = sorted(
        cc_taxa.keys(), key=lambda tx: len(cc_taxa[tx].keys()), reverse=True
    )

    # now checking non-zero cost clades/trees
    for cc in cc_clade_srt:
        if len(cc_taxa[cc].keys()) >= min_taxa:
            good_clades.append({"taxa": cc_taxa[cc], "clade": cc, "phase": "cost"})

    if debuginfo:
        eprint(" good clades: %d\n" % (len(good_clades)))

    # if we have some good clades, do some selecting

    selected_clades = []
    sel_clade_dict = {}

    # here, good_clades map to either zero_tree or cost_tree
    # they need to be re-mapped back to in_tree, and checked to remove duplicates

    for gc in good_clades:
        these_taxa = list(gc["taxa"].keys())
        selected_leaves = select_best_leaves(gc["taxa"])
        selected_names = [x.name for x in selected_leaves]

        # sel_clade_key is accessions sorted by accession

        sel_clade_key = ":".join(sorted(selected_names))
        if sel_clade_key not in sel_clade_dict:
            # get parent of selected leaves in in_tree
            in_selected_leaves = [in_clade_names[x.name] for x in selected_leaves]
            in_parent = in_tree.common_ancestor(in_selected_leaves)

            gc_cost = selected_clade_cost(in_tree, in_parent, in_selected_leaves)
            db_score = db_rank(selected_leaves, n_taxa)

            sel_info = {
                "leaves": in_selected_leaves,
                "cost": gc_cost,
                "clade": in_parent,
                "db_rank": db_score,
                "taxa": these_taxa,
                "sel_hash": sel_clade_key,
                "phase": gc["phase"],
            }

            selected_clades.append(sel_info)
            sel_clade_dict[sel_clade_key] = sel_info
        else:
            # this should not happen, as the sel_clade_key is the set of leaf.names
            if gc["phase"] == sel_clade_dict[sel_clade_key]["phase"]:
                sys.stderr.write(
                    "*** %s duplicate sel_clade_key: %s [%s] %s\n"
                    % (target_id, sel_clade_key, gc["phase"], repr(in_parent))
                )

    # check for N-taxa top_parent
    these_taxa = {}

    # check for the possibility of a much better cost tree with one less taxon
    # len(selected_clades) is the number of candidate solutions
    if drop_taxa and len(selected_clades) > 0:
        new_selected = []
        old_tested = {}

        for t_clade in selected_clades:
            if (t_clade['cost'] <= zero_eps):
                continue
            # recursively updated new_selected, old_tested
            dropped = []
            new_selected = remove_to_improve(
                in_tree,
                t_clade,
                0,
                4,
                min_taxa,
                new_selected,
                old_tested,
                dropped,
                tree_drop_fact=tree_drop_fact,
                verbose=verbose,
            )

        if len(new_selected) > 0:
            for new in new_selected:
                in_selected_names = [x.name for x in new["leaves"]]
                sel_clade_key = ":".join(sorted(in_selected_names))

                if sel_clade_key not in sel_clade_dict:
                    new["phase"] = "remove"
                    new["sel_hash"] = sel_clade_key
                    selected_clades.append(new)
                elif verbose:
                    old_phase = sel_clade_dict[sel_clade_key]["phase"]
                    sys.stderr.write(
                        "*** %s duplicate sel_clade_key: %s [%s/%s]\n"
                        % (target_id, sel_clade_key, old_phase, "remove")
                    )

    if len(selected_clades) > 0:
        selected_sort = sorted(
            selected_clades, key=lambda x: (x["cost"], -len(x["leaves"]), x["db_rank"])
        )

        # ndata output
        n_data_in = n_data_info = ">{}\n".format(target_id)
        if create_ndata_file:
            nd_fd = open(ndata_fname, "w")
            nd_fd.write(n_data_info)
        for select in selected_sort:  # for each selected clade
            # to extract original costs rather than collapsed ones
            gc_clade_canonicals = get_gc_clade_canonicals(select, orthologs_df)
            newselect = tree_clade_cost(
                in_tree,
                sel_clade=select,
                canonicals=gc_clade_canonicals,
                orthoid=target_id,
                debuginfo=debuginfo,
            )
            if newselect != "NA":
                n_data_info = print_ndata(newselect)
                n_data_in += n_data_info
                if create_ndata_file:
                    nd_fd.write(n_data_info)
        if create_ndata_file:
            nd_fd.close()

        if create_lab_file:
            lfd = open(lab2_fname, "w")
            lfd.write("label\trank\n")

        # lab2 output
        l_rank = 0
        for sel in selected_sort:
            l_rank += 1
            for leaf in sel["leaves"]:
                label = " % s\t%d\n" % (leaf.name, l_rank)
                labels += label
                if create_lab_file:
                    lfd.write(label)
        if create_lab_file:
            lfd.close()

    return in_tree, n_data_in, labels


# helper functions of scan_clade_taxa
def wrp_is_preterminal(zc):
    for c in zc:
        if c.is_terminal:
            return True


def eval_clade(zc, zc_taxa, min_taxa, this_cost):
    these_taxa = get_leaf_taxa_set(zc, {})

    zc_taxa_set = {}

    # if not enough taxa, quit
    if len(these_taxa) < min_taxa:
        return zc_taxa

    # there are enough taxa, but is the cost low enough?
    for zt in these_taxa.keys():
        for ztc in sorted(these_taxa[zt], key=lambda x: x.branch_length):
            if ztc.branch_length < this_cost:
                if zt not in zc_taxa_set:
                    zc_taxa_set[zt] = [ztc]
                else:
                    zc_taxa_set[zt].append(ztc)

    if len(zc_taxa_set) >= min_taxa:
        zc_taxa[zc.name] = zc_taxa_set

    return zc_taxa


def scan_clade_taxa(this_tree, this_cost, min_taxa):
    """
    given a clade, extract the taxa with branch lengths < this_cost
    needs to check also check taxa immediately below root
    modified to check for min_taxa during scan
    """
    zc_taxa = {}
    for zc in this_tree.find_clades():
        if zc.is_terminal():
            continue
        elif wrp_is_preterminal(zc):
            zc_taxa = eval_clade(zc, zc_taxa, min_taxa, this_cost)
        elif zc.name == "Root":
            zc_taxa = eval_clade(zc, zc_taxa, min_taxa, this_cost)

    return zc_taxa
