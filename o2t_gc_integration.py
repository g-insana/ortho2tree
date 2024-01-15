#!/usr/bin/env python
# coding: utf-8
"""
module providing functions involved in integration with the genecentric pipeline
"""
import os
import sys
import pandas as pd
from o2t_utils import eprint


prev_changes_columns = [
    "pantherid",
    "org",
    "canon_acc",
    "canon_len",
    "canon_md5",
    "prop_acc",
    "prop_len",
    "prop_md5",
    "rank_score",
    "canon_cost",
    "prop_cost",
    "n_sp",
    "n_tax",
    "n_canon",
    "wn_canon",
    "scaled_prop_f",
    "scaled_p_cost",
    "p_len_diff",
    "clade",
    "clade_members",
    "MANEstatus",
    "release",
    "dataset",
]


def check_suggestion_against_prev(
    prevgc_group_df, taxon, canon_acc, p_acc, orthoid, verbose=True
):
    """
    prev_sugg (if it exists) and new_sugg get merged with the following rules:
    1) if new_sugg and prev_sugg are same, prev_sugg will be printed
    2) if new_sugg contains a "flip flop", i.e. a pair where we suggest A->B where prev_sugg had "B->A",
       prev_sugg will be printed and new_sugg ignored
    3) if new_sugg contains a novel suggestion for the same panther_id and same taxon involving one of the accessions from prev_sugg (e.g. B->C),
       then the suggestion "A->B" from prev_sugg is removed and the new one is printed instead
    4-5) otherwise they don't conflict and they both get printed (e.g. different clade even if same taxon)

    returns:
      True/False (whether the new_sugg should get printed out to gc_output)
      conflict_text (in case we need to remove the prev_sugg) (due to parallel processing, it's not enough to mark it as conflict in the df)
    """
    # first check if there is prev_sugg for this orthogroup and taxon
    prev_sugg_taxon = prevgc_group_df[prevgc_group_df["org"] == taxon]
    if len(
        prev_sugg_taxon
    ):  # GCTODO check what happens if there are multiple prev and new suggestions for same taxon?
        prevgc_match = prevgc_group_df[
            prevgc_group_df[["oldcanon", "replacement"]]
            .isin([canon_acc, p_acc])
            .any(axis=1)
        ]
        if len(
            prevgc_match
        ):  # if there is one (or more?) suggestion mentioning one of canon_acc or p_acc
            if len(
                prevgc_group_df[
                    (prevgc_group_df["oldcanon"] == canon_acc)
                    & (prevgc_group_df["replacement"] == p_acc)
                ]
            ):
                # case 1: prev_sugg identical to old_sugg -> old_sugg is printed
                # if verbose: eprint("INT new_sugg seen before: {}->{}".format(canon_acc, p_acc))
                return (
                    False,
                    "",
                )  # the old_sugg continues to get printed instead of the identical new_sugg
            elif len(
                prevgc_group_df[
                    (prevgc_group_df["oldcanon"] == p_acc)
                    & (prevgc_group_df["replacement"] == canon_acc)
                ]
            ):
                # case 2: new_sugg is a flip flop of a previous suggestion
                if verbose:
                    eprint(
                        "INT flip flop: {}<>{} [{}]".format(canon_acc, p_acc, orthoid)
                    )
                return (
                    False,
                    "",
                )  # old_sugg continues to get printed instead of the flipflop new_sugg
            else:  # one (or more) prev_sugg conflict with new_sugg: remove prev_sugg(s)
                if verbose:
                    eprint(
                        "INT new_sugg {}->{} conflicts with prev_gc: {} [{}]".format(
                            canon_acc,
                            p_acc,
                            prevgc_match[["oldcanon", "replacement"]].values,
                            orthoid,
                        )
                    )
                # case 3: new_sugg contains novel suggestion conflicting with prev_sugg, so we need to remove prev_sugg
                # we remove prev_sugg while new_sugg is printed instead
                prevgc_group_df.loc[prevgc_match.index, "conflict"] = True
                # remove prev_sugg in multithread processing
                return True, "{}".format(
                    prevgc_group_df.loc[
                        prevgc_match.index, ["oldcanon", "replacement"]
                    ].to_csv(index=False, header=False, sep="\t")
                )
        else:
            if verbose:
                eprint("INT case 4 for {} [{}]".format(canon_acc, orthoid))
            # case 4: prev_sugg and new_sugg both printed as they don't conflict (for now, maybe next sugg for same taxon will cause a conflict of case 3)
            return (
                True,
                "",
            )  # we print new_sugg as it contains not conflicting novel suggestion
    else:  # no prev_sugg for this taxon
        if verbose:
            eprint("INT case 5 for {} [{}]".format(canon_acc, orthoid))
        # the new_sugg is printed (totally new as for this orthoid+taxon there was no prev_sugg)
        return True, ""


def read_prev_changes(ortho_df, config=None):
    prevgc_file = config["prevgc_file"]
    eprint("\nINTEGRATION: Reading prevgc output file '{}'".format(prevgc_file))
    prevgc_df = pd.read_csv(
        prevgc_file,
        sep="\t",
        header=None,
        comment="#",
        skip_blank_lines=True,
        names=prev_changes_columns,
    )

    # remove duplicates, caused by a change in pthr family assignment across analyses, e.g.:
    # PTHR11461:SF13	HUMAN	sp|P01019
    # PTHR11461:SF379	HUMAN	sp|P01019
    prevgc_df.drop_duplicates(
        subset=["org", "canon_acc", "prop_acc"], keep="last", inplace=True
    )

    eprint("  prevgc file contains {} unique suggestions".format(len(prevgc_df)))
    prevgc_df[["oldcanon_type", "oldcanon"]] = prevgc_df["canon_acc"].str.split(
        "|", expand=True
    )
    prevgc_df[["replacement_type", "replacement"]] = prevgc_df["prop_acc"].str.split(
        "|", expand=True
    )

    # CLEANING:
    # 1) entries where either cano or suggestion are no more present in genecentric are removed
    prevgc_notfound = prevgc_df[
        ~prevgc_df[["oldcanon", "replacement"]].isin(ortho_df.index).all(axis=1)
    ]
    prevgc_notfound_count = len(prevgc_notfound)
    eprint(
        "INT There are at least {} accessions from prevgc file that cannot be found in current genecentric data. Dropping them from prevgc file".format(
            prevgc_notfound_count
        )
    )
    prevgc_df.drop(
        prevgc_notfound.index, inplace=True
    )  # remove those which do not pass check #1

    # 2) also if the oldcanon and its replacement are no more in the same gc group
    prevgc_df.set_index("oldcanon", inplace=True)
    prevgc_df["can_groupid"] = ortho_df.loc[
        prevgc_df[prevgc_df.index.isin(ortho_df.index)].index, "groupid"
    ]
    prevgc_df.reset_index(inplace=True)
    prevgc_df.set_index("replacement", inplace=True)
    prevgc_df["prop_groupid"] = ortho_df.loc[
        prevgc_df[prevgc_df.index.isin(ortho_df.index)].index, "groupid"
    ]
    prevgc_df.reset_index(inplace=True)
    groupid_check = prevgc_df.can_groupid == prevgc_df.prop_groupid
    prevgc_notfound = prevgc_df[~groupid_check]
    if len(prevgc_notfound):
        eprint(
            "INT There are at least {} accessions from prevgc file where the genecentric groupid is different between old_canon and its replacement".format(
                len(prevgc_notfound)
            )
        )

    prevgc_df = prevgc_df[groupid_check]  # remove those that do not pass check #2

    # 3) also if the md5sum information changed (i.e. sequence changed)
    canon_md5_check = prevgc_df[["oldcanon", "canon_md5"]].rename(
        columns={"oldcanon": "acc", "canon_md5": "md5sum"}
    ).set_index("acc") == pd.DataFrame(ortho_df.loc[prevgc_df["oldcanon"]]["md5sum"])
    canon_md5_changed = canon_md5_check[~canon_md5_check.md5sum].index.values
    prop_md5_check = prevgc_df[["replacement", "prop_md5"]].rename(
        columns={"replacement": "acc", "prop_md5": "md5sum"}
    ).set_index("acc") == pd.DataFrame(ortho_df.loc[prevgc_df["replacement"]]["md5sum"])
    prop_md5_changed = prop_md5_check[~prop_md5_check.md5sum].index.values
    all_md5_changed = set(canon_md5_changed).union(set(prop_md5_changed))
    if len(all_md5_changed):
        eprint(
            "INT Notice: there are {} accessions from prevgc file where the sequence is different than what it was when originally added: {}. Dropping them from prevgc file.".format(
                len(all_md5_changed), all_md5_changed
            )
        )

    # combine 1&2 with 3:
    prevgc_notfound_or_changed = pd.concat(
        [
            prevgc_notfound,
            prevgc_df[prevgc_df["oldcanon"].isin(canon_md5_changed)],
            prevgc_df[prevgc_df["replacement"].isin(prop_md5_changed)],
        ]
    )

    # remove those which didn't pass condition 3:
    prevgc_df = prevgc_df[~prevgc_df["oldcanon"].isin(canon_md5_changed)]
    prevgc_df = prevgc_df[~prevgc_df["replacement"].isin(prop_md5_changed)]

    # dump to file those we will remove
    if len(prevgc_notfound) + prevgc_notfound_count + len(all_md5_changed):
        eprint(
            "INT All dropped lines from prevgc file are dumped to {}".format(
                config["prevgc_notfound_file"]
            )
        )
        prevgc_notfound_or_changed.to_csv(
            config["prevgc_notfound_file"],
            index=False,
            columns=[
                "pantherid",
                "org",
                "canon_acc",
                "canon_md5",
                "can_groupid",
                "prop_acc",
                "prop_md5",
                "prop_groupid",
                "release",
                "dataset",
            ],
            compression="gzip",
        )

    eprint("INT prevgc now contains {} suggestions".format(len(prevgc_df)))

    prevgc_df[
        "conflict"
    ] = False  # we init saying all is good; we'll mark it True if it has to be removed from output

    return prevgc_df


def dump_prev_changes(fh, prevgc_df, config=None):
    conflicts_file = config.get("conflicts_file", None)
    if conflicts_file is not None:
        if os.path.isfile(conflicts_file):
            # read in the conflict pairs written by multithread workers and collated at the end of the parallel work
            conflicts_df = pd.read_csv(
                conflicts_file,
                sep="\t",
                header=0,
                skip_blank_lines=True,
                names=["oldcanon", "replacement"],
            )
            if len(conflicts_df):
                # mark the conflicts in prevgc_df
                conflicts_df.set_index(["oldcanon", "replacement"], inplace=True)
                prevgc_df.set_index(["oldcanon", "replacement"], inplace=True)
                prevgc_match = prevgc_df.loc[conflicts_df.index]
                if len(conflicts_df) != len(prevgc_match):
                    eprint(
                        "    => ERROR: For some reason we could not find in prevgc all the conflict pairs listed in '{}'. Only found {} out of {}.".format(
                            conflicts_file, len(prevgc_match), len(conflicts_df)
                        )
                    )
                prevgc_df.loc[conflicts_df.index, "conflict"] = True
                prevgc_df.reset_index(inplace=True)
            else:
                eprint(
                    "    => ERROR: There are no conflict pairs specified in file '{}'".format(
                        conflicts_file
                    )
                )
        else:
            eprint(
                "    => NOTICE: Cannot mark multithread conflicts as specified conflicts file '{}' is not available for reading. No conflicts??".format(
                    conflicts_file
                )
            )
    prevgc_conflicts = prevgc_df[prevgc_df["conflict"]]
    if len(prevgc_conflicts):
        eprint(
            "INT Notice: there are {} suggestions from prevgc file that conflict with new suggestions and will be removed".format(
                len(prevgc_conflicts)
            )
        )
        eprint(
            "INT All removed lines from prevgc file are dumped to {}".format(
                config["prevgc_conflict_file"]
            )
        )
        prevgc_conflicts.to_csv(
            config["prevgc_conflict_file"],
            index=False,
            columns=[
                "pantherid",
                "org",
                "canon_acc",
                "can_groupid",
                "prop_acc",
                "prop_groupid",
                "release",
                "dataset",
            ],
            compression="gzip",
        )
    prevgc_ok = prevgc_df[~prevgc_df["conflict"]].sort_values(
        by=["release", "pantherid", "org", "canon_acc"]
    )  # sort before dump
    eprint("INT dumping {} prev suggestions as gc output".format(len(prevgc_ok)))
    prevgc_ok.to_csv(
        fh, index=False, header=False, sep="\t", columns=prev_changes_columns
    )  # dump with same format as original file


def simulate_changes(ortho_df, config=None):
    """
    To test result of integration into gene-centric pipeline,
    the simulate option (specifying a sugg_file, previous suggestions file)
    will read in the specified file and apply all the changes therein found,
    i.e. change the canonicals from the canon_acc to the prop_acc as per the file.
    """
    if config["skip_reprocessing_orthogroups"]:
        eprint(
            "    => ERROR: simulation prev_suggestions file has been specified but we are not re-processing groups, so this is meaningless, ignored!"
        )
        sys.exit(22)  # invalid
    if os.path.isfile(config["sugg_file"]):
        eprint(
            "SIMULATION: Reading prev_suggestions file '{}'".format(config["sugg_file"])
        )
        prev_suggestions_df = pd.read_csv(
            config["sugg_file"],
            sep="\t",
            header=None,
            comment="#",
            skip_blank_lines=True,
            names=[
                "pantherid",
                "org",
                "canon_acc",
                "canon_len",
                "prop_acc",
                "prop_len",
                "rank_score",
                "canon_cost",
                "prop_cost",
                "n_sp",
                "n_tax",
                "n_canon",
                "wn_canon",
                "scaled_prop_f",
                "scaled_p_cost",
                "p_len_diff",
                "clade",
                "clade_members",
                "MANEstatus",
            ],
            usecols=["pantherid", "org", "canon_acc", "prop_acc"],
        )
        prev_suggestions_df[["oldcanon_type", "oldcanon"]] = prev_suggestions_df[
            "canon_acc"
        ].str.split("|", expand=True)
        prev_suggestions_df[["replacement_type", "replacement"]] = prev_suggestions_df[
            "prop_acc"
        ].str.split("|", expand=True)
        prev_suggestions_df.drop(columns=["canon_acc", "prop_acc"], inplace=True)

        # minimal consistency check: all the oldcanon we are changing should actually be canonicals and vice versa
        oldcanon_check = ortho_df.loc[
            prev_suggestions_df["oldcanon"]
        ].is_canonical.unique()
        if oldcanon_check != [True]:
            eprint(
                "    => ERROR: Not all the canon_acc from the prev_suggestions file are actually canonicals in the ortho_df! Please check!"
            )
            sys.exit(33)  # numerical error
        replacement_check = ortho_df.loc[
            prev_suggestions_df["replacement"]
        ].is_canonical.unique()
        if replacement_check != [False]:
            eprint(
                "    => ERROR: Not all the prop_acc from the prev_suggestions file are actually isoforms in the ortho_df! Please check!"
            )
            sys.exit(33)  # numerical error

        # further consistency check: replacement and oldcanon should be in same genecentric group
        prev_suggestions_df.set_index("oldcanon", inplace=True)
        prev_suggestions_df["can_groupid"] = ortho_df.loc[
            prev_suggestions_df.index, "groupid"
        ]
        prev_suggestions_df.reset_index(inplace=True)
        prev_suggestions_df.set_index("replacement", inplace=True)
        prev_suggestions_df["prop_groupid"] = ortho_df.loc[
            prev_suggestions_df.index, "groupid"
        ]
        prev_suggestions_df.reset_index(inplace=True)
        groupid_check = (
            prev_suggestions_df.can_groupid == prev_suggestions_df.prop_groupid
        )
        if groupid_check.unique() != [True]:
            eprint(
                "    => ERROR: Not all the canon_acc/prop_acc from the prev_suggestions file are in the same genecentric group!"
            )
            eprint("    => CHECK the following inconsistencies")
            eprint(
                prev_suggestions_df.loc[groupid_check[~groupid_check].index][
                    [
                        "pantherid",
                        "org",
                        "oldcanon",
                        "can_groupid",
                        "replacement",
                        "prop_groupid",
                    ]
                ].to_csv(index=False, sep="\t")
            )
            sys.exit(33)  # numerical error

        # apply the suggestions
        ortho_df.loc[prev_suggestions_df["oldcanon"], "is_canonical"] = False
        ortho_df.loc[prev_suggestions_df["replacement"], "is_canonical"] = True
        eprint(
            "SIMULATION: Applied {} prev_suggestions to the dataframe".format(
                len(ortho_df.loc[prev_suggestions_df["oldcanon"]].index)
            )
        )
    else:
        eprint(
            "    => ERROR: Cannot open specified prev_suggestions file '{}' for reading".format(
                config["sugg_file"]
            )
        )
