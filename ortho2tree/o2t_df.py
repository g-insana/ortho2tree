#!/usr/bin/env python
# coding: utf-8
"""
module providing functions
* create_ortho_df(), ex-novo creation of ortho_df
* read_ortho_df(), reading previously created ortho_df from file
* ortho_df_stats(), calculating stats and filtering low taxa groups
"""
import os
import sys
import pandas as pd
from .o2t_utils import eprint
from .o2t_gc_integration import read_prev_changes, simulate_changes
from .o2t_refseq import integrate_missing_geneid, integrate_refseq
from .o2t_outliers import flag_maxseqlen, flag_outliers_in_df

sys.stderr.write("INFO: Panda version: {}\n".format(pd.__version__))


def read_ortho_df(config=None):
    """
    read a previously created dataframe from file
    return ortho_df and prevgc_df dataframes
    """
    eprint(
        "\nReading orthogroup data from cached file {}".format(
            config["orthogroup_df_cachefile"]
        )
    )
    ortho_df = pd.read_csv(config["orthogroup_df_cachefile"], index_col="acc")

    if config.get("prevgc_file", False):
        # read in previous changes (for integration in final cumulative output_gc)
        prevgc_df = read_prev_changes(ortho_df, config=config)
    else:
        prevgc_df = pd.DataFrame()

    # if we wanted to see the trembl fragments:
    # tr_fragments_df = pd.read_csv(config['tr_fragments_df_cachefile'], index_col='acc')
    # eprint("Also loaded {} trembl fragments".format(len(tr_fragments_df)))

    eprint(
        "\nThere are {} accessions flagged as outliers".format(
            ortho_df["outlier"].sum()
        )
    )

    return ortho_df, prevgc_df


def create_ortho_df(panther_df, gc_df, config=None, get_sequences_fn=None):
    """
    create ortholog groups mapping panther and genecentric data
    return ortho_df and prevgc_df dataframes
    """
    if config["add_refseq"] and os.path.exists(
        os.path.join(config["dataset_maindir"], config["missing_geneid_mapfile"])
    ):
        # optional additional geneid mapping for refseq integration
        integrate_missing_geneid(gc_df, config=config)

    ortho_df = pd.merge(
        panther_df,
        gc_df,
        how="outer",  # full outer join: we keep accessions only found in panther and not in gc
        left_index=True,
        right_index=True,
        sort=False,
        suffixes=("_x", "_y"),
        copy=False,
        indicator=False,
        validate="1:1",
    )

    if len(
        ortho_df[
            (ortho_df["org_x"].notnull())
            & (ortho_df["org_y"].notnull())
            & (ortho_df["org_x"] != ortho_df["org_y"])
        ]
    ):
        eprint(
            "ERROR! Conflicting organism data between panther and genecentric!!",
            ortho_df[
                (ortho_df["org_x"].notnull())
                & (ortho_df["org_y"].notnull())
                & (ortho_df["org_x"] != ortho_df["org_y"])
            ],
        )
        sys.exit(1)

    ortho_df["org_x"] = ortho_df["org_x"].fillna(ortho_df["org_y"])
    ortho_df.drop(columns=["org_y"], inplace=True)
    ortho_df.rename(columns={"org_x": "org"}, inplace=True)

    if len(ortho_df[(ortho_df["org"].isnull())]):
        eprint(
            "ERROR! Accession with missing organism???"
        )  # xkcd-2200 unreachable state
        sys.exit(1)

    not_in_gc = ortho_df[ortho_df["groupid"].isnull()].index.to_list()
    if len(not_in_gc):  # from 2023_03
        eprint(
            "There are {} accessions in panther data but not in genecentric.. removing them".format(
                len(not_in_gc)
            )
        )
        ortho_df.drop(not_in_gc, inplace=True)

    # set defaults and convert types changed after merge
    ortho_df["groupid"] = ortho_df["groupid"].astype(int)
    ortho_df["seqlen"] = ortho_df["seqlen"].astype(int)
    ortho_df["geneid"] = ortho_df["geneid"].astype(int)

    # find isoforms in panther:
    iso_in_panther = ortho_df[
        (ortho_df["pantherid"].notnull())
        & (ortho_df["groupid"].notnull())
        & (~(ortho_df["is_canonical"]))
    ].copy(deep=False)
    if len(iso_in_panther):
        eprint(
            "There are {} genecentric accessions in panther data; e.g.: ".format(
                len(iso_in_panther)
            )
        )
        eprint(iso_in_panther.head().index.to_list())

    # gc integration
    prevgc_df = pd.DataFrame()
    if config.get("prevgc_file", False):
        if os.path.isfile(config["prevgc_file"]):
            # read in previous changes (for integration in final cumulative output_gc)
            prevgc_df = read_prev_changes(ortho_df, config=config)
        else:
            eprint(
                "    => ERROR: Cannot open specified prev_suggestions file '{}' for reading".format(
                    config["prevgc_file"]
                )
            )
            config["prevgc_file"] = False

    # simulation
    if config.get("sugg_file", False):
        # optional simulation of applying changes (previous suggestions)
        simulate_changes(ortho_df, config=config)

    # refseq integration
    if config["add_refseq"]:
        # optionally add (novel, not in UP) refseq entries to the dataframe via geneid mapping
        refseq_df = integrate_refseq(
            ortho_df, config=config, get_sequences_fn=get_sequences_fn
        )
        # add the good refseq to the main dataframe
        if not refseq_df.empty:
            ortho_df = pd.concat(
                [ortho_df, refseq_df]
            )  # append the new refseq entries to ortho_df

    # now we'll add pantherid to all isoforms based on their genecentric group
    # and panther information known for that group

    # prepare mapping between groupid and pantherid and use it to fill based on groupid
    groupid2pantherid = ortho_df.loc[
        (ortho_df["pantherid"].notnull()) & (ortho_df["groupid"].notnull()),
        ["pantherid", "groupid"],
    ].drop_duplicates()
    groupid2pantherid.set_index("groupid", inplace=True)
    if -1 in groupid2pantherid.index:
        # drop where missing group
        groupid2pantherid.drop(-1, inplace=True)

    missing_pantherid_count = len(ortho_df[ortho_df.pantherid.isnull()])
    eprint(
        "\nGenecentric accessions without pantherid: {}".format(missing_pantherid_count)
    )
    assignable_missing_pantherid = ortho_df[
        ortho_df.pantherid.isnull() & ortho_df.groupid.isin(groupid2pantherid.index)
    ]
    assignable_pantherid_count = len(assignable_missing_pantherid)
    ortho_df.loc[assignable_missing_pantherid.index, "pantherid"] = ortho_df.loc[
        assignable_missing_pantherid.index, "groupid"
    ].map(groupid2pantherid.to_dict()["pantherid"], na_action="ignore")

    eprint(
        "Assigned pantherid to {} genecentric accessions. {} cannot be mapped".format(
            assignable_pantherid_count,
            missing_pantherid_count - assignable_pantherid_count,
        )
    )
    # unassignable_accs = ortho_df[ortho_df.pantherid.isnull()].index.to_list() #if we'd like to look at these

    if config["superfamily_level"]:
        # we keep the lower SF subdivision as a separate lo_pantherid column
        groupid2lo_pantherid = ortho_df.loc[
            (ortho_df["lo_pantherid"].notnull()) & (ortho_df["groupid"].notnull()),
            ["lo_pantherid", "groupid"],
        ].drop_duplicates()
        groupid2lo_pantherid.set_index("groupid", inplace=True)
        if -1 in groupid2lo_pantherid.index:
            # drop where missing group
            groupid2lo_pantherid.drop(-1, inplace=True)
        ortho_df.loc[assignable_missing_pantherid.index, "lo_pantherid"] = ortho_df.loc[
            assignable_missing_pantherid.index, "groupid"
        ].map(groupid2lo_pantherid.to_dict()["lo_pantherid"], na_action="ignore")

    # cleanup of dataframe
    # find accessions in panther which are isoform in current gc (they were canonical when panther was last updated):
    # these isoforms in panther could create orphan iso issue: they could split genecentric groups in different panther groups
    canonicals_of_panther_isoforms = (
        ortho_df[
            (ortho_df.groupid.isin(ortho_df.loc[iso_in_panther.index]["groupid"]))
            & ortho_df.is_canonical
        ][["pantherid", "groupid"]]
        .reset_index()
        .set_index("groupid")
    )
    canonicals_of_panther_isoforms.rename(
        columns={"pantherid": "can_pantherid"}, inplace=True
    )
    iso_in_panther.reset_index(inplace=True)
    iso_in_panther.set_index("groupid", inplace=True)
    iso_in_panther["can"] = canonicals_of_panther_isoforms["acc"]
    iso_in_panther["can_pantherid"] = canonicals_of_panther_isoforms["can_pantherid"]
    iso_in_panther.set_index("acc", inplace=True)
    orphan_isoforms = iso_in_panther[
        (iso_in_panther["pantherid"] != iso_in_panther["can_pantherid"])
    ]
    eprint(
        "There are {} panther orphan isoforms: with pantherid different than that of their canonical. Rectifying.".format(
            len(orphan_isoforms)
        )
    )
    ortho_df.loc[orphan_isoforms.index, "pantherid"] = orphan_isoforms["can_pantherid"]

    # how many canonicals now missing pantherid mapping?
    unmapped_canonicals = ortho_df[
        (ortho_df["is_canonical"]) & (ortho_df["pantherid"].isnull())
    ]
    if len(unmapped_canonicals):
        eprint(
            "There are {} unmapped canonicals. E.g.:".format(len(unmapped_canonicals))
        )
        eprint(unmapped_canonicals.head().index.to_list())
        eprint(
            "Removing them from df but writing them to file '{}'".format(
                config["unmapped_canonicals_file"]
            )
        )
        with open(
            config["unmapped_canonicals_file"], "w", encoding="utf-8"
        ) as outputfh:
            for can, groupid in unmapped_canonicals["groupid"].to_dict().items():
                outputfh.write("{}\t{}\n".format(can, groupid))
        ortho_df.dropna(
            subset=["pantherid"], inplace=True
        )  # remove unmapped entries from dataframe

    # are there groups without a single canonical?
    groups_without_canonicals = ortho_df.groupby("pantherid")["is_canonical"].any()
    groups_without_canonicals = groups_without_canonicals[
        ~groups_without_canonicals
    ].index
    if len(groups_without_canonicals):
        eprint(
            "There are {} groups without a single canonical. E.g.: {}".format(
                len(groups_without_canonicals),
                groups_without_canonicals.to_list()[0:10],
            )
        )
        eprint("Removing them")
        ortho_df.drop(
            ortho_df[ortho_df["pantherid"].isin(groups_without_canonicals)].index,
            inplace=True,
        )

    # check no missing values in the dataframe
    if ortho_df.isnull().any().any():
        eprint("Error, we have a column with missing values:")
        eprint(ortho_df.isnull().any())
        sys.exit(34)

    # dropping trembl fragments from the df
    tr_fragments_df = ortho_df[
        (ortho_df["entry_type"] == "tr") & (ortho_df["is_fragment"] == 1)
    ]
    eprint(
        "\nDropping {} trembl fragments from the dataframe".format(len(tr_fragments_df))
    )
    ortho_df.drop(tr_fragments_df.index, inplace=True)
    tr_fragments_df.to_csv(
        config["tr_fragments_df_dumpfile"],
        index=True,
        index_label="acc",
        compression="gzip",
    )
    eprint("TR fragments dumped as {}".format(config["tr_fragments_df_dumpfile"]))

    # identify and mark outliers
    eprint("\nMarking outliers, please wait..")

    flag_outliers_in_df(ortho_df, config)

    # set all canonicals and original panthers to non-outlier, regardless of what just computed above
    panther_potential_outliers = ortho_df.loc[
        ortho_df.index.intersection(panther_df.index).intersection(
            ortho_df[ortho_df["outlier"]].index
        )
    ].index
    if len(panther_potential_outliers):
        eprint(
            "Avoiding setting as outliers {} acc coming from panther data, e.g. {}".format(
                len(panther_potential_outliers),
                panther_potential_outliers.to_list()[0:10],
            )
        )
        ortho_df.loc[panther_potential_outliers, "outlier"] = False
    # gc_can_potential_outliers = ortho_df[(ortho_df['is_canonical']) & (ortho_df['outlier'])].index
    # if len(gc_can_potential_outliers):
    #    print("Avoiding setting as outliers {} gene centric canonicals {}".format(len(gc_can_potential_outliers), gc_can_potential_outliers.to_list()))
    #    ortho_df.loc[gc_can_potential_outliers, 'outlier'] = False

    eprint(
        "\nIdentified {} outliers based on sequence length".format(
            len(ortho_df[ortho_df["outlier"]])
        )
    )

    # remove all those beyond max_seqlen (too long isoforms plus too long canonicals with their isoforms):
    flag_maxseqlen(ortho_df, config["max_seqlen"])

    eprint(
        "\nThere are now {} accessions flagged as outliers".format(
            ortho_df["outlier"].sum()
        )
    )
    if config["detect_outliers_with_median_and_can_lengths"]:
        eprint(
            "{} flagged due to distance from median, {} flagged due to distance from canonicals max/min seqlengths".format(
                ortho_df["outlier_by_median"].sum(),
                ortho_df["outlier_by_can_len"].sum(),
            )
        )
    elif config["detect_outliers_with_quart_and_can_lengths"]:
        eprint(
            "{} flagged due to distance from q1/q3, {} flagged due to distance from canonicals max/min seqlengths".format(
                ortho_df["outlier_by_quart"].sum(), ortho_df["outlier_by_can_len"].sum()
            )
        )

    # create a dump of the orthogroup data
    if config["dump_orthogroup_data"]:
        ortho_df.to_csv(
            config["orthogroup_df_dumpfile"],
            index=True,
            index_label="acc",
            compression="gzip",
        )
        eprint("\nDataframe dumped as {}".format(config["orthogroup_df_dumpfile"]))

    return ortho_df, prevgc_df


def ortho_df_stats(ortho_df, config=None):
    """
    Print a series of stats and filter out groups lower than min_taxa
    """
    if config.get("print_stats", True):
        # statistics on sequence length per species
        eprint(
            "\nStats on seqlen per species:\n",
            ortho_df.groupby("org")["seqlen"].describe()[
                ["count", "min", "max", "std", "mean"]
            ],
        )

        # statistics on outlier per species
        eprint(
            "\nStats on outliers per species:\n",
            ortho_df.groupby(["org", "outlier"])["seqlen"].describe()["count"],
        )

        # statistics on number of different taxa per group
        eprint(
            "\nStats on number of different taxa per group:\n",
            ortho_df.groupby(["pantherid"])["org"].nunique().describe(),
        )

        # list of groups by number of different taxa:
        eprint(
            "\nPrinted list of groups and number of different taxa to file '{}', e.g.:".format(
                config["groups_by_taxa_count_file"]
            )
        )
        ortho_df.groupby(["pantherid"])["org"].nunique().sort_values(
            ascending=False
        ).to_csv(config["groups_by_taxa_count_file"], header=["n_taxa"])
        eprint(
            ortho_df.groupby(["pantherid"])["org"]
            .nunique()
            .sort_values(ascending=False)
            .head()
        )

    # remove groups lower than min_taxa
    orthogroups = set(ortho_df["pantherid"].drop_duplicates().values)
    all_groups_count = len(
        orthogroups
    )  # before removing the ones with n_taxa < min_taxa_threshold
    low_taxa_groups = (
        ortho_df[
            ortho_df.groupby("pantherid")["org"].transform("nunique")
            < config["min_taxa_threshold"]
        ]["pantherid"]
        .drop_duplicates()
        .values
    )  # if we want to keep track of them
    ortho_df = ortho_df[
        ortho_df.groupby("pantherid")["org"].transform("nunique")
        >= config["min_taxa_threshold"]
    ]
    orthogroups = set(ortho_df["pantherid"].drop_duplicates().values)
    all_canonicals = set(ortho_df[ortho_df["is_canonical"]].index.to_list())
    eprint(
        "\nRemoved {} groups which had number of different taxa lower than threshold ({}). Total number of groups now: {}, total number of canonicals: {}, total number of accessions: {}".format(
            all_groups_count - len(orthogroups),
            config["min_taxa_threshold"],
            len(orthogroups),
            len(all_canonicals),
            len(ortho_df),
        )
    )

    # write list of low taxa groups
    if len(low_taxa_groups):
        with open(config["low_taxa_groups_file"], "w", encoding="utf-8") as outputfh:
            for orthoid in low_taxa_groups:
                outputfh.write("{}\n".format(orthoid))
        eprint(
            "Written list of low taxa groups to file '{}'".format(
                config["low_taxa_groups_file"]
            )
        )

    if config.get("print_stats", True):
        # statistics on outlier per species
        eprint(
            "\nStats on sequences and outliers per species after min_taxa threshold filter:\n",
            ortho_df.groupby("org")["seqlen"].describe()["count"],
            "\n\n",
            ortho_df.groupby(["org", "outlier"])["seqlen"].describe()["count"],
        )
        eprint(
            "\nTotal sequences, by status: ",
            ortho_df.groupby("outlier")["pantherid"].count(),
        )
        # list of groups by number of entries (excluding lowtaxa and outliers)
        eprint(
            "\nPrinted list of groups and number of entries to file '{}'".format(
                config["groups_by_entry_count_file"]
            )
        )
        pd.DataFrame(ortho_df[~ortho_df["outlier"]]["pantherid"].value_counts()).to_csv(
            config["groups_by_entry_count_file"]
        )

        # stats on "longest sequence"
        ortho_df_sp = ortho_df[ortho_df["entry_type"] == "sp"]
        ortho_df_onlytr = ortho_df[~ortho_df["groupid"].isin(ortho_df_sp["groupid"])]

        orthogroups_sp = ortho_df_sp.groupby("groupid")
        orthogroups_onlytr = ortho_df_onlytr.groupby("groupid")

        # calculate statistics for groups with sp (curated) entries
        a = pd.DataFrame(orthogroups_sp["seqlen"].max())
        a.columns = ["longest_sp_seqlen"]
        a["cano_seqlen"] = ortho_df_sp[ortho_df_sp["is_canonical"]].set_index(
            "groupid"
        )["seqlen"]
        a["cano_is_longest"] = a["cano_seqlen"] >= a["longest_sp_seqlen"]

        # calculate statistics for groups with only trembl (unreviewed) entries
        b = pd.DataFrame(orthogroups_onlytr["seqlen"].max())
        b.columns = ["longest_tr_seqlen"]
        b["cano_seqlen"] = ortho_df_onlytr[ortho_df_onlytr["is_canonical"]].set_index(
            "groupid"
        )["seqlen"]
        b["cano_is_longest"] = b["cano_seqlen"] >= b["longest_tr_seqlen"]

        eprint("\nNumber of total gc: {}".format(len(ortho_df.groupby("groupid"))))
        eprint("Number of gc which have sp entries: {}".format(len(orthogroups_sp)))
        eprint(
            "Number of gc where sp cano is the longest sp sequence in the group: {}".format(
                a["cano_is_longest"].sum()
            )
        )
        eprint("Number of gc with only tr entries: {}".format(len(orthogroups_onlytr)))
        eprint(
            "Number of gc where tr cano is the longest tr sequence in the group: {}".format(
                b["cano_is_longest"].sum()
            )
        )

    return ortho_df
