#!/usr/bin/env python
# coding: utf-8
"""
module providing functions needed for refseq integration
"""
import os
import pandas as pd
from o2t_utils import eprint


def integrate_missing_geneid(gc_df, config=None):
    """
    integrate missing geneid information from a provided mapping file
    """
    # find number of missing geneid accessions
    missing_geneid = gc_df.loc[(gc_df["geneid"] == 0)]
    eprint(
        "Note: we have {} accessions without geneid information and a mapping file is present".format(
            len(missing_geneid)
        )
    )

    gc_df.loc[missing_geneid.index, "geneid"] = None  # set to NaN to use fillna

    # load mapping file
    missing_geneid_map = pd.read_csv(
        os.path.join(config["dataset_maindir"], config["missing_geneid_mapfile"]),
        header=0,
        sep="\t",
        names=["taxon_id", "acc", "geneid"],
        dtype={"taxon_id": int, "acc": str, "geneid": int},
    )
    if missing_geneid_map.duplicated("acc").any():
        eprint(
            "WARNING! {} duplicated acc in mapping file! Dropping duplicates".format(
                missing_geneid_map.duplicated("acc").sum()
            )
        )
        missing_geneid_map.drop_duplicates("acc", inplace=True)  # in case there are...
    missing_geneid_map.set_index("acc", inplace=True)

    # use mapping file to fill in missing information
    eprint(
        "Filling the missing geneid from mapping file '{}', which contains {} mappings".format(
            config["missing_geneid_mapfile"], len(missing_geneid_map)
        )
    )
    gc_df.geneid.fillna(missing_geneid_map.geneid, inplace=True)
    gc_df.geneid.fillna(0, inplace=True)  # back to 0 those we could not fill
    gc_df["geneid"] = gc_df["geneid"].astype(int)  # back to integer

    # compare results
    eprint(
        "We now have {} accessions without geneid information".format(
            len(gc_df.loc[(gc_df["geneid"] == 0)])
        )
    )


def integrate_refseq(ortho_df, config=None, get_sequences_fn=None):
    """
    to add refseq entries to the ortho_df based on geneid2refseq mapping
    """
    g2r_df = None
    refseq_df = pd.DataFrame()
    if os.path.isfile(config["geneid2refseq_mapfile"]):
        eprint(
            "processing geneid2refseq file {}".format(config["geneid2refseq_mapfile"])
        )
        g2r_df = pd.read_csv(
            config["geneid2refseq_mapfile"],
            sep="\t",
            header=None,
            names=["tax_id", "geneid", "transcript", "protein"],
            dtype={"tax_id": int, "geneid": int, "transcript": str, "protein": str},
            usecols=["tax_id", "geneid", "protein"],
        )
        g2r_df = g2r_df.loc[
            g2r_df["tax_id"].isin(config["tax_ids"])
        ]  # dropping taxa not in the analysis
        g2r_df = g2r_df.loc[g2r_df["protein"] != "-"]  # keeping only cds
        g2r_df.set_index("geneid", inplace=True)
        eprint(
            " Refseq mappings loaded by tax_id:",
            g2r_df.groupby("tax_id")["protein"].describe()["count"],
        )
    else:
        eprint(" NOTICE: no geneid2refseq mapping file found")

    if g2r_df is not None:  # if we have a mapping between geneid and refseq
        has_geneid = ortho_df.loc[
            (ortho_df["geneid"] != 0)
        ]  # subset of ortho_df with geneid information
        refseq_df = pd.merge(
            has_geneid,
            g2r_df,  # table join
            how="inner",
            on="geneid",
            sort=False,
            suffixes=("_x", "_y"),
            copy=False,
            indicator=False,
        ).set_index("protein")

        refseq_df.drop(columns=["tax_id"], inplace=True)  # no need
        refseq_df["entry_type"] = "ref"  # refseq
        refseq_df["is_canonical"] = False
        refseq_df["md5sum"] = None
        refseq_df["seqlen"] = None

        # keep only one refseq in case of multiple (e.g. XP_003640289 would be duplicated in two gc groups 4674006 and 4681279):"
        refseq_df = (
            refseq_df.reset_index()
            .sort_values(by=["org", "geneid", "is_canonical", "pantherid"])
            .drop_duplicates(subset=["protein", "org", "geneid"], keep="last")
            .set_index("protein")
        )
        # (by sorting first and keeping 'last' we preferentially keep the one associated with a canonical)

        # fetching and filling information:
        eprint(
            "Fetching sequence information for {} refseq entries".format(len(refseq_df))
        )
        refseq_seqinfo = get_sequences_fn(
            refseq_df.index.unique(), format="DataFrameNoSeq"
        )  # retrieve information of those sequences
        refseq_df.seqlen.fillna(
            refseq_seqinfo.seqlen, inplace=True
        )  # fill in seqlen info
        refseq_df.md5sum.fillna(
            refseq_seqinfo.md5sum, inplace=True
        )  # fill in md5sum info

        # drop refseq that we could not retrieve
        no_sequence_data_accs = refseq_df[refseq_df["seqlen"].isnull()].index.to_list()
        eprint(
            "We could not retrieve sequence information for {} refseq accessions: {}...".format(
                len(no_sequence_data_accs),
                no_sequence_data_accs[0 : min(10, len(no_sequence_data_accs))],
            )
        )
        refseq_df.drop(no_sequence_data_accs, inplace=True)

        # drop refseq duplicated sequences (same md5 and same geneid):
        refseq_df.drop_duplicates(subset=["geneid", "md5sum"], inplace=True)
        refseq_count_before = refseq_df.shape[0]

        # drop entries that we already have in uniprot (by md5sum being equal):
        refseq_df = refseq_df[~(refseq_df.md5sum.isin(ortho_df.md5sum.unique()))]
        refseq_count_after = refseq_df.shape[0]

        refseq_df["seqlen"] = refseq_df["seqlen"].astype(int)  # seqlen back to int

        refseq_df.index.rename(
            "acc", inplace=True
        )  # rename the index from protein to acc
        eprint(
            "Adding {} new refseq sequences to the analysis (excluding {} entries already in UP)".format(
                refseq_count_after, refseq_count_before - refseq_count_after
            )
        )

    return refseq_df
