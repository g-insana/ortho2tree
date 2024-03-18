#!/usr/bin/env python
# coding: utf-8
"""
module providing functions for loading panther and genecentric data
"""
import os
import gzip
import shutil
import pandas as pd
from .o2t_utils import eprint, download_file
import urllib.request
from urllib.error import URLError
from contextlib import closing


GCURL = "https://SPECIFYGCURL"


def get_panther_seq_classification(organism, config=None):
    file_prefix = "{}{}_".format(config["panther_data_dir"], config["panther_version"])
    url_prefix = "ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/{}_".format(
        config["panther_version"]
    )
    tab_file_gz = file_prefix + organism + ".gz"

    if not os.path.isfile(
        tab_file_gz
    ):  # if we haven't yet downloaded it: cache for next time
        try:
            with closing(urllib.request.urlopen(url_prefix + organism)) as r:
                with gzip.open(tab_file_gz, mode="wb", compresslevel=9) as f_out:
                    shutil.copyfileobj(r, f_out)
        except URLError as e:
            if "ftp error: error_perm" in str(e.reason):
                eprint(
                    "File not found or permission issue when fetching {}: {}".format(
                        url_prefix + organism, e.reason
                    )
                )
                return None
            else:
                eprint("URLError occurred: {}".format(e.reason))
                return None
        except Exception as e:
            eprint(f"An unexpected error occurred: {e}")
            return None
    return pd.read_csv(
        tab_file_gz, sep="\t", header=None, usecols=[0, 3], names=["gene", "pantherid"]
    )


def get_panther_df(config=None):
    tax2oscode = {}
    organisms = set(config["tax2org"].values())
    org2tax = {v: k for k, v in config["tax2org"].items()}
    df = pd.DataFrame(columns=["pantherid", "org", "acc"])
    for organism in sorted(organisms):
        eprint("Fetching mapping for {}".format(organism))
        org_df = get_panther_seq_classification(organism, config=config)
        if org_df is None:
            eprint(
                "    => ERROR: could not get panther mapping for {}".format(organism)
            )
            continue
        org_df = org_df.join(
            org_df["gene"]
            .str.split("|", expand=True)
            .rename(columns={0: "org", 1: "alt", 2: "acc"})
        )
        org_df["org"] = [x for x in org_df["org"]]
        taxid = org2tax[organism]
        oscode = org_df["org"][0]
        tax2oscode[taxid] = oscode
        del org_df["gene"]
        del org_df["alt"]  # unless needed..
        org_df["acc"] = [x[10:] for x in org_df["acc"]]
        eprint(
            "Loaded {} rows for {} aka {}, taxid {}".format(
                len(org_df), organism, oscode, taxid
            )
        )
        df = pd.concat([df, org_df], ignore_index=True)
    df.set_index("acc", inplace=True)
    eprint(" total number of entries from panther: {}".format(len(df)))

    if config["superfamily_level"]:
        df["lo_pantherid"] = [x.split(":", 1)[1] for x in df["pantherid"]]
        df["pantherid"] = [x.split(":", 1)[0] for x in df["pantherid"]]
    return df, tax2oscode


def get_gc_df(config={}, db_conn=None):
    GC_GROUPS_SQL = """
with /*+ PARALLEL */ mane_query as (
select distinct nvl(d2d.isoform, gce.accession) accession--, d2d.primary_id, d2d.secondary_id
from sptr.gene_centric_entry gce
inner join dbentry d on gce.accession=d.accession
inner join dbentry_2_database d2d on d.dbentry_id=d2d.dbentry_id
where release='{0}'
and gce.tax_id=9606 --only human has MANE
and database_id='MAS'
)
select /*+ PARALLEL */ distinct gce.tax_id, gce.accession, gce.entry_type, gce.group_id,
nvl(gce.is_canonical,0) is_canonical, s.length seqlength, s.md5 md5sum, nvl(geneid.primary_id,0) geneid,
is_fragment, CASE WHEN m.accession is NULL THEN 0 ELSE 1 END has_mane
from sptr.gene_centric_entry gce
inner join dbentry d on gce.accession=d.accession
inner join sequence s on d.dbentry_id=s.dbentry_id
left outer join dbentry_2_database partition("GeneID") geneid on d.dbentry_id=geneid.dbentry_id
left outer join mane_query m on m.accession=gce.accession
where release='{0}'
and gce.tax_id in ({1})
and exists (select 1 from sptr.gene_centric_entry gce3 where gce.group_id=gce3.group_id and release='{0}' and gce3.is_canonical=1) --needs a canonical
and not exists (select 1 from sptr.gene_centric_entry gce2 where gce.group_id=gce2.group_id and release='{0}' and gce.accession=gce2.accession and gce2.is_canonical=1 and gce.is_canonical is null)
--and (is_fragment=0 or gce.entry_type=0) --we get all including fragments
--order by gce.tax_id,gce.group_id,nvl(gce.is_canonical,0) desc,gce.accession --optional
""".format(
        config["up_release"], ",".join([str(x) for x in config["tax_ids"]])
    )

    gc_dtypes = {
        "tax_id": int,
        "acc": str,
        "trembl": bool,
        "groupid": int,
        "is_canonical": bool,
        "seqlen": int,
        "md5sum": str,
        "geneid": int,
        "is_fragment": bool,
        "has_mane": bool,
    }
    gc_columns = [
        "tax_id",
        "acc",
        "trembl",
        "groupid",
        "is_canonical",
        "seqlen",
        "md5sum",
        "geneid",
        "is_fragment",
        "has_mane",
    ]

    genecentric_tabfile = "genecentric_groups.{}{}.all.out.gz".format(
        config["up_release"], config["dataset_maindir"]
    )
    if os.path.exists(
        os.path.join(config["dataset_maindir"], genecentric_tabfile)
    ):  # for faster execution, we cache it
        eprint("Loading genecentric data from file {}".format(genecentric_tabfile))

        gc_df = pd.read_csv(
            os.path.join(config["dataset_maindir"], genecentric_tabfile),
            header=0,
            names=gc_columns,
            dtype=gc_dtypes,
        )
    elif config["gc_from_sql"]:
        eprint("executing GC_GROUPS_SQL query")
        gc_df = pd.read_sql_query(GC_GROUPS_SQL, con=db_conn)
        gc_df.columns = gc_columns
        gc_df = gc_df.astype(dtype=gc_dtypes)
        eprint("saving gc_data as csv file")
        gc_df.to_csv(
            os.path.join(config["dataset_maindir"], genecentric_tabfile),
            index=False,
            compression="gzip",
        )  # cache for next time
    else:  # download pre-generated file
        url = "{}{}".format(GCURL, genecentric_tabfile)
        eprint("Downloading genecentric data from {}".format(url))
        download_file(url, os.path.join(config["dataset_maindir"], genecentric_tabfile))
        gc_df = pd.read_csv(
            os.path.join(config["dataset_maindir"], genecentric_tabfile),
            sep="\t",
            header=None,
            names=gc_columns,
            dtype=gc_dtypes,
        )

    gc_df.drop_duplicates(
        inplace=True,
        subset=["acc", "groupid", "is_canonical", "seqlen", "md5sum", "trembl"],
        keep="first",
    )  # remove multiple geneid, if any
    gc_df = gc_df.sort_values(by="is_canonical", ascending=False).drop_duplicates(
        subset=["groupid", "md5sum"], keep="first"
    )  # remove multiple identical sequences in same gcgroup, keeping the canonical
    gc_df.set_index("acc", inplace=True)
    gc_df = gc_df.loc[gc_df["tax_id"].isin(config["tax_ids"])]  # drop any extra taxa

    # oscode from tax_id
    gc_df["org"] = [config["tax2oscode"][x] for x in gc_df["tax_id"]]
    gc_df.drop(columns=["tax_id"], inplace=True)

    # mark entry_type as sp or tr
    gc_df["entry_type"] = gc_df["trembl"].map({True: "tr", False: "sp"})
    gc_df.drop(columns=["trembl"], inplace=True)
    return gc_df
