#!/usr/bin/env python
# coding: utf-8
"""
module providing several general purpose utility functions
"""
import os
import re
import sys
import time
import pandas as pd
from hashlib import md5  # for md5sum calculation
from requests import get  # for web retrieval

REFASTAHEADER_STRIP = re.compile(r"^>.*?\n")  # strip existing fasta header
REFASTAHEADER_NOID = re.compile(r"^>(\S*\|\S*\|)\S*\n")  # strip ID from UP fasta header


def eprint(*myargs, **kwargs):
    """
    print to stderr, useful for error messages and to not clobber stdout
    """
    print(*myargs, file=sys.stderr, **kwargs)


def batch(iterable, batchsize=1):
    """
    split iterable in constant-size chunks
    """
    length = len(iterable)
    for index in range(0, length, batchsize):
        yield iterable[index : min(index + batchsize, length)]


def download_file(url, filename):
    """
    download file. parameters: url, filename
    """
    with get(url, stream=True) as r:
        r.raise_for_status()
        with open(filename, "wb") as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    return filename


def print_subfiles(filenames, fh=None, offset=0):
    """
    read and print a series of temporary sub-files writing them sequentially
    to a given output filehandle or print them to stdout
    """
    BUFFERSIZE = 1048576
    mode = "r"
    if fh is None:
        fh = sys.stdout  # print to stdout if no fh given
        mode = "r"
    for filename in filenames:
        with open(filename, mode) as infile:
            ##straight, unbuffered, will read whole file in memory
            # fh.write(infile.read())
            ##buffered:
            filesize = os.path.getsize(filename)
            if offset > 0:
                infile.seek(offset)
                filesize -= offset
            blockcount, remainder = divmod(filesize, BUFFERSIZE)
            for _ in range(blockcount):
                buffered = infile.read(BUFFERSIZE)
                fh.write(buffered)
            fh.write(infile.read(remainder))


def delete_files(filenames, path=None):
    """
    delete the temporary sub-files
    optionally specify a path where to look for the files
    """
    for filename in filenames:
        if path is not None:
            filename = os.path.join(path, filename)
        if os.path.isfile(filename):
            os.remove(filename)


def list_files_with_extension(path, extension):
    """
    Lists files in a directory with the specified extension.
    """
    return set(
        filename for filename in os.listdir(path) if filename.endswith(extension)
    )


def check_all_files_exist(path, names, extension):
    """
    Ensure that all the files specified (as name + extension for name in names)
    are present under path
    """
    filenames_found = list_files_with_extension(path, extension)
    filenames_expected = set(name + "." + extension for name in names)

    filenames_missing = filenames_expected.difference(filenames_found)

    if len(filenames_missing):
        eprint(
            "WARNING: there are {} expected but missing files, e.g. {}".format(
                len(filenames_missing), list(filenames_missing)[0:10]
            )
        )
        return False

    return True


def secs2time(secs):
    """
    human readable printout of seconds elapsed
    e.g.:
    secs2time(3663)
    """
    minutes, seconds = divmod(secs, 60)
    hours, minutes = divmod(minutes, 60)
    return "{:02.0f}h {:02.0f}m {:02.0f}s".format(hours, minutes, seconds)


def elapsed_time(start_time, work_done=None):
    """
    compute elapsed time from given start_time (in seconds) and return
    well formatted string. If work_done argument given, compute speed of the process
    e.g.:
    start_secs = time.time(); iterations_done = 10; time.sleep(2); print(" '-- Elapsed: {}, {} it/s --'".format(*elapsed_time(start_secs, iterations_done)))
    print(" '-- Elapsed: {} --'".format(elapsed_time(start_secs)))
    """
    process_time = time.time() - start_time
    if work_done is None:
        return secs2time(process_time)
    process_speed = round(work_done / process_time, 2)
    return secs2time(process_time), process_speed


# Sequence utils
def sequences2fasta(sequences, organisms={}):
    """
    fasta will be created with sequences ordered by organism (in alphabetic order) and then alphabetically by ACC
    this for consistency of alignment
    e.g.:
    sequences2fasta(get_sequences(['Q8BZ60','A0A075B784','XP_011514185']))
    sequences2fasta(get_sequences(['G3UYE7', 'Q8N9L1','A0A075B784','XP_011514185','A0A2I2Y9W8']), organisms={'Q8N9L1':'HUMAN', 'A0A075B784': 'MOUSE', 'XP_011514185': 'HUMAN', 'A0A2I2Y9W8': 'GORGO', 'G3UYE7': 'MOUSE'})
    """
    org2acc = {}  # (keys: organisms, values: accs)
    if (
        organisms
    ):  # if we are given organisms, group sequences by org and sort them alphabetically inside each group
        for org in organisms.values():
            org2acc[org] = []
        for acc in sequences:
            org = organisms[acc]
            org2acc[org].append(acc)
        for org in org2acc.keys():
            org2acc[org] = sorted(
                org2acc[org]
            )  # for each organism, sequences sorted alphabetically by acc
    else:  # otherwise, simply sort sequences alphabetically
        sequences_acc = sorted(sequences.keys())
        organisms = {acc: "" for acc in sequences_acc}
        org2acc = {"": sequences_acc}
    fasta_txt = ""
    for org in sorted(org2acc.keys()):  # organisms in alphabetic order
        # eprint("printing for {}".format(org)) #debug
        for acc in org2acc[org]:
            # eprint("  {}".format(acc)) #debug
            if acc.find("_") != -1:  # refseq
                if org:
                    fasta_txt += ">ref_iso|{}|{}\n".format(acc, org.upper())
                else:
                    fasta_txt += ">ref_iso|{}|UNKNOWN\n".format(acc)
                fasta_txt += REFASTAHEADER_STRIP.sub("", sequences[acc])
            else:
                if org:
                    fasta_txt += REFASTAHEADER_NOID.sub(
                        r">\1{}\n".format(org.upper()), sequences[acc]
                    )
                else:
                    fasta_txt += sequences[acc]  # keep as it is

    return fasta_txt


def get_sequences_swpread(
    accessions, orthoid="", format="Dict", config={}, db_conn=None
):
    """
    input: a list of accessions and optionally an orthoid (for error msgs only) and a format
    returns: depending on the format parameter, the sequence information in that format
    formats: DataFrame, DataFrameNoSeq (same as DataFrame but without the sequence, only sequence info), Dict
    """
    MAXLIST = 900  # to avoid ORA-01795: maximum number of expressions in a list is 1000

    SEQUENCE_FULL = """
        select d.accession as acc,
        (case when d.entry_type in (0,13) then 'sp' else 'tr' end) as entry_type,
        s.length as seqlen,
        s.sequence as sequence,
        s.md5 as md5sum
        """

    SEQUENCE_INFO = """
        select d.accession as acc,
        (case when d.entry_type in (0,13) then 'sp' else 'tr' end) as entry_type,
        s.length as seqlen,
        s.md5 as md5sum
        """

    SEQUENCE_FASTA = """
        select d.accession as acc,
        '>'||(case when d.entry_type in (0,13) then 'sp'  else 'tr' end)||'|'||d.accession||'|'||d.name||'\n'||s.sequence||'\n' as sequence
        """

    QUERY_TAIL = """
    from dbentry d
    inner join sequence s on d.dbentry_id=s.dbentry_id
    where d.accession in ({})
    and d.deleted='N' and d.entry_type in (0, 1, 13) and d.merge_status<>'R'
    and d.first_public is not null
    """

    SEQUENCE_FULL += QUERY_TAIL
    SEQUENCE_INFO += QUERY_TAIL
    SEQUENCE_FASTA += QUERY_TAIL

    def _get_df(db_conn, accessions, format):
        acclist = ["'" + item + "'" for item in accessions]
        accstring = ",".join(acclist)
        if format == "DataFrame":
            query = SEQUENCE_FULL.format(accstring)
        elif format == "DataFrameNoSeq":
            query = SEQUENCE_INFO.format(accstring)
        else:  # format=dict
            query = SEQUENCE_FASTA.format(accstring)

        return pd.read_sql_query(query, index_col="ACC", con=db_conn)

    if len(accessions) > MAXLIST:  # work in batches
        df = pd.DataFrame()
        for accbatch in batch(accessions, MAXLIST):
            df = df.append(_get_df(db_conn, accbatch, format))
    else:
        df = _get_df(db_conn, accessions, format)

    if format == "Dict":
        df["SEQUENCE"] = df["SEQUENCE"].astype(str)  # convert from oracle LOB to str
        return df.to_dict()["SEQUENCE"]
    else:
        if format == "DataFrame":
            df["SEQUENCE"] = df["SEQUENCE"].astype(
                str
            )  # convert from oracle LOB to str
            df.columns = ["entry_type", "seqlen", "sequence", "md5sum"]
        elif format == "DataFrameNoSeq":
            df.columns = ["entry_type", "seqlen", "md5sum"]
        else:
            eprint("UNKNOWN FORMAT SPECIFIED")
        return df


def get_sequences_uniparc(
    accessions, orthoid="", format="Dict", config={}, db_conn=None
):
    """
    input: a list of accessions and optionally an orthoid (for error msgs only) and a format
    returns: depending on the format parameter, the sequence information in that format
    formats: DataFrame, DataFrameNoSeq (same as DataFrame but without the sequence, only sequence info), Dict
    """
    MAXLIST = 900  # to avoid ORA-01795: maximum number of expressions in a list is 1000

    SEQUENCE_FULL = """
        select
        ac as acc,
        (case when dbid in (2,24) then 'sp' when dbid=3 then 'tr' else 'ref' end) as entry_type,
        len as seqlen,
        nvl(seq_long,seq_short) as sequence,
        md5 as md5sum
        """

    SEQUENCE_INFO = """
        select
        ac as acc,
        (case when dbid in (2,24) then 'sp' when dbid=3 then 'tr' else 'ref' end) as entry_type,
        len as seqlen,
        md5 as md5sum
        """

    SEQUENCE_FASTA = """
        select ac as acc,
        '>'||(case when dbid in (2,24) then 'sp' when dbid=3 then 'tr' else 'ref' end)||'|'||ac||'|'||x.upi||'\n'||nvl(seq_long,seq_short)||'\n'  as sequence
        """

    QUERY_TAIL = """
        from xref x inner join protein p on p.upi=x.upi
        WHERE deleted='N' --only active sequences, ignoring obsolete
        and dbid in (2, 3, 24, 41) --sp (2:SWISSPROT) trembl (3:TREMBL) sp isoforms (24:SWISSPROT_VARSPLIC) refseq (41)
        and ac in ({})
        """

    SEQUENCE_FULL += QUERY_TAIL
    SEQUENCE_INFO += QUERY_TAIL
    SEQUENCE_FASTA += QUERY_TAIL

    def _get_df(db_conn, accessions, format):
        acclist = ["'" + item + "'" for item in accessions]
        accstring = ",".join(acclist)
        if format == "DataFrame":
            query = SEQUENCE_FULL.format(accstring)
        elif format == "DataFrameNoSeq":
            query = SEQUENCE_INFO.format(accstring)
        else:  # format=dict
            query = SEQUENCE_FASTA.format(accstring)

        return pd.read_sql_query(query, index_col="ACC", con=db_conn)

    if len(accessions) > MAXLIST:  # work in batches
        df = pd.DataFrame()
        for accbatch in batch(accessions, MAXLIST):
            df = df.append(_get_df(db_conn, accbatch, format))
    else:
        df = _get_df(db_conn, accessions, format)

    if format == "Dict":
        df["SEQUENCE"] = df["SEQUENCE"].astype(str)  # convert from oracle LOB to str
        return df.to_dict()["SEQUENCE"]
    else:
        if format == "DataFrame":
            df["SEQUENCE"] = df["SEQUENCE"].astype(
                str
            )  # convert from oracle LOB to str
            df.columns = ["entry_type", "seqlen", "sequence", "md5sum"]
        elif format == "DataFrameNoSeq":
            df.columns = ["entry_type", "seqlen", "md5sum"]
        else:
            eprint("UNKNOWN FORMAT SPECIFIED")
        return df


def get_sequences_web(accessions, orthoid="", format="Dict", config={}):
    """
    input: a list of accessions and optionally an orthoid (for error msgs only) and a format
    returns: depending on the format parameter, the sequence information in that format
    formats: DataFrame, DataFrameNoSeq (same as DataFrame but without the sequence, only sequence info), Dict
    """
    REFASTAHEADER_SIMPLE = re.compile(
        r"^>(\S*) .*?\n"
    )  # strip header portion after first blank
    REFASTAHEADER_UPTYPE = re.compile(r"^>(sp|tr)|.*\n")  # get sp or tr

    method = "api"  # ebi protein api

    if method == "api":
        headers = {"user-agent": "ortho/1.1.0", "Accept": "text/x-fasta"}
    else:  # rest
        headers = {"user-agent": "ortho/1.1.0", "Accept": "application/json"}
    sequences = []
    seqlens = []
    md5sums = []
    geneids = []
    accs_found = []
    entry_types = []
    for acc in accessions:
        seq = None
        entry_type = ""
        seqraw = ""
        geneid = 0  # default if none found
        read_from_cache = False
        if config.get("cache_sequences_flag", False):
            fastaseq_fname = os.path.join(config["fasta_dir"], acc + ".fasta")
            if os.path.isfile(fastaseq_fname):
                # eprint("Reading {} from sequence cache".format(acc)) #debug
                with open(fastaseq_fname) as fastaseq_fd:
                    seq = fastaseq_fd.read()
                # sequences[acc] = seq #OLD
                read_from_cache = True
        if seq is None or seq == "":
            # eprint("Fetching {} from api".format(acc)) #debug
            # sequence retrieval from api
            if acc.find("_") != -1:  # ncbi refseq, get from uniparc via query
                # website (slower)
                # url = "https://www.uniprot.org/uniparc/?query=database%3ARefSeq+{}&sort=score&format=fasta".format(acc)
                # api (faster)
                if method == "api":
                    url = "https://www.ebi.ac.uk/proteins/api/uniparc?dbtype=RefSeq&dbid={}&rfActive=true&size=1".format(
                        acc
                    )
                else:  # rest
                    url = "https://rest.uniprot.org/beta/uniparc/search?query={}&fields=sequence".format(
                        acc
                    )
                    # TODO
                entry_type = "ref"

            else:  # uniprot
                if method == "api":
                    url = "https://www.ebi.ac.uk/proteins/api/proteins/{}".format(acc)
                else:
                    url = "https://rest.uniprot.org/beta/uniprotkb/{}?fields=sequence,reviewed,xref_geneid".format(
                        acc
                    )

            r = get(url, headers=headers)
            # if not r.ok:
            #    r.raise_for_status()
            if r.status_code in (200, 206):
                if method == "api":
                    seq = REFASTAHEADER_SIMPLE.sub(
                        r">\1\n", r.text
                    )  # keep only identifier
                else:
                    if entry_type != "ref":
                        entry_type = r.json()["entryType"]
                        if entry_type == "UniProtKB reviewed (Swiss-Prot)":
                            entry_type = "sp"
                        else:
                            entry_type = "tr"
                    seqjson = r.json()["sequence"]
                    seqraw = seqjson["value"]
                    xrefs = r.json()["uniProtKBCrossReferences"]
                    if type(xrefs) == list and len(xrefs) > 0:
                        geneid = xrefs[0]["id"]
                        # eprint('found geneid: {}'.format(geneid)) #debug
                    seq = ">{}|{}\n{}".format(entry_type, acc, seqraw)
                if config.get("cache_sequences_flag", False):
                    with open(fastaseq_fname, "w") as fastaseq_fd:
                        fastaseq_fd.write(seq)
            else:
                eprint(
                    "{} ERROR: problems retrieving acc {}, http status code is '{}' ({})".format(
                        orthoid, acc, r.status_code, r.reason
                    )
                )
                continue

        if seq is not None:
            # annotation/labelling of sequence headers
            if method == "api" or read_from_cache:  # api fasta
                seqraw = (
                    REFASTAHEADER_STRIP.sub("", seq).replace("\n", "").upper()
                )  # only the aa residues
                seqlens.append(len(seqraw))
                md5sums.append(md5(seqraw.encode()).hexdigest().upper())
                if entry_type == "ref" or (acc.find("_") != -1):
                    entry_type = "ref"
                    seq = ">ref|" + seq[1:]
                else:
                    entry_type = REFASTAHEADER_UPTYPE.sub(r"\1", seq)
            else:  # rest json
                seqlens.append(seqjson["length"])
                md5sums.append(seqjson["md5"])
            accs_found.append(acc)
            sequences.append(seq)
            geneids.append(geneid)
            entry_types.append(entry_type)
    if format == "DataFrame":
        df = pd.DataFrame(
            {
                "acc": accs_found,
                "entry_type": entry_types,
                "seqlen": seqlens,
                "sequence": sequences,
                "md5sum": md5sums,
                "geneid": geneids,
            },
            columns=["acc", "entry_type", "seqlen", "sequence", "md5sum", "geneid"],
        )
        df.set_index("acc", inplace=True)
        return df
    elif format == "DataFrameNoSeq":
        df = pd.DataFrame(
            {
                "acc": accs_found,
                "entry_type": entry_types,
                "seqlen": seqlens,
                "md5sum": md5sums,
                "geneid": geneids,
            },
            columns=["acc", "entry_type", "seqlen", "md5sum", "geneid"],
        )
        df.set_index("acc", inplace=True)
        return df
    else:  # dict
        return dict(zip(accs_found, sequences))


def get_orthologs_df_from_pantherid(orthoid, ortho_df, remove_outliers=False):
    if remove_outliers:
        return ortho_df[(ortho_df["pantherid"] == orthoid) & ~(ortho_df["outlier"])]
    else:
        return ortho_df[ortho_df["pantherid"] == orthoid]
