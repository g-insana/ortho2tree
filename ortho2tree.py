#!/usr/bin/env python
# coding: utf-8

# # Ortho2tree pipeline
# ### by Giuseppe Insana and William R Pearson
# ### https://github.com/g-insana/ortho2tree

# In[ ]:


import argparse
import os
import sys
import time  # elapsed times
import yaml  # for configuration files
from multiprocessing import Pool, cpu_count, current_process, set_start_method
from platform import platform

from ortho2tree.o2t_buildtree import build_tree_for_orthogroup
from ortho2tree.o2t_df import create_ortho_df, ortho_df_stats, read_ortho_df
from ortho2tree.o2t_gc_integration import dump_prev_changes
from ortho2tree.o2t_load import get_gc_df, get_panther_df
from ortho2tree.o2t_output import (
    clean_up_tempfiles,
    combine_and_print_output,
    output_headers,
)
from ortho2tree.o2t_scan_ndata import scan_ndata_file
from ortho2tree.o2t_utils import (
    check_all_files_exist,
    elapsed_time,
    eprint,
    get_orthologs_df_from_pantherid,
    get_sequences_uniparc,
    get_sequences_swpread,
    get_sequences_web,
)

sys.stderr.write("*** Ortho2tree pipeline ***\n")

# for multithreading
if sys.version_info >= (3, 8, 10) and platform().find("macOS") != -1:
    sys.stderr.write("INFO: Using fork MP context macos and py >=3.8.10\n")
    try:
        # NOTE: because spawn fails due to freeze_support
        set_start_method("fork")
    except RuntimeError:
        sys.stderr.write("NOTICE: context already set\n")

# optional, for progressbar
tqdm_available = True
try:
    from tqdm.auto import tqdm
except ImportError:
    tqdm_available = False

# for oracle database connection
try:
    import cx_Oracle
except ImportError:
    sys.stderr.write("NOTICE: Oracle module not available\n")

# optional, for connection engine
sqlalchemy_available = True
try:
    import sqlalchemy
except ImportError:
    sqlalchemy_available = False
    sys.stderr.write("NOTICE: sqlalchemy module not available\n")


# ## Part 1) Configuration and setup

# In[ ]:


# Note: any or all of the following parameters can be overridden with the yaml configuration file config = { # # # Performance and UI options # # #
config = {
    "threads": 7,  # how many parallel threads?
    "progressbar": True,  # show tqdm progressbar (if tqdm module is present)
    "print_stats": True,  # print statistics about the df
    # # # DB, API and Files options # # #
    "gc_from_sql": True,  # use oracle db sql query for gc groups? (if false, get data from tsv file)
    "seq_from_sql": True,  # use oracle db sql queries to retrieve sequences? (if false, get data via protein API)
    "use_uniparc_for_seq_retrieval": True,  # otherwise use swpread
    "cache_sequences_flag": True,  # if getting sequence data via protein API, should we cache the files?
    "cache_alignments_flag": True,  # should we store and cache alignments?
    "cache_trees_flag": True,  # should we store and cache trees?
    "create_pdf_files_flag": True,  # should we create pdf files for suggestions? implies create lab_lt and faX files
    "create_pdf_files4confirmed_flag": True,  # should we create pdf files also for confirmed canonical solutions?
    "panther_data_dir": "PANTHER_Sequence_Classification_files/",  # where to store local copy of panther data?
    "fasta_dir": "fasta",  # directory name where sequences cache can be stored
    "geneid2refseq_mapfile": "refprots_geneid2refseq.csv.gz",  # geneid 2 refseq mapping file
    "missing_geneid_mapfile": "up_refseq_geneid_map95_90.tab",  # missing geneid 2 up acc mapping file
    # Note for following paths: dataset_name/ is the same name as the cfg file. e.g. 5taxa/ for 5taxa.cfg
    "n_data_dir": "n_data",  # directory name where output n_data will be written to (under dataset/)
    "lab_data_dir": "lab_data",  # directory name where output .lab_lt files will be written to (under dataset/)
    "aln_data_dir": "aln_data",  # directory name where output alignments will optionally be written to (under dataset/)
    "tree_data_dir": "tree_data",  # directory name where trees will optionally be written to (under dataset/)
    "pdf_data_dir": "pdf_data",  # directory name where pdf files for suggestions will optionally be written to (under dataset/)
    "semaphores_dir": "processed",  # directory name where to create semaphore files (under dataset/)
    # # # DATA version options # # #
    "panther_version": "PTHR17.0",  # Panther version
    "up_release": "2024_01",  # UniProt Release version
    # # # Parameters for the ANALYSIS # # #
    "min_taxa_threshold": 3,  # minimum number of different organisms for building an alignment
    "taxa_threshold": 3,  # number of different organisms that should be in the low-cost clade for acceptance (if taxa in tree more than this, otherwise use min_taxa_threshold)
    "tree_max_cost": 0.05,  # exclude tree_to_ndata solutions with costs higher than this (default=0.05)
    "tree_drop_fact": 1.5,  # improvement required to drop a taxon
    "superfamily_level": False,  # set to true to work at superfamily level
    "add_refseq": False,  # whether we want to add refseq sequences to orthogroups
    # # parameters used when scanning n_data2 files to select suggestions # #
    "min_delta_threshold": 0.005,  # minimum cost difference for individual taxa
    "min_delta_sp_threshold": 0.02,  # minimum cost difference for individual sp taxa
    "suggestion_score_difference": -0.001,  # cost difference
    "suggestion_taxa_threshold": 3,  # minimum number of taxa
    "suggestion_min_canon": 1,  # minimum number of canonicals
    "suggestion_max_clade_cost": 0.02,  # max clade cost
    "suggestion_only_zero_cost": False,  # only consider suggestions with zero cost
    "suggestion_ranking_weights": {
        "n_sp": 16.0,
        "n_tax": 4.0,
        "n_canon": 0.0,
        "wn_canon": 12.0,
        "scaled_prop_f": 4.0,
        "scaled_p_cost": 8.0,
        "p_len_diff": 1.0,
    },
    "suggestion_taxon_weights": {"HUMAN": 2.0, "MOUSE": 2.0, "RAT": 1.0},
    "suggestion_taxon_weight_default": 0.5,
    "skip_reprocessing_orthogroups": False,  # set this to skip (re)processing of orthogroups, useful if only want to re-score existing n_data files
    # # # OUTLIERS options # # #
    "remove_outliers": True,  # should we remove sequences that have been flagged as outliers?
    "max_seqlen": 10000,  # maximum sequence length, sequences beyond this length (and their isoforms) will not be considered
    # std/mean based outlier detection:
    # outliers are identified as those under /threshold_low/ or beyond /threshold_hi/ standard deviations from mean of sequence lengths
    "detect_outliers_with_mean_std": False,
    "outliers_detection_threshold_std_lo": 2,
    "outliers_detection_threshold_std_hi": 2,
    # median based outlier detection:
    # outliers are identified as those with seqlen under /threshold_median_lo/ or beyond /threshold_median_hi/ times the median
    "detect_outliers_with_median": False,
    "outliers_detection_threshold_median_lo": 0.75,
    "outliers_detection_threshold_median_hi": 2,
    # quartile based outlier detection:
    # outliers are identified as those with seqlen under /threshold_quart_lo/ times Q1 or beyond /threshold_quart_hi/ times the Q3 value
    "detect_outliers_with_quart": False,
    "outliers_detection_threshold_quart_lo": 0.5,
    "outliers_detection_threshold_quart_hi": 2,
    # canonical seqlen outlier detection:
    # outliers are identified as those with seqlen under /threshold_can_lo/ times the can_min_len (min seqlen of canonicals in group) or beyond /threshold_can_hi/ times the can_max_len (max seqlen of canonicals in group)
    "outliers_detection_threshold_can_lo": 0.5,
    "outliers_detection_threshold_can_hi": 2,
    # median + canonical seqlen outlier detection:
    "detect_outliers_with_median_and_can_lengths": False,
    # quartile + canonical seqlen outlier detection:
    "detect_outliers_with_quart_and_can_lengths": True,
    # # # DATAFRAME reading and caching # # #
    "dump_orthogroup_data": True,  # do we want to dump the data into output files marked with a certain /outstamp/?
    "outstamp": "240102",  # this will be appended to output file names
    "use_cached_orthogroup_data": False,  # do we want to read cached dataframe from previously dumped files, marked with a certain /instamp/?
    "instamp": "240101",  # this will be appended to input file names
    # # # GCINTEGRATION # # #
    "sugg_file": False,  # if a filename is specified, it will be read to simulate genecentric application of ortho2tree using previously created suggestions
    "prevgc_file": False,  # if a filename is specified, it will be read and used for cumulative integrated changes for genecentric pipeline
    # # # ORGANISMS definition # # #
    # use panther organism names as values in the following dict, with tax_id as keys
    "tax2org": {9606: "human", 10090: "mouse", 10116: "rat", 9913: "cow", 9615: "dog"},
}


# ### Configuration override from a yaml configuration file, if present and parsing of arguments

# In[ ]:


config["groups2run"] = []
config["debug_mode"] = False
config["ignore_cache"] = False

# if we are running in the notebook (ortho2tree.ipynb)
if sys.argv[0].find("pykernel_launcher") != -1:
    # default set for notebook; override to test other sets
    config["dataset_name"] = "pln2024_01"
else:  # we are running from the shell as ortho2tree.py
    USAGE_EXAMPLE = """
    Examples:
       -set=5taxa #will do the analysis on the whole set
       -set=5taxa -id=PTHR19918:SF1               #only for one orthogroup
       -set=5taxa -id=PTHR19918:SF1 PTHR40139:SF1 #only for two orthogroups
       -set=mamRS -file=list_of_ids.txt           #for a series of groups listed in a file
    """
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=USAGE_EXAMPLE
    )
    parser.add_argument(
        "-set",
        dest="dataset_name",
        required=True,
        type=str,
        help="set for the analysis. a file SET.cfg should be present",
    )
    parser.add_argument(
        "-d",
        dest="debug_mode",
        required=False,
        action="store_true",
        help="print verbose/debug messages",
    )
    parser.add_argument(
        "-nocache",
        dest="ignore_cache",
        required=False,
        action="store_true",
        help="do not use cache, re-create aln/trees and do not save them",
    )
    parser.add_argument(
        "-no_stats",
        dest="no_stats",
        required=False,
        action="store_true",
        help="do not print any stats on the dataframe",
    )
    parser.add_argument(
        "-id",
        dest="single_group",
        required=False,
        type=str,
        nargs="+",
        help="to only work on one or few group(s)",
    )
    parser.add_argument(
        "-file",
        dest="list_filename",
        required=False,
        type=argparse.FileType("r"),
        help="to work on a series of groups, from a file",
    )
    parser.add_argument(
        "-sugg",
        dest="sugg_file",
        required=False,
        type=str,
        help="to simulate integration of canonical suggestions reading a previosly generated changes file; note file should be placed in the set main dir",
    )
    parser.add_argument(
        "-prevgc",
        dest="prevgc_file",
        required=False,
        type=str,
        help="to integrate previosly generated changes file; note file should be placed in the set main dir",
    )
    parser.add_argument(
        "-outstamp",
        dest="outstamp",
        required=False,
        type=str,
        help="to name and timestamp the output files and the dumps; this overrides the outstamp parameter from the config",
    )
    args = parser.parse_args()
    if args.single_group is not None and args.list_filename is not None:
        eprint("    => ERROR: either pass -id or -file, not both")
        sys.exit(22)
    if args.ignore_cache:
        config["ignore_cache"] = True
    if args.debug_mode:
        config["debug_mode"] = True
    if args.single_group is not None:
        config["groups2run"] = list(args.single_group)
    if args.list_filename is not None:
        for line in args.list_filename.readlines():
            config["groups2run"].append(line.rstrip())
    if args.sugg_file is not None:
        if os.path.isfile(os.path.join(args.dataset_name, args.sugg_file)):
            eprint(
                "NOTICE: Suggestion file override from command line: {}".format(
                    args.sugg_file
                )
            )
            # override default from config
            config["sugg_file"] = args.sugg_file
        else:
            eprint("    => ERROR: file specified as -sugg does not exist!")
            sys.exit(2)
    if args.prevgc_file is not None:
        if os.path.isfile(os.path.join(args.dataset_name, args.prevgc_file)):
            eprint(
                "NOTICE: prevgc_file override from command line: {}".format(
                    args.prevgc_file
                )
            )
            # override default from config
            config["prevgc_file"] = args.prevgc_file
        else:
            eprint("    => ERROR: file specified as -prevgc does not exist!")
            sys.exit(2)

    config["dataset_name"] = args.dataset_name

config["dataset_maindir"] = config["dataset_name"]
dataset_configfile = config["dataset_name"] + ".cfg"
if os.path.isfile(dataset_configfile):
    eprint(
        "NOTICE: Default configuration overridden by YAML file {}".format(
            dataset_configfile
        )
    )
    with open(dataset_configfile, encoding="utf-8") as fh:
        config.update(yaml.safe_load(fh))

# if we are running from script with argparse
if sys.argv[0].find("pykernel_launcher") == -1:
    if args.outstamp is not None:
        config["outstamp"] = args.outstamp
    if args.no_stats:
        config["print_stats"] = False
else:  # we are in jupyter
    if sys.version_info >= (3, 8, 10) and platform().find("macOS") != -1:
        eprint(
            "NOTICE: Forcing single thread mode due to multiprocessing issue with jupyter macOS py>=3.8.10"
        )
        config["threads"] = 1

if not tqdm_available:
    config["progressbar"] = False

eprint("\nINFO: Configured parameters: {}".format(config))


# ### File definitions, creation of paths

# In[ ]:


# timestamped files:
file_keys = [
    "orthogroup_df_cachefile",
    "orthogroup_df_dumpfile",  # to read/write the mapped and processed dataframe
    "prevgc_notfound_file",  # entries no more present in gc but that were previous gc suggestions
    "prevgc_conflict_file",  # conflicting suggestions being removed
    "groups_by_taxa_count_file",  # list of groups and their sizes in number of taxa
    "groups_by_entry_count_file",  # list of groups and their sizes in number of entries, excluding lowtaxa and outliers
    "low_taxa_groups_file",  # list of groups with taxa size lower than threshold
    "tr_fragments_df_cachefile",
    "tr_fragments_df_dumpfile",  # to read or write trembl fragments
    "unmapped_canonicals_file",  # list of canonicals not mapped
]
for key in file_keys:
    if key.endswith("_cachefile"):
        config[key] = f"{key[:-10]}_dump{config['instamp']}.gz"
    elif key.endswith("_dumpfile"):
        config[key] = f"{key[:-9]}_dump{config['outstamp']}.gz"
    elif key in ["unmapped_canonicals_file", "low_taxa_groups_file"]:
        config[key] = f"{key[:-5]}{config['outstamp']}"
    else:  # compressed
        config[key] = f"{key[:-5]}{config['outstamp']}.gz"


# output intermediate files:
config["output_keys"] = ["conflict", "changes", "confirm", "skipped", "gc"]
for key in config["output_keys"]:
    config[key + "_outfile"] = "output_" + key + config["outstamp"]

# headers for output files
headers = output_headers(config=config)

# taxa and tax_ids
if config["taxa_threshold"] > len(config["tax2org"]):
    eprint(
        "WARNING: the specified level taxa_threshold ({}) was higher than the number of species defined for the analysis ({}). Using the latter for threshold".format(
            config["taxa_threshold"], len(config["tax2org"])
        )
    )
    config["taxa_threshold"] = len(config["tax2org"])
eprint(
    "\nINFO: Number of species defined for the analysis: {}, threshold for n_data printout: {}".format(
        len(config["tax2org"]), config["taxa_threshold"]
    )
)
config["tax_ids"] = set(config["tax2org"].keys())

# file paths creation
if not os.path.exists(config["dataset_maindir"]):
    os.mkdir(config["dataset_maindir"])

# path join
paths_to_join = [
    "n_data_dir",
    "lab_data_dir",
    "aln_data_dir",
    "tree_data_dir",
    "pdf_data_dir",
    "semaphores_dir",
    "unmapped_canonicals_file",
    "orthogroup_df_dumpfile",
    "tr_fragments_df_dumpfile",
    "prevgc_notfound_file",
    "prevgc_conflict_file",
    "groups_by_taxa_count_file",
    "low_taxa_groups_file",
    "orthogroup_df_cachefile",
    "tr_fragments_df_cachefile",
    "gc_outfile",
    "changes_outfile",
    "confirm_outfile",
    "skipped_outfile",
    "conflict_outfile",
]
for key in paths_to_join:
    config[key] = os.path.join(config["dataset_maindir"], config[key])

# cached dataframe
if config["use_cached_orthogroup_data"] and not os.path.isfile(
    config["orthogroup_df_cachefile"]
):
    eprint("    => ERROR: no such file {}".format(config["orthogroup_df_cachefile"]))
    sys.exit(2)


# gc integration
files_to_check = ["sugg_file", "prevgc_file"]
for key in files_to_check:
    file_path = config[key]
    if file_path:
        if file_path == "False":
            config[key] = False
        else:
            config[key] = os.path.join(config["dataset_maindir"], file_path)

if config["ignore_cache"]:
    config["cache_alignments_flag"] = False
    config["cache_trees_flag"] = False
    config["cache_sequences_flag"] = False

directories_to_create = [
    "panther_data_dir",
    "n_data_dir",
    "lab_data_dir",
    "aln_data_dir",
    "tree_data_dir",
    "semaphores_dir",
]

config["create_lablt_files_flag"] = False
config["create_faX_files_flag"] = False
if config["create_pdf_files_flag"]:
    # the following files are needed to create the pdf files for the suggestions
    config["create_lablt_files_flag"] = True
    config["create_faX_files_flag"] = True
    directories_to_create.append("pdf_data_dir")

if config["cache_sequences_flag"]:
    directories_to_create.append("fasta_dir")

for key in directories_to_create:
    directory_path = config[key]
    if not os.path.exists(directory_path):
        os.mkdir(directory_path)


# ### Checks for db access (UniProt only)

# In[ ]:


if "cx_Oracle" not in sys.modules:
    config["gc_from_sql"] = False
    config["seq_from_sql"] = False

config["db_connection"] = ""
if config["gc_from_sql"]:
    if not os.path.isfile("swpread_connection.pass"):
        config["gc_from_sql"] = False
    else:  # read connection details
        with open(
            "swpread_connection.pass",
            "r",
            encoding="utf-8",  # for uniprotkb (format: username/password@dbname)
        ) as f_in:
            config["db_connection"] = f_in.read().rstrip()
    try:
        cx_Oracle.connect(config["db_connection"])
    except cx_Oracle.DatabaseError as e:
        print("WARNING: No db connection to uniprot for Genecentric data: {}".format(e))
        config["gc_from_sql"] = False

config["db_connection_uniparc"] = ""
if config["seq_from_sql"]:
    if not os.path.isfile("uatst_connection.pass"):
        config["seq_from_sql"] = False
    else:
        with open(
            "uatst_connection.pass", "r", encoding="utf-8"  # for uniparc
        ) as f_in:
            config["db_connection_uniparc"] = f_in.read().rstrip()
        try:
            cx_Oracle.connect(config["db_connection_uniparc"])
        except cx_Oracle.DatabaseError as e:
            print(
                "WARNING: no db connection to uniparc: '{}'; reverting to API for sequence retrieval".format(
                    e
                )
            )
            config["seq_from_sql"] = False


def dbconnect(dbconnection):
    """
    Connect to db and returns cx_Oracle instance.
    """
    try:
        return cx_Oracle.connect(dbconnection)
    except PermissionError as excep:
        eprint(
            "    => ERROR: You have not specified the correct username/password@dbname",
            file=sys.stderr,
        )
        eprint(str(excep), file=sys.stderr)
        sys.exit(11)


# ### Wrapper functions

# In[ ]:


def get_sequences_wrapper(config={}):
    """
    wrapper to have get_sequences appropriate to current setup (via db or via web api)
    """

    if config["seq_from_sql"]:  # get sequences via database access
        if config["use_uniparc_for_seq_retrieval"]:
            db_connection_details = config["db_connection_uniparc"]
            get_sequences_db_function = get_sequences_uniparc
            db_database_name = "UNIPARC"
        else:
            db_connection_details = config["db_connection"]
            get_sequences_db_function = get_sequences_swpread
            db_database_name = "SWPREAD"

        if sqlalchemy_available:
            eprint(
                "We will retrieve sequences using {} database via sqlalchemy".format(
                    db_database_name
                )
            )
            engine = sqlalchemy.create_engine(
                "oracle+cx_oracle://{}".format(db_connection_details.replace("/", ":")),
                pool_size=config["threads"],
                max_overflow=2,
            )

            def _wrapper(accessions, orthoid="", format="Dict", config={}):
                try:
                    db_conn = engine.pool.connect()
                    return get_sequences_db_function(
                        accessions,
                        orthoid=orthoid,
                        format=format,
                        config=config,
                        db_conn=db_conn.driver_connection,
                    )
                except Exception as e:
                    eprint("ERROR fetching sequences: {}".format(e))
                finally:
                    db_conn.close()

        else:
            eprint(
                "We will retrieve sequences using {} database via cx_Oracle".format(
                    db_database_name
                )
            )

            db_user, db_tmp = db_connection_details.split("/", 1)
            db_pwd, db_name = db_tmp.split("@", 1)
            oracle_pool = cx_Oracle.SessionPool(
                user=db_user,
                password=db_pwd,
                dsn=db_name,
                min=1,
                max=config["threads"] + 2,
            )
            db_conn = oracle_pool.acquire()
            oracle_pool.release(db_conn)

            def _wrapper(accessions, orthoid="", format="Dict", config={}):
                try:
                    db_conn = oracle_pool.acquire()
                    return get_sequences_uniparc(
                        accessions,
                        orthoid=orthoid,
                        format=format,
                        config=config,
                        db_conn=db_conn,
                    )
                except cx_Oracle.DatabaseError as exc:
                    err = exc.args
                    eprint("Oracle-Error-Code:", err.code)
                    eprint("Oracle-Error-Message:", err.message)
                except Exception as e:
                    eprint("ERROR fetching sequences: {}".format(e))
                finally:
                    oracle_pool.release(db_conn)

    else:  # get sequences via online API retrieval
        eprint("We will retrieve sequences using WEB API")

        def _wrapper(accessions, orthoid="", format="Dict", config={}):
            return get_sequences_web(
                accessions, orthoid=orthoid, format=format, config=config
            )

    return _wrapper


get_sequences = get_sequences_wrapper(config=config)


def build_tree_for_orthogroup_wrapper(orthoid):
    """
    wrapper to avoid reprocessing or to force reprocessing
    """
    semaphore_file = os.path.join(config["semaphores_dir"], orthoid + ".done")
    if FORCE_REPROCESS:
        build_tree_for_orthogroup(
            orthoid,
            ortho_df=ortho_df,
            verbose=VERBOSE_PROCESS,
            config=config,
            get_sequences_fn=get_sequences,
        )
        if not os.path.isfile(semaphore_file):
            open(semaphore_file, "a").close()  # touch semaphore file
    else:
        if not os.path.isfile(semaphore_file):  # to avoid recreating
            build_tree_for_orthogroup(
                orthoid,
                ortho_df=ortho_df,
                verbose=VERBOSE_PROCESS,
                config=config,
                get_sequences_fn=get_sequences,
            )
            open(semaphore_file, "a").close()  # touch semaphore file


def scan_ndata_file_wrapper(orthoid):
    """
    multithread wrapper in order to write into temporary files, one for each thread
    """
    output_data = {}
    (
        output_data["gc"],
        output_data["changes"],
        output_data["confirm"],
        output_data["skipped"],
        output_data["conflict"],
    ) = scan_ndata_file(
        orthoid,
        ortho_df=ortho_df,
        prevgc_df=prevgc_df,
        config=config,
        get_sequences_fn=get_sequences,
    )

    if current_process().name == "MainProcess":
        workerid = "XX"
    else:
        workerid = "{}".format(current_process().name.split("-")[1])

    for key in config["output_keys"]:
        if output_data[key] != "":
            tempfile = "{}..{}".format(config[key + "_outfile"], workerid)
            with open(tempfile, "a", encoding="utf-8") as outputfh:
                outputfh.write(output_data[key])


# # Part 2.1: Load of Panther Data

# In[ ]:


if not config["use_cached_orthogroup_data"]:
    start_secs = time.time()
    panther_df, config["tax2oscode"] = get_panther_df(config=config)
    eprint(
        "\n*** Part 2.1 workflow completed {} -- Elapsed: {} --\n".format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            elapsed_time(start_secs),
        )
    )


# # Part 2.2: Load of Genecentric Data

# In[ ]:


if not config["use_cached_orthogroup_data"]:
    start_secs = time.time()
    db_conn = None
    if config["gc_from_sql"]:
        db_conn = dbconnect(config["db_connection"])
    gc_df = get_gc_df(config=config, db_conn=db_conn)

    eprint(
        "Loaded {} gc accessions, of which {} canonicals".format(
            len(gc_df), len(gc_df[gc_df["is_canonical"]])
        )
    )
    eprint(
        "\n*** Part 2.2 workflow completed {} -- Elapsed: {} --\n".format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            elapsed_time(start_secs),
        )
    )


# # Part 3: Data processing

# In[ ]:


start_secs = time.time()
if config["use_cached_orthogroup_data"]:
    ortho_df, prevgc_df = read_ortho_df(config=config)
else:
    ortho_df, prevgc_df = create_ortho_df(
        panther_df, gc_df, config=config, get_sequences_fn=get_sequences
    )
    # optionally, free memory by deleting df no more needed
    del panther_df
    del gc_df

if len(config["groups2run"]) == 0:  # unless we are doing single groups
    # print stats and filter out groups lower than min_taxa
    ortho_df = ortho_df_stats(ortho_df, config=config)

orthogroups = set(ortho_df["pantherid"].drop_duplicates().values)
all_canonicals = set(ortho_df[ortho_df["is_canonical"]].index.to_list())

eprint(
    "We'll work on {} canonicals in {} orthogroups with total {} accessions".format(
        len(all_canonicals), len(orthogroups), len(ortho_df)
    )
)
eprint(
    "\n*** Part 3 workflow completed {} -- Elapsed: {} --\n".format(
        time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()), elapsed_time(start_secs)
    )
)


# # Part 4: Build alignments and trees, identify clades

# In[ ]:


# set to True to reprocess those already done
FORCE_REPROCESS = False
# set to True when debugging for high verbosity of each group
VERBOSE_PROCESS = False
# set to False to skip checking for all groups to be done
CHECK_ALL_PROCESSED = True


# In[ ]:


# if we want to skip or remove any group (e.g. with issues or bad data)
problematic_groups = []
for problematic_group in problematic_groups:
    if problematic_group in orthogroups:
        eprint("removed problematic group {}".format(problematic_group))
        orthogroups.remove(problematic_group)


# In[ ]:


if config["skip_reprocessing_orthogroups"] or len(config["groups2run"]):
    # we won't (re)process orthogroups and instead simply re-score, parsing existing n_data files
    eprint(
        "\n*** Part 4 workflow skipped due to skip_reprocessing_orthogroups config option or specified list of groups\n"
    )
else:
    start_secs = time.time()
    if config["threads"] > 1:
        eprint(
            "Working in {} parallel threads; your OS reports {} cpus.".format(
                config["threads"], cpu_count()
            )
        )
        pool = Pool(config["threads"])
        iterator = pool.imap(build_tree_for_orthogroup_wrapper, orthogroups)
        if config["progressbar"]:
            _ = list(tqdm(iterator, total=len(orthogroups)))
        else:
            iterator
        pool.close()  # no more work to submit
        pool.join()  # wait workers to finish
    else:
        iterator = tqdm(orthogroups) if config["progressbar"] else orthogroups
        for orthoid in iterator:
            build_tree_for_orthogroup_wrapper(orthoid)
    eprint(
        "\n*** Part 4 workflow completed {} -- Elapsed: {}, {} g/s --\n".format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            *elapsed_time(start_secs, len(orthogroups)),
        )
    )
    if CHECK_ALL_PROCESSED:
        check_all_files_exist(config["semaphores_dir"], orthogroups, "done")


if len(config["groups2run"]):
    groups_not_found = []
    for single_group in config["groups2run"]:
        if single_group not in orthogroups:
            eprint("  ERROR: '{}' not in the df".format(single_group))
            groups_not_found.append(single_group)
            continue
        if config["debug_mode"]:
            eprint("  Running group '{}'".format(single_group))
            eprint(
                get_orthologs_df_from_pantherid(single_group, ortho_df)
                .reset_index()
                .sort_values(
                    by=["org", "is_canonical", "outlier", "entry_type", "acc"],
                    ascending=[True, False, False, True, True],
                )[
                    [
                        "entry_type",
                        "acc",
                        "org",
                        "seqlen",
                        "groupid",
                        "is_canonical",
                        "outlier",
                    ]
                ]
                .set_index("acc")
            )
            build_tree_for_orthogroup(
                single_group,
                ortho_df=ortho_df,
                verbose=config["debug_mode"],
                debuginfo=config["debug_mode"],
                config=config,
                get_sequences_fn=get_sequences,
            )
        else:
            build_tree_for_orthogroup(
                single_group,
                ortho_df=ortho_df,
                verbose=config["debug_mode"],
                config=config,
                get_sequences_fn=get_sequences,
            )
    config["groups2run"] = [
        x for x in config["groups2run"] if x not in groups_not_found
    ]


# # Part 5: parse clades and create output files

# In[ ]:


if len(config["groups2run"]):
    output_all = {key: "" for key in config["output_keys"]}
    output_data = {}
    for key in config["output_keys"]:
        output_all[key] = ""

    for orthoid in config["groups2run"]:
        (
            output_data["gc"],
            output_data["changes"],
            output_data["confirm"],
            output_data["skipped"],
            output_data["conflict"],
        ) = scan_ndata_file(
            orthoid,
            ortho_df=ortho_df,
            prevgc_df=prevgc_df,
            config=config,
            get_sequences_fn=get_sequences,
        )
        for key in config["output_keys"]:
            output_all[key] += output_data[key]

    for key in config["output_keys"]:
        if output_all[key] != "":
            print(headers[key], output_all[key], sep="")
        else:
            eprint("NOTICE: no output for {}".format(key))
else:
    start_secs = time.time()

    if config["threads"] > 1:  # multithread
        for key in config["output_keys"]:
            # remove leftover tempfiles
            clean_up_tempfiles(config[key + "_outfile"])
        eprint(
            "* Working in {} parallel threads; your OS reports {} CPUs.".format(
                config["threads"], cpu_count()
            )
        )

        pool = Pool(config["threads"])
        if config["progressbar"]:
            _ = list(
                tqdm(
                    pool.imap(scan_ndata_file_wrapper, orthogroups),
                    total=len(orthogroups),
                )
            )
        else:
            pool.imap(scan_ndata_file_wrapper, orthogroups)
        pool.close()  # no more work to submit
        pool.join()  # wait workers to finish

        for key in config["output_keys"]:
            combine_and_print_output(
                config[key + "_outfile"],
                key,
                headers[key],
                prevgc_df=prevgc_df,
                config=config,
            )
    else:  # single thread
        output_fh = {}  # output filehandles
        output_data = {}
        for key in config["output_keys"]:
            eprint("* Creating output file: {}".format(config[key + "_outfile"]))
            fh = open(config[key + "_outfile"], "w")
            fh.write(headers[key])
            output_fh[key] = fh

        if config["progressbar"]:
            iterator = tqdm(orthogroups)
        else:
            iterator = orthogroups

        if not prevgc_df.empty:  # prepend the previous suggestions to gc output
            dump_prev_changes(output_fh["gc"], prevgc_df, config=config)

        for orthoid in iterator:
            (
                output_data["gc"],
                output_data["changes"],
                output_data["confirm"],
                output_data["skipped"],
                output_data["conflict"],
            ) = scan_ndata_file(
                orthoid,
                ortho_df=ortho_df,
                prevgc_df=prevgc_df,
                config=config,
                get_sequences_fn=get_sequences,
            )
            for key in config["output_keys"]:
                if output_data[key] != "":
                    output_fh[key].write(output_data[key])

        # close files
        for fh in output_fh.values():
            fh.close()

    eprint(
        "\n*** Part 5 workflow completed {} -- Elapsed: {}, {} g/s --\n".format(
            time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
            *elapsed_time(start_secs, len(orthogroups)),
        )
    )


# #
