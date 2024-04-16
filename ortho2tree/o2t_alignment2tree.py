#!/usr/bin/env python
# coding: utf-8
"""module providing function alignment2tree() to create a tree from an alignment"""
import itertools
import re
import sys

# import os  # timing
# import time  # timing

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from .o2t_utils import eprint

# PDIST_CALC = DistanceCalculator("pam30")
CONSTRUCTOR = DistanceTreeConstructor()
# Note: uncomment all the code lines marked #timing to perform analyis of pipeline times


def _matrixcompare(m1, m2):
    """
    comparison of two distance matrices
    """
    total_zero = True
    max_delta = 0.0
    for i in range(0, len(m1.matrix)):
        del_line = [
            "%.5f" % (m1.matrix[i][j] - m2.matrix[i][j])
            for j in range(0, len(m1.matrix[i]))
        ]
        all_zero = True
        for j in range(0, i):
            delta = m1.matrix[i][j] - m2.matrix[i][j]
            if delta > 2.0 * sys.float_info.epsilon:
                all_zero = False
                total_zero = False
            if delta > max_delta:
                max_delta = delta

        if not all_zero:
            eprint(" ".join(del_line))
    if not total_zero:
        eprint("some non-zero deltas: %g" % (max_delta))
        return False
    return True


def _pairwise(seq1, seq2, g_thresh):
    """
    Calculate pairwise distance from two sequences
    Returns a value between 0 (identical sequences) and 1 (completely different, or seq1 is an empty string.)
    Only consider presence of gaps for scoring
    If parameter g_thresh is passed, only gaps longer than that threshold will be considered for the cost.
    """
    score = alen = gap_len = 0
    in_gap = False

    for l1, l2 in zip(seq1, seq2):
        if l1 != l2:  # they differ
            alen += 1
            if l1 == "-" or l2 == "-":
                in_gap = True
                gap_len += 1
        else:  # they are the same
            if l1 != "-":  # but not both gap
                alen += 1
                if in_gap:
                    if gap_len > g_thresh:
                        score += gap_len
                    in_gap = False
                    gap_len = 0

    if in_gap and gap_len > g_thresh:
        score += gap_len

    if alen > 0:
        return float(score) / alen
    else:
        return 1.0


def _get_distance(msa, g_thresh):
    """
    given a Bio.Align.MultipleSeqAlignment return a DistanceMatrix
    Only consider presence of gaps for scoring
    Optionally only considers gaps longer than specified g_thresh
    NOTE: barebone version of Bio.Phylo.TreeConstruction.DistanceCalculator.get_distance
    """
    names = [s.id for s in msa]
    dm = DistanceMatrix(names)
    for seq1, seq2 in itertools.combinations(msa, 2):
        dm[seq1.id, seq2.id] = _pairwise(seq1, seq2, g_thresh)
    return dm


def bt_gap_dist(s1, s2):
    if len(s1) != len(s2):
        sys.stderr.write(
            "*** simple_dist lengths differ %s: %d :: %s: %d\n"
            % (s1.id, len(s1), s2.id, len(s2))
        )
        return 1.0

    score = 0.0
    alen = 0
    for ix, s1c in enumerate(s1):
        s2c = s2[ix]
        if s1c == s2c:
            if s1c != "-":
                score += 1.0
                alen += 1
        else:
            alen += 1
            if s1c != "-" and s2c != "-":
                score += 1.0

    if alen != 0:
        # prevent division by zero
        d_score = 1.0 - score / alen
    else:
        d_score = 1

    return d_score


def bt_add2set(entry, entry_set, value):
    if entry not in entry_set:
        entry_set[entry] = value


def bt_extract_acc(entry):
    acc = entry
    entry_fields = entry.split("|")
    if len(entry_fields) == 3:
        acc = entry_fields[1]
    elif len(entry_fields) == 2:
        acc = entry_fields[0]
    elif len(entry_fields) == 1:
        acc = entry_fields[0]
    else:
        sys.stderr.write("*** ERROR *** unrecognized entry: %s\n" % (entry))

    # remove version number if present
    acc = re.sub(r"\.\d+$", "", acc)
    return acc


def alignment2tree(aln, orthoid="", verbose=False, config=None):
    # semaphore_file = os.path.join(config["semaphores_dir"], orthoid + ".done")  # timing
    g_thresh = config.get(
        "gap_threshold", 0
    )  # to only consider gaps longer than this threshold

    # distance matrix using gap-based distances
    # gdist_ori = PDIST_CALC.get_distance(aln)
    # for ix, _ in enumerate(aln):
    #    gdist_ori.matrix[ix][ix] = 0.0
    #    for jx in range(ix):
    #        gdist_ori.matrix[ix][jx] = bt_gap_dist(aln[ix], aln[jx])

    # process_start_time = time.time()  # timing
    gdist = _get_distance(aln, g_thresh)
    # process_name = "t1"  # getdist new internal code #timing
    # process_duration = time.time() - process_start_time  # timing
    # with open(semaphore_file, "a") as fh:  # timing
    #   fh.write("{}\t{}\n".format(process_name, process_duration))  # timing

    # process_start_time = time.time() # timing
    njtree_g = CONSTRUCTOR.nj(gdist)
    # process_name = "t3"  # constructor #timing
    # process_duration = time.time() - process_start_time  # timing
    # with open(semaphore_file, "a") as fh:  # timing
    #   fh.write("{}\t{}\n".format(process_name, process_duration))  # timing

    # here we need to set all negative branch lengths to zero (it
    # would be better to preserve distances, but zero is a start)
    for clade in njtree_g.find_clades():
        if (clade.branch_length < 0.0):
            clade.branch_length = 0.0

    # process_start_time = time.time()  # timing
    try:  # if it fails, it may mean worst_canon_dist = 0
        njtree_g.root_at_midpoint()
    except Exception as e:
        eprint(f"The following error occurred: {e}")
        # find the worst canonical distances, exclude isoforms with worse differences
        # (remove alignments worse than canonical)
        keep_set = {}
        for entry in aln:
            bt_add2set(bt_extract_acc(entry.id), keep_set, True)

        worst_canon_dist = 0.0
        for ix, name in enumerate(gdist.names):
            i_acc = bt_extract_acc(name)
            for jx in range(ix + 1):
                j_name = gdist.names[jx]
                j_acc = bt_extract_acc(j_name)
                if j_acc in keep_set and i_acc in keep_set:
                    if gdist.matrix[ix][jx] > worst_canon_dist:
                        worst_canon_dist = gdist.matrix[ix][jx]

        if worst_canon_dist == 0.0:
            # we are done: the whole tree can be written as single clade
            # to n_data and it is an output_confirm
            eprint(
                "ALN2TREE: write n_data for {} as worst_canon_dist is zero".format(
                    orthoid
                )
            )
            return njtree_g
        else:
            # maybe it's a bad tree, warn as a possible output_skipped
            eprint(
                "ALN2TREE: {} to be skipped? worst_canon_dist is {}".format(
                    orthoid, worst_canon_dist
                )
            )
            return None
    # process_name = "t4"  # rooting #timing
    # process_duration = time.time() - process_start_time  # timing
    # with open(semaphore_file, "a") as fh:  # timing
    #   fh.write("{}\t{}\n".format(process_name, process_duration))  # timing

    return njtree_g
