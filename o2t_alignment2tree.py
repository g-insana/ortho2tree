#!/usr/bin/env python
# coding: utf-8
"""module providing function alignment2tree() to create a tree from an alignment"""
import re
import sys
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Align import MultipleSeqAlignment


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


def alignment2tree(aln, dist_only=False, better=False, verbose=False):
    pdist_calc = DistanceCalculator("pam30")
    keep_set = {}
    for entry in aln:
        # entry_id = entry.id #UNUSED
        bt_add2set(bt_extract_acc(entry.id), keep_set, True)

    pdist = pdist_calc.get_distance(aln)

    # make a dummy holder for gap-based distances
    gdist = pdist_calc.get_distance(aln)

    for ix, _ in enumerate(aln):
        gdist.matrix[ix][ix] = 0.0
        for jx in range(ix):
            gdist.matrix[ix][jx] = bt_gap_dist(aln[ix], aln[jx])

    max_name_len = 0

    if dist_only and not better:
        for name in pdist.names:
            if len(name) > max_name_len:
                max_name_len = len(name)

        print("gap distances:")
        for ix, name in enumerate(gdist.names):
            print(" % -*s" % (max_name_len, name), end="")
            for jx in range(ix + 1):
                print("  %.4f" % (gdist.matrix[ix][jx]), end="")
            print()
        print()

    # if better -- find the worst canonical distances, exclude isoforms with worse differences
    # (remove alignments worse than canonical)
    worst_name = ""
    worst_jname = ""
    worst_canon_dist = 0.0
    for ix, name in enumerate(gdist.names):
        i_acc = bt_extract_acc(name)
        for jx in range(ix + 1):
            j_name = gdist.names[jx]
            j_acc = bt_extract_acc(j_name)
            if j_acc in keep_set and i_acc in keep_set:
                if gdist.matrix[ix][jx] > worst_canon_dist:
                    worst_canon_dist = gdist.matrix[ix][jx]
                    worst_name = name
                    worst_jname = j_name

    if better:
        sys.stderr.write(
            "worst: [%s][%s]: %.3f\n" % (worst_name, worst_jname, worst_canon_dist)
        )

        good_set = {}
        good_list = []
        good_ind = []

        if verbose:
            sys.stderr.write("keep_set:\n")
            for acc in keep_set:
                sys.stderr.write(acc + "\n")

        ndim = len(gdist.names)
        for ix, name in enumerate(gdist.names):
            i_name_acc = bt_extract_acc(name)

            if i_name_acc not in keep_set:
                is_good = False
                min_cost = 1000.0

                # go across matrix to diagonal
                for jx in range(ix):
                    j_name = gdist.names[jx]
                    j_name_acc = bt_extract_acc(j_name)
                    if j_name_acc not in keep_set:
                        continue
                    this_cost = gdist.matrix[ix][jx]
                    if this_cost < min_cost:
                        min_cost = this_cost
                    if verbose:
                        sys.stderr.write(
                            " % s v %s [%s]: %.3f < %.3f "
                            % (
                                i_name_acc,
                                j_name_acc,
                                (j_name_acc in keep_set),
                                this_cost,
                                worst_canon_dist,
                            )
                        )
                    if j_name_acc in keep_set and this_cost < worst_canon_dist:
                        is_good = True
                        if verbose:
                            sys.stderr.write("True\n")
                    else:
                        if verbose:
                            sys.stderr.write("\n")

                # go down matrix from diagonal
                for jx in range(ix + 1, ndim):
                    j_name_acc = bt_extract_acc(gdist.names[jx])
                    if j_name_acc not in keep_set:
                        continue
                    this_cost = gdist.matrix[jx][ix]
                    if this_cost < min_cost:
                        min_cost = this_cost
                    if verbose:
                        sys.stderr.write(
                            " % s v %s [%s]: %.3f < %.3f "
                            % (
                                i_name_acc,
                                j_name_acc,
                                (j_name_acc in keep_set),
                                this_cost,
                                worst_canon_dist,
                            )
                        )
                    if j_name_acc in keep_set and this_cost < worst_canon_dist:
                        is_good = True
                        if verbose:
                            sys.stderr.write("True\n")
                    else:
                        if verbose:
                            sys.stderr.write("\n")

                if is_good:
                    good_set[name] = True
                    good_list.append(name)
                else:
                    sys.stderr.write(
                        " %s excluding: %s  %.4f\n" % (sys.argv[0], name, min_cost)
                    )

            else:
                good_set[name] = True
                good_list.append(name)

        good_matrix = []

        good_align = MultipleSeqAlignment([])

        good_cnt = 0
        for ix, name in enumerate(gdist.names):
            if name not in good_set:
                continue

            good_cnt += 1

            good_row = []

            good_align.add_sequence(aln[ix].id, str(aln[ix].seq))

            for jx in range(ix + 1):
                if gdist.names[jx] not in good_set:
                    continue
                good_row.append(gdist.matrix[ix][jx])

            good_matrix.append(good_row)
            good_ind.append(ix)

        sys.stderr.write(
            "** good cnt: %d; len(good_ind): %d\n" % (good_cnt, len(good_ind))
        )

        if dist_only:
            for name in good_list:
                if len(name) > max_name_len:
                    max_name_len = len(name)

            for ix, name in enumerate(good_list):
                print(" % -*s" % (max_name_len, name), end="")
                for jx in range(ix + 1):
                    print("  %.4f" % (good_matrix[ix][jx]), end="")
                print()

            print()

    constructor = DistanceTreeConstructor()

    njtree_g = constructor.nj(gdist)
    if worst_canon_dist > 0.0:
        njtree_g.root_at_midpoint()

    if not dist_only:
        return njtree_g
