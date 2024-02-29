#!/usr/bin/env python
# coding: utf-8
"""
module providing functions related with identification of outliers
"""
import numpy as np
from .o2t_utils import eprint


def _detect_outliers_std(data, threshold=3):
    """
    NOT USED
    detect outliers in a collection based on mean&std
    outliers are identified as those beyond /threshold/ standard deviations from mean
    default threshold is 3 (3 std)
    set, list or dict can all be given as argument
    if dict is given, values are checked for outliers and corresponding keys are returned
    return also mean and std

    examples:
      detect_outliers([3, 5, 6, 8, 9, 1000], threshold = 2)
      detect_outliers({'a':3, 'b':5, 'd': 6, 'c':1000}, threshold = 1)
    """
    if isinstance(data, dict):
        datavalues = list(data.values())
    else:
        datavalues = data
    # eprint(datavalues)
    mean = np.mean(datavalues)
    std = np.std(datavalues)
    outliers = []
    for x in datavalues:
        z_score = (x - mean) / std
        if np.abs(z_score) > threshold:
            outliers.append(x)
    if isinstance(data, dict):
        return [k for (k, v) in data.items() if v in set(outliers)], mean, std
    return outliers, mean, std


def detect_outliers_median(data, threshold_lo=0.75, threshold_hi=1.25):
    """
    detect outliers in a collection based on median
    outliers are identified as those under /threshold_lo/ times the median and
    those above /threshold_hi/ times the median
    default thresholds are 0.75 and 1.25
    set, list or dict can all be given as argument
    if dict is given, values are checked for outliers and corresponding keys are returned
    returns also median value

    examples:
      detect_outliers_median([3, 5, 6, 8, 9, 1000], threshold_lo=0.1, threshold_hi=1.9)
      detect_outliers_median([30, 27, 23, 24, 33])
      detect_outliers_median({'a':7, 'b':5, 'd': 6, 'c':1000})
    """
    if isinstance(data, dict):
        datavalues = list(data.values())
    else:
        datavalues = data
    # eprint(datavalues)
    median = np.median(datavalues)
    outliers = []
    for x in datavalues:
        if x < threshold_lo * median or x > threshold_hi * median:
            outliers.append(x)
    if isinstance(data, dict):
        return [k for (k, v) in data.items() if v in set(outliers)], median
    return outliers, median


def flag_maxseqlen(df, max_seqlen):
    """
    mark as outliers all those with sequence length beyond a threshold
    in case of canonicals, mark also their isoforms as outliers
    """
    # mark as outlier entries > max_seqlen
    seqs_too_long = df[df["seqlen"] > max_seqlen].index
    eprint(
        "Marking as outliers {} sequences whose length is > {}".format(
            len(seqs_too_long), max_seqlen
        )
    )
    df.loc[seqs_too_long, "outlier"] = True
    # identify canonicals marked as outliers and making their isoforms outliers too
    groups_cano_outlier = df[df["outlier"] & df["is_canonical"]]["groupid"].to_list()
    eprint(
        "Identified {} outlier canonicals; marking their isoforms as outliers too. GC group ids involved: {}".format(
            len(groups_cano_outlier), groups_cano_outlier
        )
    )
    df.loc[df["groupid"].isin(groups_cano_outlier), "outlier"] = True


def flag_outliers_in_df(df, config=None):
    """
    detect outliers and mark them as such in the dataframe
    """

    orthogroups = df.groupby("pantherid")

    if (
        config["detect_outliers_with_median"]
        or config["detect_outliers_with_median_and_can_lengths"]
    ):
        seqlen_medians = orthogroups.seqlen.transform(
            "median"
        )  # find median of sequence lengths in orthogroup
        df["outlier_by_median"] = np.where(
            df[
                "is_canonical"
            ],  # create column flagging outliers below or beyond specified thresholds times median sequence length of orthogroup sequences
            False,  # canonical never outlier
            ~(
                df.seqlen.between(
                    seqlen_medians * config["outliers_detection_threshold_median_lo"],
                    seqlen_medians * config["outliers_detection_threshold_median_hi"],
                )
            ),
        )

    if (
        config["detect_outliers_with_quart"]
        or config["detect_outliers_with_quart_and_can_lengths"]
    ):  # for quartile based
        seqlen_q1 = orthogroups.seqlen.transform("quantile", 0.25)  # 1st quartile
        seqlen_q3 = orthogroups.seqlen.transform("quantile", 0.75)  # 3rd quartile
        df["outlier_by_quart"] = np.where(
            df[
                "is_canonical"
            ],  # create column flagging outliers below or beyond specified thresholds times q1 and q3 of sequence lengths of orthogroup sequences
            False,  # canonical never outlier
            ~(
                df.seqlen.between(
                    seqlen_q1 * config["outliers_detection_threshold_quart_lo"],
                    seqlen_q3 * config["outliers_detection_threshold_quart_hi"],
                )
            ),
        )

    if (
        config["detect_outliers_with_median_and_can_lengths"]
        or config["detect_outliers_with_quart_and_can_lengths"]
    ):  # for quartile based
        orthogroups_canonicals = df[df["is_canonical"]].groupby(
            "pantherid"
        )  # canonicals grouped by orthogroup
        max_can_len_dict = (
            orthogroups_canonicals.seqlen.max()
        )  # find max sequence length of canonicals in each orthogroup
        min_can_len_dict = (
            orthogroups_canonicals.seqlen.min()
        )  # find max sequence length of canonicals in each orthogroup
        df["min_can_len"] = df["pantherid"].map(
            min_can_len_dict, na_action="ignore"
        )  # create column with max sequence length of canonicals/orthogroup
        df["max_can_len"] = df["pantherid"].map(
            max_can_len_dict, na_action="ignore"
        )  # create column with min sequence length of canonicals/orthogroup
        df["outlier_by_can_len"] = np.where(
            df[
                "is_canonical"
            ],  # create column flagging outliers below or beyond specified thresholds times min and max sequence length of orthogroup canonicals
            False,  # canonical never outlier
            ~(
                df.seqlen.between(
                    df.min_can_len * config["outliers_detection_threshold_can_lo"],
                    df.max_can_len * config["outliers_detection_threshold_can_hi"],
                )
            ),
        )

    if config["detect_outliers_with_mean_std"]:
        seqlen_means = orthogroups.seqlen.transform("mean")
        seqlen_stds = orthogroups.seqlen.transform("std")

        df["outlier"] = np.where(
            df["is_canonical"],
            False,  # canonical never outlier
            ~(
                df.seqlen.between(
                    seqlen_means
                    - seqlen_stds * config["outliers_detection_threshold_std_lo"],
                    seqlen_means
                    + seqlen_stds * config["outliers_detection_threshold_std_hi"],
                )
            ),
        )
    elif config["detect_outliers_with_median"]:
        df.rename(columns={"outlier_by_median": "outlier"}, inplace=True)
    elif config["detect_outliers_with_quart"]:
        df.rename(columns={"outlier_by_quart": "outlier"}, inplace=True)
    elif config["detect_outliers_with_median_and_can_lengths"]:
        df["outlier"] = np.where(
            df["is_canonical"], False, df.outlier_by_median & df.outlier_by_can_len
        )  # flag as outliers only those that satisfy both conditions
    elif config["detect_outliers_with_quart_and_can_lengths"]:
        df["outlier"] = np.where(
            df["is_canonical"], False, df.outlier_by_quart & df.outlier_by_can_len
        )  # flag as outliers only those that satisfy both conditions
    else:
        eprint("Notice: no outlier detection will be performed")
        df["outlier"] = False
