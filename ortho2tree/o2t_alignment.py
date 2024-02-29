#!/usr/bin/env python
# coding: utf-8
"""
module providing function get_alignment() which
makes exteral call to muscle in order to align sequences
"""
import subprocess
from Bio import AlignIO
from .o2t_utils import eprint
from .config_muscle import MUSCLE_CALL, ALIGN_FORMAT


def get_alignment(fasta_sequences, textoutput=False):
    """
    argument: a string containing fasta sequences
    returns:
        if textoutput: string with alignment in aligned fasta format
        otherwise: Bio.Align.MultipleSeqAlignment instance
    """
    with subprocess.Popen(
        MUSCLE_CALL,
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    ) as p:
        if textoutput:
            align, error = p.communicate(
                input=fasta_sequences
            )  # returns stdout, stderr
            if error != "":
                eprint(error)
        else:
            p.stdin.write(fasta_sequences)
            p.stdin.close()
            if p.returncode is not None:
                eprint("Error! return code: {}".format(p.returncode))
                return None
            try:
                # get an align object
                align = AlignIO.read(p.stdout, ALIGN_FORMAT)
            except ValueError:
                eprint("Error with alignment")
                return None
    return align
