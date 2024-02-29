#!/usr/bin/env python
# coding: utf-8
"""configuration for muscle"""
# muscle configuration
MUSCLE_EXE = "muscle3.8"  # path or name of muscle executable

# if using muscle3.8:
MUSCLE_CALL = [MUSCLE_EXE, "-clwstrict", "-quiet", "-maxiters", "2", "-diags"]
ALIGN_FORMAT = "clustal"

# if using muscle5.1:
# MUSCLE_CALL = [MUSCLE_EXE, "-align", "/dev/stdin", "-output", "/dev/stdout"]
# ALIGN_FORMAT = "fasta"
