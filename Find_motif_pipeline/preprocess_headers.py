#!/usr/bin/env python3
"""
===============================================================================
Script: preprocess_headers.py
Purpose: Clean FASTA headers for uniformity and traceability.

NEW RULE:
    - For each sequence name, rename to everything before the second '_'
      (eg 's1c1_g1610.t1_Chr1 XP_xxx stuff' → 's1c1_g1610.t1').
    - Applies to all FASTA files: full sequences, alignment, N-terminal and RCL motifs.
      This ensures robust cross-matching for downstream motif context extraction.

Usage:
  python preprocess_headers.py <input_fasta> <output_fasta>
===============================================================================
"""

import sys
from Bio import SeqIO

if len(sys.argv) != 3:
    sys.exit("Usage: python preprocess_headers.py <input_fasta> <output_fasta>")

in_fasta, out_fasta = sys.argv[1:3]

def rename_id(header):
    """
    Given a header line, returns the part up to the second '_'
    Example: s1c1_g1610.t1_Chr1 XP_xxx stuff => s1c1_g1610.t1
    """
    parts = header.split('_', 2) # limit splits in case there are more underscores
    if len(parts) >= 2:
        return parts[0] + '_' + parts[1].split()[0]  # In case 2nd _ is just before carriage/space
    else:
        return header.split()[0]

with open(out_fasta, "w") as out_f:
    for rec in SeqIO.parse(in_fasta, "fasta"):
        rec.id = rename_id(rec.id)
        rec.name = rec.id
        rec.description = ""  # No extra info
        SeqIO.write(rec, out_f, "fasta")