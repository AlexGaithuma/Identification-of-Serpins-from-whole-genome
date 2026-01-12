#!/usr/bin/env python3
"""
===============================================================================
Script: build_motif_matrix.py
Purpose: Build a matrix (sequence x motif/context) for motif occurrence using FIMO.

Biological basis:
 - Each row: sequence; Each column: motif "context:motifid" (from MEME discovered motifs,
   using `meme.txt` as the per-context motif file)
 - Matrix values: 1 if motif detected in sequence, 0 else; columns can be weighted later.

Pipeline requirements and context:
 - For each context-specific directory (e.g. RCL_context), expects `meme.txt` in MEME motif
   format, as output by MEME Suite.
 - All motifs in every `meme.txt` for every context are searched with FIMO over the input
   sequences, and the results are collated into a presence/absence matrix across all
   contexts+motifs (uniquely labeled as context:motif).
 - The result is a "fingerprint" matrix for clustering and downstream analysis.

===============================================================================
"""

import os
import glob
import argparse
import sys
import pandas as pd
from io import StringIO
from Bio import SeqIO
import subprocess

def run_fimo_on_meme_motif(meme_file, seq_fasta):
    """
    Run FIMO for given MEME-format motif file and input sequence FASTA.
    Returns a DataFrame with motif hits, or empty DataFrame if no hits.
    """
    # FIMO expects MEME format for motifs.
    cmd = [
        "fimo",
        "--verbosity", "1",
        "--text",
        meme_file,
        seq_fasta
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[WARNING] FIMO failed for {meme_file} on {seq_fasta}: {result.stderr}", file=sys.stderr)
        return pd.DataFrame()
    # Remove comment lines, parse as tsv
    lines = [l for l in result.stdout.splitlines() if l.strip() and not l.startswith("#")]
    if len(lines) < 2:
        return pd.DataFrame()  # No real hits
    return pd.read_csv(StringIO('\n'.join(lines)), sep='\t')

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--context_dir", required=True,
        help="Directory with per-context subdirs each containing a meme.txt motif file (MEME Suite format)")
    parser.add_argument("--fasta", required=True,
        help="FASTA file with sequences/samples to scan for motif matches")
    parser.add_argument("--out", required=True,
        help="Tab-delimited output matrix; rows=sequences, columns=context:motif_id, values=1/0 for hit")
    args = parser.parse_args()

    # Discover (non-empty) meme.txt motif files in all context subdirs
    context_dirs = sorted(glob.glob(os.path.join(args.context_dir, "*")))
    context_meme_tuples = []
    for ctxdir in context_dirs:
        meme_file = os.path.join(ctxdir, "meme.txt")
        if os.path.isfile(meme_file) and os.path.getsize(meme_file) > 0:
            ctxname = os.path.basename(ctxdir)
            context_meme_tuples.append((ctxname, meme_file))

    if not context_meme_tuples:
        sys.exit("[ERROR] No motif files found (no meme.txt present or not non-empty) in any context subdir!")

    print(f"[INFO] Found {len(context_meme_tuples)} context motif files.")

    # Get all target sequence IDs up front
    seqs = [rec.id for rec in SeqIO.parse(args.fasta, 'fasta')]
    if not seqs:
        sys.exit(f"[ERROR] No sequence IDs found in your FASTA ({args.fasta})")

    # For each context, scan with FIMO and aggregate all motif hits.
    motif_columns = set()
    presence = []  # will be: list of (seqid, colname)
    for ctxname, meme_file in context_meme_tuples:
        print(f"[INFO] Scanning context: {ctxname}")
        fimo_hits = run_fimo_on_meme_motif(meme_file, args.fasta)
        if fimo_hits.empty:
            print(f"  [INFO] No motif hits in context {ctxname}.")
            continue
        # Use 'motif_id' for unique motif; 'sequence_name' for sample/sequence
        if 'motif_id' in fimo_hits.columns and 'sequence_name' in fimo_hits.columns:
            for _, row in fimo_hits.iterrows():
                motif = f"{ctxname}:{row['motif_id']}"
                motif_columns.add(motif)
                presence.append((row['sequence_name'], motif))
        else:
            print(f"  [WARNING] FIMO output columns missing 'motif_id' or 'sequence_name' for {ctxname}.")

    motif_columns = sorted(motif_columns)
    if not motif_columns:
        sys.exit("[ERROR] No motif features detected by FIMO for any context. "
                 "Check your meme.txt content and FIMO compatibility.")

    # Build presence/absence matrix
    mat = pd.DataFrame(0, index=seqs, columns=motif_columns)
    for seqid, motifcol in presence:
        if seqid in mat.index and motifcol in mat.columns:
            mat.loc[seqid, motifcol] = 1

    mat.index.name = "SeqID"
    mat.to_csv(args.out, sep="\t")
    print(f"[SUCCESS] Wrote motif matrix: {args.out} (shape={mat.shape})")