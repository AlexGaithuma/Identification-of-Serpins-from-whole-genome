#!/usr/bin/env python3
"""
Identify proteins containing user-defined motifs near the C-terminus (e.g., RCL motif search)
with options for motif fasta input and detailed verbose reporting.

Author: Alex Gaithuma
Date: September 3rd 2025

This script:
- Reads motifs from a FASTA file (--motif_file)
- Reads protein sequences from a protein FASTA file (--proteins_fasta)
- Searches each protein's C-terminal window for best motif matches (allowing mismatches)
- Outputs FASTA and GFF of sequences with at least one motif hit

Verbose reporting at every key step!
"""

import argparse
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
import sys

def parse_motifs(fasta):
    """
    Parse a motif FASTA file into a list of motifs.

    Motifs can be on multiple lines, so lines are joined.
    Empty lines and headers (starting with '>') are handled.
    """
    motifs = []
    current = ''
    print(f"[INFO] Reading motifs from: {fasta}")
    try:
        with open(fasta) as fin:
            for line in fin:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('>'):
                    if current:
                        motifs.append(current)
                        current = ''
                else:
                    current += line
            if current:
                motifs.append(current)
    except Exception as e:
        print(f"[ERROR] Failed to read motifs: {e}")
        sys.exit(1)
    print(f"[INFO] Parsed {len(motifs)} motifs from motif file.")
    if len(motifs) == 0:
        print("[WARNING] No motif sequences parsed. Check motif FASTA file!")
    return motifs

def mismatches(seq1, seq2):
    """Count number of non-matching residues between two sequences of equal length."""
    return sum(a != b for a, b in zip(seq1, seq2))

def max_consecutive_mismatches_ok(seq1, seq2, max_consec=1):
    """
    Return True if no more than max_consec consecutive mismatches occur.

    max_consec=1 means: No two (or more) mismatches in a row allowed.
    """
    consec = 0
    for a, b in zip(seq1, seq2):
        if a != b:
            consec += 1
            if consec > max_consec:
                return False
        else:
            consec = 0
    return True

def process_record(args):
    """
    Process a single record for motif search (to be run in parallel).
    Returns best motif hit info or None.
    """
    rec_str, motif_infos, window_len, verbose = args
    import ast
    rec_dict = ast.literal_eval(rec_str)
    seq_id = rec_dict['id']
    seq = rec_dict['seq']
    if len(seq) < window_len:
        search_start = 0
        region = seq
        if verbose: print(f"[WARNING] Sequence {seq_id} shorter than window_len; searching full sequence.")
    else:
        search_start = max(0, len(seq)-window_len)
        region = seq[search_start:]
    region_len = len(region)
    best = None
    for motif, L, max_mis in motif_infos:
        if region_len < L:
            continue
        for i in range(region_len-L+1):
            subseq = region[i:i+L]
            mis = mismatches(subseq, motif)
            if mis <= max_mis:
                if not max_consecutive_mismatches_ok(subseq, motif, max_consec=1):
                    continue
                hit_start = search_start + i + 1  # 1-based
                hit_end = hit_start + L - 1
                if (best is None) or (mis < best[3]) or (mis == best[3] and hit_start < best[0]):
                    best = (hit_start, hit_end, L, mis)
    if best:
        if verbose:
            print(f"[MATCH] {seq_id}: Best motif found at {best[0]}-{best[1]} (Length {best[2]}, {best[3]} mismatches)")
        return (seq_id, rec_dict, best)
    else:
        if verbose:
            print(f"[INFO] {seq_id}: No motif found in search region.")
        return None

def main():
    parser = argparse.ArgumentParser(
        description="Identify proteins with C-terminal motif hits (e.g., RCL motif)."
    )
    parser.add_argument("--proteins_fasta", required=True, help="Input proteins fasta file (required)")
    parser.add_argument("--motif_file", required=True, help="Motif (or RCL motif) fasta file (required)")
    parser.add_argument("--prefix", required=True, help="Prefix for output files")
    parser.add_argument("--window_len", type=int, default=100, help="C-terminal window length for motif search (default: 100)")
    parser.add_argument("--verbose", action="store_true", help="Print detailed progress and matching info")

    args = parser.parse_args()

    proteins_fasta = args.proteins_fasta
    motif_file = args.motif_file
    output_fasta = f"{args.prefix}_Identified_serpins.fasta"
    output_gff   = f"{args.prefix}_Identified_serpins.gff"
    window_len = args.window_len
    verbose = args.verbose

    print(f"[INFO] Starting motif search on protein FASTA: {proteins_fasta}")
    print(f"[INFO] Using motif FASTA: {motif_file}")
    print(f"[INFO] Output FASTA: {output_fasta}")
    print(f"[INFO] Output GFF: {output_gff}")
    print(f"[INFO] Window size for C-terminal search: {window_len}")

    # Step 1: Parse motifs
    motifs = parse_motifs(motif_file)
    # Only keep motifs of reasonable length (10 to 30); set max mismatch at 30%
    motif_infos = []
    for motif in motifs:
        L = len(motif)
        if 10 <= L <= 30:
            max_mis = int(L * 0.3)
            motif_infos.append((motif, L, max_mis))
        else:
            print(f"[WARNING] Skipping motif (length {L}) outside 10-30 AA: {motif[:30]}...")

    if not motif_infos:
        print("[ERROR] No suitable motifs for analysis after filtering on length. Exiting.")
        sys.exit(2)
    print(f"[INFO] Using {len(motif_infos)} motifs for search.")

    # Step 2: Read proteins and store as list for parallel processing
    protein_records = []
    print(f"[INFO] Reading protein FASTA: {proteins_fasta}")
    for rec in SeqIO.parse(proteins_fasta, "fasta"):
        protein_records.append(str({'id': rec.id, 'desc': rec.description, 'seq': str(rec.seq)}))
    total_records = len(protein_records)
    print(f"[INFO] Loaded {total_records} protein records.")

    if total_records == 0:
        print("[ERROR] No protein records found in input. Exiting.")
        sys.exit(3)

    # Step 3: Parallel motif search (multiprocessing)
    hits = {}
    records_with_hits = {}
    num_workers = min(98, multiprocessing.cpu_count(), 80)
    print(f"[INFO] Using {num_workers} worker processes for motif search.")

    args_for_pool = [
        (rec_str, motif_infos, window_len, verbose)
        for rec_str in protein_records
    ]

    # Use ProcessPoolExecutor for parallel record searching
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(process_record, tuple_args) for tuple_args in args_for_pool]
        completed = 0
        for future in as_completed(futures):
            completed += 1
            if not verbose and (completed % 1000 == 0 or completed == total_records):
                print(f"[INFO] Processed {completed}/{total_records} records...", end="\r", file=sys.stderr)
            result = future.result()
            if result:
                seq_id, rec_dict, best = result
                hits[seq_id] = best
                records_with_hits[seq_id] = rec_dict
    print()
    print(f"[INFO] Finished processing. Hits found in {len(hits)} of {total_records} proteins.")

    # Step 4: Write output FASTA and GFF files
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq

    print(f"[INFO] Writing {len(records_with_hits)} hit protein(s) to FASTA: {output_fasta}")
    with open(output_fasta, "w") as outfa:
        for rec_dict in records_with_hits.values():
            record = SeqRecord(Seq(rec_dict['seq']), id=rec_dict['id'], description=rec_dict['desc'])
            SeqIO.write(record, outfa, "fasta")

    print(f"[INFO] Writing GFF output (hit coordinates, mismatches, annotation) to: {output_gff}")
    with open(output_gff, "w") as gff:
        gff.write("##gff-version 3\n")
        for sid in hits:
            s, e, L, mis = hits[sid]
            gff.write(f"{sid}\tmotif_search\tMotif\t{s}\t{e}\t.\t+\t.\tNote=MotifHit;Length={L};Mismatches={mis}\n")

    print("[INFO] Done! Check output files.")

if __name__ == "__main__":
    main()