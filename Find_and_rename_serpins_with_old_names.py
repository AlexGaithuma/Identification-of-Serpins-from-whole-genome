#!/usr/bin/env python3
"""
Motif-to-Protein Unique Assignment Script
=========================================
* Each motif from --in_motifs FASTA is assigned to its best-matching protein in --in_prot,
  using mismatches or 1 allowed gap.
* Each protein is represented ONCE in the output, and is assigned to the motif where it 
  is the best hit (one-to-one assignment).
* If >2 mismatches, use -like postfix in name.
* Any unassigned motif after all processed is forcibly matched to a not-yet-assigned protein,
  or (if all used) to the best of all proteins, with -like if needed.
* Outputs renamed FASTA, detailed GTF, and RCL.log with all verbose explanations and summaries.

USAGE:
    python motif_unique_assignment.py \
        --in_prot proteins.fasta --in_motifs motifs.fasta \
        --out_protein renamed.fa --out_gtf motifs.gtf

Author: Alex Gaithuma
Date: September 3rd 2025
"""

import sys
import argparse

# ------------------ Step 1: FASTA Parsing Utility ---------------------------------
def parse_fasta(fpath):
    """
    Parse a multi-FASTA file.

    Yields:
        (header(str), sequence(str))
    """
    h, seq = None, []
    with open(fpath) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if h:
                    yield (h, "".join(seq))
                h, seq = line, []
            else:
                seq.append(line)
        if h:
            yield (h, "".join(seq))

# ------------------ Step 2: Hamming Distance --------------------------------------
def hamming(seq1, seq2):
    """Count mismatched positions (strings must be same length)."""
    return sum(a != b for a, b in zip(seq1, seq2))

# ------------------ Step 3: Motif scan - mismatches only --------------------------
def scan_motifs(seq, motif_seq, max_mismatches):
    """
    Scan sliding windows of motif size in seq, allowing up to max_mismatches.
    Returns: [(mismatches, start, end, matched_window)]
    """
    mlen, out = len(motif_seq), []
    for i in range(len(seq)-mlen+1):
        window = seq[i:i+mlen]
        mismatches = hamming(window, motif_seq)
        if mismatches <= max_mismatches:
            out.append((mismatches, i+1, i+mlen, window))
    return out

# ------------------ Step 4: Motif scan with single gap allowed ---------------------
def scan_motifs_gap(seq, motif_seq, max_mismatches):
    """
    Allow a single gap (deletion) in motif or sequence. 
    Returns: [(mismatches,start,end,window,gap_pos,gap_type)]
    """
    mlen, res = len(motif_seq), []
    # (1) Gap in motif: remove each motif position
    for gap in range(mlen):
        motif_gapped = motif_seq[:gap] + motif_seq[gap+1:]
        winlen = mlen - 1
        for i in range(len(seq) - winlen + 1):
            window = seq[i:i+winlen]
            mismatches = hamming(window, motif_gapped)
            if mismatches <= max_mismatches:
                res.append((mismatches, i+1, i+winlen, window, gap+1, "motif_gap"))
    # (2) Gap in sequence: remove each pos from window
    for gap in range(mlen):
        winlen = mlen
        for i in range(len(seq) - winlen + 1):
            window_full = seq[i:i+winlen]
            if len(window_full) < winlen: continue
            window_gapped = window_full[:gap] + window_full[gap+1:]
            if len(window_gapped) != mlen-1: continue
            mismatches = hamming(window_gapped, motif_seq[:gap] + motif_seq[gap+1:])
            if mismatches <= max_mismatches:
                res.append((mismatches, i+1, i+winlen-1, window_gapped, gap+1, "seq_gap"))
    return res

# --------------------------- Step 5: Main Logic -----------------------------------
def main():
    parser = argparse.ArgumentParser(
        description='One-to-one motif/protein matching; all motifs represented; -like for >2 mismatches; protein appears only once in output.'
    )
    parser.add_argument('--in_prot', required=True, help="Input protein sequences (FASTA)")
    parser.add_argument('--in_motifs', required=True, help="Motif sequences (FASTA)")
    parser.add_argument('--out_protein', required=True, help="Renamed protein sequences (FASTA)")
    parser.add_argument('--out_gtf', required=True, help="GTF features for motif assignments")
    parser.add_argument('--max_mismatches', type=int, default=2, help="Max mismatches per match window (default: 2)")
    parser.add_argument('--verbose', action='store_true', help="All output in RCL.log (always on)")
    args = parser.parse_args()

    ### --- Step 6: Load motifs and proteins --- ###
    motifs = []
    motif_seqs = {}
    for h, seq in parse_fasta(args.in_motifs):
        mname = h[1:].split()[0]
        motifs.append(mname)
        motif_seqs[mname] = seq

    protein_ids = []
    protein_headers = {}
    protein_seqs = {}
    for h, seq in parse_fasta(args.in_prot):
        pid = h[1:].split()[0]
        protein_ids.append(pid)
        protein_headers[pid] = h
        protein_seqs[pid] = seq

    # Open outputs
    fout = open(args.out_protein, "w")
    fgtf = open(args.out_gtf, "w")
    log = open("RCL.log", "w")

    #### --- Step 7: Find best motif matches for every protein === ###
    log.write("# Phase 1: Collect all motif/protein best-hits and evidence\n")
    motif2matches = {}  # motif -> [ (protein, match_type, mismatches, start, end, window, gapinfo) ]
    prot2motif_matchinfo = {}  # protein -> list of all motif matches, chosen below

    for mname in motifs:
        motif_seq = motif_seqs[mname]
        matches = []
        for pid in protein_ids:
            seq = protein_seqs[pid]
            # 1) Mismatch hits
            for mismatches, start, end, window in scan_motifs(seq, motif_seq, args.max_mismatches):
                matches.append((pid, "mismatch", mismatches, start, end, window, None))
            # 2) Gap hits (0 mismatch)
            for mismatches, start, end, window, gap_pos, gap_type in scan_motifs_gap(seq, motif_seq, max_mismatches=0):
                matches.append((pid, "gap0", mismatches, start, end, window, (gap_pos, gap_type)))
            # 3) Gap hits (+mismatches)
            for mismatches, start, end, window, gap_pos, gap_type in scan_motifs_gap(seq, motif_seq, args.max_mismatches):
                matches.append((pid, "gapmm", mismatches, start, end, window, (gap_pos, gap_type)))
        motif2matches[mname] = matches
        log.write(f"Motif {mname}: {len(matches)} total candidate motif-matches found.\n")
        for m in matches:
            pid, mtype, mm, start, end, window, gap = m
            desc = f"type={mtype}, mismatches={mm}, window={start}-{end}, {window}"
            if gap: desc += f", gap_type={gap[1]}, gap_pos={gap[0]}"
            log.write(f"    {pid}: {desc}\n")
    log.write("\n")

    # ---------- Step 8: Assign motifs to proteins "Best match matrix" -------------#
    # Each motif selects its best protein (lowest mismatches, prefer mismatch>gap0>gapmm, then earliest, then protein order)
    # But: Each protein can only be assigned once.
    protein_assignment = {} # protein_id -> (assigned motif, match params)
    motif_assignment = {}   # motif -> (assigned protein_id, params)
    already_used_proteins = set()

    log.write("# Phase 2: Main assignment - motifs to proteins, each protein once\n")
    # Motif order guarantees motif input tiebreak
    for mname in motifs:
        # Only consider proteins not already used
        matches = [m for m in motif2matches[mname] if m[0] not in already_used_proteins]
        if not matches:
            log.write(f"Motif {mname}: No eligible proteins remain! Forcing assignment in next phase.\n")
            continue
        # Ranking preferred: least mismatches, prefer mi>gap0>gapmm, earliest start, protein input order
        matches = sorted(matches, key=lambda m: (m[2], {"mismatch":0,"gap0":1,"gapmm":2}[m[1]], m[3], protein_ids.index(m[0])))
        best = matches[0]
        pid, match_type, mismatches, start, end, window, gapinfo = best
        protein_assignment[pid] = (mname, match_type, mismatches, start, end, window, gapinfo)
        motif_assignment[mname] = (pid, match_type, mismatches, start, end, window, gapinfo)
        already_used_proteins.add(pid)
        log.write(f"Assigning motif {mname} to protein {pid} ({match_type}, mismatches={mismatches}, region={start}-{end})\n")

    # ---------- Step 9: For any unassigned motif, forcibly pick best protein --- ###
    forced_assignments = {}
    for mname in motifs:
        if mname in motif_assignment:
            continue
        # Try for unused proteins first
        candidates = [pid for pid in protein_ids if pid not in already_used_proteins]
        matches = [m for m in motif2matches[mname] if m[0] in candidates]
        if not matches:
            # All proteins assigned, use best-overall, even if results in multiple motifs on one protein
            matches = motif2matches[mname]
        if not matches:
            log.write(f"Motif {mname}: Could NOT be matched to anything (should not happen!)\n")
            continue
        best = sorted(matches, key=lambda m: (m[2], {"mismatch":0,"gap0":1,"gapmm":2}[m[1]], m[3], protein_ids.index(m[0])))[0]
        pid, match_type, mismatches, start, end, window, gapinfo = best
        forced_assignments[mname] = (pid, match_type, mismatches, start, end, window, gapinfo)
        already_used_proteins.add(pid)
        # Assign, but mark forced
        protein_assignment[pid] = (mname, match_type, mismatches, start, end, window, gapinfo, "forced")
        motif_assignment[mname] = (pid, match_type, mismatches, start, end, window, gapinfo, "forced")
        log.write(f"FORCED assign motif {mname} to protein {pid} ({match_type}, mismatches={mismatches}, {start}-{end})\n")

    # -------------- Step 10: Write output 1x per input protein, with renamed header ---------#
    log.write("# Phase 3: Writing renamed FASTA and GTF (one sequence per protein)\n")
    for pid in protein_ids:
        entry = protein_assignment.get(pid)
        if not entry:
            # No motif found or forced assignment. Can drop, or output as original.
            h = protein_headers[pid]
            seq = protein_seqs[pid]
            print(h, file=fout)
            for i in range(0, len(seq), 60):
                print(seq[i:i+60], file=fout)
            log.write(f"{pid}: No motif assigned; copied as-is.\n")
            continue
        # Unpack
        if len(entry) == 8:
            mname, match_type, mismatches, start, end, window, gapinfo, forced = entry
        else:
            mname, match_type, mismatches, start, end, window, gapinfo = entry
            forced = None
        # Rename
        if mismatches > 2:
            motif_label = mname + "-like"
        else:
            motif_label = mname
        h = protein_headers[pid]
        orig_tail = h[1+len(pid):]
        newhead = f">{pid}_{motif_label}{' [forced]' if forced else ''}{orig_tail}"
        seq = protein_seqs[pid]
        print(newhead, file=fout)
        for i in range(0, len(seq), 60):
            print(seq[i:i+60], file=fout)
        # GTF
        attrs = (f'gene_id "{pid}"; motif "{motif_label}"; matched "{window}"; motif_seq "{motif_seqs[mname]}"; mismatches "{mismatches}";')
        if gapinfo:
            attrs += f' gap_type "{gapinfo[1]}"; gap_pos "{gapinfo[0]}";'
        if forced: attrs += ' forced "true";'
        gtfrow = [pid, "MotifFinder", motif_label, str(start), str(end), ".", "+", ".", attrs]
        print('\t'.join(gtfrow), file=fgtf)
        # Log
        log.write(f"{pid}: {motif_label} assigned. Motif={mname}, Type={match_type}, mismatches={mismatches}, region={start}-{end}\n")
        log.write(f"       Window: {window}\n")
        if gapinfo:
            log.write(f"       Gap: type={gapinfo[1]}, pos={gapinfo[0]}\n")
        log.write('\n')

    # -------- Step 11: Log motif summary and ensure all motifs represented ----------
    log.write("# Summary: Assignment table motif->protein\n")
    for mname in motifs:
        assign = motif_assignment.get(mname)
        if assign:
            pid = assign[0]
            mismatches = assign[2]
            forced = " [FORCED]" if (len(assign) == 8 and assign[7] == "forced") else ""
            label = mname+"-like" if mismatches>2 else mname
            log.write(f"{mname} -> {pid} ({label}){forced}\n")
        else:
            log.write(f"{mname} -> NOT ASSIGNED\n")
    log.write("# All motifs and proteins processed.\nDONE.\n")

    fout.close(); fgtf.close(); log.close()
    print("All outputs written. See RCL.log for match details and summary.")

if __name__ == "__main__":
    main()