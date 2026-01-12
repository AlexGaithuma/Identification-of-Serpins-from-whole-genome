#!/usr/bin/env python3
"""
===============================================================================
Script: extract_contexts.py

Purpose: This script extracts unbiased, result-driven, and biologically/practically relevant 
context regions for motif/grouping analyses in serpin proteins (or homologs).

**Strategic Summary:**
    - All sequence regions are treated with an agnostic approach: regions are determined 
      by observed variability, motif calls, or biological boundaries—without preconceived
      structural anchor constraints. This prevents bias toward known motif/active-site frames.
    - Biological focus is on capturing both conserved (e.g., motif) and adaptive/variable 
      (e.g., hypervariable loop or insert) functional signatures, reflecting both ancient
      evolutionary pressures and recent functional innovation.

**For each input set, the output is tailored for downstream clustering and summary, and
all steps are robust to gaps/missing data if inputs are correctly prepared.**

Inputs:
    <input_fasta>:      FASTA file of all current (primary) serpin protein sequences.
    <input_rcl>:        FASTA containing all RCL (reactive center loop) region(s); should be pre-extracted or annotated.
    <input_nterm>:      FASTA containing all N-terminal region(s) of interest.
    <meme_xml>:         The XML output file from a global MEME motif discovery run.
    <out_context_dir>:  Directory where all derived context FASTA files will be written.

Outputs:
    * N-terminal regions, RCL regions, full-length sequences, top de novo "variable" 
      windows (entropy-driven, evidence of adaptive evolution/diversification), and 
      all motif-local windows (from global motif search).
===============================================================================
"""

import sys, os
import subprocess
import tempfile
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
import numpy as np
import xml.etree.ElementTree as ET

def die(msg, code=1):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(code)

def write_fasta(path, entries, overwrite=True):
    """
    Outputs a FASTA file given a list of (name, seq) tuples.
    Overwrite prevention option for pipeline safety.
    """
    if not overwrite and os.path.isfile(path):
        die(f"Refusing to overwrite existing file: {path}")
    with open(path, "w") as out:
        for name, seq in entries:
            out.write(f">{name}\n{seq}\n")

def extract_nterm_context(seqrecords, length=35):
    """
    For each input protein sequence, extract the N-terminal segment of specified length.
    - **Biological rationale:** The N-terminus can harbor localization signals, regulatory regions, 
      or even sites for co- or post-translational modification.
    - **Conserved/variable balance:** The minimum length filter ensures robustness if inputs have 
      partial N/C truncation.
    """
    contexts = []
    for rec in seqrecords:
        subseq = str(rec.seq)[:length]
        if len(subseq) >= length//2:
            contexts.append((rec.id + "_Nterm", subseq))
    return contexts

def extract_rcl_contexts_from_fasta(rcl_records, minlen=8):
    """
    Use all RCL (reactive center loop) regions as is from the provided FASTA.
    - **Biological rationale:** The RCL is the key substrate recognition region in serpins.
    - **Why filter very short?** To avoid spurious (truncated) fragments that don't represent 
      the full functional loop.
    """
    return [(rec.id + "_RCL", str(rec.seq)) for rec in rcl_records if len(str(rec.seq)) >= minlen]

def extract_full_seqs(seqrecords):
    """
    Uses full-length protein sequence for each entry after stripping alignment gaps.
    - **Biological rationale:** Full sequence contexts are important for downstream 
      "all-vs-all" motif analysis and for general-purpose clustering.
    """
    return [(rec.id + "_Full", str(rec.seq).replace("-", "")) for rec in seqrecords]

def check_homogeneous_length(seqrecords):
    """
    Checks whether all SeqRecords have the same sequence length.
    - **Required for:** Entropy matrix and sliding window-based analyses.
    """
    lengths = [len(rec.seq) for rec in seqrecords]
    if len(set(lengths)) != 1:
        return False
    return True

def run_mafft(input_fasta, thread=4, threadtb=5, threadit=0):
    """
    Runs MAFFT alignment externally and returns aligned SeqRecords.
    - **Biological rationale:** Sequence alignment is critical for meaningful entropy-based
      site-by-site variability analysis. Without an alignment, 'positions' are not comparable.
    - **Pipeline integration:** Automates MAFFT invocation so users with only unaligned
      input do not need to run it manually.
    """
    with tempfile.NamedTemporaryFile("w", delete=False, suffix=".fasta") as tmp_in, \
         tempfile.NamedTemporaryFile("r", delete=False, suffix=".fasta") as tmp_out:
        SeqIO.write(input_fasta, tmp_in, "fasta")
        tmp_in_name = tmp_in.name
        tmp_out_name = tmp_out.name

    # Command as per user specification
    cmd = [
        "mafft",
        "--thread", str(thread),
        "--threadtb", str(threadtb),
        "--threadit", str(threadit),
        "--reorder",
        "--maxiterate", "100000",
        "--retree", "1",
        "--globalpair",
        tmp_in_name
    ]
    try:
        with open(tmp_out_name, "w") as outf:
            subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, check=True)
        aligned_records = list(SeqIO.parse(tmp_out_name, "fasta"))
    except Exception as e:
        die(f"MAFFT alignment failed: {e}")
    finally:
        os.unlink(tmp_in_name)
        os.unlink(tmp_out_name)
    return aligned_records

def compute_entropy_matrix(seqrecords):
    """
    Computes Shannon entropy at each residue position across a set of equally-long sequences.
    - **Biological rationale:** Highlights positions under diversifying selection, adaptive
      function, or rapid evolution, serving as seed regions for de novo motif search.
    - **Gap filtering:** Gaps and ambiguous sites are omitted per position.
    """
    arr = np.array([list(str(x.seq)) for x in seqrecords])
    ncol = arr.shape[1]
    entropy = []
    for i in range(ncol):
        cols = [aa for aa in arr[:,i] if aa not in "-Xx.?*"]
        total = len(cols)
        if total == 0:
            entropy.append(0.)
        else:
            cnts = Counter(cols)
            ent = -sum((cnt/total) * np.log2(cnt/total) for cnt in cnts.values())
            entropy.append(ent)
    return arr, np.array(entropy)

def extract_top_variable_windows(seqrecords, winlen=20, min_entropy=1.2, top_n=2, verbose=True):
    """
    Slides a window of given length over the aligned input sequences and extracts the most
    diverse region(s) as indicated by local mean Shannon entropy.
    - **Biological rationale:** These represent adaptively evolving "exposed" or "loop" regions,
      new substrate binding sites, or the result of functional divergence.
    - **Practicality:** Redundancy avoidance: do not select overlapping regions.
    - **Robustness:** Ignores windows that fail the minimum mean-entropy cutoff.
    """
    if not check_homogeneous_length(seqrecords):
        if verbose:
            print("INFO: Aligning sequences with MAFFT for high-entropy window extraction.", file=sys.stderr)
        seqrecords = run_mafft(seqrecords)
    arr, entropy = compute_entropy_matrix(seqrecords)
    ncol = arr.shape[1]
    winstats = []
    for i in range(ncol - winlen + 1):
        mean_ent = np.mean(entropy[i:i+winlen])
        winstats.append((mean_ent,i,i+winlen))
    winstats.sort(reverse=True)
    selected = []
    used = np.zeros(ncol, dtype=bool)
    for ent, start, stop in winstats:
        if ent < min_entropy: continue
        if np.any(used[start:stop]): continue
        selected.append((start,stop))
        used[start:stop] = True
        if len(selected) == top_n: break
    out = []
    for widx, (start, stop) in enumerate(selected, 1):
        for rec, arrseq in zip(seqrecords, arr):
            subseq = "".join([c for c in arrseq[start:stop] if c not in "-Xx.?*"])
            if len(subseq) >= (stop-start)//2:   # At least half non-gap
                out.append((f"{rec.id}_VarWin{widx}", subseq))
    return out

def extract_all_meme_motif_windows(seqrecords, meme_xml, context_len=21, verbose=True):
    """
    Extracts a local context window for each instance of each MEME motif in the global motif XML file.
    - **Biological rationale:** Provides both motif-core and its functional surroundings,
      crucial for context-aware clustering or evolutionary footprinting.
    - **Robust prefix matching:** Handles differences in sequence naming (removes pipeline header artifacts).
    """
    motif_windows = defaultdict(list)
    try:
        tree = ET.parse(meme_xml)
    except Exception as e:
        die(f"Could not parse MEME XML at '{meme_xml}': {e}")

    root = tree.getroot()
    motifs_elements = [m for m in root.findall('.//motif')]
    for mi, motif_elem in enumerate(motifs_elements):
        motif_id = motif_elem.get('id', f'Motif_{mi+1}')
        motif_name = motif_elem.get('name', motif_id)
        motif_width = int(motif_elem.get('width', "0"))
        motifkey = f"Motif{motif_id}_{motif_name}".replace(" ", "_")
        insts = []
        for instance in motif_elem.findall('.//instance'):
            seqid = instance.get('sequence_name') or instance.get('seq_id') or instance.get('seq') or None
            pos = int(float(instance.get('position', instance.get('pos', "0"))))
            width = int(instance.get('width', motif_width))
            insts.append((seqid, pos, width))
        for seqid, pos, width in insts:
            matched_seq, recid = None, None
            for rec in seqrecords:
                # If the ID matches exactly or prefix (handle _sep cases)
                if rec.id == seqid or rec.id.split('_')[0] == str(seqid).split('_')[0]:
                    matched_seq = str(rec.seq)
                    recid = rec.id
                    break
            if matched_seq and recid:
                start = max(0, pos - (context_len//2))
                end = min(len(matched_seq), pos + width + (context_len//2))
                region = matched_seq[start:end]
                if len(region) >= width:
                    motif_windows[motifkey].append( (f"{recid}_{motifkey}_at_{pos}", region) )
    return motif_windows

def parse_fasta_file(path):
    """
    Fast and safe parse of a FASTA file.
    """
    if not os.path.isfile(path):
        die(f"Input file not found: {path}")
    try:
        return list(SeqIO.parse(path, 'fasta'))
    except Exception as e:
        die(f"Failed to parse FASTA file '{path}': {e}")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print(__doc__, file=sys.stderr)
        sys.exit(1)
    in_fasta, in_rcl, in_nterm, meme_xml, outdir = sys.argv[1:6]
    os.makedirs(outdir, exist_ok=True)

    seqrecords = parse_fasta_file(in_fasta)
    rcl_records = parse_fasta_file(in_rcl)
    nterm_records = parse_fasta_file(in_nterm)

    # 1. All N-terminal region contexts (each is a region from the user-supplied N-term FASTA)
    nterm_contexts = extract_rcl_contexts_from_fasta(nterm_records, minlen=4)
    write_fasta(f"{outdir}/Nterm_context.fasta", nterm_contexts)

    # 2. All RCL motif/context windows (biologically-critical for serpin substrate specificity)
    rcl_contexts = extract_rcl_contexts_from_fasta(rcl_records, minlen=8)
    write_fasta(f"{outdir}/RCL_context.fasta", rcl_contexts)

    # 3. Full sequence context (stripped of gaps; provides background for global clustering)
    full_contexts = extract_full_seqs(seqrecords)
    write_fasta(f"{outdir}/Full_context.fasta", full_contexts)

    # 4. De novo/high-entropy "variable" windows—useful for investigating adaptive sites/clusters.
    variable_contexts = extract_top_variable_windows(seqrecords, winlen=21, min_entropy=1.2, top_n=2)
    if variable_contexts:
        write_fasta(f"{outdir}/VarWin_context.fasta", variable_contexts)

    # 5. All MEME/global motif windows (core motif and context, for every instance/site/sequence)
    meme_contexts = extract_all_meme_motif_windows(seqrecords, meme_xml, context_len=21)
    for motifname, instances in meme_contexts.items():
        write_fasta(f"{outdir}/{motifname}_context.fasta", instances)