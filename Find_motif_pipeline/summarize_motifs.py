#!/usr/bin/env python3
"""
===============================================================================
Script: summarize_motifs.py

Purpose:
    Comprehensive, weighted summary of motifs discovered in all context windows.
    For each motif, this script:
      - Extracts motif statistics (width, info content, enrichment, e-value, etc)
      - Calculates and outputs a biologically meaningful consensus sequence
      - Outputs all instance sites as a per-motif FASTA (for logos or instance analysis)
      - Optionally creates motif logos (PNG) with WebLogo
      - Summarizes cluster, sequence and instance associations

Biological rationale:
    Key motifs have specificity for one/several clusters, high occurrence, and high info content.
    Sequence logos and consensus sequence are essential for rapid downstream analysis and scanning.

Output:
    motifs_summary.tsv: per-motif statistics and consensus/instance mapping
    motifs.fasta: consensus per motif (for scanning/tools)
    motifs_logo/: motif logos (PNG, if WebLogo is installed)
    CONTEXT_MOTIFNAME_context.fasta: instance matches (per motif)
    motifs_network.gml: motif-motif overlap/graph for visualization
===============================================================================
"""

import os, glob, argparse, sys
import pandas as pd
import numpy as np
from collections import Counter
import xml.etree.ElementTree as ET
from Bio import SeqIO

def shannon_entropy(seq):
    """Compute Shannon entropy of cluster assignments or site distribution."""
    freqs = Counter(seq)
    total = sum(freqs.values())
    return -sum((c/total)*np.log2(c/total) for c in freqs.values() if c>0) if total > 0 else 0

def parse_all_motifs(meme_context_dir, clusters_df):
    """
    Parses all context motif XMLs:
      - Extracts all stats for each motif.
      - Maps motif site `sequence_id` to sequence name using <sequence> entries.
      - Computes a consensus sequence from PWM/PPM.
      - Writes per-motif instance FASTA (real instance sequence if possible, otherwise consensus).
    """
    all_entries = []
    all_fasta = []
    xmls = sorted(glob.glob(os.path.join(meme_context_dir,"*/meme.xml")))

    # For entropy calculations, index clusters by their ID column (SeqID or SampleID)
    if not clusters_df.empty:
        col = "SeqID" if "SeqID" in clusters_df.columns else ("SampleID" if "SampleID" in clusters_df.columns else clusters_df.columns[0])
        clusters = clusters_df.set_index(col)
    else:
        clusters = None

    motif_idx = 1
    for xml in xmls:
        context = os.path.basename(os.path.dirname(xml))
        try:
            tree = ET.parse(xml)
            root = tree.getroot()

            # Mapping sequence_id -> sequence_name for context window
            seqid_to_name = {seq.attrib['id']: seq.attrib['name'] for seq in root.findall('.//sequence')}

            # Guess original input context FASTA to fetch original substrings
            for trydir in ["../tmp/contexts", "../../tmp/contexts"]:
                fasta_guess = os.path.normpath(os.path.join(meme_context_dir, trydir, f"{context}.fasta"))
                if os.path.exists(fasta_guess):
                    break
            else:
                fasta_guess = None
            seq_records = {}
            if fasta_guess and os.path.exists(fasta_guess):
                seq_records = {rec.id: str(rec.seq) for rec in SeqIO.parse(fasta_guess, "fasta")}
            else:
                print(f"[WARNING] Could not find FASTA for context {context}; will use consensus as fallback for instance FASTA.", file=sys.stderr)

            for motif_elem in root.findall('.//motif'):
                motif_id = motif_elem.attrib.get('id', f'M{motif_idx}')
                motif_name = motif_elem.attrib.get('name', motif_id)
                width = int(motif_elem.attrib.get('width', "0"))
                evalue = motif_elem.attrib.get('e_value') or motif_elem.attrib.get('evalue') or None
                try: evalue = float(evalue)
                except Exception: evalue = np.nan

                # Compute info content (sum for all columns)
                info_content = None
                pwmrows = []
                probs = motif_elem.find('probabilities')
                aa_order = []
                if probs is not None:
                    for arridx, arr in enumerate(probs.findall('./alphabet_matrix/alphabet_array')):
                        values = [float(val.text) for val in arr.findall('./value')]
                        labels = [val.attrib['letter_id'] for val in arr.findall('./value')]
                        if arridx == 0: aa_order = labels
                        pwmrows.append(values)
                    if pwmrows:
                        arr = np.array(pwmrows).T
                        with np.errstate(divide='ignore'):
                            entropy = -np.nansum(arr*np.log2(arr+1e-10),axis=1)
                        info_content = float(np.sum(np.log2(arr.shape[1])-entropy))

                # Extract all site sequences (for per-motif FASTA): prefer <contributing_site>, fallback <instance>
                seqs = []
                cluster_list = []
                for_site_context = {}
                contrib_sites = motif_elem.find('contributing_sites')
                if contrib_sites is not None:
                    for cs in contrib_sites.findall('contributing_site'):
                        seqid = cs.attrib.get('sequence_id')
                        seqname = seqid_to_name.get(seqid, seqid)
                        pos = int(float(cs.attrib.get('position', 0)))
                        instance_seq = ""
                        site_elem = cs.find('./site')
                        if site_elem is not None and list(site_elem):
                            instance_seq = ''.join([val.attrib.get('letter_id','X') for val in site_elem.findall('./letter_ref')])
                        elif seqname in seq_records:
                            instance_seq = seq_records[seqname][pos:pos+width]
                        else:
                            instance_seq = "X"*width
                        seqs.append(seqname)
                        for_site_context[seqname] = instance_seq
                        if clusters is not None and seqname in clusters.index:
                            cluster_list.append(clusters.at[seqname,"GroupID"])
                # Legacy/compat: <instance>
                for inst in motif_elem.findall('.//instance'):
                    seqid = inst.attrib.get('sequence_name') or inst.attrib.get('seq_id') or inst.attrib.get('seq')
                    pos = int(inst.attrib.get('position',0)) if 'position' in inst.attrib else 0
                    seqname = seqid_to_name.get(seqid, seqid)
                    if seqname not in seqs:
                        seqs.append(seqname)
                        if clusters is not None and seqname in clusters.index:
                            cluster_list.append(clusters.at[seqname,"GroupID"])
                    if seqname in seq_records:
                        instance_seq = seq_records[seqname][pos:pos+width]
                        for_site_context[seqname] = instance_seq
                    else:
                        for_site_context[seqname] = "X"*width

                seq_entropy = shannon_entropy(seqs)
                group_entropy = shannon_entropy(cluster_list)
                count = len(seqs)

                # Compute actual consensus from PPM
                consensus = ""
                if probs is not None and pwmrows:
                    for j in range(width):
                        if len(pwmrows[j]) == len(aa_order):
                            letter_probs = list(zip(pwmrows[j], aa_order))
                            consensus += max(letter_probs)[1]
                        else:
                            consensus += "X"
                if len(consensus) != width or not consensus:
                    consensus = "N"*width

                # --- Write per-motif instance FASTA ---
                instfasta_path = os.path.join(meme_context_dir, f"{context}_{motif_id}_context.fasta")
                with open(instfasta_path, "w") as fh:
                    for seqname in seqs:
                        site_string = for_site_context.get(seqname, consensus)
                        fh.write(f">{seqname}\n{site_string}\n")

                # --- Build summary row (Sequence is column 2 after Motif_ID) ---
                all_entries.append({
                    "Motif_ID": motif_idx,
                    "Sequence": consensus,  # New: biologically useful for display/comparison
                    "Motif_Name": f"{context}_{motif_id}",
                    "Context": context,
                    "Width": width,
                    "InfoContent": info_content,
                    "Evalue": evalue,
                    "ClusterEntropy": group_entropy,
                    "SeqEntropy": seq_entropy,
                    "Instances": ";".join(str(s) for s in seqs),
                    "Occurrence": count,
                })
                all_fasta.append(f">{context}_{motif_id}\n{consensus}\n")
                motif_idx += 1
        except Exception as ex:
            print(f"[summarize_motifs.py] Skipping motif file {xml} ({ex})", file=sys.stderr)
    return pd.DataFrame(all_entries), all_fasta

def create_motif_logos(summarydf, context_dir, outdir):
    """
    For each motif, runs WebLogo on its instance FASTA file to create a PNG logo.
    Detects peptide logo by default.
    """
    os.makedirs(f"{outdir}/motifs_logo", exist_ok=True)
    for name in summarydf["Motif_Name"].unique():
        cfile = os.path.join(context_dir, f"{name}_context.fasta")
        ofile = os.path.join(outdir, "motifs_logo", f"{name}.png")
        if os.path.exists(cfile):
            fmt_flag = '-A protein'  # Change to '-A dna' for DNA motifs
            os.system(f"weblogo --resolution 600 -F png {fmt_flag} -o {ofile} < {cfile} 2>/dev/null")

def create_motif_network(summarydf, out_path):
    """
    Create a motif-motif (or motif-group) network graph for downstream visualization.
    Nodes: motifs; Edge: Jaccard (normalized overlap) of shared instance sequences.
    """
    import networkx as nx
    G = nx.Graph()
    for idx, row1 in summarydf.iterrows():
        for jdx, row2 in summarydf.iterrows():
            if idx >= jdx: continue
            set1 = set(row1["Instances"].split(";"))
            set2 = set(row2["Instances"].split(";"))
            overlap = len(set1 & set2) / max(len(set1|set2),1)
            if overlap > 0:
                G.add_edge(row1["Motif_Name"], row2["Motif_Name"], weight=overlap)
    nx.write_gml(G, out_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--meme_dir", required=True,
                        help="(Unused, for legacy compat—provide any value)")
    parser.add_argument("--context_dir", required=True,
                        help="Directory containing context motif subfolders (each with meme.xml)")
    parser.add_argument("--clusters", required=True,
                        help="CSV file with sample/sequence IDs and cluster assignments")
    parser.add_argument("--out_dir", required=True,
                        help="Directory to write summary, fasta, logos, and network")
    args = parser.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    clustersdf = pd.read_csv(args.clusters) if os.path.isfile(args.clusters) else pd.DataFrame()
    motifsumm, motiffasta = parse_all_motifs(args.context_dir, clustersdf)

    # Ensure Sequence is the second column (after Motif_ID)
    summary_cols = ["Motif_ID", "Sequence", "Motif_Name", "Context", "Width",
                    "InfoContent", "Evalue", "ClusterEntropy", "SeqEntropy", "Instances", "Occurrence"]
    motifsumm = motifsumm[summary_cols]

    motifsumm.to_csv(os.path.join(args.out_dir, "motifs_summary.tsv"), sep="\t", index=False)
    with open(os.path.join(args.out_dir, "motifs.fasta"), "w") as out:
        out.writelines(motiffasta)

    try: create_motif_network(motifsumm, os.path.join(args.out_dir, "motifs_network.gml"))
    except Exception: pass

    try: create_motif_logos(motifsumm, args.context_dir, args.out_dir)
    except Exception: pass

    print(f"[summarize_motifs.py] Motif summary, consensus FASTA, per-motif instance FASTA, logos, and network written.")