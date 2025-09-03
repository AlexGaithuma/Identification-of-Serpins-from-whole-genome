#!/usr/bin/env python3

"""
Advanced RCL Serpin Clustering Pipeline with ESM embeddings + motif features + supervised grouping.

This script performs a sophisticated clustering that:

1. Loads all RCL serpin sequences and cleans them.
2. Deduplicates the sequences, embeds them via ESM-1b transformer embedding.
3. Identifies sequence motifs using the MEME suite for motif discovery.
4. Extracts motif presence vectors per sequence (binary or scores).
5. Combines ESM embeddings + motif presence vectors into joint feature vectors.
6. Performs semi-supervised clustering:
   - Uses sequences in `--in_group1` and `--in_group2` as known clusters.
   - Leverages these sets to guide clustering by training a semi-supervised classifier
     (e.g. using embeddings+motif features to build a feature-space).
   - Applies advanced unsupervised clustering (HDBSCAN) on all sequences afterward.
7. Labels previously classified sequences accordingly and assigns others based on clustering.
8. Produces original outputs:
   - `rcl_dendrogram.png`
   - `rcl_tsne.png`
   - `group_sizes.png`
9. Outputs a cluster assignment file including every sequence (duplicates included) mapped to clusters.
10. Optionally, performs Multiple Sequence Alignments (MSA) per cluster and produces sequence logos.

Dependencies (install via pip):

    pip install biopython torch fair-esm hdbscan umap-learn logomaker matplotlib seaborn pandas
    # Install MEME Suite separately from https://meme-suite.org/meme/ and add meme-tools to PATH.
    # MAFFT must be installed and accessible from PATH.

Usage example:

    python cluster_rcl_serpins.py --in_RCLs all_sequences.fasta --in_group1 known_group1.fasta --in_group2 known_group2.fasta \
        -out_clusters clusters_output.txt --fig_dir figures --tmp_dir tmp_dir

NOTES:

- MEME motif discovery is computationally expensive. Consider limiting input size or pre-run motifs when scaling up.
- This approach synergizes sequence deep embeddings (ESM) reflecting global sequence context,
  and motif features capturing conserved local patterns critical for function (e.g. RCL specificity).
- Semi-supervised classification helps leverage biologically validated groups to guide partitioning,
  improving biologically meaningful clusters.
- Advanced alternatives: incorporate structure prediction embeddings, use graph-based protein embeddings,
  or pretrained motif representations.
- Future work could integrate phylogenetic methods or combine sequence+structure clustering for evolutionary insight.

"""

import os
import sys
import argparse
import subprocess
import tempfile
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import hdbscan
import umap.umap_ as umap

from Bio import SeqIO, AlignIO
from Bio.Seq import Seq

import torch

try:
    import esm
except ImportError:
    raise ImportError("Please install 'fair-esm' via 'pip install fair-esm' to run this script.")

try:
    import logomaker
except ImportError:
    raise ImportError("Please install 'logomaker' via 'pip install logomaker' for logo generation.")

from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE


def clean_aa_sequence(seq):
    """
    Clean protein sequence by uppercasing and replacing non-standard AAs with 'X'.
    Essential for consistent embedding input.

    Args:
        seq (str): Raw AA sequence.
    Returns:
        str: Cleaned AA sequence.
    """
    standard_aas = set("ACDEFGHIKLMNPQRSTVWY")
    seq = seq.upper()
    cleaned = "".join(aa if aa in standard_aas else "X" for aa in seq)
    return cleaned


def load_fasta_sequences(fasta_file):
    """
    Load sequences from a fasta file and clean them.

    Args:
        fasta_file (str): Path to fasta file.
    Returns:
        list of SeqRecord: cleaned sequences.
    """
    records = list(SeqIO.parse(fasta_file, "fasta"))
    for rec in records:
        rec.seq = Seq(clean_aa_sequence(str(rec.seq)))
    return records


def deduplicate_sequences(records):
    """
    Remove exact duplicates (cleaned sequences), keeping first occurrence.

    Args:
        records (list of SeqRecord)
    Returns:
        list of unique SeqRecord
    """
    seen = {}
    unique_records = []
    for rec in records:
        seq_str = str(rec.seq)
        if seq_str not in seen:
            seen[seq_str] = rec.id
            unique_records.append(rec)
    return unique_records


def embed_sequences_esm(records, batch_size=16, device=None):
    """
    Embed sequences using ESM1b model averages token embeddings for representation.

    Args:
        records (list of SeqRecord): sequences to embed.
        batch_size (int)
        device (str or None): 'cuda'/'cpu', autodetected if None.

    Returns:
        np.ndarray: embeddings array (num_seqs x embedding_dim)
    """
    import esm

    if device is None:
        device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"[INFO] Using device '{device}' for ESM embedding.")

    model, alphabet = esm.pretrained.esm1b_t33_650M_UR50S()
    model.eval()
    model = model.to(device)
    batch_converter = alphabet.get_batch_converter()

    data = [(rec.id, str(rec.seq)) for rec in records]

    all_embeddings = []

    with torch.no_grad():
        for i in range(0, len(data), batch_size):
            batch = data[i : i + batch_size]
            labels, strs, tokens = batch_converter(batch)
            tokens = tokens.to(device)
            results = model(tokens, repr_layers=[33], return_contacts=False)
            token_embeddings = results["representations"][33]
            sequence_representations = token_embeddings[:, 1:-1].mean(1)
            all_embeddings.append(sequence_representations.cpu().numpy())

    embeddings = np.vstack(all_embeddings)
    print(f"[INFO] Computed embeddings for {embeddings.shape[0]} sequences with dim {embeddings.shape[1]}.")
    return embeddings


def run_meme(input_fasta, output_dir, nmotifs=5, minw=6, maxw=12):
    """
    Run MEME motif discovery on input fasta sequences.

    Args:
        input_fasta (str): Input fasta path
        output_dir (str): Output directory
        nmotifs (int): number of motifs to find
        minw (int): minimum motif width
        maxw (int): maximum motif width

    Returns:
        str: path to MEME motif XML output
    """
    os.makedirs(output_dir, exist_ok=True)
    meme_out_dir = os.path.join(output_dir, "meme_out")
    if not os.path.exists(meme_out_dir):
        os.makedirs(meme_out_dir)

    # FIXED MEME CALL:
    # 1. Remove -dna
    # 2. Keep only -protein
    # 3. Replace -quiet by -nostatus (to suppress progress) or omit for verbose
    # 4. Add -maxsize with a large number to avoid size limit stopping
    # 5. Use correct flags for output directory (-oc) or (-o)

    cmd = [
        "meme",
        input_fasta,
        "-oc",
        meme_out_dir,
        "-protein",
        "-mod",
        "zoops",
        "-nmotifs",
        str(nmotifs),
        "-minw",
        str(minw),
        "-maxw",
        str(maxw),
        "-nostatus",
        "-maxsize",
        "1000000",  # 1 million characters max size to handle larger inputs
    ]

    print(f"[INFO] Running MEME motif discovery on {input_fasta}...")

    try:
        subprocess.run(cmd, check=True)
        print(f"[INFO] MEME completed; output dir: {meme_out_dir}")
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] MEME failed: {e}")
        sys.exit(1)

    meme_xml = os.path.join(meme_out_dir, "meme.xml")
    if not os.path.exists(meme_xml):
        print(f"[ERROR] MEME XML output not found at {meme_xml}")
        sys.exit(1)

    return meme_xml


def parse_meme_motif_presence(meme_xml, sequences):
    """
    Parse MEME XML output and determine per sequence motif presence matrix.

    Uses MAST or direct parsing.

    Args:
        meme_xml (str): Path to MEME XML output with motifs.
        sequences (list of SeqRecord): sequences in same order as input to MEME.

    Returns:
        np.ndarray of shape (num_sequences, num_motifs): binary or float scores indicating motif presence.
    """
    import xml.etree.ElementTree as ET
    import tempfile

    meme_out_dir = os.path.dirname(meme_xml)
    motifs_file = os.path.join(meme_out_dir, "meme.txt")
    if not os.path.exists(motifs_file):
        print(f"[ERROR] MEME motifs text output not found (expected {motifs_file})")
        sys.exit(1)

    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp_seq_fa:
        SeqIO.write(sequences, tmp_seq_fa.name, "fasta")
        query_fasta = tmp_seq_fa.name

    mast_out = os.path.join(meme_out_dir, "mast_out")
    if not os.path.exists(mast_out):
        os.makedirs(mast_out)

    mast_cmd = [
        "mast",
        motifs_file,
        query_fasta,
        "-oc",
        mast_out,
        "-nohtml",
        "-nostatus",
        "-sep",
    ]

    print("[INFO] Running MAST to identify motif hits per sequence...")
    try:
        subprocess.run(mast_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] MAST failed: {e}")
        sys.exit(1)

    mast_xml = os.path.join(mast_out, "mast.xml")
    if not os.path.exists(mast_xml):
        print(f"[ERROR] Expected MAST XML output missing at {mast_xml}")
        sys.exit(1)

    tree = ET.parse(mast_xml)
    root = tree.getroot()

    meme_tree = ET.parse(meme_xml)
    meme_root = meme_tree.getroot()
    motif_ids = [motif.attrib['id'] for motif in meme_root.findall("motifs/motif")]
    n_motifs = len(motif_ids)
    s_to_idx = {seq.id: idx for idx, seq in enumerate(sequences)}
    n_seqs = len(sequences)

    motif_matrix = np.zeros((n_seqs, n_motifs), dtype=float)

    for seq_elem in root.findall("sequence"):
        seq_id = seq_elem.attrib.get("label")
        if seq_id not in s_to_idx:
            continue
        seq_idx = s_to_idx[seq_id]
        for hit in seq_elem.findall("hit"):
            motif_id = hit.attrib.get("motif")
            if motif_id not in motif_ids:
                continue
            motif_idx = motif_ids.index(motif_id)
            motif_matrix[seq_idx, motif_idx] = 1.0

    os.remove(query_fasta)
    return motif_matrix


def combine_features(esm_embeddings, motif_features):
    scaler = StandardScaler()
    esm_scaled = scaler.fit_transform(esm_embeddings)
    motif_scaled = StandardScaler().fit_transform(motif_features)
    combined = np.hstack([esm_scaled, motif_scaled])
    return combined


def semi_supervised_group_classifier(
    combined_features,
    unique_seq_records,
    known_group1_records,
    known_group2_records,
):
    seq_to_idx = {str(rec.seq): i for i, rec in enumerate(unique_seq_records)}
    labels = np.full(shape=(len(unique_seq_records),), fill_value=-2, dtype=int)
    for rec in known_group1_records:
        idx = seq_to_idx.get(str(rec.seq))
        if idx is not None:
            labels[idx] = 0
    for rec in known_group2_records:
        idx = seq_to_idx.get(str(rec.seq))
        if idx is not None:
            labels[idx] = 1

    labeled_mask = labels >= 0
    X_train = combined_features[labeled_mask]
    y_train = labels[labeled_mask]

    X_unlabeled = combined_features[~labeled_mask]

    if len(np.unique(y_train)) < 2 or X_train.shape[0] < 10:
        print("[WARNING] Insufficient labeled sequences to train supervised classifier; skipping classification.")
        return None

    print(f"[INFO] Training supervised classifier on {X_train.shape[0]} sequences in known groups.")

    clf = LogisticRegression(max_iter=1000, class_weight="balanced", random_state=42)
    clf.fit(X_train, y_train)

    pred_labels = clf.predict(X_unlabeled)

    labels[~labeled_mask] = pred_labels

    print("[INFO] Semi-supervised classification complete.")
    return labels


def cluster_unknowns_with_hdbscan(embeddings, preassigned_labels, min_cluster_size=5, min_samples=5):
    unknown_mask = preassigned_labels == -2
    unknown_embeddings = embeddings[unknown_mask]

    if unknown_embeddings.shape[0] == 0:
        print("[INFO] No unknown sequences to cluster.")
        return preassigned_labels

    print(f"[INFO] Clustering {unknown_embeddings.shape[0]} unknown sequences with HDBSCAN...")

    clusterer = hdbscan.HDBSCAN(
        min_cluster_size=min_cluster_size,
        min_samples=min_samples,
        metric="euclidean",
        cluster_selection_method="eom",
    )
    cluster_labels = clusterer.fit_predict(unknown_embeddings)

    print(f"[INFO] HDBSCAN found {len(set(cluster_labels)) - (1 if -1 in cluster_labels else 0)} clusters among unknowns.")

    mapped_labels = np.full_like(preassigned_labels, fill_value=-1)
    mapped_labels[preassigned_labels == 0] = 0
    mapped_labels[preassigned_labels == 1] = 1

    max_known_label = 1
    unique_clusters = set(cluster_labels)
    if -1 in unique_clusters:
        unique_clusters.remove(-1)

    cluster_id_map = {}
    next_label = max_known_label + 1
    for cl in sorted(unique_clusters):
        cluster_id_map[cl] = next_label
        next_label += 1

    unk_indices = np.where(unknown_mask)[0]
    for idx_input, cl in zip(unk_indices, cluster_labels):
        if cl == -1:
            mapped_labels[idx_input] = -1
        else:
            mapped_labels[idx_input] = cluster_id_map[cl]

    return mapped_labels


def map_clusters_to_all_sequences(all_records, unique_records, unique_labels):
    seq_to_label = {str(rec.seq): label for rec, label in zip(unique_records, unique_labels)}
    all_assignments = []
    for rec in all_records:
        cleaned_seq = clean_aa_sequence(str(rec.seq))
        label = seq_to_label.get(cleaned_seq, -1)
        all_assignments.append((rec, label))
    return all_assignments


def plot_dendrogram(embeddings, labels, seq_ids, outdir):
    import scipy.cluster.hierarchy as sch
    print("[INFO] Generating dendrogram...")
    linkage = sch.linkage(embeddings, method="ward")
    plt.figure(figsize=(14, 6))
    dendro = sch.dendrogram(linkage, labels=seq_ids, leaf_rotation=90)
    ax = plt.gca()
    xlbls = ax.get_xmajorticklabels()

    label_dict = dict(zip(seq_ids, labels))
    max_label = max(labels) if len(labels) else 0
    palette = sns.color_palette("tab10", n_colors=max_label + 2)

    for lbl in xlbls:
        seq_id = lbl.get_text()
        cluster_id = label_dict.get(seq_id, -1)
        if cluster_id == -1:
            lbl.set_color("black")
        else:
            lbl.set_color(palette[cluster_id % len(palette)])

    plt.title("Hierarchical Clustering Dendrogram of RCL sequences")
    plt.tight_layout()
    fn = os.path.join(outdir, "rcl_dendrogram.png")
    plt.savefig(fn, dpi=300)
    plt.close()
    print(f"[INFO] Saved dendrogram to {fn}")


def plot_tsne(embeddings, labels, outdir):
    print("[INFO] Running t-SNE...")
    tsne_model = TSNE(
        n_components=2,
        perplexity=15,
        metric="euclidean",
        random_state=42,
        n_jobs=-1,
        init="pca",
        learning_rate="auto",
    )
    tsne_results = tsne_model.fit_transform(embeddings)

    plt.figure(figsize=(8, 6))
    unique_labels = sorted(set(labels))
    if -1 in unique_labels:
        unique_labels = [l for l in unique_labels if l != -1] + [-1]

    palette = sns.color_palette("tab10", n_colors=len(unique_labels))

    for idx, clust_id in enumerate(unique_labels):
        idxs = np.where(labels == clust_id)[0]
        label_str = "Noise" if clust_id == -1 else f"Group {clust_id + 1}"
        plt.scatter(
            tsne_results[idxs, 0],
            tsne_results[idxs, 1],
            label=label_str,
            s=50,
            edgecolors="k",
            linewidth=0.4,
            alpha=0.7,
            color=palette[idx],
        )

    plt.legend()
    plt.title("t-SNE of RCL Protein Embeddings Colored by Cluster")
    plt.xlabel("Dimension 1")
    plt.ylabel("Dimension 2")
    plt.tight_layout()

    fn = os.path.join(outdir, "rcl_tsne.png")
    plt.savefig(fn, dpi=300)
    plt.close()
    print(f"[INFO] Saved t-SNE plot to {fn}")


def plot_cluster_sizes(labels, outdir):
    from collections import Counter
    counts = Counter(labels)
    sorted_clusters = sorted(counts.keys())
    sizes = [counts[k] for k in sorted_clusters]
    display_labels = ["Noise" if k == -1 else f"Group {k + 1}" for k in sorted_clusters]

    plt.figure(figsize=(8, 5))
    palette = sns.color_palette("tab10", n_colors=len(sorted_clusters))
    plt.bar(display_labels, sizes, color=palette)
    plt.xticks(rotation=45)
    plt.ylabel("Number of sequences")
    plt.title("Sequences per Cluster / Group")
    plt.tight_layout()
    fn = os.path.join(outdir, "group_sizes.png")
    plt.savefig(fn, dpi=300)
    plt.close()
    print(f"[INFO] Saved cluster sizes plot to {fn}")


def run_mafft(records, tmp_dir, out_alignment):
    input_fp = os.path.join(tmp_dir, "to_align.fasta")
    SeqIO.write(records, input_fp, "fasta")

    cmd = ["mafft", "--auto", "--quiet", input_fp]

    try:
        with open(out_alignment, "w") as out_f:
            subprocess.run(cmd, check=True, stdout=out_f, stderr=subprocess.DEVNULL)
    except subprocess.CalledProcessError:
        print(f"[ERROR] MAFFT failed on {input_fp}")
        return None

    alignment = AlignIO.read(out_alignment, "fasta")
    return alignment


def alignment_to_logomatrix(alignment):
    aas = list("ACDEFGHIKLMNPQRSTVWY-")
    length = alignment.get_alignment_length()

    rows = []
    for pos in range(length):
        counts = dict.fromkeys(aas, 0)
        for rec in alignment:
            aa = rec.seq[pos].upper()
            if aa not in aas:
                aa = "-"
            counts[aa] += 1
        total = sum(counts.values())
        freqs = {aa: counts[aa] / total if total else 0 for aa in aas}
        rows.append(freqs)

    df = pd.DataFrame(rows)
    return df


def plot_sequence_logo(alignment, group_label, fig_dir):
    freq_matrix = alignment_to_logomatrix(alignment)
    width = max(10, freq_matrix.shape[0] * 0.5)

    plt.figure(figsize=(width, 4))
    logo = logomaker.Logo(
        freq_matrix,
        shade_below=0.5,
        fade_below=0.5,
        font_name="Arial Rounded MT Bold",
    )
    plt.title(f"Sequence Logo - {group_label}")
    plt.ylabel("Frequency")
    plt.xlabel("Position")
    plt.tight_layout()

    out_fp = os.path.join(fig_dir, f"seq_logo_{group_label.replace(' ', '_')}.png")
    plt.savefig(out_fp, dpi=300)
    plt.close()
    print(f"[INFO] Saved sequence logo for {group_label} to {out_fp}")


def main():
    parser = argparse.ArgumentParser(description="Cluster RCL serpins with embeddings, motifs, and supervised guidance.")

    parser.add_argument("--in_RCLs", required=True, help="Input FASTA with RCL sequences")
    parser.add_argument("--in_group1", required=True, help="FASTA of biologically confirmed Group 1 sequences")
    parser.add_argument("--in_group2", required=True, help="FASTA of biologically confirmed Group 2 sequences")
    parser.add_argument("-out_clusters", required=True, help="Output clusters assignment file")
    parser.add_argument("--fig_dir", default="figures", help="Directory to save figures")
    parser.add_argument("--tmp_dir", default="tmp", help="Temporary directory for intermediate files")
    parser.add_argument("--min_cluster_size", type=int, default=5, help="HDBSCAN minimum cluster size")
    parser.add_argument("--min_samples", type=int, default=5, help="HDBSCAN min samples")
    parser.add_argument("--batch_size", type=int, default=16, help="Batch size for ESM embedding")

    args = parser.parse_args()

    os.makedirs(args.fig_dir, exist_ok=True)
    os.makedirs(args.tmp_dir, exist_ok=True)

    print(f"[INFO] Loading all sequences from {args.in_RCLs} ...")
    all_records = load_fasta_sequences(args.in_RCLs)
    print(f"[INFO] Loaded {len(all_records)} sequences.")

    print(f"[INFO] Loading known biologically confirmed groups.")
    group1_recs = load_fasta_sequences(args.in_group1)
    group2_recs = load_fasta_sequences(args.in_group2)

    print(f"[INFO] Deduplicating sequences for embedding...")
    unique_records = deduplicate_sequences(all_records)
    print(f"[INFO] {len(unique_records)} unique sequences found.")

    print("[INFO] Generating ESM embeddings for unique sequences...")
    esm_embeddings = embed_sequences_esm(unique_records, batch_size=args.batch_size)

    unique_fasta_path = os.path.join(args.tmp_dir, "unique_sequences.fasta")
    SeqIO.write(unique_records, unique_fasta_path, "fasta")

    # run the fixed MEME command here
    meme_xml = run_meme(unique_fasta_path, args.tmp_dir, nmotifs=5, minw=6, maxw=12)

    motif_features = parse_meme_motif_presence(meme_xml, unique_records)
    print(f"[INFO] Obtained motif presence matrix with shape {motif_features.shape}.")

    combined_features = combine_features(esm_embeddings, motif_features)
    print(f"[INFO] Combined feature matrix shape: {combined_features.shape}")

    known_g1 = deduplicate_sequences(group1_recs)
    known_g2 = deduplicate_sequences(group2_recs)

    labels_semisupervised = semi_supervised_group_classifier(
        combined_features,
        unique_records,
        known_g1,
        known_g2,
    )

    if labels_semisupervised is None:
        labels_semisupervised = np.full((len(unique_records),), -2, dtype=int)

    final_labels = cluster_unknowns_with_hdbscan(
        combined_features,
        labels_semisupervised,
        min_cluster_size=args.min_cluster_size,
        min_samples=args.min_samples,
    )

    all_assignments = map_clusters_to_all_sequences(all_records, unique_records, final_labels)

    with open(args.out_clusters, "w") as out_f:
        cluster_groups = defaultdict(list)
        for rec, label in all_assignments:
            cluster_groups[label].append(rec)

        for clust_label in sorted(cluster_groups.keys()):
            group_name = "Noise" if clust_label == -1 else f"Group_{clust_label + 1}"
            out_f.write(f">{group_name}\n")
            for rec in cluster_groups[clust_label]:
                out_f.write(f"{rec.id}\t{str(rec.seq)}\n")
            out_f.write("\n")

    print(f"[INFO] Saved cluster assignments to {args.out_clusters}.")

    unique_ids = [rec.id for rec in unique_records]

    plot_dendrogram(combined_features, final_labels, unique_ids, args.fig_dir)
    plot_tsne(combined_features, final_labels, args.fig_dir)
    plot_cluster_sizes(final_labels, args.fig_dir)

    print("[INFO] Generating MSAs and sequence logos per cluster for unique sequences...")

    for clust in sorted(set(final_labels)):
        if clust == -1:
            print("[INFO] Skipping noise group for logos.")
            continue
        group_name = f"Group_{clust + 1}"
        cluster_recs = [rec for rec, lbl in zip(unique_records, final_labels) if lbl == clust]
        if len(cluster_recs) < 2:
            print(f"[INFO] Skipping {group_name}, too few sequences (<2) for alignment.")
            continue
        aln_fp = os.path.join(args.tmp_dir, f"{group_name}_aligned.fasta")
        alignment = run_mafft(cluster_recs, args.tmp_dir, aln_fp)
        if alignment:
            plot_sequence_logo(alignment, group_name, args.fig_dir)
        else:
            print(f"[WARNING] Failed to align {group_name}, skipping logo.")

    print("[INFO] Pipeline completed successfully.")


if __name__ == "__main__":
    main()