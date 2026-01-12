#!/bin/bash

###############################################################################
# serpin_motif_pipeline.sh
#
# PIPELINE: Adaptive, Motif-Based Functional and Structural Clustering of Serpins
#
# PURPOSE & BIOLOGICAL BASIS:
#   This pipeline discovers motifs from all biologically meaningful contexts
#   (N-terminal, RCL, variable, de novo, and full sequence) using MEME,
#   weighs and summarizes them, and clusters sequences using only the most
#   discriminative features. Region/feature selection is data-driven; pipeline
#   is checkpointed to allow resuming on failure or extension.
#
# CONTRIBUTION OF EACH STEP:
#   - Preprocessing: Uniform, robust headers for all relevant input files
#   - Full-length MEME: Discovers all de novo motifs (objective data)
#   - Context Extraction: Generates distinct context windows (RCL, N-term, full, variable,
#       all MEME sites) to maximize relevant motif discovery for grouping
#   - Context Motif Discovery: Runs MEME on every context window FASTA
#   - Motif Matrix: All sequence x motif/context profiles, for global analysis
#   - Motif Clustering: Uses motif/region importance (info content, group entropy, etc.)
#       for data-driven sequence grouping
#   - Motif Summary: Sequence logos, network, group stats, fully interpretable output
#   - Visualization: Heatmaps, motif/group plots, and tSNE/UMAP motif landscapes
#
# USAGE:
#   ./serpin_motif_pipeline.sh \
#      --in_serpins serpins.fasta \
#      --in_aln serpins.aln.fasta \
#      --in_RCL all_rcl_contexts.fasta \
#      --in_nterm all_nterm_contexts.fasta \
#      --out_dir results_dir
###############################################################################

set -euo pipefail

if [[ $# -eq 0 ]]; then
  echo "USAGE:"
  echo "  $0 \\"
  echo "    --in_serpins <full_sequences.fa>"
  echo "    --in_aln <alignment.fa>"
  echo "    --in_RCL <all_rcl_contexts.fa>"
  echo "    --in_nterm <all_nterm_contexts.fa>"
  echo "    --out_dir <result_dir>"
  exit 1
fi

IN_SERPINS=""
IN_RCL=""
IN_NTERM=""
IN_ALN=""
OUT_DIR=""

while [[ $# -gt 0 ]]; do
  key="$1"
  case $key in
    --in_serpins)  IN_SERPINS="$2"; shift ;;
    --in_RCL)      IN_RCL="$2"; shift ;;
    --in_nterm)    IN_NTERM="$2"; shift ;;
    --in_aln)      IN_ALN="$2"; shift ;;
    --out_dir)     OUT_DIR="$2"; shift ;;
    *) echo "Unknown option $1"; exit 1 ;;
  esac
  shift
done

for var in IN_SERPINS IN_RCL IN_NTERM IN_ALN OUT_DIR; do
  [[ -z "${!var}" ]] && echo "Missing required parameter: $var" >&2 && exit 1
done

mkdir -p "$OUT_DIR/tmp" "$OUT_DIR/meme_whole" "$OUT_DIR/meme_context" "$OUT_DIR/tmp/contexts"

#-------------------------------------------------------------------------------
# STEP 1. Preprocessing: Standardize headers for all major files
# Biological purpose:
#    Robust, comparable names (e.g., s1c1_g1610.t1) across all input sources.
#    Prevents mismatch in context or MEME instance tracking downstream.
#-------------------------------------------------------------------------------
step1a="$OUT_DIR/tmp/serpins_short.fasta"
step1b="$OUT_DIR/tmp/rcl_short.fasta"
step1c="$OUT_DIR/tmp/nterm_short.fasta"
step1d="$OUT_DIR/tmp/serpins_short_aln.fasta"

if [[ -s "$step1a" && -s "$step1b" && -s "$step1c" && -s "$step1d" ]]; then
    echo "[1] Preprocessing outputs found. Skipping."
else
    echo "[1] Preprocessing and unifying sequence headers for all inputs..."
    python3 preprocess_headers.py "$IN_SERPINS" "$OUT_DIR/tmp/serpins_short.fasta"
    python3 preprocess_headers.py "$IN_RCL"     "$OUT_DIR/tmp/rcl_short.fasta"
    python3 preprocess_headers.py "$IN_NTERM"   "$OUT_DIR/tmp/nterm_short.fasta"
    python3 preprocess_headers.py "$IN_ALN"     "$OUT_DIR/tmp/serpins_short_aln.fasta"
    echo "    > Outputs: serpins_short.fasta, rcl_short.fasta, nterm_short.fasta, serpins_short_aln.fasta"
fi

#-------------------------------------------------------------------------------
# STEP 2. Motif Discovery (FULL SEQ): Already provided as global MEME results
# Biological purpose:
#    Use all-discovered motifs from unbiased MEME discovery as feature pool for
#    clustering and context selection downstream. Maximizes coverage, minimizes bias.
#-------------------------------------------------------------------------------
step2="$OUT_DIR/meme_whole/meme.xml"
if [[ -s "$step2" ]]; then
    echo "[2] Global MEME xml output found. Skipping redundant MEME run."
else
    echo "[2] Running MEME for full-length motifs (whole sequence set)..."
    meme "$OUT_DIR/tmp/serpins_short.fasta" \
      -oc "$OUT_DIR/meme_whole" \
      -nmotifs 1000 -minw 6 -maxw 30 \
      > "$OUT_DIR/meme_whole/meme.log" 2>&1
    echo "    > Output: $OUT_DIR/meme_whole/meme.xml"
fi

#-------------------------------------------------------------------------------
# STEP 3. Extract rich, unbiased, and statistically driven motif contexts
# Biological purpose:
#    ALL N-terminal, RCL, high-entropy (variable), full-seq, and MEME motif
#    windows (all sites of all MEME motifs) are included. Avoids bias, maximizes
#    motif coverage for potential functional and structural grouping.
#-------------------------------------------------------------------------------
CTX_OUTPUT="$OUT_DIR/tmp/contexts/Nterm_context.fasta"
if [[ -s "$CTX_OUTPUT" ]]; then
    echo "[3] All context windows already extracted. Skipping."
else
    echo "[3] Extracting all context and variable regions from all sources..."
    python3 extract_contexts.py \
        "$OUT_DIR/tmp/serpins_short.fasta" \
        "$OUT_DIR/tmp/rcl_short.fasta" \
        "$OUT_DIR/tmp/nterm_short.fasta" \
        "$OUT_DIR/meme_whole/meme.xml" \
        "$OUT_DIR/tmp/contexts"
fi

#-------------------------------------------------------------------------------
# STEP 4. Run MEME for context-specific motif discovery, and export motifs in MEME format
# Biological purpose:
#    For each context window (e.g., domain or segment), discover motifs using MEME.
#    Each motif file (meme.txt, same as MEME format) will be used as FIMO input
#    in downstream motif presence matrix construction.
#-------------------------------------------------------------------------------

SKIP="yes"
for ctxfasta in "$OUT_DIR"/tmp/contexts/*.fasta; do
    ctxname=$(basename "$ctxfasta" .fasta)
    outdir="$OUT_DIR/meme_context/$ctxname"
    motifmeme="$outdir/${ctxname}_context.meme"
    # Check for complete and non-empty MEME motif file, not just meme.xml
    if [[ ! -s "$outdir/meme.txt" ]]; then
        SKIP="no"
        break
    fi
    # Optionally, also ensure the context.meme exists (for legacy compatibility)
    if [[ ! -s "$motifmeme" ]]; then
        SKIP="no"
        break
    fi
done

if [[ "$SKIP" == "yes" ]]; then
    echo "[4] MEME context motif searches and motif exports already present. Skipping."
else
    echo "[4] Motif discovery on all context windows using MEME and export motifs per context (all *_context.fasta)..."
    for ctxfasta in "$OUT_DIR"/tmp/contexts/*.fasta; do
        ctxname=$(basename "$ctxfasta" .fasta)
        outdir="$OUT_DIR/meme_context/$ctxname"
        mkdir -p "$outdir"
        motifmeme="$outdir/${ctxname}_context.meme"
        # If meme.txt exists and is not empty, skip MEME run
        if [[ -s "$outdir/meme.txt" ]]; then
            echo "    Context $ctxname: meme.txt detected, skipping MEME run."
        else
            echo "    Discovering motifs: $ctxname"
            meme "$ctxfasta" \
                 -oc "$outdir" \
                 -nmotifs 1000 -minw 5 -maxw 20 \
                 > "$outdir/meme.log" 2>&1
        fi
        # Always (re-)export the motif MEME-format file for downstream steps if not present
        if [[ -s "$outdir/meme.txt" ]]; then
            # Use the full meme.txt file as motif source for FIMO
            if [[ ! -s "$motifmeme" ]]; then
                cp "$outdir/meme.txt" "$motifmeme"
                echo "    Exported motifs (MEME format) for context $ctxname: $motifmeme"
            fi
        else
            echo "    WARNING: meme.txt not found or empty for context $ctxname; cannot export motifs!"
        fi
    done
fi

#-------------------------------------------------------------------------------
# STEP 5. Construct sequence-by-motif fingerprint matrix
# Biological purpose:
#    Each sequence is fingerprinted by its context motif content (all, not just classic).
#    The "weight" of each feature will be learned in clustering and summarization.
#-------------------------------------------------------------------------------
step5="$OUT_DIR/motif_matrix.tsv"
if [[ -s "$step5" ]]; then
    echo "[5] Motif occurrence matrix found. Skipping."
else
    echo "[5] Building all-sequence x all-motif occurrence matrix (using FIMO)..."
    python3 build_motif_matrix.py \
        --context_dir "$OUT_DIR/meme_context" \
        --fasta "$OUT_DIR/tmp/serpins_short.fasta" \
        --out "$OUT_DIR/motif_matrix.tsv"
    echo "    > Output: $OUT_DIR/motif_matrix.tsv"
fi

#-------------------------------------------------------------------------------
# STEP 6. Cluster by weighted motif fingerprint
# Biological purpose:
#    Hierarchical clustering, with motif importance dynamically computed
#    from info content, MEME stats, and group-purity, so groupings are
#    optimally functional/structural.
#-------------------------------------------------------------------------------
step6="$OUT_DIR/motif_clusters.csv"
if [[ -s "$step6" ]]; then
    echo "[6] Motif clustering output found. Skipping."
else
    echo "[6] Weighted clustering of motif fingerprints (by info content, diversity)..."
    python3 cluster_by_motifs.py \
        --matrix "$OUT_DIR/motif_matrix.tsv" \
        --out "$OUT_DIR/motif_clusters.csv"
    echo "    > Output: $OUT_DIR/motif_clusters.csv"
fi

#-------------------------------------------------------------------------------
# STEP 7. Motif weighting, summary, and sequence logos/network
# Biological purpose:
#    Quantify each motif's "power" for discriminating groups, by
#    information content (entropy), MEME E-values/width, and group-purity;
#    generate sequence logos and motif-motif network for deep mechanistic/marker interpretation.
#-------------------------------------------------------------------------------
step7a="$OUT_DIR/motifs_summary.tsv"
step7b="$OUT_DIR/motifs.fasta"
if [[ -s "$step7a" && -s "$step7b" ]]; then
    echo "[7] Motif summary, FASTA and logos found. Skipping."
else
    echo "[7] Motif summary: context, info content, group entropy, sequence logos, and motif network..."
    python3 summarize_motifs.py \
        --meme_dir "$OUT_DIR/meme_whole" \
        --context_dir "$OUT_DIR/meme_context" \
        --clusters "$OUT_DIR/motif_clusters.csv" \
        --out_dir "$OUT_DIR"
    echo "    > Outputs: $OUT_DIR/motifs_summary.tsv, $OUT_DIR/motifs.fasta, motifs_logo/"
fi

#-------------------------------------------------------------------------------
# STEP 8. Visualization and interpretation
# Biological purpose:
#    Motif/group matrix-weighted heatmaps, logos, importance, and network:
#    visual evidence of discriminative/marker motifs and their sequence/group associations.
#-------------------------------------------------------------------------------
step8a="$OUT_DIR/motif_heatmap.png"
step8b="$OUT_DIR/motif_info_content.png"
if [[ -s "$step8a" && -s "$step8b" ]]; then
    echo "[8] Visualization outputs found. Skipping."
else
    echo "[8] Motif/group visualization (heatmap, info, logos, tSNE)..."
    #Rscript visualize_motifs.R \
    #    --summary "$OUT_DIR/motifs_summary.tsv" \
    #    --matrix "$OUT_DIR/motif_matrix.tsv" \
    #    --clusters "$OUT_DIR/motif_clusters.csv" \
    #    --out_dir "$OUT_DIR"
    echo "    > Outputs: $OUT_DIR/motif_heatmap.png, $OUT_DIR/motif_info_content.png, ... "
fi

echo
echo "========================================================================="
echo " Motif-driven, adaptive serpin analysis pipeline finished!"
echo "   All output in:   $OUT_DIR"
echo "   Steps/progress checkpointed for restartability."
echo "   See motifs_summary.tsv, motif_heatmap.png, and logos/network in motifs_logo/"
echo "========================================================================="