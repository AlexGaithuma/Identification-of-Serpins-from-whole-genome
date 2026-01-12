#!/usr/bin/env python3
"""
===============================================================================
Script: cluster_by_motifs.py

Purpose:
    Cluster biological sequences (samples/rows) according to motif occurrence profiles
    (motif features/columns). Motif columns are weighted by informativeness, using
    feature importance from classification between initial clusters, and re-clustered
    with emphasis on most group-distinguishing motifs.

Biological basis:
 - Motif content (sequence features, e.g., k-mers, Pfam signatures, species markers)
   is a fundamental biological fingerprint in genomics and metagenomics.
 - Functionally, certain motifs better distinguish biological classes (gene families,
   subpopulations, or taxa). Motifs with high variance/diversity or high
   cluster-separability are enriched for distinguishing biology.
 - This method discovers latent structure: (1) clusters samples unsupervised; (2)
   identifies motifs which explain these cluster differences (via supervised learning but
   using unsupervised cluster labels); (3) re-clusters with motifs reweighted to emphasize
   informative signals. The final groups should be better aligned with underlying biology.
 - This is *unsupervised*: requires no prior group/label information.

Pipeline Integration:
 - Designed for automated workflows. Checks input/output. Integrates with bash/cwl
   pipelines (step 6). Produces output for downstream analysis or visualization.

===============================================================================
"""

import argparse
import pandas as pd
import numpy as np
import sys
import os

from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering

# ----------------------------- Argument Parsing ------------------------------

parser = argparse.ArgumentParser(
    description="Cluster biological sequences based on motif occurrence profiles, "
                "using motif informativeness scores (by Random Forest) to weight "
                "motif columns for final cluster assignments."
)
parser.add_argument("--matrix", required=True,
    help="Input motif occurrence matrix (tab-delimited; rows: samples/sequences; "
         "columns: motifs; first column: sequence/sample IDs).")
parser.add_argument("--out", required=True,
    help="Output CSV: cluster/group assignments for all samples.")
args = parser.parse_args()

# ----------------------------- Data Loading ----------------------------------

# Expects:
# - tab-delimited input file
# - row 0: header (motif column names)
# - column 0: sample/sequence identifiers

if not os.path.exists(args.matrix):
    sys.exit(f"[ERROR] Input file '{args.matrix}' does not exist. Please check the path.")

try:
    mat = pd.read_csv(args.matrix, sep='\t', index_col=0)
except Exception as e:
    sys.exit(f"[ERROR] Failed to read the motif matrix: {e}\n"
             "Is the file properly tab-delimited with header and sequence/sample IDs?")

# Validate shape: at least 2 samples and 1 motif column
if mat.shape[0] < 2:
    sys.exit(f"[ERROR] Input matrix must have >=2 sequences (rows). Rows found: {mat.shape[0]}")
if mat.shape[1] < 1:
    sys.exit(f"[ERROR] Input matrix contains no motif features (columns). "
             "Check for proper header, sample ID column, and tab separation.")

print(f"[INFO] Loaded motif occurrence matrix: {mat.shape[0]} sequences, {mat.shape[1]} motif columns.")

# ----------------------------- Feature Scaling -------------------------------

# Motif abundance can span different value ranges (e.g., rare vs. frequent motifs).
# Standardization (zero mean, unit variance) ensures fair comparison and prevents
# highly-abundant motifs biasing clustering outcome.
scaler = StandardScaler()
X_std = scaler.fit_transform(mat.values)

# ----------------------------- Initial Clustering ----------------------------

# First, assign clusters in an unbiased, unsupervised manner.
# Using hierarchical clustering (Agglomerative), with number of clusters heuristic:
# - n_clusters = 2 + log(number of samples), giving moderate granularity.
n_clusters = max(2, 2 + int(np.floor(np.log(len(mat)))))
clustering = AgglomerativeClustering(n_clusters=n_clusters)
provisional_labels = clustering.fit_predict(X_std)

print(f"[INFO] Step 1: Initial clustering complete ({n_clusters} clusters assigned).")

# ----------------------------- Motif Informativeness Scoring -----------------

# Next, measure each motif's ability to distinctively separate the above clusters.
# Biological rationale: motifs best explaining observed clusters are more likely to
# reflect real biological discrimination (families, groups, classes, etc).
# We quantify each motif's cluster-separating power via feature importance in
# a Random Forest classifier (classification: motif profile => cluster label).

rf = RandomForestClassifier(
    n_estimators=200,
    random_state=1,
    class_weight='balanced_subsample', # Guards against cluster size imbalance
    n_jobs=-1
)
rf.fit(mat.values, provisional_labels)
motif_importances = pd.Series(rf.feature_importances_, index=mat.columns)

print("[INFO] Step 2: Top 5 most informative motifs (for separating clusters):")
print(motif_importances.sort_values(ascending=False).head())

# ----------------------------- Weighted Clustering ---------------------------

# To emphasize informative motifs, apply their informativeness scores as weights:
# Multiply each column by its feature importance, upweighting discriminative patterns,
# downweighting background/noise motifs. This sharpens biological signal in
# group assignment.
X_weighted = mat.values * motif_importances.values

# Redo clustering with the motif-weighted matrix, optionally using 'average' linkage
# (helps robustness to noise/outliers).
final_clustering = AgglomerativeClustering(
    n_clusters=n_clusters,
    linkage='average'
)
final_labels = final_clustering.fit_predict(X_weighted)

# Silhouette score provides an unsupervised quality assessment (<1 best, >0 expected)
try:
    silh = silhouette_score(X_weighted, final_labels)
    print(f"[INFO] Step 3: Final clustering silhouette score = {silh:.3f} (higher: better cohesion/separation)")
except Exception as e:
    print(f"[WARNING] Could not compute silhouette score: {e}")

# ----------------------------- Output ----------------------------------------

# Output: CSV with columns [SampleID, GroupID, GroupName]
# - SampleID: as input
# - GroupID: numeric cluster index (for downstream filtering/naming)
# - GroupName: user-friendly cluster (e.g., "Grp1", "Grp2" ...)
output_df = pd.DataFrame({
    "SampleID": mat.index,
    "GroupID": final_labels,
    "GroupName": [f"Grp{gid+1}" for gid in final_labels]
})

try:
    output_df.to_csv(args.out, index=False)
    print(f"[SUCCESS] Step 4: Motif cluster assignments written to: {args.out}")
except Exception as e:
    sys.exit(f"[ERROR] Could not write output file '{args.out}': {e}")

# ----------------------------- End of Script ---------------------------------