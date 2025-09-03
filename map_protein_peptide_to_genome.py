#!/usr/bin/env python3

"""
=================================================================================
  Script:   protein_peptide_genome_mapping.py

  PURPOSE:
    1. Map a full LENGTH protein sequence to genomic DNA using Exonerate.
    2. For all exons found, check whether a PEPTIDE (sub-protein region) occurs in coding region(s)
       (allowing up to 3 AA mismatches).
    3. If not, try mapping the peptide as an independent region on DNA.
    4. Print out step-by-step what is attempted, what is found, and print warnings if not found.

  USAGE:
    python protein_peptide_genome_mapping.py --Genomic_DNA <genomic_dna.fa> --in_protein <protein.fa> --in_peptide_region <peptide.fa>

  REQUIREMENTS:
    * BioPython (for reading FASTA and basic sequence types)
    * Exonerate on your $PATH

=================================================================================
"""

import os
import sys
import tempfile
import subprocess
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def verbose(msg):
    """Print detailed step-by-step messages."""
    print(f"[INFO] {msg}", flush=True)

def warn(msg):
    """Print warnings to the user."""
    print(f"[WARNING] {msg}", flush=True)

def load_single_fasta(filepath, seqtype_descr):
    """
    Load the FIRST record from a fasta as a string; exit if file/record missing.
    """
    if not os.path.exists(filepath):
        warn(f"Input file not found: {filepath}")
        sys.exit(1)
    try:
        with open(filepath) as f:
            records = list(SeqIO.parse(f, "fasta"))
            if not records:
                warn(f"File '{filepath}' appears empty or non-FASTA.")
                sys.exit(1)
            if len(records) > 1:
                warn(f"File '{filepath}' contains more than one sequence record; ONLY the first will be used ({seqtype_descr}).")
            seq = str(records[0].seq)
            verbose(f"Loaded {seqtype_descr} of length {len(seq)} from '{filepath}'")
            return seq
    except Exception as e:
        warn(f"Failed to parse '{filepath}' as FASTA: {e}")
        sys.exit(1)

def write_fasta(seq, outfile, name):
    """
    Write a sequence (string) to FASTA file.
    """
    with open(outfile, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")
    verbose(f"Wrote sequence '{name}' to FASTA file {outfile}")

def run_exonerate_query(dna_seq, protein_seq, name, min_percent=30, max_percent=90, output_gtf=None):
    """
    STEP 1: Map protein or peptide to DNA with Exonerate, searching for exons.

    Args:
        dna_seq, protein_seq: strings
        name: label for mapping
        min_percent: lowest identity threshold
        max_percent: highest identity threshold to try
        output_gtf: optional, output file for GTF features
    Returns:
        percent: which threshold mapping succeeded at, or None
        exons: list of (start,end,strand) for coding exons found
        gene_start, gene_end, strand
    -------------------------------------------------------------
    * Tries a range of identities low-high
    * Reports mapping step, and writes any found mapping GTF if requested
    """
    verbose("STEP 1: Running Exonerate to map (protein/peptide) to genomic DNA.")
    with tempfile.TemporaryDirectory() as tmpdir:
        query_fasta = os.path.join(tmpdir, "query.fa")
        target_fasta = os.path.join(tmpdir, "target.fa")
        out_file = output_gtf if output_gtf else os.path.join(tmpdir, "out.gtf")
        write_fasta(protein_seq, query_fasta, name)
        write_fasta(dna_seq, target_fasta, "target")
        for percent in range(max_percent, min_percent-1, -5):
            verbose(f"Trying Exonerate for mapping at minimum identity {percent}%...")
            exo_cmd = [
                "exonerate", "--model", "protein2genome",
                "--showtargetgff", "yes",
                "--showalignment", "no",
                "--showvulgar", "no",
                "--showcigar", "no",
                "--ryo", "",
                "--query", query_fasta,
                "--target", target_fasta,
                "--percent", str(percent)
            ]
            try:
                result = subprocess.run(
                    exo_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
                )
            except Exception as e:
                warn(f"Exonerate failed at identity {percent}%: {e}")
                continue
            exons = []
            gene_start, gene_end, strand = None, None, "+"
            for line in result.stdout.splitlines():
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                if fields[2] == "gene":
                    gene_start = int(fields[3])
                    gene_end = int(fields[4])
                    strand = fields[6]
                elif fields[2].lower() == "cds":
                    exons.append((int(fields[3]), int(fields[4]), fields[6]))
            if not exons:
                verbose(f"No exons mapped at {percent}%, will try a lower threshold.")
                continue
            if output_gtf:
                with open(output_gtf, "w") as outf:
                    outf.write(f"chr3\tExonerate\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tgene_id \"serpin_gene\";\n")
                    for i, (start, end, st) in enumerate(exons):
                        outf.write(f"chr3\tExonerate\texon\t{start}\t{end}\t.\t{st}\t.\tgene_id \"serpin_gene\"; exon_number \"{i+1}\";\n")
                        outf.write(f"chr3\tExonerate\tCDS\t{start}\t{end}\t.\t{st}\t0\tgene_id \"serpin_gene\";\n")
                verbose(f"Wrote mapped gene and exon data to '{output_gtf}'")
            verbose(f"Mapping SUCCEEDED at {percent}%: {len(exons)} exons found ({gene_start}-{gene_end}, {strand})")
            return percent, exons, gene_start, gene_end, strand
        warn("Exonerate mapping failed for all identity thresholds. No mapping produced.")
        return None, [], None, None, "+"

def find_peptide_region_in_mapping(exons, gene_start, dna_seq, region_prot):
    """
    STEP 2: Search for the peptide region in translated exons from mapped protein

    Args:
        exons: list of (start, end, strand)
        gene_start: unused, for reference only
        dna_seq: string, entire genomic sequence
        region_prot: string, peptide
    Returns:
        matches: list of ("exact"/"fuzzy", exon#, nt_start, nt_end, aa_pos, <extra>)
    """
    verbose("STEP 2: Scanning all mapped exons for your peptide region (exact or ≤3 mismatches).")
    def aa_mismatches(a, b):
        n = min(len(a), len(b))
        return sum(c1 != c2 for c1, c2 in zip(a[:n], b[:n])) + abs(len(a)-len(b))
    matches = []
    for i, (start, end, st) in enumerate(exons):
        # Extract & translate exon sequence according to strand
        seq = dna_seq[start-1:end] if st == "+" else str(Seq(dna_seq[start-1:end]).reverse_complement())
        aa = str(Seq(seq).translate())
        verbose(f"Examining exon {i+1} ({start}-{end}, strand {st}, {len(aa)} AA residues).")
        idx = aa.find(region_prot)
        if idx >= 0:
            verbose(f"  Peptide region found EXACTLY in exon {i+1} at AA pos {idx}.")
            matches.append(("exact", i+1, start, end, idx, aa))
            continue
        # Fuzzy: check for windows with ≤3 mismatches
        for s in range(len(aa) - len(region_prot) + 1):
            frag = aa[s:s+len(region_prot)]
            mm = aa_mismatches(frag, region_prot)
            if mm <= 3:
                verbose(f"  Fuzzy match in exon {i+1} at AA pos {s}: {frag} (mismatches={mm})")
                matches.append(("fuzzy", i+1, start, end, s, aa, mm, frag))
    if matches:
        verbose("Peptide region (or similar) FOUND in mapped exons.")
        return matches
    else:
        verbose("Peptide region NOT found in mapped exons (even with ≤3 mismatches).")
        return None

def main():
    # -------------------------------
    # Command-line argument handling
    # -------------------------------
    parser = argparse.ArgumentParser(
        description="""
  protein_peptide_genome_mapping.py

  - Map full protein to genome and check presence/location of peptide region within coding mapping.
  - Peptide and protein must be in single-record FASTA files.
  - All input files required as options.

  Example run:
    python protein_peptide_genome_mapping.py --Genomic_DNA genome.fa --in_protein protein.fa --in_peptide_region peptide.fa
        """
    )
    parser.add_argument('--Genomic_DNA', required=True, help="Input: genomic DNA, FASTA (single record)")
    parser.add_argument('--in_protein', required=True, help="Input: protein sequence, FASTA (single record)")
    parser.add_argument('--in_peptide_region', required=True, help="Input: peptide region, FASTA (single record)")
    args = parser.parse_args()

    # -----------------------------
    # STEP 0: Sequence Loading
    # -----------------------------
    GENOMIC_DNA = load_single_fasta(args.Genomic_DNA, "genomic DNA")
    PROTEIN     = load_single_fasta(args.in_protein, "protein")
    PEPTIDE_REGION = load_single_fasta(args.in_peptide_region, "peptide region")

    # ---------------------------------------------------------------
    # STEP 1: Map full protein to genome
    # ---------------------------------------------------------------
    percent, exons, gene_start, gene_end, strand = run_exonerate_query(
            GENOMIC_DNA, PROTEIN, "full_protein", min_percent=30, max_percent=90, output_gtf="Ch3_serpin.gtf")
    if percent is not None and exons:
        verbose(f"Full-protein mapping identity: {percent}%. {len(exons)} exons found, gene region {gene_start}-{gene_end} ({strand})")
    else:
        warn("No exons or gene region mapped for the protein! Check input or try lowering identity required.")
        sys.exit(2)

    # ---------------------------------------------------------------
    # STEP 2: Search all mapped exons for peptide region
    # ---------------------------------------------------------------
    region_hits = find_peptide_region_in_mapping(exons, gene_start, GENOMIC_DNA, PEPTIDE_REGION)

    # Output results
    pep_desc = "-- Your peptide region"
    if region_hits:
        print(f"{pep_desc} found INSIDE mapped exonic regions of the protein:")
        for hit in region_hits:
            if hit[0] == "exact":
                print(f"  [EXACT] In exon {hit[1]} (nucleotides {hit[2]}..{hit[3]}) at AA pos {hit[4]}")
            elif hit[0] == "fuzzy":
                print(f"  [FUZZY] ({hit[6]} mismatches) in exon {hit[1]} ({hit[2]}..{hit[3]}) at AA pos {hit[4]}: {hit[7]}")
    else:
        print(f"{pep_desc} NOT found in protein mapped exons (even with up to 3 mismatches).")
        # ---------------------------------------------------------
        # STEP 3: Try peptide region alone as a standalone mapping
        # ---------------------------------------------------------
        verbose("Peptide not found; mapping peptide as a stand-alone protein query...")
        percent2, region_exons, _, _, _ = run_exonerate_query(
            GENOMIC_DNA, PEPTIDE_REGION, "peptide_query", min_percent=30, max_percent=90, output_gtf=None)
        if region_exons and percent2:
            print(f"Peptide region independently mapped (as separate exon) at {percent2}% identity: {region_exons}")
        else:
            print("Peptide region NOT mapped anywhere on genome, even as stand-alone, even at low identity.")

    verbose("Step-by-step mapping and checking complete. Exiting.")

if __name__ == "__main__":
    main()