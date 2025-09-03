#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict

def parse_protein_gff(gff_file):
    """
    Parses protein GFF to:
      - collect transcript IDs
      - collect RCL protein residue coordinates per transcript
    """
    transcript_ids = set()
    rcl_coords = {}

    with open(gff_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue
            seqid = fields[0]
            ftype = fields[2]
            start = int(fields[3])
            end = int(fields[4])

            transcript_ids.add(seqid)
            if ftype == "RCL":
                rcl_coords[seqid] = (start, end)

    return transcript_ids, rcl_coords

def parse_gtf_and_check_strand(gtf_file, transcripts_of_interest):
    """
    Parses GTF into transcript-wise features with strand checks.
    Returns:
      transcripts: dict transcript_id -> dict:
        {
          'chrom': str,
          'strand': '+' or '-',
          'gene_id': str,
          'features': list of full feature lines (list of fields)
        }
    Filters out transcripts with inconsistent strand information.
    """
    transcripts = {}
    strand_mismatch = set()
    transcript_strands = defaultdict(set)

    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            chrom = fields[0]
            feat_type = fields[2]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attr_text = fields[8]

            transcript_id = None
            gene_id = None
            attrs = [a.strip() for a in attr_text.strip().split(";") if a.strip()]
            for a in attrs:
                if a.startswith("transcript_id"):
                    parts = a.split()
                    if len(parts) > 1:
                        transcript_id = parts[1].strip('"')
                elif a.startswith("gene_id"):
                    parts = a.split()
                    if len(parts) > 1:
                        gene_id = parts[1].strip('"')

            if transcript_id is None or transcript_id not in transcripts_of_interest:
                continue

            transcript_strands[transcript_id].add(strand)
            if transcript_id in strand_mismatch:
                continue

            if len(transcript_strands[transcript_id]) > 1:
                strand_mismatch.add(transcript_id)
                sys.stderr.write(f"Warning: transcript {transcript_id} has features on multiple strands. Skipping.\n")
                continue

            if transcript_id not in transcripts:
                transcripts[transcript_id] = {
                    "chrom": chrom,
                    "strand": strand,
                    "gene_id": gene_id if gene_id else transcript_id,
                    "features": []
                }
            transcripts[transcript_id]["features"].append(fields)

    for tx in strand_mismatch:
        if tx in transcripts:
            del transcripts[tx]

    return transcripts

def protein_residue_to_genomic(protein_start, protein_end, cds_exons, strand):
    cds_nt_start = (protein_start - 1) * 3 + 1
    cds_nt_end = protein_end * 3

    curr_len = 0
    genomic_start = None
    genomic_end = None

    for exon_start, exon_end in cds_exons:
        exon_len = exon_end - exon_start + 1
        next_len = curr_len + exon_len

        if genomic_start is None and (cds_nt_start > curr_len and cds_nt_start <= next_len):
            offset = cds_nt_start - curr_len - 1
            genomic_start = exon_start + offset if strand == '+' else exon_end - offset

        if genomic_end is None and (cds_nt_end > curr_len and cds_nt_end <= next_len):
            offset = cds_nt_end - curr_len - 1
            genomic_end = exon_start + offset if strand == '+' else exon_end - offset

        curr_len = next_len
        if genomic_start is not None and genomic_end is not None:
            break

    if genomic_start is None or genomic_end is None:
        raise ValueError(f"Protein residue coords {protein_start}-{protein_end} don't map in CDS exons.")

    return (genomic_start, genomic_end) if genomic_start < genomic_end else (genomic_end, genomic_start)

def extract_cds_exons(transcript_features):
    cds_exons = []
    for fields in transcript_features:
        if fields[2].lower() == "cds":
            cds_exons.append((int(fields[3]), int(fields[4])))
    if not cds_exons:
        raise ValueError("No CDS features found for transcript")
    strand = transcript_features[0][6]
    if strand == '+':
        cds_exons.sort(key=lambda x: x[0])
    else:
        cds_exons.sort(key=lambda x: x[0], reverse=True)
    return cds_exons

def attributes_to_dict(attr_text):
    d = {}
    for attr in attr_text.strip().split(";"):
        attr = attr.strip()
        if not attr:
            continue
        if " " in attr:
            k, v = attr.split(" ", 1)
            d[k] = v.strip('"')
        else:
            d[attr] = ""
    return d

def dict_to_attributes(d):
    out = []
    for k, v in d.items():
        if v == "":
            out.append(k)
        else:
            out.append(f'{k} "{v}"')
    return "; ".join(out) + ";"

def write_transcript(transcript_data, rcl_coords):
    chrom = transcript_data["chrom"]
    strand = transcript_data["strand"]
    gene_id = transcript_data["gene_id"]
    transcript_id = None
    features = transcript_data["features"]

    for f in features:
        attr_dict = attributes_to_dict(f[8])
        if "transcript_id" in attr_dict:
            transcript_id = attr_dict["transcript_id"]
            break
    if not transcript_id:
        transcript_id = "undefined"

    # Group features by type with ordering priority
    # Desired order: start_codon, mRNA, CDS, exon, motif_search (RCL), stop_codon

    # Initialize buckets
    buckets = {
        "start_codon": [],
        "mRNA": [],
        "CDS": [],
        "exon": [],
        "stop_codon": [],
    }
    others = []

    for f in features:
        ftype = f[2]
        if ftype == "motif_search":
            # skip motif_search as it corresponds to RCL feature, handled later
            continue
        elif ftype in buckets:
            buckets[ftype].append(f)
        else:
            others.append(f)

    # Sort each bucket by start coordinate ascending
    for k in buckets:
        buckets[k].sort(key=lambda x: int(x[3]))
    others.sort(key=lambda x: int(x[3]))

    out_lines = []

    def fix_attrs(fields):
        attr_dict = attributes_to_dict(fields[8])
        ftype = fields[2]

        attr_dict["gene_id"] = gene_id
        # Modify transcript_id based on type:
        if ftype == "start_codon":
            attr_dict["transcript_id"] = f"start_codon_{transcript_id}"
        elif ftype == "mRNA":
            attr_dict["transcript_id"] = f"mRNA_{transcript_id}"
        elif ftype == "CDS":
            attr_dict["transcript_id"] = f"CDS_{transcript_id}"
        elif ftype == "exon":
            attr_dict["transcript_id"] = f"exon_{transcript_id}"
        elif ftype == "stop_codon":
            attr_dict["transcript_id"] = f"stop_codon_{transcript_id}"
        else:
            if "transcript_id" not in attr_dict or not attr_dict["transcript_id"]:
                attr_dict["transcript_id"] = transcript_id

        fields[8] = dict_to_attributes(attr_dict)

    # Write start_codon, mRNA, CDS, exon in order
    for ftype in ["start_codon", "mRNA", "CDS", "exon"]:
        for f in buckets[ftype]:
            fix_attrs(f)
            out_lines.append("\t".join(f) + "\n")

    # Append motif_search (RCL) feature line if available after CDS and exon but before stop_codon
    if transcript_id in rcl_coords:
        prot_start, prot_end = rcl_coords[transcript_id]
        try:
            cds_exons = extract_cds_exons(features)
            g_start, g_end = protein_residue_to_genomic(prot_start, prot_end, cds_exons, strand)
        except Exception as e:
            sys.stderr.write(f"Warning: could not add RCL feature for {transcript_id}: {e}\n")
        else:
            rcl_attr = dict_to_attributes({
                "gene_id": f"RCL_{transcript_id}",
                "transcript_id": f"RCL_{transcript_id}"
            })
            rcl_fields = [
                chrom,
                "motif_search",
                "RCL",
                str(g_start),
                str(g_end),
                ".",
                strand,
                ".",
                rcl_attr
            ]
            out_lines.append("\t".join(rcl_fields) + "\n")

    # Write stop_codon last
    for f in buckets["stop_codon"]:
        fix_attrs(f)
        out_lines.append("\t".join(f) + "\n")

    # Write any other features at the end
    for f in others:
        attr_dict = attributes_to_dict(f[8])
        attr_dict["gene_id"] = gene_id
        if "transcript_id" not in attr_dict or not attr_dict["transcript_id"]:
            attr_dict["transcript_id"] = transcript_id
        f[8] = dict_to_attributes(attr_dict)
        out_lines.append("\t".join(f) + "\n")

    return out_lines

def main():
    parser = argparse.ArgumentParser(description="Filter and process Serpin transcripts and add RCL features with proper strand consistency and standard GTF formatting.")
    parser.add_argument("--in_GTF", required=True, help="Input GTF file with chromosome features")
    parser.add_argument("--in_PROT_GFF", required=True, help="Input protein GFF with RCL feature protein coordinates")
    parser.add_argument("--out", required=True, help="Output filtered and corrected GTF file")

    args = parser.parse_args()

    transcript_ids, rcl_coords = parse_protein_gff(args.in_PROT_GFF)
    if not transcript_ids:
        sys.exit("No transcripts found in protein GFF.")

    transcripts = parse_gtf_and_check_strand(args.in_GTF, transcript_ids)
    if not transcripts:
        sys.stderr.write("No valid transcripts with consistent strand found in GTF.\n")
        sys.exit(1)

    def min_start(trans):
        starts = [int(f[3]) for f in trans['features']]
        return min(starts)

    sorted_tx = sorted(transcripts.values(), key=lambda x: (x["chrom"], min_start(x)))

    with open(args.out, "w") as out_f:
        out_f.write('##gtf-version 2\n')
        for transcript_data in sorted_tx:
            lines = write_transcript(transcript_data, rcl_coords)
            for l in lines:
                out_f.write(l)

if __name__ == "__main__":
    main()