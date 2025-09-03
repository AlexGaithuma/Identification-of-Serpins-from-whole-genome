import os
import tempfile
import subprocess

GENOMIC_DNA = ("ATGAAGCTCGCCCGTGCCCAGAACCTTTTTTCCCTCCAGCTTCTCAAAAAACTCTCCAGCGAAAAGCCAGAATCGAAC"
    "ATCTTCTTCTCGCCGACCAGCATCTCTGTCGCACTCGCCATGGTCTACGCCGGCGCCAGGGGCAAGTCTGAAACGGAA"
    "ATTTCCACTGCTCTCGGCCACACCGCAGCGGGACTGTCCAGCAGAGAATCTATCCTGGAGTCTTACAAAAAGATTCTA"
    "GCAGAACAGCAGACTGATGACAATGTAAGTTTGATGATAGCCAATGCAGTGTTTGTAAAGAAAACACTGAACGTCCTC"
    "GAAAGTTAGCAGAAAGAGCCCGTTGACATCTTTGCGGCAATGTTCCGGCCGGTGGATGTTGGCGCAGGAAATTCTTCC"
    "ATGGAATTGGAGGTTAACGAGTGGCTGAAGAACAAGACCAGAGGCAAGATTTCGGGCTTCAAAATTCCCCCAAACACC"
    "ATAATGGCGCTTCTTAATGACATTTACTTCTAGGGGCTATGGAACACTCTTTCGATCCCAACGATACCTGCATTCTTT"
    "CTTTTCATAACAAAGGATTTGAGGTGGTTATGGTAGAAACCATGACGCGATTGGCTAAGGTACCGTTCACCTCTGAAC"
    "CCGAGTTTGAAGCCATAGAACTATCTTACAAAGGGGACCAACACTGCATGGTAATCATACTTCCAAGTGAAAAAAAGG"
    "AACTTCCCAAACTCCGCGACGCCATGACTGTCGAGAGCATAAAACAAATCCAGAAAAGTTTGAAAGATGAGACTGTCA"
    "AGATTCAACTTCCCAAGTTCGACTTGAAAACAGAGTACGGCCTCGTCCCTGCCGTGAAGAAGATGGGTGTGCGTTCAG"
    "TCTTCTCTGACGCCGATCTGTCAAGAGTCACTGGTGATGGAGTGCTGAGAGTAACCGAAGTCCAGCACAAGGCAGCCA"
    "TTGAGGTCAACGAAGAAGGAACGGTTGCCGCCGGCGCAACAGTAGTGGTAATAGAAAAACGTGCTCCACATTCCTTTG"
    "CCCTGGATAGACCTTTCTTCTTTTACATCCGCGAGAAAGCTACCGGCCGTATGCTATTTTTGGGAGAAGTGCACGCACT"
    "CCCAGCCGCTAAGCCCACCCTTGGTTAG"
)
PROTEIN = ("MKLARAQNLFSLQLLKELSSEKPESNIFFSPTSISVALAMVYAGARGKSETEISTALGHTAAGLSSRESILESYKKIL"
           "AEQQTDDNGLWNASFDPNDTCILSFHNKGFEVVMVETMTRLAKVPFTSEPEFEAIELSYKGDQHCMVIILPSEKKEL"
           "PKLRDAMTIESIKQIQKSLKDETVKIQLPKFDLKTEYGLVPAVKKMGVRSVFSDADLSRVTGDGVLRVTEVHHKAAI"
           "EVNEEGTVAAGATVVVIEKRAPHSFALDRPFLFYIREKATGRMLFLGEVHALPAAKPTLG")

PEPTIDE_REGION = "GLWNASFDPNDTCILSFHNK"

def write_fasta(seq, outfile, name):
    with open(outfile, "w") as f:
        f.write(f">{name}\n")
        for i in range(0, len(seq), 60):
            f.write(seq[i:i+60] + "\n")

def run_exonerate_query(dna_seq, protein_seq, name, min_percent=30, max_percent=90, output_gtf=None):
    with tempfile.TemporaryDirectory() as tmpdir:
        query_fasta = os.path.join(tmpdir, "query.fa")
        target_fasta = os.path.join(tmpdir, "target.fa")
        out_file = output_gtf if output_gtf else os.path.join(tmpdir, "out.gtf")

        write_fasta(protein_seq, query_fasta, name)
        write_fasta(dna_seq, target_fasta, "target")

        for percent in range(max_percent, min_percent-1, -5):  # Lower successively, but also try high first
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
                result = subprocess.run(exo_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except Exception:
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
                continue
            if output_gtf:
                with open(output_gtf, "w") as outf:
                    outf.write(f"chr3\tExonerate\tgene\t{gene_start}\t{gene_end}\t.\t{strand}\t.\tgene_id \"serpin_gene\";\n")
                    for i, (start, end, st) in enumerate(exons):
                        outf.write(f"chr3\tExonerate\texon\t{start}\t{end}\t.\t{st}\t.\tgene_id \"serpin_gene\"; exon_number \"{i+1}\";\n")
                        outf.write(f"chr3\tExonerate\tCDS\t{start}\t{end}\t.\t{st}\t0\tgene_id \"serpin_gene\";\n")
            return percent, exons, gene_start, gene_end, strand
        return None, [], None, None, "+"

def find_peptide_region_in_mapping(exons, gene_start, dna_seq, region_prot):
    """
    For each exon, translate its sequence and check if region_prot is present
    Allow for up to 3 mismatches (edit distance).
    """
    from Bio.Seq import Seq
    def aa_mismatches(a, b):
        # edit/mismatch
        n = min(len(a), len(b))
        return sum(c1 != c2 for c1, c2 in zip(a[:n], b[:n])) + abs(len(a)-len(b))

    # exons: list of (start, end, strand)
    matches = []
    for i, (start, end, st) in enumerate(exons):
        if st == "+":
            seq = dna_seq[start-1:end]
        else:
            seq = str(Seq(dna_seq[start-1:end]).reverse_complement())
        aa = str(Seq(seq).translate())
        # Exact match
        idx = aa.find(region_prot)
        if idx >= 0:
            matches.append(("exact", i+1, start, end, idx, aa))
            continue
        # Fuzzy: try all substrings for low mismatch
        min_mm = len(region_prot)
        for s in range(len(aa) - len(region_prot) + 1):
            frag = aa[s:s+len(region_prot)]
            mm = aa_mismatches(frag, region_prot)
            if mm <= 3: # Allow up to 3 mismatches (tuneable)
                matches.append(("fuzzy", i+1, start, end, s, aa, mm, frag))
    if matches:
        return matches
    return None

def main():
    # Exon map for full protein
    percent, exons, gene_start, gene_end, strand = run_exonerate_query(
            GENOMIC_DNA, PROTEIN, "full_protein", min_percent=30, max_percent=90, output_gtf="Ch3_serpin.gtf")

    print(f"Full protein mapping: identity ~{percent}%, exons: {exons}, gene {gene_start}-{gene_end}, strand {strand}")
    region_hits = find_peptide_region_in_mapping(exons, gene_start, GENOMIC_DNA, PEPTIDE_REGION)
    if region_hits:
        print("-- peptide region GLWNASFDPNDTCILSFHNK found within exon-mapping:")
        for hit in region_hits:
            if hit[0] == "exact":
                print(f"  exact in exon {hit[1]} ({hit[2]}..{hit[3]}) at AA pos {hit[4]}")
            elif hit[0] == "fuzzy":
                print(f"  fuzzy ({hit[6]} mismatches) in exon {hit[1]} ({hit[2]}..{hit[3]}) at AA pos {hit[4]}: {hit[7]}")
    else:
        print("-- peptide region NOT found in any mapped exon, even with up to 3 AA mismatches.")
        # Try mapping peptide region as stand-alone protein
        percent, region_exons, _, _, _ = run_exonerate_query(
            GENOMIC_DNA, PEPTIDE_REGION, "middle_pept", min_percent=30, max_percent=90, output_gtf=None)
        if region_exons and percent:
            print("Middle peptide could be mapped (as separate exon) at %s%% id: %s" % (percent, region_exons))
        else:
            print("No mapping for peptide region found as separate exon at all.")

if __name__ == "__main__":
    main()