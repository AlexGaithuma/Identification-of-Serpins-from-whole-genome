# Identification-of-Serpins-from-whole-genome
A collection of identification of serpin sequences from whole genome using Reactive center loop sequences from closely identical serpins and matching the names of previously characterized serpins with the new genome identified serpins.

# Usage
1. Create a file with all the RCL sequences from previous studies. Example

       cat Serpin_RCL.fasta \
       45_RCLs.fasta \
       > All_RCL.fasta

2. Use Serpin_annotation.py to Identify and annotate serpin sequences from protein sequences (from genome or other source).
   Example: here we use Chromosome specific input protein sequences. (Optional: You can add the newly identified RCLs to the old and re-run the script. see if you pick any new sequences: Play around with the --window_len option as well)

       for i in {1..16}
       do
       mkdir -p Genome/Chr${i}
       cd Genome/Chr${i}

       python Serpin_annotation.py \
       --proteins_fasta Chr${i}_Annotation/braker.aa \
       --motif_file 45_RCLs.fasta
       --prefix Chr${i} \
       --window_len 100 \
       --verbose
       done

 2. Use Serpin_gtf_filter.py to parse the braker GTF file to filter and extract only the identified serpins.
   Example for 16 chromosomes: we name the sequences with the Chr1, Chr2, Chr3...prefix

       for i in {1..16}
       do
       python3 Serpin_gtf_filter.py \
       --in_GTF Chr${i}_Annotation/braker.gtf \
       --in_PROT_GFF Chr${i}_Identified_serpins.gff \
       --out Chr${i}_filtered_merged_Serpins.gtf
       done

3. Use Find_and_rename_serpins_with_old_names.py to match the previously identified serpins with newly identified serpins.

       cat Chr*_Identified_serpins.fasta >All_serpins_all_chromosomes.fasta
   
       python Match_and_rename_serpins_with_old_names.py \
       --verbose \
       --max_mismatches 6 \
       --in_prot All_serpins_all_chromosomes.fasta \
       --in_motifs 45_RCLs.fasta \
       --out_gtf All_serpins_RCLs_output.gtf \
       --out_protein All_serpins_RCLs_matched_seqs.fasta
