#!/bin/bash
###############################################################################
# Trichoderma Protein Homology Evidence Preparation Pipeline
# Author: Dr. Alfiky
# Description: Shell commands used to download, curate, and deduplicate 
#              protein sequences from NCBI RefSeq for Trichoderma species.
# Requirements: NCBI Entrez Direct, seqkit, CD-HIT, BLAST+
# Tested on: Ubuntu Linux
###############################################################################

# ---------------------------------------------------------------------------
# 1. Install NCBI Entrez Direct toolkit
# ---------------------------------------------------------------------------
sudo apt update
sudo apt install -y ncbi-entrez-direct

# ---------------------------------------------------------------------------
# 2. Download curated RefSeq protein sequences for a specific Trichoderma species
#    Example shown for Trichoderma harzianum
# ---------------------------------------------------------------------------
esearch -db protein -query "Trichoderma harzianum[Organism] AND srcdb_refseq[PROP]" | \
  efetch -format fasta > T_harzianum_refseq.faa

# ---------------------------------------------------------------------------
# 3. Confirm retrieval success and check sequence count
#    Replace 'XX' with the species name abbreviation
# ---------------------------------------------------------------------------
cd ~/path/
head T_XX_refseq.faa
grep -c "^>" T_XX_refseq.faa

# ---------------------------------------------------------------------------
# 4. Filter sequences by minimum length (100 aa) and remove fragment annotations
# ---------------------------------------------------------------------------
seqkit seq -m 100 T_XX_refseq.faa | \
  seqkit grep -v -r -p "fragment" > T_XX.cleaned.faa

# ---------------------------------------------------------------------------
# 5. Deduplicate sequences using CD-HIT at 99% identity
# ---------------------------------------------------------------------------
cd-hit -i T_XX.cleaned.faa -o T_XX.nr.faa -c 0.99 -n 5 -T 4 -M 16000

# ---------------------------------------------------------------------------
# 6. Summarize final curated protein dataset
# ---------------------------------------------------------------------------
seqkit stats T_species.nr.faa

# ---------------------------------------------------------------------------
# 7. Example for Trichoderma longibrachiatum filtering
# ---------------------------------------------------------------------------
seqkit seq -m 101 -r -p /path/T_longibrachiatum_protein_filtered_cdhit.faa | \
  seqkit grep -v -r -p "fragment" \
  > /path/T_longibrachiatum_filtered.faa

# ---------------------------------------------------------------------------
# 8. Download all Trichoderma longibrachiatum proteins from NCBI (non-RefSeq)
# ---------------------------------------------------------------------------
esearch -db protein -query "Trichoderma longibrachiatum[Organism]" | \
  efetch -format fasta > ~/path/T_longibrachiatum_all.faa

# ---------------------------------------------------------------------------
# 9. Remove redundancy and retain high-confidence sequences
# ---------------------------------------------------------------------------
cd-hit -i T_longibrachiatum_all.faa -o T_longibrachiatum_nr.faa -c 0.99 -n 5 -d 0 -M 16000 -T 4

# ---------------------------------------------------------------------------
# 10. Count resulting non-redundant sequences
# ---------------------------------------------------------------------------
grep -c '^>' ~/path/T_longibrachiatum_nr.faa

# ---------------------------------------------------------------------------
# 11. Retrieve all Trichoderma RefSeq proteins into one FASTA file
# ---------------------------------------------------------------------------
esearch -db protein -query "Trichoderma[Organism] AND srcdb_refseq[PROP]" | \
  efetch -format fasta > /path/RefSeq_Trichoderma_proteins.faa

# ---------------------------------------------------------------------------
# 12. Build a BLAST protein database
# ---------------------------------------------------------------------------
makeblastdb -in /path/RefSeq_Trichoderma_proteins.faa \
  -dbtype prot \
  -out /path/RefSeq_Trichoderma_proteins_db \
  -parse_seqids

# ---------------------------------------------------------------------------
# 13. Run BLASTP to identify supported sequences for T. longibrachiatum
# ---------------------------------------------------------------------------
blastp -query /path/T_longibrachiatum_filtered.faa \
  -db /path/RefSeq_Trichoderma_proteins_db \
  -outfmt 6 -evalue 1e-5 -max_target_seqs 1 -num_threads 10 \
  -out /path/T_longibrachiatum_vs_RefSeq.tsv

# ---------------------------------------------------------------------------
# 14. Extract supported sequence IDs
# ---------------------------------------------------------------------------
cut -f1 /path/T_longibrachiatum_vs_RefSeq.tsv | \
  sort | uniq > /path/supported_ids.txt

# ---------------------------------------------------------------------------
# 15. Retrieve supported sequences using seqkit
# ---------------------------------------------------------------------------
seqkit grep -f /path/supported_ids.txt \
  /path/T_longibrachiatum_filtered.faa \
  > /path/T_longibrachiatum_supported.faa

# ---------------------------------------------------------------------------
# 16. Count number of retained proteins
# ---------------------------------------------------------------------------
grep -c '^>' /path/T_longibrachiatum_supported.faa

# ---------------------------------------------------------------------------
# 17. Concatenate all FASTA protein sets (adjust filenames as needed)
# ---------------------------------------------------------------------------
cat *.faa > pooled_proteins_raw.faa

# ---------------------------------------------------------------------------
# 18. Deduplicate the pooled dataset at 99% identity
# ---------------------------------------------------------------------------
cd-hit -i pooled_proteins_raw.faa -o pooled_proteins_cdhit99.faa -c 0.99 -n 5 -T 10 -M 16000

# ---------------------------------------------------------------------------
# 19. Count sequences before and after deduplication
# ---------------------------------------------------------------------------
grep -c '^>' /path/pooled_Trichoderma_proteins.faa
grep -c '^>' /path/pooled_Trichoderma_cdhit99.faa

###############################################################################
# End of Script
###############################################################################
