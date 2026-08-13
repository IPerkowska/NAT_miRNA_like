DOXC21-A NAT-derived small RNA analysis

This repository contains datasets and custom scripts generated and used in the study:

Perkowska et al., “In silico evidence for small RNA candidates derived from the DOXC21-A NAT locus in Arabidopsis thaliana.”

The study investigates the potential regulatory function of the natural antisense transcript (NAT; AT3G19002) overlapping DOXC21-A (AT3G19000) in Arabidopsis thaliana. The analyses integrate sequence conservation, RNA secondary structure prediction, reanalysis of publicly available small RNA-seq data, target prediction, functional enrichment, and phylogenetic analysis of DOXC21-related proteins.

Sequence data

* NAT_extended.fa — Full-length sequence of the AT3G19002 NAT
* nat_mirna_like_candidates.fasta — Candidate NAT-derived small RNA sequences
* DOXC21_590_sequences_MAFFT_alignment.fasta — Multiple sequence alignment of DOXC21-related proteins used for phylogenetic analysis

Small RNA and target analysis

* NAT_siRNA_psRNATarget_full_predictions.xlsx — Full psRNATarget output for predicted targets of NAT-derived small RNA candidates
* NAT_noncano_candidates.xlsx — Filtered list of NAT-derived small RNA candidates
* GO_enrichment_results.csv — Gene Ontology enrichment results for predicted target genes

Custom scripts

* detect_noncano_miRNA.py — Identification and filtering of non-canonical NAT-derived small RNA candidates
* generate_RNA_structures_VARNA.py — Automated visualization of predicted RNA secondary structures using VARNA

The small RNA analysis included reanalysis of publicly available A. thaliana small RNA-seq libraries from Phytophthora capsici-infected plants to examine reads mapping to the AT3G19002 NAT locus.

Custom analyses were performed in Python using standard scientific libraries, including pandas and matplotlib.
