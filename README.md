This repository contains datasets and custom scripts generated and used in the study:

Perkowska et al., “In silico evidence for DOXC21-A NAT-derived small RNA candidates and qPCR of putative target genes in Arabidopsis thaliana.”

The project investigates the regulatory potential of the natural antisense transcript (AT3G19002) overlapping DOXC21-A, 
including structural prediction of miRNA-like candidates, target prediction, Gene Ontology enrichment analysis, and phylogenetic reconstruction of DOXC homologs.

Sequence data
	•	NAT_extended.fa — Full-length NAT sequence (AT3G19002)
	•	nat_mirna_like_candidates.fasta — Candidate NAT-derived miRNA-like sequences
	•	DOXC21_590_sequences_MAFFT_alignment.fasta — Multiple sequence alignment used for phylogenetic analysis

Target prediction and functional analysis
	•	NAT_siRNA_psRNATarget_full_predictions.xlsx — Full psRNATarget output for predicted NAT-derived small RNA targets
	•	NAT_noncano_candidates.xlsx — Filtered miRNA-like candidate list
	•	GO_enrichment_results.csv — Gene Ontology (Biological Process) enrichment analysis results

Custom scripts
	•	detect_noncano_miRNA.py — Script for identifying non-canonical miRNA-like hairpin structures
	•	generate_RNA_structures_VARNA.py — Script for automated RNA secondary structure visualization using VARNA

All analyses were performed in Python (v3.13) using standard scientific libraries (pandas, matplotlib).
