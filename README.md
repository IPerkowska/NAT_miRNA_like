# DOXC21-A NAT-derived small RNA analysis

This repository contains datasets and custom scripts generated and used in the study:

**Perkowska et al., “Integrated in silico and experimental analysis of a stress-responsive DOXC21-A natural antisense transcript suggests small RNA-associated regulation in Arabidopsis thaliana.”**

The study investigates the potential regulatory role of the natural antisense transcript (NAT; AT3G19002) overlapping DOXC21-A (AT3G19000) in *Arabidopsis thaliana*. The analyses include sequence analysis, RNA secondary structure prediction, reanalysis of publicly available small RNA sequencing datasets, target prediction, functional enrichment, and phylogenetic analysis of DOXC21-related proteins.

## Sequence and analysis data

- `NAT_extended.fa` — Extended sequence containing the AT3G19002 NAT region used for structural analysis.
- `nat_mirna_like_candidates.fasta` — Candidate NAT-derived small RNA sequences.
- `DOXC21_590_sequences_MAFFT_alignment.fasta` — Multiple sequence alignment of DOXC21-related proteins used for phylogenetic analysis.
- `NAT_noncano_candidates.xlsx` — Filtered list of candidate NAT-derived non-canonical small RNA precursor regions.
- `NAT_siRNA_psRNATarget_full_predictions.xlsx` — Full psRNATarget output for predicted targets of NAT-derived small RNA candidates.
- `GO_enrichment_results.csv` — Gene Ontology enrichment results for predicted target genes.

## Custom scripts

- `detect_noncano_miRNA.py` — Identification and filtering of candidate non-canonical small RNA precursor structures.
- `generate_RNA_structures_VARNA.py` — Automated visualization of predicted RNA secondary structures using VARNA.

## Small RNA-seq reanalysis

Publicly available *Arabidopsis thaliana* small RNA sequencing datasets from two independent NCBI SRA BioProjects were reanalyzed to investigate small RNA reads associated with the AT3G19002 NAT locus.

### PRJNA607881 — *Phytophthora capsici* infection time course

The `PRJNA607881_Pcapsici` directory contains files generated during reanalysis of NCBI SRA BioProject PRJNA607881.

The dataset comprises 15 *A. thaliana* leaf small RNA-seq libraries (SRR11128435–SRR11128449) representing a *Phytophthora capsici* infection time course (0, 3, 6, 12, and 24 h post-inoculation; three biological replicates per time point).

The analysis quantified small RNA reads mapping to the AT3G19002 NAT region extended by 500 bp on both sides and examined their read-length distribution. No strand-specific filtering was applied in this analysis.

#### Files

- `AT3G19002_sRNAseq_mapping_pipeline.sh` — Script used to quantify reads overlapping the analyzed region and determine their read-length distribution.
- `nat.bed` — BED file defining the AT3G19002 NAT region extended by 500 bp on both sides.
- `AT3G19002_sRNAseq_read_counts.tsv` — Number of reads overlapping the analyzed region for individual libraries.
- `AT3G19002_sRNAseq_read_length_distribution.tsv` — Read-length distribution of reads mapping to the analyzed region for individual libraries.
- `AT3G19002_sRNAseq_read_length_summary.tsv` — Aggregated read-length summary across the analyzed libraries.
The mapping script expects preprocessed, genome-aligned and coordinate-sorted BAM files matching the pattern `trimmed_*_sorted.bam`. These intermediate files are not included in the repository.

### PRJNA1236841 — DAB-associated dataset

The `PRJNA1236841_DAB` directory contains files generated during reanalysis of NCBI SRA BioProject PRJNA1236841.

The analyzed ncRNA-seq libraries were derived from *A. thaliana* seedlings representing genotypes included in the DAB-associated dataset.

The analysis examined small RNA reads mapping to the AT3G19002 NAT region, including read abundance, read-length distribution, and strand-resolved mapping.

#### Files

- `AT3G19002_sRNAseq_read_counts.tsv` — Number of reads mapping to the analyzed AT3G19002 NAT region for individual libraries.
- `AT3G19002_sRNAseq_read_length_distribution.tsv` — Read-length distribution of reads mapping to the analyzed NAT region.
- `AT3G19002_strand_specific_counting.sh` — Script used for post hoc strand-resolved counting of NAT-region-mapping reads based on SAM flags.
- `AT3G19002_strand_specific_read_counts.tsv` — Plus- and minus-strand read counts for individual libraries.
- `NAT_region.bed` — BED file defining the AT3G19002 NAT region used for locus and strand-resolved analyses.
The strand-specific counting script expects preprocessed, genome-aligned and coordinate-sorted BAM files matching the pattern `*_sorted.bam`. These intermediate files are not included in the repository.

## Raw sequencing data

Raw sequencing reads are not redistributed in this repository and are publicly available from the NCBI Sequence Read Archive under:

- PRJNA607881
- PRJNA1236841

## Software

Custom analyses were performed using Python and standard scientific libraries, including pandas and matplotlib. Additional tools used in the analyses included SAMtools, Bowtie, ViennaRNA/RNAfold, VARNA, and R.

Detailed preprocessing procedures, software versions, and analysis parameters are described in the associated manuscript.
