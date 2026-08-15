#!/bin/bash
set -euo pipefail

# ============================================================
# AT3G19002 NAT locus small RNA-seq mapping analysis
# Dataset: NCBI SRA BioProject PRJNA607881
#
# This script quantifies small RNA reads overlapping the
# AT3G19002 NAT locus and calculates their length distribution.
# No strand-specific filtering is applied.
#
# Input:
#   nat.bed                  - AT3G19002 NAT locus coordinates
#   trimmed_*_sorted.bam     - preprocessed, genome-aligned,
#                              sorted sRNA-seq BAM files
#
# Output:
#   nat_counts/AT3G19002_sRNAseq_read_counts.tsv
#   nat_lengths/AT3G19002_sRNAseq_read_length_distribution.tsv
#   nat_sam/*_nat.sam
#
# Raw sequencing data: NCBI SRA BioProject PRJNA607881
# (SRR11128435-SRR11128449)
# ============================================================

# Input region
BED="nat.bed"

# Output folders
mkdir -p nat_counts
mkdir -p nat_lengths
mkdir -p nat_sam

# Summary files
echo -e "sample\treads_in_NAT" > nat_counts/AT3G19002_sRNAseq_read_counts.tsv
echo -e "sample\tlength\tcount" > nat_lengths/AT3G19002_sRNAseq_read_length_distribution.tsv

# Loop through all sorted BAM files
for bam in trimmed_*_sorted.bam; do
    sample=$(basename "$bam" _sorted.bam)

    echo "Processing $sample ..."

    # Count total reads in NAT region
    count=$(samtools view -L "$BED" "$bam" | wc -l | tr -d ' ')
    echo -e "${sample}\t${count}" >> nat_counts/AT3G19002_sRNAseq_read_counts.tsv

    # Save SAM reads from NAT region
    samtools view -L "$BED" "$bam" > "nat_sam/${sample}_nat.sam"

    # Count read lengths
    awk '{
        if($1 !~ /^@/) {
            print length($10)
        }
    }' "nat_sam/${sample}_nat.sam" | sort | uniq -c | awk -v s="$sample" '{
        print s "\t" $2 "\t" $1
    }' >> nat_lengths/AT3G19002_sRNAseq_read_length_distribution.tsv

done

echo "Done."
echo "Read counts: nat_counts/AT3G19002_sRNAseq_read_counts.tsv"
echo "Length distribution: nat_lengths/AT3G19002_sRNAseq_read_length_distribution.tsv"