#!/bin/bash
set -euo pipefail

# ============================================================
# AT3G19002 strand-specific small RNA read counting
# Dataset: NCBI SRA BioProject PRJNA1236841
#
# This script counts reads mapping to the plus and minus strands
# within the AT3G19002 NAT genomic region.
#
# Input:
#   NAT_region.bed
#   *_sorted.bam
#
# Output:
#   AT3G19002_strand_specific_read_counts.tsv
# ============================================================

echo -e "sample\tplus\tminus" > AT3G19002_strand_specific_read_counts.tsv

for bam in *_sorted.bam; do
    plus=$(samtools view -c -F 16 -L NAT_region.bed "$bam")
    minus=$(samtools view -c -f 16 -L NAT_region.bed "$bam")

    echo -e "$bam\t$plus\t$minus" >> AT3G19002_strand_specific_read_counts.tsv
done

echo "Done."
echo "Strand-specific read counts: AT3G19002_strand_specific_read_counts.tsv"