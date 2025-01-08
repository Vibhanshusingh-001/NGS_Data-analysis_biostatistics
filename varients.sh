#!/bin/bash

# Define input files and reference genome
BAM_FILE="fixed_1.bam"
REF_GENOME="GCA_000001405.29_GRCh38.p14_genomic.fna"
VCF_OUTPUT="variants.vcf"
FILTERED_VCF="filtered_variants.vcf"

# Step 1: Generate Pileup File
samtools mpileup -f "$REF_GENOME" "$BAM_FILE" > "sample.pileup"

# Step 2: Call Variants with bcftools
bcftools call -c -v -o "$VCF_OUTPUT" "sample.pileup"

# Step 3: Filter Variants (e.g., minimum depth of 10, QUAL > 20)
bcftools filter -i 'QUAL>20 && DP>10' -o "$FILTERED_VCF" "$VCF_OUTPUT"

echo "Variant calling and filtering completed. Results saved in $FILTERED_VCF"
