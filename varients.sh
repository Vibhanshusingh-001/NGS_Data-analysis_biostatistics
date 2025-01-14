
BAM_FILE="fixed_1.bam"
REF_GENOME="GCA_000001405.29_GRCh38.p14_genomic.fna"
VCF_OUTPUT="variants.vcf"
FILTERED_VCF="filtered_variants.vcf"

#  Generate Pileup File
samtools mpileup -f "$REF_GENOME" "$BAM_FILE" > "sample.pileup"

#  Call Variants with bcftools
bcftools call -c -v -o "$VCF_OUTPUT" "sample.pileup"

#   Filter Variants 
bcftools filter -i 'QUAL>20 && DP>10' -o "$FILTERED_VCF" "$VCF_OUTPUT"

echo "Variant calling and filtering completed. Results saved in $FILTERED_VCF"
