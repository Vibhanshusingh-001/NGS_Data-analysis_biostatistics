#!/bin/bash

# Directories and files
input_dir="input"                 # Directory containing BAM files
output_dir="output_mutect"               # Directory to save VCF files
reference="GCA_000001405.29_GRCh38.p14_genomic.fna"            # Reference genome file

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Loop through BAM files in the input directory
for bam_file in "$input_dir"/*.bam; do
    if [ -f "$bam_file" ]; then
        # Get the base name of the BAM file
        base_name=$(basename "$bam_file" .bam)
        output_vcf="$output_dir/${base_name}.vcf"
        
        echo "Processing $bam_file..."
        
        # Main Mutect2 command
        gatk Mutect2 \
            -R "$reference" \
            -I "$bam_file" \
            -O "$output_vcf" \

        # Check for errors
        if [ $? -eq 0 ]; then
            echo "Completed $bam_file -> $output_vcf"
        else
            echo "Error processing $bam_file" >&2
        fi
    else
        echo "No BAM files found in $input_dir."
    fi
done

echo "All samples processed."
