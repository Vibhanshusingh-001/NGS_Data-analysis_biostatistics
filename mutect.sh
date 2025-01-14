
input_dir="input"                
output_dir="output_mutect"               
reference="GCA_000001405.29_GRCh38.p14_genomic.fna"          

mkdir -p "$output_dir"

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

        # Check for error
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
