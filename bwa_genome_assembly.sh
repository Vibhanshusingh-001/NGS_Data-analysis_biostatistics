#!/bin/bash

# Load BWA and Samtools modules (uncomment if needed)
# module load bwa
# module load samtools

# Define reference genome
REF="GCA_000001405.29_GRCh38.p14_genomic.fna"

# Define specific samples with their prefixes
SAMPLES=("PA220KH-lib09-P19-Tumor_S2_L001" "PA221MH-lib09-P19-Norm_S1_L001")

# Process each sample
for SAMPLE in "${SAMPLES[@]}"; do
    # Define file names for paired-end reads
    READ1="${SAMPLE}_R1_001.fastq.gz"
    READ2="${SAMPLE}_R2_001.fastq.gz"
    SAM_OUTPUT="${SAMPLE}.sam"
    BAM_OUTPUT="${SAMPLE}.bam"
    SORTED_BAM_OUTPUT="${SAMPLE}_sorted.bam"

    # Check if input files exist
    if [[ -f "$READ1" && -f "$READ2" ]]; then
        echo "Processing $SAMPLE..."

        # Align paired-end reads using BWA MEM
        bwa mem "$REF" "$READ1" "$READ2" > "$SAM_OUTPUT"

        # Convert SAM to BAM
        samtools view -Sb "$SAM_OUTPUT" > "$BAM_OUTPUT"

        # Sort the BAM file
        samtools sort "$BAM_OUTPUT" -o "$SORTED_BAM_OUTPUT"

        # Optionally index the sorted BAM file
        samtools index "$SORTED_BAM_OUTPUT"

        echo "Finished processing $SAMPLE."
    else
        echo "Warning: Input files for $SAMPLE not found."
    fi
done

