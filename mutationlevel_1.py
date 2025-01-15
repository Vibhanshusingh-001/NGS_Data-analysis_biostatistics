import numpy as np
import pysam
def calculate_background_mutation_rate(vcf_file, bam_file):
    """
    Calculate the median background mutation rate and required reads for confident mutation calls.
    Parameters:
    - vcf_file (str): fixed_1.vcf
    - bam_file (str): fixed_1.bam
    Returns:
    - dict: Median mutation rate, total mutations, total reads, mutations per million reads.
    """
    # Count variants from the VCF file
    with open(vcf_file, 'r') as vcf:
        variants = [line for line in vcf if not line.startswith("#")]
    total_mutations = len(variants) 
    # Calculate total reads from the BAM file
    bam = pysam.AlignmentFile(bam_file, "rb")
    total_reads = sum([stats.mapped for stats in bam.get_index_statistics()])
    bam.close() 
    #  Calculate mutations per million reads
    mutations_per_million_reads = (total_mutations / total_reads) * 1e6  
    # Calculate confidence threshold (optional scaling for higher confidence)
    confidence_threshold = mutations_per_million_reads + 1.5 * (np.median(np.abs(total_mutations - mutations_per_million_reads)))
    return {
        "Total Mutations": total_mutations,
        "Total Reads": total_reads,
        "Mutations Per Million Reads": mutations_per_million_reads,
        "Confidence Threshold": confidence_threshold,
    }
vcf_file = "fixed_1.vcf"  
bam_file = "fixed_1.bam"  
results = calculate_background_mutation_rate(vcf_file, bam_file)
print("=== Background Mutation Analysis ===")
print(f"Total Mutations: {results['Total Mutations']}")
print(f"Total Reads: {results['Total Reads']}")
print(f"Mutations Per Million Reads: {results['Mutations Per Million Reads']:.2f}")
print(f"Confidence Threshold: {results['Confidence Threshold']:.2f}")
