# NGS_Data-analysis_biostatistics
Analyzing complex datasets
### Task 1: Data Handling and Statistical Analysis

#### Objective: 
Handle complex data and apply statistical methods effectively

#### Background: 
CpG methylation is an epigenetic marker that varies across tissue types. However, the methylation status of a single CpG site is unreliable as a biomarker due to
errors introduced by bisulfite sequencing, sampling techniques, and biological variability.

#### Definition: 
Phased Methylation Pattern (PMP) is a unique set of coordinates that includes the DNA strand (‘f’ for forward (+) or ‘r’ for reverse (-)), the relative positions of three CpG
sites on the same strand (e.g., x:y:z), and their methylation status (e.g., ‘000’ for all unmethylated or ‘111’ for all methylated). It represents a combined epigenetic signature 

#### Hypothesis: 
Phased methylation patterns (PMPs) can act as reliable biomarkers to differentiate tissue types, providing higher specificity compared to individual CpG sites.

#### Dataset: 
The dataset (Link to Data) summarizes phased methylation patterns from NGS results across two tissues. Key columns include:

	#### Strand:Indicates the DNA strand (‘f’ or ‘r’).
	#### CpG Coordinates: Relative positions of three CpG sites (x:y:z).
	#### Methylation Status: Eight possible patterns (‘000’ to ‘111’).
	#### Sample ID: Unique identifier for each sample.
	#### Replicate: Indicates technical replicates.
	#### Tissue: Tissue type (Tissue #1 or Tissue #2).

 # Coverage Analysis
 ### Calculate the median and coefficient of variation (CV) for single CpG coverage in each tissue
	import pandas as pd		
	df = pd.read_csv('PupilBioTest.csv')
	
	cpg_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	# list to store metrics for each CpG column
	results = []
	
	# loop through each CpG column 
	for cpg in cpg_columns:
	    # group by tissue
	    grouped = df.groupby('Tissue')[cpg].agg(
	        Median='median',
	        Mean='mean',
	        StdDev='std'
	    ).reset_index()
	    
	 
	    grouped['CV'] = grouped['StdDev'] / grouped['Mean']
	    
	    # add a column to indicate which CpG column this analysis belongs to
	    grouped['CpG_Column'] = cpg
	    
	    # append the result for this CpG column to the results list
	    results.append(grouped)
	
	# merge all results 
	final_results = pd.concat(results, ignore_index=True)
	

	final_results.to_csv('single_cpg_coverage_statistics.csv', index=False)
	print(final_results)



 ### Generate plots summarizing the coverage statistics
	import pandas as pd
	import seaborn as sns
	import matplotlib.pyplot as plt
	
	results = pd.read_csv('single_cpg_coverage_statistics.csv')
	
	# Plot 
	sns.set(style="whitegrid")
	plt.figure(figsize=(16, 10))
	
	# 1. Median plot
	plt.subplot(2, 2, 1)
	sns.barplot(data=results, x='CpG_Column', y='Median', hue='Tissue', palette='muted')
	plt.title('Median Coverage per CpG Column')
	plt.xlabel('CpG Column')
	plt.ylabel('Median Coverage')
	
	# 2. Mean plot
	plt.subplot(2, 2, 2)
	sns.barplot(data=results, x='CpG_Column', y='Mean', hue='Tissue', palette='pastel')
	plt.title('Mean Coverage per CpG Column')
	plt.xlabel('CpG Column')
	plt.ylabel('Mean Coverage')
	
	# 3. Standard Deviation plot
	plt.subplot(2, 2, 3)
	sns.barplot(data=results, x='CpG_Column', y='StdDev', hue='Tissue', palette='dark')
	plt.title('Standard Deviation of Coverage per CpG coverage Column')
	plt.xlabel('CpG Column')
	plt.ylabel('Standard Deviation')
	
	# 4. CV plot
	plt.subplot(2, 2, 4)
	sns.barplot(data=results, x='CpG_Column', y='CV', hue='Tissue', palette='colorblind')
	plt.title('Coefficient of Variation (CV) per CpG coverage Column')
	plt.xlabel('CpG Column')
	plt.ylabel('Coefficient of Variation')
	
	plt.tight_layout()
	plt.savefig('coverage_statistics_plots.png', dpi=300)
	plt.show()

 ![coverage_statistics_plots](https://github.com/user-attachments/assets/1908bc71-d701-4f50-88d4-207e02d01292)



 # Biomarker Identification
 ### Identify PMPs with high specificity for tissue differentiation, minimizing false positives for Tissue #1 while allowing some false negatives. Use statistical or machine learning approaches to assign confidence (e.g., p-values) to each PMP

	 
	
		
	
	import pandas as pd
	from scipy.stats import fisher_exact
	from statsmodels.stats.multitest import multipletests
	chunk_size = 500
	file_path = "PupilBioTest.csv"  
	methylation_cols = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	# initialize an empty list to store melted chunks
	melted_chunks = []
	
	# process the dataset in chunks
	for chunk in pd.read_csv(file_path, chunksize=chunk_size):
	    
	    melted_chunk = chunk.melt(
	        id_vars=['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'],
	        value_vars=methylation_cols,
	        var_name='Methylation_Status',
	        value_name='Count'
	    )
	    melted_chunks.append(melted_chunk)
	
	# combine into a single DataFrame
	df_melted = pd.concat(melted_chunks, ignore_index=True)
	
	# group data by tissue, strand, coordinates, and methylation status
	grouped = df_melted.groupby(['Tissue', 'strand', 'CpG_Coordinates', 'Methylation_Status']).agg(
	    Total_Count=('Count', 'sum')
	).reset_index()
	
	# Normalize counts by total counts per tissue
	grouped['Normalized_Count'] = grouped['Total_Count'] / grouped.groupby('Tissue')['Total_Count'].transform('sum')
	
	#  comparison table
	tissue1 = 'cfDNA'  
	tissue2 = 'islet'  
	tissue1_data = grouped[grouped['Tissue'] == tissue1]
	tissue2_data = grouped[grouped['Tissue'] == tissue2]
	
	comparison = pd.merge(
	    tissue1_data[['strand', 'CpG_Coordinates', 'Methylation_Status', 'Total_Count']],
	    tissue2_data[['strand', 'CpG_Coordinates', 'Methylation_Status', 'Total_Count']],
	    on=['strand', 'CpG_Coordinates', 'Methylation_Status'],
	    how='outer',
	    suffixes=('_T1', '_T2')
	).fillna(0)
	
	# Fisher's Exact Test for each PMP
	p_values = []
	for _, row in comparison.iterrows():
	    contingency_table = [
	        [row['Total_Count_T1'], sum(tissue1_data['Total_Count']) - row['Total_Count_T1']],
	        [row['Total_Count_T2'], sum(tissue2_data['Total_Count']) - row['Total_Count_T2']],
	    ]
	    _, p = fisher_exact(contingency_table)
	    p_values.append(p)
	
	comparison['P_Value'] = p_values
	
	# adusted p-values for multiple testing with error handling
	if not comparison['P_Value'].isnull().all() and len(comparison['P_Value']) > 0:
	    # Remove NaN P-values and apply multiple testing correction
	    valid_p_values = comparison['P_Value'].dropna()
	    if len(valid_p_values) > 0:
	        try:
	            comparison['Adjusted_P_Value'] = multipletests(valid_p_values, method='fdr_bh')[1]
	        except Exception as e:
	            print(f"Error during multiple testing correction: {e}")
	            comparison['Adjusted_P_Value'] = None
	    else:
	        print("Error: No valid P-values found for multiple testing correction.")
	        comparison['Adjusted_P_Value'] = None
	else:
	    print("Error: No valid P-values found for multiple testing correction.")
	    comparison['Adjusted_P_Value'] = None
	
	# Filter PMPs
	significant_pmps = comparison[comparison['Adjusted_P_Value'] < 0.05]
	
	significant_pmps.to_csv('statistical_significant_pmps.csv', index=False)
	print("Significant PMPs based on statistical approach saved to 'statistical_significant_pmps.csv'")


 ### Calculate the mean variant read fraction (VRF) for each PMP in both tissues
	
	
	import pandas as pd
	
	file_path = "PupilBioTest.csv"  # Replace with your actual file name
	methylation_cols = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	df = pd.read_csv(file_path)
	df.columns = df.columns.str.strip()  
	
	required_columns = ['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'] + methylation_cols
	for col in required_columns:
	    if col not in df.columns:
	        raise ValueError(f"Column '{col}' not found in the dataset. Please check your input file.")
	
	# melt the dataset to create a 'Methylation_Status' column
	df_melted = df.melt(
	    id_vars=['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'],
	    value_vars=methylation_cols,
	    var_name='Methylation_Status',
	    value_name='Count'
	)
	
	# Create PMP
	df_melted['PMP'] = (
	    df_melted['strand'] + ':' + df_melted['CpG_Coordinates'] + ':' + df_melted['Methylation_Status']
	)
	
	# group data by Tissue and PMP to calculate total counts
	grouped = df_melted.groupby(['Tissue', 'PMP']).agg(
	    Total_Count=('Count', 'sum')
	).reset_index()
	
	# Calculate the total counts for each tissue
	grouped['Total_Tissue_Count'] = grouped.groupby('Tissue')['Total_Count'].transform('sum')
	
	# Calculate VRF
	grouped['VRF'] = grouped['Total_Count'] / grouped['Total_Tissue_Count']
	
	output_file = "vrf_results.csv"
	grouped.to_csv(output_file, index=False)
	
	print(f"VRF results saved to {output_file}")
	


 # Address the following questions 
 ### How does sequencing depth affect specificity confidence?
Sequencing depth directly affects specificity confidence by improving accuracy and reducing errors. Higher depth ensures reliable detection of true variants or methylation patterns, minimizes false positives/negatives, and increases statistical power.
#### High Sequencing Depth: Benefits
 	Increased Accuracy
 	Minimized False Positives/Negatives
 	Rare Variant Detection
 	Improved Statistical Power
 	Consistent Normalization
#### Low Sequencing Depth: Drawbacks
 	Higher Error Rates
 	Ambiguous Results
 	Reduced Sensitivitys
### For the top 10 PMPs, estimate the threshold of reads required to confidently call Tissue #2 at a sequencing depth of 1 million reads.
	
         import pandas as pd
	
	file_path = "statistical_significant_pmps.csv"  # Use the provided data file
	data = pd.read_csv(file_path)
	
	# Sort by Adjusted P-Value (column 'Adjusted_P_Value') to get the top 10 PMPs
	top_pmps = data.sort_values(by='Adjusted_P_Value').head(10)
	
	# Assuming a sequencing depth of 1 million reads for Tissue #2
	sequencing_depth_t2 = 1_000_000
	
	# Estimate the threshold of reads required for each PMP
	top_pmps['Threshold_Reads_T2'] = top_pmps['Total_Count_T2'] / top_pmps['Total_Count_T2'].sum() * sequencing_depth_t2
	
	print(top_pmps[['CpG_Coordinates', 'Methylation_Status', 'Total_Count_T2', 'Threshold_Reads_T2']])

 ![Screenshot from 2025-01-14 15-16-40](https://github.com/user-attachments/assets/0c28fed0-d93d-482d-8d07-3a05d15a6985)

	
		
		
	
	
	
	

















### Validate the hypothesis by comparing the specificity of the top 10 PMPs against individual CpG sites


	import pandas as pd
	from scipy.stats import fisher_exact
	from statsmodels.stats.multitest import multipletests
	import matplotlib.pyplot as plt
	
	file_path = "PupilBioTest.csv"
	
	methylation_cols = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	# Initialize result storage
	pmp_results = []
	cpg_results = []
	
	# Process data in chunks
	chunk_size = 200
	for chunk in pd.read_csv(file_path, chunksize=chunk_size):
	    melted_chunk = chunk.melt(
	        id_vars=['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'],
	        value_vars=methylation_cols,
	        var_name='Methylation_Status',
	        value_name='Count'
	    )
	
	    grouped = melted_chunk.groupby(['Tissue', 'strand', 'CpG_Coordinates', 'Methylation_Status']).agg(
	        Total_Count=('Count', 'sum')
	    ).reset_index()
	
	    tissue1 = 'cfDNA'
	    tissue2 = 'Islet'
	
	    tissue1_data = grouped[grouped['Tissue'] == tissue1]
	    tissue2_data = grouped[grouped['Tissue'] == tissue2]
	
	    # PMPs Analysis
	    comparison = pd.merge(
	        tissue1_data[['strand', 'CpG_Coordinates', 'Methylation_Status', 'Total_Count']],
	        tissue2_data[['strand', 'CpG_Coordinates', 'Methylation_Status', 'Total_Count']],
	        on=['strand', 'CpG_Coordinates', 'Methylation_Status'],
	        how='outer',
	        suffixes=('_T1', '_T2')
	    ).fillna(0)
	
	    for _, row in comparison.iterrows():
	        contingency_table = [
	            [int(row['Total_Count_T1']), int(sum(tissue1_data['Total_Count']) - row['Total_Count_T1'])],
	            [int(row['Total_Count_T2']), int(sum(tissue2_data['Total_Count']) - row['Total_Count_T2'])],
	        ]
	        if any(x < 0 for contingency_row in contingency_table for x in contingency_row) or sum(sum(contingency_table, [])) == 0:
	            continue
	
	        _, p_value = fisher_exact(contingency_table)
	
	        true_negatives = contingency_table[1][1]
	        total_negatives = contingency_table[1][0] + contingency_table[1][1]
	        specificity = true_negatives / total_negatives if total_negatives > 0 else 0
	
	        pmp_results.append((row['CpG_Coordinates'], row['Methylation_Status'], p_value, specificity))
	
	    # CpG Site Analysis
	    individual_cpg = melted_chunk.groupby(['CpG_Coordinates', 'Tissue']).agg(
	        Total_Count=('Count', 'sum')
	    ).reset_index()
	
	    for cpg, group in individual_cpg.groupby('CpG_Coordinates'):
	        tissue_counts = group.set_index('Tissue')['Total_Count'].to_dict()
	        count_t1 = tissue_counts.get(tissue1, 0)
	        count_t2 = tissue_counts.get(tissue2, 0)
	        contingency_table = [
	            [count_t1, sum(group['Total_Count']) - count_t1],
	            [count_t2, sum(group['Total_Count']) - count_t2],
	        ]
	        if any(x < 0 for contingency_row in contingency_table for x in contingency_row) or sum(sum(contingency_table, [])) == 0:
	            continue
	
	        _, p_value = fisher_exact(contingency_table)
	
	        true_negatives = contingency_table[1][1]
	        total_negatives = contingency_table[1][0] + contingency_table[1][1]
	        specificity = true_negatives / total_negatives if total_negatives > 0 else 0
	
	        cpg_results.append((cpg, p_value, specificity))
	
	# Convert results to DataFrame
	pmp_df = pd.DataFrame(pmp_results, columns=['CpG_Coordinates', 'Methylation_Status', 'P_Value', 'Specificity'])
	cpg_df = pd.DataFrame(cpg_results, columns=['CpG_Coordinates', 'P_Value', 'Specificity'])
	
	# Adjust p-values
	pmp_df['Adjusted_P_Value'] = multipletests(pmp_df['P_Value'], method='fdr_bh')[1]
	cpg_df['Adjusted_P_Value'] = multipletests(cpg_df['P_Value'], method='fdr_bh')[1]
	
	# Filter top 10 PMPs and CpG sites
	top_pmps = pmp_df.sort_values(by='Specificity', ascending=False).head(10)
	top_cpgs = cpg_df.sort_values(by='Specificity', ascending=False).head(10)
	
	# Compare Specificity Distributions
	plt.boxplot([top_pmps['Specificity'], top_cpgs['Specificity']], labels=['PMPs', 'CpG Sites'])
	plt.ylabel('Specificity')
	plt.title('Specificity Comparison: PMPs vs CpG Sites')
	plt.show()
	
	top_pmps.to_csv("top_10_pmps.csv", index=False)
	top_cpgs.to_csv("top_10_cpg_sites.csv", index=False)
	
	print("Top 10 PMPs and CpG sites saved. Specificity comparison completed.")
	
The analysis shows that individual CpG sites have perfect specificity (1.0), while PMPs have very low specificity (0.0). Although CpG sites demonstrate high reliability in this study,their use as biomarkers is limited by errors from bisulfite sequencing, sampling techniques, and biological variability. PMPs aim to overcome these issues but currently fail to achieve sufficient specificity. Future work should focus on improving PMP design or integrating CpG site reliability with PMP robustness for better epigenetic biomarkers.	
	
	
	

# Task 2: NGS Data Analysis
## Objective: Evaluate candidates’ ability to process and analyze raw sequencing data.
### 1) Quality Control 
### (a)Perform quality checks using tools like FastQC and summarize quality metrics (e.g., sequence counts, per-base quality, read duplication levels)
	# command to run fastqc
	fastqc <filename.fastq> -o <output_directory>

#### Per Base Sequence Quality
The graph in this section shows the quality scores for each base position in the sequencing reads. The quality score is shown by the y-axis, while the read’s location is represented by the x-axis. The graph’s line or bars show how the quality scores vary with the length of the reads.
![image](https://github.com/user-attachments/assets/4e121559-4399-463e-82f9-93f9c529221e)
### Sequence Duplication level
This metric reveals if specific sequences are overrepresented as a result of biases in the library preparation or sequencing process by giving insights into the degree of duplication within a sequencing dataset. Increased duplication levels may indicate that particular sequences were amplified or enriched more frequently than others during library preparation. The results of subsequent investigations, including variant calling or differential gene expression analysis, may become skewed as a result of biassed representation.
![image](https://github.com/user-attachments/assets/f8139186-d3b7-4125-a2db-a43511714464)
### Per Sequence Quality Scores
The graph in this section shows the average quality score distribution over all sequences. The number or percentage of sequences with that average quality score is shown on the y-axis, which is parallel to the x-axis and displays quality scores.
![image](https://github.com/user-attachments/assets/a8b3bf4b-2115-4327-b6f4-97dfe996e58a)
### 2) Alignment and Mutation Calling
#### (a)Align the samples to the human genome using tools like Bowtie2 or BWA.
Alignment score is greater than 95% 

	REF="GCA_000001405.29_GRCh38.p14_genomic.fna"
	
	SAMPLES=("PA220KH-lib09-P19-Tumor_S2_L001" "PA221MH-lib09-P19-Norm_S1_L001")
	
	for SAMPLE in "${SAMPLES[@]}"; do
	    # Define file names for paired-end reads
	    READ1="${SAMPLE}_R1_001.fastq.gz"
	    READ2="${SAMPLE}_R2_001.fastq.gz"
	    SAM_OUTPUT="${SAMPLE}.sam"
	    BAM_OUTPUT="${SAMPLE}.bam"
	    SORTED_BAM_OUTPUT="${SAMPLE}_sorted.bam"
	
	    
	    if [[ -f "$READ1" && -f "$READ2" ]]; then
	        echo "Processing $SAMPLE..."
	
	        # align paired-end reads using BWA MEM
	        bwa mem "$REF" "$READ1" "$READ2" > "$SAM_OUTPUT"
	
	        #  SAM to BAM
	        samtools view -Sb "$SAM_OUTPUT" > "$BAM_OUTPUT"
	
	        samtools sort "$BAM_OUTPUT" -o "$SORTED_BAM_OUTPUT"
	
	        samtools index "$SORTED_BAM_OUTPUT"
	
	        echo "Finished processing $SAMPLE."
	    else
	        echo "Warning: Input files for $SAMPLE not found."
	    fi
	done
	
#### (b)identify somatic mutations present in the cancer sample but absent in the normal tissue.
#### (i)Benchmark Software: Use established tools such as Mutect2,Strelka2, or VarScan2 for somatic mutation identification and background mutation estimation. 
## mutect2
	
	
	
	input_dir="input"                      
	output_dir="output_mutect"              
	reference="GCA_000001405.29_GRCh38.p14_genomic.fna" # Reference genome file
	
	mkdir -p "$output_dir"
	
	for bam_file in "$input_dir"/*.bam; do
	    if [ -f "$bam_file" ]; then
	        # Get the base name of the BAM file
	        base_name=$(basename "$bam_file" .bam)
	        output_vcf="$output_dir/${base_name}.vcf"
	        
	        echo "Processing $bam_file..."
	        
	        #  Mutect2 command
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

#### (ii)Custom Code Development: Write your own scripts,leveraging tools like Samtools, bcftools, or Python/Rlibraries, to perform mutation detection and calculate therequired metrics.

### levraging tools for mutation detection
	
	
	
	BAM_FILE="fixed_1.bam"
	REF_GENOME="GCA_000001405.29_GRCh38.p14_genomic.fna"
	VCF_OUTPUT="variants.vcf"
	FILTERED_VCF="filtered_variants.vcf"
	
	#  generate Pileup File
	samtools mpileup -f "$REF_GENOME" "$BAM_FILE" > "sample.pileup"
	
	# call Variants with bcftools
	bcftools call -c -v -o "$VCF_OUTPUT" "sample.pileup"
	
	#  filter Variants 
	bcftools filter -i 'QUAL>20 && DP>10' -o "$FILTERED_VCF" "$VCF_OUTPUT"
	
	echo "Variant calling and filtering completed. Results saved in $FILTERED_VCF"

### calculation of the required metrics

	import vcfpy
	import pandas as pd
	
	vcf_file = "fixed_1.vcf"
	
	variant_data = []
	
	# Parse VCF file
	vcf_reader = vcfpy.Reader(open(vcf_file, 'r'))
	for record in vcf_reader:
	    chrom = record.CHROM
	    pos = record.POS
	    ref = record.REF
	    alt = ','.join([str(a) for a in record.ALT])
	    qual = record.QUAL
	    dp = record.INFO.get('DP', 0)  # Read Depth
	    af = record.INFO.get('AF', [0])[0]  # Variant Allele Frequency (if available)
	
	    ti_tv = 'Transition' if {ref, alt}.issubset({'A', 'G', 'C', 'T'}) and \
	        (ref + alt in ['AG', 'GA', 'CT', 'TC']) else 'Transversion'
	
	    # Append to results
	    variant_data.append([chrom, pos, ref, alt, qual, dp, af, ti_tv])
	
	df = pd.DataFrame(variant_data, columns=['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'DP', 'AF', 'Ti/Tv'])
	df.to_csv('variant_metrics.csv', index=False)
	


### visualisation 
	
	library(ggplot2)
	
	variant_data <- read.csv("variant_metrics.csv")
	
	ti_count <- sum(variant_data$Ti.Tv == "Transition")
	tv_count <- sum(variant_data$Ti.Tv == "Transversion")
	ti_tv_ratio <- ti_count / tv_count
	cat("Transition/Transversion Ratio:", ti_tv_ratio, "\n")
	# Plot VAF distribution
	ggplot(variant_data, aes(x = AF)) +
	  geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
	  labs(title = "Variant Allele Frequency Distribution", x = "VAF", y = "Count")
   ![Screenshot from 2025-01-08 23-43-30](https://github.com/user-attachments/assets/50e8fcf7-e13b-4dd4-9a26-dca224622797)

			
#### (c) Use the normal tissue to calculate the median background mutation level.The background mutation level accounts for sequencing errors or biasesthat can mimic true mutations. Determine how many reads per million are required to confidently call a given mutation.
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
	
	vcf_file = "fixed_1.vcf"  # same for second sample 
	bam_file = "fixed_1.bam"  
	
	results = calculate_background_mutation_rate(vcf_file, bam_file)
	
	print("=== Background Mutation Analysis ===")
	print(f"Total Mutations: {results['Total Mutations']}")
	print(f"Total Reads: {results['Total Reads']}")
	print(f"Mutations Per Million Reads: {results['Mutations Per Million Reads']:.2f}")
	print(f"Confidence Threshold: {results['Confidence Threshold']:.2f}")
   ![Screenshot from 2025-01-14 21-52-34](https://github.com/user-attachments/assets/56ded55b-ffec-466e-bb07-e60165c42036)

		
	

