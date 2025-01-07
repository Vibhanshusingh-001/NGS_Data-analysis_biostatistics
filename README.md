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
	
	# Columns for single CpG coverage
	cpg_columns = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	# Initialize a results dictionary to store metrics for each CpG column
	results = []
	
	# Loop through each CpG column and calculate statistics
	for cpg in cpg_columns:
	    # Group by Tissue and calculate Median, Mean, StdDev, and CV for the current CpG column
	    grouped = df.groupby('Tissue')[cpg].agg(
	        Median='median',
	        Mean='mean',
	        StdDev='std'
	    ).reset_index()
	    
	    # Add CV (Coefficient of Variation) to the grouped data
	    grouped['CV'] = grouped['StdDev'] / grouped['Mean']
	    
	    # Add a column to indicate which CpG column this analysis belongs to
	    grouped['CpG_Column'] = cpg
	    
	    # Append the result for this CpG column to the results list
	    results.append(grouped)
	
	# Concatenate all results into a single DataFrame
	final_results = pd.concat(results, ignore_index=True)
	
	# Save the results to a CSV file or print them
	final_results.to_csv('single_cpg_coverage_statistics.csv', index=False)
	print(final_results)



 ### Generate plots summarizing the coverage statistics
	import pandas as pd
	import seaborn as sns
	import matplotlib.pyplot as plt
	
	# Load your results 
	results = pd.read_csv('single_cpg_coverage_statistics.csv')
	
	# Plot settings
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
	
	# Adjust layout and save the figure
	plt.tight_layout()
	plt.savefig('coverage_statistics_plots.png', dpi=300)
	plt.show()

 ![coverage_statistics_plots](https://github.com/user-attachments/assets/1908bc71-d701-4f50-88d4-207e02d01292)



 # Biomarker Identification
 ### Identify PMPs with high specificity for tissue differentiation, minimizing false positives for Tissue #1 while allowing some false negatives. Use statistical or machine learning approaches to assign confidence (e.g., p-values) to each PMP
 ## The analysis was conducted on a small subset of the data rather than the entire dataset due to system limitations.
	 #!/usr/bin/env python3
	
	import pandas as pd
	from scipy.stats import fisher_exact
	from statsmodels.stats.multitest import multipletests
	
	# Load dataset in chunks
	chunk_size = 50
	file_path = "PupilBioTest.csv"  # Replace with your file
	methylation_cols = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	# Initialize an output file
	output_file = "statistical_results_temp.csv"
	with open(output_file, 'w') as f:
	    f.write("strand,CpG_Coordinates,Methylation_Status,Total_Count_T1,Total_Count_T2,P_Value\n")
	
	# Process dataset in chunks
	for chunk in pd.read_csv(file_path, chunksize=chunk_size):
	    # Melt the chunk
	    melted_chunk = chunk.melt(
	        id_vars=['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'],
	        value_vars=methylation_cols,
	        var_name='Methylation_Status',
	        value_name='Count'
	    )
	    
	    # Group by tissue, strand, coordinates, and methylation status
	    grouped = melted_chunk.groupby(['Tissue', 'strand', 'CpG_Coordinates', 'Methylation_Status']).agg(
	        Total_Count=('Count', 'sum')
	    ).reset_index()
	    
	    # Separate tissues
	    tissue1 = 'cfDNA'  # Replace with actual tissue name
	    tissue2 = 'Islet'  # Replace with actual tissue name
	    
	    tissue1_data = grouped[grouped['Tissue'] == tissue1]
	    tissue2_data = grouped[grouped['Tissue'] == tissue2]
	    
	    # Merge for comparison
	    comparison = pd.merge(
	        tissue1_data[['strand', 'CpG_Coordinates', 'Methylation_Status', 'Total_Count']],
	        tissue2_data[['strand', 'CpG_Coordinates', 'Methylation_Status', 'Total_Count']],
	        on=['strand', 'CpG_Coordinates', 'Methylation_Status'],
	        how='outer',
	        suffixes=('_T1', '_T2')
	    ).fillna(0)
	    
	    # Perform Fisher's Exact Test for each row
	    results = []
	    for _, row in comparison.iterrows():
	        contingency_table = [
	            [row['Total_Count_T1'], sum(tissue1_data['Total_Count']) - row['Total_Count_T1']],
	            [row['Total_Count_T2'], sum(tissue2_data['Total_Count']) - row['Total_Count_T2']],
	        ]
	        _, p_value = fisher_exact(contingency_table)
	        results.append((row['strand'], row['CpG_Coordinates'], row['Methylation_Status'],
	                        row['Total_Count_T1'], row['Total_Count_T2'], p_value))
	    
	    # Write intermediate results to a file
	    with open(output_file, 'a') as f:
	        for res in results:
	            f.write(",".join(map(str, res)) + "\n")
	
	# Load results for multiple testing correction
	df_results = pd.read_csv(output_file)
	valid_p_values = df_results['P_Value'].dropna()
	if len(valid_p_values) > 0:
	    df_results['Adjusted_P_Value'] = multipletests(valid_p_values, method='fdr_bh')[1]
	else:
	    df_results['Adjusted_P_Value'] = None
	
	# Filter significant PMPs
	significant_pmps = df_results[df_results['Adjusted_P_Value'] < 0.05]
	
	# Save final significant results
	significant_pmps.to_csv('statistical_significant_pmps.csv', index=False)
	print("Significant PMPs saved to 'statistical_significant_pmps.csv'")


 ### Calculate the mean variant read fraction (VRF) for each PMP in both tissues
	 #!/usr/bin/env python3
	
	import pandas as pd
	
	# Load dataset
	file_path = "PupilBioTest.csv"  # Replace with your actual file name
	methylation_cols = ['`000', '`001', '`010', '`011', '`100', '`101', '`110', '`111']
	
	# Load the data and clean column names
	df = pd.read_csv(file_path)
	df.columns = df.columns.str.strip()  # Remove leading/trailing spaces from column names
	
	# Verify columns exist
	required_columns = ['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'] + methylation_cols
	for col in required_columns:
	    if col not in df.columns:
	        raise ValueError(f"Column '{col}' not found in the dataset. Please check your input file.")
	
	# Melt the dataset to create a 'Methylation_Status' column
	df_melted = df.melt(
	    id_vars=['strand', 'CpG_Coordinates', 'Sample_ID', 'Replicate', 'Tissue'],
	    value_vars=methylation_cols,
	    var_name='Methylation_Status',
	    value_name='Count'
	)
	
	# Create Phased Methylation Pattern (PMP)
	df_melted['PMP'] = (
	    df_melted['strand'] + ':' + df_melted['CpG_Coordinates'] + ':' + df_melted['Methylation_Status']
	)
	
	# Group data by Tissue and PMP to calculate total counts
	grouped = df_melted.groupby(['Tissue', 'PMP']).agg(
	    Total_Count=('Count', 'sum')
	).reset_index()
	
	# Calculate the total counts per tissue
	grouped['Total_Tissue_Count'] = grouped.groupby('Tissue')['Total_Count'].transform('sum')
	
	# Calculate Variant Read Fraction (VRF)
	grouped['VRF'] = grouped['Total_Count'] / grouped['Total_Tissue_Count']
	
	# Save the VRF results to a file
	output_file = "vrf_results.csv"
	grouped.to_csv(output_file, index=False)
	
	print(f"VRF results saved to {output_file}")
	
		




