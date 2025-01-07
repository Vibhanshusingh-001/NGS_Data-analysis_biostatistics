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

