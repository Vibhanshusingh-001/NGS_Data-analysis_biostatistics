import pandas as pd

file_path = "statistical_significant_pmps.csv"  # Use the provided data file
data = pd.read_csv(file_path)

# sort by Adjusted P-Value (column 'Adjusted_P_Value') to get the top 10 PMPs
top_pmps = data.sort_values(by='Adjusted_P_Value').head(10)

# assuming a sequencing depth of 1 million reads for Tissue #2
sequencing_depth_t2 = 1_000_000

# estimate the threshold of reads required for each PMP
top_pmps['Threshold_Reads_T2'] = top_pmps['Total_Count_T2'] / top_pmps['Total_Count_T2'].sum() * sequencing_depth_t2

print(top_pmps[['CpG_Coordinates', 'Methylation_Status', 'Total_Count_T2', 'Threshold_Reads_T2']])
