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

