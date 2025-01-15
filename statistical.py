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
