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

