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
