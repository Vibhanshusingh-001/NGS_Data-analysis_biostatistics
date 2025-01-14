import vcfpy
import pandas as pd

vcf_file = "fixed_1.vcf"

variant_data = []

vcf_reader = vcfpy.Reader(open(vcf_file, 'r'))
for record in vcf_reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = ','.join([str(a) for a in record.ALT])
    qual = record.QUAL
    dp = record.INFO.get('DP', 0)  # Read Depth
    af = record.INFO.get('AF', [0])[0]  # Variant Allele Frequency (if available)

    # Calculate metrics (VAF, transitions/transversions)
    ti_tv = 'Transition' if {ref, alt}.issubset({'A', 'G', 'C', 'T'}) and \
        (ref + alt in ['AG', 'GA', 'CT', 'TC']) else 'Transversion'

    variant_data.append([chrom, pos, ref, alt, qual, dp, af, ti_tv])

df = pd.DataFrame(variant_data, columns=['CHROM', 'POS', 'REF', 'ALT', 'QUAL', 'DP', 'AF', 'Ti/Tv'])
df.to_csv('variant_metrics.csv', index=False)

print("Metrics saved to 'variant_metrics.csv'")
