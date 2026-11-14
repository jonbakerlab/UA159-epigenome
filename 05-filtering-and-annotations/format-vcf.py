import os
import pandas as pd
import argparse

# %% Argparse preamble
parser = argparse.ArgumentParser(description="Convert the filtered modkit table TSV files (derived from the bedMethyl file) into a VCF compatible format.")
parser.add_argument("input1", help="Provide the path to a filtered modkit table TSV file.")
args = parser.parse_args()

# %% Import the TSV file into a DataFrames
vcf_format = pd.read_csv(args.input1, sep='\t', usecols=['Contig', 'Position'])

# rename the first two columns according to the format of a VCF
vcf_format.columns.values[0] = '#CHROM'
vcf_format.columns.values[1] = 'POS'

# Add empty columns as placeholders to format the VCF file
for column in ['ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE']:
    vcf_format[column] = None

# Make a new directory (if not made already) and export the VCF DataFrame into it
os.makedirs('snippy', exist_ok=True)
vcf_format.to_csv(os.path.join('snippy', os.path.splitext(os.path.basename(args.input1))[0] + '.vcf'), sep='\t', index=False)