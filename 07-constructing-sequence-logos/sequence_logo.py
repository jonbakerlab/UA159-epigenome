import os
import pandas as pd
import re
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import logomaker
import matplotlib.pyplot as plt

# %% --- Argeparse preamble --- %%

parser = argparse.ArgumentParser(description="Generates a sequence logo of all of the bases surrounding modification sites identified")
parser.add_argument("input1", help="Provide the path to the annotation table TSV file, e.g., annotated_table_6mA.tsv")
parser.add_argument("input2", help="Provide the path to the FASTA file (FA or FNA file okay)")
parser.add_argument("input3", help="Provide the strain name, e.g., UA159")
args = parser.parse_args()

# Import the TSV input as a DataFrame and only keep the 'Position' and 'Strand' column
tsv_filename = args.input1
df = pd.read_csv(tsv_filename, sep="\t", usecols=['Contig', 'Position', 'Strand'])

# Import the FASTA file path
fna_filename =  args.input2

# Strain name
strain = args.input3

# %% --- Extract the sequence around the each of the modification positions --- %%

# No. bases before and after modification position
window_size = 12 # no. bases before and after modification position

# Load the genome
records = SeqIO.to_dict(SeqIO.parse(fna_filename, "fasta"))

# For loop cycling through all the rows in the data frame
for row_index in range(0, df.shape[0]):

    # Extract the position for that row
    position = df.loc[row_index, 'Position']
    strand = df.loc[row_index, 'Strand']
    contig = df.loc[row_index, 'Contig']

    record = records[contig]

    # Locate the start and end positions for the given window
    start = max(0, position - window_size)
    end = position + window_size - 1

    # If the strand is '-'...
    if strand == '-':
        window_seq = str(record.seq[start:end].reverse_complement())
    else:
        window_seq = str(record.seq[start:end])

    # Insert sequence in the 'Sequence' column of the DataFrame
    df.loc[row_index, 'Sequence'] = window_seq

# Filter out all the Data Frame such that all of the positions have the same sequence string length
sequence_df = df[(df["Sequence"].apply(len) == ((window_size * 2) - 1) )].reset_index(drop=True)

# %% --- Create the sequence logo of every methylation site --- %%
sequences = sequence_df['Sequence'].values

# Create a count matrix
counts_df = logomaker.alignment_to_matrix(sequences)

# Plot title
title = "Sequence Logo: " + strain

# Plot logo
plt.figure(figsize=(10, 4))
logo = logomaker.Logo(counts_df)
plt.title(title)
plt.savefig(strain + "_sequence_logo.pdf")
print(f"✅ Done! {title} has been successfully generated.")
